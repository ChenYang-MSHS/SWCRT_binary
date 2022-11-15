require(boot)
require(lme4)
require(expm)

calc.power <- function(m, method,
                       reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                       sigmas = c(1, 0.5, 0.75), 
                       prevalences = 0.5, 
                       I = 8, 
                       sig.level = 0.05){
  J = length(reg.coef) - 4
  cluster.index <- 1:I
  alloc <- rep(1:J, each = (I/J))
  S <- length(unique(alloc))
  T.period <- S + 1
  state.ini <- rep(0, I)
  SW.design <- sapply(1:T.period, function(s){
    state.ini[which(alloc < s)] <- 1
    return(state.ini)
  })
  W.list <- lapply(cluster.index, function(i){
    rep(SW.design[i,], each = m)
  })
  W <- unlist(W.list)
  e.list <- lapply(cluster.index, function(i){
    rep(0:S, each = m)
  })
  e <- unlist(e.list)
  cluster.label <- rep(cluster.index, (rep(m, I) * T.period))
  period.label <- unlist(lapply(rep(m, I), function(i){
    rep(1:T.period, each = i)
  }))
  individual.label <- unlist(lapply(rep(m, I), function(i){
    rep(1:i, T.period)
  }))
  sigma.alpha <- sigmas[1]
  sigma.nu <- sigmas[2]
  sigma.psi <- sigmas[3]
  I.matrix <- function(n) Matrix(diag(rep(1, n)),sparse = TRUE)
  J.matrix <- function(n) Matrix(matrix(1, nrow = n, ncol = n), sparse = TRUE)
  D.list <- lapply(rep(m, I), function(i){
    sigma.alpha^2 * J.matrix(i * T.period) + sigma.nu^2 * kronecker(I.matrix(T.period), J.matrix(i)) + sigma.psi^2 * kronecker(J.matrix(T.period), I.matrix(i))
  })
  X.list <- lapply(cluster.index, function(i){
    rep(rep(c(1,0), c(floor(m * prevalences), (m - floor(m * prevalences)))), T.period)
  })
  X <- unlist(X.list)
  M.list <- lapply(cluster.index, function(i){
    model.matrix(~ factor(e.list[[i]]) + W.list[[i]] * X.list[[i]])
  })
  M <- model.matrix(~ factor(e) + W * X)
  L.list <- lapply(M.list, function(i){
    u0 <- inv.logit(as.numeric(i %*% reg.coef))
    Matrix(diag(u0 * (1 - u0)), sparse = TRUE)
  })
  if(method == "AML"){
    Sigma.AML <- lapply(cluster.index, function(i){
      D.list[[i]] + solve(L.list[[i]])
    })
    method.list <- lapply(cluster.index, function(i){
      t(M.list[[i]]) %*% solve(Sigma.AML[[i]]) %*% M.list[[i]]
    })
    FI <- Reduce("+", method.list)
    var.interact <- solve(FI)[length(reg.coef), length(reg.coef)]
  }
  else if(method == "GEE"){
    A.list <- lapply(L.list, function(i){
      i - (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2) * i %*% i
    })
    Vtilde.list <- lapply(cluster.index, function(i){
      L.list[[i]] %*% D.list[[i]] %*% L.list[[i]] + A.list[[i]]
    })
    constant <- 16 * sqrt(3) / (15 * pi)
    cluster.constant <- 1/sqrt(1 + constant^2 * (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2))
    Ltilde.list <- lapply(M.list, function(i){
      u <- as.numeric(inv.logit(cluster.constant * i %*% reg.coef))
      diag(u * (1 - u))
    })
    Sigma.GEE.inv <- lapply(cluster.index, function(i){
      cluster.constant^2 * Ltilde.list[[i]] %*% solve(Vtilde.list[[i]]) %*% Ltilde.list[[i]]
    })
    method.list <- lapply(cluster.index, function(i){
      t(M.list[[i]]) %*% Sigma.GEE.inv[[i]] %*% M.list[[i]]
    })
    FI <- Reduce("+", method.list)
    var.interact <- solve(FI)[length(reg.coef), length(reg.coef)]
  }
  else if(method == "GEE-KC"){
    A.list <- lapply(L.list, function(i){
      i - (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2) * i %*% i
    })
    Vtilde.list <- lapply(cluster.index, function(i){
      L.list[[i]] %*% D.list[[i]] %*% L.list[[i]] + A.list[[i]]
    })
    rm(D.list)
    constant <- 16 * sqrt(3) / (15 * pi)
    cluster.constant <- 1/sqrt(1 + constant^2 * (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2))
    Ltilde.list <- lapply(M.list, function(i){
      u <- inv.logit(cluster.constant * as.numeric(i %*% reg.coef))
      Matrix(diag(u * (1 - u)), sparse = TRUE)
    })
    GEE.list <- lapply(cluster.index, function(i){
      cluster.constant^2 * t(M.list[[i]]) %*% Ltilde.list[[i]] %*% solve(Vtilde.list[[i]]) %*% Ltilde.list[[i]] %*% M.list[[i]]
    })
    FI.GEE <- Reduce("+", GEE.list)
    method.list <- lapply(cluster.index, function(i){
      grad <- cluster.constant * Ltilde.list[[i]] %*% M.list[[i]]
      H <- grad %*% solve(FI.GEE) %*% t(grad) %*% solve(Vtilde.list[[i]])
      Fa <- sqrtm(solve(I.matrix(T.period * m) - H))
      B <- cluster.constant * Ltilde.list[[i]] %*% solve(Vtilde.list[[i]]) %*% Fa
      t(M.list[[i]]) %*% B %*% Vtilde.list[[i]] %*% t(B) %*% M.list[[i]]
    })
    FI.GEE.RBC <- Reduce("+", method.list)
    V.ultra <- solve(FI.GEE) %*% FI.GEE.RBC %*% solve(FI.GEE)
    var.interact <- V.ultra[length(reg.coef), length(reg.coef)]
  }
  else if(method == "GEE-MD"){
    A.list <- lapply(L.list, function(i){
      i - (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2) * i %*% i
    })
    Vtilde.list <- lapply(cluster.index, function(i){
      L.list[[i]] %*% D.list[[i]] %*% L.list[[i]] + A.list[[i]]
    })
    rm(D.list)
    constant <- 16 * sqrt(3) / (15 * pi)
    cluster.constant <- 1/sqrt(1 + constant^2 * (sigma.alpha^2 + sigma.nu^2 + sigma.psi^2))
    Ltilde.list <- lapply(M.list, function(i){
      u <- inv.logit(cluster.constant * as.numeric(i %*% reg.coef))
      diag(u * (1 - u))
    })
    GEE.list <- lapply(cluster.index, function(i){
      cluster.constant^2 * t(M.list[[i]]) %*% Ltilde.list[[i]] %*% solve(Vtilde.list[[i]]) %*% Ltilde.list[[i]] %*% M.list[[i]]
    })
    FI.GEE <- Reduce("+", GEE.list)
    method.list <- lapply(cluster.index, function(i){
      grad <- cluster.constant * Ltilde.list[[i]] %*% M.list[[i]]
      H <- grad %*% solve(FI.GEE) %*% t(grad) %*% solve(Vtilde.list[[i]])
      Fa <- solve(I.matrix(T.period * m) - H)
      B <- cluster.constant * Ltilde.list[[i]] %*% solve(Vtilde.list[[i]]) %*% Fa
      t(M.list[[i]]) %*% B %*% Vtilde.list[[i]] %*% t(B) %*% M.list[[i]]
    })
    FI.GEE.RBC <- Reduce("+", method.list)
    V.ultra <- solve(FI.GEE) %*% FI.GEE.RBC %*% solve(FI.GEE)
    var.interact <- V.ultra[length(reg.coef), length(reg.coef)]
  }
  else stop("The input method is unavailable.")
  calc.power <- pnorm((reg.coef[length(reg.coef)] 
                       / sqrt(var.interact) - qnorm(1 - (sig.level / 2)))) + pnorm((-reg.coef[length(reg.coef)] 
                                                                                    / sqrt(var.interact) - qnorm(1 - (sig.level / 2))))
  return(calc.power)
}

cluster.size <- function(size.range, method,
                         reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                         sigmas = c(1, 0.5, 0.75), 
                         prevalences = 0.5, 
                         I = 8, 
                         sig.level = 0.05,
                         req.power = 0.8){
  calc.powers <- sapply(size.range, function(m){
    calc.power(m, 
               method = method,
               reg.coef = reg.coef, 
               sigmas = sigmas, 
               prevalences = prevalences, 
               I = I,
               sig.level = sig.level)
  })
  if(max(calc.powers) < req.power){
    max.power <- max(calc.powers)
    warning.message <- paste("The largest power", max.power, "is smaller than the required power! please increase the cluster size.", sep = " ")
    warning(warning.message)
  }
  else return(list("size"=size.range[min(which(calc.powers >= req.power))], "power"=min(calc.powers[calc.powers >= req.power])))
}

number.cluster <- function(m=46, method,
                         reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                         sigmas = c(1, 0.5, 0.75), 
                         prevalences = 0.5, 
                         number.range,
                         sig.level = 0.05,
                         req.power = 0.8){
  J = length(reg.coef) - 4
  if(!all(number.range %% J == 0)) stop("Not all numbers of clusters are multiple of the sequence length!")
  else{
    calc.powers <- sapply(number.range, function(I){
      calc.power(m, 
                 method = method,
                 reg.coef = reg.coef, 
                 sigmas = sigmas, 
                 prevalences = prevalences, 
                 I = I, 
                 sig.level = sig.level)
    })
    if(max(calc.powers) < req.power){
      max.power <- max(calc.powers)
      warning.message <- paste("The largest power", max.power, "is smaller than the required power! please increase the cluster size.", sep = " ")
      warning(warning.message)
    }
    else return(list("Number"=number.range[min(which(calc.powers >= req.power))], "power"=min(calc.powers[calc.powers >= req.power])))
    
  }
  }
