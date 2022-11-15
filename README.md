# Power calculation for SW-CRT design based on GLMM with binary outcomes

The `PowerSWCRT.R` file contains three R functions: `calc.power` (for power calculation), `cluster.size` (cluster size calculation for given power), and `number.cluster` (determine the required number of clusters given the power and cluster size).

## Power calculation

In order to calculate the power, one needs to specify the following parameters:
- cluser size `m` (number of participants in each cluster per sequential period);
- `method` (one of the following four methods: "AML", "GEE", "GEE-KC", "GEE-MD");
- all fixed effects `reg.coef` (a vector of all fixed effects, including the baseline effect, sequential period fixed effects, overall treatment effects, single covariate fixed effect, and heterogenity of treatment effect);
- all standard deviations for random effects `sigmas` (a vector of length 3 including the cluster random effect, sequential period random effect nested in clusters, and individual repeated measurement random effect nested in clusters);
- proportion of the single binary covariate within each cluster `prevalence`;
- number of clusters `I`;
- significance level (type I error) `sig.level`.

For example, for an SW-CRT design with 8 clusters, 4 sequential periods, half of the individuals with single binary covariate equal to 1 within each cluster, the power can be calculated as follows: 
```r
source("PowerSWCRT.R")
(power <- calc.power(m = 20, method = "GEE-MD",
                       reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                       sigmas = c(1, 0.5, 0.75), 
                       prevalences = 0.5, 
                       I = 8, 
                       sig.level = 0.05))
```
using the method "GEE-MD" under the 0.05 significance level.
## Sample size calculation

The sample size calculation is based on the power calculated by `calc.power`. Hence in addition to the parameter specification of `calc.power`, one has to provide:
- the required power `req.power` and;
- a range of candidate cluster sizes `size.range` instead of a cluser size value `m`.

For example, given the required power 0.8, for the abovementioned specification the sample size for each cluster can be determined from candidate sizes 40, 42, 44,..., 60 as follows:
```r
cluster.size(size.range = seq(40,60,2), method = "GEE-MD",
                         reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                         sigmas = c(1, 0.5, 0.75), 
                         prevalences = 0.5, 
                         I = 8, 
                         sig.level = 0.05,
                         req.power = 0.8)
```
The result looks like:
```
$size
[1] 46

$power
[1] 0.8006238
```
## Number of clusters
In some cases, one may consider varying the number of clusters instead of cluster sizes. For these cases, `number.cluster` can help determine the desirable number of clusters for the SW-CRT design. The idea is similar to `cluster.size`, i.e. finding out the smallest number of clusters such that the required power can be achieved. Hence in addition to the parameter specification of `calc.power`, one has to provide:
- the required power `req.power` and;
- a range of candidate number of clusters `number.range` which must be multiple of the number of sequential periods, instead of a single value of number of clusters `I`.

For example, given the required power 0.8 and fixed cluster size 46, for the abovementioned specification the number of clusters can be determined from candidate sizes 4, 8, 12,..., 40 as follows:
```r
number.cluster(m = 46, method = "GEE-MD",
                           reg.coef = c(log(0.15/0.85), 0.1, 0.2, 0.3, 0.4, log(1.68), log(1.5), log(2)), 
                           sigmas = c(1, 0.5, 0.75), 
                           prevalences = 0.5, 
                           number.range = seq(4,40,4),
                           sig.level = 0.05,
                           req.power = 0.8)
```
The result looks like:
```
$Number
[1] 8

$power
[1] 0.8006238
```
