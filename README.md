# USLM (Small area estimation via unmatched sampling and linking models)
This package implements small area estimation via unmatched sampling and linking models, as proposed by the following papers.

Sugasawa, S., Kubokawa, T. and Rao, J. N. K. (2018).  Small area estimation via unmatched sampling and linking models.  *TEST* 27, 407-427. https://link.springer.com/article/10.1007/s11749-017-0551-5

Functions are implemented in `USLM-function.R` available in the repository.
```{r}
source("USLM-function.R")
```
 
We use the following synthetic data using logistic link function.
```{r}
set.seed(1)
m=50
x1=runif(m,1,3); x2=runif(m,-1,1)
X=cbind(1,x1,x2)
Di=runif(m,0,0.05)

beta=c(-2,1,0.5); A=0.7
mu=as.vector(X%*%beta)

logit=function(x){ exp(x)/(1+exp(x)) }
th=logit(mu+sqrt(A)*rnorm(m))
y=th+sqrt(Di)*rnorm(m)
y
```

Apply the proposed method with logistic link.

INPUT (`USLM`):
- `y`: A vector of responce variables.
- `X`: A matrix of covariates whose column vectors have the same dimension as y.
- `Di`: A vector of sampling variances.
- `link`: A link function chosen from "exp" (exponetial link) or "logit" (logit link). The default is "exp".
- `n`: The number of GH quadrature. The default value is 100.
- `dd`: A small positive value. The EM algorithm is terminated if the relative difference is smaller than this value. Default value is 0.0001
- `Init`: A vector of initial parameter values. The default is NA in which the initial vlaues are set to A=1 and 0 for regression coeffieicnts.

OUTPUT: 
- `Est`: A vector of estimates of model parameters.
- `EB`: A vector of empirical Bayes estimates.
- `AEB`: A vector of approximate empirical Bayes estimates.

```{r}
fit1=USLM(y=y,X=X,Di=Di,link="logit")
fit1$Est
fit1$EB
fit1$AEB
```

Next, calcualte mean squared error estimates.

INPUT (`mseUSLM`): Additionally to the input in the function `USLM`
- `mc`: The number of Monte Carlo samples used for computing the leading term of MSE. The default value is 1000
- `B`: The number of bootstrap samples for MSE

OUTPUT: A vector of MSE estimates.
```{r}
mseUSLM(y=y,X=X,Di=Di,link="logit",B=100)
```

Finally, fit the USLM with P-spline.

INPUT (`spUSLM`): Additionally to the input in the function `USLM`
- `Z`: A matrix for P-spline.
- `mc`: The initial number of Monte Carlo samples used in E-step. The default value is 3000.
- `burn`: The number of burn-in period in E-step. The default value is 2000.
- `band`: The value of a stap-size in MH update in E-step. The default value is 0.3.

OUTPUT: 
- `Est`: A vector of estimates of model parameters.
- `EB`: A vector of empirical Bayes estimates.
- `AEB`: A vector of approximate empirical Bayes estimates.
- `Gamma`: A vector of estimates of gamma (coefficients of spline terms).
```{r}
K=10   # number of knots
Kappa=kmeans(X[,-1],K)$centers
Z=matrix(NA,m,K)
for(j in 1:K){
  for(i in 1:m){ 
    r=sqrt(sum((X[i,-1]-Kappa[j,])^2))
    Z[i,j]=exp(-0.5*r^2/(0.5)^2)
  }
}

fit2=spUSLM(y=y,X=X,Z,Di,link="logit")
fit2$Est
fit2$EB
plot(fit1$EB,fit2$EB)
```
