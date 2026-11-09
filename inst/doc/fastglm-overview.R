## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", message=FALSE)

## -----------------------------------------------------------------------------
library(fastglm)
set.seed(41)
d <- data.frame(x=rnorm(500), group=factor(rep(letters[1:4],125)),
                exposure=runif(500,.5,2))
d$y <- rpois(500, d$exposure*exp(.2+.4*d$x))
fit <- fastglm(y~x+group+offset(log(exposure)), data=d, family=poisson(),
               method="auto")
summary(fit)
predict(fit, newdata=d[1:4,], type="response", se.fit=TRUE)

## -----------------------------------------------------------------------------
fastglm_diagnostics(fit)

## -----------------------------------------------------------------------------
vcov(fit)
confint(fit, parm="x")
vcov(fit, type="HC1")

## -----------------------------------------------------------------------------
cluster <- rep(seq_len(100), each=5)
robust <- fastglm_vcov(fit, type="HC1", cluster=cluster)
fastglm_average(fit, vcov.=robust)

## -----------------------------------------------------------------------------
fastglm_wald_test(fit,terms="group")
fastglm_wald_test(fit,terms="group",vcov.=robust)

## -----------------------------------------------------------------------------
group_comparison <- fastglm_wald_test(fit,contrast=c(groupc=1,groupb=-1))
group_comparison
as.data.frame(group_comparison)

## -----------------------------------------------------------------------------
X <- model.matrix(~x+group, d)
Y <- cbind(outcome1=rnorm(nrow(d))+.4*d$x,
           outcome2=rnorm(nrow(d))-.2*d$x)
design <- fastglm_prepare(X)
models <- fastglm_batch(design, Y)
vapply(models, function(z) coef(z)["x"], numeric(1))

## -----------------------------------------------------------------------------
set.seed(148)
mu1 <- exp(.2+.3*d$x)
mu2 <- exp(.5-.2*d$x)
Y_positive <- cbind(outcome1=rgamma(nrow(d),shape=3,scale=mu1/3),
                    outcome2=rgamma(nrow(d),shape=3,scale=mu2/3))
design_firth <- fastglm_prepare(X,firth=TRUE)
models_firth <- fastglm_batch(design_firth,Y_positive,family=Gamma("log"))
vapply(models_firth,function(z) coef(z)["x"],numeric(1))
vapply(models_firth,function(z) z$dispersion,numeric(1))

## -----------------------------------------------------------------------------
compact <- fastglm(X, d$y, poisson(), method="auto",
                   inference="none", retain.data=FALSE)
coef(compact)

