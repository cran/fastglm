## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", message=FALSE)

## -----------------------------------------------------------------------------
library(fastglm)
set.seed(544)
d <- data.frame(x=rnorm(500), exposure=runif(500,.8,1.2))
d$y <- rnbinom(500,mu=d$exposure*exp(.3+.5*d$x),size=1.5)

## -----------------------------------------------------------------------------
X <- model.matrix(~x,d)
nb <- fastglm_nb(X,d$y,offset=log(d$exposure))
c(theta=nb$theta,SE.theta=nb$SE.theta)
summary(nb)

## -----------------------------------------------------------------------------
fastglm_confint_theta(nb)

## -----------------------------------------------------------------------------
positive <- subset(d,y>0)
truncated <- fastglm_truncated(y~x+offset(log(exposure)),data=positive,
                               dist="negbin")
head(predict(truncated,se.fit=TRUE)$fit)

## -----------------------------------------------------------------------------
hurdle <- fastglm_hurdle(y~x+offset(log(exposure))|x,data=d,dist="negbin")
summary(hurdle)
predict(hurdle,newdata=d[1:4,],type="response",se.fit=TRUE)
predict(hurdle,newdata=d[1:4,],type="zero")
predict(hurdle,newdata=d[1:4,],type="count")

## -----------------------------------------------------------------------------
set.seed(545)
d$inflated <- ifelse(rbinom(nrow(d),1,plogis(-.8+.4*d$x)),0,d$y)
zi <- fastglm_zi(inflated~x+offset(log(exposure))|x,data=d,
                 dist="poisson",em.maxit=300)
predict(zi,newdata=d[1:4,],type="prob",at=0:4)
predict(zi,newdata=d[1:4,],type="response",se.fit=TRUE)

## -----------------------------------------------------------------------------
fastglm_wald_test(hurdle,terms="x")
fastglm_wald_test(hurdle,terms="x",component="count")
fastglm_wald_test(zi,terms="x",component="zero")
fastglm_wald_test(truncated,terms="x")

## -----------------------------------------------------------------------------
regression <- names(coef(nb))
joint_covariance <- vcov(nb,full=TRUE)[regression,regression,drop=FALSE]
fastglm_wald_test(nb,terms="x",vcov.=joint_covariance)

## -----------------------------------------------------------------------------
replicates <- simulate(hurdle,nsim=3,seed=546)
c(observed=mean(d$y==0),vapply(replicates,function(z) mean(z==0),numeric(1)))

## -----------------------------------------------------------------------------
Xs <- as(model.matrix(~x,d),"dgCMatrix")
nb_sparse <- fastglm_nb(Xs,d$y,offset=log(d$exposure))
hurdle_sparse <- fastglm_hurdle(Xs,response=d$y,z=Xs,
                                offset=log(d$exposure),dist="negbin")
max(abs(unname(coef(nb_sparse))-unname(coef(nb))))

## -----------------------------------------------------------------------------
read_block <- function(k) d[((k-1)*100+1):(k*100),]
count_source <- fastglm_source(read_block,5,
    formula=y~x+offset(log(exposure))|x,prototype=d)
hurdle_stream <- fastglm_hurdle(count_source,dist="negbin")
max(abs(unname(coef(hurdle_stream))-unname(coef(hurdle))))
summary(hurdle_stream)

## -----------------------------------------------------------------------------
cluster <- rep(seq_len(100),each=5)
V <- fastglm_vcov(hurdle,type="HC1",cluster=cluster,full=TRUE)
fastglm_average(hurdle,vcov.=V)

## -----------------------------------------------------------------------------
boot <- fastglm_bootstrap(hurdle,R=49,seed=547,type="cluster",cluster=cluster,
    statistic=function(model) mean(predict(model,newdata=d,type="response")))
boot

