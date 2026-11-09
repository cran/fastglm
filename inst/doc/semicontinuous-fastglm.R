## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", message=FALSE)

## -----------------------------------------------------------------------------
library(fastglm)
set.seed(619)
d <- data.frame(age=rnorm(800,50,10), severity=rnorm(800))
probability <- plogis(-.4+.02*(d$age-50)+.5*d$severity)
positive_mean <- exp(5+.015*(d$age-50)+.4*d$severity)
d$cost <- rbinom(800,1,probability)*rgamma(800,shape=2,scale=positive_mean/2)
fit <- fastglm_twopart(cost~age+severity|age+severity,data=d)
summary(fit)

## -----------------------------------------------------------------------------
fastglm_wald_test(fit,terms="severity",component="positive")
fastglm_wald_test(fit,terms="severity",component="binary")
fastglm_wald_test(fit,terms="severity")

## -----------------------------------------------------------------------------
new <- data.frame(age=c(40,60),severity=c(0,0))
cbind(probability_positive=1-predict(fit,new,type="zero"),
      positive_mean=predict(fit,new,type="count"),
      marginal_mean=predict(fit,new,type="response"))
predict(fit,new,type="response",se.fit=TRUE)

## -----------------------------------------------------------------------------
replicates <- simulate(fit,nsim=2,seed=620)
vapply(replicates,mean,numeric(1))
mean(d$cost)

## -----------------------------------------------------------------------------
change_severity <- function(design) {
    design$count[,"severity"] <- design$count[,"severity"]+1
    design$zero[,"severity"] <- design$zero[,"severity"]+1
    design
}
fastglm_average(fit,contrast=change_severity)

## -----------------------------------------------------------------------------
cluster <- rep(seq_len(160),each=5)
V <- fastglm_vcov(fit,type="HC1",cluster=cluster,full=TRUE)
fastglm_average(fit,contrast=change_severity,vcov.=V)
fastglm_wald_test(fit,terms="severity",vcov.=V)

## -----------------------------------------------------------------------------
read_cost <- function(k) d[((k-1)*200+1):(k*200),]
cost_source <- fastglm_source(read_cost,4,
    formula=cost~age+severity|age+severity,prototype=d)
stream_fit <- fastglm_twopart(cost_source)
max(abs(unname(coef(stream_fit))-unname(coef(fit))))
fastglm_average(stream_fit)

