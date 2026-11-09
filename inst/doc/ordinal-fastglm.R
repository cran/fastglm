## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,comment="#>",message=FALSE,
                      eval=requireNamespace("MASS",quietly=TRUE))
library(fastglm)

## -----------------------------------------------------------------------------
housing <- MASS::housing
levels(housing$Sat)
c(rows=nrow(housing),residents=sum(housing$Freq))

fit_ml <- fastglm_ordinal(Sat~Infl+Type+Cont,data=housing,
                          weights=Freq,tol=1e-10)
fit_br <- update(fit_ml,firth=TRUE)
profiles <- unique(housing[,c("Infl","Type","Cont")])[1:4,]
summary(fit_br)

## -----------------------------------------------------------------------------
cbind(maximum_likelihood=coef(fit_ml,component="location"),
      mean_bias_reduction=coef(fit_br,component="location"))
exp(coef(fit_br,component="location"))
fastglm_diagnostics(fit_br)

## -----------------------------------------------------------------------------
fastglm_wald_test(fit_br,terms="Type")
fastglm_wald_test(fit_br,contrast=c(InflHigh=1,InflMedium=-1))

## -----------------------------------------------------------------------------
fit_probit <- update(fit_br,link="probit")
fit_cloglog <- update(fit_br,link="cloglog")
link_models <- list(logit=fit_br,probit=fit_probit,cloglog=fit_cloglog)
link_averages <- do.call(rbind,lapply(names(link_models),function(link) {
    data.frame(link=link,fastglm_average(link_models[[link]],weights=housing$Freq))
}))
rownames(link_averages) <- NULL
link_averages

## -----------------------------------------------------------------------------
fit_partial <- fastglm_ordinal(
    Sat~Infl+Type+Cont,data=housing,weights=Freq,
    nonparallel=~Cont,tol=1e-10)
summary(fit_partial)
exp(fit_partial$gamma)

## -----------------------------------------------------------------------------
fastglm_wald_test(fit_partial,terms="Cont")
predict(fit_partial,newdata=profiles,se.fit=TRUE)
fastglm_average(fit_partial,weights=housing$Freq)

## -----------------------------------------------------------------------------
fastglm_wald_test(fit_partial,contrast=c(
    "nonparallel:ContHigh:Low|Medium"=1,
    "nonparallel:ContHigh:Medium|High"=-1))

## -----------------------------------------------------------------------------
fastglm_ordinal_assess(fit_partial)

## -----------------------------------------------------------------------------
wine <- data.frame(contact=c(0,1,0,1),temperature=c(0,0,1,1))
wine$rating <- I(rbind(c(4,9,5,0,0),c(1,7,8,2,0),
                       c(0,5,8,3,2),c(0,1,5,7,5)))
colnames(wine$rating) <- paste0("rating",1:5)
wine_br <- fastglm_ordinal(
    rating~contact+temperature,data=wine,
    nonparallel=~temperature,firth=TRUE,tol=1e-9)
summary(wine_br)
fastglm_diagnostics(wine_br)

## -----------------------------------------------------------------------------
pred <- predict(fit_br,newdata=profiles,se.fit=TRUE)
pred
predict(fit_br,newdata=profiles,type="cumulative")

## -----------------------------------------------------------------------------
fastglm_average(fit_br,weights=housing$Freq)

low_contact <- housing
low_contact$Cont <- factor(rep("Low",nrow(housing)),
                           levels=levels(housing$Cont))
contact_difference <- fastglm_average(
    fit_br,newdata=low_contact,weights=housing$Freq,
    contrast=function(x) {
        x[,"ContHigh"] <- 1
        x
    })
contact_difference

## ----ordinal-marginal-effects-------------------------------------------------
set.seed(110303)
severity <- data.frame(exposure=runif(500,0,3),
                       group=factor(rep(c("A","B"),250)))
location <- with(severity,.7*exposure-.15*exposure^2+
                  .35*exposure*(group=="B"))
severity$outcome <- ordered(cut(location+rlogis(nrow(severity)),
    c(-Inf,.3,1.6,Inf),labels=c("low","moderate","high")))
severity_fit <- fastglm_ordinal(
    outcome~exposure+I(exposure^2)+exposure:group,
    data=severity,firth=TRUE)
exposure_effect <- fastglm_marginal_effects(
    severity_fit,"exposure",newdata=severity)
exposure_effect

## ----ordinal-partial-marginal-effects-----------------------------------------
severity_partial <- update(severity_fit,nonparallel=~exposure)
fastglm_marginal_effects(severity_partial,"exposure",newdata=severity)

## -----------------------------------------------------------------------------
confint(fit_ml,parm="ContHigh",method="profile")
fastglm_profile_test(fit_ml,
    parm=c("TypeApartment","TypeAtrium","TypeTerrace"))

## -----------------------------------------------------------------------------
type_score <- fastglm_score_test(fit_br,terms="Type")
type_score
as.data.frame(type_score)

## -----------------------------------------------------------------------------
fastglm_score_test(fit_br,terms="Type",null.estimation="meanBR")
fastglm_score_test(fit_br,
    contrast=c(InflHigh=1,InflMedium=-1),null.estimation="meanBR")

## ----ordinal-partial-term-test------------------------------------------------
fastglm_po_test(fit_partial,terms="Cont")
fastglm_po_test(fit_partial,terms="Cont",method="lr")

## ----ordinal-partial-term-tests-----------------------------------------------
fit_varying <- update(fit_partial,nonparallel=~Type+Cont)
fastglm_po_test(fit_varying,terms=c("Type","Cont"),by.term=TRUE)

## ----ordinal-partial-constraint-----------------------------------------------
contact_equality <- c(
    "nonparallel:ContHigh:Low|Medium"=1,
    "nonparallel:ContHigh:Medium|High"=-1)
fastglm_profile_test(fit_partial,contrast=contact_equality)
fastglm_score_test(fit_partial,contrast=contact_equality)
confint(fit_partial,parm="nonparallel:ContHigh:Low|Medium",
        method="profile")

## ----ordinal-partial-adjusted-score-------------------------------------------
fit_partial_br <- update(fit_partial,firth=TRUE)
fastglm_score_test(fit_partial_br,contrast=contact_equality,
                   null.estimation="meanBR")

## -----------------------------------------------------------------------------
V_robust <- fastglm_vcov(fit_br,type="HC0")
cbind(model_based=sqrt(diag(vcov(fit_br))),
      sandwich=sqrt(diag(V_robust)))
fastglm_average(fit_br,weights=housing$Freq,vcov.=V_robust)

## ----ordinal-partial-robust---------------------------------------------------
V_partial <- fastglm_vcov(fit_partial_br,type="HC0")
fastglm_wald_test(fit_partial_br,contrast=contact_equality,
                  vcov.=V_partial)
fastglm_average(fit_partial_br,weights=housing$Freq,vcov.=V_partial)

## -----------------------------------------------------------------------------
boot_location <- fastglm_bootstrap(
    fit_br,R=199,type="parametric",seed=61103,
    statistic=function(fit) coef(fit,component="location"))
boot_location

## -----------------------------------------------------------------------------
probability_residual <- residuals(fit_br)
head(cbind(probability_scale=probability_residual,
           residuals(fit_br,type="cumulative")))
by(seq_len(nrow(housing)),housing$Infl,function(i)
    weighted.mean(probability_residual[i],housing$Freq[i],na.rm=TRUE))

## -----------------------------------------------------------------------------
if(requireNamespace("broom",quietly=TRUE)) {
    broom::tidy(fit_br,component="location",conf.int=TRUE,exponentiate=TRUE)
    broom::glance(fit_br)
    broom::augment(fit_br,newdata=profiles,se_fit=TRUE)
}

## -----------------------------------------------------------------------------
set.seed(61105)
training_frequency <- rbinom(nrow(housing),housing$Freq,.7)
evaluation_frequency <- housing$Freq-training_frequency
training_fit <- fastglm_ordinal(
    Sat~Infl+Type+Cont,data=housing,weights=training_frequency,
    firth=TRUE,tol=1e-10)
assessment <- fastglm_ordinal_assess(
    training_fit,newdata=housing,y=housing$Sat,
    weights=evaluation_frequency,bins=5)
assessment$scores

## ----fig.width=7, fig.height=3.5----------------------------------------------
plot(assessment,type="calibration")

## ----fig.width=8, fig.height=8------------------------------------------------
assumption_check <- fastglm_ordinal_assess(fit_br,proportional.odds=TRUE)
head(assumption_check$proportional.odds[
    ,c("category","term","estimate","common.estimate","converged")])
plot(assumption_check,type="proportional.odds")

## -----------------------------------------------------------------------------
po_test <- fastglm_po_test(fit_br)
po_test

## -----------------------------------------------------------------------------
type_test <- fastglm_po_test(fit_br,terms="Type")
type_test
as.data.frame(type_test)

## -----------------------------------------------------------------------------
po_terms <- fastglm_po_test(fit_br,by.term=TRUE)
po_terms
as.data.frame(po_terms)

## -----------------------------------------------------------------------------
exposure <- matrix(c(0,1),ncol=1,dimnames=list(NULL,"exposed"))
sparse_counts <- rbind(c(12,18,20),c(0,0,4))
colnames(sparse_counts) <- c("low","middle","high")
sparse_fit <- fastglm_ordinal(exposure,sparse_counts,firth=TRUE)
fastglm_po_test(sparse_fit)

## -----------------------------------------------------------------------------
fastglm_po_test(sparse_fit,null.estimation="meanBR")

## -----------------------------------------------------------------------------
simulated <- simulate(fit_br,nsim=2,seed=61104)
head(simulated[[1]])
all(rowSums(simulated[[1]])==housing$Freq)

## -----------------------------------------------------------------------------
prepared <- fastglm_ordinal_prepare(
    ~Infl+Type+Cont,data=housing,
    levels=levels(housing$Sat),firth=TRUE)
prepared_original <- fastglm_ordinal(
    prepared,housing$Sat,weights=housing$Freq,tol=1e-10)
all.equal(coef(prepared_original),coef(fit_br),tolerance=1e-8)

repeated <- fastglm_ordinal_batch(
    prepared,simulated,weights=rep(1,nrow(housing)),tol=1e-10)
vapply(repeated,function(fit) fit$converged,logical(1))
predict(repeated[[1]],newdata=profiles)

## -----------------------------------------------------------------------------
restored <- unserialize(serialize(prepared,NULL))
restored_fit <- fastglm_ordinal(
    restored,simulated[[1]],weights=rep(1,nrow(housing)),tol=1e-10)
all.equal(coef(restored_fit),coef(repeated[[1]]),tolerance=1e-10)

## ----ordinal-partial-prepared-------------------------------------------------
prepared_partial <- fastglm_ordinal_prepare(
    ~Infl+Type+Cont,data=housing,levels=levels(housing$Sat),
    nonparallel=~Cont,firth=TRUE)
partial_repeated <- fastglm_ordinal_batch(
    prepared_partial,simulated,weights=rep(1,nrow(housing)),tol=1e-10)
summary(partial_repeated[[1]])

## -----------------------------------------------------------------------------
housing_source <- fastglm_source(
    function(k) housing[((k-1L)*18L+1L):(k*18L),,drop=FALSE],
    n_chunks=4L,formula=Sat~Infl+Type+Cont,
    prototype=housing,weights="Freq")
fit_stream <- fastglm_ordinal_streaming(housing_source,firth=TRUE)
fastglm_diagnostics(fit_stream)

## ----ordinal-partial-streaming------------------------------------------------
partial_stream <- fastglm_ordinal_streaming(
    housing_source,nonparallel=~Cont,firth=TRUE,tol=1e-10)
summary(partial_stream)
predict(partial_stream,newdata=profiles,se.fit=TRUE)

## ----ordinal-source-robust----------------------------------------------------
V_source <- fastglm_vcov(partial_stream,type="HC0")
fastglm_average(partial_stream,newdata=housing,weights=housing$Freq,
                vcov.=V_source)

## ----ordinal-source-assessment------------------------------------------------
source_assessment <- fastglm_ordinal_assess(partial_stream,bins=5)
source_assessment
plot(source_assessment,type="calibration")

## ----ordinal-source-common-slope-tests----------------------------------------
fastglm_po_test(partial_stream,terms="Cont")
fastglm_po_test(partial_stream,terms="Cont",method="lr")

## ----ordinal-boundary-row-diagnostics-----------------------------------------
fastglm_diagnostics(wine_br,boundary=TRUE,max.rows=3)

## ----ordinal-pattern-boundary-------------------------------------------------
pattern_x <- matrix(c(0,1),ncol=1,
                    dimnames=list(NULL,"group"))
pattern_counts <- rbind(c(low=6,middle=0,high=4),
                        c(low=3,middle=5,high=2))
pattern_message <- tryCatch(
    fastglm_ordinal(pattern_x,pattern_counts,nonparallel="group",firth=TRUE),
    error=conditionMessage)
pattern_message

## -----------------------------------------------------------------------------
no_predictors <- matrix(numeric(),nrow=1,ncol=0)
endpoint_counts <- matrix(c(0,8,2),nrow=1,
                           dimnames=list(NULL,c("low","middle","high")))
endpoint_fit <- fastglm_ordinal(no_predictors,endpoint_counts,firth=TRUE)
coef(endpoint_fit,component="threshold")
predict(endpoint_fit)

## -----------------------------------------------------------------------------
endpoint_boot <- fastglm_bootstrap(
    endpoint_fit,R=199,type="parametric",seed=61105)
endpoint_boot
data.frame(
  `Category pattern`=c("Empty endpoint","Empty interior category"),
  `Proportion of refits`=c(mean(lengths(
    endpoint_boot$diagnostics$empty.end.categories)>0),
    mean(endpoint_boot$diagnostics$boundary)),check.names=FALSE)

## -----------------------------------------------------------------------------
interior_counts <- matrix(c(4,0,7,3),nrow=1,
                           dimnames=list(NULL,c("low","mild","moderate","high")))
interior_fit <- fastglm_ordinal(no_predictors,interior_counts,firth=TRUE)
coef(interior_fit,component="threshold")
predict(interior_fit,se.fit=TRUE)
diag(vcov(interior_fit))

