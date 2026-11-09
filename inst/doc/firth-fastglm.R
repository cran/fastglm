## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  echo = TRUE,
  eval = all(vapply(c("MASS", "brglm2", "logistf", "microbenchmark", "statmod"),
                    requireNamespace, logical(1), quietly = TRUE))
)

## ----eval = FALSE-------------------------------------------------------------
# fastglm(x, y, family = binomial(), firth = TRUE)

## -----------------------------------------------------------------------------
library(fastglm)

## -----------------------------------------------------------------------------
data(sex2, package = "logistf")
X <- model.matrix(~ age + oc + vic + vicl + vis + dia, data = sex2)
y <- sex2$case

fit_f <- fastglm(X, y, family = binomial(), firth = TRUE)
fit_l <- logistf::logistf(case ~ age + oc + vic + vicl + vis + dia,
                          data = sex2, pl = FALSE, plconf = NULL,
                          control = logistf::logistf.control(
                              lconv = 1e-12, gconv = 1e-12, xconv = 1e-12,
                              maxit = 200L))

cbind(fastglm = unname(coef(fit_f)), logistf = unname(coef(fit_l)))
max(abs(unname(coef(fit_f)) - unname(coef(fit_l))))

## -----------------------------------------------------------------------------
mb <- microbenchmark::microbenchmark(
    fastglm = fastglm(X, y, family = binomial(), firth = TRUE),
    logistf = logistf::logistf(case ~ age + oc + vic + vicl + vis + dia,
                               data = sex2, pl = FALSE, plconf = NULL),
    times = 3L
)
knitr::kable(summary(mb, unit = "ms")[, c("expr", "min", "median", "max")], digits = 2,
             col.names = c("Method", "Minimum (ms)", "Median (ms)", "Maximum (ms)"), row.names = FALSE)

## -----------------------------------------------------------------------------
library(brglm2)
set.seed(123)
n <- 500
x1 <- rnorm(n);  x2 <- rnorm(n)
X <- cbind(1, x1, x2)
eta_true <- 0.5 + 0.3 * x1 - 0.2 * x2

## -----------------------------------------------------------------------------
# Largest coefficient difference from brglm2, by link.
difference_table <- function(differences) {
    knitr::kable(data.frame(Link = names(differences),
                            Difference = formatC(differences, format = "e", digits = 1)),
                 col.names = c("Link", "Largest coefficient difference"), row.names = FALSE)
}
y_bin <- rbinom(n, 1, plogis(eta_true))
df <- data.frame(y = y_bin, x1 = x1, x2 = x2)

differences <- sapply(c("logit", "probit", "cloglog"), function(lnk) {
    fam <- binomial(lnk)
    f <- fastglm(X, y_bin, family = fam, firth = TRUE, tol = 1e-10)
    b <- glm(y ~ x1 + x2, data = df, family = fam, method = "brglmFit",
             type = "AS_mean", control = list(epsilon = 1e-10, maxit = 200))
    max(abs(unname(coef(f)) - unname(coef(b))))
})
difference_table(differences)

## -----------------------------------------------------------------------------
y_pois <- rpois(n, exp(eta_true))

differences <- sapply(c("log", "sqrt"), function(lnk) {
    fam <- poisson(lnk)
    f <- fastglm(X, y_pois, family = fam, firth = TRUE, tol = 1e-10)
    b <- glm(y ~ x1 + x2, data = data.frame(y = y_pois, x1 = x1, x2 = x2),
             family = fam, method = "brglmFit", type = "AS_mean",
             control = list(epsilon = 1e-10, maxit = 200))
    max(abs(unname(coef(f)) - unname(coef(b))))
})
difference_table(differences)

## -----------------------------------------------------------------------------
y_gam <- rgamma(n, shape = 5, rate = 5 / exp(eta_true))

differences <- sapply(c("log", "inverse"), function(lnk) {
    fam <- Gamma(lnk)
    f <- fastglm(X, y_gam, family = fam, firth = TRUE, tol = 1e-10)
    b <- glm(y ~ x1 + x2, data = data.frame(y = y_gam, x1 = x1, x2 = x2),
             family = fam, method = "brglmFit", type = "AS_mean",
             control = list(epsilon = 1e-10, maxit = 200))
    max(abs(unname(coef(f)) - unname(coef(b))))
})
difference_table(differences)

## -----------------------------------------------------------------------------
y_gauss <- eta_true + rnorm(n, 0, 0.5)

differences <- sapply(c("identity", "log"), function(lnk) {
    fam <- gaussian(lnk)
    y0 <- if (lnk == "log") pmax(y_gauss, 0.01) else y_gauss
    f <- fastglm(X, y0, family = fam, firth = TRUE, tol = 1e-10)
    b <- glm(y ~ x1 + x2, data = data.frame(y = y0, x1 = x1, x2 = x2),
             family = fam, method = "brglmFit", type = "AS_mean",
             control = list(epsilon = 1e-10, maxit = 200))
    max(abs(unname(coef(f)) - unname(coef(b))))
})
difference_table(differences)

## -----------------------------------------------------------------------------
mu_ig <- exp(eta_true)
y_ig <- statmod::rinvgauss(n, mean = mu_ig, shape = 5)
y_ig[y_ig <= 0] <- 0.01

f <- fastglm(X, y_ig, family = inverse.gaussian("log"), firth = TRUE, tol = 1e-10)
b <- glm(y ~ x1 + x2, data = data.frame(y = y_ig, x1 = x1, x2 = x2),
         family = inverse.gaussian("log"), method = "brglmFit", type = "AS_mean",
         control = list(epsilon = 1e-10, maxit = 200))
difference_table(c(log = max(abs(unname(coef(f)) - unname(coef(b))))))

## -----------------------------------------------------------------------------
f <- fastglm(X, y_bin, family = binomial(), firth = TRUE, tol = 1e-10)
b <- glm(y ~ x1 + x2, data = df, family = binomial(), method = "brglmFit",
         type = "AS_mean", control = list(epsilon = 1e-10, maxit = 200))

cbind(fastglm_SE = f$se, brglm2_SE = summary(b)$coefficients[, "Std. Error"])
max(abs(f$se - summary(b)$coefficients[, "Std. Error"]))

## -----------------------------------------------------------------------------
fit_gamma <- fastglm(X, y_gam, family=Gamma("log"), firth=TRUE,
                    tol=1e-10)
fastglm_diagnostics(fit_gamma)

## -----------------------------------------------------------------------------
housing <- MASS::housing
housing_design <- model.matrix(~Infl+Type+Cont,housing)
housing_response <- with(housing,cbind(Freq*(Sat=="Low"),Freq*(Sat!="Low")))
fit_housing <- fastglm(housing_design,housing_response,
    family=binomial("cloglog"),firth=TRUE,tol=1e-10)
fastglm_diagnostics(fit_housing)

## -----------------------------------------------------------------------------
set.seed(112101)
separation_design <- cbind("(Intercept)"=1,
    exposure=seq(-2,2,length.out=240),nuisance=rnorm(240,sd=.25))
separated_fit <- fastglm(separation_design,
    as.numeric(separation_design[,"exposure"]>0),
    family=binomial("cloglog"),firth=TRUE,
    weights=rep(c(0,1,2,1),60),offset=cos(seq_len(240))/8,
    tol=1e-10,maxit=500)
summary(separated_fit)
fastglm_diagnostics(separated_fit)

## -----------------------------------------------------------------------------
gamma_blocks <- split(seq_len(nrow(X)),ceiling(seq_len(nrow(X))/250L))
gamma_reader <- local({
    design <- X
    response <- y_gam
    blocks <- gamma_blocks
    function(k) {
        i <- blocks[[k]]
        list(X=design[i,,drop=FALSE],y=response[i])
    }
})
qr_gamma <- fastglm_streaming(gamma_reader,length(gamma_blocks),Gamma("log"),
                              method="qr",firth=TRUE,tol=1e-8)
fastglm_diagnostics(qr_gamma)

## -----------------------------------------------------------------------------
f <- fastglm(X, y_bin, family=binomial(), firth=TRUE)
c(deviance=f$deviance, penalized.deviance=f$penalized.deviance)

## -----------------------------------------------------------------------------
fit_profile <- fastglm(case ~ age + oc + vic + vicl + vis + dia,
                       data=sex2, family=binomial(), firth=TRUE,
                       tol=1e-10, maxit=200)
profile_ci <- confint(fit_profile, parm=c("age","dia"), method="profile")
cbind(estimate=coef(fit_profile)[rownames(profile_ci)], profile_ci)
attr(profile_ci,"diagnostics")[,c("parameter","side","converged","lr.error")]

## -----------------------------------------------------------------------------
fastglm_profile_test(fit_profile, parm=c("oc","vic"), null=0)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 10000;  p <- 5
x1 <- rnorm(n);  x2 <- rnorm(n); x3 <- rnorm(n)
x4 <- rnorm(n);  x5 <- rnorm(n)
X <- cbind("(Intercept)"=1, x1, x2, x3, x4, x5)
eta <- 0.5 + 0.3 * x1 - 0.2 * x2 + 0.1 * x3

families <- list(
    `binomial logit`   = list(fam = binomial("logit"),
                              y = rbinom(n, 1, plogis(eta))),
    `binomial probit`  = list(fam = binomial("probit"),
                              y = rbinom(n, 1, plogis(eta))),
    `binomial cloglog` = list(fam = binomial("cloglog"),
                              y = rbinom(n, 1, plogis(eta))),
    `poisson log`      = list(fam = poisson("log"),
                              y = rpois(n, exp(eta))),
    `Gamma log`        = list(fam = Gamma("log"),
                              y = rgamma(n, shape = 5, rate = 5 / exp(eta))),
    `gaussian identity`= list(fam = gaussian("identity"),
                              y = eta + rnorm(n, 0, 0.5))
)

results <- list()
for (nm in names(families)) {
    fi <- families[[nm]]
    fit_fast <- function() fastglm(X,fi$y,family=fi$fam,firth=TRUE,
                                   tol=1e-8,maxit=100)
    fit_reference <- function() brglm2::brglmFit(X,fi$y,family=fi$fam,
        control=brglm2::brglmControl(type="AS_mean",epsilon=1e-8,maxit=100))
    f <- fit_fast(); b <- fit_reference()
    error <- max(abs(unname(coef(f))-unname(coef(b))))
    stopifnot(f$converged,b$converged,error<1e-5)
    mb <- microbenchmark::microbenchmark(
        fastglm=fit_fast(),brglm2=fit_reference(),times=5L)
    timing <- summary(mb,unit="ms")
    ms <- setNames(timing$median,as.character(timing$expr))
    results[[nm]] <- data.frame(family=nm,fastglm_ms=ms[["fastglm"]],
        brglm2_ms=ms[["brglm2"]],ratio=ms[["brglm2"]]/ms[["fastglm"]],
        fastglm_iterations=f$iter,brglm2_iterations=b$iter,
        max_coefficient_difference=error)
}
timing_table <- do.call(rbind,results)
timing_table$max_coefficient_difference <- formatC(
    timing_table$max_coefficient_difference,format="e",digits=2)
knitr::kable(timing_table,row.names=FALSE,digits=3,
    col.names=c("Family and link","fastglm (ms)","brglm2 (ms)",
                "Time ratio","fastglm iterations","brglm2 iterations",
                "Largest coefficient difference"))

