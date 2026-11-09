## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>", message=FALSE,
                      eval=requireNamespace("MASS",quietly=TRUE))

## -----------------------------------------------------------------------------
library(fastglm)
set.seed(82)
n <- 1200
raw <- data.frame(x=rnorm(n),g=factor(rep(letters[1:4],length.out=n)))
raw$y <- rpois(n,exp(.2+.3*raw$x))
X <- model.matrix(~x+g,raw)

## -----------------------------------------------------------------------------
Xs <- as(X,"dgCMatrix")
sparse_fit <- fastglm(Xs,raw$y,poisson(),method="llt")
dense_fit <- fastglm(X,raw$y,poisson(),method="llt")
max(abs(coef(sparse_fit)-coef(dense_fit)))

## -----------------------------------------------------------------------------
Xb <- bigmemory::as.big.matrix(X)
mapped_fit <- fastglm(Xb,raw$y,poisson(),method="llt")
max(abs(coef(mapped_fit)-coef(dense_fit)))

## -----------------------------------------------------------------------------
chunk_rows <- 300L
get_data <- function(k) raw[((k-1L)*chunk_rows+1L):(k*chunk_rows),]
source <- fastglm_source(get_data, n_chunks=4L,
                         formula=y~x+g, prototype=raw)
stream_fit <- fastglm_streaming(source, family=poisson())
max(abs(coef(stream_fit)-coef(dense_fit)))
c(passes=stream_fit$passes, chunks=stream_fit$chunk.reads)
summary(stream_fit)

## -----------------------------------------------------------------------------
qr_stream <- fastglm_streaming(source,family=poisson(),method="qr")
max(abs(coef(qr_stream)-coef(stream_fit)))
head(predict(qr_stream,type="response",se.fit=TRUE)$se.fit)

## -----------------------------------------------------------------------------
set.seed(149)
positive <- raw
positive$y <- rgamma(n,shape=2,scale=exp(.2+.3*positive$x)/2)
positive_source <- fastglm_source(
    function(k) positive[((k-1L)*chunk_rows+1L):(k*chunk_rows),],
    4L,formula=y~x+g,prototype=positive)
bias_stream <- fastglm_streaming(positive_source,Gamma("log"),
                                 firth=TRUE,tol=1e-9)
bias_dense <- fastglm(X,positive$y,Gamma("log"),firth=TRUE,tol=1e-9)
max(abs(coef(bias_stream)-coef(bias_dense)))
c(passes=bias_stream$passes,blocks_read=bias_stream$chunk.reads)

## -----------------------------------------------------------------------------
chunk_means <- numeric(4)
predict(stream_fit,type="response",
        sink=function(value,k) chunk_means[k] <<- mean(value))
chunk_means
fastglm_average(stream_fit)

## -----------------------------------------------------------------------------
cluster_ids <- rep(seq_len(120),each=10)
get_cluster <- function(k) cluster_ids[((k-1L)*chunk_rows+1L):(k*chunk_rows)]
V <- fastglm_vcov(stream_fit,type="HC1",cluster=get_cluster)
sqrt(diag(V))
fastglm_average(stream_fit,vcov.=V)

## -----------------------------------------------------------------------------
housing <- MASS::housing
housing_source <- fastglm_source(
    function(k) housing[((k-1L)*18L+1L):(k*18L),,drop=FALSE],
    n_chunks=4L,formula=Sat~Infl+Type+Cont,
    prototype=housing,weights="Freq")
ordinal_stream <- fastglm_ordinal_streaming(housing_source,firth=TRUE)
summary(ordinal_stream)
fastglm_diagnostics(ordinal_stream)

## -----------------------------------------------------------------------------
ordinal_partial_stream <- fastglm_ordinal_streaming(
    housing_source,nonparallel=~Cont,firth=TRUE)
summary(ordinal_partial_stream)

## -----------------------------------------------------------------------------
category_totals <- matrix(0,4L,nlevels(housing$Sat),
    dimnames=list(paste0("block",1:4),levels(housing$Sat)))
predict(ordinal_stream,sink=function(value,k) {
    i <- ((k-1L)*18L+1L):(k*18L)
    category_totals[k,] <<- colSums(value*housing$Freq[i])
})
category_totals

## ----ordinal-source-covariance------------------------------------------------
V_ordinal <- fastglm_vcov(ordinal_partial_stream,type="HC0")
fastglm_average(ordinal_partial_stream,newdata=housing,
                weights=housing$Freq,vcov.=V_ordinal)

## ----ordinal-source-clusters--------------------------------------------------
set.seed(110305)
ordered_data <- raw
ordered_data$rating <- ordered(cut(.5*ordered_data$x+rlogis(n),
    c(-Inf,-.5,.75,Inf),labels=c("low","middle","high")))
ordinal_clusters <- rep(seq_len(120),length.out=n)
ordinal_source <- fastglm_source(
    function(k) ordered_data[((k-1L)*chunk_rows+1L):(k*chunk_rows),],
    n_chunks=4L,formula=rating~x+g,prototype=ordered_data)
ordinal_cluster_fit <- fastglm_ordinal_streaming(
    ordinal_source,nonparallel=~x,firth=TRUE)
ordinal_cluster_ids <- function(k)
    ordinal_clusters[((k-1L)*chunk_rows+1L):(k*chunk_rows)]
V_ordinal_cluster <- fastglm_vcov(
    ordinal_cluster_fit,type="HC1",cluster=ordinal_cluster_ids)
summary(ordinal_cluster_fit,vcov.=V_ordinal_cluster)

## ----ordinal-source-assessment------------------------------------------------
ordinal_assessment <- fastglm_ordinal_assess(ordinal_partial_stream,bins=5)
ordinal_assessment

## ----ordinal-source-assessment-output-----------------------------------------
score_directory <- tempfile("ordinal-scores-")
dir.create(score_directory)
save_scores <- function(value,k)
    saveRDS(value,file.path(score_directory,paste0("block-",k,".rds")))
fastglm_ordinal_assess(ordinal_partial_stream,sink=save_scores)

## ----ordinal-source-assessment-cleanup, include=FALSE-------------------------
unlink(score_directory,recursive=TRUE)

## ----ordinal-source-po-test---------------------------------------------------
fastglm_po_test(ordinal_partial_stream,terms="Cont")

## ----ordinal-source-boundary-diagnostics--------------------------------------
fastglm_diagnostics(ordinal_partial_stream,boundary=TRUE,max.rows=3)

## ----eval=FALSE---------------------------------------------------------------
# parquet_source <- fastglm_source_arrow(
#     files=c("part-1.parquet","part-2.parquet"),
#     formula=y~x+g, prototype=prototype)
# fit <- fastglm_streaming(parquet_source, family=poisson())
# 
# queries <- c("SELECT y, x, g FROM outcomes WHERE id <= 10000 ORDER BY id",
#              "SELECT y, x, g FROM outcomes WHERE id > 10000 ORDER BY id")
# database_source <- fastglm_source_dbi(connection,queries,
#                                       formula=y~x+g,prototype=prototype)

## -----------------------------------------------------------------------------
hurdle_source <- fastglm_source(get_data,4,
                                 formula=y~x+g|x,prototype=raw)
hurdle_stream <- fastglm_hurdle(hurdle_source,dist="poisson")
fastglm_average(hurdle_stream)

