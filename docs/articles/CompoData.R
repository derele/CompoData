## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )
library(ggplot2)
theme_set(theme_bw())

## ---- fig.height=4, fig.show='hold', fig.width=4-------------------------
library(reshape)
## initialize the random number generator (-> reproducibility)
set.seed(1209)

n <- 3 ## number of rows
m <- 3 ## number of columns

means <- c(8, 28, 330) ## per row (species) mean

countData <- matrix(rnbinom(m * n, mu = means, size = 10), 
                    ncol = m)

colnames(countData) <- c("venus", "earth", "mars")
rownames(countData) <- c("tigers", "unicorns", "ladybugs")

ggplot(melt(countData),
       aes(x=X2, y=value, fill=X1)) +
    geom_col()

## ---- fig.show='hold', fig.height=4, fig.width=4-------------------------
depth <- c(1, 0.1, 0.5) ## per column (planet) depth
mu <- as.matrix(means)%*%t(as.matrix(depth))

countData <- matrix(rnbinom(m * n, mu = mu, size = 10), 
                   ncol = m)
colnames(countData) <- c("venus", "earth", "mars")
rownames(countData) <- c("tigers", "unicorns", "ladybugs")

ggplot(melt(countData),
       aes(x=X2, y=value, fill=X1)) +
    geom_col()

## ---- results='as.is', fig.height=4, fig.show='hold', fig.width=4--------
scaling.factor <- colSums(countData)/max(colSums(countData))

ggplot(melt(t(t(countData)*1/scaling.factor)),
       aes(x=X2, y=value, fill=X1)) +
    geom_col()

## ---- fig.height=4, fig.show='hold', fig.width=4-------------------------
library(pheatmap)
n <- 300 ## number of rows (species)
m <- 36 ## number of columns (samples)

## the relationship between mean and dispersion
dispMeanRel <- function(x) 4/x + 0.1

## low mean and high standard deviation to get sparse counts (lots of
## zeros)
interceptMean <- 0.1
interceptSD <- 2

## the depth as influenced by our sequencing
## depth <- runif(m, min=0.4, max=1)

##  for now assume we have even depths
depth <- rep(1, times=m)

## beta will be the mean of the negative binomial it is here still on
## the log scale we make two more columns wich give a deviation from
## our intercept SD (see below)
beta <- cbind(rnorm(n, interceptMean, interceptSD),
              rnorm(n, 0, 1),
              rnorm(n, 0, 2))

## the sample different ecosystems
colData <- data.frame(planet =
                          factor(rep(c("venus", "earth", "mars"),
                                     each = m/3)))
rownames(colData) <- paste(colData$planet, 1:m, sep="_")

## a model matrix telling us which sample belongs to which condition
planet.design <- stats::model.matrix.default(~colData$planet)

## the final sample means based on design, log-means and the depth of
## sequencing. By matrix multiplication mars gets 
mu <- t(2^(planet.design %*% t(beta)) * depth)

## the dispersion via the mean relationship 
dispersion <- dispMeanRel(2^(beta[, 1]))

## and based of on these the negative binomial counts
countData <- matrix(rnbinom(m * n, mu = mu, size = dispersion), 
                    ncol = m)

## we now remove species for which we have only zero counts (we
## wouldn't have this in a real dataset)
countData <- countData[rowSums(countData)>0, ]

colnames(countData) <- rownames(colData)
rownames(countData) <- paste("species", 1:nrow(countData), sep="_")

pheatmap(log10(countData+1),
         annotation_col=colData,
         show_rownames=FALSE,
         show_colnames=FALSE)

## ---- results='as.is', fig.height=4, fig.show='hold', fig.width=4--------

ColS <- data.frame(sample.sum=colSums(countData), num=1:ncol(countData),
                   planet=colData$planet)

ggplot(ColS, aes(x=reorder(num, sample.sum), y=sample.sum, fill=planet))+
    geom_bar(stat = "identity")

scaling.factor <- colSums(countData)/max(colSums(countData))

countData.scaled <- t(t(countData)*1/scaling.factor)

pheatmap(log10(countData.scaled+1),
         annotation_col=colData,
         show_rownames=FALSE,
         show_colnames=FALSE)

## ---- results='as.is', fig.height=4, fig.show='hold', fig.width=4--------
library(vegan)
countData.rare <- lapply(1:100, function (i){
    t(rrarefy(t(countData),
              sample=min(ColS$sample.sum)))
})

pheatmap(log10(countData.rare[[1]]+1),
         annotation_col=colData,
         show_rownames=FALSE,
         show_colnames=FALSE)

richness <- lapply(countData.rare, function(y){
    apply(y, 2, function(x) length(x[x>0]))
})

rich.mat <- matrix(unlist(richness), nrow=m)

richness <- data.frame(planet=colData$planet,
                       richness.mean=rowMeans(rich.mat),
                       richness.sd=apply(rich.mat, 1, sd),
                       n=1:nrow(rich.mat))

ggplot(richness, aes(as.character(n), richness.mean)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=richness.mean-richness.sd,
                      ymax=richness.mean+richness.sd))+
    facet_wrap(~planet, scales="free_x")

## ---- results='as.is', fig.height=4, fig.show='hold', fig.width=4, cache=TRUE----
library(ALDEx2)


g.mean <- function (x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

countData.clr <- apply(countData, 2, function(x){
    log(x+0.0001/g.mean(x+0.0001))
})

pheatmap(countData.clr,
         annotation_col=colData,
         show_rownames=FALSE,
         show_colnames=FALSE)

countData.clrP <- aldex.clr(countData,
                            conds=as.character(colData$planet))

res.aldex <- aldex.glm(countData.clrP, conditions=colData$planet)

## ---- results='as.is', fig.height=4, fig.show='hold', fig.width=4--------
res <- merge(res.aldex, as.data.frame(res.deseq), by=0)

ggplot(res, aes(glm.eBH, padj)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10()


