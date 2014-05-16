library(OCplus)
A = MAsim.smyth(ng=5000, n=50, p0=0.2)
A = (A - min(A) + runif(1,0,1))/10

#rownames(A) = c(1:dim(A)[1])
#colnames(A) = c(1:dim(A)[2])
#write.gct(A, "50+50x5000.gct")

plotPCA = function(my.exprs, my.title="PCA plot", my.classes=NULL) {
    # my.exprs = row-major (gct is column major)
    my.svd = svd(cov(my.exprs))
    my.iloads = solve(t(my.svd$v))
    my.weights = my.svd$d
    my.scores = my.exprs %*% my.iloads

    my.title = paste("PCA of different leukemia types, capturing\n", 
        round(100*sum(my.weights[1:2])/sum(my.weights), digits=1), "% of variance")
    my.pc1 = paste("PC1 (", round(100*my.weights[1]/sum(my.weights), 1), "%)")
    my.pc2 = paste("PC2 (", round(100*my.weights[2]/sum(my.weights), 1), "%)")

    plot(my.scores[,1], my.scores[,2], main=my.title, xlab=my.pc1, ylab=my.pc2)#, pch=my.classes, col=my.classes)
#    legend(x="bottomright", legend=my.types, col=c(1:length(my.types)), pch=c(1:length(my.types)))
}

source("nmf.r")
A = read.gct("20+20x1000.gct")
system.time(runNMFinJobs(A, k=c(2:5), num.clusterings=5, maxniter=10000, seed=123, njobs=1))

#system.time(nmfconsensus("50+50x5000.gct", k.init=2, k.final=5, num.clusterings=5, maxniter=10000, error.function="euclidean"))
