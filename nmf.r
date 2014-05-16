#!/usr/bin/env Rscript

library(BatchJobs)
dyn.load("libnmf.so") # make sure library is available

runExample = function() {
#    library(OCplus)
#    A = MAsim.smyth(ng=5000, n=50, p0=0.2)
#    A = (A - min(A) + runif(1,0,1))/10

    A=read.dataset("20+20x1000.gct")
#    A=read.dataset("50+50x5000.gct")
    runNMFinJobs(A, k=c(2:5), num.clusterings=10, maxniter=10000, seed=123, njobs=4)
}

#####################################################################################
#
# Low-level cluster/thread controlling functions and calling of NMF C library
#
#####################################################################################

#FIXME: change num.clusterings to seed that every run gets
doNMF = function(A, k, maxniter, seed=seed, tolerance=1e-4, num.clusterings=NULL) {
    # wrapper function that takes k, num.clusterings and performs maxniter iterations
    #
    # takes:   A - data matrix
    #          k - the number of clusters (matrix rows)
    #          maxniter - number of NMF iterations to perform
    #          seed - the random seed
    #          tolerance - tolerance for convergence checks in updates
    #          num.clusterings - ignored, BatchJobs needs this
    # returns: a list of the results matrices W and H, with iterations performed

    A = as.matrix(A)
    m = dim(A)[1]
    n = dim(A)[2]
    W = as.double(runif(m*k, 0, 1)) # standard random init in (0,1)
    H = as.double(runif(k*n, 0 ,1)) # could use libnmf's generatematrix for it

    # calls to add: nmf_als, mu, neals, alspg, pg
    dyn.load("libnmf.so")
    result = .C("nmf_mu", a=as.double(A), w0=as.double(W), h0=as.double(H),
                   pm=as.integer(m), pn=as.integer(n), pk=as.integer(k), 
                   maxiter=as.integer(maxniter), pTolX=as.double(tolerance), 
                   pTolFun=as.double(tolerance), DUP=F, PACKAGE='libnmf')

    dim(W) = c(m,k)
    dim(H) = c(k,n)

    return(list(W=W, H=H, iter=result$maxiter)) # calc membership in C?; dnorm return value??
}

createJobArray = function(A, k, num.clusterings, maxniter, seed, tolerance=1e-4) {
    # create job array with batchExpandGrid(k, num.clusterings)
    #
    # takes:   A - data matrix
    #          k - number of clusters (matrix rows)
    #          num.clusterings - clusterings to perform for each k
    #          maxniter - iterations to perform in NMF updates
    #          seed - the random seed to use
    # returns: the job registry object

    reg = makeRegistry(id="NMFJobs", file.dir=tempfile(), seed=seed)
    batchExpandGrid(reg, 
                    doNMF, 
                    k=k, 
                    num.clusterings=c(1:num.clusterings),
                    more.args=list(A=A, maxniter=maxniter, seed=seed, tolerance=tolerance))
    return(reg)
}

reduceGridBy = function(reg, by="x", fun) {
    # summarise data over num.clusterings to get consensus clusters
    #
    # takes:   reg - the job registry
    #          by  - the grid axis to *keep* when summarising data
    #          fun - the function to apply along 'by'
    # returns: a list for each k with the consensus results

    # summarise the job ids for each 'by' as matrix rows
    grid = reduceResultsMatrix(reg, fun=function(job, res) job$pars)
#    stopifnot() # number of result=number of submitted TODO

    reduceLevels = grid[,colnames(grid) == by]
    if (nrow(grid) == 1) { # names are only kept if grid has > 1 row, fix that
        names(reduceLevels) = rownames(grid)
    }
    reduceMatrix = do.call(rbind, split(as.integer(names(reduceLevels)), reduceLevels))
#IMPROVEMENT: coordinate with runNMFinJobs to be able to work on finished rows of matrix

    # apply fun to each row of reduceMatrix and return results as a list
    result = list()
    for (i in 1:dim(reduceMatrix)[1]) {
        res = reduceResultsList(reg, ids=reduceMatrix[i,], fun=function(job, res) res)
        result[[rownames(reduceMatrix)[i]]] = fun(res)
    }
    return(result)
}

#####################################################################################
#
# Medium-level calculation functions in R
#
#####################################################################################

runNMFinJobs = function(A, k, num.clusterings, maxniter, seed, njobs=1) {
    if (1 %in% k)
        stop("Need at least two clusters to compute standard deviation")

    reg = createJobArray(A=A, k=k, num.clusterings=num.clusterings, maxniter=maxniter, seed=seed)
    chunked = chunk(getJobIds(reg), n.chunks=njobs, shuffle=TRUE)
    submitJobs(reg, chunked)
    waitForJobs(reg)
# IMPROVEMENT:
# could process done jobs here as soon as they are available
# instead of waiting for all and then processing all
    result = reduceGridBy(reg, by="k", fun=computeConsensusMatrixFromClusterings)
    computeConsensusAndSaveFiles(result)
}

computeConsensusMatrixFromClusterings = function(listOfResults) {
    # summarises clusterings with same k to consensus
    #
    # takes:   listOfResults - list with each item being doNMF result, k=constant
    # returns: summarised data across clusterings 

    # get a list of cluster assignments (row of highest value in H)
    clusterAssignments = lapply(listOfResults, function(x) apply(x$H, 2, order)[1,])

# this saves memory but uses for loops
#    n = length(clusterAssignments[[1]])
#    connect.matrix = matrix(0, n, n)
#    for(l in clusterAssignments)
#        connect.matrix = connect.matrix + outer(l,l, function(x,y) as.integer(x==y))
##       for(e1 in 1:n)
##           for(e2 in 1:n)
##               connect.matrix[e1,e2] = connect.matrix[e1,e2] + as.integer(l[e1]==l[e2]) 

    # create a similarity matrix of each sample with each sample across clusterings
    connect.matrix = Reduce('+', lapply(clusterAssignments, function(l) 
                                        outer(l,l, function(x,y) as.integer(x==y))))

    return(connect.matrix/length(clusterAssignments)) # do I need more here? if so, return a list
}

computeConsensusAndSaveFiles = function(resultList) {
    k.vec = names(resultList)
    num.k = length(k.vec)
    k.init = as.integer(k.vec[1])
    k.final = as.integer(k.vec[num.k])
    cols = dim(resultList[[1]])[1]
    connect.matrix.ordered <- array(0, c(num.k, cols, cols))

    rho <- vector(mode = "numeric", length = num.k) #TODO: BROAD code
    k.vector <- vector(mode = "numeric", length = num.k) #TODO: BROAD code
    col.names <- as.character(c(1:cols)) #names(D) #???
    doc.string = ""
    directory = "./temp"
    input.ds = "myfile"

    for(k.index in 1:num.k) {
        k = k.vec[k.index]
        connect.matrix = resultList[[as.character(k)]]

        dist.matrix = as.dist(1 - connect.matrix)
        HC = hclust(dist.matrix, method="average")
        dist.coph = cophenetic(HC)

################################ TODO: below is broad code and it's pretty bad
        k.vector[k.index] <- k
        rho[k.index] <- cor(dist.matrix, dist.coph)
        rho[k.index] <- signif(rho[k.index], digits = 4)

        connect.matrix.ordered[k.index,,] <- connect.matrix[HC$order, HC$order]

        # compute consensus clustering membership
        membership <- cutree(HC, k = k)
        max.k <- max(membership)
        items.names.ordered <- col.names[HC$order]
        membership.ordered <- membership[HC$order]
        results <- data.frame(cbind(membership.ordered, items.names.ordered))

        if (k > k.init) {
            all.membership <- cbind(all.membership, membership)
        }
        else {
            all.membership <- cbind(membership)
        }

        sub.string <- paste(doc.string, " k=", k, sep="")
        matrix.abs.plot(connect.matrix.ordered[k.index,,], sub=sub.string, log = F, main = "Ordered Consensus Matrix", ylab = "samples", xlab ="samples")
        plot(HC, xlab="samples", cex = 0.75, labels = col.names, sub = sub.string, col = "blue", main = paste("Ordered Linkage Tree. Coph=", rho[k.index]))
        # dev.off()???

        resultsGct <- data.frame(membership.ordered)
        row.names(resultsGct) <- items.names.ordered
        filename <- paste(directory, doc.string, ".", "consensus.k.",k, ".gct", sep="", collapse="")
        write.gct(resultsGct, filename)

#        H.sorted <- H.saved[,HC$order]
#        sub.string <- paste(doc.string, " k=", k, sep="")
#        matrix.abs.plot(H.sorted, sub = sub.string, log = F, main = "Example H matrix (ordered)", ylab = "metagenes", xlab ="samples")
#        metagene.plot(H = H.sorted, sub = sub.string, main = "Example metagenes (ordered)", xlab = "samples", ylab = "metagenes")
#        dev.off()

        filename <- paste(directory, doc.string, ".", "consensus.plot.k", k, ".pdf", sep="", collapse="")
        pdf(file=filename, width = 8.5, height = 11)
        nf <- layout(matrix(c(1), 1, 1, byrow=T), c(1, 1), c(1, 1), TRUE)
        conlabel <- paste("Consensus k =", k, sep=" ", collapse="")
        sub.string <- paste("Consensus matrix k=", k, "; dataset= ", input.ds, sep="")
        ConsPlot(connect.matrix.ordered[k.index,,], col.labels = membership.ordered, col.names = items.names.ordered, main = " ", sub=sub.string, xlab=" ", ylab=" ")
        dev.off()
    } # end of loop over k


    # Save consensus matrices in one file
    filename <- paste(directory, doc.string, ".", "consensus.all.k.plot.pdf", sep="")
    pdf(file=filename, width = 8.5, height = 11)

    nf <- layout(matrix(c(1:16), 4, 4, byrow=T), c(1, 1, 1, 1), c(1, 1, 1, 1), TRUE)

    for (k in 1:num.k) { 
        matrix.abs.plot(connect.matrix.ordered[k,,], log = F, main = paste("k=", k.vector[k]), 
            sub = paste("Cophenetic coef.=", rho[k]), ylab = "samples", xlab ="samples")
    }

    y.range <- c(1 - 2*(1 - min(rho)), 1)
    plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
         xlab = "k", ylab="Cophenetic correlation", type = "n")
    lines(k.vector, rho, type = "l", col = "black")
    points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
    dev.off()

    filename <- paste(directory, doc.string, ".", "cophenetic.plot.pdf", sep="")
    pdf(file=filename, width = 8.5, height = 11)


    # Write the membership matrix
    resultsmembership <- data.frame(all.membership)
    row.names(resultsmembership) <- col.names
    filename <- paste(directory, doc.string, ".", "membership", ".gct", sep="", collapse="")
    write.gct(resultsmembership , filename)

    y.range <- c(1 - 2*(1 - min(rho)), 1)
    plot(k.vector, rho, main ="Cophenetic Coefficient", xlim=c(k.init, k.final), ylim=y.range, 
         xlab = "k", ylab="Cophenetic correlation", type = "n")
    lines(k.vector, rho, type = "l", col = "black")
    points(k.vector, rho, pch=22, type = "p", cex = 1.25, bg = "black", col = "black")
    dev.off()

    xx <- cbind(k.vector, rho)
    write(xx, file= paste(directory, doc.string, ".", "cophenetic.txt", sep=""))
}

#####################################################################################
#
# Plotting routines
#
#####################################################################################

read.dataset <- function(file) {
	result <- regexpr(paste(".gct","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(read.gct(file))
	result <- regexpr(paste(".res","$",sep=""), tolower(file))
	if(result[[1]] != -1)
		return(read.res(file))
	stop("Input is not a res or gct file.")	
}

matrix.abs.plot <- function(V, axes = F, log = F, transpose = T, matrix.order = T, max.v = 1, min.v = 0, main = " ", sub = " ", xlab = " ", ylab = "  ") {
    if (log == T)
        V = log(V)

    if (max.v == 1 && min.v == 0) {
        max.v = max(V)
        min.v = min(V)
    }

    B = max.v - V + min.v
    if (matrix.order == T)
        B = apply(B, 2, rev)

	if (transpose == T)
	    B <- t(B)
    
    image(z = B, zlim = c(min.v, max.v), axes = axes, 
          col = rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), 
          main = main, sub = sub, xlab = xlab, ylab = ylab)

    return(list(B, max.v, min.v))
}

metagene.plot <- function(H, main = " ", sub = " ", xlab = "samples ", ylab = "amplitude") {
	k <- length(H[,1])
	S <- length(H[1,])
	index <- 1:S

	plot(index, H[1,], xlim=c(1, S), ylim=c(min(H), max(H)), main = main, 
         sub = sub, ylab = ylab, xlab = xlab, type="n")

	for (i in 1:k)
	    lines(index, H[i,], type="l", col = i, lwd=2)
}

ConsPlot <- function(V, col.labels, col.names, main = " ", sub = " ", xlab=" ", ylab=" ") {
    # Plots a heatmap plot of a consensus matrix
    cols <- length(V[1,])
	B = max(V) - V + min(V)
    B = apply(B, 1, rev)

    col.names2 <- rev(col.names)
    col.labels2 <- rev(col.labels)
    D <- matrix(0, nrow=(cols + 1), ncol=(cols + 1))

    col.tag <- vector(length=cols, mode="numeric")
    current.tag <- 0
    col.tag[1] <- current.tag
    for (i in 2:cols) {
        if (col.labels[i] != col.labels[i - 1]) {
            current.tag <- 1 - current.tag
        }
        col.tag[i] <- current.tag
    }

    col.tag2 <- rev(col.tag)
    D[(cols + 1), 2:(cols + 1)] <- ifelse(col.tag %% 2 == 0, 1.02, 1.01)
    D[1:cols, 1] <- ifelse(col.tag2 %% 2 == 0, 1.02, 1.01)
    D[(cols + 1), 1] <- 1.03
    D[1:cols, 2:(cols + 1)] <- B[1:cols, 1:cols]

    col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75), 
                 "#BBBBBB", "#333333", "#FFFFFF")

    image(1:(cols + 1), 1:(cols + 1), t(D), col = col.map, axes=FALSE, 
          main=main, sub=sub, xlab= xlab, ylab=ylab)

    col.names <- sapply(col.names, function(x) paste("      ", substr(x, 1, 12), sep=""))
    col.names2 <- sapply(col.names2, function(x) paste(substr(x, 1, 12), "     ", sep=""))

    axis(2, at=1:cols, labels=col.names2, adj= 0.5, tick=FALSE, las = 1, 
         cex.axis=0.50, font.axis=1, line=-1)
    axis(2, at=1:cols, labels=col.labels2, adj= 0.5, tick=FALSE, las = 1, 
         cex.axis=0.65, font.axis=1, line=-1)
    axis(3, at=2:(cols + 1), labels=col.names, adj= 1, tick=FALSE, las = 3, 
         cex.axis=0.50, font.axis=1, line=-1)
    axis(3, at=2:(cols + 1), labels=as.character(col.labels), adj = 1, 
         tick=FALSE, las = 1, cex.axis=0.65, font.axis=1, line=-1)
}

read.res <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in RES format and converts it into an R data frame
#
   header.cont <- readLines(filename, n = 1)
   temp <- unlist(strsplit(header.cont, "\t"))
   colst <- length(temp)
   header.labels <- temp[seq(3, colst, 2)]
   ds <- read.delim(filename, header=F, row.names = 2, sep="\t", quote="", 
                    skip=3, blank.lines.skip=T, comment.char="", as.is=T)
   colst <- length(ds[1,])
   cols <- (colst - 1)/2
   rows <- length(ds[,1])
   A <- matrix(nrow=rows - 1, ncol=cols)
   A <- ds[1:rows, seq(2, colst, 2)]
   table1 <- data.frame(A)
   names(table1) <- header.labels
   return(table1)
}

read.gct <- function(filename = "NULL") { 
    # Reads a gene expression dataset in GCT format and converts it into an R data frame
    ds <- read.delim(filename, header=T, sep="\t", quote="", skip=2, row.names=1, 
                     blank.lines.skip=T, comment.char="", as.is=T)
    ds <- ds[-1]
    return(ds)
}

write.gct <- function (gct, filename) {
    f <- file(filename, "w")
    cat("#1.2", "\n", file = f, append = TRUE, sep = "")
    cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")

    cat("Name", "\t", file = f, append = TRUE, sep = "")
    cat("Description", "\t", file = f, append = TRUE, sep = "")
    cat(c(1:dim(gct)[2]), file = f, append = TRUE, sep = "\t")

    names <- names(gct)
    for (j in 1:length(names)) {
        cat("\t", names[j], file = f, append = TRUE, sep = "")
    }
    cat("\n", file = f, append = TRUE, sep = "")
    oldWarn <- options(warn = -1)

    m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  2)
    m[, 1] <- rownames(gct)
    m[, 2] <- rownames(gct)
    index <- 3
    for (i in 1:dim(gct)[2]) {
        m[, index] <- gct[, i]
        index <- index + 1
    }
    write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", 
                col.names = FALSE, row.names = FALSE)
    close(f)
    options(warn = 0)
    return(gct)
}

