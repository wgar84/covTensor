ReTransETD <-
    function (etd, mean.matrix, at.sd = c (-1.96, 1.96), n.dim = NULL, log.matrix = TRUE)
    ### reverts eigenmatrices back to the manifold of Sym-PD matrices,
    ### at points defined with respect to the SD of each direction
    ### for n.dim principal eigenmatrices
    {
        if (is.null (n.dim))
            n.dim <- length (etd $ values)
        at.actual.values <- outer (sqrt (etd $ values [1:n.dim]), at.sd)
        sampled.matrices <-
            aaply (1:n.dim, 1, function (i)
                   etd $ matrices [, , i] %o% at.actual.values [i, ])
        mean.sqrt <- sqrtm (mean.matrix)
        if (log.matrix)
            sampled.matrices <-
                aaply (sampled.matrices, c (1, 4), function (x)
                       mean.sqrt %*% expm (x) %*% mean.sqrt)
        else
            sampled.matrices <-
                aaply (sampled.matrices, c (1, 4), function (x)
                       mean.sqrt %*% x %*% mean.sqrt)
        sampled.matrices
    }

CenterMatrices <-
    function (matrix.array, mean.matrix = NULL, log.matrix = TRUE, parallel = TRUE)
    ### center matrices to the mean matrix in the manifold of Sym-PD matrices,
    ### and (optionally) log-transforms then
    {
        if (is.null (mean.matrix))
            mean.matrix <- MeanMatrix (matrix.array)
        center.matrix <- solve (sqrtm (mean.matrix))
        out <- aaply (matrix.array, 3, function (x)
                      center.matrix %*% x %*% center.matrix,
                      .parallel = parallel)
        if (log.matrix)
            out <- aaply (out, 1, logm, .parallel = parallel)
        out <- aperm (out, c (2, 3, 1), resize = TRUE)
        out
    }


