RetransETD <-
    function (etd, mean.matrix, at.sd = c (-1.96, 1.96), n.dim = NULL)
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
        sampled.matrices <-
            aaply (sampled.matrices, c (1, 4), function (x)
                   mean.sqrt %*% expm (x) %*% mean.sqrt)
        sampled.matrices
    }
