ReTransETD <-
    function (etd.values, etd.matrices, mean.matrix,
              at.sd = c (-1.96, 1.96), n.dim = NULL, log.matrix = TRUE)
    ### reverts eigenmatrices back to the manifold of Sym-PD matrices,
    ### at points defined with respect to the SD of each direction
    ### for n.dim principal eigenmatrices
    {
        if (is.null (n.dim))
            n.dim <- length (etd.values)
        at.actual.values <- outer (sqrt (etd.values [1:n.dim]), at.sd)
        sampled.matrices <-
            aaply (1:n.dim, 1, function (i)
                   etd.matrices [, , i] %o% at.actual.values [i, ])
        mean.sqrt <- sqrtm (mean.matrix)
        if (log.matrix)
            sampled.matrices <-
                aaply (sampled.matrices, c (1, 4), function (x)
                       mean.sqrt %*% expm (x) %*% mean.sqrt)
        else
            sampled.matrices <-
                aaply (sampled.matrices, c (1, 4), function (x)
                       mean.sqrt %*% x %*% mean.sqrt)
        dim ()
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
        dimnames (out) <- dimnames (matrix.array)
        out
    }

BuildMatrix <-
    function (scores, etd, mean.matrix, standardized = FALSE, log.matrix = TRUE)
    ### builds matrix using etd as a basis
    {
        n.pm <- length (scores)
        if (standardized)
            scores <- scores * etd $ values [1:n.pm]
        if (n.pm == 1)
            {
                projected <- etd $ matrices [, , 1] * scores
                matrix.sum <- projected
            }
        else
            {
                projected <-
                    aaply (1:n.pm, 1, function (i) etd $ matrices [, , i] * scores [i])
                matrix.sum <- aaply (projected, c (2, 3), sum)
            }
        if (log.matrix)
          matrix.sum <- expm (matrix.sum)
        recenter.matrix <- sqrtm (mean.matrix)
        out <- recenter.matrix %*% matrix.sum %*% recenter.matrix
        dimnames (out) <- dimnames (etd $ matrices) [1:2]
        out
    }

BuildMatrix.NC <-
    function (scores, etd, standardized = FALSE, log.matrix = TRUE)
    ### builds matrix using etd as a basis
    {
        n.pm <- length (scores)
        if (standardized)
            scores <- scores * etd $ values [1:n.pm]
        if (n.pm == 1)
            {
                projected <- etd $ matrices [, , 1] * scores
                matrix.sum <- projected
            }
        else
            {
                projected <-
                    aaply (1:n.pm, 1, function (i) etd $ matrices [, , i] * scores [i])
                matrix.sum <- aaply (projected, c (2, 3), sum)
            }
        if (log.matrix)
          out <- expm (matrix.sum)
        dimnames (out) <- dimnames (etd $ matrices) [1:2]
        out
    }


ProjectMatrix <-
  function (matrix, etd, mean.matrix, inv.sqrt = FALSE)
  {
    if (!inv.sqrt)
      mean.matrix <- solve (sqrtm (mean.matrix))
    log.center.mat <- logm (mean.matrix %*% matrix %*% mean.matrix)
    aaply (etd $ matrices, 3, FrobInner, B = log.center.mat)
  }
