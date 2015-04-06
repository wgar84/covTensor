EigenTensorDecomposition <-
    function (matrix.array, return.projection = FALSE)
    ### performs eigentensor decomposition on a set of matrices in array form,
    ### following Basser \& Pajevic, 2007
    {
        n.traits <- dim (matrix.array) [1]
        n.matrix <- dim (matrix.array) [3]
        Sigma <- BuildSigma (matrix.array)
        n.tensor <- min (n.matrix - 1, dim (Sigma) [1])
        eigen.dec <- eigen (Sigma)
        eigen.matrices <-
            aaply (eigen.dec $ vectors [, 1:n.tensor], 2, ### ematrices with non-zero evals
                   function (x)
                   {
                       eigen.mat <- diag (x [1:n.traits])
                       eigen.mat [upper.tri (eigen.mat)] <-
                           x [(n.traits+1):length (x)] / sqrt (2)
                       eigen.mat <- eigen.mat + t(eigen.mat)
                       diag (eigen.mat) <- diag (eigen.mat) / 2
                       eigen.mat
                   })
        eigen.matrices <- aperm (eigen.matrices, c(2, 3, 1), resize = TRUE)
        dimnames (eigen.matrices) <-
            list (rownames (matrix.array),
                  colnames (matrix.array),
                  paste ('PM', 1:n.tensor, sep = ''))
        out <- list ('values' = eigen.dec $ values [1:n.tensor],
                     'matrices' = eigen.matrices,
                     'Sigma' = Sigma)
        if (return.projection)
            {
                project <- aaply (matrix.array, 3, function (A, B)
                                  aaply (B, 3, FrobInner, B = A),
                                  B = eigen.matrices)
                out $ projection <- project
                out $ projection <- scale (out $ projection, center = FALSE, scale = TRUE)
            }
        return (out)
    }

