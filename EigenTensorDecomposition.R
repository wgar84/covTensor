EigenTensorDecomposition <- function (matrix.array)
    ### performs eigentensor decomposition on a set of matrices in array form, following Basser \& Pajevic, 2007
    {
        n.traits <- dim (matrix.array) [1]
        n.matrix <- dim (matrix.array) [3]
        Sigma <- BuildSigma (matrix.array)
        n.tensor <- min (n.matrix, dim (Sigma) [1])
        eigen.dec <- eigen (Sigma)
        eigen.matrices <-
            aaply (eigen.dec $ vectors [, 1:n.tensor], 2, ### only 
                   function (x)
                   {
                       eigen.mat <- diag (x [1:n.traits])
                       eigen.mat [upper.tri (eigen.mat)] <- x [(n.traits+1):length (x)] / sqrt (2)
                       eigen.mat <- eigen.mat + t(eigen.mat)
                       diag (eigen.mat) <- diag (eigen.mat) / 2
                       eigen.mat
                   })
        return (list ('values' = eigen.dec $ values [1:n.tensor], 'matrices' = eigen.matrices))
    }
