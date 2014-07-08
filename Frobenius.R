FrobDist <- function (A, B)
    ### Computes Frobenius distance between two matrices
    {
        difference <- A - B
        return (sqrt (sum (diag (t (difference) %*% difference))))
    }

FrobNorm <- function (A)
    ### Computes Frobenius norm of a matrix
    sqrt (sum (diag (t (A) %*% A)))
    
