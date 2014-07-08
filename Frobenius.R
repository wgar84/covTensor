FrobDist <- function (A, B)
    ### Computes Frobenius distance between two matrices
    {
        difference <- A - B
        return (sqrt (sum (diag (t (difference) %*% difference))))
    }

FrobNorm <- function (A)
    ### Computes Frobenius norm of a matrix
    sqrt (sum (diag (t (A) %*% A)))
    
FrobInner <- function (A, B)
  ### Computes Frobenius inner product between two matrices
  sum (diag (t (A) %*% B))

FrobCos <- function (A, B)
  ### Computes cosine of the angle between two matrices
  FrobInner (A, B) / FrobNorm (A) * FrobNorm (B)
