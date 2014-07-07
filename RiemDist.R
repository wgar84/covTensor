RiemDist <- function (A, B)
    ### Computes Riemannian distance between two matrices, following Moakher, 2005.
    sqrt (sum (log (eigen (solve (A) %*% B) $ values) ^ 2))
    
