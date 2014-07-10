MatrixDist <-
  function (A, B, method = c('Frobenius', 'Riemann'), ...)
  UseMethod ('MatrixDist')

FrobeniusDist <-
  function (A, B)
    ### Computes Frobenius distance between two matrices
    {
        difference <- A - B
        return (sqrt (sum (diag (t (difference) %*% difference))))
    }

RiemannDist <-
  function (A, B)
    ### Computes Riemannian distance between two matrices, following Moakher, 2005.
    sqrt (sum (log (eigen (solve (A) %*% B) $ values) ^ 2))

MatrixDist.matrix <-
  function (A, B, method = c('Frobenius', 'Riemann'))
  {
    DistFunc <- eval (parse (text = paste (method, 'Dist', sep = '')))
    DistFunc (A, B)
  }

MatrixDist.array <-
  function (A, B = NULL, method = c('Frobenius', 'Riemann'), parallel = F)
  {
    DistFunc <-  eval (parse (text = paste (method, 'Dist', sep = '')))
    if (!is.null (B))
      aaply (A, 3, DistFunc, B = B)
    else
      Pairwise (A, DistFunc, parallel = parallel)
  }

MatrixDist.list <-
  function (A, B = NULL, method = c('Frobenius', 'Riemann'), parallel = F)
  {
    DistFunc <-  eval (parse (text = paste (method, 'Dist', sep = '')))
    if (!is.null (B))
      laply (A, 3, DistFunc, B = B)
    else
      Pairwise (A, DistFunc, parallel = parallel)
  }

