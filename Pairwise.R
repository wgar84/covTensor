Pairwise <-
  function (A, func, parallel = FALSE, ...)
  UseMethod ('Pairwise')

Pairwise.array <-
  function (A, func, parallel = FALSE, ...)
  {
    coerce2Dist <-
      function (x)
        {
          pair.mat <- array (0, c (n.mat, n.mat))
          dimnames (pair.mat) <- list (dimnames (A) [[3]], dimnames (A) [[3]])
          pair.mat [lower.tri (pair.mat)] <- x
          pair.mat <- as.dist (pair.mat)
          return (pair.mat)
        }
    n.mat <- dim (A) [3]
    index <- which (lower.tri (diag (n.mat)), arr.ind = TRUE)
    pair.comp <-
      aaply (index, 1, function (ij) func (A [, , ij [1]], A [, , ij [2]], ...),
             .parallel = parallel)
    if (is.null (dim (pair.comp)))
      return (coerce2Dist (pair.comp))
    else
      return (alply (pair.comp, 2, coerce2Dist, .dims = TRUE))
  }

Pairwise.list <-
  function (A, func, parallel = FALSE, ...)
  {
    coerce2Dist <-
      function (x)
        {
          pair.mat <- array (0, c (n.mat, n.mat))
          dimnames (pair.mat) <- list (names (A), names (A))
          pair.mat [lower.tri (pair.mat)] <- x
          pair.mat <- as.dist (pair.mat)
          return (pair.mat)
        }
    n.mat <- length (A)
    index <- which (lower.tri (diag (n.mat)), arr.ind = TRUE)
    pair.comp <-
      aaply (index, 1, function (ij) func (A [[ij [1]]], A [[ij [2]]], ...),
             .parallel = parallel)
    if (is.null (dim (pair.comp)))
      return (coerce2Dist (pair.comp))
    else
      return (alply (pair.comp, 2, coerce2Dist, .dims = TRUE))
  }
