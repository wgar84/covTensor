MeanMatrix <- function (matArray, tol = 1e-14, max.steps = 100, parallel = FALSE)
  {
      ### Computes (geometric) mean matrix using gradient optimization
      ### in the manifold of PD-Sym matrices (from Woods, 2003)
      logm.single <- function (Ai, inv.Mk) return (logm (Ai %*% inv.Mk))
      A <- matArray
      N <- dim (A) [3]
      Mk <- diag (nrow (A))
      i <- 1
      repeat
          {
              inv.Mk <- solve (Mk)
              centered.now <- aaply (A, 3, logm.single, inv.Mk = inv.Mk,
                                     .parallel = parallel)
              centered.now <- aperm(centered.now, c(2, 3, 1))
              o <- array (- rowMeans (centered.now), c(nrow (A), nrow (A)))
              frob.norm.o <- FrobNorm (o)
              print (frob.norm.o)
              if (frob.norm.o < tol)
                  break
              Mk <- expm (- o) %*% Mk
              if (i == max.steps)
                  stop ('Convergence has not been achieved in number of steps.')
              i <- i + 1
          }
      return (Mk)
  }

CheapMean <- function (mat.array, parallel = FALSE, tol = 1e-08)
  {
    ## DA Bini
    logm.single <- function (Aj, Ai.is) return (logm (Ai.is %*% Aj %*% Ai.is))
    n.mat <- dim (mat.array) [3]
    A.nu <- mat.array
    repeat
      {
        A.nu <-
          aaply (1:n.mat, 1, function (i)
                 {
                   sqA.nu <- sqrtm(A.nu[, , i])
                   arr.log <-
                     aaply (A.nu, 3, logm.single, Ai.is = solve (sqA.nu))
                   arr.log <- aperm (arr.log, c(2, 3, 1))
                   arr.sum <- expm (aaply (arr.log, c(1, 2), mean))
                   sqA.nu %*% arr.sum %*% sqA.nu
                 }, .parallel = parallel)
        A.nu <- aperm (A.nu, c(2, 3, 1))
        dist.mat <- Pairwise(A.nu, RiemannDist, parallel = parallel)
        print (mean (dist.mat))
        if (mean (dist.mat) < tol)
          break
      }
    A.nu [, , 1]
  }
