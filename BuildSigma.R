BuildSigma <- function (matrix.array)
    ### constructs Sigma matrix (second-order tensor representation of a fourth-order covariance
    ### tensor between symmetric matrices), following Hine et al., 2009
    {
        variances <- aaply (matrix.array, 3, diag)
        covariances <- aaply (matrix.array, 3, function (x) x [lower.tri (x)])
        block.var <- var (variances)
        block.off <- cov (variances, covariances)
        block.cov <- var (covariances)
        upper.block <- cbind (block.var, sqrt (2) * block.off)
        lower.block <- cbind (sqrt (2) * t (block.off), 2 * block.cov)
        Sigma <- rbind (upper.block, lower.block)
        return (Sigma)
    }
