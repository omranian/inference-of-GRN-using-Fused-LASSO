
# this function should be substituted with the fused.lasso function in lqa in the following manner:
# 1: trace(fused.lasso,edit=T)
# 2: substitute the following function with the existing  one
function (lambda = NULL, ...) 
{
  argList = list(...)
  w <- argList$...
  lambda.check(lambda)
  if (length(lambda) != 2) 
    stop("The fused.lasso penalty must consist on two parameters! \n")
  names(lambda) <- c("lambda1", "lambda2")
  first.derivative <- function(beta = NULL, ...) {
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    p <- length(beta)
    if (p < 2) 
      stop("There must be at least two regressors! \n")
    vec1 <- c(rep(lambda[1], p), rep(lambda[2], (3/4) * p))
    vec2 <- abs(drop(t(a.coefs(beta, ... = w)) %*% beta))
    return(vec1 * vec2)
  }
  a.coefs <- function(beta = NULL, ...) {
    argList = list(...)
    w <- argList$...
    if (is.null(beta)) 
      stop("'beta' must be the current coefficient vector \n")
    p <- length(beta)
    if (p < 2) 
      stop("There must be at least two regressors! \n")
    if (p > 2) {
      h1 <- cbind(-diag((3/4) * p), matrix(0, (3/4) * p, (1/4) * p))
      h2 <- cbind(matrix(0, (3/4) * p, (1/4) * p), diag((3/4) * p))
      mat1 <- h1 + h2
      mat2 <- diag(w)
      a.coefs.mat <- cbind(mat2, t(mat1))
    }
    else a.coefs.mat <- cbind(diag(2), c(-1, 1))
    return(a.coefs.mat)
  }
  structure(list(penalty = "fused.lasso", lambda = lambda, 
                 first.derivative = first.derivative, a.coefs = a.coefs), class = "penalty")
}
