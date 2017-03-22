# Some tensor functions snagged from Peter via David's Github
mat <- function(A, k) {
  Ak <- t(apply(A, k, "c"))
  if (nrow(Ak) != dim(A)[k]){
    Ak <- t(Ak)
  }
  return(Ak)
}
amprod <- function(A, M, k) {
  K <- length(dim(A))
  AM <-crossprod(t(M), mat(A, k))
  AMA <- array(AM, dim = c(dim(M)[1], dim(A)[-k]))
  return(aperm(AMA, match(1:K, c(k, (1:K)[-k]))))
}
atrans <- function(A, B) {
  X <- A
  for (k in 1:length(B)) {
    X <- amprod(X, B[[k]], k)
  }
  return(X)
}