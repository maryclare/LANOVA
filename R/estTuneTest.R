estTuneTest <- function(Y, lowpen = FALSE) {
  
  p <- dim(Y)
  
  H <- vector("list", length(p))
  for (i in 1:length(p)) {
    H[[i]] <- diag(p[i]) - rep(1, p[i])%*%t(rep(1, p[i]))/p[i]
  }
  
  
  r <- atrans(Y, H) 
  bar.r.2 <- mean(r^2)
  bar.r.4 <- mean(r^4)
  
  sigma.4.c <- prod(p^3/((p - 1)*(p^2 - 3*p + 3)))*(bar.r.4/3 - bar.r.2^2)
  sigma.2.cz <- prod(p/(p - 1))*bar.r.2
  sigma.2.c <- sqrt(ifelse(sigma.4.c >= 0, sigma.4.c, 0))
  sigma.2.z <- ifelse(sigma.2.cz - sigma.2.c >= 0, sigma.2.cz - sigma.2.c, 0)
  
  lambda.c <- ifelse(sigma.2.c != 0, sqrt(2/sigma.2.c), Inf)
  
  test.stat <- sqrt(prod(p))*(sigma.4.c/(sqrt(8/3)*(sigma.2.cz)^2))  
  
  out <- list("lambda.c" = lambda.c,
              "sigma.2.z" = sigma.2.z,
              "test.stat" = test.stat)
  
  if (lowpen & length(p) == 2) {
    
    m.1 <- t(rep(1, p[1])/p[1])
    m.2 <- t(rep(1, p[2])/p[2])
    
    H.1 <- H[[1]]
    H.2 <- H[[2]]
    
    ahat <- atrans(Y, list(H.1, m.2))[ , 1]
    bhat <- atrans(Y, list(m.1, H.2))[ , 1]
    
    s2a <- sum(ahat^2)/(p[1] - 1) - p[1]*bar.r.2/((p[1] - 1)*(p[2] - 1))
    s2b <- sum(bhat^2)/(p[2] - 1) - p[2]*bar.r.2/((p[1] - 1)*(p[2] - 1))
    
    out <- c(out, 
             "lambda.a" = ifelse(s2a > 0, sqrt(2/s2a), Inf),
             "lambda.b" = ifelse(s2b > 0, sqrt(2/s2b), Inf))
    
    
    
  } else if (lowpen & length(p) == 3) {
    
    
    m.1 <- t(rep(1, p[1])/p[1])
    m.2 <- t(rep(1, p[2])/p[2])
    m.3 <- t(rep(1, p[3])/p[3])
    
    H.1 <- H[[1]]
    H.2 <- H[[2]]
    H.3 <- H[[3]]
    
    ahat <- atrans(Y, list(H.1, m.2, m.3))[, , 1]
    bhat <- atrans(Y, list(m.1, H.2, m.3))[, , 1]
    dhat <- atrans(Y, list(m.1, m.2, H.3))[1, , ]
    
    ehat <- atrans(Y, list(H.1, H.2, m.3))[, , 1]
    fhat <- atrans(Y, list(H.1, m.2, H.3))[, 1, ]
    ghat <- atrans(Y, list(m.1, H.2, H.3))[1, , ]
    
    s2e <- (mean(ehat^2) - bar.r.2/(p[3] - 1))*(p[1]*p[2])/((p[1] - 1)*(p[2] - 1))
    s2f <- (mean(fhat^2) - bar.r.2/(p[2] - 1))*(p[1]*p[3])/((p[1] - 1)*(p[3] - 1))
    s2g <- (mean(ghat^2) - bar.r.2/(p[1] - 1))*(p[2]*p[3])/((p[2] - 1)*(p[3] - 1))
    s2a <- (p[1]/(p[1] - 1))*(mean(ahat^2) - mean(ehat^2)/(p[2] - 1) - mean(fhat^2)/(p[3] - 1) + bar.r.2/((p[2] - 1)*(p[3] - 1)))
    s2b <- (p[2]/(p[2] - 1))*(mean(bhat^2) - mean(ehat^2)/(p[1] - 1) - mean(ghat^2)/(p[3] - 1) + bar.r.2/((p[1] - 1)*(p[3] - 1)))
    s2d <- (p[3]/(p[3] - 1))*(mean(dhat^2) - mean(fhat^2)/(p[1] - 1) - mean(ghat^2)/(p[2] - 1) + bar.r.2/((p[1] - 1)*(p[2] - 1)))
    
    out <- c(out, 
             "lambda.a" = ifelse(s2a > 0, sqrt(2/s2a), Inf),
             "lambda.b" = ifelse(s2b > 0, sqrt(2/s2b), Inf),
             "lambda.d" = ifelse(s2d > 0, sqrt(2/s2d), Inf),
             "lambda.e" = ifelse(s2e > 0, sqrt(2/s2e), Inf),
             "lambda.f" = ifelse(s2f > 0, sqrt(2/s2f), Inf),
             "lambda.g" = ifelse(s2g > 0, sqrt(2/s2g), Inf))
    
  } else if (lowpen) {
    cat("Penalties for lower-order parameters for tensors with more than three modes not supported by this function.\n")
  } 
  
  return(out)
  
}