#' Function for soft thresholding
#' @export
softt<-function(x,lambda){ sign(x)*pmax(abs(x)-lambda/2,0) }


#' Function for estimating unknown mean parameters
#'
#' \code{estM} uses block coordinate descent to obtain a
#' LANOVA penalized estimate of a mean matrix or three-way tensor.
#'
#' @param \code{Y} a matrix or three-way tensor of observed data
#' @param \code{ett} a list returned by the \code{estTune} function that contains the tuning parameter estimates to be used
#' @return A list of unknown mean parameter estimates.
#' @export
estM <-function(Y, ett, tol=1e-7) {

  if (length(dim(Y)) == 2) {

    n<-nrow(Y) ; p<-ncol(Y)

    sigma.2.z <- ett$sigma.2.z
    lambda.c <- ifelse(!is.infinite(ett$lambda.c), ett$lambda.c, 9999)
    lambda.a <- lambda.b <- 0

    if ("lambda.a" %in% names(ett)) {
      lambda.a <- ifelse(!is.infinite(ett$lambda.a),
                         ett$lambda.a, 9999)
      lambda.b <- ifelse(!is.infinite(ett$lambda.b),
                         ett$lambda.b, 9999)
    }

    mu.f <- mean(Y)
    a.f <- apply(Y, 1, mean) - mu.f
    b.f <- apply(Y, 2, mean) - mu.f
    C.f <- Y - (mu.f + outer(a.f, b.f, "+"))

    delta <- dev.old <-Inf; iter<-0
    while(delta > tol) {

    ## update a
    yr<-apply( sweep(Y-C.f,2,b.f,"-"),1,mean) - mu.f
    a.f<-softt(yr, lambda.a*sigma.2.z*2/p)

    ## update b
    yc<-apply( sweep(Y-C.f,1,a.f,"-"),2,mean) - mu.f
    b.f<-softt(yc, lambda.b*sigma.2.z*2/n)

    ## update C
    C.f<-softt( Y- mu.f - outer(a.f,b.f,"+"), lambda.c*sigma.2.z*2)

    ## update mu
    mu.f <-mean( Y - outer(a.f,b.f,"+") - C.f )

    dev.new<-sum(( Y-(mu.f + outer(a.f,b.f,"+") + C.f))^2)/(2*sigma.2.z) +
      lambda.a*sum(abs(a.f)) + lambda.b*sum(abs(b.f)) + lambda.c*sum(abs(C.f))
    delta <- dev.old - dev.new
    dev.old <- dev.new
    iter <- iter + 1
    if(iter%%100 == 0) { cat(iter,dev.old, delta,"\n") }
  }
  out <- list(mu.f=mu.f,a.f=a.f,b.f=b.f,C.f=C.f,
              M.f = mu.f + outer(a.f,b.f,"+") + C.f # , dev=2*dev.old
              )
  }

  if (length(dim(Y)) == 3) {

    n<-dim(Y)[1] ; p<-dim(Y)[2] ; q <- dim(Y)[3]

    m.1 <- t(rep(1, n)/n)
    m.2 <- t(rep(1, p)/p)
    m.3 <- t(rep(1, q)/q)

    H.1 <- diag(n) - rep(1, n)%*%t(rep(1, n))/n
    H.2 <- diag(p) - rep(1, p)%*%t(rep(1, p))/p
    H.3 <- diag(q) - rep(1, q)%*%t(rep(1, q))/q

    sigma.2.z <- ett$sigma.2.z
    lambda.c <- ifelse(!is.infinite(ett$lambda.c), ett$lambda.c, 9999)
    lambda.a <- lambda.b <- lambda.d <- lambda.e <- lambda.f <- lambda.g <- 0

    if ("lambda.a" %in% names(ett)) {
      lambda.a <- ifelse(!is.infinite(ett$lambda.a), ett$lambda.a, 9999)
      lambda.b <- ifelse(!is.infinite(ett$lambda.b), ett$lambda.b, 9999)
      lambda.d <- ifelse(!is.infinite(ett$lambda.d), ett$lambda.d, 9999)
      lambda.e <- ifelse(!is.infinite(ett$lambda.e), ett$lambda.e, 9999)
      lambda.f <- ifelse(!is.infinite(ett$lambda.f), ett$lambda.f, 9999)
      lambda.g <- ifelse(!is.infinite(ett$lambda.g), ett$lambda.g, 9999)
    }


    mu.f <- atrans(Y, list(m.1, m.2, m.3))[, , 1]
    a.f <- atrans(Y, list(H.1, m.2, m.3))[, , 1]
    b.f <- atrans(Y, list(m.1, H.2, m.3))[, , 1]
    d.f <- atrans(Y, list(m.1, m.2, H.3))[1, , ]
    E.f <- atrans(Y, list(H.1, H.2, m.3))[, , 1]
    F.f <- atrans(Y, list(H.1, m.2, H.3))[, 1, ]
    G.f <- atrans(Y, list(m.1, H.2, H.3))[1, , ]

    C.f <- atrans(Y, list(H.1, H.2, H.3))

    M.f <- sweep(C.f, c(2, 3), G.f, "+")
    M.f <- sweep(M.f, c(1, 3), F.f, "+")
    M.f <- sweep(M.f, c(1, 2), E.f, "+")
    M.f <- sweep(M.f, c(3), d.f, "+")
    M.f <- sweep(M.f, c(2), b.f, "+")
    M.f <- sweep(M.f, c(1), a.f, "+")
    M.f <- M.f + mu.f

    delta <- dev.old <-Inf; iter<-0
    while(delta > tol) {

      ## update a
      yr <- atrans(Y - sweep(M.f, 1, a.f, "-"), list(H.1, m.2, m.3))[, , 1]
      a<-softt(yr, lambda.a*sigma.2.z*2/(p*q))
      M.f <- sweep(M.f, 1, a - a.f, "+")
      a.f <- a

      ## update b
      yr <- atrans(Y - sweep(M.f, 2, b.f, "-"), list(m.1, H.2, m.3))[, , 1]
      b<-softt(yr, lambda.b*sigma.2.z*2/(n*q))
      M.f <- sweep(M.f, 2, b - b.f, "+")
      b.f <- b

      ## update d
      yr <- atrans(Y - sweep(M.f, 3, d.f, "-"), list(m.1, m.2, H.3))[1, , ]
      d<-softt(yr, lambda.d*sigma.2.z*2/(n*p))
      M.f <- sweep(M.f, 3, d - d.f, "+")
      d.f <- d

      ## update e
      yr <- atrans(Y - sweep(M.f, c(1, 2), E.f, "-"), list(H.1, H.2, m.3))[, , 1]
      E<-softt(yr, lambda.e*sigma.2.z*2/(q))
      M.f <- sweep(M.f, c(1, 2), E - E.f, "+")
      E.f <- E

      ## update f
      yr <- atrans(Y - sweep(M.f, c(1, 3), F.f, "-"), list(H.1, m.2, H.3))[, 1, ]
      FF<-softt(yr, lambda.f*sigma.2.z*2/(p))
      M.f <- sweep(M.f, c(1, 3),  FF - F.f, "+")
      F.f <- FF

      ## update g
      yr <- atrans(Y - sweep(M.f, c(2, 3), G.f, "-"), list(m.1, H.2, H.3))[1, , ]
      G<-softt(yr, lambda.g*sigma.2.z*2/(n))
      M.f <- sweep(M.f, c(2, 3), G - G.f, "+")
      G.f <- G

      ## update C
      C<-softt( Y - (M.f - C.f), lambda.c*sigma.2.z*2)
      M.f <- (M.f - C.f) + C
      C.f <- C

      ## update mu
      mu<-mean( Y - (M.f -mu.f))
      M.f <- (M.f - mu.f) + mu
      mu.f <- mu

      dev.new<-sum(( Y-M.f)^2)/(2*sigma.2.z) +
        lambda.a*sum(abs(a.f)) + lambda.b*sum(abs(b.f)) + lambda.d*sum(abs(d.f)) +
        lambda.e*sum(abs(E.f)) + lambda.f*sum(abs(F.f)) + lambda.g*sum(abs(G.f)) +
        lambda.c*sum(abs(C))
      delta <- dev.old - dev.new
      dev.old <- dev.new
      iter <- iter + 1
      if(iter%%1 == 0) { cat(iter,dev.old, delta,"\n") }
    }
    out <- list(mu.f=mu.f,a.f=a.f,b.f=b.f,d.f = d.f, E.f = E.f,
                F.f = F.f,
                G.f = G.f,
                C.f=C.f,
                M.f = M.f #,  dev=2*dev.old
                )
  }

  return(out)
}


