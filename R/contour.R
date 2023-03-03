
.contour.curv.coord <- function (L, N = NULL, tol = 0.1) {
  
  if (is.null (N)) {
    try.once <- FALSE
    N <- 8
  } else {
    try.once <- TRUE
  }
  
  repeat {
    M <- cbind (L$x, L$y)
    M <- rbind (M, M[1, ])
    l <- sapply (1:(nrow(M) - 1), function (i) sqrt ((M[i, 1] - M[i+1, 1])^2 + (M[i, 2] - M[i+1, 2])^2))
    L.tot <- c (0, cumsum (l))
    L.1 <- L.tot / max (L.tot)
    
    suppressWarnings (x.i <- approx (L.1, M[, 1], rev (seq (1, 0, length.out = N + 1)[-1]))$y)
    suppressWarnings (y.i <- approx (L.1, M[, 2], rev (seq (1, 0, length.out = N + 1)[-1]))$y)
    
    if (try.once) break
    
    tol_ <- max (.contour.distance (L, list (x=x.i, y=y.i)))
    
    if (tol_ < tol) break
    N <- N * 2
  }
  
  return (list (x=x.i, y=y.i, x0=L$x, y0=L$y, s=L.1[-length (L.1)], l=max (L.tot)))
}
################################################################################

.contour.distance <- function (ref, L) {
  M.ref <- cbind (ref$x, ref$y)
  apply (M.ref, 1, function (Mi) {
    sqrt (min ((L$x - Mi[1])^2 + (L$y - Mi[2])^2))
  })
}

################################################################################

.contour.filtering <- function (L, bw = 0.1) {
  X <- fft (L$x) / length (L$x)
  Y <- fft (L$y) / length (L$y)
  
  N <- length (L$x) * bw / 2
  
  X[(N + 1):(length(X) + 1 - N)] <- 0
  Y[(N + 1):(length(X) + 1 - N)] <- 0
  
  x.f <- Re (fft (X, inverse=TRUE))
  y.f <- Re (fft (Y, inverse=TRUE))
  
  L$x <- x.f
  L$y <- y.f
  return (L)
}

################################################################################

.contour.simplify <- function (L.i, tol = 0.5) {
  L <- L.i
  L$x <- c(L$x, L$x[1])
  L$y <- c(L$y, L$y[1])
  splits <- c(1, length (L.i$x))
  
  err.splits <- c()
  
  repeat {
    x.s <- approx (splits, L$x[splits], 1:length (L$x))$y
    y.s <- approx (splits, L$y[splits], 1:length (L$y))$y
    
    err <- (L$x - x.s)^2 + (L$y - y.s)^2
    err.mean <- sqrt (mean (err, na.rm=TRUE))
    new.split <- which.max (err)
    splits <- sort (c(splits, new.split))
    # err.splits <- c(err.splits, sqrt (err[new.split]))
    err.splits <- c(err.splits, err.mean)
    
    if (err.mean < tol) break
  }
  return (list (x=L$x[splits], y=L$y[splits], err.splits=err.splits))
}