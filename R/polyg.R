
#' @importFrom stats approx
polyg.UniformCurvCoord <- function (pol, N = NULL, tol = 0.1) {
  
  .get.err <- function (pol, pol.ref) {
    if (length(pol)==0) return(pol)
    M <- as.matrix(pol)
    apply (pol, 1, function (Mi) {
      sqrt (min ((pol.ref[,1] - Mi[1])^2 + (pol.ref[,2] - Mi[2])^2))
    })
  }
  
  if (length(pol)==0) return(pol)
  if (!(ncol(pol)==2|ncol(pol)==3)) stop("pol must have 2 or 3 columns")
  
  m.f <- is.matrix(pol)

  
  if (is.null (N)) {
    try.once <- FALSE
    N <- 8
  } else {
    try.once <- TRUE
  }
  
  M <- round (as.matrix(pol[,1:2]),6)
  closed <- all(M[1, ] == M[nrow(M), ])
  if (!closed) M <- rbind (M, M[1, ])
  
  l <- sapply (1:(nrow(M) - 1), function (i) sqrt ((M[i, 1] - M[i+1, 1])^2 + (M[i, 2] - M[i+1, 2])^2))
  L.tot <- c (0, cumsum (l))
  L.1 <- L.tot / max (L.tot)
  
  repeat {
    suppressWarnings (x.i <- approx (L.1, M[, 1], rev (seq (1, 0, length.out = N + 1)[-1]))$y)
    suppressWarnings (y.i <- approx (L.1, M[, 2], rev (seq (1, 0, length.out = N + 1)[-1]))$y)
    if (try.once) break
    tol_ <- max (.get.err (M, pol.ref = data.frame(x=x.i, y=y.i)))
    if (tol_ < tol) break
    N <- N * 2
  }
  pol_ <-data.frame(x=x.i, y=y.i)
  if (ncol(pol)==3) pol_$z  <- pol[1,3] 
  if (closed) pol_ <- pol_[c(1:nrow(pol_), 1), ]
  if(m.f) return (as.matrix(pol_))
  return(pol_)
}


################################################################################
#' @importFrom stats fft
.polyg.filtering <- function (pol, bw = 0.1) {
  if (length(pol)==0) return(pol)
  if (!(ncol(pol)==2|ncol(pol)==3)) stop("pol must have 2 or 3 columns")
  
  m.f <- is.matrix(pol)
  le <- nrow(pol)
  X <- fft (pol[,1]) / le
  Y <- fft (pol[,2]) / le
  
  N <- le * bw / 2
  
  LE <- length(X)
  X[(N + 1):(LE + 1 - N)] <- 0
  Y[(N + 1):(LE + 1 - N)] <- 0
  
  x.f <- Re (fft (X, inverse=TRUE))
  y.f <- Re (fft (Y, inverse=TRUE))
  
  pol_ <-data.frame(x=x.f, y=y.f)
  if (ncol(pol)==3) pol_$z  <- pol[1,3] 
  if(m.f) return(as.matrix(pol_))
  return(pol_)
}

################################################################################

#' @importFrom stats approx
polyg.simplify <- function (pol, tol = 0.5) {
  le <- nrow (pol)
  if (le ==0) return(pol)
  if (!(ncol(pol)==2|ncol(pol)==3)) stop("pol must have 2 or 3 columns")
  # pol <- rbind(pol, pol[1,])
  
  splits <- c(1, nrow(pol))
  err.splits <- c()
  
  repeat {
    x.s <- approx (splits, pol[splits,1], 1:le)$y
    y.s <- approx (splits, pol[splits,2], 1:le)$y
    
    err <- (pol[,1] - x.s)^2 + (pol[,2]  - y.s)^2
    err.mean <- sqrt (mean (err, na.rm=TRUE))
    new.split <- which.max (err)
    splits <- sort (c(splits, new.split))
    # err.splits <- c(err.splits, sqrt (err[new.split]))
    err.splits <- c(err.splits, err.mean)
    
    if (err.mean < tol) break
  }
  return (pol[splits,])
}

################################################################################

polyg.sort <- function(pol, clockwise = TRUE){
  mf <- is.matrix(pol)
  df <- is.data.frame(pol)
  if (!(mf|df)) stop("must be a matrix or dataframe with 2 or 3 column")
  ncol_pol <- ncol(pol)
  if (ncol_pol>3 & ncol_pol<2) stop("must be a matrix or dataframe with 2 or 3 column")
  pt <- matrix(.polygcleanC(as.vector(t(as.matrix(pol))),ncol=ncol_pol, 
                                      sort = TRUE,clockwise=clockwise), byrow = TRUE,
               ncol=ncol_pol,
               dimnames=list(NULL,c("x","y","z")[1:ncol_pol]))
  
  if (mf) return(pt)
  return(as.data.frame(pt))
#   le <- nrow(pol)
#   if (le ==0) return(pol)
#   if (!(ncol(pol)==2|ncol(pol)==3)) stop("pol must have 2 or 3 columns")
# 
#   if (!all(round(pol[1,],6)==round(pol[le,],6))) return (pol)
#   
#   m.f <- is.matrix(pol)
#   pol <- as.data.frame(pol)
#   first <- order(pol[1:(le-1),1],pol[1:(le-1),2])[1]
#   pol <- rbind(pol[first:(le-1),], pol[-(first:le),], pol[first[1],])
#   
#   if (polyg.is.clockwise(pol)!= clockwise) pol <- pol[le:1, ]
#   rownames(pol) <- NULL
#   if (m.f) return(as.matrix(pol))
#   return(pol)
}

################################################################################

polyg.is.clockwise <- function(pol) {
  le <- nrow(pol)
  if (le==0) return(FALSE)
  pol <- round(pol,6)
  pol <- pol[c(TRUE, !((diff(pol[,1])==0) & (diff(pol[,2])==0))),]
  first <- order(pol[1:(le-1),1],pol[1:(le-1),2])[1]
  pol <- rbind(pol[first:(le-1),], pol[-(first:le),])
  if(ncol(pol)==2) pol <- cbind(pol,rep(0,le-1))
  r <- vector.product(pol[1,]-pol[le-1,],pol[2,]-pol[1,])
  while (r[3]==0) {
    pol <- pol[-1,]
    le <- length(pol)
    if (length(pol)==0) return(FALSE)
    first <- order(pol[,1],pol[,2])[1]
    pol <- rbind(pol[first:le,], pol[-(first:le),])
    r <- vector.product(pol[1,]-pol[le-1,],pol[2,]-pol[1,])
  }
  return(r[3]<0)
}
################################################################################

polyg.area <- function(pol) {
  
  pol_le <- nrow(pol)
  pol <- round (pol,6)
  if (!all(pol[1,]==pol[pol_le,])) return(0.0)
  idx <- 2:pol_le
  return(abs(sum (pol[idx-1,1] * pol[idx,2] - pol[idx,1] * pol[idx-1,2]) / 2))
}

################################################################################

.polyg.add.common.pt <- function(pol1, pol2,...) {
  eps=1e-9
  m.f1 <- is.matrix(pol1)
  m.f2 <- is.matrix(pol2) 
  pol1 <- as.data.frame(pol1)
  pol2 <- as.data.frame(pol2)
  
  pol1$idx <- 1:nrow(pol1)
  pol2$idx <- 1:nrow(pol2)
  u1_x <- diff(pol1[,1])
  u1_y <- diff(pol1[,2])
  u2_x <- diff(pol2[,1])
  u2_y <- diff(pol2[,2])
  d1 <- sqrt(u1_x^2 + u1_y^2)
  d2 <- sqrt(u2_x^2 + u2_y^2)
  u1_x <- u1_x/d1; u1_y <- u1_y/d1
  u2_x <- u2_x/d2; u2_y <- u2_y/d2
  
  L <- .addcommonptC (pt1_x=pol1[,1], pt1_y=pol1[,2], pt2_x=pol2[,1], pt2_y=pol2[,2],
                      u1_x, u1_y,u2_x, u2_y, d1, d2, eps = eps)
  add.pol <- matrix(NA, nrow =length(L$k1), ncol =5, 
                    dimnames = list(NULL,c("x","y","z","idx1","idx2")))
  
  pol1 <- pol1[-nrow(pol1),]
  pol2 <- pol2[-nrow(pol2),]
  if (!is.null(L$k1)){
    add.pol[ ,1:2] <- cbind(L$k1 * u1_x[L$pt1.index] + pol1[L$pt1.index,1],
                            L$k1 * u1_y[L$pt1.index] + pol1[L$pt1.index,2])
    add.pol[ ,3] <- pol1[1,3]
    delta1 <- L$k1/d1[L$pt1.index]
    delta2 <- L$k2/d2[L$pt2.index]
    add.pol[ ,4] <- L$pt1.index + delta1
    add.pol[ ,5] <- L$pt2.index + delta2
    # add.pol <- add.pol[order(add.pol[,4]), ]
    f <- !duplicated(round(add.pol[,4],-log10(eps)))
    
    add.pol <- add.pol [f , ]
    if (sum(f)==1) add.pol <- matrix(add.pol, ncol=5,dimnames = list(NULL,c("x","y","z","idx1","idx2")))
    rdelta1 <- round(delta1[f],-log10(eps))
    suppress.f <-  rdelta1 == 0 |  rdelta1 == 1
    if (any(suppress.f)) {
      suppress.idx <- unique(L$pt1.index[f][suppress.f]+ rdelta1[suppress.f])
      pol1 <- pol1[-suppress.idx,]
    }
    
    rdelta2 <- round(delta2[f],-log10(eps))
    suppress.f <-  rdelta2 == 0 |  rdelta2 == 1
    if (any(suppress.f)) {
      suppress.idx <- unique(L$pt2.index[f][suppress.f] + rdelta2[suppress.f])
      pol2 <- pol2[-suppress.idx,]
    }
    
  }
  
  
  pol1 <- as.data.frame(rbind(as.matrix(pol1), add.pol[,c(1:3,4)]))
  pol1 <- pol1 [!duplicated (pol1$idx),]
  pol1 <- pol1[order(pol1$idx),]
  pol1 <- pol1[c(1:nrow(pol1),1),]
  row.names(pol1) <- NULL
  
  add.pol[,3] <-  pol2[1,3]
  pol2 <-  as.data.frame(rbind(as.matrix(pol2), add.pol[,c(1:3,5)]))
  pol2 <- pol2 [!duplicated (pol2$idx),]
  pol2 <- pol2[order(pol2$idx),]
  pol2 <- pol2[c(1:nrow(pol2),1),]
  row.names(pol2) <- NULL
  
  pol1$idx <- NULL
  pol2$idx <- NULL
  if (m.f1) pol1 <- as.matrix(pol1)
  if (m.f2) pol2 <- as.matrix(pol2)
  return(list(pol1 = pol1, pol2 = pol2))
}

