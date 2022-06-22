#' @importFrom stats pnorm 
#' @export
.toy.rtdose <- function (CT.vol, PTV.radius, PTV.G, beam.nb) {
  
  sim.r <- PTV.radius + 5
  D.tot <- CT.vol
  D.tot$vol3D.data[] <- 0
  D <- CT.vol
  D$vol3D.data[] <- 0
  N <- CT.vol
  N$vol3D.data[] <- 0
  DSA <- 600
  theta <- rev (seq (360, 0, length.out = beam.nb + 1)[-1] * pi / 180)
  
  cone.b.a <- sim.r / DSA #angle du cone (= 1cm @ 1m)
  
  # distances min et max depuis le point de foc
  DMAX <- DSA + 256 * sqrt (2)
  d.theta <- min (CT.vol$dxyz) / DMAX / 1.5 # angle de simulation pour 1mm en distal
  uv <- as.matrix (expand.grid (seq (-cone.b.a, cone.b.a, by=d.theta),
                                seq (-cone.b.a, cone.b.a, by=d.theta)))
  r <- sqrt(uv[, 1]^2 + uv[, 2]^2) * DSA
  transverse.int <- 1 - pnorm (r, (PTV.radius + sim.r)/2, 1)
  
  uvt <- cbind (uv[r <= cone.b.a * DSA, ], transverse.int[r <= cone.b.a * DSA])
  colnames (uvt) <- c("u", "v", "t")
  
  
  for (angle in theta) {
    
    # point de foc
    foc <- PTV.G + DSA * c(sin (angle), -cos (angle), 0)
    
    PTV2foc <- foc - PTV.G
    PTV2foc <- PTV2foc / sqrt (sum (PTV2foc^2))
    
    n.b <- -PTV2foc # base du faisceau
    v.b <- c(n.b[2], n.b[3], n.b[1])
    u.b <- vector.product (n.b, v.b)
    u.b <- u.b / sqrt (sum (u.b^2))
    v.b <- vector.product (n.b, u.b)
    
    d.mM <- Re (polyroot (c (sum (foc^2) - 2 * 128^2, 2 * sum (foc * n.b), sum (n.b^2))))
    d.min <- d.mM[1]
    d.max <- d.mM[2]
    
    b.dirs <- t (apply (uvt, 1, function (uvt) {
      n.beamlet <- n.b + uvt[1] * u.b + uvt[2] * v.b
      c (n.beamlet / sqrt (sum (n.beamlet^2)), uvt)
    }))
    
    beamlets <- apply (b.dirs, 1, function (bl) {
      ds <- .95 * min (CT.vol$dxyz)
      L <- get.line (CT.vol, foc, bl[1:3], seq (d.min, d.max, by=ds))
      idx <- range (which (L$value > -900))
      L <- L[idx[1]:idx[2], ]
      L$mu <- L$value / 1000 + 1
      L$s_ <- cumsum (L$mu * ds)
      D <- pnorm(L$s_, 2, 5) * exp (-L$s_ / 150)
      L$D <- D / max (D) * bl[6]
      L$PTV <- (L$x - PTV.G[1])^2 + (L$y - PTV.G[2])^2 +
        (L$z - PTV.G[3])^2 <= PTV.radius^2
      return (L)
    })
    
    D$vol3D.data[] <- 0
    N$vol3D.data[] <- 0
    for (bl in beamlets) {
      ijk <- round (get.ijk.from.xyz (cbind (bl$x, bl$y, bl$z), CT.vol)) + 1
      D$vol3D.data[ijk] <- D$vol3D.data[ijk] + bl$D
      N$vol3D.data[ijk] <- N$vol3D.data[ijk] + 1
    }
    avg <- which (N$vol3D.data != 0)
    D$vol3D.data[avg] <- D$vol3D.data[avg] / N$vol3D.data[avg]
    D.tot$vol3D.data[avg] <- D.tot$vol3D.data[avg] + D$vol3D.data[avg]
    
  }
  
  D.tot$vol3D.data <- D.tot$vol3D.data / max (D.tot$vol3D.data) * 52
  D.tot$min.pixel <- 0
  D.tot$max.pixel <- 52
  
  D.tot$modality <- "rtdose"
  D.tot$description <- "IMRT|PTV52"
  
  return (D.tot)
}