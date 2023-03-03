#' Load a toy patient for test
#' @description The \code{toy.load.patient} creates a dummy "patient" class object.
#' It is used for the test.
#' @param modality String vector, whose elements are chosen among the modalities 
#' "ct", "mr", "rtstruct" and "rtdose".
#' @param roi.name String vector, whose elements are chosen among the regions of 
#' interest (RoI) "eye", "optical nerve", "brain", "labyrinth processing unit",
#' "energy unit", "gizzard", "ghost container" and "exhaust valve". Note that the
#' RoI "couch", "patient" and "ptv" are still present.
#' @param dxyz Vector of length 3, representing the x, y, z steps in mm, between 
#' ct, mr and rtdose voxels.
#' @param beam.nb Positive integer. Number of radiotherapy beams in rtdose modality.
#' @return Returns an toy object of "patient" class, containing the modalities 
#' defined in \code{modality}. See \link[espadon]{espadon.class} for class definitions.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' step <- 5
#' pat <- toy.load.patient (dxyz = rep (step, 3), beam.nb = 2)
#' str (pat, max.level = 2)

#' @importFrom stats runif spline
#' @export
toy.load.patient <- function ( 
    modality = c("ct", "mr", "rtdose", "rtstruct"),
    roi.name = c("eye", "optical nerve", "brain","labyrinth processing unit",
                "energy unit", "gizzard","ghost container", "exhaust valve" ), 
    dxyz = c(1,1,1), beam.nb = 7){
  m <- match(modality,c("ct", "mr", "rtdose", "rtstruct"))
  
  if (all(is.na(m))) {
    warning ("modality must contain at least one known object.")
    return (NULL)
  }
  
  modality <- modality[!is.na(m)]
  ref1.do.nb <- rep(0,3)
  f <-!is.na(match( c("ct", "rtdose", "rtstruct"), modality))
  ref1.do.nb[f] <- 1:sum(f)
  
  
  S <- NULL
  rtstruct.f <- "rtstruct" %in% modality
  ct.f <- "ct" %in% modality
  mr.f <- "mr" %in% modality
  rtdose.f<- "rtdose" %in% modality
  
  #initialisation
  
  MR.date <- "20220523"
  CT.date <- "20220609"
  rt.date <- "20220610"
  
  mr.offset <- c(50,100,-70)
  pat.name <- "PMan"
  pat.pseudo <- "PM"
  bd <- "19800522"
  sex <- c("M","O")
  if (!rtstruct.f) sex <- sex[1]
  
  PM.radius <- 120
  ptv.radius <- 10
  ptv.G <- c (-30.2,50.5,-30.1)
  
  
  
 

  
  pat <- list()
  pat$patient <- data.frame(PIN = pat.pseudo, name =pat.name,  birth.date=bd, sex=sex)
  pat$pat.pseudo <- pat$patient[1,1]
  pat$description <- data.frame(PIN=character(0), modality=character(0),
                                obj=character(0), ref.pseudo = character(0),
                                nb.of.suboject = character(0), description=character(0),
                                nb = character(0),max=character(0),object.alias = character(0))
  pat$description.by.reg <- list()
  ########
  pat$T.MAT = list ()
  pat$T.MAT$ref.info <- list ()
  pat$T.MAT$ref.info <- data.frame(ref.pseudo=character(0), ref=character(0))
  
  if (any (!is.na(match(modality, c("ct", "rtdose")))))
    pat$T.MAT$ref.info <- rbind(pat$T.MAT$ref.info, 
                                data.frame(ref.pseudo = "ref1", 
                                           ref = "2.16.840.1.114357.58485835081.6250.1551498438.0.3"))
  if (any (!is.na(match(modality, c("rtstruct")))))
    pat$T.MAT$ref.info <- rbind(pat$T.MAT$ref.info, 
                                data.frame(ref.pseudo = "ref1", 
                                           ref = "2.16.840.1.114357.58485835081.6250.1551498438.0.3"))
  
  
  if (any (!is.na(match(modality, c("mr")))))
    pat$T.MAT$ref.info <- rbind(pat$T.MAT$ref.info, 
                                data.frame(ref.pseudo = "ref2", 
                                           ref = "2.16.840.1.114357.58094835081.62350.15267498438.0.1"))
  pat$T.MAT$ref.info <- unique(pat$T.MAT$ref.info[,1:2])
  pat$T.MAT$reg.info <- list()
  pat$T.MAT$reg.info$patient <- data.frame(patient = pat.pseudo, patient.name = pat.name, patient.bd=bd, patient.sex=sex)
  pat$T.MAT$reg.info$file <- data.frame(t = "ref1<-ref2", path="local")
  
  pat$T.MAT$matrix.description <- data.frame(t=character(0), src=character(0),
                                             dest=character(0),type=character(0))
  pat$T.MAT$matrix.list <- list()
  ref.pseudo <- unique (pat$T.MAT$ref.info$ref.pseudo)
  if (length(ref.pseudo) == 2) {
    pat$T.MAT$matrix.description <- data.frame (t=c("ref1<-ref1","ref2<-ref1","ref1<-ref2", "ref2<-ref2"), 
                                                src = c("ref1","ref1", "ref2", "ref2"),
                                                dest = c("ref1","ref2", "ref1", "ref2"),
                                                type = rep("RIGID", 4))
    pat$T.MAT$matrix.list[["ref1<-ref1"]] <- diag(4)
    pat$T.MAT$matrix.list[["ref2<-ref1"]] <- cbind(diag(1,4,3), c(mr.offset,1))
    pat$T.MAT$matrix.list[["ref1<-ref2"]] <- cbind(diag(1,4,3), c(-mr.offset,1))
    pat$T.MAT$matrix.list[["ref2<-ref2"]] <- diag(4)
  } else if (ref.pseudo == "ref1"){ 
    pat$T.MAT$matrix.description <- data.frame (t=c("ref1<-ref1"), 
                                                src = c("ref1"),
                                                dest = c("ref1"),
                                                type = "RIGID")
    pat$T.MAT$matrix.list[["ref1<-ref1"]] <- diag(4)
  } else if (ref.pseudo == "ref2") {
    pat$T.MAT$matrix.description <- data.frame (t=c("ref2<-ref2"), 
                                                src = c("ref2"),
                                                dest = c("ref2"),
                                                type = "RIGID")
    pat$T.MAT$matrix.list[["ref2<-ref2"]] <- diag(4)
  }
  class (pat$T.MAT ) <- "t.mat"
  ###########
  
  fat <- function (vol, ijk, radius, value) {
    dijk <- as.matrix (expand.grid (-radius:radius, -radius:radius, -radius:radius))
    keep <- apply (dijk, 1, function (v) sum(v^2) <= radius^2)
    dijk <- dijk[keep, ]
    dummy <- apply (dijk, 1, function (d) vol$vol3D.data[t (t (ijk) + d)] <<- value)
    return (vol)
  }
  
  n.ijk  <- ceiling(c(128, 128, 128) / dxyz)*2
  pt000 <- -(n.ijk/2)*dxyz + c(1,1,1)
  b <- vol.create (n.ijk = n.ijk, dxyz = dxyz, pt000=pt000, modality = "binary",
                   default.value = FALSE, ref.pseudo = "ref1", 
                   frame.of.reference = "2.16.840.1.114357.58485835081.6250.1551498438.0.3",
                   number = 1)
  b$max.pixel <- TRUE
  
  seqx <- seq(pt000[1], (n.ijk[1]-1)*dxyz[1] + pt000[1],dxyz[1])
  seqy <- seq(pt000[2], (n.ijk[2]-1)*dxyz[2] + pt000[2],dxyz[2])
  seqz <- seq(pt000[3], (n.ijk[3]-1)*dxyz[3] + pt000[3],dxyz[3])
  xyz<-expand.grid(seqx, seqy,seqz)
  
  
  
  
  
  
  xyz.L.eye.G <- 0.95 * PM.radius * c (sin (40/180*pi) * sin (45/180*pi),
                                       -cos (40/180*pi) * sin (45/180*pi),
                                       cos (45/180*pi))
  xyz.R.eye.G <- 0.95 * PM.radius * c (sin (-40/180*pi) * sin (45/180*pi),
                                       -cos (-40/180*pi) * sin (45/180*pi),
                                       cos (45/180*pi))
  
  
  # Pacman external contour
  
  PM1 <- xyz[, 1]^2 + xyz[, 2]^2 + xyz[, 3]^2 <= PM.radius^2 &
    (xyz[, 2] - 2 * xyz[, 3] > 0 | -xyz[, 2] - 2 * xyz[, 3] < 0)
  b.PM.ext <- b
  b.PM.ext$vol3D.data[PM1] <- TRUE
  # display.plane (b.PM.ext, view.type = "sagi", main = "PM.ext",view.coord = xyz.L.eye.G[3])
  if (rtstruct.f) {
    S <- struct.merge(S,struct.from.bin(b.PM.ext, roi.name = "patient" , roi.nb = 2,
                                        roi.color = "#ffffff", roi.type = "EXTERNAL",
                                        external.only = TRUE))
    # display.plane (struct = S, view.type = "sagi", main = "PM.ext")
  }
  
  # Pacman interior
  PM.shielding <- 6
  PM1 <- xyz[, 1]^2 + xyz[, 2]^2 + xyz[, 3]^2 <= (PM.radius - PM.shielding)^2 &
    (xyz[, 2] - 2 * xyz[, 3] > 2 * PM.shielding | -xyz[, 2] - 2 * xyz[, 3] < -2 * PM.shielding)
  b.PM.int <- b
  b.PM.int$vol3D.data[PM1] <- TRUE
  # display.plane (b.PM.int, view.coord = xyz.L.eye.G[3], main = "PM int")
  
  # b.bone <- bin.subtraction (b.PM.ext, b.PM.int)
  # display.plane (b.bone, view.coord = xyz.L.eye.G[3], main = "b.bone", view.type = "sagi")
  
  # couch

  PM1 <- xyz[, 2] >= PM.radius |
    (xyz[, 2] >= 0.3 * PM.radius & xyz[, 3] >= 0.4 * PM.radius & xyz[, 3] <= 0.45 * PM.radius) |
    (xyz[, 2] >= 0.3 * PM.radius & xyz[, 3] <= -0.4 * PM.radius & xyz[, 3] >= -0.45 * PM.radius)
  b.couch <- b
  b.couch$vol3D.data[PM1] <- TRUE
  b.couch <- bin.subtraction (b.couch, b.PM.ext)
  
  if (rtstruct.f) {
    S <- struct.merge(S,struct.from.bin(b.couch, roi.name = "couch" , roi.nb = 1,
                                        roi.color = "#808080", roi.type = "SUPPORT",
                                        external.only = TRUE))
  }
  # display.plane (b.couch, view.type = "sagi", main = "couch")
  
  # external mask
  b.ext.mask <- bin.inversion (bin.sum (b.PM.ext, b.couch))
  # display.plane (b.ext.mask, view.coord = xyz.L.eye.G[3], main = "ext mask")
   

  
  # eyes
  f.eye <- "eye" %in% roi.name
  b.eye <- NULL
  b.eyes_shield <- NULL
  if  (f.eye){
    PM.eyes.radius <- 20

    b.eye.L <- b
    b.eye.R <- b
    
    PM1 <- ((xyz[, 1] - xyz.L.eye.G[1]) / PM.eyes.radius)^2 +
      ((xyz[, 2] - xyz.L.eye.G[2]) / PM.eyes.radius)^2 +
      ((xyz[, 3] - xyz.L.eye.G[3]) / PM.eyes.radius)^2 <= 1
    b.eye.L$vol3D.data[PM1] <- TRUE
    
    PM1 <- ((xyz[, 1] - xyz.R.eye.G[1]) / PM.eyes.radius)^2 +
      ((xyz[, 2] - xyz.R.eye.G[2]) / PM.eyes.radius)^2 +
      ((xyz[, 3] - xyz.R.eye.G[3]) / PM.eyes.radius)^2 <= 1
    b.eye.R$vol3D.data[PM1] <- TRUE
    
    b.eye.L <- bin.intersection (b.eye.L,b.PM.ext)
    b.eye.R <- bin.intersection (b.eye.R,b.PM.ext)
    
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.eye.L, roi.name = "left eye" , roi.nb = 3,
                                          roi.color = "#ff0000", roi.type = "ORGAN", external.only = TRUE))
      S <- struct.merge(S,struct.from.bin(b.eye.R, roi.name = "right eye" , roi.nb = 4,
                                          roi.color = "#00ff00", roi.type = "ORGAN", external.only = TRUE))
    }
   
    # display.plane (bin.sum (b.eye.L, b.eye.R), view.coord = xyz.L.eye.G[3], main = "eye L|R")
   
    
    # eyes shielding
    PM.eyes.radius.s <- 26
    b.eye_shield.L <- b
    b.eye_shield.R <- b
    
    PM1 <- ((xyz[, 1] - xyz.L.eye.G[1]) / PM.eyes.radius.s)^2 +
      ((xyz[, 2] - xyz.L.eye.G[2]) / PM.eyes.radius.s)^2 +
      ((xyz[, 3] - xyz.L.eye.G[3]) / PM.eyes.radius.s)^2 <= 1
    b.eye_shield.L$vol3D.data[PM1] <- TRUE
    
    PM1 <- ((xyz[, 1] - xyz.R.eye.G[1]) / PM.eyes.radius.s)^2 +
      ((xyz[, 2] - xyz.R.eye.G[2]) / PM.eyes.radius.s)^2 +
      ((xyz[, 3] - xyz.R.eye.G[3]) / PM.eyes.radius.s)^2 <= 1
    b.eye_shield.R$vol3D.data[PM1] <- TRUE
    
    b.eyes_shield <- bin.intersection (bin.sum (b.eye_shield.L, b.eye_shield.R),
                                       b.PM.ext)
    b.eye <- bin.sum (b.eye.L, b.eye.R)
    
    rm(list = c("b.eye_shield.R","b.eye_shield.L","b.eye.L","b.eye.R"))
    # display.plane (b.eye, view.coord = xyz.L.eye.G[3], main = "eye")
    # display.plane (b.eyes_shield, view.coord = xyz.L.eye.G[3], main = "ext eye_shield")
   
    b.ext.mask <- bin.sum (b.ext.mask, b.eye)
    #display.plane (b.ext.mask, view.coord = xyz.L.eye.G[3], main = "ext mask")
    b.PM.int <- bin.subtraction (b.PM.int, b.eyes_shield)
    # display.plane (b.PM.int, view.coord = xyz.L.eye.G[3], main = "PM int")
  }
  
  f.brain <- "brain" %in% roi.name
  b.brain.ext <- NULL
  b.brain.int <- NULL
  if (f.brain){
    # brain exterior
    PM.brain.radius <- 30
    PM1 <- (xyz[, 1] / (2.2 * PM.brain.radius))^2 +
      ((xyz[, 2] - 0.4 * PM.radius) / PM.brain.radius)^2 +
      ((xyz[, 3] - 0.5 * PM.radius) / PM.brain.radius)^2 <= 1
    b.brain.ext <- b
    b.brain.ext$vol3D.data[PM1] <- TRUE
    
    # display.plane (b.brain.ext, view.type  = "sagi", main = "brain ext")
    
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.brain.ext, roi.name = "brain" , roi.nb = 5,
                                          roi.color = "#0080f0", roi.type = "ORGAN", external.only = TRUE))}
    
    # brain interior
    PM1 <- (xyz[, 1] / (2.2 * PM.brain.radius - 3))^2 +
      ((xyz[, 2] - 0.4 * PM.radius) / (PM.brain.radius - 3))^2 +
      ((xyz[, 3] - 0.5 * PM.radius) / (PM.brain.radius - 3))^2 <= 1
    b.brain.int <- b
    b.brain.int$vol3D.data[PM1] <- TRUE
    # display.plane (b.brain.int, view.type  = "sagi", main = "brain int")
  }
  
  # optical nerves
  f.on <- "optical nerve" %in% roi.name
  b.on <- NULL
  b.on.s <- NULL
  if (f.on){
    Sp.way <- function(cp){
      d <- 0
      for (i in 2:nrow (cp)) d <- c (d, sqrt (sum ((cp[i, ] - cp[i-1, ])^2)))
      d.tot <- cumsum (d)
      df <- as.data.frame (cbind (cp, d.tot))
      x <- spline (df$d.tot, df$V1, n=max (df$d.tot) * 10)$y
      y <- spline (df$d.tot, df$V2, n=max (df$d.tot) * 10)$y
      z <- spline (df$d.tot, df$V3, n=max (df$d.tot) * 10)$y
      # lines3d (cp)
      # points3d (x, y, z, col="red")
      xyz.on <- cbind (x, y, z)
      xyz.on
    }
    cp <- matrix (c(xyz.R.eye.G, xyz.R.eye.G + c(0, 20, 0),
                    c(0, 0.4 * PM.radius, 0.5 * PM.radius)), ncol=3, byrow=TRUE)
    ijk <- get.ijk.from.xyz(Sp.way(cp), b) + 1
    b.on.s.R <- b
    b.on.s.R <- fat (b.on.s.R, ijk, 3, TRUE)
    b.on.R <- b
    b.on.R <- fat (b.on.R, ijk, 1, TRUE)
    
    cp <- matrix (c(c(0, 0.4 * PM.radius, 0.5 * PM.radius),
                    xyz.L.eye.G + c(0, 20, 0), xyz.L.eye.G), ncol=3, byrow=TRUE)
    ijk <- get.ijk.from.xyz(Sp.way(cp), b) + 1
    b.on.s.L <- b
    b.on.s.L <- fat (b.on.s.L, ijk, 3, TRUE)
    b.on.L <- b
    b.on.L<- fat (b.on.R, ijk, 1, TRUE)
    
    
    b.on.s <- bin.sum (b.on.s.L,b.on.s.R)
    
    if (f.brain) {
      b.on.s.L <- bin.subtraction(b.on.s.L, b.brain.int)
      b.on.s.R <- bin.subtraction(b.on.s.R, b.brain.int)
      b.brain.ext <- bin.subtraction(b.brain.ext, b.on.s)
    }
    if (f.eye) {
      b.on.s.L <- bin.subtraction(b.on.s.L, b.eye)
      b.on.s.R <- bin.subtraction(b.on.s.R, b.eye)
      b.eyes_shield <- bin.subtraction(b.eyes_shield, b.on.s)
    }
    b.on.s.R <- bin.subtraction(b.on.s.R, b.on.R)
    b.on.s.L <- bin.subtraction(b.on.s.L, b.on.L)
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.on.s.L, roi.name = "left optical nerve" , roi.nb = 6,
                                          roi.color = "#8f0000", roi.type = "ORGAN", external.only = FALSE))
      S <- struct.merge(S,struct.from.bin(b.on.s.R, roi.name = "right optical nerve" , roi.nb = 7,
                                          roi.color = "#008f00", roi.type = "ORGAN", external.only = FALSE))
    }
    b.on <- bin.sum (b.on.L,b.on.R)
    rm(list = c("b.on.L","b.on.R","b.on.s.L","b.on.s.R"))
    # display.plane (b.on.s, view.coord = xyz.L.eye.G[3]-3, interpolate = FALSE, main = "ON shield R")
    # display.plane (b.on.s, view.coord = xyz.L.eye.G[3]-3, interpolate = FALSE, main = "ON shield L")
    # display.plane (b.on, view.coord = xyz.L.eye.G[3]-3, interpolate = FALSE, main = "ON")
  }
  
  f.gc <- "ghost container" %in% roi.name
  b.gc.ext <- NULL
  b.gc.int <- NULL
  if (f.gc){
    # ghost container exterior
    PM.gc.radius <- 20
    PM1 <- (xyz[, 1] / (2.2 * PM.gc.radius))^2 +
      (xyz[, 2] / (2.2 * PM.gc.radius))^2 +
      ((xyz[, 3] + 0.7 * PM.radius) / PM.gc.radius)^2 <= 1
    b.gc.ext <- b
    b.gc.ext$vol3D.data[PM1] <- TRUE
    # display.plane (b.gc.ext, view.type  = "sagi", main = "gc ext")
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.gc.ext, roi.name = "ghost container" , roi.nb = 8,
                                          roi.color = "#ff8800", roi.type = "ORGAN", external.only = TRUE))}
    
    ## ghost container interior
    PM1 <- (xyz[, 1] / (2.2 * PM.gc.radius - 3))^2 +
      (xyz[, 2] / (2.2 * PM.gc.radius - 3))^2 +
      ((xyz[, 3] + 0.7 * PM.radius) / (PM.gc.radius - 3))^2 <= 1
    b.gc.int <- b
    b.gc.int$vol3D.data[PM1] <- TRUE
    # display.plane (b.gc.int, view.type  = "sagi", main = "gc int")
  }
  
  f.lpu <- "labyrinth processing unit"  %in% roi.name
  b.lpu <- NULL
  if (f.lpu){
    # Labyrinth processing unit
    t <- seq (0, 2 * pi, by=0.01)
    x <- 10 * (cos (t) + 2 * cos (2 * t))
    y <- 6 * (sin (t) - 2 * sin (2 * t)) + 0.4 * PM.radius
    z <- 6 * 2 * sin (3 * t) + 0.5 * PM.radius
    xyz.lpu <- cbind (x, y, z)
    ijk <- get.ijk.from.xyz(xyz.lpu, b) + 1
    b.lpu <- b
    b.lpu <- fat (b.lpu, ijk, 2, TRUE)
    if (rtstruct.f) {
      S <- struct.merge(S, struct.from.bin(b.lpu, roi.name = "labyrinth processing unit" , 
                                           roi.nb = 9,
                                           roi.color = "#f759cb", roi.type = "ORGAN", external.only = TRUE))}
    # display.plane (b.lpu, view.type  = "sagi", main = "lab unit")
  }
  
  # gizzard shield
  f.gizzard <- "gizzard" %in% roi.name
  b.gz.ext <- NULL
  b.gz.int <- NULL
  if  (f.gizzard){
    cp <- matrix (c(0, -20, 0,
                    0, 10, 0,
                    0, 20, 0,
                    0, 20, -30,
                    20, 20, -30,
                    20, -10, -30,
                    -20, -10, -30,
                    -20, -10, -50,
                    -20, 20, -50,
                    20, 20, -50,
                    0, 0, -30,
                    0, 0, -55,
                    0, 0, -80), ncol=3, byrow=TRUE)
    d <- 0
    for (i in 2:nrow (cp)) d <- c (d, sqrt (sum ((cp[i, ] - cp[i-1, ])^2)))
    d.tot <- cumsum (d)
    df <- as.data.frame (cbind (cp, d.tot))
    x <- spline (df$d.tot, df$V1, n=max (df$d.tot) * 10)$y
    y <- spline (df$d.tot, df$V2, n=max (df$d.tot) * 10)$y
    z <- spline (df$d.tot, df$V3, n=max (df$d.tot) * 10)$y
    xyz.gz <- cbind (x, y, z)
    ijk <- get.ijk.from.xyz(xyz.gz, b) + 1
    b.gz.ext <- b
    b.gz.ext <- fat (b, ijk, 5, TRUE)
    # display.plane (b.gz.ext, view.type  = "sagi", main = "gz shield")
    
    b.gz.int <- b
    b.gz.int <- fat (b, ijk, 2, TRUE)
    
    b.gz.int <- bin.intersection (b.gz.int, b.PM.ext)
    b.gz.ext <- bin.intersection (b.gz.ext,  b.PM.ext)
    if (f.gc){ 
      b.gc.ext <- bin.subtraction (b.gc.ext, b.gz.int)
      b.gz.int <- bin.subtraction (b.gz.int, b.gc.int)
      b.gz.ext <- bin.subtraction (b.gz.ext, b.gc.int)
    }
    
    
    
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.gz.ext, roi.name = "gizzard", roi.nb = 10,
                                          roi.color = "#ffc914", roi.type = "ORGAN", external.only = TRUE))}
  }
  # display.plane (b.gz.int, view.type  = "sagi", main = "gz int")
  
  ## exhaust valve
  f.ev <- "exhaust valve" %in% roi.name
  b.ev <- NULL
  if (f.ev){
    PM1 <- xyz[, 3] <= -0.7 * PM.radius & xyz[, 3] >= -PM.radius &
      xyz[, 1]^2 + xyz[, 2]^2 <= PM.shielding^2
    b.ev <- b
    b.ev$vol3D.data[PM1] <- TRUE
    
    b.ev <- bin.intersection (b.ev, b.PM.ext)
    if (f.gc){ 
      b.gc.ext <-  bin.subtraction (b.gc.ext, b.ev)
      b.ev <- bin.subtraction (b.ev, b.gc.int)
    }
    
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.ev, roi.name = "exhaust valve", roi.nb = 11,
                                          roi.color = "#80ffff", roi.type = "ORGAN", external.only = TRUE))}
  }
  # display.plane (b.ev, view.type  = "sagi", main = "gz int")
  ## energy unit
  f.eu <- "energy unit" %in% roi.name
  b.eu <- NULL
  
  if (f.eu){
    PM.eu.G <- 0.6 * PM.radius * c(0, 1, -0.4)
    PM1 <- abs (xyz[, 1] - PM.eu.G[1]) < 0.1 * PM.radius &
      abs (xyz[, 2] - PM.eu.G[2]) < 0.1 * PM.radius &
      abs (xyz[, 3] - PM.eu.G[3]) < 0.2 * PM.radius
    b.eu <- b
    b.eu$vol3D.data[PM1] <- TRUE  
    
    if (rtstruct.f) {
      S <- struct.merge(S,struct.from.bin(b.eu, roi.name = "energy unit", roi.nb = 12,
                                          roi.color = "#fa4305", roi.type = "ORGAN", external.only = TRUE))}
    
    ## chaos energy flow eu to ghost container
    if (f.gc){
      for (i in 1:20) {
        src <- PM.eu.G - c (0, 0, 0.2 * PM.radius)
        t <- runif (1, -pi, pi)
        dst <- c (0, 0, -0.7 * PM.radius) + c(runif (1, 0, 50) * c(cos (t), sin (t)), 0)
        n <- dst - src
        l <- sqrt (sum (n^2))
        n <- n / l
        u_ <- c(n[2], n[3], n[1])
        u <- vector.product (n, u_)
        v <- vector.product (n, u)
        l.s <- seq (0, l, by=.5)
        r <- cumsum (c (0, sample (seq (-l/20, l/20, length.out = length (l.s)-1))))
        t <- cumsum (runif (length (r), -pi/30, pi/30))
        M <- t (apply (cbind (r, l.s, t), 1, function (rst) src + rst[2] * n + rst[1] * cos (rst[3]) * u + rst[1] * sin (rst[3]) * v))
        keep <- apply (M, 1, function (v) sum(v^2) <= (PM.radius - PM.shielding)^2)
        M <- M[keep, ]
        ijk <- get.ijk.from.xyz (M, b) + 1
        b.eu$vol3D.data <- b.eu$vol3D.data | fat (b, ijk, 2, TRUE)$vol3D.data
      }
    }
    
    ## chaos energy flow eu to brain
    if (f.brain){
      for (i in 1:20) {
        src <- PM.eu.G + c (0, 0, 0.2 * PM.radius)
        t <- runif (1, -pi, pi)
        dst <- c (0, 0.4 * PM.radius, 0.5 * PM.radius) + c(runif (1, 0, PM.brain.radius) *
                                                             c(2 * cos (t), sin (t)), 0)
        n <- dst - src
        l <- sqrt (sum (n^2))
        n <- n / l
        u_ <- c(n[2], n[3], n[1])
        u <- vector.product (n, u_)
        v <- vector.product (n, u)
        l.s <- seq (0, l, by=.5)
        r <- cumsum (c (0, sample (seq (-l/20, l/20, length.out = length (l.s)-1))))
        t <- cumsum (runif (length (r), -pi/30, pi/30))
        M <- t (apply (cbind (r, l.s, t), 1, function (rst) src + rst[2] * n + rst[1] * cos (rst[3]) * u + rst[1] * sin (rst[3]) * v))
        keep <- apply (M, 1, function (v) sum(v^2) <= (PM.radius - PM.shielding)^2)
        M <- M[keep, ]
        ijk <- get.ijk.from.xyz (M, b) + 1
        b.eu$vol3D.data <- b.eu$vol3D.data | fat (b, ijk, 2, TRUE)$vol3D.data
      }
    }
  }
  
  
  # PTV
  b.ptv <- NULL
  # if  (rtdose.f | rtstruct.f){
  PM1 <- ((xyz[, 1] - ptv.G[1]) / ptv.radius)^2 +
    ((xyz[, 2] - ptv.G[2]) / ptv.radius)^2 +
    ((xyz[, 3] - ptv.G[3]) / ptv.radius)^2 <= 1
  b.ptv <- b
  b.ptv$vol3D.data[PM1] <- TRUE
  if (rtstruct.f) {
    S <- struct.merge(S,struct.from.bin(b.ptv, roi.name = "ptv" , roi.nb = 13,
                                        roi.color = "#ff80ff", roi.type = "PTV",
                                        external.only = TRUE))
  }
  # }
  
  # display.plane (b.eu, view.type  = "sagi", view.coord = PM.eu.G[1], main = "energy")
  b.bone <- bin.subtraction (bin.inversion (b.ext.mask), b.PM.int)
  ## CT computation
  
  fill <- function (src, bin, value, sd) {
    src$vol3D.data[bin$vol3D.data] <- rnorm (sum (bin$vol3D.data), value, sd)
    return (src)
  }
  
  ct.n <- paste0(CT.date,"_ref1_do",ref1.do.nb[1],"_ct")
  if (ct.f | rtdose.f){
    n <-ct.n
    pat$ct <- list()
    pat$ct[[paste0(n,1)]] <-
      vol.create (n.ijk = n.ijk, dxyz = dxyz, pt000=pt000, modality = "ct",
                  default.value = -1000, ref.pseudo = "ref1", 
                  description = "FULL^BODY|initial",
                  frame.of.reference = "2.16.840.1.114357.58485835081.6250.1551498438.0.3",
                  number = 1)
    
    pat$ct[[1]]$patient <- pat$pat.pseudo
	pat$ct[[1]]$patient.name <- pat$patient$name[1]
    pat$ct[[1]]$patient.bd <- pat$patient$birth.date[1]
    pat$ct[[1]]$patient.sex <- pat$patient$sex[1]
    pat$ct[[1]]$object.name <- n
    pat$ct[[1]]$object.alias <- paste0(n,1)
    pat$ct[[1]]$acq.date = CT.date
    pat$ct[[1]]$study.date = CT.date
    pat$ct[[1]]$creation.date = CT.date
    
    pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.PM.int, -50, 10)
    pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.bone, 800, 10)
    
    if (!is.null(b.eu)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.eu, 80, 10)
    if (!is.null(b.gc.ext)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.gc.ext, 800, 10)
    if (!is.null(b.brain.ext)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.brain.ext, 800, 10)
    if (!is.null(b.gz.ext)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.gz.ext, 800, 10)
    if (!is.null(b.gz.int)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.gz.int, -1000, 0)
    if (!is.null(b.ev)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.ev, 800, 10)
    if (!is.null(b.on.s)) pat$ct[[1]] <- fill ( pat$ct[[1]] , b.on.s, 800, 10)
    if (!is.null(b.on)) pat$ct[[1]] <- fill ( pat$ct[[1]] , b.on, 10, 10)
    if (!is.null(b.gc.int)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.gc.int, -1000, 0)
    if (!is.null(b.brain.int)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.brain.int, 10, 10)
    if (!is.null(b.lpu)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.lpu, 150, 30)
    if (!is.null(b.eye)) pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.eye, -10, 10)
    if (!is.null(b.ptv)) pat$ct[[1]]$vol3D.data <- pat$ct[[1]]$vol3D.data + b.ptv$vol3D.data*10
    
    pat$ct[[1]]  <- fill ( pat$ct[[1]] , bin.inversion (b.PM.ext), -1000, 0)
    pat$ct[[1]]  <- fill ( pat$ct[[1]] , b.couch, 0, 5)
    
    pat$ct[[1]]$max.pixel <- max (pat$ct[[1]] $vol3D.data)
    pat$ct[[1]]$min.pixel <- min (pat$ct[[1]] $vol3D.data)
    
    if (rtdose.f) {
      pat$rtdose <- list ()
      n <-paste0(rt.date,"_ref1_do",ref1.do.nb[2],"_rtdose")
      # .toy.rtdose  (pat$ct[[1]], ptv.radius, ptv.G, beam.nb)
      

      pat$rtdose[[paste0(n,1)]] <- .toy.rtdose.from.bin(D.max = 52, DSA= 600, 
                                                        beam.nb = beam.nb, 
                                                        vol = pat$ct[[1]], 
                                                        bin = b.ptv, 
                                                        description =  "IMRT|PTV52")
        
      pat$rtdose[[1]]$patient <- pat$pat.pseudo
	  pat$rtdose[[1]]$patient.name <- pat$patient$name[1]
      pat$rtdose[[1]]$patient.bd <- pat$patient$birth.date[1]
      pat$rtdose[[1]]$patient.sex <- pat$patient$sex[1]
      pat$rtdose[[1]]$object.name <- n
      pat$rtdose[[1]]$object.alias <- paste0(n,1)
      if (rtstruct.f) pat$rtdose[[1]]$ref.object.alias <- ""
      pat$rtdose[[1]]$acq.date = CT.date
      pat$rtdose[[1]]$study.date = CT.date
      pat$rtdose[[1]]$creation.date = CT.date
      pat$rtdose[[1]]$unit = "GY"
      
    }
    if (!ct.f) pat$ct <- NULL
  }
  
  if (mr.f){
    
    pat$mr <- list()
    
    b <- vol.create (n.ijk = n.ijk, dxyz = dxyz, pt000=pt000 + mr.offset,default.value = 0)
    b <- fill (b, b.PM.int, 100, 30)
    b <- fill (b, b.bone, 20, 3)
    
    if (!is.null(b.eu)) b <- fill (b, b.eu, 50, 10)
    if (!is.null(b.gc.ext)) b <- fill (b, b.gc.ext, 20, 3)
    if (!is.null(b.brain.ext)) b <- fill (b, b.brain.ext, 20, 3)
    if (!is.null(b.gz.ext)) b <- fill (b, b.gz.ext, 20, 3)
    if (!is.null(b.gz.int)) b <- fill (b, b.gz.int, 5, 1)
    if (!is.null(b.ev)) b <- fill (b, b.ev, 20, 3)
    if (!is.null(b.on.s)) b <- fill (b, b.on.s, 20, 3)
    if (!is.null(b.on)) b <- fill (b, b.on, 200, 50)
    if (!is.null(b.gc.int)) b <- fill (b, b.gc.int, 5, 1)
    if (!is.null(b.brain.int)) b <- fill (b, b.brain.int, 100, 30)
    if (!is.null(b.lpu)) b <- fill (b, b.lpu, 200, 30)
    if (!is.null(b.eye)) b <- fill (b, b.eye, 80, 30)
    if (!is.null(b.ptv)) b$vol3D.data <- b$vol3D.data + b.ptv$vol3D.data*200
    
    b <- fill (b, bin.inversion (b.PM.ext), 5, 1)
    b <- fill (b, b.couch, 20, 5)
    b$vol3D.data[b$vol3D.data < 0] <- 0
    b$max.pixel <- max (b$vol3D.data)
    b$min.pixel <- min (b$vol3D.data)
    
    perm <- c(3,1,2)
    m <- as.matrix (expand.grid(1:b$n.ijk[1],
                                1:b$n.ijk[2],
                                1:b$n.ijk[3]))
    n <-paste0(MR.date,"_ref2_do1_mr")
    pat$mr[[paste0(n,1)]] <- 
      vol.create (n.ijk = n.ijk[perm], dxyz = dxyz[perm], pt000=b$xyz0[1,], 
                  modality = "mr", default.value = 0, ref.pseudo = "ref2", 
                  description = "FULL^BODY|T1 GADO",
                  frame.of.reference = "2.16.840.1.114357.58094835081.62350.15267498438.0.1",
                  number = 1)
    pat$mr[[1]]$orientation <- as.vector(diag(3)[,perm])[1:6]
    pat$mr[[1]]$xyz.from.ijk[,1:3] <-b$xyz.from.ijk[,perm]   
    pat$mr[[1]]$xyz0  <- matrix((as.matrix (expand.grid (0, 0, pat$mr[[1]]$k.idx,1))%*% 
                                           t(pat$mr[[1]]$xyz.from.ijk))[ ,1:3],ncol=3)
    m <- m[order(m[,perm[3]],m[,perm[2]],m[,perm[1]]),]
    pat$mr[[1]]$vol3D.data <-array(b$vol3D.data[m], dim = pat$mr[[1]]$n.ijk)
    pat$mr[[1]]$max.pixel <- b$max.pixel
    pat$mr[[1]]$min.pixel <- b$min.pixel
    
    pat$mr[[1]]$patient <- pat$pat.pseudo
	pat$mr[[1]]$patient.name <- pat$patient$name[1]
    pat$mr[[1]]$patient.bd <- pat$patient$birth.date[1]
    pat$mr[[1]]$patient.sex <- pat$patient$sex[1]
    pat$mr[[1]]$object.name <- n
    pat$mr[[1]]$object.alias <- paste0(n,1)
    pat$mr[[1]]$acq.date = MR.date
    pat$mr[[1]]$study.date = MR.date
    pat$mr[[1]]$creation.date = MR.date
    
    rm(b)
    
  }
  
  if (rtstruct.f) {
    pat$rtstruct <- list ()
    n <-paste0(rt.date,"_ref1_do",ref1.do.nb[3],"_rtstruct")
    pat$rtstruct[[paste0(n,1)]] <- S
    pat$rtstruct[[1]]$patient <- pat$pat.pseudo
	pat$rtstruct[[1]]$patient.name <- pat$patient$name[1]
    pat$rtstruct[[1]]$patient.bd <- pat$patient$birth.date[1]
    pat$rtstruct[[1]]$patient.sex <- pat$patient$sex[2]
    pat$rtstruct[[1]]$object.name <- n
    pat$rtstruct[[1]]$object.alias <- paste0(n,1)
    if (rtdose.f) pat$rtdose[[1]]$ref.object.alias <- pat$rtstruct[[1]]$object.alias
    if (ct.f) pat$rtstruct[[1]]$ref.object.alias <- pat$ct[[1]]$object.alias
    pat$rtstruct[[1]]$description <- "TREATMENT|RS: Approved Structure Set"
    pat$rtstruct[[1]]$ref.object.alias <-  ct.n
    pat$rtstruct[[1]]$approval.status <- "APPROVED"
    pat$rtstruct[[1]]$study.date = rt.date
    pat$rtstruct[[1]]$creation.date = rt.date
    
    l <- length(pat$rtstruct[[1]])
    pat$rtstruct[[1]] <- pat$rtstruct[[1]] [c (1:7,l-1, 8:15,l,16:(l-2))]
    class(pat$rtstruct[[1]]) <- "struct"
  }
  
  pat$description <- do.call(rbind.data.frame, 
                             lapply(do.call(c, pat[c("ct", "mr", "rtdose", "rtstruct")]), function(l) {
                               nb <- switch(l$modality, "rtstruct" =  l$nb.of.roi,  
                                            "rtdose" = l$n.ijk[3], "ct" = l$n.ijk[3],
                                            "mr" = l$n.ijk[3],NA)
                               max.pix <- NA
                               idx <- grep ("max[.]pixel", names(l))
                               if (length(idx)>0)  max.pix <- l[[idx]]
                               c(l$patient, l$modality, l$object.name, l$ref.pseudo,
                                 1, l$description, nb,  max.pix, l$object.alias)
                             }))
  
  if (nrow(pat$T.MAT$ref.info)==2)
    pat$description <- rbind(pat$description, c(pat.pseudo, "reg", "ref1_from_ref2", 
                                                "ref1",1,"ref1 from ref2",2,NA,
                                                "ref1_from_ref2"))
  colnames(pat$description) <- c ("PIN", "modality", "obj", "ref.pseudo", "nb.of.subobject" ,"description", "nb", "max","object.alias")
  pat$description <- pat$description[order(pat$description$PIN,
                                           pat$description$ref.pseudo,
                                           pat$description$modality),]
  pat$description$max <- suppressWarnings(as.character(round(as.numeric(pat$description$max),3)))
  pat$description$nb <- suppressWarnings(as.numeric(pat$description$nb))
  row.names(pat$description) <- NULL
  
  pat$description.by.reg[[1]] <- pat$description 
  return (pat)
}


