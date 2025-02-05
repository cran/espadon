## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6, fig.height = 4, out.width = "90%", fig.align = 'center'
)

## -----------------------------------------------------------------------------
par(cex = 0.5, title.cex=0.8)

## ----eval = FALSE-------------------------------------------------------------
# library(espadon)
# vol <- load.obj.from.dicom (image_path)

## ----eval = FALSE-------------------------------------------------------------
# library(espadon)
# dicom.to.Rdcm.converter("D:/dcm/patient001", "D:/Rdcm/patient001“)

## ----eval = FALSE-------------------------------------------------------------
# library(espadon)
# dcm.filename <- file.choose ()
# df <- dicom.parser (dicom.raw.data.loader (dcm.filename), as.txt = TRUE,
#                     try.parse = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# library(espadon)
# pat.dir <- choose.dir () # patient folder containing mr, ct, rt-struct, rt-dose, rt-plan and reg…
# pat <- load.patient.from.dicom (pat.dir)
# display.obj.links (pat)

## ----eval = FALSE-------------------------------------------------------------
# S <- load.obj.data (pat$rtstruct[[1]])
# CT <- load.obj.data (pat$ct[[1]])
# MR <- load.obj.data (pat$mr[[1]])
# D <- load.obj.data (pat$rtdose[[1]])

## ----eval = TRUE--------------------------------------------------------------
library (espadon)
pat <- toy.load.patient(modality = c("ct", "mr", "rtdose", "rtstruct"),
                        dxyz = c(2, 2, 2),beam.nb = 3)
S <- pat$rtstruct[[1]]
CT <- pat$ct[[1]]
MR <- pat$mr[[1]]
D <- pat$rtdose[[1]]
display.obj.links (pat)

## ----echo = TRUE, eval = TRUE-------------------------------------------------
par(mar = c(4, 4, 2, 1), cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)
display.kplane (CT, k =10, col = pal.RVV(1000), 
                flop = TRUE, interpolate = TRUE)

ptv.center <- as.numeric(S$roi.info[S$roi.info$name == "ptv", c("Gx","Gy","Gz")])


display.plane(CT, D, S, main ="transverse", view.coord = ptv.center[3], 
              legend.shift = -100)
display.plane(CT, D, S, main ="frontal", view.type = "front", 
              view.coord = ptv.center[2], legend.shift = -100)
display.plane(CT, D, S, main ="sagittal", view.type = "sagi",
              view.coord = ptv.center[1], legend.shift = -100)

## ----echo = TRUE, eval = TRUE-------------------------------------------------
par(mar = c(4, 4, 2, 1), cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.9)
plot(CT,view.type = "yz", view.coord = ptv.center, col = pal.RVV(100), las = 1)
plot(D, add = TRUE, view.type = "yz", view.coord = ptv.center, 
     col= pal.rainbow((100)))
Scut <- plot(S, add = TRUE, view.type = "yz", view.coord = ptv.center, lwd = 2)
legend("topleft",legend = Scut$roi.info$name, text.col = Scut$roi.info$color,
       cex = 0.8, bty="n")

## ----eval = FALSE-------------------------------------------------------------
# 
# library(rgl)
# bg3d ("black")
# par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), ncol = 4),
#        windowRect = c(0, 50, 300, 300), zoom = 0.7)
# 
# display.3D.contour (S, roi.sname = "eye", display.ref = CT$ref.pseudo,
#                     T.MAT = pat$T.MAT)
# display.3D.stack (MR, k.idx = round(seq(3, length(MR$k.idx)-3,length.out=10)),
#                   display.ref = CT$ref.pseudo,
#                   T.MAT = pat$T.MAT,
#                   ktext = FALSE)
# display.3D.sections (CT, cross.pt = c(0, 110, 0),
#                      col = pal.RVV (200, alpha = c(rep (0, 90), rep (1, 110))),
#                      breaks = seq(-1000, 1000, length.out = 201),
#                      border.col = "#844A39")
# 
# rglwidget()
# close3d()

## ----eval = FALSE-------------------------------------------------------------
# nesting.MR <- nesting.roi (MR, S, roi.sname = "brain", T.MAT = pat$T.MAT,
#                            vol.restrict = TRUE, xyz.margin = c(2, 2, 2))
# display.3D.stack (nesting.MR, ktext = FALSE, line.lwd= 2)
# display.3D.contour (S, roi.sname = "brain", display.ref = MR$ref.pseudo,
#                     T.MAT = pat$T.MAT)
# 

## ----eval = FALSE-------------------------------------------------------------
# MR.on.CT <- vol.regrid (MR, CT, T.MAT=pat$T.MAT, interpolate = TRUE)
# display.3D.stack (MR.on.CT, display.ref = CT$ref.pseudo , T.MAT = pat$T.MAT,
#                   ktext = FALSE, line.lwd = 2, line.col = "#844A39")

## ----eval = TRUE--------------------------------------------------------------
# Points determination : centers of the eyes and the chiasm found in rtstruct data
roi.idx <- select.names (S$roi.info$roi.pseudo, 
                         roi.name = c("left eye" ,"right eye","ptv"))
pt <- S$roi.info[roi.idx, c("Gx", "Gy", "Gz")]
origin.pt <- as.numeric (apply (pt, 2, mean))

# Plane creation in the CT frame of reference 
new.plane <- get.plane (CT, origin = origin.pt, 
                        plane.orientation =  orientation.create (pt), 
                        interpolate = TRUE)

# Display of the new plane, in the cut plane frame of reference
tmat <- ref.cutplane.add (new.plane, origin = origin.pt, 
                          ref.cutplane = "ref.plane", T.MAT = pat$T.MAT)
plan <- vol.in.new.ref (new.plane, "ref.plane", tmat)
display.plane (plan, struct= S, roi.idx = roi.idx, T.MAT = tmat, 
               view.coord = plan$xyz0[3],
               bottom.col = pal.RVV(200), 
               bottom.breaks = seq(-1000, 1000, length.out = 201), 
               main = "cut plane", bg = "#379DA2", sat.transp = T, 
               legend.plot = FALSE, interpolate = TRUE)

## ----eval = TRUE--------------------------------------------------------------
nesting.MR <- nesting.roi (MR, S, roi.sname = "patient", T.MAT = pat$T.MAT, 
                           vol.restrict = TRUE, 
                           xyz.margin = c(1, 1, 1)) # to reduce the computation time 
bin.patient <- bin.from.roi (nesting.MR, struct = S, roi.sname = "patient", 
                           T.MAT = pat$T.MAT)
MR.patient <- vol.from.bin (nesting.MR, bin.patient)
plot(MR.patient)

## ----eval = TRUE--------------------------------------------------------------
bin.patient_ <- bin.from.roi (nesting.MR, struct = S, roi.sname = "patient", 
                           T.MAT = pat$T.MAT)
MR.patient_ <- vol.from.bin (nesting.MR, bin.patient_)
par(mfrow=c(1,2))
plot(bin.patient, main = "binary", cut.interpolate  = FALSE)
plot(bin.patient_, main="weight", cut.interpolate  = FALSE)
par(mfrow=c(1,1))

## ----eval = TRUE--------------------------------------------------------------
bin.min.150 <- bin.from.vol (MR.patient, min = 70)

display.plane (bin.min.150, view.coord = 0, view.type = "sagi",
               main = "bin.min.150", bg = "#379DA2", legend.plot = FALSE, 
               sat.transp = TRUE, interpolate = TRUE)

## ----eval = TRUE--------------------------------------------------------------
bin.smooth <- bin.closing (bin.min.150, radius = 2 * max(bin.min.150$dxyz))
bin.cluster <- bin.clustering (bin.smooth)

n=5
col <- c ("#00000000", rainbow (n))
breaks <- seq (-0.5, n + 0.5, length.out = n+2)
display.plane (MR, top = bin.cluster, main = "After clustering", 
               view.coord = 0, view.type = "sagi",
               top.col = col, top.breaks = breaks, 
               interpolate = FALSE)

## ----eval = FALSE-------------------------------------------------------------
# # Creation of new patient contours
# new.bin.patient <- bin.from.vol (bin.cluster, min = 1, max = 1)
# S.patient <- struct.from.bin (new.bin.patient , roi.name = "patient",
#                                 roi.color = "#FF0000" )
# 
# # Display of the first 3 largest sub-volumes, and ventricle contours
# library (rgl)
# bg3d ("black")
# par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
#                             ncol = 4),
#        windowRect = c(0, 50, 300, 300), zoom = 0.7)
# display.3D.stack(bin.cluster, k.idx = bin.cluster$k.idx,
#                  col = c ("#00000000", rainbow (3, alpha = 1)),
#                  breaks = seq (0, 3, length.out = 5), border = FALSE,
#                  ktext = FALSE)
# play3d (spin3d (rpm = 4), duration = 15)
# 
# bg3d ("black")
# par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
#                             ncol = 4),
#        windowRect = c(0, 50, 300, 300), zoom = 1)
# display.3D.contour (S.patient, roi.col="yellow")
# play3d (spin3d (rpm = 4), duration = 15)
# 

## ----eval = TRUE--------------------------------------------------------------

bin.leye <- bin.sum (bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]], 
                                   roi.name = "left eye"),
                     bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]], 
                                   roi.name = "left optical nerve"))

bin.reye <- bin.sum (bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]], 
                                   roi.name = "right eye"),
                     bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]], 
                                   roi.name = "right optical nerve"))

mesh.leye <- mesh.from.bin (bin.leye)
mesh.reye <- mesh.from.bin (bin.reye)

## ----eval = FALSE-------------------------------------------------------------
# display.3D.mesh (mesh.eye, color = "burlywood2", specular = "black",
#                  alpha = 0.8)
# display.3D.mesh (mesh.proc, color = "red", specular = "black",
#                  alpha = 1)

## ----eval = TRUE--------------------------------------------------------------
# resample D on the cut planes on which the chiasm has been outlined
D.on.CT <- vol.regrid (D, CT)

# Calculation of histogram with simulation of random displacements according to 
# a normal distribution of standard deviation 1 mm in all directions:
h <- histo.from.roi (D.on.CT, S, 
                     roi.name="ptv", 
                     breaks = seq(0, 80, 0.1),
                     alias = "ptv", MC = 10)

DVH <- histo.DVH (h, alias = "PTV")

display.DVH (DVH, MC.plot = TRUE, col = "#ff0000", lwd = 2)

## ----eval = TRUE--------------------------------------------------------------
indices <- rt.indices.from.roi (D.on.CT, S, target.roi.name = "ptv", 
                     healthy.roi.sname = c("eye", 'proc', 'ghost'),
                     presc.dose = round(0.9 * D$max.pixel),
                     conformity.indice = "", homogeneity.indices="",
                     gradient.indices = "",
                     V.xpc = c(5, 50, 75), DVH = FALSE)
head(indices$dosimetry)
head(indices$volume)

## ----eval = TRUE--------------------------------------------------------------


# The binary object for brain is :
bin.brain <- bin.from.roi (CT, S, roi.sname = "brain", modality ="weight")

# The mesh objects for brain:
# FYI: the smooth.iteration argument avoids staircase meshing, due to z-axis binning.
#      As a result, ROI contours are also smoothed in XY.
mesh.brain <- mesh.from.bin (bin.brain, smooth.iteration = 10, alias = "brain")


# spatial similarity between lest ocular system and brain
metrics.le <- c(list(sp.similarity.from.bin (bin.leye , bin.brain)),
                sp.similarity.from.mesh (mesh.leye , mesh.brain, 
                                         hausdorff.quantile = c (0.5, 0.95),
                                         surface.tol = 1:3))

metrics.le[[1]]
metrics.le[[2]]
metrics.le[[3]]

