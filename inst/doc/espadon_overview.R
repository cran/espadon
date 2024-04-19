## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  vol <- load.obj.from.dicom (image_path)

## ----eval = FALSE-------------------------------------------------------------
#  dicom.to.Rdcm.converter("D:/dcm/patient001", "D:/Rdcm/patient001“)

## ----eval = FALSE-------------------------------------------------------------
#  dcm.filename <- file.choose ()
#  df <- dicom.parser (dicom.raw.data.loader (dcm.filename), as.txt = TRUE, try.parse = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  library(espadon)
#  pat.dir <- choose.dir () # patient folder containing mr, ct, rt-struct, rt-dose, rt-plan and reg…
#  pat <- load.patient.from.dicom (pat.dir)
#  display.obj.links (pat)

## ----eval = FALSE-------------------------------------------------------------
#  S <- load.obj.data (pat$rtstruct[[1]])
#  CT <- load.obj.data (pat$ct[[1]])
#  MR <- load.obj.data (pat$mr[[1]])
#  D <- load.obj.data (pat$rtdose[[1]])

## ----eval = FALSE-------------------------------------------------------------
#  CT <- load.obj.from.dicom(ct_dicom_file_path)
#  display.kplane (CT, k =10, col = pal.RVV(1000),
#                  breaks = seq(-1000, 1000, length.out = 1001),
#                  ord.flip = TRUE, interpolate = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  library(espadon)
#  library(rgl)
#  pat.dir <- choose.dir () # patient folder containing mr, ct, rt-struct, and reg
#  pat <- load.patient.from.dicom (pat.dir)
#  S <- load.obj.data(pat$rtstruct[[1]])
#  CT <- load.obj.data(pat$ct[[1]])
#  MR <- load.obj.data(pat$mr[[1]])
#  
#  bg3d ("black")
#  par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1), ncol = 4),
#         windowRect = c(0, 50, 300, 300), zoom = 0.7)
#  
#  display.3D.contour (S, roi.sname = "eye", display.ref = CT$ref.pseudo,
#                      T.MAT = pat$T.MAT)
#  display.3D.stack (MR, k.idx = c(10,35,70, 105, 140,165),
#                    display.ref = CT$ref.pseudo,
#                    T.MAT = pat$T.MAT,
#                    ktext = FALSE)
#  display.3D.sections (CT, cross.pt = c(0, 150, 0),
#                       col = pal.RVV (200, alpha = c(rep (0, 90), rep (1, 110))),
#                       breaks = seq(-1000, 1000, length.out = 201),
#                       border.col = "#844A39")
#  
#  play3d (spin3d (rpm = 4), duration = 15)

## ----eval = FALSE-------------------------------------------------------------
#  nesting.MR <- nesting.roi (MR, S, roi.sname = "brain", T.MAT = pat$T.MAT,
#                             vol.restrict = TRUE, mar)
#  display.3D.stack (nesting.MR, ktext = FALSE, line.lw d= 2)
#  display.3D.contour (S, roi.sname = "brain", display.ref = MR$ref.pseudo,
#                      T.MAT = pat$T.MAT)

## ----eval = FALSE-------------------------------------------------------------
#  MR.on.CT <- vol.regrid (MR, CT, T.MAT=pat$T.MAT, interpolate = TRUE)
#  display.3D.stack (MR.on.CT, display.ref = CT$ref.pseudo , T.MAT = pat$T.MAT,
#                    ktext = FALSE, line.lwd = 2, line.col = "#844A39")

## ----eval = FALSE-------------------------------------------------------------
#  # Points determination : centers of the eyes and the chiasm found in rtstruct data
#  roi.idx <- select.names (S$roi.info$roi.pseudo, roi.name c("eye", "chiasm"))
#  pt <- S$roi.info[roi.idx, c("Gx", "Gy", "Gz")]
#  origin.pt <- as.numeric (apply (pt, 2, mean))
#  
#  # Plane creation in the CT frame of reference
#  new.plane <- get.plane (CT, origin = origin.pt,
#                          plane.orientation =  orientation.create (pt),
#                          interpolate = TRUE)
#  
#  # Display of the new plane, in the cut plane frame of reference
#  tmat <- ref.cutplane.add (new.plane, origin.pt = origin.pt,
#                            ref.cutplane = "ref.plane", T.MAT = pat$T.MAT)
#  plan <- vol.in.new.ref (new.plane, "ref.plane", tmat)
#  display.plane (plan, struct= S, roi.idx = idx, T.MAT = tmat,
#                 view.coord = plan$patient.xyz0[3],
#                 bottom.col = pal.RVV(200),
#                 bottom.breaks = seq(-1000, 1000, length.out = 201),
#                 main = "cut plane", bg = "#379DA2", sat.transp = T,
#                 legend.plot = FALSE, interpolate = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  nesting.MR <- nesting.roi (MR, S, roi.sname = "brain", T.MAT = pat$T.MAT,
#                             vol.restrict = TRUE,
#                             xyz.margin = c(20, 20, 20)) # to reduce the computation time
#  bin.brain <- bin.from.roi (nesting.MR, struct = S, roi.sname = "brain",
#                             T.MAT = pat$T.MAT)
#  MR. brain <- vol.from.bin (nesting.MR, bin.brain)

## ----eval = FALSE-------------------------------------------------------------
#  bin.min.150 <- bin.from.vol (MR.brain, min = 150)
#  
#  display.plane (bin.min.150, view.coord = 0, view.type = "sagi",
#                 main = "bin.min.150", bg = "#379DA2", legend.plot = FALSE,
#                 sat.transp = TRUE, interpolate = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  bin.smooth <- bin.opening (bin.min.150, radius = 1.5))
#  bin.cluster <- bin.clustering (bin.smooth)
#  
#  # Creation of ventricle contours
#  bin.ventricle <- bin.from.vol (bin.cluster, min = 1, max = 1)
#  S.ventricle <- struct.from.bin (bin.ventricle , roi.name = "ventricle",
#                                  roi.color = "#FF0000" )
#  
#  # Display of the first 3 largest sub-volumes, and ventricle contours
#  library (rgl)
#  bg3d ("black")
#  par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
#                              ncol = 4),
#         windowRect = c(0, 50, 300, 300), zoom = 0.7)
#  display.3D.stack(bin.cluster, k.idx = bin.cluster$k.idx,
#                   col = c ("#00000000", rainbow (3, alpha = 1)),
#                   breaks = seq (0, 3, length.out = 5), border = FALSE,
#                   ktext = FALSE)
#  play3d (spin3d (rpm = 4), duration = 15)
#  
#  bg3d ("black")
#  par3d (userMatrix = matrix (c(0, 0, -1, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1),
#                              ncol = 4),
#         windowRect = c(0, 50, 300, 300), zoom = 1)
#  display.3D.contour (S.ventricle)
#  play3d (spin3d (rpm = 4), duration = 15)
#  

## ----eval = FALSE-------------------------------------------------------------
#  library (espadon)
#  patient.folder <- choose.dir ()
#  pat <- load.patient.from.dicom (patient.folder, data = TRUE)
#  
#  bin.lung <- bin.sum (bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]],
#                                     roi.name = "lungl"),
#                       bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]],
#                                     roi.name = "lungr"))
#  bin.heart <- bin.from.roi (pat$ct[[1]], struct = pat$rtstruct[[1]],
#                             roi.name = "heart")
#  
#  mesh.lung <- mesh.from.bin (bin.lung)
#  mesh.heart <- mesh.from.bin (bin.heart)
#  
#  display.3D.mesh (mesh.lung, color = "burlywood2", specular = "black",
#                   alpha = 0.8)
#  display.3D.mesh (mesh.heart, color = "red", specular = "black",
#                   alpha = 1)

## ----eval = FALSE-------------------------------------------------------------
#  library (espadon)
#  patient.folder <- choose.dir ()
#  pat <- load.patient.from.dicom (patient.folder, data = FALSE)
#  
#  S <- load.obj.data (pat$rtstruct[[1]])
#  CT <- load.obj.data (pat$ct[[1]])
#  D <- load.obj.data (pat$rtdose[[1]])
#  
#  # resample D on the cut planes on which the chiasm has been outlined
#  D.on.CT <- vol.regrid (D, CT)
#  
#  # Calculation of histogram with simulation of random displacements according to
#  # a normal distribution of standard deviation 1 mm in all directions:
#  h <- histo.from.roi (D.on.CT, S, roi.name="chiasm", breaks = seq(0, 80, 0.1),
#                       alias = "chiasm", MC = 100 )
#  
#  DVH <- histo.DVH (h, alias = "Chiasm DVH")
#  
#  display.DVH (DVH, MC.plot = TRUE, col = "#ff0000", lwd = 2)

## ----eval = FALSE-------------------------------------------------------------
#  library (espadon)
#  patient.folder <- choose.dir ()
#  pat <- load.patient.from.dicom (patient.folder, data = FALSE)
#  
#  S <- load.obj.data (pat$rtstruct[[1]])
#  CT <- load.obj.data (pat$ct[[1]])
#  D <- load.obj.data (pat$rtdose[[1]])
#  
#  # resample D on the cut planes on which the chiasm has been outlined
#  D.on.CT <- vol.regrid (D, CT)
#  
#  rt.indices.from.roi (D.on.CT, S, target.roi.name = "ptv",
#                       healthy.roi.name = "chiasm",
#                       presc.dose = 0.9 * 70)

## ----eval = FALSE-------------------------------------------------------------
#  library (espadon)
#  patient.folder <- choose.dir ()
#  pat <- load.patient.from.dicom (patient.folder, data = FALSE)
#  
#  S <- load.obj.data (pat$rtstruct[[1]])
#  CT <- load.obj.data (pat$ct[[1]])
#  
#  # The binary objects for 2 ROIs named "roi A" and "roi B" respectively :
#  bin.A <- bin.from.roi (CT, S, roi.sname = "roi A", alias = "ROI A")
#  bin.B <- bin.from.roi (CT, S, roi.sname = "roi B", alias = "ROI B")
#  
#  # The mesh objects for these 2 ROIs:
#  # FYI: the smooth.iteration argument avoids staircase meshing, thanks to z-axis binning.
#  #      As a result, ROI contours are also smoothed in XY.
#  mesh.A <- mesh.from.bin(bin.A, smooth.iteration = 10, alias = "ROI A")
#  mesh.B <- mesh.from.bin(bin.B, smooth.iteration = 10, alias = "ROI B")
#  
#  # spatial similarity
#  sp.similarity.from.bin (bin.A , bin.B)
#  sp.similarity.from.mesh (mesh.A , mesh.B, hausdorff.quantile = c (0.5, 0.95), surface.DSC.tol = 1:3)

