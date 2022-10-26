#' Restrict volume to a binary selection
#' @description The \code{nesting.bin} function restricts a "volume" class 
#' object to the rectangular parallelepiped circumscribed to the selected voxels.
#' @param vol "volume" class object, containing data to restrict.
#' @param sel.bin "volume" class object, of "binary" modality, specifying the selected voxels.
#' @param xyz.margin Vector of length 3, representing the distances in mm to be added
#' to the x, y and z directions of the rectangular parallelepiped circumscribed
#' to the voxels selected in \code{sel.bin}, in the cutting planes frame of reference. 
#' By default xyz.margin = c (0, 0, 0).
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be 
#' \code{paste (vol$description,"restricted to", sel.bin$description)}.
#' @return Returns a "volume" class object, in which 3D volume is limited to the
#' rectangular parallelepiped circumscribed to the voxels selected by \code{sel.bin}, increased by the
#' requested margins.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for  better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name = "brain", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' b <- bin.from.vol (CT, min = 0, max =200)
#' 
#' CT.restricted <- nesting.bin (CT, b, xyz.margin =  rep (step, 3))
#' display.plane (bottom = CT.restricted, top = b, view.type = "sagi",
#'              bottom.col = pal.RVV (1000),
#'              bottom.breaks = seq (-1000, 1000, length.out = 1001),
#'              bg = "#00ff00",  interpolate  = FALSE)
#' @export              
nesting.bin <- function (vol, sel.bin, alias = "", description = NULL, xyz.margin = c (0, 0, 0)) {
  
  if (!is (vol, "volume")) {
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  if (!is (sel.bin, "volume")) {
    warning ("sel.bin should be a volume class object.")
    return (NULL)
  }
  if ((sel.bin$modality!="binary")) {
    warning ("vol must be of binary modality.")
    return (NULL)
  }
  
  if(is.null(vol$vol3D.data)){
    warning ("empty vol$vol3D.data.")
    return (NULL)
  }
  if(is.null(sel.bin$vol3D.data)){
    warning ("empty sel.bin$vol3D.data.")
    return (NULL)
  }
  
  #verifier que les volumes ont le même support
  if (suppressWarnings(!grid.equal (vol, sel.bin))) {
    warning ("vol and sel.bin volumes must share the same grid.")
    return (NULL)
  }
  
  t.mat <- ref.cutplane.add(vol, T.MAT = NULL)
  vol_ <- vol.in.new.ref(vol, new.ref.pseudo=paste0(vol$ref.pseudo,"m"), t.mat)
  bin_ <- vol.in.new.ref(sel.bin, new.ref.pseudo=paste0(vol$ref.pseudo,"m"), t.mat)
  
  ijk.margin <- abs(xyz.margin / bin_$dxyz)
  ijk.margin[is.infinite(ijk.margin)] <- 0
  ijk.margin <- ceiling(ijk.margin)
  if (is.null(description)) description <-  paste (vol$description, "restricted to", 
                                                   sel.bin$description)
  new.vol <- vol.copy (vol_, alias = alias, description = description)
  range.ijk <- apply(get.ijk.from.index(which(bin_$vol3D.data==TRUE), bin_),2,range)
  new.vol <- .vol.border.tuning (new.vol, pre.nijk = -as.numeric(range.ijk[1,]) + ijk.margin, 
                                 post.nijk = -c(vol_$n.ijk[1:2]-1, vol_$k.idx[sel.bin$n.ijk[3]]) +  
                                   as.numeric(range.ijk[2,]) + ijk.margin)
  
  return (vol.in.new.ref(new.vol,new.ref.pseudo=vol$ref.pseudo, t.mat,alias = alias, 
                         description=description))
}
  