#' Restrict volume to a binary selection
#' @description The \code{nesting.bin} function restricts a "volume" class 
#' object to the rectangular parallelepiped circumscribed to the selected voxels.
#' @param obj "volume" class object, containing data to restrict.
#' @param sel.bin "volume" class object, of "binary" modality, specifying the selected voxels.
#' @param xyz.margin Numeric vector of length 3, by default set to \code{c (0, 0, 0)}.
#' See details.
#' @param obj.restrict Boolean. Used if \code{obj} is of class "volume". If 
#' \code{obj.restrict = TRUE}, the rectangular parallelepiped circumscribed to 
#' the selected voxels, enlarged by xyz.margin cannot exceed the initial volume.
#' @param alias Character string, \code{$alias} of the created object.
#' @param description Character string, describing the the created object. 
#' If \code{description = NULL}, it will be 
#' \code{paste (obj$description,"restricted to", sel.bin$description)}.
#' @param ... Argument such as T.MAT, or for deprecated arguments
#' @return Returns a "volume" class object, in which 3D volume is limited to the
#' rectangular parallelepiped circumscribed to the voxels selected by \code{sel.bin}, increased by the
#' requested margins.

#' @details If \code{obj} is of class "volume", \code{xyz.margin} represents the 
#' distances in mm to be added to the x, y and z directions of the rectangular 
#' parallelepiped circumscribed to the voxels selected in \code{sel.bin}, in the 
#' \code{obj} frame of reference.
#' 
#' @details If \code{obj} is of class “mesh”, \code{sel.bin} will undergo a 
#' margin expansion \code{xyz.margin} before the mesh points are selected..

#' @seealso \link[espadon]{add.margin}, \link[espadon]{nesting.cube} and  
#' \link[espadon]{nesting.roi}.
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
nesting.bin <- function (obj, sel.bin, alias = "", description = NULL, 
                         xyz.margin = c (0, 0, 0), obj.restrict = FALSE,...) {
  
  passed <- names(as.list(match.call())[-1])
  args <- list(...)
  if (!("obj" %in% passed)){
    if (is.null(args[['vol']])) stop('argument "obj" is missing, with no default')
    obj <- args[['vol']]
  }

  if (!is.null(args[['vol.restrict']])) obj.restrict <- args[['vol.restrict']]
  eps <- args[["eps"]]
  if (is.null(eps)) eps <- 1e-9
  
  T.MAT <- args[['T.MAT']]
  
  if (!is (sel.bin, "volume")) stop ("sel.bin should be a volume class object.")
  if (sel.bin$modality!="binary") stop ("sel.bin must be of binary modality.")
  if (is.null(sel.bin$vol3D.data)) stop ("empty sel.bin$vol3D.data.")
  
  if (is.null(description)) description <-  paste (obj$description, "restricted to", 
                                                   sel.bin$description)
  
  alias.sel <- sel.bin$object.alias
  sel.bin <- vol.in.new.ref(sel.bin,obj$ref.pseudo, T.MAT)
  sel.bin$object.alias <- alias.sel
  
  if (is (obj, "volume")) {
  #-------------------------
    if(is.null(obj$vol3D.data)) stop ("empty obj$vol3D.data.")
  
    ext.pt <- get.extreme.pt(sel.bin, pt.min = eps)
    ext.pt[1,] <- ext.pt[1,] + abs(xyz.margin[1])*c(-1,1)
    ext.pt[2,] <- ext.pt[2,] + abs(xyz.margin[2])*c(-1,1)
    ext.pt[3,] <- ext.pt[3,] + abs(xyz.margin[3])*c(-1,1)
    
    ext.pt.obj <- get.extreme.pt(obj)
    if (obj.restrict) {
      ext.pt[1,] <- .range.intersection(ext.pt[1,],ext.pt.obj[1,])
      ext.pt[2,] <- .range.intersection(ext.pt[2,],ext.pt.obj[2,])
      ext.pt[3,] <- .range.intersection(ext.pt[3,],ext.pt.obj[3,])
    }
    if (any(is.na(ext.pt))) stop("obj and the selection of sel.bin have no common zones")
    obj_ <- nesting.cube(obj,pt.min = ext.pt[,1], pt.max = ext.pt[,2],
                         alias=alias, description = description)
    # #verifier que les volumes ont le même support
    # if (suppressWarnings(!grid.equal (obj, sel.bin))) {
    #   warning ("obj and sel.bin volumes must share the same grid.")
    #   return (NULL)
    # }
    # 
    # t.mat <- ref.cutplane.add(obj, ref.cutplane="rcp", T.MAT = NULL)
    # vol_ <- vol.in.new.ref(obj, new.ref.pseudo="rcp", t.mat)
    # bin_ <- vol.in.new.ref(sel.bin, new.ref.pseudo="rcp", t.mat)
    # 
    # xyz.margin_ <- c(xyz.margin, 0) %*% t(get.rigid.M(t.mat, obj$ref.pseudo, "rcp"))
    # 
    # idx.c <- which(apply(abs(bin_$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
    # idx.r <-  which(apply(abs(bin_$xyz.from.ijk[1:3,1:3]),1,sum)==0)
    # if (length(idx.c)>0) {
    #   #2D
    #   u <- bin_$xyz.from.ijk 
    #   u[idx.r,idx.c]<- 1
    #   ijk.from.xyz <- solve(u)
    #   ijk.from.xyz[idx.r,idx.c] <- 0
    #   ijk.margin <- xyz.margin_ %*% t(ijk.from.xyz)
    # } else {
    #   ijk.margin <- xyz.margin_ %*% t(solve(bin_$xyz.from.ijk))
    # }
    # ijk.margin[1:3][bin_$n.ijk < 2] <- 0
    # ijk.margin <- ceiling(abs(ijk.margin)[1:3])
    # 
    # new.vol <- vol.copy (vol_, alias = alias, description = description)
    # range.ijk <- apply(get.ijk.from.index(which(bin_$vol3D.data==TRUE), bin_),2,range)
    # 
    # range.ijk[1,] <- range.ijk[1,]-ijk.margin
    # range.ijk[2,] <- range.ijk[2,]+ijk.margin
    # if (vol.restrict){
    #   range.ijk[1,1]<- max(range.ijk[1,1],0)
    #   range.ijk[1,2]<- max(range.ijk[1,2],0)
    #   range.ijk[1,3]<- max(range.ijk[1,3],vol_$k.idx[1])
    #   range.ijk[2,1]<- min(range.ijk[2,1],vol_$n.ijk[1]-1)
    #   range.ijk[2,2]<- min(range.ijk[2,2],vol_$n.ijk[2]-1)
    #   range.ijk[2,3]<- min(range.ijk[2,3],vol_$k.idx[vol_$n.ijk[3]])
    # }
    # new.vol <- .vol.border.tuning (new.vol, pre.nijk = -as.numeric(range.ijk[1,]), 
    #                                post.nijk = -c(vol_$n.ijk[1:2]-1, vol_$k.idx[sel.bin$n.ijk[3]]) +  
    #                                  as.numeric(range.ijk[2,]))
    # new.vol <- vol.in.new.ref(new.vol,new.ref.pseudo=vol$ref.pseudo, t.mat,alias = alias, 
    #                           description=description) 
    # 
  } else if (is (obj, "mesh")){
    
    #on aggrandit sel.bin de la marge: 
  
    if (any(abs(xyz.margin)!=0)) sel.bin <- bin.dilation(sel.bin, radius =abs(xyz.margin))
    sel.bin$object.alias <- alias.sel
    
    obj_ <- obj
    obj_$description <- description
    obj_$object.alias <- alias
    if (obj_$nb.faces==0)  return(obj_)
    
    if ((obj$ref.pseudo != sel.bin$ref.pseudo) & is.null(T.MAT)) 
      stop("obj and sel.bin have different ref.pseudo. Set T.MAT argument")
    f <- get.value.from.xyz(t(obj_$mesh$vb[1:3,]), sel.bin, T.MAT = T.MAT)>=0.5
    f[is.na(f)] <- FALSE
    to.suppress <- which(!f)
    obj_$mesh$normals <- obj_$mesh$normals[,f]
    obj_$mesh$vb <- obj_$mesh$vb[,f]
    to.keep <- which(f)
    face.keep <- is.na(match(obj_$mesh$it[1,],to.suppress)) & 
      is.na(match(obj_$mesh$it[2,],to.suppress)) & 
      is.na(match(obj_$mesh$it[3,],to.suppress))
    obj_$mesh$it <- obj_$mesh$it[ ,face.keep]
    obj_$mesh$it[1,] <- match(obj_$mesh$it[1,],to.keep)
    obj_$mesh$it[2,] <- match(obj_$mesh$it[2,],to.keep)
    obj_$mesh$it[3,] <- match(obj_$mesh$it[3,],to.keep)
    if (!is.null(obj_$mesh$remface)) obj_$mesh$remface<- obj$mesh$remface[face.keep]
   
  }
  if (alias=="") return(obj_)
  return(.set.ref.obj(obj_,list(obj,sel.bin), add=FALSE))
}
  