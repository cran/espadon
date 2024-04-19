#' Display the transverse frontal or sagittal view in the patient reference system
#' @description The \code{display.plane} function displays an overlay of images and RoI
#' closed planar contours on a plane defined by the equations x = constant (sagittal 
#' view), or y = constant (frontal view) or z = constant (transverse view) in a 
#' frame of reference chosen by the user.
#' @param bottom "volume" class object, displayed using \code{bottom.col}
#' palette. If \code{bottom = NULL}, no bottom image is displayed.
#' @param top "volume" class object, displayed as an overlay, using \code{top.col}
#' palette. If \code{top = NULL}, no overlay image is displayed.
#' @param struct "struct" class object. If \code{NULL}, no RoI is displayed. Only
#' RoI of closed planar or point type are displayed.
#' @param roi.name Vector of exact names of the RoI in the \code{struct} object.
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the \code{struct} 
#' object. By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @param struct.dxyz 3D vector. Used in case of \code{bottom} and 
#' \code{top} are set to \code{NULL}. It represents the voxel size in the \code{display.ref} 
#' frame of reference, used to calculate contours in frontal or sagittal view.
#' @param display.ref Character string. Pseudonym of the frame of reference used 
#' for display. If \code{NULL} (default), the bottom image FoR, or top image FoR 
#' (when no bottom image), or struct FoR (when no volume displayed).
#' @param T.MAT "t.mat" class object, created by \link[espadon]{load.patient.from.Rdcm} 
#' or \link[espadon]{load.T.MAT}. If \code{T.MAT} is \code{NULL}, \code{bottom}, 
#' \code{top} and \code{struct} must have the same frame of reference.
#' @param interpolate Boolean, indicating whether to apply linear interpolation, 
#' when calculating the bottom and top cuts,and then when displaying them.
#' If \code{interpolate = FALSE}, the values of the nearest voxels are used. When \code{TRUE} (by delfault), 
#' trilinear interpolation is used.
#' @param view.type Character string, defining the view to display. It must be set to
#' \itemize{
#' \item \code{"trans"} for a transverse view,
#' \item \code{"front"} for a frontal view or,
#' \item \code{"sagi"} for a sagittal view.
#' }
#' @param view.coord Numeric vector of the coordinates along the normal vector of 
#' the selected view.
#' @param bg Background color of the image. By default, this color is black.
#' @param abs.rng Vector of 2 elements indicating the minimum and maximum abscissa
#' to display on the background image.
#' @param ord.rng Vector of 2 elements indicating the minimum and maximum ordinate
#' to display on the background image.
#' @param bottom.col,top.col Vectors, representing the palette color of 
#' \code{bottom} and \code{top}.
#' @param bottom.breaks,top.breaks One of :
#' \itemize{
#' \item \code{NULL} : the minimum and the maximum value of \code{bottom} or
#' \code{top} define the range.
#' \item Vector giving the breakpoints of each color. Outside values are transparent,
#' leaving the background visible, depending on \code{sat.transp}.
#' }
#' When breaks are specified, the number of breaks must be one unit more then the number of colors.
#' @param sat.transp Boolean. If \code{TRUE}, outside values are transparent, else set
#' to \code{bottom.breaks} or \code{top.breaks} limits.
#' @param struct.lwd Line thickness of the RoI contours.
#' @param main Character string. When \code{main} different from \code{NULL}, 
#' it replaces the title, and removes the subtitle and the maximum dose indication 
#' if \code{top} is of modality rtdose.
#' @param legend.plot Boolean, that indicates whether the RoI legend should be 
#' displayed on the image. It is displayed by default.
#' @param legend.shift Numeric. It shifts (in mm) the display of the RoI legend 
#' on x-axis.
#' @param legend.roi.pseudo Boolean. If \code{TRUE}, the name used 
#' for a RoI in the legend comes from the \code{struct$roi.info$roi.pseudo} 
#' column, otherwise the \code{struct$roi.info$name} column.
#' @param ... others parameters of plot function
#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are
#' all set to \code{NULL}, all closed planar or point RoI are selected.
#' If a RoI is not present in the requested plane, the RoI legend won't mention it.
#' @note 1- The main title is given by \code{bottom}, the
#' subtitle by \code{top}.
#' @note 2- When \code{top} is in the "rtdose" modality, the maximum dose is 
#' written on the image.
#' @seealso \link[espadon]{display.kplane}, \link[espadon]{plot.volume}, 
#' \link[espadon]{plot.struct}, \link[espadon]{plot.mesh}.
#' @return Returns a display of the  transverse, sagittal or frontal plane. This plane
#' has the coordinate z = view.coord (transverse), y = view.coord (sagittal) pr 
#' x = view.coord (frontal). The display is an overlay of:
#' \itemize{
#' \item a background image of uniform color \code{bg}
#' \item the bottom image if it exists
#' \item the top image if it exists
#' \item the contours of the regions of interest if they exist in the plane considered.
#' } 
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "mr", "rtstruct", "rtdose"), 
#'                              roi.name  = "",
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' CT <- patient$ct[[1]]
#' MR <- patient$mr[[1]]
#' D <- patient$rtdose[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' display.plane (bottom = CT, top = D, struct = S, view.coord = -30, 
#'                interpolate = FALSE, legend.shift = -80)
#' # Display of CT in reference frame "ref1" and  MR in "ref2"               
#' display.plane (bottom = CT, top = MR, interpolate = FALSE)
#' 
#' # Display of CT and MR in reference frame "ref2"
#' display.plane (bottom = CT, top = MR, interpolate = FALSE, display.ref ="ref2",
#'                T.MAT = patient$T.MAT)
#' @export
#' @importFrom grDevices rainbow grey.colors
#' @importFrom methods is
display.plane <- function (bottom = NULL, top = NULL, struct = NULL, 
                           roi.name = NULL, roi.sname = NULL, roi.idx = NULL, 
                           struct.dxyz = c (0.5, 0.5, struct$thickness), 
                           display.ref = NULL, 
                           T.MAT = NULL, interpolate = TRUE,
                           view.type = c("trans", "front", "sagi"), 
                           view.coord = 0,
                           bg="#000000", abs.rng = NULL, ord.rng = NULL,
                           bottom.col = grey.colors (255, start = 0, end = 1),
                           top.col = pal.rainbow (255),
                           bottom.breaks = NULL, top.breaks = NULL, 
                           sat.transp = FALSE,
                           struct.lwd=2, main = NULL, 
                           legend.plot = TRUE, legend.shift = 0,
                           legend.roi.pseudo = TRUE,...) {
  
  args <- tryCatch(list(...), error = function(e)list())
  if(!is.null(abs.rng)) args[["xlim"]] <- abs.rng
  if(!is.null(ord.rng)) args[["ylim"]] <- ord.rng
  if(!is.null(main)) args[["main"]] <- main
  args[["bg"]] <- bg
  args[["view.type"]] <-view.type[1]
  add <- FALSE
  if(!is.null(args[["add"]]))  add <- args[["add"]]
    
  xpd <- NULL
  on.exit(
    expr = {
      if (!is.null(xpd)) par(xpd = xpd)
    })
  
  view.type <- view.type[1]
  list.roi.idx <- NULL
  
  if (!is.null(struct) & !is (struct, "struct")) stop ("struct should be a struct class object.")

  if (!is.null(bottom)) {
    if (!is (bottom, "volume")){
      stop ("bottom should be a volume class object.")
    } else if  (is.null(bottom$vol3D.data)){
        message ("bottom should have vol3D.data.")
        bottom <- NULL
    }
  }
  if (!is.null(top)) {
    if (!is (top, "volume")) {
      stop ("top should be a volume class object.")
    } else if  (is.null(top$vol3D.data)){
      message ("top should have vol3D.data.")
      top <- NULL
    }
  }
    
  if (length(view.coord)==0) stop ("view.coord length is 0.")
  if (!is.null (struct)) list.roi.idx <- select.names (struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  
  warn.ref <- FALSE
  warn.ref.struct <- FALSE
  selected.ref <- display.ref
  if (is.null(selected.ref)) {
    if (!is.null(bottom)) selected.ref <- bottom$ref.pseudo
    else if (!is.null(top)) selected.ref <- top$ref.pseudo
    else if (!is.null(list.roi.idx) & length(list.roi.idx)!=0) selected.ref <- struct$ref.pseudo
    else return(NULL)
  }
  
  if (length(unique(c(display.ref,bottom$ref.pseudo, 
                      top$ref.pseudo, struct$ref.pseudo)))!=1 & 
      is.null(T.MAT)) warning("objects have different ref.pseudo. Load T.MAT for correct display")
  if (!is.null(bottom)) {
    dum <- suppressWarnings(vol.in.new.ref(bottom, selected.ref, T.MAT))
    if (is.null(dum)) {
      warning(paste("bottom is displayed in the ref.pseudo", bottom$ref.pseudo, "instead of", selected.ref))
      warn.ref <- TRUE
    } else {bottom <- dum}}
  if (!is.null(top)) {
    dum <- suppressWarnings(vol.in.new.ref(top, selected.ref, T.MAT))
    if (is.null(dum)) {
      warning(paste("top is displayed in the ref.pseudo", top$ref.pseudo, "instead of", selected.ref))
      warn.ref <- TRUE
    } else {top <- dum}}
  
  if (is.null(bottom) & is.null(top) & is.null(list.roi.idx)) {
    stop ("nothing to display")
  }
  
  
  if (!is.null(list.roi.idx)) {
    dum <- suppressWarnings(struct.in.new.ref (struct,new.ref.pseudo= selected.ref, T.MAT))
    if (is.null(dum)) {
      warning(paste("struct is displayed in the ref.pseudo", struct$ref.pseudo, "instead of", selected.ref))
      warn.ref.struct <- TRUE
    } else {struct <- dum}}
  if (is.null (struct))  {
    list.roi.idx <- NULL
  }else{
    back.dxyz <- struct.dxyz
    if (!is.null(bottom)) {
      back.dxyz[1:2] <- c(min(bottom$dxyz[1],back.dxyz[1]), min(bottom$dxyz[2],back.dxyz[2]))
    } else if (!is.null(top)) {
      back.dxyz <- c(min(bottom$dxyz[1],back.dxyz[1]), min(bottom$dxyz[2],back.dxyz[2]))
    }
  }
  
 
  
  #############################################################################
  for (coord.idx in 1:length (view.coord)){
    args[["view.coord"]] <- view.coord[coord.idx]
    
    coord.lab <- switch (view.type, "trans" = "z", "front"="y", "sagi" = "x")
    
    top.p <-bottom.p <- NULL
    if (! is.null(bottom)){
      args_ <- args
 
      args_[["x"]] <- bottom
      args_[["col"]] <- bottom.col
      args_[["breaks"]] <- bottom.breaks 
      args_[["sat.transp"]] <- sat.transp 
      args_[["add"]] <-  add
      args_[["cut.interpolate"]] <- interpolate
      args_[["display.interpolate"]] <- interpolate
      bottom.p <- suppressMessages(do.call(plot,args_))
      if (is.na(bottom.p$max.pixel)) {
        bottom.p <- NULL; 
        message("no bottom view @", coord.lab," = ", round(view.coord[coord.idx],2), " mm")
        }
    }
    
    if (!is.null(top)) {
      args_ <- args
      args_[["x"]] <- top
      args_[["col"]] <- top.col
      args_[["breaks"]] <- top.breaks 
      args_[["sat.transp"]] <- sat.transp 
      args_[["add"]] <- !is.null(bottom.p$max.pixel) | add

      top.p <- suppressMessages(do.call(plot,args_))
      if (is.na(top.p$max.pixel)) {
        top.p <- NULL; 
        message("no top view @", coord.lab," = ", round(view.coord[coord.idx],2), " mm")
      } else if (is.null (main)) {
        if (args_[["add"]]){
          idx <- which(top.p$xyz.from.ijk[,3]!=0)
          mtext (paste (top$modality, " (", top$description,") @ ", c("x","y","z")[idx],
                        " = ",round (top.p$xyz0[1,idx],3)," mm",sep=""),
                 side=3, line=0.4, col='gray32', cex=0.8)
        }
        if (top$modality =="rtdose")
          text (par("usr")[1], par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1,
                paste("  Dose max : ",round (top.p$max.pixel, 3)," Gy",sep=""), cex=1, col="red",adj = c(0,0))
      }
      
    }
    if (length(list.roi.idx)!=0){
      args_ <- args
      args_[["x"]] <- struct
      args_[["col"]] <- NULL
      args_[["breaks"]] <- NULL
      args_[["sat.transp"]] <- NULL 
      args_[["add"]] <- !(is.null(bottom.p) & is.null(top.p)) | add
      args_[["lwd"]] <- !is.null(struct.lwd)
      args_[["interpolate"]] <- NULL
      args_[["roi.idx"]] <- list.roi.idx
      args_[["lwd"]] <- struct.lwd
      args_[["back.dxyz"]] <- 
      S  <- do.call(plot,args_)
      
      if (S$nb.of.roi>0){
        legendlabel <- S$roi.info$roi.pseudo
        legendcol <- S$roi.info$color
        type <- do.call(rbind.data.frame,lapply(S$roi.data, function(L) {
          v <- castlow.str(unique(sapply(L,function(l) l$type)))
          return(c(any(grepl("planar$",v)),any(grepl("point$",v))))}))
        legendlty <- rep(0,nrow(type)); legendlty[type[,1]] <- 1
        legendpch <- rep(" ",nrow(type)); legendpch[type[,2]] <- "+"
        if (length (legendlabel)>0 & legend.plot) {
          xpd <- par()$xpd
          par(xpd=TRUE)
          legend(par("usr")[2]+legend.shift ,par("usr")[4],
                 legend = legendlabel, col = legendcol,
                 ncol=1, lty =legendlty ,lwd=struct.lwd, pch = legendpch, bty="o", 
                 cex=0.6, text.col="white",bg="black")
        }
      }
    }
  }
  ##############################################################################
  # if ((is.null(bottom) & is.null(top)) | (warn.ref.struct)){
  #   #on construit un support pour les contours
  #   rng.x <- c (floor (min(struct$roi.info[list.roi.idx,]$min.x)), max(struct$roi.info[list.roi.idx,]$max.x))
  #   rng.y <- c (floor (min(struct$roi.info[list.roi.idx,]$min.y)), max(struct$roi.info[list.roi.idx,]$max.y))
  #   rng.z <- c (floor (min(struct$roi.info[list.roi.idx,]$min.z)), max(struct$roi.info[list.roi.idx,]$max.z))
  #   nxyz <- c(ceiling((rng.x[2] - rng.x[1])/struct.dxyz[1])+11,
  #             ceiling((rng.y[2] - rng.y[1])/struct.dxyz[2])+11,
  #             ceiling((rng.z[2] - rng.z[1])/struct.dxyz[3])+11)
  #   struct.vol3D <- vol.create (n.ijk =nxyz, pt000= c(rng.x[1]-5*struct.dxyz[1],
  #                                                     rng.y[1]-5*struct.dxyz[2],
  #                                                     rng.z[1]-5*struct.dxyz[3]),
  #                               dxyz = struct.dxyz,
  #                               ref.pseudo = struct$ref.pseudo,
  #                               frame.of.reference = struct$frame.of.reference,
  #                               alias = struct$object.alias, number = 0,
  #                               modality = struct$modality,  description = "")
  #   
  #   if (!warn.ref.struct) {
  #     struct.vol3D <- vol.in.new.ref (struct.vol3D, selected.ref, T.MAT)
  #   }
  # }
  # 
  # #centre image
  # if (!is.null(bottom)) {center.pt <- apply (get.extreme.pt (bottom),1,mean)
  # } else if (!is.null(top)) {center.pt <- apply (get.extreme.pt (top),1,mean)
  # } else center.pt <- apply (get.extreme.pt (struct.vol3D),1,mean)
  # 
  # 
  # lab <- c("x", "y", "z")
  # 
  # 
  # 
  # 
  # for (coord.idx in 1:length (view.coord)) {
  # 
  #   if (view.type=="sagi") {
  #     plane.orientation= c (0, 0, 1, 0, 1, 0, 1, 0, 0)
  #     lab.idx <- c(3,2,1)
  #     ord.flip <- TRUE
  #     w.idx <- 1
  #   } else if (view.type=="front") {
  #     plane.orientation= c (1, 0, 0, 0, 0, 1, 0, 1, 0)
  #     lab.idx <- c(1,3,2)
  #     ord.flip <- FALSE
  #     w.idx <- 2
  #   } else {
  #     plane.orientation= c(1, 0, 0, 0, 1, 0, 0, 0, 1)
  #     lab.idx <- c(1,2,3)
  #     ord.flip <- TRUE
  #     w.idx <- 3
  #   }
  #   plane.pt <- center.pt
  #   p.idx <- (1:3)[-w.idx]
  #   plane.pt[w.idx] <- view.coord[coord.idx]
  #   
  #   process.ori <- function(pt, vol,p.idx){
  #     # if(all(apply (abs(vol$xyz.from.ijk[1:3, p.idx])< 1e-4 ,2,sum)==2)){
  #     #   center.ijk <- get.ijk.from.xyz(pt, vol)
  #     #   center.ijk[p.idx] <- round(center.ijk[p.idx])
  #     #   pt <- (vol$xyz.from.ijk %*% c(center.ijk,1))[1:3]
  #     # }
  #     pt
  #   }
  #   
  #   if (!is.null(bottom)){
  #     #bottom.p <-get.plane(bottom, origin = plane.pt, plane.orientation= plane.orientation, rev.k=rev.k, interpolate =interpolate)
  #     bottom.p <- get.plane(bottom, origin = process.ori (plane.pt, bottom, p.idx), 
  #                           plane.orientation= plane.orientation,
  #                           interpolate = interpolate)
  #     if (!is.null (bottom.p)){
  #       pt000 <- c(0, 0, 0, 1) %*% t(bottom.p$xyz.from.ijk)
  #       if (is.null(bottom.breaks)){
  #         b <- .pixel.scale (bottom$min.pixel,bottom$max.pixel,length(bottom.col))
  #       } else { b <- bottom.breaks}
  #       
  #       if (is.null (main)){ 
  #         main.title <- paste (bottom$modality, " (",bottom$description,") @ ",
  #                              lab[lab.idx[3]], " = ",round (pt000[lab.idx[3]],3)," mm",sep="")
  #       } else {
  #         main.title <- main
  #       }
  #       display.kplane (vol = bottom.p, pt00= pt000[lab.idx[1:2]], dxy= bottom.p$dxyz[1:2],
  #                       col = bottom.col, breaks = b, sat.transp = sat.transp,
  #                       abs.lab = lab [lab.idx[1]],
  #                       ord.lab = lab [lab.idx[2]], ord.flip = ord.flip,
  #                       main = main.title,
  #                       bg=bg, abs.rng = abs.rng, ord.rng = ord.rng, interpolate=interpolate)
  #       if (warn.ref | warn.ref.struct) mtext ("warning : different frames of reference",side=1, line=2, col='red', cex=0.8)
  #     }
  #     if (!is.null (top)) {
  #       plane.pt[w.idx] <- pt000[w.idx]
  #       #top.p <-get.plane(top, origin = plane.pt, plane.orientation= plane.orientation, rev.k=rev.k, interpolate =interpolate)
  #       top.p <- get.plane(top, origin = process.ori (plane.pt, top,  p.idx), 
  #                          plane.orientation= plane.orientation,
  #                           interpolate = interpolate)
  #       if (!is.null (top.p)){
  #         pt000 <- c(0, 0, 0, 1) %*% t(top.p$xyz.from.ijk)
  #         if (is.null(top.breaks)){
  #           b <- .pixel.scale (top$min.pixel,top$max.pixel,length(top.col))
  #         } else { b <- top.breaks}
  #         display.kplane (vol=top.p, pt00= pt000[lab.idx[1:2]], dxy= top.p$dxyz[1:2],
  #                          col = top.col, breaks = b, sat.transp = sat.transp,
  #                          add=TRUE,  interpolate=interpolate)
  #         if (is.null (main)) {
  #           mtext (paste (top$modality, " (", top$description,") @ ", lab[lab.idx[3]], 
  #                         " = ",round (pt000[lab.idx[3]],3)," mm",sep=""),
  #                  side=3, line=0.4, col='gray32', cex=0.8)
  #           if (top$modality =="rtdose")
  #             text (par("usr")[1], par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1, 
  #                   paste("  Dose max : ",round (top.p$max.pixel, 3)," Gy",sep=""), cex=1, col="red",adj = c(0,0))
  #         }
  #       }
  #     }
  #   } else if (!is.null(top)){
  #     bottom.p <- get.plane(top, origin = process.ori (plane.pt, top, p.idx), 
  #                           plane.orientation= plane.orientation,
  #                           interpolate = interpolate)
  #     if (!is.null (bottom.p)){
  #       pt000 <- c(0, 0, 0, 1) %*% t(bottom.p$xyz.from.ijk)
  #       if (is.null(top.breaks)){
  #         b <- .pixel.scale (top$min.pixel,top$max.pixel,length(top.col))
  #       } else { b <- top.breaks}
  #       
  #       if (is.null (main)){ 
  #         main.title <- paste (top$modality, " (",top$description,") @ ", lab[lab.idx[3]],
  #                              " = ",round (pt000[lab.idx[3]],3)," mm",sep="")
  #       } else {
  #         main.title <- main
  #       }
  #       
  #       display.kplane (vol=bottom.p, pt00= pt000[lab.idx[1:2]], dxy= bottom.p$dxyz[1:2],
  #                        col = top.col, breaks = b, sat.transp = sat.transp,
  #                        abs.lab = lab [lab.idx[1]], ord.lab = lab [lab.idx[2]], ord.flip = ord.flip,
  #                        main = main.title,
  #                        bg=bg, abs.rng = abs.rng, ord.rng = ord.rng, interpolate=interpolate)
  #       if (is.null (main) & (top$modality =="rtdose"))
  #         text (par("usr")[1], par("usr")[4] - (par("usr")[4]-par("usr")[3])*0.1, paste("  Dose max : ",round (bottom.p$max.pixel, 3)," Gy",sep=""), cex=1, col="red",adj = c(0,0))
  #     
  #       if (warn.ref | warn.ref.struct) mtext ("warning : different frames of reference",side=1, line=2, col='red', cex=0.8)
  #       }
  #     
  #   } else {
  #     bottom.p <- get.plane(struct.vol3D, origin = process.ori (plane.pt, struct.vol3D,  p.idx),
  #                            plane.orientation= plane.orientation)
  #     
  #     if (!is.null (bottom.p)){
  #       pt000 <- c(0, 0, 0, 1) %*% t(bottom.p$xyz.from.ijk)
  #       if (is.null (main)){ 
  #         main.title <- paste (struct$modality, " (",top$description,") @ ", 
  #                              lab[lab.idx[3]], " = ",round (pt000[lab.idx[3]],3)," mm",sep="")
  #       } else {
  #         main.title <- main
  #       }
  #       display.kplane (vol=bottom.p, pt00= pt000[lab.idx[1:2]], dxy= bottom.p$dxyz[1:2],
  #                        abs.lab = lab [lab.idx[1]],
  #                        ord.lab = lab [lab.idx[2]], ord.flip = ord.flip,
  #                        main =  main.title,
  #                        bg=bg, abs.rng = abs.rng, ord.rng = ord.rng, interpolate=interpolate)
  #       if (warn.ref | warn.ref.struct) mtext ("warning : different frames of reference",side=1, line=2, col='red', cex=0.8)
  #     }
  #   }
  #   
  #   if (length(list.roi.idx)>0) {
  #     legendcol <- list()
  #     legendlabel <- list()
  #     legendlty <- list()
  #     legendpch <- list()
  #     label.index <- 1
  #     if (legend.roi.pseudo) {legend.name <- struct$roi.info$roi.pseudo} else {legend.name <- struct$roi.info$name}
  #     if (lab[w.idx]=="z" & 
  #         (bottom.p$ref.pseudo==struct$ref.pseudo | warn.ref) & 
  #         all(round(as.numeric(struct$ref.from.contour),6)== as.numeric(diag(4)))) {
  #       new.struct <-  .display.select.struct.by.z (struct=struct, list.roi.idx= list.roi.idx, z =bottom.p$xyz0[1,w.idx], dz = struct$thickness)
  #       
  #     } else {
  # 
  #       if (warn.ref.struct) {
  #         bottom.p <-  get.plane(struct.vol3D, origin = process.ori (plane.pt, struct.vol3D,  p.idx),
  #                                plane.orientation= plane.orientation)
  #       }
  #    
  #       t.mat <- ref.cutplane.add(bottom.p, ref.cutplane = "intern", origin = c(0,0,0))
  #       
  #       new.struct <- lapply(1:struct$nb.of.roi, function(r.idx){
  #         if (!(r.idx %in% list.roi.idx)) return (NULL)
  #         if ((length(struct$roi.data[[r.idx]]) ==1) &   
  #             (castlow.str (struct$roi.data[[r.idx]][[1]]$type) =="point")) return(struct$roi.data[[r.idx]])
  #         roi.nesting <- suppressWarnings (nesting.roi (obj=bottom.p, struct=struct, roi.idx=r.idx, 
  #                                     T.MAT=T.MAT, xyz.margin=c(1,1,1), vol.restrict=TRUE))
  #         if (is.null(roi.nesting)) return (NULL)
  #         bin <-bin.from.roi (vol=roi.nesting, struct=struct, roi.idx=r.idx, T.MAT=T.MAT)
  #         bin_ <- vol.in.new.ref(bin, new.ref.pseudo="intern", t.mat)
  #        return (.display.roi.data.from.bin (bin_))})
  #       names(new.struct) <- struct$roi.info$roi.pseudo
  #     }
  #     for (j in list.roi.idx) {
  #       if (length(new.struct[[j]])>0) {
  #         for (nb in 1:length(new.struct[[j]])){
  #           type <- castlow.str (new.struct[[j]][[nb]]$type)
  #           test.pt <-FALSE 
  #           if (type=="closedplanar" | type=="openplanar"){
  #             test.pt <- TRUE
  #             lines (new.struct[[j]][[nb]]$pt$x, new.struct[[j]][[nb]]$pt$y, col= struct$roi.info$color[j], lwd=struct.lwd)
  #             legendlty[[label.index]]<-1
  #             legendpch[[label.index]]<-" "
  #           } else if ((type=="point") &
  #             (round(new.struct[[j]][[nb]]$pt[1,w.idx],6) == round(view.coord[coord.idx],6))) {
  #               test.pt <- TRUE
  #               points(new.struct[[j]][[nb]]$pt[lab.idx[1]], new.struct[[j]][[nb]]$pt[lab.idx[2]], col= struct$roi.info$color[j], pch="+",cex=1)
  #               legendlty[[label.index]]<-0
  #               legendpch[[label.index]]<-"+"
  #           }
  #         }
  #         if (test.pt){
  #           legendcol[[label.index]]<- struct$roi.info$color[j]
  #           legendlabel[[label.index]]<- legend.name[j]
  #           label.index<-label.index+1
  #         }
  #       }
  #     }
  #     if (length (legendlabel)>0 && legend.plot) {
  #       
  #       # 	par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4]
  #       xpd <- par()$xpd
  #       par(xpd=TRUE)
  #       
  #       legend(par("usr")[2]+legend.shift ,par("usr")[4],
  #              legend = unlist (legendlabel), col = unlist (legendcol),
  #              ncol=1, lty = unlist (legendlty) ,lwd=struct.lwd, pch = unlist (legendpch), bty="o", cex=0.6, text.col="white",bg="black")
  #       
  # 
  #     }
  #   }
  # 
  # }
  
}