#' Merging of structures into a new structure
#' @description The \code{struct.merge} function merges two structures 
#' into a new one. It is useful for comparing contours, for example.
#' @param ref.struct struct class object. All RoI of this structure are kept.
#' @param add.struct struct class object. Only the selected RoI are kept for merging.
#' @param roi.name Vector of exact names of the RoI in the \code{add.struct} object. 
#' By default \code{roi.name = NULL}. See Details.
#' @param roi.sname Vector of names or parts of names of the RoI in the 
#' \code{add.struct} object. By default \code{roi.sname = NULL}. See Details.
#' @param roi.idx Vector of indices of the RoI that belong to the \code{add.struct} 
#' object. By default \code{roi.idx = NULL}. See Details.
#' @param suffix Character string. \code{'-suffix'} is added to RoI name.
#' @param alias Character string, \code{$alias} of the resulted object.
#' @param description Character string, describing the the resulted object.

#' @details If \code{roi.name}, \code{roi.sname}, and \code{roi.idx} are all 
#' \code{NULL}, then all RoI of add.struct are selected.
#' @note Beware that, when merging structures, some RoI may have same name or 
#' \code{roi.info$roi.pseudo}. In this case \code{struct.merge} prints a warning 
#' message. Consider changing \code{suffix} to avoid the ambiguity.
#' @return Returns a struct class object. See \link[espadon]{espadon.class} for 
#' class definitions.
#' @seealso \link[espadon]{struct.from.bin}.
#' @examples
#' # loading of toy-patient objects (decrease dxyz and increase beam.nb for 
#' # better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("rtdose"),
#'                              dxyz = rep (step, 3), beam.nb = 3)
#' D <- patient$rtdose[[1]]  
#' 
#' # isodose 50% Dmax Gy and 90% Dmax
#' bin50 <- bin.from.vol (D, min = 0.5 * D$max.pixel)
#' bin90 <- bin.from.vol (D, min = 0.9 * D$max.pixel)
#' S.isodose50 <- struct.from.bin (bin50, roi.name = "50pc" , 
#'                                 roi.color = "#00FFFF") 
#' S.isodose90 <- struct.from.bin (bin90, roi.name = "90pc" , 
#'                                 roi.color = "#FFFF00")
#' S <- struct.merge (S.isodose50, S.isodose90, alias = "isodose",
#'                    description = paste ("isodose of", D$object.alias)) 
#' # Dmax location :
#' z.dmax <- get.xyz.from.index(which (D$vol3D.data == D$max.pixel), D)[1,3]                 
#' display.plane(top = D, struct = S, view.coord = z.dmax, legend.shift = -50)                               


#' @export
#' @importFrom methods is

struct.merge <- function (ref.struct, add.struct, roi.name = NULL, 
                          roi.sname = NULL, 
                          roi.idx=NULL,
                          suffix = "", alias = "", description = "") {
  
  if (!is.null(ref.struct) & !is (ref.struct, "struct")) {
    warning ("ref.struct should be a struct class object.")
    return (NULL)
  }
  
  if (!is.null(add.struct) & !is (add.struct, "struct")) {
    warning ("add.struct should be a struct class object.")
    return (NULL)
  }
  
  if (is.null(ref.struct) & is.null(add.struct)) return(NULL)
  
  if  (is.null(ref.struct)) {
    roi2 <- select.names (add.struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
    if (is.null (roi2)) return(NULL)
    
    add.struct$roi.data <- add.struct$roi.data[roi2]
    add.struct$roi.info <-  add.struct$roi.info[roi2, ]
    
    add.struct$nb.of.roi <- length (add.struct$roi.data)
    add.struct$file.basename <- ""
    add.struct$file.dirname <- ""
    add.struct$object.name <- alias
    add.struct$object.alias <- alias
    add.struct$description <- description
    add.struct$object.info <- NULL
    add.struct$ref.object.name <- NULL
    add.struct$ref.object.info <- NULL
    if (!is.null(add.struct$error)) add.struct$error<-NULL
    return (add.struct)
  }
  
  ref.struct$file.basename <- ""
  ref.struct$file.dirname <- ""
  ref.struct$object.name <- alias
  ref.struct$object.alias <- alias
  ref.struct$description <- description
  ref.struct$object.info <- NULL
  ref.struct$ref.object.name <- NULL
  ref.struct$ref.object.info <- NULL
  if (!is.null(ref.struct$error)) ref.struct$error<-NULL
  
  roi2 <- select.names (add.struct$roi.info$roi.pseudo, roi.name, roi.sname, roi.idx)
  if (is.null (roi2)) return(ref.struct)

  

  if (ref.struct$thickness != add.struct$thickness | 
      any (ref.struct$ref.from.contour != add.struct$ref.from.contour) |
      ref.struct$ref.pseudo !=  add.struct$ref.pseudo) {
    warning ("structures are not compatible (thickness or ref).")
    return (NULL)
  }
  
  roi1 <- select.names (ref.struct$roi.info$roi.pseudo)
  to.add <- add.struct$roi.info[roi2, ]
  
  
  if (suffix!=""){ 
    to.add$name <- paste(to.add$name, suffix, sep="-")
    to.add$roi.pseudo <- paste(to.add$roi.pseudo, suffix, sep="-")
  }
  if (any (ref.struct$roi.info$RoI.pseudo[roi1] %in% to.add$roi.pseudo) |
      any (ref.struct$roi.info$name[roi1] %in% to.add$name))
    stop ("some RoI in ref.struct and add.struct share same name or roi.pseudo. Change suffix")
  

  ref.struct$roi.info <- rbind (ref.struct$roi.info, to.add)
  row.names(ref.struct$roi.info) <- NULL
  
  to.add_ <- add.struct$roi.data[roi2]
  ref.struct$roi.data <- c(ref.struct$roi.data, to.add_)
  # names(ref.struct$roi.data) <- ref.struct$roi.info$name
  
  keep <- !sapply(ref.struct$roi.data,is.null)
  ref.struct$roi.data <- ref.struct$roi.data[keep]
  ref.struct$roi.info <- ref.struct$roi.info[keep, ]
  ref.struct$nb.of.roi <- length (ref.struct$roi.data)
  return (ref.struct)
}