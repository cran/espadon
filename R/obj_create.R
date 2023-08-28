#' Espadon object creating
#' @description The \code{obj.create} function creates an espadon object with 
#' the essential properties it must have.
#' @param class Character string, representing an espadon class from among 
#' "volume", "struct" or "mesh".
#' @param alias Character string, \code{$alias} of the created object.
#' @return Returns a espadon class object (see \link[espadon]{espadon.class} 
#' for class definitions).

#' @examples
#' # Creation of an espadon mesh of a cube
#' M <- obj.create (class = "mesh")
#' M$mesh <- Rvcg::vcgIsotropicRemeshing (Rvcg::vcgBox(),0.5) 
#' M$nb.faces <- ncol (M$mesh$it)
#' rgl::wire3d (M$mesh)
#' @export
obj.create <- function (class = c ("","volume", "struct", "mesh"), alias = ""){
  obj <- list()
  obj$patient <- ""
  obj$patient.name <- ""
  obj$patient.bd <- ""
  obj$patient.sex <- ""
  obj$file.basename <- ""
  obj$file.dirname <- ""
  obj$object.name <- ""
  obj$object.alias <- alias
  obj$frame.of.reference <- ""
  obj$ref.pseudo <- "ref1"
  obj$modality <- ""
  obj$description <- ""
  obj$acq.date <- ""
  obj$study.date <- ""
  obj$creation.date <- format(Sys.Date(), "%Y%m%d")
  obj$number <- 1
  
  switch(class[1],
         "struct" = {
           obj$modality <- "rtstruct"
           obj$nb.of.roi <- 0
           obj$thickness <- NA
           obj$ref.from.contour <- diag(4)
           obj$roi.info <- data.frame(number=character(0),name=character(0),description=character(0),
                                      generation.algorithm = character(0), color=character(0), 
                                      roi.pseudo = character(0))
           obj$roi.obs <- data.frame(nb=character(0),roi.nb=character(0),label=character(0),
                                     code.value = character(0), code.scheme=character(0), 
                                     code.scheme.v=character(0), code.meaning=character(0),
                                     type=character(0), interpreter =character(0))
           
           obj$roi.data <- list()
           class(obj) <- "struct"
         },
         
         "mesh" = {
           obj$nb.faces <- 0
           obj$mesh <- list()
           class(obj) <- "mesh"
         },
         
         "vol" = {
           obj$unit = ""
           obj$n.ijk <- c(10,10,1)
           obj$slice.thickness = 1
           obj$min.pixel <- 0
           obj$max.pixel <- 0
           obj$dxyz <- c(1,1,1)
           obj$orientation <- c(1,0,0,0,1,0)
           obj$xyz0 <- matrix(c(0,0,0), ncol =3, dimnames= list(NULL, c("x0","y0","z0")))
           obj$xyz.from.ijk <- diag(4)
           obj$k.idx <- 0
           obj$missing.k.idx <- FALSE
           obj$cube.idx <- matrix(c(0,9,9,0,0,9,9,0,0,0,9,9,0,0,9,
                                    9,rep(0,8),rep(1,8)),ncol=8,byrow = TRUE,
                                  dimnames = list(c("i","j","k","t"), NULL))
           obj$vol3D.data <- array(0,dim=c(10,10,1))
           class(obj) <- "volume"
         }
  )

  return(obj)
}