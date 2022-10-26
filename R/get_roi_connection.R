#' List of inter-connections between RoI
#' @description The \code{get.roi.connection} function describes the interconnections 
#' between Regions of Interest (RoI), from an imaging volume of "cluster" modality,
#' created by struct.clustering.
#' @param vol "volume" class object of "cluster" modality, created by 
#' \link[espadon]{struct.clustering}
#' @return Returns the list of regions of interest (RoI), where each item in the 
#' list gives the inter-connections with other RoI.
#' @seealso \link[espadon]{struct.clustering}
#' @export
#'
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c ("mr", "rtstruct"),  
#'                              dxyz = rep (step, 3))
#' MR <- patient$mr[[1]]
#' S <- patient$rtstruct[[1]]
#' cluster.vol <- struct.clustering (MR, S, T.MAT = patient$T.MAT, verbose = FALSE)
#' 
#' get.roi.connection (cluster.vol)
#' 
get.roi.connection <- function (vol) {
  #, clust.name = vol$cluster.info$label[2]
  if (!is (vol, "volume")){
    warning ("vol should be a volume class object.")
    return (NULL)
  }
  
  if ((vol$modality!="cluster") | is.null(vol$cluster.info)) {
    warning ("vol must be of cluster modality.")
    return (NULL)
  }
  

  clust.names <- unique(unlist(sapply (vol$cluster.info$label, function (s) {
    unlist(strsplit (as.character(s), "[|]")[[1]])
  })))
  clust.names <- clust.names[-(clust.names=="NA")]
  
  if (length(clust.names)==0) return(NULL)
  
  L <- lapply(clust.names, function(clust.name){
    idx <- sapply (vol$cluster.info$label, function (s) {
      s <- any (strsplit (as.character(s), "[|]")[[1]] == clust.name)
    })
    db <- vol$cluster.info[idx, ]
    row.names(db) <- NULL
    return (db)
  })
  names(L) <- clust.names
  L
}

