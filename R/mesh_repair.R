#' Repair of a mesh
#' @description The \code{mesh.repair} function repairs holes in a \code{mesh} 
#' class object.
#' @param mesh "mesh" class object.
#' @param verbose Boolean, by default set to \code{FALSE}. Allows you to inhibit 
#' comments.
#' @return Returns a mesh, repaired by removing degenerated triangles and filling holes.
#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 4
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), 
#'                              roi.name = "gizzard", dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#' 
#' # creation of the gizzard mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "gizzard")
#' mesh.gizzard <- mesh.from.bin (bin, alias = "gizzard", verbose = FALSE)
#' 
#' repair.mesh.gizzard <- mesh.repair (mesh.gizzard, verbose = FALSE)
#' str (mesh.gizzard)
#' str (repair.mesh.gizzard)

#' @importFrom Rvcg vcgBorder vcgIsolated vcgSmooth vcgVFadj
#' @importFrom methods is
#' @export
mesh.repair <- function (mesh, verbose = TRUE) {
  
  
  if (!is (mesh, "mesh")) {
    warning ("mesh should be a mesh class object.")
    return (NULL)
  }
  mesh_ <- mesh
  mesh <- mesh_$mesh
  ##
  ## Elimination des triangles dégénérés
  ##
  
  if (verbose) cat ("mesh.repair INFO : removing degenerated triangles\n")
  repeat {
    ## listes de connexions
    f.list <- vcgVFadj (mesh) ## liste des triangles connectés à chaque vertex
    v.list <- lapply (1:length (f.list), function (v) { ## liste des vertex connectés à chaque vertex
      sel <- unique (as.vector (mesh$it[,f.list[[v]]]))
      return (list (v=v, n.v = sel[sel != v]))
    })
    v <- sort (which (sapply (v.list, function (v) length (v$n.v)) == 2), decreasing = TRUE)[1]
    if (is.na (v)) break
    
    if (verbose) cat ("mesh.repair INFO : cleaning", v, "\n")
    mesh$vb <- mesh$vb[, -v]
    t <- which (mesh$it == v, arr.ind = TRUE)[,2]
    mesh$it <- mesh$it[, -t]
    sel <- which (mesh$it > v)
    mesh$it[sel] <- mesh$it[sel] - 1
    mesh$normals <- mesh$remvert <- NULL
  }
  mesh <- vcgIsolated(mesh, silent = !verbose)
  
  if (verbose) cat ("mesh.repair INFO : removing 3 borders triangles\n")
  borders <- vcgBorder(mesh)$bordervb
  t.b <- apply (mesh$it, 2, function (t) !all(borders[t]))
  if (verbose) cat ("mesh.repair INFO : number of tr :", ncol (mesh$it) - sum (t.b) ,"\n")
  mesh$it <- mesh$it[, t.b]
  mesh$normals <- mesh$remvert <- NULL
  mesh <- vcgIsolated(mesh, silent = !verbose)
  
  
  
  if (verbose) cat ("mesh.repair INFO : removing >2 borders border vertex\n")
  
  repeat {
    f.list <- vcgVFadj (mesh) ## liste des triangles connectés à chaque vertex
    v.list <- lapply (1:length (f.list), function (v) { ## liste des vertex connectés à chaque vertex
      sel <- unique (as.vector (mesh$it[,f.list[[v]]]))
      return (list (v=v, n.v = sel[sel != v]))
    })
    
    borders <- vcgBorder(mesh)$bordervb
    b.count <- sapply (which (borders), function (v) {
      sum (borders[v.list[[v]]$n.v])
    })
    
    v <- sort (which (borders)[b.count > 2], decreasing = TRUE)[1]
    if (is.na (v)) break
    
    if (verbose)cat ("mesh.repair INFO : cleaning", v, "\n")
    mesh$vb <- mesh$vb[, -v]
    t <- which (mesh$it == v, arr.ind = TRUE)[,2]
    mesh$it <- mesh$it[, -t]
    sel <- which (mesh$it > v)
    mesh$it[sel] <- mesh$it[sel] - 1
    mesh$normals <- mesh$remvert <- NULL
  }
  mesh <- vcgIsolated(mesh, silent = !verbose)
  
  
  ##
  ## réparation des trous
  ##
  
  if (verbose) cat ("mesh.repair INFO : mesh holes filling\n")
  
  ## listes de connexions
  f.list <- vcgVFadj (mesh) ## liste des triangles connectés à chaque vertex
  v.list <- lapply (1:length (f.list), function (v) { ## liste des vertex connectés à chaque vertex
    sel <- unique (as.vector (mesh$it[,f.list[[v]]]))
    return (list (v=v, n.v = sel[sel != v]))
  })
  # table (sapply (f.list, function (v) length (v)))
  # table (sapply (v.list, function (v) length (v$n.v)))
  
  repeat {
    borders <- vcgBorder(mesh)$bordervb
    if (verbose) cat ("mesh.repair INFO : remaining borders :", sum (borders), "\n")
    b.list <- which (borders)[1]
    if (is.na (b.list)) break
    repeat {
      if (rev(b.list)[1] %in% v.list[[b.list[1]]]$n.v & length (b.list) > 2) break
      T1 <- which (mesh$it == b.list[1], arr.ind = TRUE)
      for (Tr in T1[, 2]) {
        id1 <- which (mesh$it[, Tr] == b.list[1])
        id2 <- which (mesh$it[, Tr] != b.list[1] & borders[mesh$it[, Tr]])
        if (length (id2) == 1) {
          if (id1 == 1 & id2 == 3) break
          if (id1 == 2 & id2 == 1) break
          if (id1 == 3 & id2 == 2) break
        } else if (length (id2) == 2) {
          stop ()
        }
      }
      b.list <- c(mesh$it[id2, Tr], b.list)
    }
    G <- apply (mesh$vb[, b.list], 1, mean)
    mesh$vb <- cbind (mesh$vb, G)
    dimnames(mesh$vb) <- NULL
    tr <- matrix (c (rep (ncol (mesh$vb), length (b.list)), b.list[-1], b.list[1], b.list), nrow = 3, byrow=TRUE)
    mesh$it <- cbind (mesh$it, tr)
    mesh$normals <- mesh$remvert <- NULL
    
    f.list <- vcgVFadj (mesh) ## liste des triangles connectés à chaque vertex
    v.list <- lapply (1:length (f.list), function (v) { ## liste des vertex connectés à chaque vertex
      sel <- unique (as.vector (mesh$it[,f.list[[v]]]))
      return (list (v=v, n.v = sel[sel != v]))
    })
  }
  
  mesh_$mesh <- mesh
  mesh_$nb.faces <- dim(mesh$it)[2]
  class(mesh_)
  return (mesh_)
}