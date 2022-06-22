#' Adding spherical coordinates to a mesh
#' @description The \code{mesh.spheric.proj} function adds latitude and longitude
#' coordinates to a mesh. These features map the mesh surface to a sphere. 
#' Latitude and longitude are computed using the heat diffusion approach explained by
#' \emph{Brechbühler and al} \strong{\[1\]}. 
#' @param mesh "mesh" class object.
#' @param verbose Boolean, by default set to \code{FALSE}. Allows you to inhibit 
#' comments.

#' @return
#' returns a "mesh" class object in which \code{$mesh} contains \code{Lat} 
#' and \code{lon} evaluated at vertices.
#' The function allows to have a parameterized surface  for later computations 
#' as curvature or shape index, hence, nor the surface, nor the angles are preserved.
#' In the DICOM frame of reference, latitude goes along Z axis (from feet = -1 to 
#' head = +1) and longitude turns counter clockwise (from -1 to +1).
#' @note This funtion is time consuming.

#' @importFrom Rdpack reprompt
#' @references 
#' \strong{\[1\]} \insertRef{BRECHBUHLER1995154}{espadon}

#' @examples
#' # loading of toy-patient objects (decrease dxyz for better result)
#' step <- 5
#' patient <- toy.load.patient (modality = c("ct", "rtstruct"), roi.name = "", 
#'                              dxyz = rep (step, 3))
#' CT <- patient$ct[[1]]
#' S <- patient$rtstruct[[1]]
#'
#' #creation of the patient mesh
#' bin <- bin.from.roi (CT, struct = S, roi.name = "patient")
#' m.patient <- mesh.from.bin (bin)
#' m.skin <- mesh.repair (m.patient, verbose = FALSE)
#'
#' m.proj <- mesh.spheric.proj (m.skin, verbose = FALSE)
#'
#' library (rgl)
#' col <- hcl.colors (12, "Blue-Red 3")
#' open3d()
#' shade3d (m.proj$mesh, meshColors = "vertices",
#'          color = col[round ((m.proj$mesh$Lat/2 + 0.5) * 11) + 1],
#'          specular = "#404040")
#' open3d()         
#' shade3d (m.proj$mesh, meshColors = "vertices",
#'          color = col[round ((m.proj$mesh$Lon/2 + 0.5) * 11) + 1],
#'          specular = "#404040")

#'
#' @importFrom Rvcg vcgVFadj
#' @import Matrix
#' @importFrom methods is
#' @export

mesh.spheric.proj <- function (mesh, verbose = TRUE) {
  
# ○n utilise rvcg et Matrix
  
  if (!is (mesh, "mesh")) {
    warning ("mesh should be a mesh class object.")
    return (NULL)
  }  
  
  m <- mesh$mesh
  
  if (verbose) cat ("mesh.spheric.proj INFO : computing connections\n")
  
#  calcul des conections
  
  f.list <- vcgVFadj (m) ## liste des triangles connectés à chaque vertex
  v.list <- lapply (1:length (f.list), function (v) { ## liste des vertex connectés à chaque vertex
    sel <- unique (as.vector (m$it[,f.list[[v]]]))
    return (list (v=v, n.v = sel[sel != v]))
  })
  
#  construcion de la latitude
  
  idx <- which (m$vb[3, ] == max (m$vb[3, ]))
  x.g <- mean (m$vb[1, idx])
  y.g <- mean (m$vb[2, idx])
  idx.keep <- which.min ((m$vb[1, idx] - x.g)^2 + (m$vb[2, idx] - y.g)^2)
  top <- idx[idx.keep]
  
  idx <- which (m$vb[3, ] == min (m$vb[3, ]))
  x.g <- mean (m$vb[1, idx])
  y.g <- mean (m$vb[2, idx])
  idx.keep <- which.min ((m$vb[1, idx] - x.g)^2 + (m$vb[2, idx] - y.g)^2)
  bot <- idx[idx.keep]
  
# mesh <- vcgSmooth(mesh, iteration = 3, type="HClaplace")
# shade3d (mesh, col="burlywood2", specular = "#404040")
# wire3d (mesh, col="burlywood2", specular = "#404040")
  
# spheres3d (t (mesh$vb[, top]), col="red", radius=3)
# spheres3d (t (mesh$vb[, bot]), col="green", radius=3)
  
  
#  idx <- which (sapply (f.list, function (t) length (t)) < 3)
#  mesh$it <- mesh$it[, -unlist (f.list[idx])]
#  mesh <- tmesh3d (mesh$vb, mesh$it)
  
  if (verbose) cat ("mesh.spheric.proj INFO : building Latitude\n")
  
  A <- Matrix (0, nrow=dim(m$vb)[2], ncol=dim(m$vb)[2], sparse=TRUE)
  
  diag.el.idx <- sapply (v.list, function (v) v$v)
  diag.el.val <- sapply (v.list, function (v) length (v$n.v))
  A[cbind (diag.el.idx, diag.el.idx)] <- diag.el.val
  
  out.diag.row <- unlist (sapply (v.list, function (v) rep (v$v, length (v$n.v))))
  out.diag.col <- unlist (sapply (v.list, function (v) v$n.v))
  A[cbind (out.diag.row, out.diag.col)] <- -1
  
  B <- Matrix (0, nrow=dim (m$vb)[2], ncol=1, sparse=TRUE)
  
  z.min <- min (m$vb[3, ])
  z.max <- max (m$vb[3, ])
  A[top, ] <- 0
  A[bot, ] <- 0
  A[top, top] <- 1
  A[bot, bot] <- 1
  B[top, 1] <- 1
  B[bot, 1] <- -1
  
  if (verbose) cat ("mesh.spheric.proj INFO : system inversion\n")
  
  Lat <- as.vector (solve (A, B)) # on doit pouvoir accélérer ceci
  
#  qr.A <- qr (A) ## plutôt qu'inverser la matrice, on fait une decomposition QR
#  Lat <- qr.coef (A, B) ## note, on peut aussi avoir le résidu au besois avec qr.resid
  
  ecdf(Lat)(0)
  Lat <- 2 * (ecdf (Lat)(Lat) - 0.5)
  
#  ligne de cht de date
  
  if (verbose) cat ("mesh.spheric.proj INFO : date change line\n")
  
  date.line <- top
  while (date.line[1] != bot) {
    idx <- which.min (Lat[v.list[[date.line[1]]]$n.v])
    date.line <- c(v.list[[date.line[1]]]$n.v[idx], date.line)
  }
  date.line <- rev (date.line) # du N au S
  
  candidates <- which ((m$it[1, ] %in% date.line + m$it[2, ] %in% date.line
                        + m$it[3, ] %in% date.line) == 2)
  
  idx <- sapply (candidates, function (t) {
    pos <- m$it[, t] %in% date.line
    if (pos[1] & pos[2]) {
      if (which (date.line == m$it[1, t]) < which (date.line == m$it[2, t])) return (TRUE)
    }
    if (pos[2] & pos[3]) {
      if (which (date.line == m$it[2, t]) < which (date.line == m$it[3, t])) return (TRUE)
    }
    if (pos[3] & pos[1]) {
      if (which (date.line == m$it[3, t]) < which (date.line == m$it[1, t])) return (TRUE)
    }
    return (FALSE)
  })
  
  date.line.right <- sapply (candidates[idx], function (t) {
    pos <- !m$it[, t] %in% date.line
    m$it[pos, t]
  })
  
  
  candidates <- which ((m$it[1, ] %in% date.line + m$it[2, ] %in% date.line + m$it[3, ] %in% date.line) == 1)
  
  repeat {
    add.to.date.line.right <- sapply (candidates, function (t) {
      pos <- !m$it[, t] %in% date.line[-c (1, length (date.line))]
      if (all (pos)) return (FALSE)
      any (m$it[pos, t] %in% date.line.right)
    })
    
    date.line.right_ <- unique (as.vector (sapply (candidates[add.to.date.line.right], function (t) {
      pos <- !m$it[, t] %in% date.line
      m$it[pos, t]
    })))
    
    if (length (date.line.right_) == length (date.line.right)) break
    
    date.line.right <- date.line.right_
  }
  
#  longitude
  
  if (verbose) cat ("mesh.spheric.proj INFO : building longitude\n")
  
  A <- Matrix (0, nrow=dim(m$vb)[2], ncol=dim(m$vb)[2], sparse=TRUE)
  
  diag.el.idx <- sapply (v.list, function (v) v$v)
  diag.el.val <- sapply (v.list, function (v) length (v$n.v))
  A[cbind (diag.el.idx, diag.el.idx)] <- diag.el.val
  
  out.diag.row <- unlist (sapply (v.list, function (v) rep (v$v, length (v$n.v))))
  out.diag.col <- unlist (sapply (v.list, function (v) v$n.v))
  A[cbind (out.diag.row, out.diag.col)] <- -1
  
  B <- Matrix (0, nrow=dim (m$vb)[2], ncol=1, sparse=TRUE)
  
#  on annule N & S
  A[top, ] <- 0
  A[top, top] <- 1
  A[bot, ] <- 0
  A[bot, bot] <- 1
#  on déconnecte N & S
  N.con <- v.list[[top]]$n.v
  A[N.con, top] <- 0
  A[cbind (N.con, N.con)] <- A[cbind (N.con, N.con)] - 1
  S.con <- v.list[[bot]]$n.v
  A[S.con, bot] <- 0
  A[cbind (S.con, S.con)] <- A[cbind (S.con, S.con)] - 1
#  on définit la condition cyclique
  date.line_ <- date.line[2:(length (date.line) - 1)]
  date.line.n <- sapply (date.line_, function (v) {
    sum (v.list[[v]]$n.v %in% date.line.right)
  })
  B[date.line_] <- date.line.n
  date.line.right.n <- sapply (date.line.right, function (v) {
    sum (v.list[[v]]$n.v %in% date.line_)
  })
  B[date.line.right] <- -date.line.right.n
#  On définit la température d'un vertex
  A[1, ] <- 0
  A[1, 1] <- 1
  B[1] <- 0
  
  if (verbose) cat ("mesh.spheric.proj INFO : system inversion\n")
  
  Lon <- as.vector (solve (A, B))
  
  Lon <- 2 * (Lon - min (Lon)) / (max (Lon) - min (Lon)) - 1
  Lon <- 2 * ecdf (Lon)(Lon) - 1
  
#  fin du meshage
  
  mesh$mesh$Lat <- Lat
  mesh$mesh$Lon <- Lon
  
  return (mesh)
}