#' Distance-based spatial similarity metrics calculated from the mesh.
#' @description The \code{sp.similarity.from.mesh} function computes Hausdorff 
#' distances and surface Dice similarity coefficient.
#' @param mesh1,mesh2 espadon mesh class objects
#' @param hausdorff.coeff Vector indicating the requested Hausdorff distance metrics from among 
#' 'HD.max','HD.mean'. Equal to \code{NULL} if not requested.
#' \code{NULL}, it will be added. 
#' @param hausdorff.quantile numeric vector of probabilities with values between 0 and 1,
#' representing the quantiles of the unsigned distances between \code{mesh1} and \code{mesh2}.
#' Equal to \code{NULL} if not requested.
#' @param surface.DSC.tol numeric vector representing the maximum margins of 
#' deviation which may be tolerated without penalty. Equal to \code{NULL} if not requested.
#' @return Returns a list containing (if requested): 
#' \itemize{
#' \item \code{Hausdorff} : dataframe including the maximum, mean and quantiles 
#' \item \code{surfaceDSC} : dataframe with the columns \code{tol} and \code{surfaceDSC},
#' representing respectively the requested tolerances and the surface Dice 
#' similarity coefficients defined by \emph{Nikolov et al} \strong{\[1\]}
#' }
#' @importFrom Rdpack reprompt
#' @references \strong{\[1\]} \insertRef{Nikolov2018DeepLT}{espadon}

#' @seealso \link[espadon]{sp.similarity.from.bin}
#' @examples
#' library (Rvcg)
#' # espadon mesh of two spheres of radius R1=10 and R2=11, separated by dR = 3
#' sph <- vcgSphere ()
#' mesh1 <- obj.create ("mesh")
#' mesh1$nb.faces <- ncol (sph$it)
#' mesh1$mesh <- sph
#' mesh2 <-  mesh1
#' 
#' R1 <- 10
#' R2 <- 11
#' dR <- 3
#' mesh1$mesh$vb[1:3,] <- R1 * mesh1$mesh$normals[1:3,] + mesh1$mesh$vb[1:3,]
#' mesh2$mesh$vb[1:3,] <- R2 * mesh2$mesh$normals[1:3,] + mesh2$mesh$vb[1:3,] +
#'                  matrix (c (dR, 0, 0), ncol = ncol (mesh2$mesh$vb), nrow = 3)
#'              
#' sp.similarity.from.mesh (mesh1 , mesh2, 
#'                          hausdorff.quantile = seq (0, 1, 0.05),
#'                          surface.DSC.tol = seq (0, dR + abs(R2-R1), 0.5))
#' @export
sp.similarity.from.mesh <- function (mesh1, mesh2, 
                                     hausdorff.coeff = c('HD.max', 'HD.mean'),
                                     hausdorff.quantile = c(0.5,0.95), 
                                     surface.DSC.tol = seq(0,10, 1)) {
  
  if (!all(hausdorff.quantile<=1 & hausdorff.quantile>=0)) stop("hausdorff.quantile must be between 0 and 1")
  if (!is.null(mesh1$ref.pseudo) & !is.null(mesh2$ref.pseudo))
    if (mesh1$ref.pseudo != mesh2$ref.pseudo)  
      warning("mesh1 and mesh2 do not share the same frame of reference")
  
  H1 <- dist.mesh ( mesh.test = mesh1, mesh.ref = mesh2, plot=FALSE, in.front = F)
  H2 <- dist.mesh ( mesh.test = mesh2, mesh.ref = mesh1, plot=FALSE, in.front = F)

 
  S1 <- H1$surface
  S2 <- H2$surface
  d1 <- abs(H1$distance)
  d2 <- abs(H2$distance)
  
  d <- c (d1, d2)
  nb <- length(d)
  max.d <- max(d)
  w <- c(S1/mean(S1),S2/mean(S2))
  # w <- rep(1, nb)
  ord <-order (d)
  w <- w[ord]
  d <- d[ord]
  df_ <- data.frame (d= as.factor(d), w= w)
  value<- by(df_, df_$d, FUN = function(x) sum(x$w))
  d_ <- c(0,as.numeric(names(value)))
  w_ <- c(0,as.numeric (value)/nb)
  
  mean.d <- sum(w_*d_)
  w_ <- cumsum(w_)
  
  hausdorff.coeff <- unique(hausdorff.coeff)
  hausdorff.quantile <- unique(hausdorff.quantile)
  result.quantile <- sapply(hausdorff.quantile, function(q){
    idx <- which(q<=w_)
    idx <- idx[1]
    return(d_[idx])
    # if (idx==1) return(d_[idx])
    # return(d_[idx] + (d_[idx-1]-d_[idx])*(q-w_[idx])/(w_[idx-1]- w_[idx]))
  })
  

  f.idx <- match(hausdorff.coeff, c('HD.max', 'HD.mean'))  
  f.idx <- f.idx[!is.na(f.idx)]
  label <- HD <-c()
  if (length(f.idx)!=0) { 
    label <- c('HD.max', 'HD.mean')[f.idx]
    HD <- c(max.d, mean.d)[f.idx]
    }
 
  if (length(hausdorff.quantile)!=0)   label <- c(label, paste("HD", hausdorff.quantile*100, sep="."))  
  L <- list()
  if (length(label)!=0)
    L[['Hausdorff']] <-data.frame(label = label, HD = c(HD, result.quantile))
  

  if (length(surface.DSC.tol)==0) return(L)
 
  S1.tot <- sum(S1)
  S2.tot <- sum(S2)
  surface.DSC <- data.frame(tol= surface.DSC.tol,
                            sDSC =sapply(surface.DSC.tol, 
                                         function(d) sum(S1[d1<=d]) + 
                                           sum(S2[d2<=d]))/ (S1.tot + S2.tot))
  L[['surfaceDSC']] <- surface.DSC
  return(L)
  
}


