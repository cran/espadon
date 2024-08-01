#' Indices of voxels crossed by a fan
#' @description The \code{fan.to.voxel} function computes the indices of voxels 
#' crossed by a fan. It is useful for retrieving voxel values and voxel indices 
#' of a volume (dose or ct) along the fan rays.
#' @param vol "volume" class object.
#' @param fan "fan" class object created by \link[espadon]{fan.sphere} for example.
# @param att Boolean. If \code{TRUE}, the cumulative sum  of the \code{vol} voxels 
# values along the ray path, times the distance crossed by the ray in these 
# voxels, is returned.
#' @param restrict Boolean. If \code{TRUE}, only the voxels with a value equal 
#' to \code{vol.value} are taken into account.
#' @param vol.value Value of the voxels taken into account, in case of \code{restrict = TRUE}
#' @return Returns a dataframe of 4 columns. Each line gives:
#' \itemize{
#' \item column "ray.index": the index (i.e. the row number) of the ray 
#' concerned in \code{fan$dxyz},
#' \item column "vol.index": the index of the voxel crossed in \code{vol$vol.3Ddata},
#' \item column "l.in": the distance between fan source (i.e. \code{fan$origin}) 
#' and the first face of the voxel crossed by the ray,
#' \item column "dl": the distance crossed by the ray in the voxel.
# \item if \code{att = TRUE} the cumulative sum of the voxels 
# values along the ray path, times dl (column "att").
#' }
#' If the rays do not cross any voxel, the dataframe has no row.
#' @seealso \link[espadon]{fan.beam}, \link[espadon]{fan.planar}, \link[espadon]{fan.sphere}. 
#' @export
#' @examples
#' vol <- vol.create (pt000 = c(1,10,10), dxyz = c (1 , 1, 1),
#'                    n.ijk = c(100, 100, 100)) 

#' fan.origin <- c (50,50,50)                         
#' fan <- fan.sphere (angle = 10, origin = fan.origin)

#' fan.voxel <- fan.to.voxel (vol = vol, fan = fan)
#' head (fan.voxel)
#' 
#' # Use of the 2nd column of fan.voxel to select voxels 
#' bin <- vol.copy (vol, modality = "binary")
#' bin$vol3D.data[] <- FALSE
#' bin$vol3D.data[fan.voxel[,2]] <- TRUE
#' bin$max.pixel <- TRUE
#' bin$min.pixel <- FALSE
#' display.kplane(bin, k=10)
 
fan.to.voxel <- function (vol, fan, restrict = FALSE, vol.value = 1){
  
  
  if (!is (vol, "volume")){
    stop ("vol should be a volume class object.")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
  
  if (!is (fan, "fan")){
    stop ("fan should be a fan class object.")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
  
  if (fan$ref.pseudo!=vol$ref.pseudo){
    warning ("vol ref.pseudo and fan ref.pseudo are different")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
   
  #2D
  idx.c <- which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),2,sum)==0) 
  idx.r <-  which(apply(abs(vol$xyz.from.ijk[1:3,1:3]),1,sum)==0)
  if (length(idx.c)>0) {
    if (abs(vol$xyz0[1,idx.r])>1e-6) return(NULL)
    u <- vol$xyz.from.ijk 
    u[idx.r,idx.c]<- 1
    Mat <- solve(u)
    Mat[idx.r,idx.c] <- 0
  } else {#3D
    Mat <- solve(vol$xyz.from.ijk)}
  
  
  u_ijk <- (cbind(fan$xyz,0) %*% t(Mat))[,1:3]
  O_ijk <-  as.numeric((c(fan$origin,1) %*% t(Mat))[,1:3])
  k_idx<- match(0:max(vol$k.idx),vol$k.idx)
  k_loc <- k_idx-1
  fna <- is.na(k_idx)
  k_idx[!fna] <- vol$k.idx
  k_loc[fna] <- max(vol$k.idx)+1
  k_idx[fna] <- max(vol$k.idx)+1
  n_ijk  <- as.numeric(vol$n.ijk)
  
  p <- as.numeric(t(u_ijk))
  ncol <- 4
  # if (att) ncol <- 5
  df <- as.data.frame (matrix(.fantovoxelC(p, n_ijk, k_idx, 
                                           k_loc,O_ijk, 
                                           vol_data = as.numeric(vol$vol3D.data),
                                           att = FALSE,
                                           vol_value_flag = restrict,
                                           vol_value = as.numeric(vol.value)), 
    ncol = ncol, byrow=TRUE))
  colnames(df) <- c("ray.index","vol.index","l.in","dl","att")[1:ncol]
  return(df)
}

  # df <-do.call(rbind,lapply(1:nrow(u_ijk), function(m) {
  #   l <- c(0,round(c((0:vol$n.ijk[1] - 0.5 - O_ijk[1]) / u_ijk[m,1],
  #                    (0:vol$n.ijk[2] - 0.5 - O_ijk[2]) / u_ijk[m,2], 
  #                    (c(k_idx[1]-0.5,(k_idx[1]:k_idx[length(k_idx)])+0.5)-O_ijk[3])/ u_ijk[m,3]),6))
  #   f <- which(!is.infinite(l) & l>=0)
  #   if (length(f)==0) return(NULL)
  #   l <- unique(sort(l[f]))
  #   M_mid<- floor(matrix(rep(O_ijk,length(l)-1),ncol = 3, byrow=T) +  
  #                   0.5*(l[1:(length(l)-1)] + l[2:length(l)])*u_ijk[rep(m,length(l)-1),] + 0.5)
  #   f <- which(M_mid[,1]>=0 & M_mid[,1]<vol$n.ijk[1] & M_mid[,2]>=0 & M_mid[,2]<vol$n.ijk[2] &
  #                M_mid[,3]>= k_idx[1] & M_mid[,3]<=k_idx[length(k_idx)])
  #   if (length(f)==0) return(NULL)
  #   Lin <- l[f]
  #   Lt <- (l[2:length(l)]-l[1:(length(l)-1)])[f]
  #   M_mid <- matrix(as.numeric(M_mid[f,]),nrow=length(Lt), byrow = F)
  #   f <- which(k_loc[M_mid[,3] + 1] == M_mid[,3] )
  #   if (length(f)==0) return(NULL)
  #   vol_index = M_mid[f,1] + M_mid[f,2] * vol$n.ijk[1] +
  #     (k_loc[M_mid[f,3] + 1] *vol$n.ijk[1] *  vol$n.ijk[2]) + 1
  #   data.frame(ray.index=m, vol.index = vol_index,Lin=Lin[f], Lt = Lt[f])
  # }))

.fan.to.voxel.with.att <- function (vol, fan, att,restrict = FALSE, vol.value = 1){
  
  
  if (!is (vol, "volume")){
    stop ("vol should be a volume class object.")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
  
  if (!is (fan, "fan")){
    stop ("fan should be a fan class object.")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
  
  if (fan$ref.pseudo!=vol$ref.pseudo){
    warning ("vol ref.pseudo and fan ref.pseudo are different")
    # return (data.frame(ray.index=numeric(0), vol.index=numeric(0),
    #                    l.in=numeric(0), dl=numeric(0)))
  }
  
  Mat <-solve (vol$xyz.from.ijk) 
  u_ijk <- (cbind(fan$xyz,0) %*% t(Mat))[,1:3]
  O_ijk <-  as.numeric((c(fan$origin,1) %*% t(Mat))[,1:3])
  k_idx<- match(0:max(vol$k.idx),vol$k.idx)
  k_loc <- k_idx-1
  fna <- is.na(k_idx)
  k_idx[!fna] <- vol$k.idx
  k_loc[fna] <- max(vol$k.idx)+1
  k_idx[fna] <- max(vol$k.idx)+1
  n_ijk  <- as.numeric(vol$n.ijk)
  
  p <- as.numeric(t(u_ijk))
  ncol <- 4
  if (att) ncol <- 5
  df <- as.data.frame (matrix(.fantovoxelC(p, n_ijk, k_idx, 
                                           k_loc,O_ijk, 
                                           vol_data = as.numeric(vol$vol3D.data),
                                           att = att,
                                           vol_value_flag = restrict,
                                           vol_value = as.numeric(vol.value)), 
                              ncol = ncol, byrow=TRUE))
  colnames(df) <- c("ray.index","vol.index","l.in","dl","att")[1:ncol]
  return(df)
}