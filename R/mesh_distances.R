#' @importFrom graphics boxplot
#' @importFrom rgl open3d points3d wire3d
#' @importFrom grDevices colorRampPalette
#' @importFrom Rvcg vcgClost vcgRaySearch
.mesh.distances <- function (mesh.test, mesh.ref, tol =1e-7, plot=TRUE, 
                             col.in = "blue", col.out = "green", pt.size = 10) {
  
  tol <- abs(tol)
  pt <- t(mesh.test$mesh$vb)
  it <- t(mesh.test$mesh$it)
  m <- vcgClost(mesh.test$mesh, mesh.ref$mesh)
  u <- t(m$vb)[,1:3]-pt[,1:3]
  dum <- u*u
  dum <- sqrt(dum[,1] + dum[,2] + dum[,3])
  fd <- dum!=0
  m$distance <- sign(m$quality)*dum
  in.front <- rep(TRUE,length(m$distance))
  diff  <- rep(0,length(m$distance))
  success <- in.front
  m2 <- list()
  m2$vb <- m$vb[,fd]
  m2$normals <- t(as.matrix(cbind(-u[fd,]/dum[fd],1)))
  class(m2) <- "mesh3d"
  d <- vcgRaySearch(m2,mesh.test$mesh)
  
  diff[fd] <- abs(m$distance[fd])-d$distance
  success[fd] <- d$quality ==1
  f.nok <-(diff<=-tol |  !success | (diff>=tol & diff<=1e-2))
  in.front[!f.nok] <- abs(diff[!f.nok]) < tol
  if (any(f.nok))
    in.front[f.nok]<- .meshinfront ( pt_x = as.numeric(pt[,1]),
                                     pt_y = as.numeric(pt[,2]),
                                     pt_z = as.numeric(pt[,3]),
                                     p2_x = as.numeric(t(m$vb)[f.nok,1]),
                                     p2_y = as.numeric(t(m$vb)[f.nok,2]),
                                     p2_z = as.numeric(t(m$vb)[f.nok,3]),
                                     u_x = -as.numeric(u[f.nok,1]),
                                     u_y = -as.numeric(u[f.nok,2]),
                                     u_z = -as.numeric(u[f.nok,3]),
                                     n_A = as.numeric(it[,1])-1,
                                     n_B = as.numeric(it[,2])-1,
                                     n_C = as.numeric(it[,3])-1)
  m$in.front <- in.front
  m$quality <- NULL
  
  if (plot) {
    palet <- colorRampPalette(c(col.in, "grey90",col.out))(100)
    dummy <- factor (sign(m$distance), levels=c(-1, 1), labels=c ("in", "out"))
    dev.new (width = 7, height = 5, noRStudioGD = TRUE)
    boxplot (d~in.out, data.frame (d= abs(m$distance), in.out=dummy), col=c(col.in ,col.out),
             xlab = paste(mesh.test$object.alias,"vs",mesh.ref$object.alias ), ylab = "d (mm)")
    ## affichage des mesures
    d.max <- max (abs(m$distance))
    open3d (windowRect = c (100, 100, 600, 600))
    f <- (m$in.front==1)
    wire3d(mesh.test$mesh, color="grey90", lit= FALSE)
    points3d(as.matrix(pt[f,1:3]), col=palet[1 + 99 * (m$distance[f] / d.max / 2 + 0.5)], size=pt.size)
    
    display.palette(palet, breaks=seq (-d.max, d.max, length.out=101))
    dev.new (width = 7, height = 5, noRStudioGD = TRUE)
    H <- hist (m$distance, breaks=seq (-d.max, d.max, length.out=100))
  }
  
  return(m)
}
  
  # pt.ref <- as.data.frame (t (mesh.ref$mesh$vb[1:3, ]))
  # pt <- t(mesh.test$mesh$vb)
  # it <- t(mesh.test$mesh$it)
  # tab <- data.frame(t(do.call(rbind.data.frame,vcgClost(x=pt[,1:3],
  #                                                       mesh.ref$mesh, barycentric = TRUE))))
  # row.names(tab) <- NULL
  # 
  # it.ref <-  t(mesh.ref$mesh$it)[tab$faceptr, ]
  # A <- pt.ref[it.ref[,1],]; B <- pt.ref[it.ref[,2], ]; C <- pt.ref[it.ref[,3], ]
  # P2 <-  (B-A)*tab[,c(12,12,12)] + (C-A)*tab[,c(13,13,13)] + A
  # row.names(P2) <- NULL
  # u <-P2-pt[,1:3]
  # dum <- u*u
  # dum <-sqrt(dum[,1] + dum[,2] + dum[,3])
  # tab <- data.frame(cbind(tab[,9],pt[ ,1:3],tab[,c(1:3,5:7,10)]))
  # colnames(tab) <- c("d","x.test","y.test","z.test","x.ref","y.ref","z.ref","ux","uy","uz","face.ref")
  # 
  # fd <- dum!=0 
  # tab[fd, 8:10] <- u[fd,]/dum[fd]
  # tab$d <- sign(tab$d)*dum
  # 
  # # 
  # # tab <- do.call(rbind.data.frame,lapply (1:ncol(mesh.test$mesh$vb), function (t.idx) {
  # #   pt <- as.numeric(pt[t.idx,1:3])
  # #   n.t <-  as.numeric(mesh.test$mesh$normals[1:3, t.idx])
  # #   d.vect <- cbind(ref[,1] - pt[1], ref[,2] - pt[2],ref[,3] - pt[3])
  # #   d <- sqrt (d.vect[,1]^2 + d.vect[,2]^2 + d.vect[,3]^2)
  # #   idx <- which (d==min(d))
  # #   pr <- ((ref[idx,1] - pt[1])*n.ref[idx,1] + (ref[idx,2] - pt[2])*n.ref[idx,2] + (ref[idx,3] - pt[3])*n.ref[idx,3])
  # #   idx.pr <-which.max(abs(pr))
  # #   pr <- pr[idx.pr]
  # #   idx <- idx[idx.pr]
  # #   io <- (pr>=0)
  # #   c(pt, as.numeric(ref[idx,]), d[idx], d.vect[idx,],  io)
  # # }))
  # # tab$in.front <- .meshdist ( pt_x = as.numeric(pt[,1]),
  # #                     pt_y = as.numeric(pt[,2]),
  # #                     pt_z = as.numeric(pt[,3]),
  # #                     p2_x = as.numeric(P2[,1]),
  # #                     p2_y = as.numeric(P2[,2]),
  # #                     p2_z = as.numeric(P2[,3]),
  # #                     u_x = -as.numeric(u[,1]),
  # #                     u_y = -as.numeric(u[,2]),
  # #                     u_z = -as.numeric(u[,3]),
  # #                     n_A = as.numeric(it[,1])-1,
  # #                     n_B = as.numeric(it[,2])-1,
  # #                     n_C = as.numeric(it[,3])-1)
  # in.front <- rep(TRUE,nrow(tab))
  # diff  <- rep(0,nrow(tab))
  # success <- in.front
  # m2 <- list()
  # m2$vb <- t(as.matrix(cbind(P2[fd,],1)))
  # m2$normals <- t(as.matrix(cbind(-tab[fd, c("ux","uy","uz")],1)))
  # class(m2) <- "mesh3d"
  # d <- vcgRaySearch(m2,mesh.test$mesh)
  # 
  # diff[fd] <- abs(tab$d[fd])-d$distance
  # success[fd] <- d$quality ==1
  # f.nok <-(diff<=-tol |  !success | (diff>=tol & diff<=1e-2))
  # in.front[!f.nok] <- abs(diff[!f.nok]) < tol
  # if (any(f.nok))
  #   in.front[f.nok]<- .meshdist ( pt_x = as.numeric(pt[,1]),
  #                                 pt_y = as.numeric(pt[,2]),
  #                                 pt_z = as.numeric(pt[,3]),
  #                                 p2_x = as.numeric(P2[f.nok,1]),
  #                                 p2_y = as.numeric(P2[f.nok,2]),
  #                                 p2_z = as.numeric(P2[f.nok,3]),
  #                                 u_x = -as.numeric(u[f.nok,1]),
  #                                 u_y = -as.numeric(u[f.nok,2]),
  #                                 u_z = -as.numeric(u[f.nok,3]),
  #                                 n_A = as.numeric(it[,1])-1,
  #                                 n_B = as.numeric(it[,2])-1,
  #                                 n_C = as.numeric(it[,3])-1)
  # tab$in.front <- in.front
  # 
  # 
  # if (plot) {
  #   palet <- colorRampPalette(c(col.in, "grey90",col.out))(100)
  #   
  #   dummy <- factor (sign(tab$d), levels=c(-1, 1), labels=c ("out", "in"))
  #   dev.new (width = 7, height = 5, noRStudioGD = TRUE)
  #   boxplot (d~in.out, data.frame (d= abs(tab$d), in.out=dummy), col=c(col.out, col.in),
  #            xlab = paste(mesh.test$object.alias,"vs",mesh.ref$object.alias ), ylab = "d (mm)")
  #   
  #   ## affichage des mesures
  #   d.max <- max (abs(tab$d))
  #   open3d (windowRect = c (100, 100, 600, 600))
  #   f <- (tab$in.front==1) 
  #   wire3d(mesh.test$mesh, color="grey90", lit= FALSE)
  #   points3d(as.matrix(tab[f,2:4]), col=palet[1 + 99 * (tab$d[f] / d.max / 2 + 0.5)], size=pt.size, lit = FALSE)
  #   display.palette(palet, breaks=seq (-d.max, d.max, length.out=101))
  #   dev.new (width = 7, height = 5, noRStudioGD = TRUE)
  #   H <- hist (tab$d, breaks=seq (-d.max, d.max, length.out=100))
  # }
  # return (tab)

