#' Display patient objects links
#' @description The \code{display.obj.links} function displays a graph of 
#' connections between objects of a patient.
#' The name of the objects corresponds to their modality (ct, mr, rtdose...) 
#' followed by their position in their respective lists in the patient list objects.
#' Connected objects are linked by arrows. Objects sharing the same frame of reference 
#' have the same color except for objects with warnings, errors or missing planes 
#' which are all in grey.
#' Approved objects are circled in green. 
#' By default, objects shapes are circles, except rtdose represented as squares.
#' @param pat "patient" class object, as loaded using \link[espadon]{load.patient.from.dicom}, 
#'  \link[espadon]{load.patient.from.Rdcm} or \link[espadon]{toy.load.patient}.
#' @param obj.selected Dataframe (default to NULL) containing the objects already selected, 
#' created by a previous call of \code{display.obj.links} for example.
#' @param exclusion Vector of patient file modalities that should not be displayed, 
#' as for instance \code{"mr"}...
#' @param square Vector of patient file modalities that should be enclosed by a 
#' square, as for instance \code{c ("ct", "mr")}...
#' If \code{NULL} no object name is squared.
#' @param group.by.connected.FoR Boolean. If \code{TRUE} (default), all objects 
#' sharing the same frame of reference or connected by a registration matrix have 
#' the same color. If \code{group.by.connected.FoR =FALSE}, only objects sharing 
#' the same FoR have the same color.
#' @param interactive Boolean. If \code{interactive = TRUE}, buttons are available 
#' on the graph to get information about the objects and select or remove them from 
#' the data frame of the selected objects.
#' Then simply click on the name of the object on which to apply the chosen action. 
#' If \code{interactive = FALSE} no interaction possible with the plot.
#' @param random.seed Positive Integer or \code{NULL}. If \code{random.seed = NULL},
#' the objects are laid out randomly. The layout is otherwise fixed.

#' @return The function displays all patient objects, linked by an 
#' arrow when they are connected or a line when they belongs to the same DICOM 
#' object, and with a color and a shape depending on \code{square}, \code{group.by.connected.FoR}.
#' @return When \code{interactive = TRUE}, it returns a dataframe of the selected objects, 
#' or NULL if no object is selected.
#' @examples
#' # loading of toy-patient objects
#' patient <- toy.load.patient (dxyz = c (5, 5, 5), beam.nb = 1)
#' display.obj.links (patient, group.by.connected.FoR = FALSE)
#' display.obj.links (patient, group.by.connected.FoR = TRUE)
#' display.obj.links (patient, group.by.connected.FoR = TRUE, random.seed=NULL)

#' @seealso \link[espadon]{get.obj.connection}

#' @importFrom graphics locator
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph norm_coords layout_nicely
#' @importFrom grDevices dev.cur dev.list
#' @export
display.obj.links <- function (pat, obj.selected = NULL, exclusion = NULL, 
                               square = "rtdose", group.by.connected.FoR = TRUE, 
                               interactive = FALSE, random.seed = 314) {
  # initialisation
  old.lb <- getOption("locatorBell")
  old.par <- par(no.readonly = TRUE)
  my.dev <- dev.cur()
  if (exists(".Random.seed")){ rs <- .Random.seed
  } else rs <- NA

  on.exit(
    expr = {
      if (my.dev %in% dev.list()) (par(old.par))
      options(locatorBell = old.lb)
      if (!is.na(rs[1]))  .Random.seed <- rs
    })
  

  options(locatorBell = FALSE)
  if (is.null(obj.selected)) obj.selected <- data.frame()
  if (!is.null(random.seed)) set.seed(random.seed)
  ## Plot function
  F_Plot <- function (pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, 
                      obj.list, obj.selected, choice) {
    
    if (group.by.connected.FoR) {
      M.T.MAT <- matrix (pat$T.MAT$matrix.description$type != '', nrow=sqrt (length (pat$T.MAT$matrix.description$type)))
      src.T.MAT <- matrix (pat$T.MAT$matrix.description$src, nrow=sqrt (length (pat$T.MAT$matrix.description$type)))
      codes.T.MAT <- apply (M.T.MAT, 1, function (V) {
        sum (2^(1:length (V))*V)
      })
      ncolors <- length (unique (codes.T.MAT))
      obj.code <- sapply (obj.ref, function (ref) {
        V <- M.T.MAT[src.T.MAT == ref]
        sum (2^(1:length (V))*V)
      })
      color.levels <- as.numeric (factor (obj.code))
    } else {
      color.levels <- as.numeric (factor (obj.ref))
      ncolors <- max (color.levels)
    }
    color.levels[obj.be.careful] <- max (color.levels) + 1
    
    if (!is.null (square)) {
      to.square <- sapply (square, function (keep) grepl (paste0 ("^", keep), obj.list))
    } else to.square <- matrix (rep (FALSE, length (obj.list)), nrow=length (obj.list))
    to.square <- apply (to.square, 1, function (resp) any (resp))
    
    bold <- vector("logical",length(obj.alias))
    if (dim(obj.selected)[1]!=0){
      for (i in 1:dim(obj.selected)[1]){
        bold <- bold | sapply(as.matrix(obj.alias), function (alias) (alias==obj.selected$PIN.exam[i]))
      }
    }
    par(mar=c(0,0,3,0))
    plot(range(coords[,1])*1.5, range(coords[,2])*1.5, type="n", main=pat$pat.pseudo)
    
    plot.igraph (network, layout = coords, add=TRUE,
                 vertex.color = c(hcl.colors (ncolors, "Pastel 1"), "#b0b0b0")[color.levels],
                 vertex.shape = c("circle", "square")[to.square + 1],
                 vertex.frame.color = c(NA, "green")[(obj.approv == "APPROVED") + 1],
                 vertex.label.color = "black",
                 vertex.label.font = c(1,4)[bold +1],
                 vertex.label.cex = c(0.9,1)[bold +1],
                 edge.arrow.size =0.7, edge.color="black")
    
    plot.igraph (network.same.obj, layout = coords, add=TRUE,
                 vertex.color = c(hcl.colors (ncolors, "Pastel 1"), "#b0b0b0")[color.levels],
                 vertex.shape = c("circle", "square")[to.square + 1],
                 vertex.frame.color = c(NA, "green")[(obj.approv == "APPROVED") + 1],
                 vertex.label.color = "black",
                 vertex.label.font = c(1,4)[bold +1],
                 vertex.label.cex = c(0.9,1)[bold +1],
                 edge.arrow.size = 0, edge.color="grey80")
    
    
    if (choice == "Get infos") {
      text(0, par()$usr[4]-0.08, "GET INFOS")
      text(par()$usr[1]+0.15, par()$usr[4]-0.08, "QUIT")
      text(0, par()$usr[3]+0.08, cex=0.9, font=3, "Click on an exam to display informations.")
      rect (par()$usr[1]+0.03, par()$usr[4]-0.019, par()$usr[1]+0.27, par()$usr[4]-0.141)
    } else if (choice == "Select exams") {
      text(0, par()$usr[4]-0.08, "Select the exams to be kept")
      text(par()$usr[1]+0.15, par()$usr[4]-0.08, "QUIT")
      text(0, par()$usr[3]+0.08, cex=0.9, font=3, "Click on an exam to add it at the data frame.")
      rect (par()$usr[1]+0.03, par()$usr[4]-0.019, par()$usr[1]+0.27, par()$usr[4]-0.141)
    } else if (choice == "Remove exams") {
      text(0, par()$usr[4]-0.08, "Remove from the data frame")
      text(par()$usr[1]+0.15, par()$usr[4]-0.08, "QUIT")
      text(0, par()$usr[3]+0.08, cex=0.9, font=3, "Click on an exams to remove it on the data frame.")
      rect (par()$usr[1]+0.03, par()$usr[4]-0.019, par()$usr[1]+0.27, par()$usr[4]-0.141)
    }
  }
  
  ## Action choice function
  F_Select_Action <- function(names, choice) {
    F_Bouton <- function (name, number)
    {
      number = number
      startx = 0.02
      endx = 0.18
      starty = 1.0 - number*0.06
      endy = starty - 0.04
      rect (startx, starty, endx, endy)
      text (startx+(endx-startx)/2, starty-(starty-endy)/2, name, cex=0.8, adj=c(0.5, 0.5))
      return(data.frame(name=name, x=startx+(endx-startx)/2, y=starty-(starty-endy)/2, number=number))
    }
    tryCatch(error = function(e) {choice<-"Quit";return(choice)},{
      par (usr = c(0, 1, 0, 1))
      boutons <- list()
      number <- 1:length(names)
      for (i in 1:length(names)) {boutons[[i]] <- F_Bouton(names[i], number[i])}
      boutons <- do.call(rbind, boutons)
      coords <- cbind(boutons[[2]], boutons[[3]])
      
      repeat {
        B <- unlist(locator (n=1))
        dist <- apply (coords, 1, function (Mi) sum ((B - Mi)^2))
        idx <- which.min (dist)
        dist.x <- abs(B[1] - coords[idx,1])
        dist.y <- abs(B[2] - coords[idx,2])
        if (dist.x <= 0.089 & dist.y <= 0.0265) {
          cat ("\n*****", boutons[[1]][[idx]], "*****\n")
          return(boutons[[1]][[idx]])
        }
      }
    })
  }
  
  ## Information function
  F_Display_Infos <- function (pat, coords, obj.type, obj.idx, choice) {
    tryCatch(error = function(e) {choice<-"Quit";return(choice)},{
      xy <- rbind(coords,c(par()$usr[1]+0.15, par()$usr[4]-0.08))
      repeat {
        M <- unlist (locator (n=1))
        dist <- apply (xy, 1, function (Mi) sum ((M - as.numeric (Mi))^2))
        idx <- which.min (dist)
        dist.x <- abs(M[1] - xy[idx,1])
        dist.y <- abs(M[2] - xy[idx,2])
        if (dist.x <= 0.132 & dist.y <= 0.12) {
          if (idx == dim(xy)[1]) {cat("\n***** quit *****\n"); break}
          sel.type <- names (pat)[obj.type[idx]]
          sel.idx <- obj.idx[idx]
          cat ("\n\t*****", sel.type, sel.idx, "*****\n")
          
          if (sel.type == "rtdose") {
            D <- pat$rtdose[[sel.idx]]
            cat ("object      :", D$object.alias, "\n")
            cat ("study date  :", D$study.date, "\t||\tacq date :", D$acq.date, "\t||\tcreation date :", D$creation.date, "\n")
            cat ("description :", D$description, "\n")
            cat ("status      :", D$approval.status, "\n")
            cat ("FoR         :", D$ref.pseudo, "\n")
            cat ("Dmax        :", D$max.pixel)
          } else if (sel.type == "rtstruct") {
            S <- pat$rtstruct[[sel.idx]]
            cat ("object      :", S$object.alias, "\n")
            cat ("study date  :", S$study.date, "\t||\tacq date :", S$acq.date, "\t||\tcreation date :", S$creation.date, "\n")
            cat ("description :", S$description, "\n")
            cat ("status      :", S$approval.status, "\n")
            cat ("FoR         :", S$ref.pseudo, "\n")
            cat ("roi number  :", S$nb.of.roi, "\n")
            cat (paste (sort (S$roi.info$name), collapse=" | "), "\n")
          } else if (sel.type == "ct" || sel.type == "mr" || sel.type == "pt") {
            if (sel.type == "ct") IM <- pat$ct[[sel.idx]]
            else if (sel.type == "mr") IM <- pat$mr[[sel.idx]]
            else IM <- pat$pt[[sel.idx]]
            cat ("object      :", IM$object.alias, "\n")
            cat ("study date  :", IM$study.date, "\t||\tacq date :", IM$acq.date, "\t||\tcreation date :", IM$creation.date, "\n")
            cat ("description :", IM$description, "\n")
            cat ("FoR         :", IM$ref.pseudo, "\n")
            cat ("n.ijk       :", IM$n.ijk)
          } else if (sel.type == "rtplan") {
            P <- pat$rtplan[[sel.idx]]
            cat ("object      :", P$object.alias, "\n")
            cat ("plan name   :", P$plan.info$plan.name, "\n")
            cat ("study date  :", P$study.date, "\t||\tacq date :", P$acq.date, "\t||\tcreation date :", P$creation.date, "\n")
            cat ("description :", P$description, "\n")
            cat ("status      :", P$approval.status, "\n")
            cat ("FoR         :", P$ref.pseudo, "\n")
            if (!(is.null(P$fraction.beam))) {
              cat ("beam nb     :", P$fraction.info$nb.of.beam, "\n")
              cat ("fraction nb :", P$fraction.beam$nb.of.frac.planned, "\n")
              cat ("beam dose   :", P$fraction.beam$beam.dose, "\n")
            }
            if (!(is.null(P$fraction.brachy))) {
              cat ("brachy nb   :", P$fraction.info$nb.of.brachy.app, "\n")
              cat ("fraction nb :", P$brachy.info$nb.of.frac.planned, "\n")
              cat ("brachy dose :", P$brachy.info$brachy.dose, "\n")
            }
          } else {
            # cat ("object not handled\n")
            obj <- pat[[obj.type[idx]]][[sel.idx]]
            cat ("object      :", obj$object.alias, "\n")
            cat ("study date  :", obj$study.date, "\t||\tacq date :", obj$acq.date, 
                 "\t||\tcreation date :", obj$creation.date, "\n")
            cat ("description :", obj$description, "\n")
            # cat ("status      :", obj$approval.status, "\n")
            cat ("FoR         :", obj$ref.pseudo, "\n")
          }
          cat ("\n==========================================================\n")
        }
      }
      return(choice)})
  }
  
  ## Object selection function
  F_Select_Exams <- function ( network, network.same.obj, obj.selected, pat, coords, obj.type, obj.idx, obj.alias, choice) {
    tryCatch(error = function(e) {choice<-"Quit";return(list(obj=obj.selected,ch=choice))},{
      F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, choice)
      xy <- rbind(coords,c(par()$usr[1]+0.15, par()$usr[4]-0.08))
      repeat {
        M <- unlist (locator (n=1))
        dist <- apply (xy, 1, function (Mi) sum ((M - as.numeric (Mi))^2))
        idx <- which.min (dist)
        dist.x <- abs(M[1] - xy[idx,1])
        dist.y <- abs(M[2] - xy[idx,2])
        if (dist.x <= 0.132 & dist.y <= 0.12) {
          if (idx == dim(xy)[1]) {cat("\n***** quit *****\n"); break}
          sel.type <- names (pat)[obj.type[idx]]
          sel.idx <- obj.idx[idx]
          cat ("\n*****", sel.type, sel.idx, "*****\n")
          keep <- data.frame(PIN.patient=pat$pat.pseudo, PIN.exam=obj.alias[idx], type=sel.type, idx=sel.idx)
          if (!(TRUE %in% (sapply(obj.selected$PIN.exam, function (alias) (alias==keep$PIN.exam))))){
            obj.selected <- rbind(obj.selected, keep)
          }
          F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, choice)
        }
      }
      return(list(obj=obj.selected,ch=choice))})
  }
  
  ## Object remove function
  F_Remove_Exams <- function ( network, network.same.obj, obj.selected, pat, coords, obj.type, obj.idx, obj.alias, choice) {
    tryCatch(error = function(e) {choice<-"Quit";return(list(obj=obj.selected,ch=choice))},{
      F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, choice)
      xy <- rbind(coords,c(par()$usr[1]+0.15, par()$usr[4]-0.08))
      repeat {
        M <- unlist (locator (n=1))
        dist <- apply (xy, 1, function (Mi) sum ((M - as.numeric (Mi))^2))
        idx <- which.min (dist)
        dist.x <- abs(M[1] - xy[idx,1])
        dist.y <- abs(M[2] - xy[idx,2])
        if (dist.x <= 0.132 & dist.y <= 0.12) {
          if (idx == dim(xy)[1]) {cat("\n***** quit *****\n"); break}
          sel.type <- names (pat)[obj.type[idx]]
          sel.idx <- obj.idx[idx]
          cat ("\n*****", sel.type, sel.idx, "*****\n")
          remove.PIN <- obj.alias[idx]
          remove.idx <- which(obj.selected$PIN.exam==remove.PIN)
          if (length(remove.idx)!=0) (obj.selected <- obj.selected[-remove.idx,])
          F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, choice)
        }
      }
      return(list(obj=obj.selected,ch=choice))})
  }
  
  ## Main
  ## Initialisation
  if (is.null(pat)) stop("pat is NULL")
  if (any(is.na(match(c("patient", "pat.pseudo","description","T.MAT"),names(pat))))) stop("pat is not a patient")
  obj.list <- obj.alias <- obj.name <- obj.ref <- obj.nb <- obj.approv <- c()
  obj.type <- obj.idx <- obj.be.careful <- c()
  obj.con <- list ()
  for (i in (which(names(pat)=="T.MAT")+1):length (pat)) {
    toggle = TRUE
    if (names (pat)[i] %in% exclusion) toggle <- FALSE
    if (toggle) {
      for (j in 1:length (pat[[i]])) {
        obj.type <- c(obj.type, i)
        obj.idx <- c(obj.idx, j)
        obj.list <- c (obj.list, paste(names (pat)[i], j))
        obj.alias <- c (obj.alias, pat[[i]][[j]]$object.alias)
        obj.name <- c (obj.name, pat[[i]][[j]]$object.name)
        obj.ref <- c (obj.ref, pat[[i]][[j]]$ref.pseudo)
        obj.nb <- c (obj.nb, pat[[i]][[j]]$number)
        obj.be.careful <- c(obj.be.careful,
                            !is.null (pat[[i]][[j]]$error) | !is.null (pat[[i]][[j]]$warning) |
                              ifelse (is.null (pat[[i]][[j]]$missing.k.idx), FALSE, pat[[i]][[j]]$missing.k.idx))
        if (is.null (pat[[i]][[j]]$approval.status)) obj.approv <- c (obj.approv, "")
        else obj.approv <- c (obj.approv, pat[[i]][[j]]$approval.status)
        if (is.null (pat[[i]][[j]]$ref.object.alias)) obj.con <- c(obj.con, "NA")
        else obj.con <- c(obj.con, list (pat[[i]][[j]]$ref.object.alias))
      }
    }
  }
  
  M <- matrix (0, nrow=length (obj.list), ncol=length (obj.list))
  rownames (M) <- obj.list
  colnames (M) <- obj.list
  
  for (i in 1:length (obj.list)) {
    idx <- which (sapply (obj.con, function (V) obj.alias[i] %in% V))
    M[i, idx] <- 1
    # idx <- which (sapply (obj.name, function (V) obj.name[i] %in% V))
    # M[i, idx] <- 1
  }
  # M[cbind (1:length (obj.list), 1:length (obj.list))] <- 0
  
  
  M.same.obj <- matrix (0, nrow=length (obj.list), ncol=length (obj.list))
  rownames (M.same.obj) <- obj.list
  colnames (M.same.obj) <- obj.list
  for (i in 1:length (obj.list)) {
    idx <- which (sapply (obj.name, function (V) obj.name[i] %in% V))
    M.same.obj[i, idx] <- 1
  }
  
  M.same.obj[cbind (1:length (obj.list), 1:length (obj.list))] <- 0
  
  network <- graph_from_adjacency_matrix(M)
  network.same.obj <- graph_from_adjacency_matrix(M.same.obj)
  
  coords <- norm_coords(layout_nicely(network))
  boutons <- c("Get infos","Select exams","Remove exams","Quit")
  choice <- ""
  
  ## Start repetition display until quit
  if (interactive){
    dev.new()
    my.dev <- dev.cur ()
    repeat {
      if (choice=="Quit") break
      F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, "")
      choice <- F_Select_Action(boutons, choice)
      if (choice == "Get infos")  {
        F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, choice)
        choice <- F_Display_Infos(pat, coords, obj.type, obj.idx, choice)
      } 
      if (choice == "Select exams") {
        r <-  F_Select_Exams( network, network.same.obj, obj.selected, pat, coords, obj.type, obj.idx, obj.alias, choice)
        obj.selected <- r$obj
        choice<-r$ch
      }  
      if (choice == "Remove exams") {
        r <- F_Remove_Exams( network, network.same.obj, obj.selected, pat, coords, obj.type, obj.idx, obj.alias, choice)
        obj.selected <- r$obj
        choice<-r$ch
      }
    }
    if (nrow(obj.selected)!=0){
      obj.selected <- obj.selected[order(obj.selected$PIN.patient),]
      return(obj.selected)
    }
    return (NULL)
  } 
  F_Plot(pat, network, network.same.obj, coords, obj.approv, obj.be.careful, obj.ref, obj.list, obj.selected, "")
  
}

