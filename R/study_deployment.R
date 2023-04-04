#' Deployment of DICOM files from multiple patients
#' @description The \code{study.deployment} function duplicates DICOM data from 
#' multiple patients, so that it becomes data independent of the original data. 
#' This function simplifies the analysis of multi-center or multi-expert studies 
#' in dosimetry challenges, contouring consensus searches, etc.
#' @param pats.dir Name of the directory in which all patient directories are 
#' stored, each containing the DICOM files to be duplicated.
#' @param deploy.dir Name of the directory where all patient files will be 
#' duplicated.
#' @param design.matrix Boolean matrix. See Details.
#' @param pid.prefix string vector of length 1 or string vector of length 
#' \code{ncol(design.matrix)}, representing the prefix added to the new 
#' unique identifier of the deployed patient (tag (0010,0020)).
#' @param white.list Names vector, representing a part of the DICOM tag name
#' UI value representation, other than those defined by the DICOM 
#' standard, which will be modified. By default, the UID name containing 'instance' 
#' or 'reference' will be modified.
#' @param black.list Names vector, representing a part of the DICOM tag name
#' UI value representation, other than those defined by the DICOM 
#' standard, which will not be modified. By default, the frame of reference UID 
#' will not be modified.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.

#' @details 
#' The \code{design.matrix} argument defines how patients DICOM files will be deployed.
#' The names of the lines must match the names of the directories contained in \code{pats.dir}.
#' The names of the columns are for example the different experts or hospitals 
#' who will study the patient files. These experts will only review the patients 
#' files defined by \code{rownames(design.matrix)[design.matrix[,"expert"]]}.
#' @return Creates the \code{deploy.dir} directory, containing the expert 
#' directories defined by the \code{design.matrix} column names. Each expert 
#' directory contains as many patient directories as defined by the 
#' \code{design.matrix} row names. All patients will be independent of each other.
#' The new created patients have the pats.dir as name, and expert name as first 
#' name, and an independent patient ID, with prefix \code{pid.prefix}.
#' @examples
#' # First, save toy.dicom.raw () raw data to a temporary file/pats.dir/toy_PM 
#' # for testing.
#' toy_PM.dir <- file.path (tempdir(), "pats.dir","toy_PM") 
#' dir.create (toy_PM.dir, recursive = TRUE) 
#' dcm.filename <- tempfile (pattern = "toyrtplan", tmpdir = toy_PM.dir,
#'                           fileext = ".dcm")
#' zz <- file (dcm.filename, "wb")
#' writeBin (toy.dicom.raw (), zz, size = 1)
#' close (zz)
#' 
#' # function test:
#' pats.dir <- dirname (toy_PM.dir)
#' deploy.dir <- file.path (tempdir(), "deploy.dir") 
#' design.matrix <- matrix(TRUE, nrow = length (dir (pats.dir)), ncol=3,
#'                         dimnames = list (basename (dir (pats.dir)),
#'                                          c("Dr Quinn","Dr Who","Dr House")))
#' design.matrix
#' study.deployment (pats.dir, deploy.dir, design.matrix, 
#'                  pid.prefix = c("zz_", "yy_", "xx_"))
#' 
#' # check result
#' list.files(deploy.dir, recursive = TRUE)
#' load.patient.from.dicom(deploy.dir)$patient

#' # Cleaning  temporary directory
#' unlink (pats.dir, recursive = TRUE)
#' unlink (deploy.dir, recursive = TRUE)
#' @export
#' @import progress
study.deployment <- function (pats.dir, deploy.dir, 
                              design.matrix = matrix(TRUE,nrow = length (dir(pats.dir)), ncol=1,
                                                     dimnames = list(basename (dir (pats.dir)),"expert_1")),
                              pid.prefix = "",
                              white.list = c("instance","reference"), 
                              black.list = c("frame of reference","class"), 
                              tag.dictionary = dicom.tag.dictionary()) {
  
  if (!dir.exists(pats.dir)) stop(paste0("'",pats.dir, "' does'nt exist."))
  dir.pat.list <- dir(pats.dir)
  if (length(dir.pat.list) == 0) stop(paste0("'",pats.dir, "' contains no patients to be deployed."))
  if (is.null(colnames(design.matrix))) 
    stop("design.matrix should have column names corresponding to the names of the experts who will study the patients")
  m <- is.na(match(rownames(design.matrix),basename (dir(pats.dir))))
  if (all(m)) 
    stop("Check the row names in design.matrix or the contents of pats.dir: the row names should match the patient directories in pats.dir.")
  if (any(m)) warning("Some of patients do not exist in pats.dir.")
  pat.list <- row.names (design.matrix)
  expert.list <- colnames (design.matrix) 
  

  if (any(!m)){ 
    pat.list <- pat.list[!m]
    design.matrix <-  as.matrix(do.call(rbind.data.frame,lapply(rownames(design.matrix)[!m], function(r) design.matrix[r,])))
    row.names (design.matrix) <- pat.list 
    colnames (design.matrix) <- expert.list
  }
  
  f <- design.matrix==1
  nb <- sum(f)
  if (nb==0) stop("No distribution request")
  
  
  pat.ID <- paste0( "999",as.integer (runif (1, 10000, 99999-nb)) + c(1:nb))
  pat.ID.matrix <- design.matrix
  pat.ID.matrix[] <- ""
  pat.ID.matrix[f] <-pat.ID
  
  if ((ncol(design.matrix)!= length (pid.prefix)) & length (pid.prefix)!=1) 
    stop ("pid.prefix length must be a vector of length 1 or a vector of length ncol(design.matrix)")
  
  if (length (pid.prefix)==1) pid.prefix <- rep(pid.prefix, ncol(design.matrix))
    
  pat.pseudo.matrix <- matrix(as.character(sapply(row.names(design.matrix), function(rn) 
    paste(rn,colnames(design.matrix), sep = "^"))), ncol = ncol(design.matrix), byrow = TRUE)
  pat.pseudo.matrix[!f] <- ""
  
  SOPclassUID <- function(dicom.raw.data, dicom.df,i = grep( "^[(]0008,0016[)]$",dicom.df$tag)) {
    
    return (gsub("[.]","",sub("^1[.]2[.]840[.]10008","",
                              dicom.tag.parser (dicom.df$start[i], 
                                                dicom.df$stop[i], dicom.df$VR[i], 
                                                dicom.df$endian[i], dicom.raw.data))))
  }
  
  
  for (pl.idx in 1:length(pat.list)) {
    if (any(design.matrix[pl.idx,])){
      p <- pat.list[pl.idx]
      p.dir <- file.path (pats.dir, p)
      cat ("loading", p.dir, "\n")
      
      input <- list.files (file.path (pats.dir, p), full.names = FALSE, recursive = TRUE)
      pb <- progress_bar$new(format = " downloading [:bar] :percent",
                             total = length(input),  width= 60)
      mapping.UID.df <- data.frame(Value=character(), New=character())
      
      for (input.file in input) {
        
        dicom.raw.data <- dicom.raw.data.loader (file.path(pats.dir, p,input.file))
        dicom.df <-  dicom.browser (dicom.raw.data,   full.info = TRUE, 
                                    tag.dictionary = tag.dictionary)
        
        if (!is.null(dicom.df)){
          UID.tab0 <- .dicom.get.UID (dicom.raw.data = dicom.raw.data,
                                      white.list=white.list,
                                      black.list = black.list,
                                      tag.dictionary = tag.dictionary,
                                      dicom.browser= dicom.df)
          
          # pn met tous les tag de protocole 002 à la fin 
          f0002 <- grepl("^[(]0002,", UID.tab0$tag)
          UID.tab <- rbind(UID.tab0[!f0002, ],UID.tab0[f0002, ])
          rownames(UID.tab) <- NULL
          
          #on actualise la table mapping.UID.df
          UID.tab <- UID.tab[match(unique(UID.tab$Value),UID.tab$Value),]
          f.UID <- is.na(match(unique(UID.tab$Value),mapping.UID.df$Value))
          if (any(f.UID)) {
            UID.tab.to.add <- UID.tab[f.UID,]
            UID.tab.to.add$New = ""
            
            SOPclass <- SOPclassUID(dicom.raw.data, dicom.df)
            UID.tab.to.add$New[grepl("[(]0008,0014[)]$",UID.tab.to.add$tag)] <- .espadon.UID()
            UID.tab.to.add$New[grepl("[(]0008,0018[)]$",UID.tab.to.add$tag)] <- paste(
              .espadon.UID(),"patID", SOPclass,
              as.character(floor(as.numeric(Sys.time())*100)/100), sep=".")
            
            f.1155 <- grepl("[(]0008,1155[)]$",UID.tab.to.add$tag)
            if (sum(f.1155)>0) {
              SOP.i.1155 <- sapply(UID.tab.to.add$line[f.1155]-1,
                                   function(i ) SOPclassUID (dicom.raw.data, dicom.df,i))
              dum <- as.character((floor(as.numeric(Sys.time())*100) + 1:length(SOP.i.1155))/100)
              UID.tab.to.add$New[f.1155] <- paste(.espadon.UID(),
                                                  "patID",paste(SOP.i.1155,dum, sep="."), sep=".")
            }
            f <- which(UID.tab.to.add$New=="")
            dum <- as.character((floor(as.numeric(Sys.time())*100) + 1:length(f))/100)
            UID.tab.to.add$New[f] <- paste(paste(.espadon.UID(), "patID", sep="."), dum, sep=".")
            mapping.UID.df <- rbind(mapping.UID.df,UID.tab.to.add[, c("Value", "New")])
          }
          
          m.UID <- match(UID.tab0$Value,mapping.UID.df$Value)
          tag.vect <- c("(0010,0010)", "(0010,0020)",dicom.df$tag [UID.tab0$line])
          
          for (el.idx in 1:length(expert.list)) {
            pat.ID <- pat.ID.matrix[pl.idx, el.idx]
            pat.pseudo <- pat.pseudo.matrix[pl.idx, el.idx]
            e <- expert.list[el.idx]
            e.p.dir <- file.path (deploy.dir, e, p)
            if (nchar(pat.pseudo) != 0) {
              value.vect <- c (pat.pseudo, paste0(pid.prefix[el.idx], gsub("^999","",pat.ID)), 
                               gsub("patID",pat.ID,mapping.UID.df$New[m.UID]))
              
              dr <- dicom.set.tag.value (dicom.raw.data, tag=tag.vect, 
                                         tag.value =value.vect, 
                                         tag.dictionary = tag.dictionary, 
                                         dicom.browser = dicom.df)
              
              out.fn <- file.path(e.p.dir, input.file)
              out.fn <- file.path(dirname(out.fn), paste0("E",el.idx,"P",pl.idx,"_", basename (out.fn)))
              if (!dir.exists (dirname(out.fn))) dir.create (dirname(out.fn), recursive = TRUE)
              zz <- file (out.fn, "wb")
              writeBin (dr, zz, size = 1)
              close (zz)
            }
          }
        }
        pb$tick ()
      }
    }
  }
}
