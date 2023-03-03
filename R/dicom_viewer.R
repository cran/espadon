#' DICOM content viewer
#' @description the \code{dicom.viewer} function displays the data of a DICOM file.
#' @param dcm espadon object of class "volume", "rtplan", "struct" provided by
#'  DICOM files, or DICOM filename, or Rdcm filename, or raw vector  representing 
#'  the binary extraction of the DICOM file.  
#' @param txt.sep String. Used if \code{as.txt = TRUE}. Separator of the tag value elements.
#' @param txt.length Positive integer. Used if \code{as.txt = TRUE}. Maximum number 
#' of letters in the representation of the TAG value.
#' @param tag.dictionary Dataframe, by default equal to \link[espadon]{dicom.tag.dictionary}, 
#' whose structure it must keep. This dataframe is used to parse DICOM files.

#' @param height,width Height and width in pixel of the DICOM table.
#' @param ... Additional argument \code{dicom.browser} when previously calculated by 
#' \link[espadon]{dicom.browser}. Argument \code{nb} or \code{dicom.nb} representing the 
#' number of DICOM file, when \code{dcm} contains multiple DICOM files.
#' @return Returns the DICOM file description in a browser window.
#' @seealso \link[espadon]{xlsx.from.dcm}, \link[espadon]{xlsx.from.Rdcm}, \link[espadon]{dicom.parser} 
#' @examples
#' if (interactive ()) dicom.viewer (toy.dicom.raw ())

#' @importFrom shiny fixedPanel fixedPage div reactiveVal observe observeEvent 
#' HTML stopApp isolate tags runGadget browserViewer
# shinyApp
#' @import DT
#' @importFrom shinyWidgets noUiSliderInput updateNoUiSliderInput setBackgroundColor
#' @import js 
#' @export


dicom.viewer <- function(dcm, txt.sep = "\\", txt.length = 100, 
                         tag.dictionary = dicom.tag.dictionary (),
                         height = 600, width = 900,...){
  fontsize <- 12  
  name <- ""
  # if (!is.raw(dcm) & !is.character(dcm[100])) return(NULL)
  # if (is.character(dcm)){ 
  #   name <- basename(dcm[1])
  #   dcm <- dicom.raw.data.loader (dcm[1]) }
  
  df <- dicom.parser(dcm, try.parse = TRUE,  tag.dictionary = tag.dictionary,
                     txt.sep = txt.sep, txt.length=txt.length)
  Encoding(df$Value) <- df$Value[df$TAG==("(0008,0005)")]#"UTF-8"
  
  if (is.null(df)) return (NULL)
  df <- df[!grepl("[(]FFFE,E00D[)]$|[(]FFFE,E0DD[)]$",df$TAG),]
  rownames(df) <- NULL
  
  df_ <- do.call(rbind,lapply(strsplit(df$TAG," "), function(l) {
    le <- length(l)
    return(data.frame(encaps=le-1,
                      TAG = HTML(paste(rep("  ",le-1),collapse=""),
                                 l[le])))}))
  df_$VM <- df$VM
  df_$VR <- df$VR
  
  df_$Value <- df$Value
  TAG <- paste0("^",gsub(")","[)]",gsub("(","[(]",df$TAG, fixed=T), fixed=T))
  f <- c(df_$encaps[1:(nrow(df_)-1)]+1 == df_$encaps[2:nrow(df_)], FALSE)
  df_$deploy <- " "
  df_$deploy[f] <- ">"
  
  df_ <- df_[,c(2,4,6,3,5,1)]
  df_$display <- df_$encaps==0
  df_$ttag <- df$TAG
  if ("(0008,0005)" %in% df_$ttag) 
    Encoding(df_$Value) <- df_$Value[df_$ttag==("(0008,0005)")]#"UTF-8"
  
  gl <- grepl(",0000[)]$",df$TAG)
  if (any(gl)) {
    f <- df_[gl,]$VM==""
    if (any(f))  df_[gl,]$VM[f] <- "Group Length"
  }
  gl <- df_$VR=="UN"
  if (any(gl)) df_[gl,]$VM <- "Private TAG"
  
  
  
  
  
  nbrow <- floor(height/fontsize/2.4)
  fontjs <- paste0("function(settings, json) {",
                   "$('body').css({'font-family': 'Consolas'});",
                   "$(this.api().table().header()).css({'font-size': '", fontsize ,"px'});",
                   "$(this.api().tables().body()).css({'font-size': '", fontsize ,"px'});}"
  )
  
  ro <- 1
  colwidth <-apply(sapply(df_[,c("TAG","deploy","VM","VR","Value")],nchar),2,max) 
  df_$TAG <- paste(df_$TAG,sapply(nchar(df_$TAG) ,function(i) paste(rep(" ",colwidth[1]-i),collapse="" )))
  df_$VM <- paste(df_$VM,sapply(nchar(df_$VM) ,function(i) paste(rep(" ",colwidth[3]-i),collapse="" )))
  df_$VR<- paste(df_$VR,sapply(nchar(df_$VR) ,function(i) paste(rep(" ",colwidth[4]-i),collapse="" )))
  df_$Value  <- paste(df_$Value ,sapply(nchar(df_$Value ) ,function(i) paste(rep(" ",colwidth[5]-i),collapse="" )))
  
  
  ui <-fixedPage(
    setBackgroundColor("white"),
    
    fixedPanel(
      top=25,
      left=10,
      width = 40,
      height = height,
      noUiSliderInput(inputId="sI", label=NULL,min=1, max= sum(df_$display)-nbrow,
                      # range=c(1,sum(df_$display)-nbrow),
                      value=1,  step=1,  update_on="change",
                      tooltips =FALSE, orientation="vertical",
                      color="#D0D0D0",
                      height = paste0(round((nbrow+1)*fontsize*2.1),"px"),
                      width = "5px"             # width = "300px", margin = 100,
      )),
    fixedPanel(
      left=50,
      width = width- 60,
      height = height,
      tags$style("#df_data { white-space:pre; }"),
      
      DTOutput("df_data"))
  )
  
  
  
  server <- function(input, output, session) {
    session$onSessionEnded (function() {stopApp()})
    
    disp_df <- reactiveVal(df_)
    # db <- reactiveVal(db_)
    maxrow <- reactiveVal(sum(df_$display)-nbrow)
    observeEvent(maxrow(),{
      m <- isolate(maxrow())
      updateNoUiSliderInput(
        session = session,
        inputId = "sI",
        range = c(1, m)
      )
    })
    
    
    output$df_data <-  renderDT ({
      db <- disp_df()
      ro <- input$sI
      db <- db[db$display,c("TAG","VR","deploy","VM","Value")][ro:(ro+ nbrow),]
      colnames(db) <- rep("",5)#c("TAG","VR","","VM","Value")
      
      hot <- datatable (db,  class ="display compact",
                        options = list(dom = 't',
                                       class ='white-space: pre',
                                       paging = FALSE,
                                       escape  = FALSE,
                                       ordering = FALSE,
                                       autoWidth = TRUE,
                                       scrollX = TRUE,
                                       
                                       initComplete = JS(fontjs)
                        ),
                        rownames=FALSE,
                        selection=list(mode="single", target="cell"))%>%
        formatStyle(1, fontWeight = 'bold')%>%
        formatStyle(4, target='row', backgroundColor = styleEqual("Private TAG", "#E2E2E2"))%>%
        formatStyle(5,  color = 'blue')
      
      
      return (hot)})
    
    
    observe({
      cell <- input$df_data_cell_clicked
      
      if (length(cell)!=0) {
        if (cell$col[1]==2 & cell$value[1]!=" "){
          db <- isolate (disp_df())
          ro <- isolate(input$sI)
          tg <-(db[db$display,"ttag"][ro:(ro+ nbrow)][cell$row[1]])
          idx.row <- which(db$ttag == tg )
          
          if (db$deploy[idx.row]==">"){
            flag <- grepl(paste0(TAG[idx.row],"[ ]"), db$ttag) & db$encaps ==(db$encaps[idx.row]+1)
            db$display[flag] <- TRUE
            db$deploy [idx.row] <- "v"
            
          } else if(db$deploy[idx.row]=="v"){
            flag <- grepl(paste0(TAG[idx.row],"[ ]"), db$ttag) & (db$encaps > db$encaps[idx.row])
            db$display[flag] <- FALSE
            db$deploy [flag & db$deploy!=" "] <- ">"
            db$deploy [idx.row] <- ">"
            
          }
          maxrow(sum(db$display)-nbrow)
          disp_df(db)
        }
        
      }
    })
    
  }
  
  # shinyApp (ui, server)
  # runGadget(ui, server, viewer  = dialogViewer(name, width = width, height = height))
  runGadget(ui, server, viewer = browserViewer(browser = getOption("browser")))
}