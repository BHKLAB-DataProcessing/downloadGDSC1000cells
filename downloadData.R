install.packages("readxl")
library("readxl")
require(downloader)

tmpdir <- tempdir()

cellFileURL <- "ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-6.0/"
cellFileName <- "Cell_Lines_Details.xlsx"

## TODO:: Make it download into tmpdir first

## download sample information
message("Download cell info")
myfn <- file.path(tmpdir, "gdsc1000_cellinfo.xlsx")

dwl.status <- download.file(url=sprintf("%s/%s",cellFileURL,cellFileName), destfile=file.path(tmpdir,cellFileName), quiet=TRUE)
if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }


require(gdata)
cell.info <- as.data.frame(read_excel(file.path(tmpdir, "gdsc1000_cellinfo.xlsx"),  sheet=1 ))
cell.info <- cell.info[-nrow(cell.info),]
cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
cellcuration <- cell_all[,c("CGP.cellid", "GDSC.SNP.cellid", "CGP_EMTAB3610.cellid", "unique.cellid")]
EMTAB3610_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"CGP_EMTAB3610.cellid"])))
SNP_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"GDSC.SNP.cellid"])))
CGP_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"CGP.cellid"])))
unique_matches <- match(toupper(gsub(pattern=badchars, "", x=cell.info$Sample.Name)), toupper(gsub(pattern=badchars, "", x=cellcuration[,"unique.cellid"])))

matches <- cbind(EMTAB3610_matches, SNP_matches, CGP_matches, unique_matches)

doubleMatchedCells <- apply(matches,1,function(x){
  return(ifelse(!all(is.na(x)), !all(c(x[1]==x[2], x[1]==x[3], x[1]==x[4]), na.rm=T), FALSE))
})

matchesNames <- cbind(cell_all$unique.cellid[EMTAB3610_matches], cell_all$unique.cellid[SNP_matches], cell_all$unique.cellid[CGP_matches], cell_all$unique.cellid[unique_matches])

##### for whatever reason, the SNP_matches are correct while the others are wrong except for JUKRAT which the first column in correct
matches <- apply(matches, 1, function(x){
  if(!is.na(x[4])){
    return(x[4])   
  } else {
    return(ifelse(length(unique(na.omit(x[1:3]))), unique(na.omit(x[1:3])), NA))
  }
})
cell.info$cellid <- cell.info$Sample.Name
cell.info$cellid[!is.na(matches)] <- cell_all$unique.cellid[na.omit(matches)]
cell.info$cellid[grep("KM-H2", cell.info$cellid)] <- cell.info$Sample.Name[grep("KM-H2", cell.info$cellid)]
cell.info$cellid[grep("^T[-]T$", cell.info$Sample.Name)] <- "T.T"
cell.info$cellid[grep("^TT$", cell.info$Sample.Name)] <- "TT"
save(cell.info, file="/pfs/out/cellInfo.RData")