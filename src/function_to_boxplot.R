#########################################
#   function_to_boxplot.R                #
#########################################
#
# Script to compute boxplots for the functions
# measurements
################################################
# AP-G. Spanish National Center for Biotechnology
# apascualgarcia.github.io
# Redited and cleaned, August, 2024
################################################

rm(list=ls())
# START EDITING --------------------

# ... File and directory with the functions, ASV table and metadata
fileOTU="seqtable_readyforanalysis.csv"
fileFun =  "20241205_Functions_ALL.csv" #"20151016_Functions_remainder.csv"
dirData = "6_finalfiles"
file.meta = "metadata_Time0D-7D-4M_May2022_wJSDpart-merged_ext.csv"
dirMeta =  "7.1_classes"

#... Directory for the output figures
dirFigs = "7.7_functions"

# ... A vector with the plots you want to create among the options displayed below
modes.plot = c("raw","norm","log","lognorm")
# "raw": Raw functions
# "norm": Functions normalized by cell count
# "log": Functions in logarithmic scale
# "lognorm": Logarithmic and normalized


# --- Packages & scripts

scripts <- c("clean_ASV_table.R") # "heatmap.2.mod.R") # list of functions to load

packages <- c("ggplot2", "ggpubr", "ggsignif")#,
#"GGally")# "mixtools") # list of packages to load

###### STOP EDITING -------------

# --- Set the main directory
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1]
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
#dirSrc=here::here() # src of the repository
setwd(this.dir)

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")

n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
setwd(dirSrc)
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}

# SOURCE FUNCTIONS ---------

n_scripts <- length(scripts) # count how many packages are required

for(n in 1:n_scripts){
  cat("Loading script #", n, " of ", n_scripts, "... Currently Loading: ", scripts[n], "\n", sep = "")
  lib_load <- paste("source(\"",scripts[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}

# --- Create some labels depending on your analysis
FunMean=c("CO2","ATP","G","N","X","P","C","Cells") # This is what we want to see in the labels
FunNorm=c("CO2/Cells","ATP/Cells","G/Cells","N/Cells","X/Cells", "P/Cells","C/Cells","Cells")
FunLog=c("log(CO2)","log(ATP)","log(G)","log(N)","log(X)", "log(P)","log(C)","log(Cells)")
FunLogNorm=c("log(CO2/Cells)","log(ATP/Cells)","log(G/Cells)","log(N/Cells)","log(X/Cells)", "log(P/Cells)","log(C/Cells)","log(Cells)")

# Load data ---------------

#... read in functional data
setwd(this.dir)
setwd(dirData)
dd.data=read.csv(fileFun, sep = "\t")
# ... read ASV table
ASV.table=read.table(fileOTU,sep="\t")
dim(ASV.table)
ASV.table=as.matrix(ASV.table)
setwd(this.dir)

# Load metadata and subset samples ----------
setwd(this.dir)
setwd(dirMeta)
sample_meta = read.csv(file.meta, sep = "\t")
setwd(this.dir)

# Clean data   ----------
# --- Extract a set with only time 7
exclude_exp = c("4M","0D")
clean.data.list = clean_ASV_table(ASV.table, sample_meta,
                                  match_exp = T, exclude_exp = exclude_exp)

ASV.table.final = clean.data.list$ASV.table
sample_meta.final = clean.data.list$sample_md

# --- Work on the functions
#exclude negative controls
dd.data=dd.data[dd.data$Community!="blank",]
dd.data=na.omit(dd.data)

# ... create new field for sampleid
dd.data$sampleid = paste(dd.data$Community, dd.data$Replicate, sep = ".")

# ... average few technical replicates
dd = dd.data
comm.id = paste(dd.data$Community,dd.data$Replicate, sep = '.')
dd=aggregate(dd.data[4:ncol(dd.data)],
              list(comm.id),function(x){mean(x,na.rm=T)})
colnames(dd)[1] = "sampleid"
rownames(dd) = dd$sampleid

# select the functions of interest
idxx=colnames(dd) # create a vector with the remaining names, select those containing "7" and convert it into a matrix
dd=subset(dd,select=grep("7", idxx))
cc=dd$CPM7
cc.log=log10(cc)
dd=subset(dd,select= -c(pgRPC.7,CPM7))


# Transform the data -------------------
# .... Change units
dd$mgCO2.7=dd$mgCO2.7*1000 # miligr to microgr
dd$ATP7=dd$ATP7/1000 # nanomol to micromol
newcol=dd$mG7+dd$mX7 # carbon uptake

dd=cbind(dd,newcol)

# .... Normalizations, etc.
dd.log=log10(dd)
dd.norm=dd/cc # normalize by cell count
dd.norm.log=log10(dd.norm) # Improve normality
dd=cbind(dd,cc)
dd.log=cbind(dd.log,cc.log)
dd.norm=cbind(dd.norm,cc)
dd.norm.log=cbind(dd.norm.log,cc.log)
rownames(dd.log)=rownames(dd)
rownames(dd.norm.log)=rownames(dd)
rownames(dd.norm)=rownames(dd)
colnames(dd.log)=FunLog #FunMean
colnames(dd.norm.log)=FunLogNorm
colnames(dd.norm)=FunNorm
colnames(dd)=FunMean


# ---- Identify communities with missing functions and reshape
matched = match(sample_meta.final$sampleid, rownames(dd))
not.found.final = unique(as.character(sample_meta.final$sampleid[is.na(matched)]))
#length(not.found.final) # 16 missed at least one replicate --> **fixed, labeling error in the plate reader**
#[1] "RDM01.2"  "WYC15.2"  "WYC41.2"  "WYC42.2"  "WYD09.2"  "WYM02.2" 
#[7] "WYM21.2"  "WYM25.2"  "WYM27.2"  "WYM29.2"  "WYM33.2"  "WYM39.2" 
#[13] "WYT116.2" "WYT63.2"  "WYT87.2"  "WYT94.2" 
sample_meta.final.matched = sample_meta.final[!is.na(matched),]
matched = matched[!is.na(matched)]
dd.final.matched = dd[matched, ]
dd.log.final.matched = dd.log[matched, ]
dd.norm.log.final.matched = dd.norm.log[matched, ]
dd.norm.final.matched = dd.norm[matched,]

# --- Define your partition and function dataset
Part  = as.character(sample_meta.final.matched$partition)
names(Part) = sample_meta.final.matched$sampleid
# table(Part)
# Class1 Class2 
# 794    306 

# Plot results -----------------

# .... Fix your color codes
colorCodes =c("Class1"="chocolate4","Class2"="chartreuse")
grouping=as.data.frame(Part)
colorVec=colorCodes[Part]
colorVec=as.data.frame(colorVec) # For ggplot, for plot comment
colnames(colorVec)="Color"
colnames(grouping)="Group"

# --- Select the dataset to plot
for(mode in modes.plot){
  if(mode == "norm"){
    xx=dd.norm.final.matched
    yin.array = c("CO2/Cells","ATP/Cells","X/Cells","G/Cells","N/Cells","P/Cells")
  }else if(mode == "log"){
    xx=dd.log.final.matched
    yin.array=c("log(Cells)","log(CO2)","log(ATP)",
                "log(X)","log(G)","log(N)","log(P)") # May 2024
  }else if(mode == "lognorm"){
    xx = dd.norm.log.final.matched
    yin.array = c("log(CO2/Cells)","log(ATP/Cells)",
                  "log(X/Cells)","log(G/Cells)","log(N/Cells)","log(P/Cells)")
  }else{
    xx=dd.final.matched
    yin.array = c("Cells","CO2","ATP",
                  "X","G","N","P")
  }
  ylab.array = yin.array
  nplots = length(yin.array)
  xy=cbind(xx,colorVec,grouping) # Here colorVec converts into a factor
  
  # --- Plot boxes
  xin="Group"
  xin_quo = enquo(xin)
  xlab="Community class"
  group= "Group" # "Color"#  "Group" #"Color"
  
  # --- Start plot
  setwd(this.dir)
  setwd(dirFigs)
  box=list()
  nbox = 7
  empty.plot=F
  if(nplots != nbox){
    empty.plot=TRUE #FALSE # if you want to create an empty plot for the first box
  }
  i=0
  for(k in 1:nbox){
    if(empty.plot == TRUE){
      if(k == 1){
        box[[k]]=ggplot() + theme_void()
        next
      }
    }
    i=i+1
    ylab=ylab.array[i]
    yin=yin.array[i]
    yin_quo = enquo(yin)
    xlabel = ""
    box[[k]]= ggboxplot(xy, xin, yin,
                        #geom_boxplot(xy, aes(xin, yin),
                        #combine=TRUE, # If several variables are used for y 
                        #color = "grouping", fill = "grouping",palette = Codes,
                        #color = "Color", 
                        #short.panel.labs = FALSE,
                        #panel.labs = c("1","2","3","4","5","6"),
                        #notch=TRUE,
                        ylab = ylab,
                        xlab = xlabel,
                        xlim = c(1.3,1.7),
                        width = 0.4,
                        fill = group, 
                        palette = colorCodes,#colorCodes,#c(red1","green4"), #colorCodes,# Codes,
                        #legend="none",
                        #add=c("jitter"),#,"jitter"),#median",
                        #order=OrderCodes, # not needed after reordering factor levels
                        alpha = 0.5, ggtheme = theme_bw())+
      labs(fill="Class")+
      theme(axis.text.x= element_blank(),
            axis.text.y = element_text(size=30),
            axis.title = element_text(size =40),
            legend.title = element_text(size=40),
            legend.text = element_text(size=30),
            legend.key.size = unit(2,"cm"),
            plot.margin=margin(t=0.2,r=0.5,b=0,l=1,unit = "cm"))+
      # geom_signif(comparisons=list(c("1","2"),c("2","3"),c("3","4"),
      #                              c("4","5"),c("5","6")),
      #             map_signif_level = TRUE, 
      #             textsize=6,vjust=2)# Time 0 # margin_top = -0.1=
      geom_signif(comparisons=list(c("Class1","Class2")), 
                  map_signif_level = TRUE, 
                  textsize=14,vjust=2) # Time 7
  }
  # --- Avoid mistakes
  boxCell=box[[1]]
  boxC02=box[[2]]
  boxATP=box[[3]]
  boxX=box[[4]]
  boxG=box[[5]]
  boxN=box[[6]]
  boxP=box[[7]]
  #box[[7]]
  
  fileOut=paste0("MultiPlot_BoxPlots_",mode,"_Time7_matched.pdf")
  pdf(file=fileOut,height = 8,width = 25)
  #if(nplots == 7){
  gg =  ggarrange(boxCell,boxC02,boxATP,boxG,boxX,boxN,boxP,
                   nrow=1,ncol=nbox, legend = c("bottom"),
                   common.legend = TRUE)
  # }else{
  #   gg = ggarrange(boxC02,boxATP,boxG,boxX,boxN,boxP,
  #                  nrow=1,ncol=nplots, legend = c("bottom"),
  #                  common.legend = TRUE)
  # }
  print(gg)
  
  dev.off()
}
