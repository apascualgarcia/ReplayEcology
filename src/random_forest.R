###################################
#     random_forest.R  
##################################
# Author: Alberto Pascual-García
# Copyright (c)  Alberto Pascual-García,  2024
# Web:  apascualgarcia.github.io
# 
# Date: 2024-04-03
# Script Name: random_forest.R   
# Script Description:
#
#
# Notes:
#
#
rm(list = ls())
######### START EDITING ------------
# --- Random forest parameters
optimization=0 # Should RF parameters be optimized? (=1)
run_RF=0 # Should RF be run (=1) or just read from file (=0)
ntree.min=1000 # number of trees, mandatory if optimization = 0 (and only used in that case)
mtry.min="default" # number of variables randomly selected to build the trees. Fix to "default" 
# if you don't have an informed guess
partial.plots=0 # should partial plots be generated (=1), fix to 0 otherwise. It will
# generate plots for the 10 most important variables. This is controlled by
# variable Nsel in the section plots below
select.class="partition" # determine the classification you want to predict, it should
# be present in the metadata table

# SET WORKING DIRECTORY -----------------------------

# --- Input and output files
file.ASV = "seqtable_readyforanalysis.csv"
file.taxa = "taxa_wsp_readyforanalysis.csv"
file.Meta = "metadata_Time0D-7D-4M_May2022_wJSDpart_ext.csv"

# --- Directories
this.dir=strsplit(rstudioapi::getActiveDocumentContext()$path, "/src/")[[1]][1] # don't edit, just comment it if problems...
dirSrc=paste(this.dir,"/src/",sep="") # Directory where the code is
dirASV=paste(this.dir,"/6_finalfiles/",sep="") # Dir of ASV table
dirMeta=paste(this.dir,"/7.3_phyloseq/",sep="") # Dir of metadata
dirOut=paste(this.dir,"/7.6_keystone/random_forest",sep="") # Dir of output data 

# --- Packages

packages <- c("tidyverse", "stringr", "ggplot2",
              "caret", # additional RF functions
              "randomForest", "randomForestExplainer") # list of packages to load

scripts <- c("clean_ASV_table.R") # list of functions to load

# --- Rehaping data parameters
nreads = 10000 # minimum number of reads to consider a sample
exclude_exp = c("4M") # A vector of characters with the experiments that should be excluded
match_exp = TRUE # Set to true if only starting communities that were resurrected should be included
output.label = "Time0D_7D_matched" 

###### STOP EDITING -------------

# INSTALL PACKAGES & LOAD LIBRARIES -----------------
cat("INSTALLING PACKAGES & LOADING LIBRARIES... \n\n", sep = "")

n_packages <- length(packages) # count how many packages are required

new.pkg <- packages[!(packages %in% installed.packages())] # determine which packages aren't installed

# install missing packages
if(length(new.pkg)){
  install.packages(new.pkg)
}

# load all requried libraries
for(n in 1:n_packages){
  cat("Loading Library #", n, " of ", n_packages, "... Currently Loading: ", packages[n], "\n", sep = "")
  lib_load <- paste("library(\"",packages[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}
# SOURCE FUNCTIONS ---------
setwd(dirSrc)
n_scripts <- length(scripts) # count how many packages are required

for(n in 1:n_scripts){
  cat("Loading script #", n, " of ", n_scripts, "... Currently Loading: ", scripts[n], "\n", sep = "")
  lib_load <- paste("source(\"",scripts[n],"\")", sep = "") # create string of text for loading each library
  eval(parse(text = lib_load)) # evaluate the string to load the library
}


# READ INPUT FILES ----------
# --- Read ASVs table
setwd(dirASV)
ASV.table=read.table(file = file.ASV, sep="\t")
colnames(ASV.table)[1:5]
rownames(ASV.table)[1:5]
dim(ASV.table)
head(ASV.table)[1:5,1:5]

# ..... read metadata. Samples present in metadata were those passing the filtering
setwd(dirMeta)
sample_md <-read.table(file = file.Meta, sep="\t", header=TRUE)
head(sample_md)[1:5,1:5]

# CLEAN DATA ---------- 

# Clean data   ------
clean.data.list = clean_ASV_table(ASV.table,sample_md,match.exp = T)

ASV.table = clean.data.list$ASV.table
sample_md = clean.data.list$sample_md

head(ASV.table)[1:5,1:5]

# BUILD REGRESSION DATA ---------- 

# .... Build the table of predictors, first ASVs
# Here we will have a table with parent's ids and ASVs repeated
# four times, and an additional column specifying the replicate (1 to 4)
# we will regress these variables to the response (final class)
id.t0 = which(sample_md$Experiment == "0D")
samples.t0 = as.character(sample_md$sampleid[id.t0])

xIn.tmp = ASV.table[samples.t0, ]
xIn.tmp$sampleid = rownames(xIn.tmp) # create a new column, rownames cannot be repeated
xIn = rbind(xIn.tmp, xIn.tmp, xIn.tmp, xIn.tmp) # same predictors for four different responses
rownames(xIn) = seq(1,dim(xIn)[1]) 
xIn = xIn[,c(dim(xIn)[2],1:(dim(xIn)[2]-1))] # reorder
head(xIn)[1:5,1:5]


# .... Now add the replicate as an additional envir predictor
replicate = c(rep(1, length(id.t0)), rep(2, length(id.t0)),
        rep(3, length(id.t0)),rep(4, length(id.t0)))
xIn$replicate = replicate

# .... Extract response variable
# Each starting condition will have a response, for each parent community
# we will look for the final class

#i = 1 # debug
yIn = c()
for(i in 1:4){
  rep.name = paste0("Rep",i)
  id.rep = which(sample_md$replicate == rep.name)
  matched = match(xIn$sampleid[1:dim(xIn.tmp)[1]], sample_md$parent[id.rep])
  
  yIn.tmp = as.character(sample_md[,select.class][matched])
  
  yIn = c(yIn, yIn.tmp)
}
yIn[which(yIn == "Class1")] = 1
yIn[which(yIn == "Class2")] = 2
yIn = as.factor(yIn)

length(which(yIn == 1))
length(which(yIn == 2))

# ... Finally, drop sampleid names
xIn = subset(xIn, select = -c(sampleid))
head(xIn)[1:5,1:5]

# RANDOM FOREST computation -------------
set.seed(3032024) # today

# --- First estimate the optimal RF parameters:
# the two steps (ntree and mtry) should possibly be iterated
dir.create(dirOut)
setwd(dirOut)
mtry.default=sqrt(dim(xIn)[2]) # default number of variables taken to build the trees
if(optimization == 1){
  #range="large"
  Nrand=25
  if(range == "small"){
    ntree.test=seq(from=100,to=500,by=20) # small range
  }else{
    ntree.test=c(100,250,seq(from=500,to=10000,by=500)) # large
    ntree.test=c(ntree.test,12000,14000,16000,18000,20000) # very large
  }
  mtry.test=mtry.default #300 # if different than default this was fixed after one iteration with the next step below
  #OOB=vector(mode="numeric",length=length(ntree.test))
  OOB=matrix(0,nrow=length(ntree.test),ncol=Nrand)
  i=0
  for(ntree.tmp in ntree.test){
    i=i+1
    for(k in 1:Nrand){
      RF.tmp <-randomForest(y=yIn,x=xIn,
                            importance=T, proximity = T, 
                            ntree=ntree.tmp,mtry = mtry.test)
      OOB[i,k]=mean(RF.tmp$err.rate[,1])
    }
  }
  OOB.mean=rowMeans(OOB)
  OOB.std=apply(OOB,1, sd, na.rm = TRUE)
  OOB.df=data.frame(cbind(ntree.test,OOB.mean,OOB.std))
  #rownames(OOB.df)=ntree.test
  OOB.min.id=which.min(OOB.mean)
  ntree.min=OOB.df$ntree.test[OOB.min.id]
  ymin=OOB.df$OOB.mean-OOB.df$OOB.std/sqrt(Nrand)
  ymax=OOB.df$OOB.mean+OOB.df$OOB.std/sqrt(Nrand)
  fileOut=paste("optimization_ntree_Class-",select.class,"_mtry",trunc(mtry.test),
                "_",range,".csv",sep="")
  write.table(OOB.df,file = fileOut,sep="\t",quote = FALSE,row.names = FALSE)
  plotOut=paste("Plot_optimization_ntree_Class-",select.class,"_mtry",trunc(mtry.test),
                "_",range,".pdf",sep="")
  pdf(plotOut)
  g=ggplot()+
    geom_vline(xintercept = ntree.min,linetype = 'dotted', col = 'red')+
    geom_point(data=OOB.df,aes(x=ntree.test,y=OOB.mean))+
    ylab("Mean OOB error")+xlab("Number of trees")+
    scale_y_continuous(trans='log10')+
    geom_errorbar(aes(x=ntree.test,ymin=ymin,ymax=ymax))+
    theme_bw()
  print(g)
  dev.off()
  
  # --- Now we fix the optimal ntree and look for the optimization of mtry
  ntree.in=ntree.min
  #ntree.in = 6000 # same order of magnitude than the minimum.
  mtry.test=seq(from=mtry.default/2,to= dim(xIn)[2],by=mtry.default/2)
  #OOB=vector(mode="numeric",length=length(mtry.test))
  OOB=matrix(0,nrow=length(mtry.test),ncol=Nrand)
  i=0
  for(mtry.tmp in mtry.test){
    i=i+1
    for(k in 1:Nrand){
      RF.tmp <-randomForest(y=yIn,x=xIn,
                            importance=T, proximity = T, 
                            ntree=ntree.in,mtry=mtry.tmp)
      OOB[i,k]=mean(RF.tmp$err.rate[,1])
    }
  }
  OOB.mean=rowMeans(OOB)
  OOB.std=apply(OOB,1, sd, na.rm = TRUE)
  OOB.df=data.frame(cbind(mtry.test,OOB.mean,OOB.std))
  #rownames(OOB.df)=ntree.test
  OOB.min.id=which.min(OOB.mean)
  mtry.min=OOB.df$mtry.test[OOB.min.id]
  ymin=OOB.df$OOB.mean-OOB.df$OOB.std/sqrt(Nrand)
  ymax=OOB.df$OOB.mean+OOB.df$OOB.std/sqrt(Nrand)
  fileOut=paste("optimization_mtry_Class-",select.class,"_ntree-",trunc(ntree.min),
                "_",range,".csv",sep="")
  write.table(OOB.df,file = fileOut,sep="\t",quote = FALSE,row.names = FALSE)
  plotOut=paste("Plot_optimization_mtry_Class-",select.class,
                "_ntree",ntree.in,".pdf",sep="")
  pdf(plotOut)
  g=ggplot()+
    geom_point(data=OOB.df,aes(x=mtry.test,y=OOB.mean))+
    ylab("Mean OOB error")+xlab("Number of variables")+
    geom_errorbar(aes(x=mtry.test,ymin=ymin,ymax=ymax))
  print(g)
  dev.off()
}else{ # if we do not optimize
  if(mtry.min == "default"){ # we need to give a value if the user choose a default value
    mtry.min = mtry.default
  }
}

# --- With the optimal parameters run the RF
ntree.in=ntree.min
mtry.in=mtry.min
fileOut=paste("RandForestOut_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".RDS",sep="")
if(run_RF == 1){
  RF.out <-randomForest(y=yIn,x=xIn,
                        importance=T, proximity = T, 
                        ntree=ntree.in,mtry=mtry.in)
  # --- Have a look at the output
  RF.out
  mean(RF.out$err.rate[,1]) # this is the mean OOB
  saveRDS(RF.out, file = fileOut)
}else{
  RF.out = readRDS(file = fileOut)
}

RF.out 

# ANALYSE -------------

# --- Extract important variables
setwd(dirOut)
importance.df = measure_importance(RF.out)
ASV.top = important_variables(importance.df)

fileOut=paste("varImpExt_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".tsv",sep="")

write.table(importance.df,file=fileOut,sep = "\t",quote=FALSE)

fileOut=paste("varTop_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".tsv",sep="")

write.table(ASV.top,file=fileOut,sep = "\t",quote=FALSE)

# Plots ---------

# first and overview, it will generate an html file. Takes time.
explain_forest(RF.out)

# --- Error
plotOut=paste("Plot_Error_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".pdf",sep="")
pdf(plotOut)
plot(RF.out)
dev.off()


# --- Variable importance
fileOut=paste("varImp_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".tsv",sep="")
varImp.out=varImp(RF.out)
varImp.out=varImp.out[which(rowSums(varImp.out) != 0),]
write.table(varImp.out,file=fileOut,sep = "\t",quote=FALSE)
plotOut=paste("Plot_varImp_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".pdf",sep="")
pdf(plotOut,width=10)
#varImpPlot(RF.out, type=1,main="")#,xlab="Variable importance",title="")
#dev.off()

# Get variable importance from the model fit
ImpData <- as.data.frame(importance(RF.out))
ImpData$Var.Names <- row.names(ImpData)
ImpData.sort.idx = sort(ImpData$MeanDecreaseAccuracy, 
                        decreasing = TRUE, index.return = T)
quantile(ImpData.sort.idx$x)
Imp.Data.sort = ImpData[ImpData.sort.idx$ix, ]
explained = 50
id.exp = which(Imp.Data.sort$MeanDecreaseAccuracy > explained)
Imp.Data.sort = Imp.Data.sort[id.exp,]

ggplot(Imp.Data.sort, aes(x=Var.Names, y=MeanDecreaseAccuracy)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=MeanDecreaseAccuracy), color="skyblue") +
  geom_point(aes(size = MeanDecreaseGini), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

dev.off()

# Variable importance, multiway

plotOut=paste("Plot_accuracyVsGini_Class-",select.class,
              "_ntree",trunc(ntree.in),"_mtry",trunc(mtry.in),
              ".pdf",sep="")

pdf(plotOut, height =6)
p = plot_multi_way_importance(
  importance_frame = importance.df,
  x_measure = "accuracy_decrease",
  y_measure = "gini_decrease",
  size_measure = "p_value",
  min_no_of_trees = 0,
  no_of_labels = 10 #,
  #main = "Multi-way importance plot"
)
p = p + xlab("Accuracy decrease") +ylab("Gini index decrease")
p
#print(p)
dev.off()

# --- Partial plots sorted by importance
# ... these plots may take a long time
if(partial.plots == 1){
  Nsel=10 # select only 10 vars
  imp <- importance(RF.out)
  impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
  
  impvar=impvar[1:Nsel]
  Nlev=levels(yIn)
  
  # Interpret partial plots
  # https://stats.stackexchange.com/questions/147763/meaning-of-y-axis-in-random-forest-partial-dependence-plot
  # p=seq(from=0.005, to=0.995, by=0.005) # to understand the plots
  # plot(p,log(p/(1-p))) # plot the logit function
  for (i in seq_along(impvar)) {
    var.lab=str_replace_all(impvar[i],pattern="[[:punct:]]",replacement = "")
    var.lab=str_replace_all(var.lab,pattern = " ",replacement = "_")
    plotOut=paste("Plot_PartialDependence_",trunc(ntree.in),"_mtry",trunc(mtry.in),
                  "Class-",select.class,"_Var-",var.lab,
                  ".pdf",sep="")
    pdf(file=plotOut,width=12)
    op <- par(mfrow=c(2, 3))
    for(level in levels(yIn)){
      partialPlot(RF.out, xIn, impvar[i], xlab=impvar[i],
                  which.class = level,
                  main=paste("Partial Dependence for class", level)) #ylim=c(30, 70))
    }
    par(op)
    dev.off()
  }
}


# 
# 
# data(airquality)
# airquality <- na.omit(airquality)
# set.seed(131)
# ozone.rf <- randomForest(Ozone ~ ., airquality, importance=TRUE)
# imp <- importance(ozone.rf)
# impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
# op <- par(mfrow=c(2, 3))
# for (i in seq_along(impvar)) {
#   partialPlot(ozone.rf, airquality, impvar[i], xlab=impvar[i],
#               main=paste("Partial Dependence on", impvar[i]),
#               ylim=c(30, 70))
# }
# par(op)
