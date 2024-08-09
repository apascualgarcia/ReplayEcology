##################################################
# prepare_NCBIupload_filtered.R
##################################################

# The 2843 filtered sequence files were sorted into separate folders 
# based on study and (for this study) replicate within study, in order 
# to facilitate their deposition at NCBI (given only 1000 files can be 
# deposited in each submission) and potential users' future interaction 
# with them. For reproducibility purposes, it was necessary to deposit all 
# 2843 files (even those not specific to this study) because ASV inference 
# was performed on all samples (see below). This process required the use of 
# text processing in R, given there was no explicit label for study in the filenames.

# Cardiff, July-August 2024
# European Centre for the Environment & Human Health, University of Exeter
# Matt Lloyd Jones
# https://github.com/befriendabacterium/


here::here()
library(dplyr)

# READ IN DATA ------------------------------------------------------------

#read in sample metadata
sample_md<-read.csv('4_dada2/metadata_Time0D-7D-4M_May2022.csv', sep='')

#read in geographical attributes (additional metadata)
geog_attrib<-read.csv('1_raw/geographical_attributes.csv')
geog_attrib$latlong<-paste(geog_attrib$Northings,'N',geog_attrib$Westings,'W')

#get filenames of the 2843 sequence files
filename_fullpaths<-list.files(c('3_filtered'), full.names = T, include.dirs = F, pattern = '.fastq.gz')
filename<-gsub('3_filtered/','',filename_fullpaths)

# PARSE THE DATA IN THE FILE NAMES INTO A CONVENTION ----------------------

#create vector of the various ways in which columns are labelled
columns<-c(1,2,3,4,5,6,7,8,9,'01','02','03','04','05','06','07','08','09','10','11','12')

# #create empty vector of wells
wells<-c()
#loop through column labels appending to wells to get all possible well labels
for (i in 1:length(columns)){
  wells<-c(wells,paste(LETTERS[1:8], columns[i], sep=''))
}

# # #add underscores either side so know to remove this string
wells_strings<-paste(wells,'.',sep='')
#join them with an 'or' (|) operator
wells_strings<-paste(wells_strings, collapse='|')

#generate a string for the run by matching run names
run<-unlist(stringr::str_extract_all(filename, ("052214DR16s|021216DR515F")))

#generate a string for the day by matching day names in filenames
day<-unlist(stringr::str_extract_all(filename, ("day0|day7")))

#create sample name column by removing file extension
sample_name<-gsub(".fastq.gz","",filename)

#snip out the well names (everything before first period with sub)
sample_name<-sub(".*?\\.", "", sample_name)
sample_name

#remove days from sample names
sample_name<-gsub("day0_|day7_","",sample_name)

#sample title (e.g. community name)
sample_name<-sapply(strsplit(sample_name, "_"), `[`, 1)

#make site column
site<-sub("^([[:alpha:]]*).*", "\\1", sample_name)
#make sample number per site column (N.B. some will be NA because don't have such numbers)
samplenumber_site<-as.numeric(gsub(".*?([0-9]+).*", "\\1",  sample_name))
#create latlong column in Northings and Westings
latlong<-geog_attrib$latlong[stringdist::amatch(gsub("\\..*","", sample_name), geog_attrib$SampleID, maxDist=2)] 
#create date of collection column
collection_date<-geog_attrib$Date[stringdist::amatch(gsub("\\..*","", sample_name), geog_attrib$SampleID, maxDist=2)] # add collection date column
#convert date to sensible format
collection_date<-as.character(lubridate::ymd(collection_date))

#combine into one dataframe
NCBIdf_all<-as.data.frame(
  cbind(run,
        day,
        filename,
        filename_fullpaths,
        sample_name,
        site,
        samplenumber_site,
        latlong,
        collection_date
  )
)

#remove duplicate communities (exactly same sequences so must have been duplicated by MRDNA post sequencing when batching up files since numbers of samples in batches overlap sometimes)
#identify duplicate rows based on the community and rep columns using the 'duplciated' function
#we define the lower numbered ones as duplicateds (fromLast=T) 
duprows<-which(duplicated(NCBIdf_all$sample_name, fromLast = T))

#separate duplicates into new df
NCBIdf_all_dups<-NCBIdf_all[duprows,]
NCBIdf_all_dups$inclusionexclusionreason<-'duplicate'

#remove duplicates from main df
NCBIdf_all<-NCBIdf_all[-duprows,]

#separate core communities into new df
NCBIdf_all_core<-NCBIdf_all[grep('[0-9]a|c', NCBIdf_all$sample_name, invert=F),]
NCBIdf_all_core$inclusionexclusionreason<-'core community'

#remove 96 core communities from day 0
NCBIdf_all<-NCBIdf_all[grep('[0-9]a|c', NCBIdf_all$sample_name, invert=T),]

# FILTER OUT RELEVANT SAMPLES FROM SEQUENCE TABLE ----------------------------------------------------------

#create vectors of sample IDs per experiment separated with | (or) operator for searching
stringstofind_thisstudy<-paste(sample_md$parent, collapse = '|')
stringstofind_thisstudy<-paste(stringstofind_thisstudy,'Negative',sep='|')
stringstofind_scheuerlstudy<-'Iso|Mix|Nov2|Riz|CTRL|Bag'
stringstofind_mombrikotbstudy<-c('KF|MU|PUD|PU|RU|LF|MF|NF|QF|PF|RF|SF|KU|QU|NU|SU|LU|QU|NEG|POS')
stringstonotfind_misc<-paste(stringstofind_thisstudy,stringstofind_scheuerlstudy,stringstofind_mombrikotbstudy, sep='|')

#identify the samples/rows needed using the above vector
samples_thisstudy<-grep(stringstofind_thisstudy, NCBIdf_all$sample_name)
samples_scheuerlstudy<-grep(stringstofind_scheuerlstudy, NCBIdf_all$sample_name)

#remove the scheuerl studies from the treehole dataset
samples_thisstudy<-samples_thisstudy[-which(samples_thisstudy%in%samples_scheuerlstudy)]
samples_mombrikotbstudy<-grep(stringstofind_mombrikotbstudy,NCBIdf_all$sample_name)
samples_misc<-grep(stringstonotfind_misc, NCBIdf_all$sample_name, invert = T)

# ASSEMBLE NCBI DATAFRAMES FOR UPLOADING ----------------------------------
#Here we assemble dataframes of just the treehole samples from this experiment, for uploading to NCBI. 
#So we copy the relevant data to '3_filtered' (a few steps back in the pipeline), to facilitate re-use and uploading to NCBI

#create NCBI df for this study
NCBIdf_thisstudy<-NCBIdf_all[samples_thisstudy,]
NCBIdf_thisstudy$inclusionexclusionreason<-'this study'

length(grep("day0",NCBIdf_thisstudy$day))

#check number of samples for this study
nrow(NCBIdf_thisstudy)

#add additional day column fr this experiment, based on sequencing run
#NCBIdf_thisstudy$day<-as.factor(NCBIdf_thisstudy$run)
#levels(NCBIdf_thisstudy$day)<-c(7,0)

#the day 0 negative appears to have been run on the second run instead of the first, so to keep things balanced (i.e. 1 negative per rep), relabel as day 0
NCBIdf_thisstudy$day[NCBIdf_thisstudy$sample_name=='Negative']<-'day0'

#split this study by sampling occasion/replicate because only allowed 1000 files per deposit
NCBIdf_thisstudy_day0<-NCBIdf_thisstudy[NCBIdf_thisstudy$day=='day0',]
NCBIdf_thisstudy_day7r1<-NCBIdf_thisstudy[grep('\\.1',NCBIdf_thisstudy$sample_name),]
NCBIdf_thisstudy_day7r2<-NCBIdf_thisstudy[grep('\\.2',NCBIdf_thisstudy$sample_name),]
NCBIdf_thisstudy_day7r3<-NCBIdf_thisstudy[grep('\\.3',NCBIdf_thisstudy$sample_name),]
NCBIdf_thisstudy_day7r4<-NCBIdf_thisstudy[grep('\\.4',NCBIdf_thisstudy$sample_name),]

#we are going to bind the rogue Negative that doesn't have a rep number to rep 1, just so everything is accounted for...
#NCBIdf_thisstudy_day7r1<-rbind(NCBIdf_thisstudy_day7r1,NCBIdf_thisstudy[NCBIdf_thisstudy$sample_name=='Negative',])

#check all add up to 2180
sum(nrow(NCBIdf_thisstudy_day0),
    nrow(NCBIdf_thisstudy_day7r1),
    nrow(NCBIdf_thisstudy_day7r2),
    nrow(NCBIdf_thisstudy_day7r3),
    nrow(NCBIdf_thisstudy_day7r4)
)

#create NCBI df for other studies and samples
NCBIdf_scheuerlstudy<-NCBIdf_all[samples_scheuerlstudy,]
NCBIdf_scheuerlstudy$inclusionexclusionreason<-'scheuerl study'

NCBIdf_mombrikotbstudy<-NCBIdf_all[samples_mombrikotbstudy,]
NCBIdf_mombrikotbstudy$inclusionexclusionreason<-'mombrikotb study'

NCBIdf_misc<-NCBIdf_all[samples_misc,]
NCBIdf_misc$inclusionexclusionreason<-'does not belong to any of the three studies identified'

NCBIdf_misc<-rbind(NCBIdf_misc,NCBIdf_all_dups) #put duplicates back into misc list
NCBIdf_misc<-rbind(NCBIdf_misc,NCBIdf_all_core) #put duplicates back into misc list

#check all add up to 2843
sum(nrow(NCBIdf_thisstudy_day0),
    nrow(NCBIdf_thisstudy_day7r1),
    nrow(NCBIdf_thisstudy_day7r2),
    nrow(NCBIdf_thisstudy_day7r3),
    nrow(NCBIdf_thisstudy_day7r4),
    nrow(NCBIdf_scheuerlstudy),
    nrow(NCBIdf_mombrikotbstudy),
    nrow(NCBIdf_misc)
)

#remove irrelevant columns from last 3 dataframes
irrelevant_columns<-c('site','day','samplenumber_site','latlong','collection_date')
NCBIdf_scheuerlstudy<-dplyr::select(NCBIdf_scheuerlstudy,-irrelevant_columns)
NCBIdf_mombrikotbstudy<-dplyr::select(NCBIdf_mombrikotbstudy,-irrelevant_columns)
NCBIdf_misc<-dplyr::select(NCBIdf_misc,-irrelevant_columns)

#clean up workspace
rm(list=setdiff(ls(), c("NCBIdf_thisstudy_day0",
                        "NCBIdf_thisstudy_day7r1",
                        "NCBIdf_thisstudy_day7r2",
                        "NCBIdf_thisstudy_day7r3",
                        "NCBIdf_thisstudy_day7r4",
                        "NCBIdf_scheuerlstudy",
                        "NCBIdf_mombrikotbstudy",
                        "NCBIdf_misc")))

#write csv for each timepoint/rep combo
write.csv(NCBIdf_thisstudy_day0,'3_filtered/NCBIdf_thisstudy_day0.csv', row.names=F)
write.csv(NCBIdf_thisstudy_day7r1,'3_filtered/NCBIdf_thisstudy_day7r1.csv', row.names=F)
write.csv(NCBIdf_thisstudy_day7r2,'3_filtered/NCBIdf_thisstudy_day7r2.csv', row.names=F)
write.csv(NCBIdf_thisstudy_day7r3,'3_filtered/NCBIdf_thisstudy_day7r3.csv', row.names=F)
write.csv(NCBIdf_thisstudy_day7r4,'3_filtered/NCBIdf_thisstudy_day7r4.csv', row.names=F)
write.csv(NCBIdf_scheuerlstudy,'3_filtered/NCBIdf_scheuerlstudy.csv', row.names=F)
write.csv(NCBIdf_mombrikotbstudy,'3_filtered/NCBIdf_mombrikotbstudy.csv', row.names=F)
write.csv(NCBIdf_misc,'3_filtered/NCBIdf_misc.csv', row.names=F)

# COPY FILES --------------------------------------------------------------

#copy sequence files to separate folders for uploading to NCBI
file.copy(NCBIdf_thisstudy_day0$filename_fullpaths,'3_filtered/thisstudy_day0/')
file.copy(NCBIdf_thisstudy_day7r1$filename_fullpaths,'3_filtered/thisstudy_day7r1')
file.copy(NCBIdf_thisstudy_day7r2$filename_fullpaths,'3_filtered/thisstudy_day7r2')
file.copy(NCBIdf_thisstudy_day7r3$filename_fullpaths,'3_filtered/thisstudy_day7r3')
file.copy(rev(NCBIdf_thisstudy_day7r4$filename_fullpaths),'3_filtered/thisstudy_day7r4')
file.copy(NCBIdf_scheuerlstudy$filename_fullpaths,'3_filtered/scheuerlstudy',)
file.copy(NCBIdf_mombrikotbstudy$filename_fullpaths,'3_filtered/mombrikotbstudy')
file.copy(NCBIdf_misc$filename_fullpaths,'3_filtered/misc')
