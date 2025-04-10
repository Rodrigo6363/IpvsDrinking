####Open new file and copy and paste this script into it to edit it, rest of scripts remain uneditted unless you need to uncomment a line

######### I. Define General ########
#Set these variables

#Set working directory, folder names, and file names for data and metadata
setwd("")
FolderName<-"src/data/" #folder name with / at the end
FolderList <-"src/data/protein_lists/" 
Folderfig <-"src/figures"
diann_data<-"src/data/raw_data/" 
FASTA<-"src/data/metadata_fasta"
Metadata<-"src/data/metadata_fasta/samplemetadata.xlsx"

#uncomment if using spectronaut

#Set inputs: names and group
NameCond1<- "IP"#generally minus; comes first in table; second in contrast (ex. + in +vs- contrast)
NameCond2<- "DR" #generally plus; comes second in table; first in contrast and in first contrast with NC
NameControl<- "Ctrl" #generally control; comes last in table

#Set numbers of variables
cond1<-2 #IP; comes first in table
cond2<-2 #DR; comes second in table
control<-1 #control; comes last in table

#Set number of replicates to look in
#ex. look in a out of cond1 samples (5 out of 6)
a<-cond1-1
b<-cond2-1
c<-control-1 

#Set cutoffs
pval<-0.05
qval<-0.05
foldlog2<-1.5

#run MS-DAP1 Pipeline script and then return to II

#### II. MSDAP 1 variables#####
#Set these variables
#MS-DAP1 downstream Analysis: uncomment line 37 if need to change where controls are to the end
#Don't forget to open libraries at the top of MS-DAP1 downstream Analysis

#Set time stamp of MS-DAP1 output
DateTimeStamp<-"2025-04-08_17-10-05" #fill in after running MSDAP1
#Set Column assignments from full prot.input 

col.start<-4
col.end<-8

##Don't change below; these are the dynamic labels
foldchange.colNam1<-paste0("foldchange.log2.",NameCond2,".",NameControl)
pvalue.colNam1<-paste0("pvalue.log2.",NameCond2,".",NameControl)
qvalue.colNam1<-paste0("qvalue.log2.",NameCond2,".",NameControl)
foldchange.colNam2<-paste0("foldchange.log2.",NameCond1,".",NameControl)
pvalue.colNam2<-paste0("pvalue.log2.",NameCond1,".",NameControl)
qvalue.colNam2<-paste0("qvalue.log2.",NameCond1,".",NameControl)

#Column assignments after selecting columns and making labels automated
total<-cond1+cond2+control
prot.col.end<-total+3
col.start.cond1<-4
col.end.cond1<-cond1+3
col.start.cond2<-cond1+3+1
col.end.cond2<-cond1+cond2+3
col.start.ctrl<-cond1+cond2+3+1


#foldchange and p-value/q-value columns names for dea file
foldchange.colNam<-paste0("foldchange.log2.",NameCond1,".",NameCond2)
foldchange2.colNam<-paste0("foldchange.log2.",NameCond1,".",NameCond2)
pvalue.colNam<-paste0("pvalue.log2.",NameCond2,".",NameCond1)
qvalue.colNam<-paste0("qvalue.log2.",NameCond2,".",NameCond1)



