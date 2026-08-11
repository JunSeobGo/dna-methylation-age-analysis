# Copy and paste the following R software code
# Use forward slashes /”, as R will misread filepaths with backslashes#
setwd("C:/Users/ko911/OneDrive/바탕 화면/UnsortedStudy/DNAage") # 경로 설정

#install.packages("BiocManager")
#BiocManager::install("WGCNA")
#install.packages('sqldf')

library(WGCNA)
library(sqldf)

## Install the Bioconductor installer
#install.packages("BiocManager")

## Install the "impute" package from Bioconductor
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#BiocManager::install("impute")

## Install the devtools package if not already installed
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}

#install.packages("RPMM")
#library(RPMM)
source("C:/Users/ko911/OneDrive/바탕 화면/UnsortedStudy/DNAage/AdditionalFile24NROMALIZATION.R.txt")


#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("AdditionalFile22probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("AdditionalFile21datMiniAnnotation27k.csv")
datClock=read.csv("AdditionalFile23predictor.csv")

#install.packages("sqldf")
library(sqldf)

#Read in the DNA methylation data (beta values)
# For a small file, e.g. measured on the 27k platform you could just use read.csv. 
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.
dat0=read.csv.sql("AdditionalFile26MethylationDataExample55.csv") ;
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")
#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if (   sum(selectXchromosome) >=500 )  {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  } 
}



match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) { 
  missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="LogFile.txt",append=TRUE)  } 


#STEP 2: Restrict the data to 21k probes and ensure they are numeric
match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
dat1= dat0[match1,]
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)


# 단계 3: 결과 파일 datout 생성
set.seed(1)
# 데이터 정규화 수행 여부 (권장)
normalizeData = TRUE
source("AdditionalFile25StepwiseAnalysis.txt")


# 단계 4: 결과 출력
if (sum(datout$Comment != "") == 0) { 
  cat("\n 개별 샘플은 정상적으로 처리되었습니다.", file = "LogFile.txt", append = TRUE) 
} 
if (sum(datout$Comment != "") > 0) { 
  cat(paste("\n 경고: 다음 샘플에 대해 경고가 생성되었습니다.\n", datout[, 1][datout$Comment != ""], "\n 자세한 내용은 로그 파일을 확인하세요."), file = "LogFile.txt", append = TRUE) 
}

# 결과를 디렉토리에 출력
write.table(datout, "Output.csv", row.names = FALSE, sep = ",")

######test###############

# 이 작업을 수행하기 위해 연령 데이터가 포함된 샘플 주석 데이터를 읽습니다.

######Here!!!! We can change the test data!
# 여기에서 테스트 데이터를 변경할 수 있습니다.
datSample = read.csv("AdditionalFile27SampleAnnotationExample55.csv")


DNAmAge = datout$DNAmAge
medianAbsDev = function(x, y) median(abs(x - y), na.rm = TRUE)
medianAbsDev1 = signif(medianAbsDev(DNAmAge, datSample$Age), 2)
par(mfrow = c(1, 1))
verboseScatterplot(DNAmAge, datSample$Age, xlab = "DNAm Age", ylab = "Chronological Age", main = paste("All, err=", medianAbsDev1))
abline(0, 1) 