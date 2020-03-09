
mapValues=function(x,keyForMapping) {
  x=plyr::mapvalues(x,from=names(keyForMapping),to=keyForMapping,warn_missing=FALSE)
  return(x)
}

summaryTable<-function(dataForReport) {
  reportTable<-data.frame()
  if (!is.null(dataForReport$DuplicateSnpConcordance)) {
    temp<-NROW(dataForReport$DuplicateSnpConcordance$SampleLevel)
    reportTable<-rbind(reportTable,c("Samples have low concordance in duplicate SNPs",temp),stringsAsFactors =FALSE)
  }
  if (!is.null(dataForReport$DuplicateSamplesConcordance)) {
    temp<-NROW(dataForReport$DuplicateSamplesConcordance$SnpLevel)
    reportTable<-rbind(reportTable,c("Snps have low concordance in duplicate Samples",temp),stringsAsFactors =FALSE)
  }
  if (!is.null(dataForReport$G1000SamplesConcordance)) {
    temp<-NROW(dataForReport$G1000SamplesConcordance$SnpLevel)
    reportTable<-rbind(reportTable,c("Snps have low concordance in replicate 1000 Genomes Samples",temp),stringsAsFactors =FALSE)
    temp<-NROW(dataForReport$G1000SamplesConcordance$SampleLevel)
    reportTable<-rbind(reportTable,c("Samples have low concordance in replicate 1000 Genomes Samples",temp),stringsAsFactors =FALSE)
  }
  if (!is.null(dataForReport$GenderCheckTable)) {
    temp<-NROW(dataForReport$GenderCheckTable)
    reportTable<-rbind(reportTable,c("Samples have different gender in gender check",temp),stringsAsFactors =FALSE)
  }
  reportTable[,2]<-as.integer(reportTable[,2])
  colnames(reportTable)<-c("Category","Count")
  return(reportTable)
}


aimSnpMatrixPca<-function(snpMatrix) {
  #using this method: http://bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/pca-vignette.pdf
  #Alternative solution: http://www.popgen.dk/software/index.php/Rscripts
  snpMatrixXxt <- xxt(snpMatrix, correct.for.missing=FALSE)
  evv <- eigen(snpMatrixXxt,symmetric=TRUE)
  pcs <- evv$vectors[,1:10]
  return(pcs)
}

snpMatrixToRatio<-function(snpMatrix) {
  ratio01By11<-apply(as(snpMatrix, 'numeric'),1,function(x) length(which(x==1))/length(which(x==2)))
  result<-data.frame(ratio01By11=ratio01By11)
  row.names(result)<-row.names(snpMatrix)
  return(result)
}


#temp=snpsComp(snpMatrix1,snpMatrix2,dupSampleFile=sampleToG1000File)
snpsComp<-function(snpMatrix1,snpMatrix2=NULL,dupSnpFile=NULL,dupSampleFile=NULL,snpConcordanceCutoff=NULL,sampleConcordanceCutoff=NULL) {
  if (is.null(snpMatrix2)) {
    snpMatrix2=snpMatrix1
  } else { #not the same snp, keep same probes only
    temp<-intersect(colnames(snpMatrix1),colnames(snpMatrix2))
    snpMatrix1<-snpMatrix1[,temp]
    snpMatrix2<-snpMatrix2[,temp]
  }

  if (!is.null(dupSnpFile)) { #compare duplicate snps in all samples
    dupSnps<-read.delim(dupSnpFile,header=F,as.is=T)
    #only keep snp pairs both in the data
    dupSnps[,1]<-gsub(" +","",dupSnps[,1])
    dupSnps[,2]<-gsub(" +","",dupSnps[,2])
    temp1<-dupSnps[,1] %in% colnames(snpMatrix1)
    temp2<-dupSnps[,2] %in% colnames(snpMatrix2)
    dupSnpInd<-which((as.integer(temp1)+as.integer(temp2))==2)
    print(paste0(length(dupSnpInd)," duplicate probe pairs identified in the data"))
    dupSnps<-dupSnps[dupSnpInd,]

    temp1<-snpMatrix1[,as.character(dupSnps[,1])]
    temp2<-snpMatrix2[,as.character(dupSnps[,2])]
    colnames(temp1)<-make.unique(colnames(temp1))
    colnames(temp2)<-colnames(temp1)
    snpsCompResult<-sm.compare(temp1,temp2)
  } else if (!is.null(dupSampleFile)) {
    dupSamples<-read.delim(dupSampleFile,header=F,as.is=T)
    #only keep smaple pairs both in the data
    dupSamples[,1]<-gsub(" +","",dupSamples[,1])
    dupSamples[,2]<-gsub(" +","",dupSamples[,2])
    temp1<-dupSamples[,1] %in% row.names(snpMatrix1)
    temp2<-dupSamples[,2] %in% row.names(snpMatrix2)
    dupSampleInd<-which((as.integer(temp1)+as.integer(temp2))==2)
    print(paste0(length(dupSampleInd)," duplicate sample pairs identified in the data"))
    dupSamples<-dupSamples[dupSampleInd,]

    temp1<-snpMatrix1[as.character(dupSamples[,1]),]
    temp2<-snpMatrix2[as.character(dupSamples[,2]),]
    row.names(temp1)<-make.unique(row.names(temp1))
    row.names(temp2)<-row.names(temp1)
    snpsCompResult<-sm.compare(temp1,temp2)
  } else {
    break(paste0("One of dupSnpFile or dupSampleFile should be provided."))
  }

  totalSampleNum=nrow(snpsCompResult$row.wise)
  totalSnpNum=nrow(snpsCompResult$col.wise)

  #analysis SNP level
  callSampleNum<-totalSampleNum-rowSums(snpsCompResult$col.wise[,c("NA.agree","NA.disagree")])
  homAgreeRate<-snpsCompResult$col.wise[,"Hom.agree"]/callSampleNum
  allAgreeRate<-rowSums(snpsCompResult$col.wise[,c("Hom.agree","Het.agree")])/callSampleNum
  if (!is.null(dupSnpFile)) {
	  snpsCompResultSnpLevel<-data.frame(SnpDup1=dupSnps[,1],SnpDup2=dupSnps[,2],SamplesCalled=callSampleNum,
			  HomConcordanceRate=homAgreeRate,AllConcordanceRate=allAgreeRate,
			  stringsAsFactors=FALSE)
  } else if (!is.null(dupSampleFile)) {
	  snpsCompResultSnpLevel<-data.frame(Snp=row.names(snpsCompResult$col.wise),SamplesCalled=callSampleNum,
			  HomConcordanceRate=homAgreeRate,AllConcordanceRate=allAgreeRate,
			  stringsAsFactors=FALSE)
  } else {
	  break(paste0("One of dupSnpFile or dupSampleFile should be provided."))
  }
  if (!is.null(snpConcordanceCutoff)) {
	  diffInd<-which(snpsCompResultSnpLevel$HomConcordanceRate<=snpConcordanceCutoff | snpsCompResultSnpLevel$AllConcordanceRate<=snpConcordanceCutoff)
	  snpsCompResultSnpLevel<-snpsCompResultSnpLevel[diffInd,]
	  print(paste0(length(diffInd)," Snps have less than ",snpConcordanceCutoff, " Concordance"))
  }
  if (nrow(snpsCompResultSnpLevel)!=0) {
	  snpsCompResultSnpLevel<-snpsCompResultSnpLevel[order(snpsCompResultSnpLevel$AllConcordanceRate),]
  } else {
	  snpsCompResultSnpLevel<-NULL
  }

  #analysis sample level
  callSnpNum<-totalSnpNum-rowSums(snpsCompResult$row.wise[,c("NA.agree","NA.disagree")])
  homSampleAgreeRate<-snpsCompResult$row.wise[,"Hom.agree"]/callSnpNum
  allSampleAgreeRate<-rowSums(snpsCompResult$row.wise[,c("Hom.agree","Het.agree")])/callSnpNum
  if (!is.null(dupSnpFile)) {
	  snpsCompResultSampleLevel<-data.frame(Sample=row.names(snpsCompResult$row.wise),SnpsCalled=callSnpNum,
			  HomConcordanceRate=homSampleAgreeRate,AllConcordanceRate=allSampleAgreeRate,
			  stringsAsFactors=FALSE)
  } else if (!is.null(dupSampleFile)) {
	  snpsCompResultSampleLevel<-data.frame(Sample1=dupSamples[,1],Sample2=dupSamples[,2],SnpsCalled=callSnpNum,
			  HomConcordanceRate=homSampleAgreeRate,AllConcordanceRate=allSampleAgreeRate,
			  stringsAsFactors=FALSE)
  } else {
	  break(paste0("One of dupSnpFile or dupSampleFile should be provided."))
  }
  if (!is.null(sampleConcordanceCutoff)) {
	  diffInd<-which(snpsCompResultSampleLevel$HomConcordanceRate<=sampleConcordanceCutoff | snpsCompResultSampleLevel$AllConcordanceRate<=sampleConcordanceCutoff)
	  snpsCompResultSampleLevel<-snpsCompResultSampleLevel[diffInd,]
	  print(paste0(length(diffInd)," Samples have less than ",sampleConcordanceCutoff, " Concordance"))
  }
  if (nrow(snpsCompResultSampleLevel)!=0) {
	  snpsCompResultSampleLevel<-snpsCompResultSampleLevel[order(snpsCompResultSampleLevel$AllConcordanceRate),]
  } else {
	  snpsCompResultSampleLevel<-NULL
  }

  return(list(SnpLevel=snpsCompResultSnpLevel,SampleLevel=snpsCompResultSampleLevel))
}
