#' @export
initialize_report_parameter<-
  function(
    plinkBed,
    ReportRace=FALSE,
    RaceToG1000Maf=c(A="EAS",ASIAN="EAS",B="AFR",BLACK="AFR",W="EUR",WHITE="EUR",H="AMR",I="SAS"),
    MafOutlierCut=0.3
  ) {

    dataForReport=list()
    dataForReport$ReportRace=ReportRace
    dataForReport$RaceToG1000Maf=RaceToG1000Maf
    dataForReport$MafOutlierCut=MafOutlierCut

    return(dataForReport)
}

#' @export

preparePlinkData<-function(plinkBed,checkRefAllele=TRUE,dupSnpFile=NULL,dupSampleFile=NULL,aimFile=NULL,...) {
  dataForReport<-list()
  #####################################
  #read plink data
  #####################################

  plinkRaw <- read.plink(plinkBed,...) #bim and fam will be added automaticlly
  #row.names(plinkRaw$genotypes)<-plinkRaw$fam$member

  totalSampleNum=nrow(plinkRaw$fam)
  totalSnpNum=nrow(plinkRaw$map)
  print(paste0("Plink file with ",totalSampleNum," samples and ",totalSnpNum," SNPs was loaded."))

  #####################################
  #probe summary and sample summary
  #####################################
  dataForReport$SnpSummary <- col.summary(plinkRaw$genotypes)
  dataForReport$SnpSummary<-cbind(plinkRaw$map,dataForReport$SnpSummary)
  setDT(dataForReport$SnpSummary)
  setkey(dataForReport$SnpSummary,chromosome,position)
  #head(dataForReport$SnpSummary)

  dataForReport$SampleSummary <- row.summary(plinkRaw$genotypes) #takes about 5 minutes
  #head(dataForReport$SampleSummary)

  #####################################
  #hg19 reference allele
  #####################################
  if (checkRefAllele) {
    print(paste0("Checking if the major alleles are hg19 reference alleles:"))
    chrToCheck=dataForReport$SnpSummary$chromosome
    chrToCheck[chrToCheck==23]="X"
    chrToCheck[chrToCheck==24]="Y"
    chrToCheck=paste0("chr",chrToCheck)
    genomeRefAllele<-getSeq(Hsapiens,chrToCheck,dataForReport$SnpSummary$position,width=1, as.character=TRUE)
    majorAlleleAsRef<-rep(NA,totalSnpNum)
    majorAlleleAsRef[which(dataForReport$SnpSummary$allele.1==genomeRefAllele)]<-0
    majorAlleleAsRef[which(dataForReport$SnpSummary$allele.2==genomeRefAllele)]<-1
    print(paste0("Identified ",length(which(majorAlleleAsRef==1))," probes with Major allele as Reference Allele in Genome; ",length(which(majorAlleleAsRef==0))," probes with Minor allele as Reference Allele; ",length(which(is.na(majorAlleleAsRef)))," probes with No Reference Allele;"))

    dataForReport$SnpSummary$MajorAlleleAsRef<-majorAlleleAsRef
  }

  ######################################
  #work on duplicate SNPs
  ######################################
  if (!is.null(dupSnpFile)) {
    print(paste0("Analyzing Concordance of duplicate SNPs:"))
    dataForReport$DuplicateSnpsConcordance=snpsComp(snpMatrix1=plinkRaw$genotypes,dupSnpFile=dupSnpFile,sampleConcordanceCutoff=0.5)
  }

  ######################################
  #work on duplicate Samples
  ######################################
  if (!is.null(dupSampleFile)) {
    print(paste0("Analyzing Concordance of duplicate Samples:"))
    dataForReport$DuplicateSamplesConcordance=snpsComp(snpMatrix1=plinkRaw$genotypes,dupSampleFile=dupSampleFile,snpConcordanceCutoff=0.5)
  }

  #no plink outdir at this time
  #	######################################
  #	#sex check by plink
  #	######################################
  #	genderCheckFile<-paste0(plinkOutDir,"/genderCheck.sexcheck")
  #	if (!file.exists(genderCheckFile)) {
  #		print(paste0("Gender Check by Plink:"))
  #		plinkCmd<-paste0("plink --noweb --bfile ",file_path_sans_ext(plinkBed)," --maf 0.1 --check-sex --out ",file_path_sans_ext(genderCheckFile))
  #		system(plinkCmd)
  #	}
  #	dataForReport$GenderCheckTable<-read.table(genderCheckFile,header=T,as.is=T)
  #	dataForReport$GenderCheckTable<-dataForReport$GenderCheckTable[which(dataForReport$GenderCheckTable$STATUS=="PROBLEM"),]

  ######################################
  #Hardy-weinberg equilibrium
  ######################################
  #pvalue2sided=2*pnorm(-abs(z))
  hweZCol<-grep("\\z\\.HWE$",colnames(dataForReport$SnpSummary))
  if (length(hweZCol)>0) {
    print(paste0("Analyzing Hardy-weinberg equilibrium:"))
    temp<-as.data.table(2*pnorm(-abs(as.matrix(dataForReport$SnpSummary[,..hweZCol]))))
    colnames(temp)<-gsub("\\z\\.HWE$","p.HWE",colnames(temp))
    dataForReport$SnpSummary<-cbind(dataForReport$SnpSummary,temp)
  }

  #######################################
  #Ratio
  #######################################
  print(paste0("Analyzing Snps Ratios:"))
  dataForReport$SampleRatio=snpMatrixToRatio(plinkRaw$genotypes)

  #######################################
  #aimFile
  #######################################
  if (!is.null(aimFile)) {
    aimProbes<-readLines(aimFile)
    temp<-intersect(colnames(plinkRaw$genotypes),aimProbes)
    if (length(temp)>0) {
      print(paste0("Analyzing Race related AIM probes:"))
      print(paste0(length(temp)," AIM probes identified and will be used for PCA analysis."))
      aimPcaResult<-aimSnpMatrixPca(plinkRaw$genotypes[,temp])
      colnames(aimPcaResult)<-paste0("PC",1:10)
      dataForReport$AimPcaResult=aimPcaResult

    } else {
      break(paste0("aimFile provided, but no oveelaps between probes in this data and in aimFile."))
    }
  }

  return(dataForReport)
}


#' prepareDataForReport
#'
#' @param plinkBed 1212
#' @param outFile 21212
#'
#' @export
prepareDataForReport<-function(plinkBed,outFile=paste0(plinkBed,".GenoTypingQC.html"),raceFile=NULL) {
  dataForReport<-list()
  dataForReport$ReportRace=FALSE
  dataForReport$RaceToG1000Maf=c(A="EAS",ASIAN="EAS",B="AFR",BLACK="AFR",W="EUR",WHITE="EUR",H="AMR",I="SAS")
  dataForReport$MafOutlierCut=0.3

  #####################################
  #read plink data
  #####################################

  plinkRaw <- read.plink(plinkBed) #bim and fam will be added automaticlly
  #row.names(plinkRaw$genotypes)<-plinkRaw$fam$member

  totalSampleNum=nrow(plinkRaw$fam)
  totalSnpNum=nrow(plinkRaw$map)
  print(paste0("Plink file with ",totalSampleNum," samples and ",totalSnpNum," SNPs was loaded."))

  #####################################
  #probe summary and sample summary
  #####################################
  dataForReport$SnpSummary <- col.summary(plinkRaw$genotypes)
  dataForReport$SnpSummary<-cbind(plinkRaw$map,dataForReport$SnpSummary)
  setDT(dataForReport$SnpSummary)
  setkey(dataForReport$SnpSummary,chromosome,position)
  #head(dataForReport$SnpSummary)

  dataForReport$SampleSummary <- row.summary(plinkRaw$genotypes) #takes about 5 minutes
  #head(dataForReport$SampleSummary)


  ######################################
  #work on G1000 MAF
  ######################################
  if (!is.null(G1000MafFile)) {
    table1000GMerged <- fread(G1000MafFile,header=TRUE)
    colnames(table1000GMerged)<-paste0("G1000.",colnames(table1000GMerged))
    setkey(table1000GMerged,G1000.chr,G1000.pos)

    dataForReport$SnpSummary<-MergeG1000MafToSnpSummary(dataForReport$SnpSummary,table1000GMerged)
  }





  if (dataForReport$ReportRace) { #add race information
    dataForReport$SampleSummary<-cbind(dataForReport$SampleSummary,Race=dataForReport$SampleRace[,2])
  }

  if (dataForReport$ReportRace) { #Snp Summary for each RACE
    for (raceOne in dataForReport$RaceForReport) {
      snpSummaryRaceOne<-col.summary(plinkRaw$genotypes[which(dataForReport$SampleRace[,2]==raceOne),])
      colnames(snpSummaryRaceOne)<-paste0("Race.",raceOne,".",colnames(snpSummaryRaceOne))
      dataForReport$SnpSummary<-cbind(dataForReport$SnpSummary,snpSummaryRaceOne)
    }
  }
  dataForReport$SnpSummary<-cbind(plinkRaw$map,dataForReport$SnpSummary)
  setDT(dataForReport$SnpSummary)
  setkey(dataForReport$SnpSummary,chromosome,position)
  #head(dataForReport$SnpSummary)


  dataForReport$SampleSummary <- row.summary(plinkRaw$genotypes) #takes about 5 minutes
  if (dataForReport$ReportRace) { #add race information
    dataForReport$SampleSummary<-cbind(dataForReport$SampleSummary,Race=dataForReport$SampleRace[,2])
  }
  head(dataForReport$SampleSummary)




  if (!is.null(raceFile)) {
    SampleRaceRaw<-read.delim(raceFile,header=F,row.names=2,as.is=T)
    dataForReport$SampleRace<-SampleRaceRaw[plinkRaw$fam$member,]
    row.names(dataForReport$SampleRace)<-plinkRaw$fam$member

    print("Race table:")
    raceTable<-table(dataForReport$SampleRace[,2])
    print(raceTable)
    if (any(raceTable>=10)) {
      dataForReport$ReportRace=TRUE
      dataForReport$RaceForReport<-names(which(raceTable>=10)) #only interested in RACE with at least 5 samples
      dataForReport$SampleRace[,2][which(! dataForReport$SampleRace[,2] %in% dataForReport$RaceForReport)]<-NA

      #For 1000G MAF
      temp<-intersect(dataForReport$RaceForReport,names(dataForReport$RaceToG1000Maf))
      if (length(temp)<length(dataForReport$RaceForReport)) {
        warning(paste0("Different Race names in Race file and RaceToG1000Maf parameter. You may need to change one of them:\n Race provided: ",
                       paste0(dataForReport$RaceForReport,collapse=", "),"\n RaceToG1000Maf:", paste0(names(dataForReport$RaceToG1000Maf),collapse=", ")))
      }
      dataForReport$MafReportPairs<-cbind(paste0("G1000.",dataForReport$RaceToG1000Maf[temp]),paste0("Race.",temp,".MAF"))

    } else {
      dataForReport$SampleRace<-NULL
    }

    if (length(plinkRaw$fam$member)!=nrow(SampleRaceRaw)) {
      warning(paste0("There are xx samples in Plink file but Race file only have xx sampleswas provided with ",nrow(SampleRaceRaw)," samples, but ",length(plinkRaw$fam$member)-length(temp)," samples have no Race information"))
    }
  }

  return(plinkRaw)

}
