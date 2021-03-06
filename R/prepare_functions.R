#' @export
processingPlinkData=function(plinkBed,checkRefAllele=TRUE,dupSnpFile=NULL,dupSampleFile=NULL,raceFile=NULL,aimFile=NULL,
                             MafOutlierCut=0.3,
                             plinkChrToGenomeChr=c("23"="X","24"="Y","26"="MT"),
                             raceTo1000GRace=c("CHB"="EAS","JPT"="EAS","CHS"="EAS","CDX"="EAS","KHV"="EAS",
                                               "CEU"="EUR","TSI"="EUR","FIN"="EUR","GBR"="EUR","IBS"="EUR",
                                               "YRI"="AFR","LWK"="AFR","GWD"="AFR","MSL"="AFR","ESN"="AFR","ASW"="AFR","ACB"="AFR",
                                               "MXL"="AMR","PUR"="AMR","CLM"="AMR","PEL"="AMR",
                                               "GIH"="SAS","PJL"="SAS","BEB"="SAS","STU"="SAS","ITU"="SAS")
) {

  dataForReport=list(MafOutlierCut=MafOutlierCut)

  #####################################
  #read plink data
  #####################################
  #remove all "." probes/snps
  validSnps=data.table::fread(gsub("bed$","bim",plinkBed))
  validSnps=validSnps$V2[validSnps$V2!="."]

  dataForReport$plinkRaw <- read.plink(plinkBed,select.snps = validSnps) #bim and fam will be added automaticlly
  totalSampleNum=nrow(dataForReport$plinkRaw$fam)
  totalSnpNum=nrow(dataForReport$plinkRaw$map)
  print(paste0("Plink file with ",totalSampleNum," samples and ",totalSnpNum," SNPs was loaded."))

  #####################################
  #probe summary and sample summary
  #####################################
  dataForReport$SnpSummary <- col.summary(dataForReport$plinkRaw$genotypes)
  dataForReport$SnpSummary<-cbind(dataForReport$plinkRaw$map,dataForReport$SnpSummary)
  #head(dataForReport$SnpSummary)

  dataForReport$SampleSummary <- row.summary(dataForReport$plinkRaw$genotypes) #takes about 5 minutes
  #head(dataForReport$SampleSummary)

  ######################################
  #work on G1000 MAF
  ######################################
  G1000Population=GenomicScores::populations(MafDb.1Kgenomes.phase3.hs37d5::MafDb.1Kgenomes.phase3.hs37d5)
  snpGRange=dataForReport$SnpSummary[,c("chromosome","position")]
  snpGRange$chromosome=mapValues(snpGRange$chromosome,plinkChrToGenomeChr)  #change chromosome name format

  snpPositionNaInd=which(is.na(snpGRange$position))
  snpPositionNotNaInd=which(!is.na(snpGRange$position))
  if (length(snpPositionNaInd)>0) { #fill NA with chr1 and pos 1, so that can be used to get 1000G MAF
    snpGRange$chromosome[snpPositionNaInd]="1"
    snpGRange$position[snpPositionNaInd]=1
  }
  snpGRange=GenomicRanges::makeGRangesFromDataFrame(snpGRange, start.field="position",end.field="position")
  snpG1000MafTable=GenomicRanges::mcols(GenomicScores::gscores(MafDb.1Kgenomes.phase3.hs37d5::MafDb.1Kgenomes.phase3.hs37d5, snpGRange,pop=G1000Population))
  snpG1000MafTable=data.table(data.frame(snpG1000MafTable))
  dataForReport$SnpSummary=cbind(dataForReport$SnpSummary,snpG1000MafTable)

  setDT(dataForReport$SnpSummary)
  setkey(dataForReport$SnpSummary,chromosome,position)

  ######################################
  #work on Race
  ######################################
  if (!is.null(raceFile)) {
    SampleRaceRaw<-read.delim(raceFile,header=F,row.names=1,as.is=T)
    dataForReport$SampleSummary$Race=SampleRaceRaw[row.names(dataForReport$SampleSummary),1]
    dataForReport$SampleSummary$Race=mapValues(dataForReport$SampleSummary$Race,raceTo1000GRace)
    print("Race table:")
    raceTable<-table(dataForReport$SampleSummary$Race)
    print(raceTable)

    if (any(raceTable>=10)) {
      dataForReport$ReportRace=TRUE
      dataForReport$RaceForReport<-names(which(raceTable>=10)) #only interested in RACE with at least 10 samples
      dataForReport$SampleSummary$Race[which(! dataForReport$SampleSummary$Race %in% dataForReport$RaceForReport)]<-NA

      for (raceOne in dataForReport$RaceForReport) { #Snp Summary for each RACE
        snpSummaryRaceOne<-col.summary(dataForReport$plinkRaw$genotypes[which(dataForReport$SampleRace[,2]==raceOne),])
        colnames(snpSummaryRaceOne)<-paste0("Race.",raceOne,".",colnames(snpSummaryRaceOne))
        dataForReport$SnpSummary<-cbind(dataForReport$SnpSummary,snpSummaryRaceOne)
      }

      #Race MAF vs 1000G Race MAF
      temp<-intersect(dataForReport$RaceForReport,gsub("_AF","",G1000Population))
      if (length(temp)<length(dataForReport$RaceForReport)) {
        warning(paste0("Different Race names in Race file and RaceToG1000Maf parameter. You may need to change one of them:\n Race provided: ",
                       paste0(dataForReport$RaceForReport,collapse=", "),"\n RaceToG1000Maf:", paste0(names(dataForReport$RaceToG1000Maf),collapse=", ")))
      }
      dataForReport$MafReportPairs<-cbind(paste0(temp,"_AF"),paste0("Race.",temp,".MAF"))

      sampleNoRaceNum=length(setdiff(row.names(dataForReport$SampleSummary),row.names(SampleRaceRaw)))
      if (sampleNoRaceNum>0) {
        warning(paste0("There are ",nrow(dataForReport$SampleSummary)," samples in Plink file but ",sampleNoRaceNum," samples have no Race information in Race file: ",raceFile))
      }

    } else {
      dataForReport$ReportRace=FALSE
    }
  }

  #####################################
  #hg19 reference allele
  #####################################
  if (checkRefAllele) {
    print(paste0("Checking if the major alleles are hg19 reference alleles:"))
    chrToCheck=dataForReport$SnpSummary$chromosome
    snpPositionNotNaInd=which(!is.na(chrToCheck))
    #chrToCheck[chrToCheck==23]="X"
    #chrToCheck[chrToCheck==24]="Y"
    chrToCheck=mapValues(chrToCheck,plinkChrToGenomeChr)
    chrToCheck=paste0("chr",chrToCheck)

    genomeRefAllele=rep(NA,nrow(dataForReport$SnpSummary))
    genomeRefAllele[snpPositionNotNaInd]<-BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,chrToCheck[snpPositionNotNaInd],dataForReport$SnpSummary$position[snpPositionNotNaInd],width=1, as.character=TRUE)
    majorAlleleAsRef<-rep(NA,nrow(dataForReport$SnpSummary))
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
    dataForReport$DuplicateSnpsConcordance=snpsComp(snpMatrix1=dataForReport$plinkRaw$genotypes,dupSnpFile=dupSnpFile,sampleConcordanceCutoff=0.5)
  }

  ######################################
  #work on duplicate Samples
  ######################################
  if (!is.null(dupSampleFile)) {
    print(paste0("Analyzing Concordance of duplicate Samples:"))
    dataForReport$DuplicateSamplesConcordance=snpsComp(snpMatrix1=dataForReport$plinkRaw$genotypes,dupSampleFile=dupSampleFile,snpConcordanceCutoff=0.5)
  }

  ######################################
  #sex check by plink
  ######################################
  genderCheckFile<-paste0(dirname(plinkBed),"/",file_path_sans_ext(basename(plinkBed)),".sexcheck")
  if (Sys.which("plink")!="" | file.exists(genderCheckFile)) {
    if (!file.exists(genderCheckFile)) {
      print(paste0("Gender Check by Plink:"))
      genderCheckFile<-paste0(tempdir(),"/genderCheck.sexcheck")
      plinkCmd<-paste0("plink --noweb --bfile ",file_path_sans_ext(plinkBed)," --maf 0.1 --check-sex --out ",file_path_sans_ext(genderCheckFile))
      system(plinkCmd)
    }
    dataForReport$GenderCheckTable<-read.table(genderCheckFile,header=T,as.is=T)
    dataForReport$GenderCheckTable<-dataForReport$GenderCheckTable[which(dataForReport$GenderCheckTable$STATUS=="PROBLEM"),]
  } else {
    warning(paste0("Can't find plink in path nor plink sexcheck result file. Skip gender check."))
  }

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
  dataForReport$SampleRatio=snpMatrixToRatio(dataForReport$plinkRaw$genotypes)

  #######################################
  #aimFile
  #######################################
  if (!is.null(aimFile)) {
    aimProbes<-readLines(aimFile)
    aimProbes<-intersect(colnames(dataForReport$plinkRaw$genotypes),aimProbes)
    if (length(aimProbes)==0) {
      break(paste0("aimFile provided, but no oveelaps between probes in this data and in aimFile."))
    }
  } else {
    aimProbes=colnames(dataForReport$plinkRaw$genotypes)
    warning(paste0("aimFile not provided. Use all probes for PCA analysis."))
  }
  print(paste0("Analyzing Race related AIM probes:"))
  print(paste0(length(aimProbes)," AIM probes identified and will be used for PCA analysis."))
  aimPcaResult<-aimSnpMatrixPca(dataForReport$plinkRaw$genotypes[,aimProbes])
  colnames(aimPcaResult)<-paste0("PC",1:10)
  dataForReport$AimPcaResult=aimPcaResult

  #######################################
  #summaryTable
  #######################################
  dataForReport$SummaryTable=summaryTable(dataForReport)

  return(dataForReport)
}


