#' @import tools
#' @import data.table
#' @import snpStats
NULL
preparePlinkData<-function(plinkBed) {
  dataForReport<-list()
  #####################################
  #read plink data
  #####################################

  plinkRaw <- read.plink(plinkBed) #bim and fam will be added automaticlly
  row.names(plinkRaw$genotypes)<-plinkRaw$fam$member

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
  row.names(plinkRaw$genotypes)<-plinkRaw$fam$member

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
