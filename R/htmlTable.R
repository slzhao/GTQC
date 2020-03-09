# TODO: Add comment
# 
# Author: zhaos6
###############################################################################


makeFormatTable<-function(dataForTable,normalizeMin=0.2,percentColInd=NULL,makeIntoDT=FALSE,maxRow=200) {
	require(formattable)
	require(DT)
	
  if (nrow(dataForTable)>maxRow) {
    warning(paste0(nrow(dataForTable)," rows in data table. Only the first ",maxRow," were kept in the report."))
    dataForTable<-dataForTable[1:maxRow,]
  }
  
	numColInd<-which(sapply(dataForTable,class)=="numeric" | sapply(dataForTable,class)=="integer")
	
	formatList<-list()
	for (i in 1:length(numColInd)) {
		formatList[[i]]<-as.formula(paste0('area(col = c(\'',colnames(dataForTable)[numColInd[i]],'\')) ~ normalize_bar("pink",min=',normalizeMin,',na.rm=TRUE)'))
	}
	if (!is.null(percentColInd)) {
		for (i in percentColInd) {
			dataForTable[,i]<-percent(dataForTable[,i])
		}
	}
  if (makeIntoDT) {
    htmlTable<-as.datatable(formattable(dataForTable,formatList),rownames=FALSE,
                 options = list(columnDefs = list(list(className = 'dt-right', targets = "_all")))
                 )
  } else {
	  htmlTable<-formattable(dataForTable,formatList)
  }
  return(htmlTable)
}


makeSummaryTable<-function(dataForTable,valueCut=c(0.01,0.5,0.99)) {
	valueCut<-sort(valueCut)
  
  processDataForTableOne<-function(dataForTableOne,valueCut) {
    dataForTableInCategory<-cut(dataForTableOne,c(-Inf,valueCut,Inf))
    dataForTableInCategoryTable<-table(dataForTableInCategory)
    temp1<-paste0("<=",valueCut)
    temp2<-cumsum(dataForTableInCategoryTable)[-length(dataForTableInCategoryTable)]
    temp3<-temp2/sum(as.numeric(dataForTableInCategoryTable))[-length(dataForTableInCategoryTable)]
    reusltOut<-data.frame(Threshold=temp1,Below=temp2,Percentage=temp3,stringsAsFactors=FALSE)
    row.names(reusltOut)<-NULL
    return(reusltOut)
  }
  
  if (is.vector(dataForTable)) {
    reusltOut<-processDataForTableOne(dataForTable,valueCut=valueCut)
  } else { #data.frame. more than 1 column
    reusltOut<-NULL
    for (i in 1:ncol(dataForTable)) {
      dataForTableOne<-unlist(dataForTable[,..i])
      reusltOutOne<-processDataForTableOne(dataForTableOne,valueCut=valueCut)[,-1]
      colnames(reusltOutOne)<-paste0(colnames(dataForTable)[i],".",colnames(reusltOutOne))
      if (is.null(reusltOut)) {
        reusltOut<-reusltOutOne
      } else {
        reusltOut<-cbind(reusltOut,reusltOutOne)
      }
    }
    Threshold<-processDataForTableOne(dataForTableOne,valueCut=valueCut)[,1]
    reusltOut<-cbind(Threshold,reusltOut)
  }
  return(reusltOut)
}


