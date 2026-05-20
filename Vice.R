#Function to look at ordinal range change (in number of horizons) of a taxon
	#between two alternate GCMs of the same dataset
#y_1 is a taxa column of the original GCM
#y_2 is the corresponding column of a new GCM
#Negative value indicates range shrinking (number of horizons) in y_2 relative to y_1
#Positive value indicates range expansion (number of horizons) in y_2 relative to y_1
RangeChange<-function(y_1,y_2){
	rangePos<-range(which(y_1==1))
	rangePos2<-range(which(y_2==1))
	HorSpan<-rangePos[2]-rangePos[1]
	HorSpan2<-rangePos2[2]-rangePos2[1]
	#Negative value indicates range shrinking
	#Positive value indicates range expansion
	SpanDiff<-HorSpan2-HorSpan
	return(SpanDiff)
}

###########################################################################
#Function to loop through all taxa to compare range changes (in number of horizons) between two composites
#y1 is GCM 1 (original solution)
#y2 is GCM 2 (alternate solution)
#Negative value indicates range shrinking (number of horizons) in y2 relative to y1
#Positive value indicates range expansion (number of horizons) in y2 relative to y1
SumRangeChange<-function(y1,y2){
	nTaxa<-length(y1[1,])
	AllSpanDiff<-c(rep(NA,nTaxa-4))
	for(i in 5:nTaxa){
		AllSpanDiff[i-4]<-RangeChange(y1[,i],y2[,i])
	}
	print(AllSpanDiff)
	SumSpanDiff<-sum(AllSpanDiff)
	return(SumSpanDiff)
}

#############################################################################
#Function which moves each horizon, starting from the oldest horizon (=position 1 in GCM), up in the composite
	#until it either bumps up against the next horizon in its own section or it increases the solution score
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
TaxaShiftBottom_SingleHor<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nHor2<-nHor-1
	nTaxa<-length(y[1,])
	Shift_No_Change<-0
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone[InnerCount,2]!=GCM_Clone[InnerCount+1,2]){
			scoreTemp1<-GCM_Clone[InnerCount,1]
			scoreTemp2<-GCM_Clone[InnerCount+1,1]
			GCM_Clone[InnerCount,1]<-scoreTemp2
			GCM_Clone[InnerCount+1,1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			if(PenCheck>PenCur){
				Check<-1
				Control<-1
			}
			else{
				if(PenCheck<PenCur){
					print(paste("Better,",PenCheck,",",rowIndex,",",InnerCount,sep=""))
				}
				Shift_No_Change<-Shift_No_Change+1
				PenCur<-PenCheck
			}
			InnerCount<-InnerCount+1
			if(InnerCount==nHor){
				Check<-1
				Control<-0
			}
		}
		else{
			Control<-0
			Check<-1
		}
	}
	return(c(Shift_No_Change,Control))
}

#testShiftBot2<-TaxaShiftBottom_SingleHor(Dar,1)

#############################################################################
#Function to put a bottom-up vice on all horizons
	#Loop through TaxashiftBottom_SingleHor() function
#Values returned are number of ordinal positions that can be moved without
	#changing penalty or violating stratigraphy
TaxaShiftBottom_All<-function(y,PenType="Normal"){
	y<-as.data.frame(y)
	nHor<-length(y[,1])
	res<-data.frame(matrix(NA,nrow=nHor,ncol=2))
	names(res)<-c("BottomVice","Control")
	nHor2<-nHor-1
	for(i in 1:nHor2){
		#print(i)
		res[i,1:2]<-TaxaShiftBottom_SingleHor(y,i,PenType=PenType)
	}
	res[nHor,1]<-0
	return(res)
}

#############################################################################
#Function which moves each horizon, starting from the youngest horizon (=greatest position in GCM), down in the composite
	#until it either bumps up against the next horizon in its own section or it increases the solution score
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
TaxaShiftTop_SingleHor<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nTaxa<-length(y[1,])
	Shift_No_Change<-0
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone[InnerCount,2]!=GCM_Clone[InnerCount-1,2]){
			scoreTemp1<-GCM_Clone[InnerCount,1]
			scoreTemp2<-GCM_Clone[InnerCount-1,1]
			GCM_Clone[InnerCount,1]<-scoreTemp2
			GCM_Clone[InnerCount-1,1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			if(PenCheck>PenCur){
				Check<-1
				Control<-1
			}
			else{
				if(PenCheck<PenCur){
					#print(paste("Better:",PenCheck,rowIndex,InnerCount))
					print(paste("Better,",PenCheck,",",rowIndex,",",InnerCount,sep=""))
				}
				Shift_No_Change<-Shift_No_Change+1
				PenCur<-PenCheck
			}
			InnerCount<-InnerCount-1
			if(InnerCount==1){
				Check<-1
				Control<-0
			}
		}
		else{
			Check<-1
			Control<-0
		}
	}
	return(c(Shift_No_Change,Control))
}

#testShiftTop2<-TaxaShiftTop_SingleHor(Dar,1549)

#############################################################################
#Function to put a top-down vice on all horizons
	#Loop through TaxashiftTop_SingleHor() function
#Values returned are number of ordinal positions that can be moved without
	#changing penalty or violating stratigraphy
TaxaShiftTop_All<-function(y,PenType="Normal"){
	y<-as.data.frame(y)
	nHor<-length(y[,1])
	res<-data.frame(matrix(NA,nrow=nHor,ncol=2))
	names(res)<-c("TopVice","Control")
	for(i in nHor:2){
		#print(i)
		res[i,1:2]<-TaxaShiftTop_SingleHor(y,i,PenType=PenType)
	}
	res[1,1]<-0
	return(res)
}

#testShiftTop<-TaxaShiftTop_All(jTest$d,PenType="Normal")
#testShiftBot<-TaxaShiftBottom_All(jTest$d,PenType="Normal")

##############################################################################
#Function which moves each horizon, starting from the oldest horizon (=position 1 in GCM), up in the composite
	#until it either bumps up against the next horizon in its own section or it increases the solution score
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
TaxaShiftBottom_SingleHor_v2<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nHor2<-nHor-1
	nTaxa<-length(y[1,])
	Shift_No_Change<-0
	OriPos<-c()
	improvPos<-c()
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone$Section[InnerCount]!=GCM_Clone$Section[InnerCount+1]){
			scoreTemp1<-GCM_Clone$Score[InnerCount]
			scoreTemp2<-GCM_Clone$Score[InnerCount+1]
			GCM_Clone$Score[InnerCount]<-scoreTemp2
			GCM_Clone$Score[InnerCount+1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			if(PenCheck>PenCur){
				Check<-1
				Control<-1
			}
			else{
				if(PenCheck<PenCur){
					print(paste("Better,",PenCheck,",",rowIndex,",",InnerCount,sep=""))
					OriPos<-c(OriPos,rowIndex)
					improvPos<-c(improvPos,InnerCount)
				}
				Shift_No_Change<-Shift_No_Change+1
				PenCur<-PenCheck
			}
			InnerCount<-InnerCount+1
			if(InnerCount==nHor){
				Check<-1
				Control<-0
			}
		}
		else{
			Control<-0
			Check<-1
		}
	}
	return(list(Shift=c(Shift_No_Change,Control),BetterPos=c(OriPos,improvPos)))
}

#testShiftBot2_v2<-TaxaShiftBottom_SingleHor_v2(Dar,1)

#############################################################################
#Function to put a bottom-up vice on all horizons
	#Loop through TaxashiftBottom_SingleHor() function
#Values returned are number of ordinal positions that can be moved without
	#changing penalty or violating stratigraphy
TaxaShiftBottom_All_v2<-function(y,PenType="Normal"){
	y<-as.data.frame(y)
	nHor<-length(y[,1])
	res<-data.frame(matrix(NA,nrow=nHor,ncol=2))
	res2<-c(nHor,nHor)
	names(res)<-c("BottomVice","Control")
	nHor2<-nHor-1
	for(i in 1:nHor2){
		print(i)
		CurHorCalc<-TaxaShiftBottom_SingleHor_v2(y,i,PenType=PenType)
		res[i,1:2]<-CurHorCalc[[1]]
		if(is.null(CurHorCalc[[2]])==FALSE){
			res2<-as.data.frame(rbind(res2,CurHorCalc[[2]]))
		}
	}
	names(res2)<-c("OriPos","ImprovPos")
	res[nHor,1]<-0
	return(list(ShiftAll=res,PenImprov=res2))
}

#############################################################################
#Function which moves each horizon, starting from the youngest horizon (=greatest position in GCM), down in the composite
	#until it either bumps up against the next horizon in its own section or it increases the solution score
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
TaxaShiftTop_SingleHor_v2<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nTaxa<-length(y[1,])
	Shift_No_Change<-0
	OriPos<-c()
	improvPos<-c()
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone$Section[InnerCount]!=GCM_Clone$Section[InnerCount-1]){
			scoreTemp1<-GCM_Clone$Score[InnerCount]
			scoreTemp2<-GCM_Clone$Score[InnerCount-1]
			GCM_Clone$Score[InnerCount]<-scoreTemp2
			GCM_Clone$Score[InnerCount-1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			if(PenCheck>PenCur){
				Check<-1
				Control<-1
			}
			else{
				if(PenCheck<PenCur){
					#print(paste("Better:",PenCheck,rowIndex,InnerCount))
					print(paste("Better,",PenCheck,",",rowIndex,",",InnerCount,sep=""))
					OriPos<-c(OriPos,rowIndex)
					improvPos<-c(improvPos,InnerCount)
				}
				Shift_No_Change<-Shift_No_Change+1
				PenCur<-PenCheck
			}
			InnerCount<-InnerCount-1
			if(InnerCount==1){
				Check<-1
				Control<-0
			}
		}
		else{
			Check<-1
			Control<-0
		}
	}
	return(list(Shift=c(Shift_No_Change,Control),BetterPos=c(OriPos,improvPos)))
}

#testShiftTop2<-TaxaShiftTop_SingleHor_v2(Dar,1549)

#############################################################################
#Function to put a top-down vice on all horizons
	#Loop through TaxashiftTop_SingleHor() function
#Values returned are number of ordinal positions that can be moved without
	#changing penalty or violating stratigraphy
TaxaShiftTop_All_v2<-function(y,PenType="Normal"){
	y<-as.data.frame(y)
	nHor<-length(y[,1])
	res<-data.frame(matrix(NA,nrow=nHor,ncol=2))
	res2<-c(nHor,nHor)
	names(res)<-c("TopVice","Control")
	for(i in nHor:2){
		print(i)
		CurHorCalc<-TaxaShiftTop_SingleHor_v2(y,i,PenType=PenType)
		res[i,1:2]<-CurHorCalc[[1]]
		if(is.null(CurHorCalc[[2]])==FALSE){
			res2<-as.data.frame(rbind(res2,CurHorCalc[[2]]))
		}
	}
	names(res2)<-c("OriPos","ImprovPos")
	res[1,1]<-0
	return(list(ShiftAll=res,PenImprov=res2))
}

#testShiftTop2<-TaxaShiftTop_All_v2(jAll2$d,PenType="Normal")
#testShiftBot2<-TaxaShiftBottom_All_v2(jAll2$d,PenType="Normal")

##############################################################################
#Function to plot the results of TaxaShiftTop_All and TaxaShiftBot_All function
	#in ordinal position
#Plots as a moving average using the TTR package
#Top-down vice position changes are displayed as blue negative values
#Bottom-up vice position changes are displayed as red positive values
#windowSize is the moving average window size, objective choice is 2*SD of jackknife results
	#Hirn set 2SD = 40, Darr set 2SD = 16
VicePlot<-function(ViceTop,ViceBot,windowSize=20){
	library(TTR)
	smoothTop<-SMA(ViceTop,windowSize)
	smoothBot<-SMA(ViceBot,windowSize)
	nLevel<-length(ViceTop)
	xtxt<-paste("Vice Motion Allowed (",windowSize," point moving average)",sep="")
	plot(smoothTop*-1,1:nLevel,xlim=c(round(floor(min(smoothTop*-1,na.rm=TRUE)),0),round(ceiling(max(smoothBot,na.rm=TRUE)),0)),col="blue",pch=20,type="b",xlab=xtxt,ylab="Ordinal Position")
	par(new=TRUE)
	plot(smoothBot,1:nLevel,xlim=c(round(floor(min(smoothTop*-1,na.rm=TRUE)),0),round(ceiling(max(smoothBot,na.rm=TRUE)),0)),col="red",pch=20,type="b",xlab="",ylab="")
	legend("bottomleft",c("TopVice","BottomVice"),pch=c(20,20),col=c("blue","red"))
}

#VicePlot(vic[,2],vic[,3],16)
#abline(h=c(386,1387),lwd=1.5)
#text(x=-70,y=405,"Darriwilian GSSP")

#pdf("ViceAction_smooth16.pdf")
#VicePlot(testShiftTop[,1],testShiftBot[,1],windowSize=16)
#abline(h=c(386,1387),lwd=1.5)
#text(x=-60,y=410,"Darriwilian GSSP")
#text(x=60,y=1410,"Sandbian GSSP")
#dev.off()

#############################################################################
#Function in which a given horizon is moved up in the composite until it bumps up
	#against the next horizon in its own section
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#Lists the penalty score for each move
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
Bottom_SingleHor_BFI<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nHor2<-nHor-1
	nTaxa<-length(y[1,])
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	Shift_Pen<-c(PenCur)
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	HorList<-c(InnerCount)
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone[InnerCount,2]!=GCM_Clone[InnerCount+1,2]){
			scoreTemp1<-GCM_Clone[InnerCount,1]
			scoreTemp2<-GCM_Clone[InnerCount+1,1]
			GCM_Clone[InnerCount,1]<-scoreTemp2
			GCM_Clone[InnerCount+1,1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			Shift_Pen<-c(Shift_Pen,PenCheck)
			HorList<-c(HorList,InnerCount)
			InnerCount<-InnerCount+1
			if(InnerCount==nHor){
				Check<-1
				Control<-0
			}
		}
		else{
			Control<-0
			Check<-1
		}
	}
	return(cbind(HorList,Shift_Pen))
}

#testBFIBotVice<-Bottom_SingleHor_BFI(jTest$d,1368)

#############################################################################
#Function in which a given horizon is moved down in the composite until it bumps up
	#against the next horizon in its own section
#Useful for approximating the difference that is caused by not allowing ties in HA, heuristic error
	#also occasionally finds better solution but does not currently adopt those if found, equivalent of point shifts in HA
#Lists the penalty score for each move
#y is the GCM of a solution
#rowIndex is the horizon of interest
#PenType indicates whether the normal or strict penalty is to be used for checking penalty
	#this function assumes taxon records begin at column 5
Top_SingleHor_BFI<-function(y,rowIndex,PenType="Normal"){
	nHor<-length(y[,1])
	nTaxa<-length(y[1,])
	if(PenType=="Normal"){
		PenCur<-NetRangeExtension(y[,5:nTaxa])
	}
	if(PenType=="Strict"){
		PenCur<-NetRangeExtension_v2(y[,5:nTaxa])
	}
	Shift_Pen<-c(PenCur)
	GCM_Clone<-y
	Check<-0
	InnerCount<-rowIndex
	HorList<-c(InnerCount)
	#Control is a dummy code for what is preventing movement
		# 0 = stratigraphic, 1 = biology (penalty increase)
	Control<-(-1)
	while(Check==0){
		if(GCM_Clone[InnerCount,2]!=GCM_Clone[InnerCount-1,2]){
			scoreTemp1<-GCM_Clone[InnerCount,1]
			scoreTemp2<-GCM_Clone[InnerCount-1,1]
			GCM_Clone[InnerCount,1]<-scoreTemp2
			GCM_Clone[InnerCount-1,1]<-scoreTemp1
			GCM_Clone<-GCM_Clone[order(GCM_Clone[,1]),]
			if(PenType=="Normal"){
				PenCheck<-NetRangeExtension(GCM_Clone[,5:nTaxa])
			}
			if(PenType=="Strict"){
				PenCheck<-NetRangeExtension_v2(GCM_Clone[,5:nTaxa])
			}
			Shift_Pen<-c(Shift_Pen,PenCheck)
			HorList<-c(HorList,InnerCount)
			InnerCount<-InnerCount-1
			if(InnerCount==1){
				Check<-1
				Control<-0
			}
		}
		else{
			Check<-1
			Control<-0
		}
	}
	return(cbind(HorList,Shift_Pen))
}

#testBFITopVice<-Top_SingleHor_BFI(jTest$d,1368)

########################################################################################
#Function to plot a taxa shift range chart
	#Each x-value is a single horizon occurrence of a taxon
	#Excursion bars are the results of the top and bottom vices for the given horizon
		#The number of ordinal positions a horizon can shift without increasing the penalty
			#or violating stratigraphy
	#SectLabel adds the section number the horizon belongs to above each horizon entry
		#off by default
PlotTaxaViceRange<-function(GCM,myDict,taxIndex,ShiftTop,ShiftBot,SectLabel=FALSE){
	taxPos<-which(GCM[,taxIndex]==1)
	nPos<-length(taxPos)
	BotRange<-min(taxPos-ShiftTop$ShiftAll$TopVice[taxPos])
	TopRange<-max(taxPos+ShiftBot$ShiftAll$BottomVice[taxPos])
	plot(1:nPos,taxPos,pch=20,ylim=c(round(BotRange,-1)-10,round(TopRange,-1)+10),xlab="",ylab="Ordinal Shift Range")
	arrows(x0=1:nPos,y0=taxPos,x1=1:nPos,y1=taxPos-ShiftTop$ShiftAll$TopVice[taxPos],length=0.05,angle=90,code=2)
	arrows(x0=1:nPos,y0=taxPos,x1=1:nPos,y1=taxPos+ShiftBot$ShiftAll$BottomVice[taxPos],length=0.05,angle=90,code=2)
	if(SectLabel==TRUE){
		text(x=1:nPos,y=taxPos+ShiftBot$ShiftAll$BottomVice[taxPos]+10,c(GCM[taxPos,2]),cex=0.5,srt=90)
	}
	legend("topleft",myDict$TaxaName[which(myDict$GRCode==names(GCM[1,])[taxIndex])],bty="n")
}

#PlotTaxaViceRange(j$d,myDict=myDict,ShiftTop=testShiftTop,ShiftBot=testShiftBot,SectLabel=TRUE,taxIndex=125)

########################################################################################
#Function to plot a taxa shift range chart
	#Each x-value is a single horizon occurrence of a section
	#Excursion bars are the results of the top and bottom vices for the given horizon
		#The number of ordinal positions a horizon can shift without increasing the penalty
			#or violating stratigraphy
	#SectList allows legend to display the section name rather than number
PlotSectionViceRange<-function(GCM,SectNumber,ShiftTop,ShiftBot,SectList){
	SectPos<-which(GCM[,2]==SectNumber)
	nPos<-length(SectPos)
	BotRange<-min(SectPos-ShiftTop$ShiftAll$TopVice[SectPos])
	TopRange<-max(SectPos+ShiftBot$ShiftAll$BottomVice[SectPos])
	plot(1:nPos,SectPos,pch=20,ylim=c(round(BotRange,-1)-10,round(TopRange,-1)+10),xlab="",ylab="Ordinal Shift Range")
	arrows(x0=1:nPos,y0=SectPos,x1=1:nPos,y1=SectPos-ShiftTop$ShiftAll$TopVice[SectPos],length=0.05,angle=90,code=2)
	arrows(x0=1:nPos,y0=SectPos,x1=1:nPos,y1=SectPos+ShiftBot$ShiftAll$BottomVice[SectPos],length=0.05,angle=90,code=2)
	legend("topleft",SectList$SecName[which(SectList$SecNum==SectNumber)],bty="n")
}

#PlotSectionViceRange(j$d,ShiftTop=testShiftTop,ShiftBot=testShiftBot,SectNumber=1012,SectList=SectList)

########################################################################################
#Function to plot a taxa shift range chart
	#Each x-value is a single horizon occurrence of a taxon
	#Excursion bars are the results of the top and bottom vices for the given horizon
		#The number of ordinal positions a horizon can shift without increasing the penalty
			#or violating stratigraphy
	#SectLabel adds the section number the horizon belongs to above each horizon entry
		#off by default
PlotTaxaJackknifeRange<-function(GCM,taxIndex,AllJack,CI=0.95,myDict){
	taxPos<-which(GCM[,taxIndex]==1)
	AllJack<-AllJack[order(AllJack$Original.Position),]
	nPos<-length(taxPos)
	nJack<-length(AllJack[1,])
	UpperEdge<-round((nJack-5)*CI,0)
	LowerEdge<-round((nJack-5)*(1-CI),0)
	LowJack<-c(rep(NA,nPos))
	UpJack<-c(rep(NA,nPos))
	for(i in 1:nPos){
		JackTemp<-sort(AllJack[taxPos[i],5:nJack])
		LowJack[i]<-JackTemp[LowerEdge][[1]]
		UpJack[i]<-JackTemp[UpperEdge][[1]]
	}
	BotRange<-min(AllJack[taxPos,5:nJack],na.rm=TRUE)
	TopRange<-max(AllJack[taxPos,5:nJack],na.rm=TRUE)
	plot(1:nPos,taxPos,pch=20,ylim=c(round(BotRange,-1)-10,round(TopRange,-1)+10),xlab="",ylab="Ordinal Jackknife Range")
	arrows(x0=1:nPos,y0=taxPos,x1=1:nPos,y1=LowJack,length=0.05,angle=90,code=2)
	arrows(x0=1:nPos,y0=taxPos,x1=1:nPos,y1=UpJack,length=0.05,angle=90,code=2)
	legend("topleft",myDict$TaxaName[which(myDict$GRCode==names(GCM[1,])[taxIndex])],bty="n")
}

#PlotTaxaJackknifeRange(Dar,taxIndex=125,AllJack=AllJack,CI=0.95,myDict=myDict)

#######################################################################################
#Computes the largest gap in ordinal space for each taxon in a GCM
	#Assumes taxa start on column 5
BiggestGap<-function(GCM){
	nTaxa<-length(GCM[1,])
	Gaps<-c(rep(0,nTaxa-4))
	for(i in 5:nTaxa){
		taxPos<-which(GCM[,i]==1)
		nPos<-length(taxPos)
		if(nPos>1){
			for(j in 2:nPos){
				tempMax<-0
				tempDif<-taxPos[j]-taxPos[j-1]
				if(tempDif>tempMax){
					tempMax<-tempDif
				}
			}
		}
		Gaps[i-4]<-tempMax
	}
	return(Gaps)
}

#testGap<-BiggestGap(j$d)

####################################################################################
#Function to optimize a GCM given improvments found from vice output
	##NOT WORKING CURRENTLY
OptimizeViceGCM<-function(GCM,ViceOut){
	nHors<-length(GCM[,1])
	nImprov<-length(ViceOut$PenImprov[,1])
	tempScore<-GCM[,1]
	#Bottom vice
	if(ViceOut$PenImprov[2,1]-ViceOut$PenImprov[2,2]<0){
		for(i in 2:nImprov){
			print(i)
			MoveToPos<-which(GCM[,1]==tempScore[ViceOut$PenImprov[i,2]])
			MoveFromPos<-which(GCM[,1]==tempScore[ViceOut$PenImprov[i,1]])
			tempScore[MoveToPos]<-tempScore[MoveFromPos]
			Swap1<-MoveFromPos+1
			Swap2<-MoveFromPos-1
			tempScore[Swap1:MoveToPos]<-tempScore[MoveFromPos:Swap2]
		}
	}
	#Top Vice
	else{
		for(i in 2:nImprov){
			print(i)
			MoveToPos<-which(GCM[,1]==tempScore[ViceOut$PenImprov[i,2]])
			MoveFromPos<-which(GCM[,1]==tempScore[ViceOut$PenImprov[i,1]])
			tempScore[MoveToPos]<-tempScore[MoveFromPos]
			Swap1<-MoveFromPos-1
			Swap2<-MoveToPos+1
			tempScore[MoveToPos:Swap1]<-tempScore[Swap2:MoveFromPos]
		}
	}
	GCM[,1]<-tempScore
	GCM<-GCM[order(GCM[,1]),]
	return(GCM)
}

#testOpt<-OptimizeViceGCM(Dar,testShiftTop2)
#testOpt<-OptimizeViceGCM(testOpt,testShiftBot2)

####################################################################################
#Aeronian GSSP
	#Trefawr Track (#)
	#Monograptus austerus sequens (GR1815)

#Rhudannian GSSP
	#Dobs Linn (#461)
	#Akidograptus ascensus (GR1279)

#Hirnantian GSSP
	#Wangjiwan North (#153)
	#Metabolograptus extraordinarius (GR4919)

#Katian GSSP
	#Black Knob Ridge (#131)
	#Diplcanthograptus caudatus (GR3301)

#Sandbian GSSP
	#Fagelsang (#1042)
	#Nemagraptus gracilis gracilis (GR1035)

#Darriwilian GSSP
	#Huangnitang (#1012)
	#Levisograptus austrodentatus austrodentatus (GR5090)

#Dapingian GSSP
	#Huanghuachang (#1062)
	#Baltoniodus triangularis (CO1612)

#Floian GSSP
	#Diabassbrottet Quarry (#517)
	#Paratetragraptus approximatus (GR5397)

#Tremadoc GSSP
	#Green Point (#79,#100,#486)
	#Iapetognathus fluctivagus (CO1717)
