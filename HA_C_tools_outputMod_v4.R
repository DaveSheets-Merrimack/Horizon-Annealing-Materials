# HA formatting from OAL section files
# Modified by J. Boyle to cut out 9999s from GCM and to list GRCodes,
	#section numbers, and heights for each horizon into the d matrix
	#Also modified to allow cutoff for taxa occurrences to be specified
	#Added HAFile4R_v3 that allows use of biozone taxa for scaffolding starting solution

HAFile<-function(mySections)
{
	# create the HA Dmat and initial score matrix
	myTaxa=BGTaxaList(mySections)
	print(myTaxa)
	mySectionNames=BGSectionList(mySections)
	nsections=length(mySections)
	vTaxa=myTaxa[myTaxa$Occurrences>1,]	
	vTaxa=vTaxa[vTaxa$GRCode!="99999",]
	print(vTaxa)		
	#list of informative taxa
	ntaxa=length(vTaxa$GRCode)
	
	print(ntaxa)
	print(-1)
	TaxaNumbers=1:ntaxa
	names(TaxaNumbers)=vTaxa$GRCode
	
	hcount=rep(0,nsections)						
	# set up counts of the number of horizons in each section
	for(i in 1:nsections)
	{
		print(mySections[[i]]@NumHorizons)
		hcount[i]=mySections[[i]]@NumHorizons
	}
	hOffset=rep(0,nsections)					
	# set up the offsets
	for(i in 2:nsections)
	{
		hOffset[i]=hOffset[i-1]+hcount[(i-1)]
	}
	# total horizon count
	nHorizonsTotal=sum(hcount)
	# set up Dmat structure for HA
	# column 1- section number
	# column 2- horizon number
	# column 3-should be horizon height-for the moment, simply using horizon number
	# column 4 to ntaxa+3-  taxa
	
	# empty data matrix---filled with (-1) values....
    dmat=matrix(rep(-1,(ntaxa+2)*nHorizonsTotal), nrow=nHorizonsTotal,ncol=(ntaxa+3))
    
    for(i in 1:nsections)
    {
    	# fill in section number
    	print(paste(hOffset[i], " ",hcount[i]))
    	dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),1]=rep(i,hcount[i])
    	dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),2]=(1:hcount[i])
    	dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
	#Attempt to assign actual heights
	#for(t in 1:hcount){	
    		#dmat[(hOffset[t]+1),3]=mySections[[i]]@Height[t]
	#}
    	for(j in 1:mySections[[i]]@NumTaxa)
    	{
    		# check to see if taxa is valid
    		valid=1
    		if(length(grep("9999",mySections[[i]]@GRCode[j]))>0)				
    		# exclude 9999 entries
    		{
    			valid=0
    		}
    		if(sum(vTaxa$GRCode==mySections[[i]]@GRCode[j])<1)				
    		# check to see if taxa are one valid list
    		{
    			valid=0
    		}
    		print(valid)
    		print(mySections[[i]]@GRCode[j])
    		print(TaxaNumbers[mySections[[i]]@GRCode[j]])
    		if(valid)														
    		# taxa are correct, place in matrix
    		{
    			print(i)
    			print(j)
    			print(hcount[i])
    			print(hOffset[i])  
    			print(length(mySections[[i]]@Dmat[j, ]))
    			dmat[((hOffset[i]+1):(hOffset[i]+hcount[i])),(TaxaNumbers[mySections[[i]]@GRCode[j]]+3)]=mySections[[i]]@Dmat[j, ]
    		}
    	}
    }	
    
    # NOTE- as of 7/21/14, dmat is not complete- 2nd column needs to be the horizon height info.......
    
    
    # now set up starting solution
    # column 1 is the score,  column 2 is the section, column 3 is the horizon, sorted on column 1
    score=matrix(0,nrow=nHorizonsTotal,ncol=3)
    for(i in 1:nsections)
    {
    	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= sort(runif(hcount[i]))
    	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),2]=i
    	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
    }
    sc_order=order(score[,1])
    score=score[sc_order,]
    
	myBack=list(dmat=dmat,vTaxa=vTaxa,sectionList=mySectionNames,startScore=score)
}

############################################################################################
# set up the R based HA input file, from a set of OAL input sections
	#Added cutoff parameter to allow user to say how many occurrences a taxa must have to be included in composite
  
HAFile4R2<-function(mySections,cutoff=1){
	# create the HA Dmat and initial score matrix
	myTaxa=BGTaxaList(mySections)
	mySectionNames=BGSectionList(mySections)
	nsections=length(mySections)
	
	# note filtering here to taxa with more than a certain number of occurrence-saves time, but be a bit careful
	#vTaxa=myTaxa[myTaxa$Occurrences>cutoff,]
	vTaxa=myTaxa[myTaxa$Occurrences>cutoff | substring(myTaxa$GRCode,1,2)=="EB" | substring(myTaxa$GRCode,1,2)=="KB",]
	vTaxa=vTaxa[vTaxa$GRCode!="99999",]
	torder=order(vTaxa$GRCode)
	vTaxa=vTaxa[torder,]	
	#list of informative taxa
	ntaxa=length(vTaxa$GRCode)
	
	TaxaNumbers=1:ntaxa
	names(TaxaNumbers)=vTaxa$GRCode
	
	hcount=rep(0,nsections)						
	# set up counts of the number of horizons in each section
	for(i in 1:nsections)
	{
		hcount[i]=mySections[[i]]@NumHorizons
	}
	hOffset=rep(0,nsections)					
	# set up the offsets
	for(i in 2:nsections)
	{
		hOffset[i]=hOffset[i-1]+hcount[(i-1)]
	}
	# total horizon count
	nHorizonsTotal=sum(hcount)
	# set up Dmat structure for HA
	# column 1- initial score
	# column 2- section number
	# column 3- horizon number
	# column 4-should be horizon height-for the moment, simply using horizon number
	# column 5 to ntaxa+3-  taxa
	
	# empty data matrix---filled with (-1) values....
	dmat=matrix(rep(-1,(ntaxa+4)*nHorizonsTotal), nrow=nHorizonsTotal,ncol=(ntaxa+4))
    
	for(i in 1:nsections){
    		# fill in section number
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),2]=rep(mySections[[i]]@SectionNumber,hcount[i])
    		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),4]=mySections[[i]]@Height[1:mySections[[i]]@NumHorizons]
		#Output horizon heights as column 4 of dmat, added by JBoyle March 2nd 2016
    		for(j in 1:mySections[[i]]@NumTaxa){
    			# check to see if taxa is valid
    			valid=1
			# exclude 9999 entries
    			if(length(grep("9999",mySections[[i]]@GRCode[j]))>0){	
    				valid=0
    			}
			# check to see if taxa are on valid list
    			if(sum(vTaxa$GRCode==mySections[[i]]@GRCode[j])<1){	
    				valid=0
    			}
			# taxa are correct, place in matrix
    			if(valid){
    				dmat[((hOffset[i]+1):(hOffset[i]+hcount[i])),(TaxaNumbers[mySections[[i]]@GRCode[j]]+4)]=mySections[[i]]@Dmat[j, ]
    			}
		}
	}    
    
    # now set up starting solution
    # column 1 is the score,  column 2 is the section, column 3 is the horizon, colum 4 is the height sorted on column 1
    score=matrix(0,nrow=nHorizonsTotal,ncol=4)
    for(i in 1:nsections)
    {
    	temp=sort(runif(hcount[i]))
    	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    	dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp
	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),2]=mySections[[i]]@SectionNumber
    	score[(1+hOffset[i]):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
	#######Supposed to allow actual horizon heights to be included, currently not working#############
	#for(k in 1:hcount[i]){
	#	score[(1+hOffset[k]),4]=mySections[[i]]@Height[k]
	#}
    }
    #Gives index values of scores ordered from least to greatest
    sc_order=order(score[,1])
    #Reorderes scores into ascending order
    score=score[sc_order,]
    
    
    # adding column names to this matrix did not slow the calculation down
    # adding row names slowed the R version of HA by 10 to 12 percent
    # so do not use rownames
    
    addColumnNames=TRUE
    if(addColumnNames)
    {
    	temp=c("Score","Section","Horizon","Height",as.character(vTaxa$GRCode))
    	colnames(dmat)<-temp
    #	temp=as.character(mySectionNames$SecName)
    #	rownames(dmat)<-temp[dmat[,2]]
    }
    
    # set up the default style penaltry structure, assuming all taxa carry a penalty, but
    # that there are no other types of data
    # the penalty structure may have to be altered separately to allow other designs
    
    n_biostrat=ntaxa
	biostrat=1:ntaxa
	# for biostrat data, or Taxa FADs,LADs, the biostrat variable is the numbers of the columns with taxa

	n_pmag=0
	pmag=1

	# pmag is a list of the column(s) with paleomagnetic signals, or really any binary data, NA values are not counted

	n_dates=0
	dates=matrix(c(64,0,65,0,1000,65,0,66,0,1000,66,0,67,0,1000),nrow=3,byrow=TRUE)

	# each row of the dates matrix is a set of data to be entered into the passing penalty
	# the first entry on each row is the column of the lower variale
	# second entry on a row is the data type 0- singular date,  1- FAD, 2-LAD
	# third and fourth entries on each row are the column and type of the second variable
	# fifth value on each row is the weight

	n_ashes=0
	ashes=matrix(c(68,100,69,100),nrow=2,byrow=TRUE)
	n_continuous=0
	continuous=matrix(c(70,5,71,5),nrow=2,byrow=TRUE)
					PenaltySpec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,ashes=ashes,n_continuous=n_continuous,continuous=continuous)
    
    
	#Clean dmat by deleting horizons with only 9999 taxa
	displacedLast<-length(dmat[1,])
	UseList<-c()
	for (i in 1:nHorizonsTotal){
		if(max(dmat[i,5:displacedLast])==1){
			UseList<-c(UseList,i)
		}
	}
	dmat<-dmat[UseList,]

	# set up the list j expected by HorizonAnneal4, the structured form of HA
	j=list(d=dmat,TaxaName=as.character(vTaxa$GRCode),SectionName=as.character(mySectionNames$SecNum),PenaltySpec=PenaltySpec)
    
	# returned list contains many organizations of the same data
	myBack=list(j=j,dmat=dmat,vTaxa=vTaxa,sectionList=mySectionNames,startScore=score)
}

############################################################################################
# set up the R based HA input file, from a set of OAL input sections
	#Added cutoff parameter to allow user to say how many occurrences a taxa must have to be included in composite
	#keyTaxa should be a dataframe with at least two columns 1)GRCode and 2)Starting scores of each GRCode
		#the starting scores have to be the last column in the keyTaxa file

###########################################################################################
##################### DEFUNCT, keyTaxa not written well ###################################
###########################################################################################
HAFile4R_v3<-function(mySections,cutoff=1,keyTaxa="NA"){
	# create the HA Dmat and initial score matrix
	myTaxa=BGTaxaList(mySections)
	mySectionNames=BGSectionList(mySections)
	nsections=length(mySections)
	
	#note filtering here to taxa with more than a certain number of occurrence-saves time, but be a bit careful
	#vTaxa=myTaxa[myTaxa$Occurrences>cutoff,]
	vTaxa=myTaxa[myTaxa$Occurrences>cutoff | substring(myTaxa$GRCode,1,2)=="EB" | substring(myTaxa$GRCode,1,2)=="KB",]
	vTaxa=vTaxa[vTaxa$GRCode!="99999",]
	vTaxa=vTaxa[vTaxa$GRCode!="CO9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="CH9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="TR9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="EB9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="KB9999",]

	torder=order(vTaxa$GRCode)
	vTaxa=vTaxa[torder,]	
	#list of informative taxa
	ntaxa=length(vTaxa$GRCode)
	
	TaxaNumbers=1:ntaxa
	names(TaxaNumbers)=vTaxa$GRCode
	
	hcount=rep(0,nsections)						
	# set up counts of the number of horizons in each section
	for(i in 1:nsections){
		hcount[i]=mySections[[i]]@NumHorizons
	}
	hOffset=rep(0,nsections)					
	# set up the offsets
	for(i in 2:nsections){
		hOffset[i]=hOffset[i-1]+hcount[(i-1)]
	}
	# total horizon count
	nHorizonsTotal=sum(hcount)
	# set up Dmat structure for HA
	# column 1- initial score
	# column 2- section number
	# column 3- horizon number
	# column 4-should be horizon height-for the moment, simply using horizon number
	# column 5 to ntaxa+3-  taxa
	
	# empty data matrix---filled with (-1) values....
	dmat=matrix(rep(-1,(ntaxa+4)*nHorizonsTotal), nrow=nHorizonsTotal,ncol=(ntaxa+4))
    
	for(i in 1:nsections){
    		# fill in section number
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),2]=rep(mySections[[i]]@SectionNumber,hcount[i])
    		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),4]=mySections[[i]]@Height[1:mySections[[i]]@NumHorizons]
		#Output horizon heights as column 4 of dmat, added by JBoyle March 2nd 2016
    		for(j in 1:mySections[[i]]@NumTaxa){
    			# check to see if taxa is valid
    			valid=1
			# exclude 9999 entries
    			if(length(grep("9999",mySections[[i]]@GRCode[j]))>0){	
    				valid=0
    			}
			# check to see if taxa are on valid list
    			if(sum(vTaxa$GRCode==mySections[[i]]@GRCode[j])<1){	
    				valid=0
    			}
			# taxa are correct, place in matrix
    			if(valid){
    				dmat[((hOffset[i]+1):(hOffset[i]+hcount[i])),(TaxaNumbers[mySections[[i]]@GRCode[j]]+4)]=mySections[[i]]@Dmat[j, ]
    			}
		}
	}
    
	# now set up starting solution
	# column 1 is the score,  column 2 is the section, column 3 is the horizon, colum 4 is the height sorted on column 1
	score=matrix(0,nrow=nHorizonsTotal,ncol=4)
	for(i in 1:nsections){
		#If no keyTaxa list all sections are assigned scores randomly
		if(keyTaxa=="NA"){
    			temp=sort(runif(hcount[i]))
    			score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    			dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp
    		}
		if(keyTaxa!="NA"){
			#Number of biozone taxa in keyTaxa list
			nZones<-length(keyTaxa[,1])
			tempZones<-c(NULL)
			#Always have the last column of the keyTaxa list be the scores
			ZoneScoresPos<-length(keyTaxa)
			#Check whether any biozone taxa are in current section
			for(m in 1:nZones){
				#Finds whether a section has each biozone taxa in it
				zoneTaxaPres<-which(mySections[[i]]@GRCode==keyTaxa$GRCode[m])
				#Should always give a value of 0(absent) or 1(present)
				tempZones<-c(tempZones,length(zoneTaxaPres))
			}
			#If section does not contain any biozone taxa assign score randomly
			if(all(tempZones==0)){
				temp<-sort(runif(hcount[i]))
    				score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    				dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp			
			}
			if(any(tempZones==1)){
				#Identify which biozone taxa are the 
				#botBound<-min(which(tempZones>0))
				botScore<-keyTaxa[min(which(tempZones>0)),ZoneScoresPos]
				#topBound<-max(which(tempZones>0))
				topScore<-keyTaxa[max(which(tempZones>0)),ZoneScoresPos]
				#Check if scores for bounding biozones are equivalent
					#Can occurr either because there is only one biozone taxa in a section
					#or if two taxa have the same starting score
				if(botScore==topScore){
					#In cases of equivalency scores are expanded to 0.5 the standard deviation of the scores in the keyTaxa file up and down
					#Checks if min value would drop below 0
					tempMin<-botScore-0.5*sd(keyTaxa[,ZoneScoresPos])
					#Checks if max value would rise above 1
					tempMax<-topScore+0.5*sd(keyTaxa[,ZoneScoresPos])
					#Bounds min score to 0
					if(tempMin<0){
						tempMin<-0
					}
					#Bounds max score to 1
					if(tempMax>1){
						tempMax<-1
					}
					temp<-sort(runif(hcount[i],min=tempMin,max=tempMax))
					score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    					dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp				
				}
				else{
					temp<-sort(runif(hcount[i],min=botScore,max=topScore))
					score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    					dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp
				}
			}
		}
		score[(1+hOffset[i]):(hOffset[i]+hcount[i]),2]=mySections[[i]]@SectionNumber
    		score[(1+hOffset[i]):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
		#######Supposed to allow actual horizon heights to be included, currently not working#############
		#for(k in 1:hcount[i]){
		#	score[(1+hOffset[k]),4]=mySections[[i]]@Height[k]
		#}
	}
	#Gives index values of scores ordered from least to greatest
	sc_order=order(score[,1])
	#Reorderes scores into ascending order
	score=score[sc_order,]
    
	# adding column names to this matrix did not slow the calculation down
	# adding row names slowed the R version of HA by 10 to 12 percent
	# so do not use rownames
	addColumnNames=TRUE
	if(addColumnNames){
	    	temp=c("Score","Section","Horizon","Height",as.character(vTaxa$GRCode))
    		colnames(dmat)<-temp
		#temp=as.character(mySectionNames$SecName)
    		#rownames(dmat)<-temp[dmat[,2]]
	}
    
	# set up the default style penaltry structure, assuming all taxa carry a penalty, but
	# that there are no other types of data
	# the penalty structure may have to be altered separately to allow other designs
    
	n_biostrat=ntaxa
	biostrat=1:ntaxa
	# for biostrat data, or Taxa FADs,LADs, the biostrat variable is the numbers of the columns with taxa

	n_pmag=0
	pmag=1

	# pmag is a list of the column(s) with paleomagnetic signals, or really any binary data, NA values are not counted

	n_dates=0
	dates=matrix(c(64,0,65,0,1000,65,0,66,0,1000,66,0,67,0,1000),nrow=3,byrow=TRUE)

	# each row of the dates matrix is a set of data to be entered into the passing penalty
	# the first entry on each row is the column of the lower variale
	# second entry on a row is the data type 0- singular date,  1- FAD, 2-LAD
	# third and fourth entries on each row are the column and type of the second variable
	# fifth value on each row is the weight

	n_ashes=0
	ashes=matrix(c(68,100,69,100),nrow=2,byrow=TRUE)
	n_continuous=0
	continuous=matrix(c(70,5,71,5),nrow=2,byrow=TRUE)
	PenaltySpec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,ashes=ashes,n_continuous=n_continuous,continuous=continuous)
	#Clean dmat by deleting horizons with only 9999 taxa
	displacedLast<-length(dmat[1,])
	UseList<-c()
	for (i in 1:nHorizonsTotal){
		if(max(dmat[i,5:displacedLast])==1){
			UseList<-c(UseList,i)
		}
	}
	dmat<-dmat[UseList,]

	# set up the list j expected by HorizonAnneal4, the structured form of HA
	j=list(d=dmat,TaxaName=as.character(vTaxa$GRCode),SectionName=as.character(mySectionNames$SecNum),PenaltySpec=PenaltySpec)
    
	# returned list contains many organizations of the same data
	myBack=list(j=j,dmat=dmat,vTaxa=vTaxa,sectionList=mySectionNames,startScore=score)
}

#myHA_v3<-HAFile4R_v3(mySectsSyn,cutoff=4,keyTaxa=keyTaxaList)

##################################################################################################
##################################################################################################
# set up the R based HA input file, from a set of OAL input sections
	#Added cutoff parameter to allow user to say how many occurrences a taxa must have to be included in composite
	#keyTaxa should be a dataframe with at least two columns 1)GRCode and 2)Starting scores of each GRCode
		#the starting scores have to be the last column in the keyTaxa file
	#Updated keyTaxa scaffolding March 4th, 2022 (JBoyle)

HAFile4R_v4<-function(mySections,cutoff=1,keyTaxa="NA"){
	# create the HA Dmat and initial score matrix
	myTaxa=BGTaxaList(mySections)
	mySectionNames=BGSectionList(mySections)
	nsections=length(mySections)
	
	#note filtering here to taxa with more than a certain number of occurrence-saves time, but be a bit careful
	vTaxa=myTaxa[myTaxa$Occurrences>cutoff | substring(myTaxa$GRCode,1,2)=="EB" | substring(myTaxa$GRCode,1,2)=="KB",]
	vTaxa=vTaxa[vTaxa$GRCode!="99999",]
	vTaxa=vTaxa[vTaxa$GRCode!="CO9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="CH9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="TR9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="EB9999",]
	vTaxa=vTaxa[vTaxa$GRCode!="KB9999",]

	torder=order(vTaxa$GRCode)
	vTaxa=vTaxa[torder,]	
	#list of informative taxa
	ntaxa=length(vTaxa$GRCode)
	
	TaxaNumbers=1:ntaxa
	names(TaxaNumbers)=vTaxa$GRCode
	
	hcount=rep(0,nsections)						
	# set up counts of the number of horizons in each section
	for(i in 1:nsections){
		hcount[i]=mySections[[i]]@NumHorizons
	}
	hOffset=rep(0,nsections)
	if(nsections>1){					
		# set up the offsets
		for(i in 2:nsections){
			hOffset[i]=hOffset[i-1]+hcount[(i-1)]
		}
	}
	else{
		hOffset[1]=0
	}
	# total horizon count
	nHorizonsTotal=sum(hcount)
	# set up Dmat structure for HA
	# column 1- initial score
	# column 2- section number
	# column 3- horizon number
	# column 4-should be horizon height-for the moment, simply using horizon number
	# column 5 to ntaxa+3-  taxa
	
	# empty data matrix---filled with (-1) values....
	dmat=matrix(rep(-1,(ntaxa+4)*nHorizonsTotal), nrow=nHorizonsTotal,ncol=(ntaxa+4))
    
	for(i in 1:nsections){
    		# fill in section number
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),2]=rep(mySections[[i]]@SectionNumber,hcount[i])
    		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),3]=(1:hcount[i])
		dmat[(hOffset[i]+1):(hOffset[i]+hcount[i]),4]=mySections[[i]]@Height[1:mySections[[i]]@NumHorizons]
		#Output horizon heights as column 4 of dmat, added by JBoyle March 2nd 2016
    		for(j in 1:mySections[[i]]@NumTaxa){
    			# check to see if taxa is valid
    			valid=1
			# exclude 9999 entries
    			if(length(grep("9999",mySections[[i]]@GRCode[j]))>0){	
    				valid=0
    			}
			# check to see if taxa are on valid list
    			if(sum(vTaxa$GRCode==mySections[[i]]@GRCode[j])<1){	
    				valid=0
    			}
			# taxa are correct, place in matrix
    			if(valid){
    				dmat[((hOffset[i]+1):(hOffset[i]+hcount[i])),(TaxaNumbers[mySections[[i]]@GRCode[j]]+4)]=mySections[[i]]@Dmat[j, ]
    			}
		}
	}
    
	# now set up starting solution
	# column 1 is the score,  column 2 is the section, column 3 is the horizon, colum 4 is the height sorted on column 1
	score=matrix(0,nrow=nHorizonsTotal,ncol=4)
	for(i in 1:nsections){
		###############################################
		#print("Section Number")
		#print(mySections[[i]]@SectionNumber)
		###############################################
		#If no keyTaxa list all sections are assigned scores randomly
		if(keyTaxa=="NA"){
    			temp=sort(runif(hcount[i]))
    			score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    			dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp
    		}
		if(keyTaxa!="NA"){
			#Matrix with columns of KeyTaxa GRCode, Binary Local presence absence in a section, local FAD horizon index, KeyTaxaScorePosition
				#Always have the last column of the keyTaxa list be the scores
			ZonePres<-data.frame(matrix(NA,nrow=length(keyTaxa[,1]),ncol=4))
			ZonePres[,1]<-keyTaxa$GRCode
			ZonePres[,2]<-0
			ZonePres[,4]<-keyTaxa[,length(keyTaxa[1,])]
			nZones<-length(keyTaxa[,1])
			#Check whether any biozone taxa are in current section
			for(m in 1:nZones){
				#Finds whether a section has each biozone taxa in it
				if(any(mySections[[i]]@GRCode==keyTaxa$GRCode[m])){
					ZonePres[m,2]<-1
					#FAD index position in the local section Dmat of the keyTaxon of interest
					ZonePres[m,3]<-min(which(mySections[[i]]@Dmat[which(mySections[[i]]@GRCode==keyTaxa$GRCode[m]),]==1))
				}
			}
			#print(ZonePres)
			#If section does not contain any biozone taxa assign score randomly
			if(all(ZonePres[,2]==0)){
				temp<-sort(runif(hcount[i]))
    				score[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]= temp
    				dmat[(1+hOffset[i]):(hOffset[i]+hcount[i]),1]=temp			
			}
			#Stepping through sections by the number of zones they are split into by keyTaxa
			########## Fixed by JBoyle on May 18, 2022 ####################
			if(any(ZonePres[,2]==1)){
				#Number of unique zones (defined as being different assigned scores) in a section
				uniZones<-sort(unique(ZonePres[which(ZonePres[,2]==1),4]))
				#Check for reversals and ties in FAD orders relative to that expectd from keyTaxa list (i.e. local section FAD reversals)
				##Cut out key taxa/zones with a local FAD out of order
					#Logic is that we have a good constraint on the FAD orders and that reversals and ties are likely a case of an earlier taxon missing its early occurrences locally
				#Temporary table to find local zone FADs and catch ties and reversal of keytaxa zones
					#First column is the local FAD of each zone, second column is the input score for each zone, third column is the output score for each zone which differs from column 2 only where there are conflicts in zone FADs locally
				TempTable<-data.frame(matrix(NA,nrow=length(uniZones),ncol=3))
				TempTable[,2]<-uniZones
				TempTable[,3]<-uniZones
				for(m in uniZones){
					TempTable[which(round(TempTable[,2],6)==round(m,6)),1]<-min(ZonePres[which(round(ZonePres[,4],6)==round(m,6)),3],na.rm=TRUE)
				}
				#Check for ties and reversal between zones
				nUniZones<-length(uniZones)
				for(m in 1:nUniZones){
					for(n in nUniZones:m){
						RevTieCheck<-round(TempTable[n,1],6)-round(TempTable[m,1],6)
						if(RevTieCheck<1){
							TempTable[m,3]<-TempTable[n,3]
							break
						}
					}
				}
				uniZones<-sort(unique(ZonePres[which(ZonePres[,2]==1),4]))
				#print("Pre-screening unique zones")
				#print(uniZones)
				#Replaces scores of zones to reflect reversal and ties in ordering of the zones with a section
				for(p in uniZones){
					ZonePres[which(ZonePres[,4]==p),4]<-TempTable[which(TempTable[,2]==p),3]
				}
				uniZones<-sort(unique(ZonePres[which(ZonePres[,2]==1),4]))
				#print("Post-screening unique zones")
				#print(uniZones)
				#TopBound is always 1, but scores get sequentially overwritten, bottom bound of score generation shifts (assume start is 0)
				TopBound<-1
				BotBound<-0
				#Starting point of constrained zone bins has to be relative to start of section 1+hOffset[i] so need a variable to make adjustments as we move through zones
				LowerStartAdd<-0
				for(m in uniZones){
					#Lowest local horizon index for a key taxa zone
					breakHorizon<-min(ZonePres[which(round(ZonePres[,4],6)== round(m,6)),3],na.rm=TRUE)
					#print("Break Horizon")
					#print(breakHorizon)
					#Scores assigned to lower split
					nBelow<-length(which(mySections[[i]]@HorNum<breakHorizon))
					#print("nBelow")
					#print(nBelow)
					LowTemp<-(1+hOffset[i]+LowerStartAdd)
					HighTemp<-(1+hOffset[i]+nBelow-1)
					#print("Lower Split edges")
					#print(LowTemp)
					#print(HighTemp)
					TempDif<-HighTemp-LowTemp
					if(TempDif>1){
						temp<-sort(runif(HighTemp-LowTemp+1,min=BotBound,max=m))
						score[LowTemp:HighTemp,1] = temp
						dmat[LowTemp:HighTemp,1] = temp
					}
					score[(1+hOffset[i]+nBelow),1] = m
					dmat[(1+hOffset[i]+nBelow),1] = m
					#Scores assigned to upper split
					nAbove<-length(which(mySections[[i]]@HorNum>breakHorizon))
					#print("nAbove")
					#print(nAbove)
					if(nAbove==0){
						score[1+hOffset[i]+breakHorizon-1,1] = m
						dmat[1+hOffset[i]+breakHorizon-1,1] = m
					}
					else{
						temp<-sort(runif(nAbove,min=m,max=1))
						LowTemp<-(2+hOffset[i]+(nBelow-LowerStartAdd)+LowerStartAdd)
						HighTemp<-(hOffset[i]+hcount[i])
						#print("Upper Split edges")
						#print(LowTemp)
						#print(HighTemp)
						if(LowTemp==HighTemp){
							score[LowTemp,1]=temp
							dmat[LowTemp,1]=temp
						}
						else{
							score[LowTemp:HighTemp,1]= temp
   							dmat[LowTemp:HighTemp,1]=temp
						}
					}
					#New bottom bound for the next step up
					BotBound<-m
					#Adjustment to place the middle of the next bounded region
					LowerStartAdd<-nBelow+1
				}
			}
		}
		LowTemp<-(1+hOffset[i])
		HighTemp<-(hOffset[i]+hcount[i])
		score[LowTemp:HighTemp,2]=mySections[[i]]@SectionNumber
    		score[LowTemp:HighTemp,3]=(1:hcount[i])
		#######Supposed to allow actual horizon heights to be included, currently not working#############
		#for(k in 1:hcount[i]){
		#	score[(1+hOffset[k]),4]=mySections[[i]]@Height[k]
		#}
	}
	#Gives index values of scores ordered from least to greatest
	sc_order=order(score[,1])
	#Reorderes scores into ascending order
	score=score[sc_order,]
    
	# adding column names to this matrix did not slow the calculation down
	# adding row names slowed the R version of HA by 10 to 12 percent
	# so do not use rownames
	addColumnNames=TRUE
	if(addColumnNames){
	    	temp=c("Score","Section","Horizon","Height",as.character(vTaxa$GRCode))
    		colnames(dmat)<-temp
		#temp=as.character(mySectionNames$SecName)
    		#rownames(dmat)<-temp[dmat[,2]]
	}
    
	# set up the default style penaltry structure, assuming all taxa carry a penalty, but
	# that there are no other types of data
	# the penalty structure may have to be altered separately to allow other designs
    
	n_biostrat=ntaxa
	biostrat=1:ntaxa
	# for biostrat data, or Taxa FADs,LADs, the biostrat variable is the numbers of the columns with taxa

	n_pmag=0
	pmag=1

	# pmag is a list of the column(s) with paleomagnetic signals, or really any binary data, NA values are not counted

	n_dates=0
	dates=matrix(c(64,0,65,0,1000,65,0,66,0,1000,66,0,67,0,1000),nrow=3,byrow=TRUE)

	# each row of the dates matrix is a set of data to be entered into the passing penalty
	# the first entry on each row is the column of the lower variable
	# second entry on a row is the data type 0- singular date,  1- FAD, 2-LAD
	# third and fourth entries on each row are the column and type of the second variable
	# fifth value on each row is the weight
	
	#ashes maxtrix is meaningless if n_ashes is 0, same for continuous and pmag
	n_ashes=0
	ashes=matrix(c(68,100,69,100),nrow=2,byrow=TRUE)
	n_continuous=0
	continuous=matrix(c(70,5,71,5),nrow=2,byrow=TRUE)
	PenaltySpec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,ashes=ashes,n_continuous=n_continuous,continuous=continuous)
	#Clean dmat by deleting horizons with only 9999 taxa
	displacedLast<-length(dmat[1,])
	UseList<-c()
	for (i in 1:nHorizonsTotal){
		if(max(dmat[i,5:displacedLast])==1){
			UseList<-c(UseList,i)
		}
	}
	dmat<-dmat[UseList,]

	# set up the list j expected by HorizonAnneal4, the structured form of HA
	j=list(d=dmat,TaxaName=as.character(vTaxa$GRCode),SectionName=as.character(mySectionNames$SecNum),PenaltySpec=PenaltySpec)
    
	# returned list contains many organizations of the same data
	myBack=list(j=j,dmat=dmat,vTaxa=vTaxa,sectionList=mySectionNames,startScore=score)
}

#myHA_v4<-HAFile4R_v4(mySectsSyn,cutoff=4,keyTaxa=keyTaxa)


##################################################################################################
# function to zero out a given section in a data matrix
# used to set up for HA jackknifing

zeroSection=function(dmat,targetSection)
{
	# the first column of the data matrix is the section number
	# the second column of the data matrix is the horizon number within the section
	# the third column is the horizon height
	# columns 4 to ncolumn are the entries for species 1 to n
	
	ncold=dim(dmat)[2]
	print(ncold)
	print(targetSection)
	dmat[dmat[,1]==targetSection, 4:ncold] = -1
	redmat=dmat
	
}

###########################################################################################
#function to write a single formatted dmatrix form to disk, using the input dmatrix and filename

writeDmat=function(dmat,fname)
{
	nrowsx=dim(dmat)[1]
	print(nrowsx)
	ncolsx=dim(dmat)[2]
	print(ncolsx)
	cat(nrowsx,file=fname,sep="\n")
	cat((ncolsx-3),file=fname,append=TRUE,sep="\n")
	cat("Section","\t","Horizons",file=fname,append=TRUE)
	cat("",file=fname,append=TRUE,sep="\n")
	for(i in 1:nrowsx)
	{
		x=dmat[i,]
		cat(sprintf("%f",x),file=fname,sep="\t",fill=FALSE,append=TRUE)
		cat("",file=fname,append=TRUE,sep="\n")
	}
}

#################################################################################################
# function to input a dmatrix and generate a series of jackknifed dmatrices, removing one section at a time
	#need to cut off the first column (scores) of the dmat befor this will work
		#dmat=j$d[,2:nTaxa+4]
jackknifeDmat=function(dmat,rootfile)
{
	nSections=length(unique(dmat[,2]))
	for(i in 1:nSections)
	{
		redMat=zeroSection(dmat,i)
		fname=sprintf("%s_%i.txt",rootfile,i)
		writeDmat(redMat,fname)
	}
	z=0
}
