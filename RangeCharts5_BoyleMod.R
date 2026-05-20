# HA Plotting routines,  based on the HA4 data format
# this is the BoyleMod version


PlotRangeChartQ2<-function(j,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,labelactive=TRUE,barOrder="None"){
	# this version of the Range chart plot uses the j structure
	# of HA version 4, so that j$d is the data matrix, j$TaxaName has the taxa names
	# j$SectionName has the names of all the sections
	# if multiple plotting- alter hbar, try hbar=0.4
	# note that this version assumes taxa start in column 5,  so score is column 1, section is column 2, horizon # is 3, horizon height is 4
	# setting barOrder to "Fad" will sort species by FAD
	# note that this routine will return a list of the FADs and LADs of all species
	# in the data set, as positions in the composite.
	
	
	j_ch=j$d	
	pt=dim(j_ch)
	print("size of data set")
	print(pt)
	
	# get taxa range
	toffset=4;				# offset for start of taxa, Horizon anneal 4 uses offset of 4
	trange=apply(j_ch[,species_list+toffset],2,ColumnFADLAD)
	#print(trange)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    
    y=species_list
    npts=length(species_list)
	y=cbind(y,trange)
	yoff=min(y[,2])
	#y[,2]=y[,2]-yoff
	#y[,3]=y[,3]-yoff
	
	if(barOrder=='Fad')
	{
		forder=order(y[,2])
		y=y[forder,]
	}
	
	
	cat(labelactive)
	if(doplot){
		if(overplot==FALSE)
		{
	#	quartz()
		plot(c(0,npts+1),c(-(pt[1])/5,(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		  for(k in 1:(npts)){
		  	 lines(c(k,k+hbar),c(y[k,2],y[k,2]),col="black")
	         lines(c(k,k+hbar),c(y[k,3],y[k,3]),col="black")
		     lines(c(k,k),c(y[k,2],y[k,3]))	
		     lines(c(k,k)+hbar,c(y[k,2],y[k,3])) 
		     if(labelactive)
		     {
		     	text(k+0.5,y[k,2],c(j$TaxaName[y[k,1]],"   "),adj=c(1,0),cex=0.6,srt=90,font=3) 
		     	cat(j$TaxaName[y[k,1]],"\n")
		     }
		  }
		}
		else
		{	
		  hoffset=0.45
		  
		  for(k in 1:(npts)){
		     lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,2],y[k,2]),col="black",lty=2)
	         lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,3],y[k,3]),col="black",lty=2)
		     lines(c(y[k,1],y[k,1])+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		     lines(c(y[k,1],y[k,1])+hbar+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		  }
	
		}
		
	}
	if(labelactive)
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=j$TaxaName[y[,1]])
	}
	else
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=rep(NA,length(species_list)))
	}
    return(y2)	
}

#############################################


PlotRangeChartQ3<-function(j,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,labelactive=TRUE,barOrder="None",myDict=myDict){
  ###Modified by JBoyle to allow actual names to be printed on chart due to changes in HA_cTools output
  ###Also sidesteps issue of 9999_GlobalSynonym_ prefix on taxa currently causing problems
  ###Requires the name of the dictionary variable, have to specify myDict=myDict etc or causes a recurse error
  # this version of the Range chart plot uses the j structure
  # of HA version 4, so that j$d is the data matrix, j$TaxaName has the taxa names
  # j$SectionName has the names of all the sections
  # if multiple plotting- alter hbar, try hbar=0.4
  # note that this version assumes taxa start in column 5,  so score is column 1, section is column 2, horizon # is 3, horizon height is 4
  # setting barOrder to "Fad" will sort species by FAD
  # note that this routine will return a list of the FADs and LADs of all species
  # in the data set, as positions in the composite.
  # this version plots the section numbers of all finds on the range chart as well
  
  
  j_ch=j$d	
  pt=dim(j_ch)
  print("size of data set")
  print(pt)
  
  # get taxa range
  toffset=4;				# offset for start of taxa, Horizon anneal 4 uses offset of 4
  trange=apply(j_ch[,species_list+toffset],2,ColumnFADLAD)
  #print(trange)
  trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
  
  y=species_list
  npts=length(species_list)
  y=cbind(y,trange)
  yoff=min(y[,2])
  #y[,2]=y[,2]-yoff
  #y[,3]=y[,3]-yoff
  
  if(barOrder=='Fad')
  {
    forder=order(y[,2])
    y=y[forder,]
  }
  #Added by JBoyle to allow taxa to be ordered by Lad
  if(barOrder=='Lad'){
    forder=order(y[,3],decreasing=TRUE)
    y<-y[forder,]
  }
  
  
  cat(labelactive)
  if(doplot){
    if(overplot==FALSE)
    {
      #	quartz()
      plot(c(0,npts+1),c(-(pt[1]),(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
      for(k in 1:npts){
        lines(c(k,k+hbar),c(y[k,2],y[k,2]),col="black")
        lines(c(k,k+hbar),c(y[k,3],y[k,3]),col="black")
        lines(c(k,k),c(y[k,2],y[k,3]))	
        lines(c(k,k)+hbar,c(y[k,2],y[k,3])) 
        if(labelactive)
        {
          #Line below assigns names that still have 9999_GlobalSynonym glitch
          #Gives names as GRCodes
          #text(k+0.5,y[k,2],c(j$TaxaName[y[k,1]],"   "),adj=c(1,0),cex=0.6,srt=90,font=3)
          ###Altered to guarantee adding taxa name, requires the myDict variable name as a parameter
          #Also requires using the modified function HAFile4R2() in HA_C Tools.R to get GRCodes
          TaxaNamePres<-myDict$TaxaName[which(myDict$GRCode==as.character(j$TaxaName[y[k,1]]))]			
          text(k+0.5,y[k,2],c(TaxaNamePres,"   "),adj=c(1,0),cex=1.0,srt=90,font=3)
          cat(j$TaxaName[y[k,1]],"\n")
        }
      }
    }
    else
    {	
      hoffset=0.45
      
      for(k in 1:(npts)){
        lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,2],y[k,2]),col="black",lty=2)
        lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,3],y[k,3]),col="black",lty=2)
        lines(c(y[k,1],y[k,1])+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
        lines(c(y[k,1],y[k,1])+hbar+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
      }
      
    }
    
  }
  if(labelactive)
  {
    y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],GRCode=j$TaxaName[y[,1]])
  }
  else
  {
    y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=rep(NA,length(species_list)))
  }
  
  # now go through the taxa, and plot the sections they were found in
  for(k in 1:npts)
  {
    myfinds=which(j_ch[,toffset+y[k,1]]==1)
    myfinds=myfinds[!is.na(myfinds)]
    mySecs=j_ch[myfinds,2]
    for(k2 in 1:length(myfinds))
    {
      mytext=sprintf("%i",mySecs[k2])
      text(k+0.5,myfinds[k2],mytext,adj=c(1,0),cex=0.5,col='red')
    }
  }
  
  
  return(y2)	
}





###########################################
PlotRangeChartQ4<-function(j,species_list,doplot=TRUE,overplot=FALSE,hbar=0.66,labelactive=TRUE,barOrder="None",myHA){
		###Modified by JBoyle to allow actual names to be printed on chart due to changes in HA_cTools output
		###Also sidesteps issue of 9999_GlobalSynonym_ prefix on taxa currently causing problems
		###Requires the name of the dictionary variable, have to specify myDict=myDict etc or causes a recurse error
	# this version of the Range chart plot uses the j structure
	# of HA version 4, so that j$d is the data matrix, j$TaxaName has the taxa names
	# j$SectionName has the names of all the sections
	# if multiple plotting- alter hbar, try hbar=0.4
	# note that this version assumes taxa start in column 5,  so score is column 1, section is column 2, horizon # is 3, horizon height is 4
	# setting barOrder to "Fad" will sort species by FAD
	# note that this routine will return a list of the FADs and LADs of all species
	# in the data set, as positions in the composite.
	# this version plots the section numbers of all finds on the range chart as well
	
	
	j_ch=j$d	
	pt=dim(j_ch)
	print("size of data set")
	print(pt)
	
	# get taxa range
	toffset=4;				# offset for start of taxa, Horizon anneal 4 uses offset of 4
	trange=apply(j_ch[,species_list+toffset],2,ColumnFADLAD)
	#print(trange)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
    
    y=species_list
    npts=length(species_list)
	y=cbind(y,trange)
	yoff=min(y[,2])
	#y[,2]=y[,2]-yoff
	#y[,3]=y[,3]-yoff
	
	if(barOrder=='Fad')
	{
		forder=order(y[,2])
		y=y[forder,]
	}
	#Added by JBoyle to allow taxa to be ordered by Lad
	if(barOrder=='Lad'){
		forder=order(y[,3],decreasing=TRUE)
		y<-y[forder,]
	}
	
	
	#find GRCodes of all taxa in use here-find GRCodes in the Dictionary
	# j$TaxaName[y[,1]] is the list of names plotted
	
	active_Taxa_Names=j$TaxaName[y[,1]]
	active_GR_Codes=rep("NA",length(active_Taxa_Names))
	for( i in 1:length(active_Taxa_Names))
	{
	  active_Taxa_Names[i]
	  active_GR_Codes[i]=as.character(myHA$vTaxa$GRCode[trimws(myHA$vTaxa$TaxaName)==trimws(active_Taxa_Names[i])])
	}
	cat(labelactive)
	if(doplot){
		if(overplot==FALSE)
		{
	#	quartz()
		plot(c(0,npts+1),c(-(pt[1]),(pt[1]+1)),type="n",xlab="Taxa Number",ylab="Position in Composite")
		  for(k in 1:npts){
		  	lines(c(k,k+hbar),c(y[k,2],y[k,2]),col="black")
			lines(c(k,k+hbar),c(y[k,3],y[k,3]),col="black")
			lines(c(k,k),c(y[k,2],y[k,3]))	
			lines(c(k,k)+hbar,c(y[k,2],y[k,3])) 
			if(labelactive)
			{
			#Line below assigns names that still have 9999_GlobalSynonym glitch
			#Gives names 
		     	text(k+0.5,y[k,2],c(j$TaxaName[y[k,1]],"   "),adj=c(1,0),cex=0.6,srt=90,font=3)
			###Altered to guarantee adding taxa name, requires the myDict variable name as a parameter
				#Also requires using the modified function HAFile4R2() in HA_C Tools.R to get GRCodes
			#TaxaNamePres<-myDict$TaxaName[which(myDict$GRCode==as.character(j$TaxaName[y[k,1]]))]			
			#text(k+0.5,y[k,2],c(TaxaNamePres,"   "),adj=c(1,0),cex=1.0,srt=90,font=3)
		   #  	cat(j$TaxaName[y[k,1]],"\n")
		     }
		  }
		}
		else
		{	
		  hoffset=0.45
		  
		  for(k in 1:(npts)){
		     lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,2],y[k,2]),col="black",lty=2)
	           lines(c(y[k,1],y[k,1]+hbar)+hoffset,c(y[k,3],y[k,3]),col="black",lty=2)
		     lines(c(y[k,1],y[k,1])+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		     lines(c(y[k,1],y[k,1])+hbar+hoffset,c(y[k,2],y[k,3]),col="black",lty=2)	
		  }
	
		}
		
	}
	if(labelactive)
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],GRCode=active_GR_Codes,TaxaName=j$TaxaName[y[,1]],FADscore=j_ch[y[,2],1], LADscore=j_ch[y[,3],1])
	}
	else
	{
		y2=data.frame(c9code=y[,1],FAD=y[,2],LAD=y[,3],SpeciesName=rep(NA,length(species_list)))
	}
	
	# now go through the taxa, and plot the sections they were found in
	for(k in 1:npts)
	{
		myfinds=which(j_ch[,toffset+y[k,1]]==1)
		myfinds=myfinds[!is.na(myfinds)]
		mySecs=j_ch[myfinds,2]
		for(k2 in 1:length(myfinds))
		{
			mytext=sprintf("%i",mySecs[k2])
			text(k+0.5,myfinds[k2],mytext,adj=c(1,0),cex=0.5,col='red')
		}
	}
	
	
    return(y2)	
}





PlotSectionRangeQ2<-function(y,ordinal=TRUE,PositionSupport=FALSE,overplot=FALSE,SectionList=mySectList){
	par(pch='_')
	# set up range
	# plot section range charts
   # this version shows the number of stratigraphically informative
    # taxa supporting the positioning of the horizon when PositionSupport=TRUE
    # set PositionSupport=TRUE to turn on the plotting of the number of informative
    # taxa on each horizon
    # set overplot=TRUE if you want to add another solution to this plot using
	
	if(overplot)
	{
	  hbar=0.25
	}
	else
	{
		hbar=0.5
	}  
	#nsections=max(y$d[,2])
	nsections<-length(unique(mySectList$SecNum))
	npts=dim(y$d)[1]
	sectmax=rep(0,nsections)
	sectmin=rep(1e9,nsections)

	lspace=0
    
	if(ordinal)
	{
		if(TRUE)			# relic code
		{
		  plot(c(0,max(y$d[,2])+1),c(-(dim(y$d)[1]/5),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		else
		{
		  plot(c(0,max(y$d[,2])+1),c(0,dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")	
		}  
		dpt=dim(y$d)
		dord=1:dpt[1]
		for(i in 1:npts)
		{
		  #Format of function below is draw line from x=c(x1,x2),y=c(y1,y2)
			#Creates bar of each horizons position on chart
		  lines(c(y$d[i,2],y$d[i,2]+hbar),c(dord[i],dord[i]))
		  stext=sprintf("%i",as.integer(sum(y$d[i,-(1:4)]>0,na.rm=TRUE)))
		  if(PositionSupport){
		    text(y$d[i,2]+hbar*1.3, dord[i],stext,cex=1,col="red")
		  }
		  print("Hello")
		  #Next two if statements to assign section numbers = index value on chart to section tops and bottoms
		  if(dord[i]<sectmin[y$d[i,2]]){ ##################Currently breaks on 14th loop of this command for Dar set
		     sectmin[y$d[i,2]]=dord[i]
		     }
		  if(dord[i]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=dord[i]
		     }
		}
		for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(TRUE)											# add section names
		    {
			text(i+0.5,sectmin[i],paste(as.character(SectionList$SecName[i])," -"),adj=c(1,0),cex=1,srt=90,font=3)
		      #text(i+0.5,sectmin[i],paste(j$SectionName[i]," -"),adj=c(1,0),cex=1,srt=90,font=3) 
		    }
			
		}
		
	}
	else{	
		 if(TRUE)
		 {
		   plot(c(0,max(y$d[,2])+1),c(-0.2,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 else
		 {
		   	plot(c(0,max(y$d[,2])+1),c(0,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		 }
		 for(i in 1:npts)
		 {
	      lines(c(y$d[i,2],y$d[i,2]+hbar),c(y$d[i,1],y$d[i,1]))
	      if(y$d[i,1]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=y$d[i,1]
		     }
		  if(y$d[i,1]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=y$d[i,1]
		     }
	    }
	    for(i in 1:nsections)
		{
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(TRUE)
		    {
		      #text(i+0.5,sectmin[i],paste(j$SectionName[i]," -"),adj=c(1,0),cex=1,srt=90,font=3) 
			text(i+0.5,sectmin[i],paste(as.character(SectionList$SecName[i])," -"),adj=c(1,0),cex=1,srt=90,font=3)
		    }
		}
	  } 
	  
	  
}

###########################################################################################################################################
PlotSectionRangeQ3<-function(y,ordinal=TRUE,PositionSupport=FALSE,overplot=FALSE,SectionList=mySectList){
	par(pch='_')
	# set up range
	# plot section range charts
	# this version shows the number of stratigraphically informative
	# taxa supporting the positioning of the horizon when PositionSupport=TRUE
	# set PositionSupport=TRUE to turn on the plotting of the number of informative taxa on each horizon
	# set overplot=TRUE if you want to add another solution to this plot using
	#Modified by JBoyle to plot output from HAFile4R2(), places all section names at the same location currently
		#Also, any sections that have no informative taxa (not in GCM) are in red
	
	if(overplot){
		hbar=0.25
	}
	else
	{
		hbar=0.5
	}  
	SecNums<-unique(j$SectionName)
	nsections<-length(SecNums)
	npts=dim(y$d)[1]
	sectmax=NULL
	sectmin=NULL
	lspace=0
    
	if(ordinal){
		if(TRUE){
			plot(c(0,nsections+1),c(-(dim(y$d)[1]/2),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		else{
			plot(c(0,nsections+1),c(0,dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")	
		}  
		dpt=dim(y$d)
		#Creates list of 1:number of horizons
		dord=1:dpt[1]
		for(i in 1:npts){
			#Format of function below is draw line from x=c(x1,x2),y=c(y1,y2)
			#Creates bar of each horizons position on chart
			lines(c(which(SecNums==y$d[i,2]),which(SecNums==y$d[i,2])+hbar),c(dord[i],dord[i]))
			stext=sprintf("%i",as.integer(sum(y$d[i,-(1:4)]>0,na.rm=TRUE)))
			if(PositionSupport){
				text(which(SecNums==y$d[i,2])+hbar*1.3,dord[i],stext,cex=1,col="red")
			}
		}
		for(i in 1:nsections){
			sectmin<-max(which(y$d[,2]==SecNums[i]))
			sectmax<-min(which(y$d[,2]==SecNums[i]))
			#Checks for case where a section is not part of the GCM because all taxa in a section occur less than the threshold value specified in earlier steps
			if(sectmin!="-Inf"){
				print(sectmin)
				lines(c(i,i),c(sectmin,sectmax))
				lines(c(i,i)+hbar,c(sectmin,sectmax))
				# add section names
				text(i+0.5,-10,paste(as.character(SectionList$SecName[which(SectionList$SecNum==SecNums[i])])," -"),adj=c(1,0),cex=1,srt=90,font=3)
			}
			#Adds section name in red if that section has no informative taxa
			if(sectmin=="-Inf"){
				text(i+0.5,-10,paste(as.character(SectionList$SecName[which(SectionList$SecNum==SecNums[i])])," -"),adj=c(1,0),cex=1,srt=90,font=3,col="red")
			}
		}
	}
	#Will crash if this is triggered
	else{	
		if(TRUE){
			plot(c(0,max(y$d[,2])+1),c(-0.2,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		else{
			plot(c(0,max(y$d[,2])+1),c(0,max(y$d[,1]+0.1)),type="n",xlab="Section Number",ylab="Position in Composite")
		}
		for(i in 1:npts){
			lines(c(y$d[i,2],y$d[i,2]+hbar),c(y$d[i,1],y$d[i,1]))
			if(y$d[i,1]<sectmin[y$d[i,2]]){
				sectmin[y$d[i,2]]=y$d[i,1]
			}
			if(y$d[i,1]>sectmax[y$d[i,2]]){
				sectmax[y$d[i,2]]=y$d[i,1]
			}
		}
		for(i in 1:nsections){
			lines(c(i,i),c(sectmin[i],sectmax[i]))
			lines(c(i,i)+hbar,c(sectmin[i],sectmax[i]))
			if(TRUE){
				#text(i+0.5,sectmin[i],paste(j$SectionName[i]," -"),adj=c(1,0),cex=1,srt=90,font=3) 
				text(i+0.5,sectmin[i],paste(as.character(SectionList$SecName[i])," -"),adj=c(1,0),cex=1,srt=90,font=3)
			}
		}
	}
}

###########################################################################################################################################
# adding a second section for comparison
PlotSectionRange2nd<-function(j,ordinal=TRUE){
	y=j
	par(pch='_')
	# set up range
	
	
	  hbar=0.4
	  hoffset=0.45
	
	 
	nsections=max(y$d[,2])
	npts=dim(y$d)[1]
	sectmax=rep(0,nsections)
	sectmin=rep(1e9,nsections)
	labelsActive=FALSE
	lspace=0
    
	if(ordinal)
	{
		 
		dpt=dim(y$d)
		dord=1:dpt[1]
		for(i in 1:npts)
		{
		  lines(c(y$d[i,2],y$d[i,2]+hbar)+hoffset,c(dord[i],dord[i]),col='blue')
		  if(dord[i]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=dord[i]
		     }
		  if(dord[i]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=dord[i]
		     }
		}
		for(i in 1:nsections)
		{
			lines(c(i,i)+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			lines(c(i,i)+hbar+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
		
		}
		
	}
	else{	
		 
		 for(i in 1:npts)
		 {
	      lines(c(y$d[i,2],y$d[i,2]+hbar)+hoffset,c(y$d[i,1],y$d[i,1]),col='blue')
	      if(y$d[i,1]<sectmin[y$d[i,2]]){
		     sectmin[y$d[i,2]]=y$d[i,1]
		     }
		  if(y$d[i,1]>sectmax[y$d[i,2]]){
		     sectmax[y$d[i,2]]=y$d[i,1]
		     }
	    }
	    for(i in 1:nsections)
		{
			lines(c(i,i)+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			lines(c(i,i)+hbar+hoffset,c(sectmin[i],sectmax[i]),col='blue',lty=4)
			
		}
	  } 
}

######################################################################################################
#Modified section plotting routine to look at distribution of a single taxa in a composite, requires using modified HAFile4R2 command
	#due to need for GRCodes in parameters, also requires output of BGSectionList() for paramter SectionList in CountingTools.R for SectionList parameter

PlotSingleTaxon<-function(y,SectionList=mySectList,target="GR1035",OccCol="blue"){
	#List of section numbers in order they occur in the composite
	SectNums<-unique(y$d[,2])
	#Empty plot to fill
	plot(c(0,length(SectNums)+1),c(-(dim(y$d)[1]*0.75),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite",main=target)
	nSects<-length(SectNums)
	#Width of the section bars and horizons
	hbar=0.5
	#Plot and label section bars, also mark top and bottom horizons
	for(i in 1:nSects){
		#Left side vertical
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i,y1=max(which(y$d[,2]==SectNums[i])))
		#Right side vertical
		segments(x0=i+hbar,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section FAD Bar
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=min(which(y$d[,2]==SectNums[i])))
		#Section LAD Bar
		segments(x0=i,y0=max(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section Name Labels
		text(i+0.5,min(which(y$d[,2]==SectNums[i])),paste(as.character(SectionList$SecName[which(SectionList$SecNum==SectNums[i])])," -"),adj=c(1,0),cex=0.75,srt=90,font=3)
	}
	#Get column number of target taxa in d matrix
	TargetLoc<-which(colnames(y$d)==target)
	#Get indices, height values, in composite where target taxa occurs from the d matrix
	TargetOccs<-which(y$d[,TargetLoc]==1)
	NOccs<-length(TargetOccs)
	#Plot horizons where target taxa occurs
	for(i in 1:NOccs){
		segments(x0=which(y$d[TargetOccs[i],2]==SectNums),y0=TargetOccs[i],x1=which(y$d[TargetOccs[i],2]==SectNums)+hbar,y1=TargetOccs[i],col=OccCol)
	}
}

#PlotSingleTaxon(j,SectionList=mySectList,target="GR1035",OccCol="red")

###################################################################################################
#Modified section plotting routine to look at distribution of a single taxa in a composite, requires using modified HAFile4R2 command
	#due to need for GRCodes in parameters, also requires output of BGSectionList() in CountingTools.R for SectionList parameter
	#Add color-coding to section names to indicate paleoplate affinity
	#PaleoContList is uplaoded csv file with paleoplate affinities of sections
PlotSingleTaxonCol<-function(y,SectionList=mySectList,target="GR1035",OccCol="blue",PaleoContList){
	#List of section numbers in order they occur in the composite
	SectNums<-unique(y$d[,2])
	#Empty plot to fill
	plot(c(0,length(SectNums)+1),c(-(dim(y$d)[1]*0.75),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite",main=target)
	nSects<-length(SectNums)
	#Width of the section bars and horizons
	hbar=0.5
	#Create Array of Paleocontinents
	PaleoList<-unique(PaleoContList[,1])
	#Plot and label section bars, also mark top and bottom horizons
	for(i in 1:nSects){
		#Left side vertical
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i,y1=max(which(y$d[,2]==SectNums[i])))
		#Right side vertical
		segments(x0=i+hbar,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section FAD Bar
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=min(which(y$d[,2]==SectNums[i])))
		#Section LAD Bar
		segments(x0=i,y0=max(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section Name Labels
		PaleoCol<-which(PaleoList==PaleoContList[which(PaleoContList[,2]==SectNums[i]),1])
		text(i+0.5,min(which(y$d[,2]==SectNums[i])),paste(as.character(SectionList$SecName[which(SectionList$SecNum==SectNums[i])])," -"),adj=c(1,0),cex=0.75,srt=90,font=3,col=PaleoCol)
		legend(90,-700,unique(PaleoContList[,1]),pch=c(rep(NA,6)),col=c(1:6),lty=c(rep(1,6)),cex=0.66)

	}
	#Get column number of target taxa in d matrix
	TargetLoc<-which(colnames(y$d)==target)
	#Get indices, height values, in composite where target taxa occurs from the d matrix
	TargetOccs<-which(y$d[,TargetLoc]==1)
	NOccs<-length(TargetOccs)
	#Plot horizons where target taxa occurs
	for(i in 1:NOccs){
		segments(x0=which(y$d[TargetOccs[i],2]==SectNums),y0=TargetOccs[i],x1=which(y$d[TargetOccs[i],2]==SectNums)+hbar,y1=TargetOccs[i],col=OccCol)
	}
}

#PlotSingleTaxonCol(j,SectionList=mySectsList,target="GR1012",OccCol="red",PaleoContList=SectPlate)

###################################################################################################
#Function to sequentially plot and output all taxa in a solution, height and width are for adjusting the pdf outputs
PlotAllSingleTaxa<-function(y,SectionList=mySectList,OccCol="blue",height=7,width=15){
	#List of taxa GRCodes in the d matrix
	TaxaList<-colnames(y$d)[-c(1:4)]
	nTaxa<-length(TaxaList)
	#Walk through plotting occurences of each taxa, output to pdf
	for(i in 1:nTaxa){
		PlotSingleTaxon(y=y,SectionList=SectionList,target=TaxaList[i],OccCol=OccCol)
		#Next three lines send each taxa pdf file to the working directory
		WorkingDir<-getwd()
		pathFile<-paste(WorkingDir,"/",TaxaList[i],".pdf",sep="")
		dev.print(pdf,pathFile,width,height=height,width=width)
	}
}

#PlotAllSingleTaxa(j,SectionList=mySectsList,OccCol="red",height=7,width=15)


####################################################################################################
#Modified section plotting routine equivalent to PlotSectionRangeQ2(), required if using modified HAFile4R2 command
	#due to need for GRCodes in parameters, also requires output of BGSectionList() in CountingTools.R for SectionList parameter

PlotAllHors<-function(y,SectionList=mySectList,OccCol="black"){
	#List of section numbers in order they occur in the composite, presorted for plotting
	SectNums<-unique(y$d[,2])
	#Empty plot to fill
	plot(c(0,length(SectNums)+1),c(-(dim(y$d)[1]*0.75),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
	nSects<-length(SectNums)
	#Width of the section bars and horizons
	hbar=0.5
	#Plot and label section bars, also mark top and bottom horizons
	for(i in 1:nSects){
		#Left side vertical
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i,y1=max(which(y$d[,2]==SectNums[i])))
		#Right side vertical
		segments(x0=i+hbar,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section Horizon Bars
		segments(x0=i,y0=c(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=c(which(y$d[,2]==SectNums[i])),col=OccCol)
		#Section Name Labels
		text(i+0.5,min(which(y$d[,2]==SectNums[i])),paste(as.character(SectionList$SecName[which(SectionList$SecNum==SectNums[i])])," -"),adj=c(1,0),cex=1,srt=90,font=3)
	}
}

#PlotAllHors(j,SectionList=mySectsList,OccCol="red")

####################################################################################################

TaxFadOrd<-function(y){
	TaxaList<-colnames(y$d)[-c(1:4)]
	nTaxa<-length(TaxaList)
	FadOrd<-c()
	FadPos<-c()
	for(i in 1:nTaxa){
		FadPos<-c(FadPos,min(which(y$d[,i+4]==1)))
	}
	FadPos<-cbind(FadPos,1:nTaxa)
	FadOrd<-order(FadPos[,1])
	FadPos<-FadPos[FadOrd,]
	return(FadPos[,2])
}

#FADsOrd<-TaxFadOrd(j)

####################################################################################################
#Modified taxa plotting routine equivalent to PlotRangeChartQ3()
	#Plots a taxa range chart with plain bars for presence instead of section number

PlotAllTaxa<-function(y,OccCol="black",myDict=myDict){
	TaxaList<-colnames(y$d)[-c(1:4)]
	nTaxa<-length(TaxaList)
	FadOrd<-TaxFadOrd(y=y)
	FadOrd<-FadOrd+4
	#Empty plot to fill
	plot(c(0,nTaxa+1),c(-(dim(y$d)[1]*0.75),dim(y$d)[1]+2),type="n",xlab="Taxon Index",ylab="Position in Composite",main="")
	hbar<-0.5
	for(i in 1:nTaxa){
		#Left side vertical
		segments(x0=i,y0=min(which(y$d[,FadOrd[i]]==1)),x1=i,y1=max(which(y$d[,FadOrd[i]]==1)))
		#Right side vertical
		segments(x0=i+hbar,y0=min(which(y$d[,FadOrd[i]]==1)),x1=i+hbar,y1=max(which(y$d[,FadOrd[i]]==1)))
		#Section Horizon Bars
		segments(x0=i,y0=c(which(y$d[,FadOrd[i]]==1)),x1=i+hbar,y1=c(which(y$d[,FadOrd[i]]==1)),col=OccCol)
		#Taxa Name Labels
		TaxaNamePres<-myDict$TaxaName[which(myDict$GRCode==as.character(TaxaList[FadOrd[i]-4]))]
		text(i+0.5,min(which(y$d[,FadOrd[i]]==1)),paste(TaxaNamePres," -"),adj=c(1,0),cex=1,srt=90,font=3)
	}
}

#PlotAllTaxa(j,OccCol="black",myDict=myDict)

####################################################################################################
#Modified section plotting routine equivalent to PlotSectionRangeQ2(), required if using modified HAFile4R2 command
	#due to need for GRCodes in parameters, also requires output of BGSectionList() in CountingTools.R for SectionList parameter
	#Add color-coding to section names to indicate paleoplate affinity
	#PaleoContList is uplaoded csv file with paleoplate affinities of sections

PlotAllHorsCol<-function(y,SectionList=mySectList,OccCol="black",PaleoContList){
	#List of section numbers in order they occur in the composite, presorted for plotting
	SectNums<-unique(y$d[,2])
	#Empty plot to fill
	plot(c(0,length(SectNums)+1),c(-(dim(y$d)[1]*0.75),dim(y$d)[1]+2),type="n",xlab="Section Number",ylab="Position in Composite")
	nSects<-length(SectNums)
	#Width of the section bars and horizons
	hbar=0.5
	#Create Array of Paleocontinents
	PaleoList<-unique(PaleoContList[,1])
	#Plot and label section bars, also mark top and bottom horizons
	for(i in 1:nSects){
		#Left side vertical
		segments(x0=i,y0=min(which(y$d[,2]==SectNums[i])),x1=i,y1=max(which(y$d[,2]==SectNums[i])))
		#Right side vertical
		segments(x0=i+hbar,y0=min(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=max(which(y$d[,2]==SectNums[i])))
		#Section Horizon Bars
		segments(x0=i,y0=c(which(y$d[,2]==SectNums[i])),x1=i+hbar,y1=c(which(y$d[,2]==SectNums[i])),col=OccCol)
		#Section Name Labels
		PaleoCol<-which(PaleoList==PaleoContList[which(PaleoContList[,2]==SectNums[i]),1])
		text(i+0.5,min(which(y$d[,2]==SectNums[i])),paste(as.character(SectionList$SecName[which(SectionList$SecNum==SectNums[i])])," -"),adj=c(1,0),cex=0.75,srt=90,font=3,col=PaleoCol)
		legend(90,-700,unique(PaleoContList[,1]),pch=c(rep(NA,6)),col=c(1:6),lty=c(rep(1,6)),cex=0.66)
	}
}

#PlotAllHorsCol(j,SectionList=mySectsList,OccCol="red",PaleoContList=SectPlate)

####################################################################################################
#Function to calculate a smoothed running diversity curve from a HA solution
SmoothedDiversity<-function(y,windowSize=10){
	TaxaList<-colnames(y$d)[-c(1:4)]
	nTaxa<-length(TaxaList)
	nHors<-length(y$d[,1])
	FadPos<-c()
	LadPos<-c()
	Diversity<-c()
	for(i in 1:nTaxa){
		FadPos<-c(FadPos,min(which(y$d[,i+4]==1)))
		LadPos<-c(LadPos,max(which(y$d[,i+4]==1)))
	}
	if(windowSize%%2==0){
		firstStep<-(windowSize/2)+1
		lastStep<-nHors-(windowSize/2)
	}
	else{
		firstStep<-round((windowSize/2),0)
		lastStep<-nHors-round((windowSize/2),0)
	}
	for(k in 1:nHors){
		if(k<firstStep){
			Diversity<-c(Diversity,NA)
		}
		else if(k>lastStep){
			Diversity<-c(Diversity,NA)
		}
		else{
			FadLog<-FadPos>k+5
			LadLog<-LadPos<k-5
			spp<-which(FadLog==FALSE & LadLog==FALSE)
			Diversity<-c(Diversity,length(spp))
		}
	}
	return(Diversity)
}

#DivSmooth<-SmoothedDiversity(j,windowSize=10)

####################################################################################################
#Function to calculate a smoothed running number of sections curve from a HA solution
SmoothedNSects<-function(y,windowSize=10){
	Sects<-unique(y$d[,2])
	nHors<-length(y$d[,1])
	SectFad<-c()
	SectLad<-c()
	nSects<-c()
	for(i in Sects){
		SectFad<-c(SectFad,min(which(j$d[,2]==i)))
		SectLad<-c(SectLad,max(which(j$d[,2]==i)))
	}
	if(windowSize%%2==0){
		firstStep<-(windowSize/2)+1
		lastStep<-nHors-(windowSize/2)
	}
	else{
		firstStep<-round((windowSize/2),0)
		lastStep<-nHors-round((windowSize/2),0)
	}
	for(k in 1:nHors){
		if(k<firstStep){
			nSects<-c(nSects,NA)
		}
		else if(k>lastStep){
			nSects<-c(nSects,NA)
		}
		else{
			#Checks which sections have FADs above window edge
			FadLog<-SectFad>k+5
			#Checks which sections have LADs below window edge
			LadLog<-SectLad<k-5
			spanSects<-which(FadLog==FALSE & LadLog==FALSE)
			nSects<-c(nSects,length(spanSects))
		}
	}
	return(nSects)
}

#nSectSmooth<-SmoothedNSects(j,windowSize=10)

####################################################################################################
#Function to calculate a smoothed running standard deviation by taxa curve from a HA solution
	#AllJack is the compiled jackknife results from the fetchAllJacks() function
SmoothedDuration<-function(y,windowSize=10){
	TaxaList<-colnames(y$d)[-c(1:4)]
	nTaxa<-length(TaxaList)
	nHors<-length(y$d[,1])
	FadPos<-c()
	LadPos<-c()
	Duration<-c()
	DurSmooth<-c()
	SectNums<-unique(y$d[,2])
	nSects<-length(SectNums)
	#Compiles standard deviation jackknife values for each horizon
	jackSD<-apply(AllJack[,4:nSects+3],1,sd,na.rm=TRUE)
	for(i in 1:nTaxa){
		TaxPos<-which(y$d[,i+4]==1)
		FadPos<-c(FadPos,min(TaxPos))
		LadPos<-c(LadPos,max(TaxPos))
		Duration<-c(Duration,LadPos-FadPos)
	}
	if(windowSize%%2==0){
		firstStep<-(windowSize/2)+1
		lastStep<-nHors-(windowSize/2)
	}
	else{
		firstStep<-round((windowSize/2),0)
		lastStep<-nHors-round((windowSize/2),0)
	}
	for(k in 1:nHors){
		if(k<firstStep){
			DurSmooth<-c(DurSmooth,NA)
		}
		else if(k>lastStep){
			DurSmooth<-c(DurSmooth,NA)
		}
		else{
			FadLog<-FadPos>k+5
			LadLog<-LadPos<k-5
			spp<-which(FadLog==FALSE & LadLog==FALSE)
			DurSmooth<-c(DurSmooth,mean(Duration[spp]))
		}
	}
	return(DurSmooth)
}

#SmoothedDurs<-SmoothedDuration(j,windowSize=10)

####################################################################################################
#Function to look at correlation of horizons with a single target horizon
	#j is the j structure output from HorizonAnneal4
	#jack parameter is the jacknife output of j, currently if it is read in from a csv file
	#targetHor is the target horizon, using the original position in the composite = first column of jack
	#myDictionary is the Dictionary file, needed for modified version of PlotRangeChartQ3, partially to sidestep 9999Global_Synonym glitch
SingleHorPackages<-function(j,jack,targetHor=1,resDiff,resCor,myDict,TaxaList){
	resDim<-dim(jack)
	#Call the PlotRangeChartQ3 function to give background image to write on top of
	PlotRangeChartQ3(j,4:resDim[2],barOrder="Fad",myDict=myDict)
	#Next block of text copied from PlotRangeChartQ3 to get FAD-ordered taxa matrix
	j_ch=j$d
	pt=dim(j_ch)
	toffset=4
	trange=applot(j_ch[,species_list+toffset],2,ColumnFADLAD)
	trange=matrix(unlist(trange),ncol=2,byrow=TRUE)
	y=TaxaList
	npts=length(species_list)
	y=cbind(y,trange)
	forder=order(y[,2])
	y=y[forder,]

	#Plotting correlation with color gradation from red (perfect correlation) to blue (perfect negative correlation)
	library(grDevices)
	colfunc<-colorRampPalette(c("Red","Blue"))
	GCM_dim<-dim(j$d)
	#Scans through GCM matrix to replace horizons by color gradation corresponding to correlation with target horizon
	for(k in 1:npts){
		HorTaxaList<-colnames(j$d(which(j$d[,4:resDim[2]-3]==1)))
		NHorTaxa<-length(HorTaxaList)
		for(m in 1:NHorTaxa){
			lines(x0=which(y[,1]==HorTaxaList[m]),y0=resDiff[k,1],x1=which(y[,1]==HorTaxaList[m])+hbar,y1=resDiff[k,1],col=colfunc(11)[round(resCor[k,which(resCor[,1]==targetHor)],0)])
		}
	}
	return(resCor)
}
