# Version 4 of Horizon Annealing in R
# allows for a variable number of geological variables- chemostratigraphy, sequence stratigraphy, paleomagnetic
# Sections are chosen in proportion to their horizon contribution to the solution
# Added a new move type (vice) that picks a single horizon and then attempts to push it up or down until
	#it either increases the penalty or runs against an adjacent horizon in its own section
#Added a new piece to the output to keep track of the number of each move type tried and how many times it was successful

HorizonAnneal4_v3<-function(j_in,param=list(nouter=400,ninner=100,temperature=5,cooling=0.9),move=c(1:6)){
	# j_in is the j-structure used for handling the input data
	# j_in$d is the grand capture history at the core of the HA calculation
	
	# extract the dmat from the input j list and assign it to the variable d3 to be manipulated 
	# during the search
	
	# the input j list must have
	# j_in$d- which is the data matrix of presence/absence info, with initial scores for all horizons
	# j_in$TaxaName- an ordered list of the names of all taxa and geological variables in the dmatrix, ordered the
	# same way the columns of dmat are ordered.
	# j_in$SectionName- a list of section names, such that the ith item in the list 
	# corresponds to section number i
	# j_in$PenaltySpec  - this is the penalty structure for the data set as outlined in the file
	# penalty_spec_template.R
	
	d3=j_in$d
	pen_str=j_in$PenaltySpec
		
	#  simulated annnealing based on biostratigraphy+chemostratigraphy
	# d3[,1] - xpositions
	# d3[,2] -section codes
	# d3[,3]- horizon number
	# d3[,4]-horizon height
	# d3[,4+1:lastcolumn] - geological variables- of all types
	
	# param is the list of parameter values controlling the annealing
	#
	# pen_str is the penalty structure specifying how the penalty is to be calculated
	# penalty_spec=list(n_biostrat=n_biostrat,biostrat=biostrat,n_pmag=n_pmag,pmag=pmag,n_dates=n_dates,dates=dates,n_ashes=n_ashes,n_continuous=n_continuous,continuous=continuous)
	
	
	
	ptm <- proc.time()
	taxaOffset=4
	# assumes that the data has position scores in the first column, reorder the data to match	
  

	d3ord=order(d3[,1])
	d3a=d3[d3ord,]
	d3a[,1]=(d3a[,1]-min(d3a[,1]))/(max(d3a[,1])-min(d3a[,1]))
	gaprange=mean(diff(d3a[,1]))						# size of typical gap	
 
	movetrack=c(0,0,0,0,0,0)
	movetry=c(0,0,0,0,0,0)
	cat("Initial Penalty Calculation","\n")
 	# compute the biostratigraphic range extension  
	if(pen_str$n_biostrat>0){
		cv=NetRangeExtension(d3a[,taxaOffset+pen_str$biostrat])
	}
	else{
		cv=0
	}
	cat('biostrat error was ',cv,'\n')
	#compute the paleomagnetic or binary range extension
	if(pen_str$n_pmag==1){
  		cv=cv+ftransitions2(d3a[,pen_str$pmag+taxaOffset])*30    
	}	
	else{
		if(pen_str$n_pmag>1){
			cv=cv+sum(apply(d3a[,pen_str$pmag+taxaOffset],2,ftransitions2))*30
		}
	}
	# compute the AshRangeExtensionPenalty
	if(pen_str$n_ashes==1){
  		cv=cv+AshRangeExtension(d3a[,pen_str$ashes[,1]+taxaOffset])*10  
  		cat('value of penalty string on ashes',pen_str$ashes[,1]+taxaOffset,"\n")
  		cat("ash was ",AshRangeExtension(d3a[,pen_str$ashes[,1]+taxaOffset]), "\n" )  
	}
	else{
		#if(pen_str$n_ashes>1)
		#{
		#	  cv=cv+sum(apply(d3a[,pen_str$ashes[,1]+taxaOffset],2,AshRangeExtension))*10  
		#}
		#Added by JBoyle for individual ash penalties
		if(pen_str$n_ashes>1){
			for(p in 1:pen_str$n_ashes){
				cv<-cv+AshRangeExtension(d3a[,pen_str$ashes[p,1]+taxaOffset])*pen_str$ashes[p,2]
			}
		}
	}
	# compute the continuous variable penalty
	if(pen_str$n_continuous==1){
  		cv=cv+ftransitions2(d3a[,pen_str$continuous[,1]+taxaOffset])*pen_str$continuous[,2]   
	}
	else{
		if(pen_str$n_continuous>1){
			cat("in calculation of continuous, n>1")
			cv=cv+sum(apply(d3a[,pen_str$continuous[,1]+taxaOffset],2,ftransitions2))*pen_str$continuous[1,2]
			cat(cv)
		}
	}
	#compute the "Passing Penalty"
	if(pen_str$n_dates>0){
  		for(i in 1:pen_str$n_dates){
  		    cv=cv+PassingError2(d3a[,pen_str$dates[i,1]+taxaOffset],pen_str$dates[i,2],d3a[,pen_str$dates[i,3]+taxaOffset],pen_str$dates[i,4])*pen_str$dates[i,5]
  		    cat(i,pen_str$dates[i,1]+taxaOffset,pen_str$dates[i,2],pen_str$dates[i,3]+taxaOffset,pen_str$dates[i,4])
  		    cat(i,'Passing Error', PassingError2(d3a[,pen_str$dates[i,1]+taxaOffset],pen_str$dates[i,2],d3a[,pen_str$dates[i,3]+taxaOffset],pen_str$dates[i,4])*pen_str$dates[i,5],"\n")
  		}    
	}	   
	bestcv=cv
	cat("Starting penalty")
	print(bestcv)
	bestd3=d3a
	temperature=param$temperature
	#nsections=max(d3a[,2])
	nsections=length(unique(j_in$d[,2]))
	history=matrix(0,ncol=3,nrow=param$nouter)
	nhorizons=dim(d3a)[1]
	nMoves<-length(move)
  
	for(j in 1:param$nouter){
		for(i in 1:param$ninner){
			pd3=d3a
			# propose change
			#pInt<-sample.int(nsections,1)
			#psec=floor(runif(1,min=1.000001,max=nsections+0.99999))      # pick the section
			#psec<-sort(unique(j_in$d[,2]))[pInt]
			##Picks a section via picking a horizon
				#Sections will be picked proportional to the number of horizons
				#Possibly less time spent moving data-poor sections
			psec<-j_in$d[sample.int(nhorizons,1),2]
			pmove<-move[sample.int(nMoves,1)]
			#cat(j,i,"\n")
			#cat(psec,pmove,"\n")
			if(pmove==1){            #shift up or down						
				shval=runif(1,min=-0.1,max=0.1)
				ps=(pd3[,2]==psec)
				pd3[ps,1]=pd3[ps,1]+shval
				movetype=1
				movetry[1]=movetry[1]+1
			}
			if(pmove==2){                   #expand/contract
				shval=runif(1,min=-0.05,max=0.05)+1;
				ps=(pd3[,2]==psec)
				pd3min=mean(pd3[ps,1])
				pd3[ps,1]=(pd3[ps,1]-pd3min)*shval+pd3min;
				movetype=2
				movetry[2]=movetry[2]+1
			}
			#insert or remove gap
			if(pmove==3){ 					#insert or remove gap
				ps=(pd3[,2]==psec)
				while(sum(ps)<3){
					#pInt<-sample.int(nsections,1)
					#psec=floor(runif(1,min=1.000001,max=nsections+0.99999))      # pick the section
					#psec<-sort(unique(j_in$d[,2]))[pInt]
					##Picks a section via picking a horizon
						#Sections will be picked proportional to the number of horizons
						#Possibly less time spent moving data-poor sections
					psec<-j_in$d[sample.int(nhorizons,1),2]
					ps=(pd3[,2]==psec)
				}
				#cat("Section ",psec,"\n")
				breakpoint=floor(runif(1,min=2,max=sum(ps)-0.001))
				w=which(ps)	
				gap=runif(1,min=-0.59,max=5)*(pd3[w[breakpoint+1],1]-pd3[w[breakpoint],1])
				w=w[-(1:breakpoint)]
				pd3[w,1]=pd3[w,1]+gap
				movetype=3
				movetry[3]=movetry[3]+1
			}
			if(pmove==4){								#insert dogleg  (was at 0.6)
				shval=runif(1,min=-0.1,max=0.1)+1
				ps=(pd3[,2]==psec)
				while(sum(ps)<3){
					#pInt<-sample.int(nsections,1)
					#psec=floor(runif(1,min=1.000001,max=nsections+0.99999))      # pick the section
					#psec<-sort(unique(j_in$d[,2]))[pInt]
					##Picks a section via picking a horizon
						#Sections will be picked proportional to the number of horizons
						#Possibly less time spent moving data-poor sections
					psec<-j_in$d[sample.int(nhorizons,1),2]
					ps=(pd3[,2]==psec)
				}
				w=which(ps)
				breakpt=floor(runif(1,min=2,max=sum(ps)-0.001))
				gapval=diff(pd3[w,1])
				upchoice=runif(1,min=0,max=1)
				if(upchoice>0.5){
					gapval[(breakpt+1):sum(ps)-1]=gapval[(breakpt+1):sum(ps)-1]*shval
				}
				else{
					gapval[1:(breakpt)]=gapval[1:(breakpt)]*shval;
				}
				newval=cumsum(c(pd3[w[1],1],gapval))
				pd3[w,1]=newval
				movetype=4
				movetry[4]=movetry[4]+1
			}
			if(pmove==5){				# insert shuffling routine, point shift
				target=floor(runif(1,min=1.000001,max=nhorizons+.99))
				dmove=floor(runif(1,min=0.01,max=1.99))
				nmove=ceiling(abs(rnorm(1,0,4)))
				nmove=1
				movetype=5
				movetry[5]=movetry[5]+1
				if(dmove==0){					# shuffle up
					startsection=pd3[target,2]
					while(nmove>0){
						if(target==nhorizons){
							nmove=0
						}
						else{
							#cat("target ",target, pd3[target,2], pd3[target,1],"\n")
							#cat("target+1",target+1,pd3[target+1,2],pd3[target+1,1],"\n")
							if(pd3[target+1,2]==startsection){
								nmove=0
							}
							else{
								temp=pd3[target+1,1]
								pd3[target+1,1]=pd3[target,1]
								pd3[target,1]=temp
								target=target+1
								nmove=nmove-1
							}
						}     #end else on target
					}			#end while on nmove
				}			# end of if on dmove
				else{							#shuffle down
					while(nmove>0){
						if(target==1){
							nmove=0
						}
						else{
							if(pd3[target-1,2]==pd3[target,2]){
								target=target-1
							}
							else{
								temp=pd3[(target-1),1]
								pd3[(target-1),1]=pd3[target,1]
								pd3[(target),1]=temp
								target=target-1
								nmove=nmove-1
							}
						}     #end else on target
					}			#end while on nmove
				} # end of else based on dmove
			}
			if(pmove==6){	#Vice Move
				movetype=6
				movetry[movetype]=movetry[movetype]+1
				##Picks a section via picking a horizon
					#Sections will be picked proportional to the number of horizons
					#Possibly less time spent moving data-poor sections
				psec<-j_in$d[sample.int(nhorizons,1),2]
				hors<-which(pd3[,2]==psec)
				nhors<-length(hors)
				targetHor<-hors[sample.int(nhors,1)]
				#Determine if top or bottom vice
				if(targetHor==nhorizons){
					ViceType=2
				}
				if(targetHor==1){
					ViceType=1
				}
				if(targetHor!= 1 & targetHor!= nhorizons){
					ViceType<-sample.int(2,1)
				}
				#Bottom Vice
				if(ViceType==1){
					pd3<-BottomViceMove(pd3,targetHor)
				}
				#Top Vice
				if(ViceType==2){
					pd3<-TopViceMove(pd3,targetHor)
				}
			}
			# order the proposed solution
			pd3<-pd3[order(pd3[,1]),]
			# find the error for the proposed solution, starting with biostrat data
			if(pen_str$n_biostrat>0){
				pcv3=NetRangeExtension(pd3[,taxaOffset+pen_str$biostrat])
			}
			else{
				pcv3=0
			}
			#compute the paleomagnetic or binary range extension
			if(pen_str$n_pmag==1){
				pcv3=pcv3+ftransitions2(pd3[,pen_str$pmag+taxaOffset])*30    
			}	
			else{
				if(pen_str$n_pmag>1){
					pcv3=pcv3+sum(apply(pd3[,pen_str$pmag+taxaOffset],2,ftransitions2))*30
				}
			}
			# compute the AshRangeExtensionPenalty
			if(pen_str$n_ashes==1){
				pcv3=pcv3+AshRangeExtension(pd3[,pen_str$ashes[,1]+taxaOffset])*10    
			}
			else{
				if(pen_str$n_ashes>1){
					pcv3=pcv3+sum(apply(pd3[,pen_str$ashes[,1]+taxaOffset],2,AshRangeExtension))*10  
				}
			}
			# compute the continuous variable range extension    
			if(pen_str$n_continuous>0){
				for(i in pen_str$continuous[,1]){
					pcv3=pcv3+ftransitions2(pd3[,i+taxaOffset])*pen_str$continuous[1,2]
				}    
			}
			# compute the passing error penalty
			if(pen_str$n_dates>0){
				for(i in 1:pen_str$n_dates){
					pcv3=pcv3+PassingError2(pd3[,pen_str$dates[i,1]+taxaOffset],pen_str$dates[i,2],pd3[,pen_str$dates[i,3]+taxaOffset],pen_str$dates[i,4])*pen_str$dates[i,5]
				}    
			}
			# now decide whether to accept the proposed solution	
			if(pcv3<=bestcv){						# if the proposed is better than the best ever seen, accept it
				if(pcv3<bestcv){
					print(bestcv)
					movetrack[movetype]=movetrack[movetype]+1 }
				else{ 
					cat("swop","\t")
				}
					bestcv=pcv3
					bestd3=pd3
					d3a=pd3
					cv=pcv3
			}
			else{
				pch=runif(1)
				if(pch<exp(-(pcv3-cv)/temperature)){			# use the Boltzman factor to move up sometimes
					d3a=pd3
					cv=pcv3
				}
				else{
					# don't accept the change
				}
			}	
			# force rescale
			#d3a[,1]=(d3a[,1]-min(d3a[,1]))/(max(d3a[,1])-min(d3a[,1]))
			d3a[,1]=(1:length(d3a[,1]))/length(d3a[,1])
		}
		temperature=temperature*param$cooling
		history[j,1]=temperature
		history[j,2]=bestcv
		history[j,3]=pcv3
		cat("\nN outer: ",j, "T: ",temperature,"Best pen: ",bestcv,pmove,psec,"Recent prop pen: ", pcv3,"\n")
	}
	ptime=proc.time() - ptm
	
	# set up the output list,  outputting the current best grand capture history, also transfer the TaxaName, SectionName
	# and penalty spec to the output list,  by copying them from the input list
	# JBoyle edit, added variables that track the number of times each move type is tried and the number of times each move type improved the solution
	y=list(pen=bestcv,initpen=cv, d=bestd3,history=history,TaxaName=j_in$TaxaName, SectionName=j_in$SectionName,PenaltySpec=j_in$PenaltySpec,MoveTried=movetry,MoveImproved=movetrack)

	print(ptime)
	print(sprintf("Best Total Penalty: %f ", bestcv))
	if(pen_str$n_biostrat){  
		print(sprintf("Net taxa range extension %f",NetRangeExtension(bestd3[,2+pen_str$biostrat]))) 
	}
	pcv3=0
	if(pen_str$n_pmag>0){
		for(i in pen_str$pmag){
			pcv3=pcv3+ftransitions2(bestd3[,i+taxaOffset])*30
		}    
	}	
	print(sprintf("Pmag Penalty: %f",pcv3))
	print(sprintf("Pmag Reversals: %f",pcv3/30))
	acv=0
	if(pen_str$n_ashes>0){
		for(i in pen_str$ashes[,1]){
			acv=acv+AshRangeExtension(bestd3[,i+taxaOffset])*10
		}    
	}
	print(sprintf("Ash term: %i",acv))
	print(sprintf("Ash levels %i",acv/10))
	ccv=0
	if(pen_str$n_continuous>0){
		for(i in pen_str$continuous[,1]){
			ccv=ccv+ftransitions2(bestd3[,i+taxaOffset])
		}    
	}
	print(sprintf("Continuous variable term: %f ",ccv))
	pcv=0
	if(pen_str$n_dates>0){
		for(i in 1:pen_str$n_dates){
			pcv=pcv+PassingError2(bestd3[,pen_str$dates[i,1]+taxaOffset],pen_str$dates[i,2],bestd3[,pen_str$dates[i,3]+taxaOffset],pen_str$dates[i,4])*pen_str$dates[i,5]
		}    
	}
	print(sprintf("Passing Error: %i",pcv))
	movetrack=1000*movetrack/movetry
	print(sprintf("Move track values %f %f %f %f %f %f",movetrack[1],movetrack[2],movetrack[3],movetrack[4],movetrack[5],movetrack[6]))
	cat(gaprange)
	#PlotHAHistory(y)
	return(y)
}
