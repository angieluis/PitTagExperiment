###############################################################
# function to add columns to dataframe for whether in the building 
# or out, and a new time column which increases each day
###############################################################

library(schoolmath)


data.fun <- function(data, building, correct=FALSE){
  data <- as.data.frame(data) # in case problems with tibble format
  data <- data[order(data$datetime),]
  data <- data[data$Building==building,]
  #mark each entry as either in the building or out
  indiv <- sort(unique(data$EarTag[data$Building==building])) #all the animals in that building
  l <- length(indiv)
  in.or.out <- rep(NA,dim(data)[1])
  for(i in 1:l){
    in.or.out[which(data$EarTag==indiv[i])] <- rep(c("in","out"))
  }
  data$in.or.out <- in.or.out
  return(data)
}


############################################################################
## function to add columns of in or out, but with 'correction'...
############################################################################
# The detector threw out all readings that were within 15s
# of each other. So the animal could have come and gone in that time.
# That makes all the next intervals off. The time spent
# in the building at any one time is short. When I look at the 
# 'uncorrected' data, only 5% are longer than 1 hour, 54% are shorter 
# than 1 minute, 92% are shorter than 30 minutes.
# Here, you can set the max.time interval you think they had. 
# The default is 30 minutes. When a time intervals is longer
# than 30 minutes, this tries removing the 'in' and seeing if the 
# later intervals are better (less than 30 minutes). 
# Also can remove intervals over a set limit (in minutes)
############################################################################

data.fun.correct <- function(data, 
                             building, 
                             max.time.int = 30, # in minutes
                             max.duration=NULL){ # optional: remove intervals over this duration (in minutes)
  
  data <- data[order(data$EarTag, data$datetime),]
  data <- data[which(data$Building==building),]
  indiv <- sort(unique(data$EarTag[data$Building==building])) #all the animals in that building
  l <- length(indiv)

  new.dat <- data[-(1:dim(data)[1]),]
  new.dat$in.or.out <- character()
  
  for(i in 1:l){
    ind.dat <- data[which(data$EarTag==indiv[i]),]
    if(dim(ind.dat)[1]==1){next} 
     # flag any intervals longer than 30 minutes and see if removing the entrance time leads to shorter future intervals.
    dif.time <- diff(ind.dat$datetime) # time intervals in seconds
    flag <- which(dif.time>max.time.int*60 & is.odd(1:length(dif.time))) # which ones are longer than max.time minutes and are odd, because want only "in" to "out", not also "out" to "in".
      # this refers to the ind.dif$datetime that would be the start of the interval
      # (the line you would potentially remove)
    init.num.flags <- length(flag)
    counter <- 1
    
    # try to remove offending "in"s until the intervals look better (no more offending intervals) or the counter has exceeded the initial number of flags*2 - isn't going to do better
    while(length(flag)>0 & counter<(init.num.flags*2)){
      
      num.offenders <- length(flag) # how many offenders are there?
      ind.dat.try <- ind.dat[-flag[1],] #remove the "in" to the first offending interval
      new.dif.time <- diff(ind.dat.try$datetime) 
      
      #new.dif.time.odd <- new.dif.time[seq(1,length(new.dif.time),by=2)] #[flag[1]-1]# new next interval
      
      ## this [flag[1]-1] isn't right - removing too much
      
      if(new.dif.time[flag][1]<max.time.int*60 | flag[1]==length(dif.time) ){ # did removing that one improve just the next interval? Or is it the last time step?
        ind.dat <- ind.dat.try #If so, keep it
      }
      # if its the last reading that is removed, then this causes an error
      
      dif.time <- diff(ind.dat$datetime) # recalculate time intervals
      flag <- which(dif.time>max.time.int*60 & is.odd(1:length(dif.time)))  
      counter <- counter + 1 
    } # while loop   
        
   ### Add the "in" or "out" column to the data
     if(is.odd(dim(ind.dat)[1])){
       ind.dat$in.or.out <- c(rep(c("in","out"),dim(ind.dat)[1]/2),"in")
       warning(paste("Last reading is 'in', for Exp ",data$Experiment[1],", Bldg", building ,", individual ",indiv[i],sep=""))
     } else{
      ind.dat$in.or.out <- rep(c("in","out"), dim(ind.dat)[1]/2)
     }
    
    ## if removing rows of data with intervals over a max duration
    if(length(max.duration)>0){
      # calculate durations
      durs <- as.duration(interval(ind.dat$datetime[which(is.odd(1:length(ind.dat$datetime)))],ind.dat$datetime[which(is.even(1:length(ind.dat$datetime)))])) 
      rmi <- which(durs>dminutes(max.duration)) # this is interval not row
      if(length(rmi)>0){
        rmr <- (rmi[1]*2 - 1):(rmi[1]*2) # calculate what rows to remove
        if(length(rmi)>1){
          for( i in 2:length(rmi)){
            rmr <- c(rmr, (rmi[i]*2 - 1):(rmi[i]*2))
          }
        }
        ind.dat <- ind.dat[-rmr ,]
      }
      
    }
     
    new.dat <- rbind(new.dat,ind.dat)
  } # i for individuals
    
  return(new.dat)
}


############################################################################
#function to make a list of intervals the animal was in the building 
############################################################################

interval.fun=function(data,
                      building){
                      # max.duration=NULL # optional: remove intervals over a max duration
                      # works better to do this above in data.fun.corrected()

  data=data[data$Building==building,]
  indiv=sort(unique(data$EarTag))
  no.indiv=length(indiv)
  intervals=list(
  matrix(NA,nrow=length(which(data$in.or.out[data$EarTag==indiv[1]]=="out")),ncol=2))
  rm.ind <- numeric() #individuals to be removed (because <2 readings)
  for (j in 1:no.indiv){
    
    int=matrix(NA,nrow=length(which(data$in.or.out[data$EarTag==indiv[j]]=="out")),ncol=2)

    if(dim(int)[1]==0){
      rm.ind<- c(rm.ind,j)
      next}
    
    for(i in which(data$in.or.out[data$EarTag==indiv[j]]=="out")){

      int[i/2,]=as.character(data$datetime[data$EarTag==indiv[j]][(i-1):i])
	
    }
    
    int.dur <- data.frame(In=int[,1],Out=int[,2],duration=as.duration(interval(int[,1],int[,2])))
  
    intervals[[j]] <- int.dur
    intervals[[j]]$In <- as.character(intervals[[j]]$In)
    intervals[[j]]$Out <- as.character(intervals[[j]]$Out)
  
  }
  
  if(length(rm.ind)>0){
    names(intervals) <- sort(unique(data$EarTag))[-rm.ind]
  } else{
    names(intervals) <- sort(unique(data$EarTag))
  }
return(intervals)
}


############################################################################
#function to make contact matrix
############################################################################

contact.fun=function(intervallist,data,building, weighted=FALSE){
  indiv <- sort(unique(data$EarTag[data$Building==building]))
  
  contactmat=matrix(NA,length(intervallist),length(intervallist),dimnames=list(indiv,indiv))

  for(r in 1:length(intervallist)){
    for(c in 1:length(intervallist)){

      b=matrix(NA,ncol=length(intervallist[[c]][,1]),nrow=length(intervallist[[r]][,1]))

      for(i in 1:length(intervallist[[r]][,1])){	
        for(j in 1:length(intervallist[[c]][,1])){

b[i,j]=ifelse(((intervallist[[r]][i,1:2]<intervallist[[c]][j,2])[1]==1|(intervallist[[r]][i,1:2]<intervallist[[c]][j,2])[2]==1) 
&&((intervallist[[r]][i,1:2]>intervallist[[c]][j,1])[1]==1|(intervallist[[r]][i,1:2]>intervallist[[c]][j,1])[2]==1)
,1,0)

        } #i 	
      } # j

      contactmat[r,c] <- ifelse(any(b==1)==1,1,0)
      
      if(weighted==TRUE){
        contactmat[r,c] <- sum(b)
      }
    } # c
  } # r
  return(contactmat)
}




############################################################################
#function that puts it all together- just input data and get output of contact matrix
############################################################################

contactmatrix.fun=function(dataframe,building, 
                           weighted=FALSE, # if TRUE, count up number of contacts between each pair
                           corrected=FALSE,
                           max.time.int=NULL, # if corrected=TRUE, max time interval allowed
                           max.duration=NULL # optional:remove intervals over this duration (in hours)
                           ){
	
  if(corrected==TRUE){
    data <- data.fun.correct(dataframe,building,max.time.int=max.time.int,max.duration=max.duration)	
  } else{
    data <- data.fun(dataframe,building)
  }
  intervals=interval.fun(data,building)
  # get an error if an animal was only recorded once (going in but not coming out)
  # remove those individuals
  l <- lapply(intervals,length)
  if(length(which(l==0))>0){
    data <- data[-which(data$EarTag==names(intervals)[which(l==0)]),]
    intervals[[which(l==0)]] <- NULL
  }
  contactmat=contact.fun(intervals,data,building,weighted)	

  return(contactmat)
}
	
	
	

	





############################################################################
#put together the matrices from the different buildings/treatments
############################################################################



contact.all.fun=function(contact.list, # = list(contact.matrix.A,contact.matrix.B)
                         weighted=FALSE){ #add up contacts (weighted=TRUE) or just yes or no (weighted=FALSE)?
	
  indiv=character()
  for(i in 1:length(contact.list)){
    indiv <- sort(unique(c(indiv,dimnames(contact.list[[i]])[[1]])))
  }
  contact.matrix.all=matrix(NA,nrow=length(indiv), ncol=length(indiv),dimnames=list(indiv,indiv))
	
  for(r in 1:length(indiv)){
    for(c in 1:length(indiv)){
      indivr=dimnames(contact.matrix.all)[[1]][r]
      indivc=dimnames(contact.matrix.all)[[2]][c]

      contacts <- unlist(lapply(contact.list,function(x){x[which(rownames(x)==indivr),which(colnames(x)==indivc)]}))
      
      if(length(contacts)==0){
        contact.matrix.all[r,c] <- 0
      } else{
        contact.matrix.all[r,c] <- ifelse(sum(contacts)>0,1,0)
        if(weighted==TRUE){
          contact.matrix.all[r,c] <- sum(contacts)
        }
      }
    }
  }
      
      
  return(contact.matrix.all)
}	



############################################################

interval.fun2=function(data){

  indiv=sort(unique(data$EarTag))
  no.indiv=length(indiv)


  intervals=list(
    matrix(NA,nrow=length(which(data$in.or.out[data$EarTag==indiv[1]]=="out")),ncol=2))


  for (j in 1:no.indiv){
	
    int=matrix(NA,nrow=length(which(data$in.or.out[data$EarTag==indiv[j]]=="out")),ncol=2)

    for(i in which(data$in.or.out[data$EarTag==indiv[j]]=="out")){

      int[i/2,]=as.character(data$datetime[data$EarTag==indiv[j]][(i-1):i])
	
		}
    intervals[[j]]=int

  }
  return(intervals)
}


