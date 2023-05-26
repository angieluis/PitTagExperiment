##################################################################
# PIT tag experiment 
##################################################################

library(tidyverse)
library(stringr)
library(MASS)
library(maptools)
library(lme4)
library(viridis)

all.data <- read.csv("PitData.csv")

Exp.data <- all.data %>% 
  group_by(Block,Building,Treatment) %>% 
  summarise(num.readings=n())

print(Exp.data,n=Inf)

exp.use <- which(Exp.data$num.readings>3) #reader failed if less than 3 readings, so remove


#####################################################################
# First, stats on the raw data for each experiment (a few days for 
# which a building had a particular treatment)

MNAs <- read.csv("MNAs.csv")

Exp.data$num.readingsM <- numeric(dim(Exp.data)[1]) # number of readings male
Exp.data$num.readingsF <- numeric(dim(Exp.data)[1]) # numebr of readings female
Exp.data$num.indiv <- numeric(dim(Exp.data)[1]) # number of unique individuals
Exp.data$num.males <- numeric(dim(Exp.data)[1]) # number of unique males
Exp.data$num.female <- numeric(dim(Exp.data)[1]) # number of unique females
Exp.data$MNA <- numeric(dim(Exp.data)[1]) # density index for the 2-week period 
Exp.data$num.days <- numeric(dim(Exp.data)[1]) # length of the experiment (time the building had that treatment)

for(i in 1:dim(Exp.data)[1]){
  dat <- all.data[,-18] %>%
    filter(Block==Exp.data$Block[i] & Treatment==Exp.data$Treatment[i] & Building==Exp.data$Building[i])
  
  sex.info <- dat %>% group_by(EarTag,Sex) %>% summarise(n())
  
  Exp.data$num.readingsM[i] <- ifelse(length(which(sex.info$Sex=="M"))>0, sex.info$`n()`[which(sex.info$Sex=="M")],0) 
  Exp.data$num.readingsF[i] <- ifelse(length(which(sex.info$Sex=="F"))>0, sex.info$`n()`[which(sex.info$Sex=="F")],0) 
  Exp.data$num.indiv[i] <- length(unique(dat$EarTag))
  Exp.data$num.males[i] <- length(which(sex.info$Sex=="M"))
  Exp.data$num.female[i] <- length(which(sex.info$Sex=="F"))
  Exp.data$MNA[i] <- MNAs$MNA[Exp.data$Block[i]]
  Exp.data$num.days[i] <- length(unique(dat$day))
}

print(Exp.data,n=Inf)


hist(Exp.data$num.readings/Exp.data$num.days)
hist(log(Exp.data$num.readings/Exp.data$num.days))

readings.lm1 <- lm(log(num.readings/num.days) ~ Treatment + MNA , data=Exp.data)
# this one is best by stepAIC where max was Treatment * MNA * Building

summary(readings.lm1)

betas <- readings.lm1$coefficients

control.color <- rgb(62,17,81,maxColorValue=255)
bedding.color <- rgb(70,142,139,maxColorValue=255)
food.color <- rgb(249,232,85,maxColorValue=255)

control.color2 <- rgb(62,17,81,maxColorValue=255,alpha=0.7*255)
bedding.color2 <- rgb(70,142,139,maxColorValue=255,alpha=0.7*255)
food.color2 <- rgb(249,232,85,maxColorValue=255,alpha=0.7*255)


plot(Exp.data$MNA,log(Exp.data$num.readings/Exp.data$num.days), type="n",ylab="log(readings per day)",xlab="MNA")

points(Exp.data$MNA[which(Exp.data$Treatment=="control")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="control")] ,col=control.color,pch=15)
points(Exp.data$MNA[which(Exp.data$Treatment=="bedding")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="bedding")] ,col=bedding.color,pch=17)
points(Exp.data$MNA[which(Exp.data$Treatment=="food")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="food")] ,col=food.color,pch=19)


curve( betas[1]+ betas[2] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=control.color,add=TRUE,lwd=2) 
curve( betas[1] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=bedding.color,add=TRUE,lwd=2) 
curve( betas[1]+ betas[3] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=food.color,add=TRUE,lwd=2) 


#################################################################
# Activity Patterns
#################################################################

# add sunrise and sunset times. need long,lat for site
lat <- -112.790639
long <- 46.036167

all.data$datetime <- mdy_hms(paste(all.data$date,all.data$Time, sep="-"),tz="America/Denver")

Sunset.time <- sunriset(crds=matrix(c(lat,long),nrow=1),if_else(hour(ymd_hms(all.data$datetime))>12,all.data$datetime,all.data$datetime-days(1)) # subtract a day if before noon (so comparing to time of last, not next, sunset)
,POSIXct.out=TRUE,direction="sunset")$time

Sunrise.time <- sunriset(crds=matrix(c(lat,long),nrow=1),if_else(hour(ymd_hms(all.data$datetime))<12,all.data$datetime,all.data$datetime+days(1)), ## add a day if before noon (so comparing to time of next, not last, sunrise)
  POSIXct.out=TRUE,direction="sunrise")$time

Since.sunset <- as.duration(interval(Sunset.time,all.data$datetime)) # in seconds unfortunately
To.sunrise <- as.duration(interval(Sunrise.time,all.data$datetime))

all.data$Sunset.time <- Sunset.time
all.data$Sunrise.time <- Sunrise.time
all.data$Since.sunset <- Since.sunset
all.data$To.sunrise <- To.sunrise

# in relation to interval sunset-sunrise

night.duration <- as.duration(interval(all.data$Sunset.time,all.data$Sunrise.time))

all.data$Prop.night.duration <- all.data$Since.sunset/night.duration

hist(all.data$Prop.night.duration,xlab="Night Duration",main="Histogram of readings")


###################################################################
### Now think about contacts 
###################################################################

source("PitTagfunctions.r")

all.intervals <- list() 
for(i in 1:dim(Exp.data)[1]){
  dat <- all.data %>%
    filter(Block==Exp.data$Block[i] & Treatment==Exp.data$Treatment[i] & Building==Exp.data$Building[i])
  data <- data.fun(dat,Exp.data$Building[i])	
  intervals <- interval.fun(data,Exp.data$Building[i])
  #add a columns for Block info
  for(l in 1:length(intervals)){
    intervals[[l]]$Block <- Exp.data$Block[i]
    intervals[[l]]$Treatment <- Exp.data$Treatment[i]
    intervals[[l]]$Building <- Exp.data$Building[i]
  }
  all.intervals <- append(all.intervals,intervals)
}

## Examine durations in buildings

all.durations <- all.intervals[[1]]$duration
all.dur.expinfo <- all.intervals[[1]][,4:6]
for(i in 2:length(all.intervals)){
  if(is.null(dim(all.intervals[[i]])[1])){next}
  if(dim(all.intervals[[i]])[1]>0){
    all.durations <- c(all.durations,all.intervals[[i]]$duration)
    all.dur.expinfo <- rbind(all.dur.expinfo,all.intervals[[i]][,4:6])
  }
}


# undoctored data:
hist(all.durations/60,breaks=500,xlab="minutes",main="") 
hist(all.durations/60,breaks=500,xlab="minutes",main="",xlim=c(0,200)) 

length(all.durations) #13731
length(which(all.durations>dhours(1)))/length(all.durations) # 4.1% are longer than 1 hour
length(which(all.durations>dhours(2)))/length(all.durations) # 2.5% are longer than 2 hours
length(which(all.durations>dhours(5)))/length(all.durations) # 1.5% are longer than 5 hours
length(which(all.durations<dminutes(1)))/length(all.durations) # 52.1% are shorter than 1 minute
length(which(all.durations<dminutes(5)))/length(all.durations) # 80.8% are shorter than 5 minutes
length(which(all.durations<dminutes(30)))/length(all.durations) # 93.8% are shorter than 30 minutes

boxplot(all.durations/60~all.dur.expinfo$Treatment,ylab="Durations in minutes",xlab="Treatment")
boxplot(log(all.durations/60)~all.dur.expinfo$Treatment,ylab="Durations in log(minutes)",xlab="Treatment")





############ now with corrected data (before removing all over 1 hour)

all.intervals.corrected <- list() 
for(i in 1:dim(Exp.data)[1]){
  dat <- all.data %>%
    filter(Block==Exp.data$Block[i] & Treatment==Exp.data$Treatment[i] & Building==Exp.data$Building[i])
  
  data <- data.fun.correct(dat,Exp.data$Building[i],max.time.int = 30) # try to correct if >30 min duration
  intervals <- interval.fun(data,Exp.data$Building[i])

  #add columns for Block info
  for(l in 1:length(intervals)){
    intervals[[l]]$Block <- Exp.data$Block[i]
    intervals[[l]]$Treatment <- Exp.data$Treatment[i]
    intervals[[l]]$Building <- Exp.data$Building[i]
  }
  
  all.intervals.corrected <- append(all.intervals.corrected,intervals)
}
# still quite a few that end on "in"
# look at duration of intervals
all.durations.corrected <- all.intervals.corrected[[1]]$duration
all.dur.cor.expinfo <- all.intervals.corrected[[1]][,4:6]
for(i in 2:length(all.intervals.corrected)){
  if(is.null(dim(all.intervals.corrected[[i]])[1])){next}
  if(dim(all.intervals.corrected[[i]])[1]==0){next}
  all.durations.corrected <- c(all.durations.corrected,all.intervals.corrected[[i]]$duration)
  all.dur.cor.expinfo <- rbind(all.dur.cor.expinfo,all.intervals.corrected[[i]][,4:6])
}

hist(all.durations.corrected/60,breaks=500,xlab="minutes",main="") 
hist(all.durations.corrected/60,breaks=500,xlab="minutes",main="",xlim=c(0,200)) 

length(all.durations.corrected) #13448
length(which(all.durations.corrected>dhours(1)))/length(all.durations.corrected) # 1.0% > 1hr
length(which(all.durations.corrected>dhours(2)))/length(all.durations.corrected) # .7% > 2hrs
length(which(all.durations.corrected>dhours(5)))/length(all.durations.corrected) # .3% > 5hrs
length(which(all.durations.corrected<dminutes(1)))/length(all.durations.corrected) # 55.0% < 1 minute
length(which(all.durations.corrected<dminutes(5)))/length(all.durations.corrected) # 84.8% < 5 min
length(which(all.durations.corrected<dminutes(30)))/length(all.durations.corrected) # 98.5% < 30 min

boxplot(all.durations.corrected/60~all.dur.cor.expinfo$Treatment,ylab="Durations in minutes",xlab="Treatment")
boxplot(log(all.durations.corrected/60)~all.dur.cor.expinfo$Treatment,ylab="Durations in log(minutes)",xlab="Treatment",main="Corrected Data")



####################################################################
#first put together one big contact network for each treatment
###################################################################


# There are 48 intervals longer than 1 hour -  I will remove them for now and see how things look

# contactmatrix.fun is defined in "PitTagfunctions.r"
# it creates a contact matrix based on input data
# takes arguments dataframe: this needs to be cut down to just the experiment you are interested in
# building, 
# weighted: if =TRUE, then number of contacts, if =FALSE, then just 0 if no contacts and 1 if any contacts
# can also include corrections as above if corrected=TRUE, then set max.time.int, and also optionally set 
# max.duration to remove intervals over this duration (in hours)


# a few of the experiment combinations had fewer than 3 recordings, meaning that not more than 1 animal was recorded as going in and then out. So remove those.

exp.use <- which(Exp.data$num.readings>3)

for(i in exp.use){
  dat <- all.data %>%
    filter(Block==Exp.data$Block[i] & Treatment==Exp.data$Treatment[i] & Building==Exp.data$Building[i])
  
  cont.mat <- contactmatrix.fun(dat,building=Exp.data$Building[i],weighted=FALSE, corrected=TRUE,max.time.int = 30,max.duration = 60)
  assign(paste("corrected.contact.matrix",Exp.data$Block[i],Exp.data$Building[i],Exp.data$Treatment[i],sep="."),cont.mat)
  
  
}

food.treatments <- intersect(which(Exp.data$Treatment=="food"),exp.use)
food.contact.list <- list()
for(i in 1:length(food.treatments)){
  nam <- 
    paste("corrected.contact.matrix",Exp.data$Block[food.treatments[i]],Exp.data$Building[food.treatments[i]],Exp.data$Treatment[food.treatments[i]],sep=".")
  
  food.contact.list[[i]] <- get(nam)
}

Food.contact.matrix <- contact.all.fun(food.contact.list, weighted=FALSE)


########################
bedding.treatments <- intersect(which(Exp.data$Treatment=="bedding"),exp.use)
bedding.contact.list <- list()
for(i in 1:length(bedding.treatments)){
  nam <- 
    paste("corrected.contact.matrix",Exp.data$Block[bedding.treatments[i]],Exp.data$Building[bedding.treatments[i]],Exp.data$Treatment[bedding.treatments[i]],sep=".")
  
  bedding.contact.list[[i]] <- get(nam)
}

Bedding.contact.matrix <- contact.all.fun(bedding.contact.list, weighted=FALSE)


########################
control.treatments <- intersect(which(Exp.data$Treatment=="control"),exp.use)
control.contact.list <- list()
for(i in 1:length(control.treatments)){
  nam <- 
    paste("corrected.contact.matrix",Exp.data$Block[control.treatments[i]],Exp.data$Building[control.treatments[i]],Exp.data$Treatment[control.treatments[i]],sep=".")
  
  control.contact.list[[i]] <- get(nam)
}

Control.contact.matrix <- contact.all.fun(control.contact.list, weighted=FALSE)



##############################################################

library(igraph)

food.graph <- graph.adjacency(Food.contact.matrix, mode="undirected",diag=FALSE)
sex <- character()
for(i in 1:length(V(food.graph)$name)){
  sex[i] <- unique(as.character(all.data$Sex[which(all.data$EarTag==V(food.graph)$name[i])]))
}
V(food.graph)$sex <- sex
V(food.graph)$color <- ifelse(V(food.graph)$sex=="M","steelblue","tomato")
plot(food.graph,vertex.label=NA,main="Food Treatment")

bedding.graph <- graph.adjacency(Bedding.contact.matrix, mode="undirected",diag=FALSE)
sex <- character()
for(i in 1:length(V(bedding.graph)$name)){
  sex[i] <- unique(as.character(all.data$Sex[which(all.data$EarTag==V(bedding.graph)$name[i])]))
}
V(bedding.graph)$sex <- sex
V(bedding.graph)$color <- ifelse(V(bedding.graph)$sex=="M","steelblue","tomato")
plot(bedding.graph,vertex.label=NA,main="Bedding Treatment")



control.graph <- graph.adjacency(Control.contact.matrix, mode="undirected",diag=FALSE)
sex <- character()
for(i in 1:length(V(control.graph)$name)){
  sex[i] <- unique(as.character(all.data$Sex[which(all.data$EarTag==V(control.graph)$name[i])]))
}
V(control.graph)$sex <- sex
V(control.graph)$color <- ifelse(V(control.graph)$sex=="M","steelblue","tomato")
plot(control.graph,vertex.label=NA,main="Control")



par(mfrow=c(1,3))
plot(control.graph,vertex.label=NA,main="Control",vertex.color=control.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")
plot(bedding.graph,vertex.label=NA,main="Bedding Treatment",vertex.color=bedding.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")
plot(food.graph,vertex.label=NA,main="Food Treatment",vertex.color=food.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")



Treatment.sum.table <- data.frame(Network=c("food","bedding","control"), 
    Nodes=rep(NA,3),Edges=rep(NA,3),Connectance=rep(NA,3), MeanDegree=rep(NA,3),
    MeanEvCent=rep(NA,3),MeanBetweenness=rep(NA,3),Transitivity=rep(NA,3),Modularity=rep(NA,3))
for(i in 1:3){
  ig <- get(paste(Treatment.sum.table$Network[i],"graph",sep="."))
  Treatment.sum.table$Nodes[i] <-length(V(ig)$name)
  Treatment.sum.table$Edges[i] <- length(E(ig))
  Treatment.sum.table$Connectance[i] <- length(E(ig))/(length(V(ig)$name)^2)
  Treatment.sum.table$MeanDegree[i] <- mean(degree(ig))
  Treatment.sum.table$MeanEvCent[i] <- mean(evc <- unlist(evcent(ig)$vector))
  Treatment.sum.table$MeanBetweenness[i] <- mean(betweenness(ig))
  Treatment.sum.table$Transitivity[i] <- transitivity(ig)
  com <- fastgreedy.community(ig)
  Treatment.sum.table$Modularity[i] <- modularity(ig,com$membership)
}


Treatment.sum.table

##################################################################
# Now stats on each Block/building/treatment combo separately

Exp.data$Nodes  <- numeric(dim(Exp.data)[1])
Exp.data$Edges  <- numeric(dim(Exp.data)[1])
Exp.data$Connectance <- numeric(dim(Exp.data)[1])
Exp.data$MeanDegree <- numeric(dim(Exp.data)[1])
Exp.data$MeanBetweenness  <- numeric(dim(Exp.data)[1])
Exp.data$MeanEVCent  <- numeric(dim(Exp.data)[1])


for(i in exp.use){
  
  it <- paste(Exp.data$Block[i],Exp.data$Building[i],Exp.data$Treatment[i],sep=".")
  mat <- get(paste("corrected.contact.matrix.",it,sep=""))
  ig <- graph.adjacency(mat, mode="undirected",diag=FALSE)
  
  d <- degree(ig)
  betw <- betweenness(ig)
  evc <- unlist(evcent(ig)$vector)
  
  Exp.data$Nodes[i] <-length(V(ig)$name)
  Exp.data$Edges[i] <- length(E(ig))
  Exp.data$Connectance[i] <- length(E(ig))/(length(V(ig)$name)^2)
  Exp.data$MeanDegree[i] <- mean(d)
  Exp.data$MeanBetweenness[i] <- mean(betw)
  Exp.data$MeanEVCent[i] <- mean(evc) #Eigenvector centrality, includes weights
  
  assign(paste("graph",it,sep="."),ig)
  
}

print(Exp.data,n=Inf)

trt.table <- Exp.data %>%
  group_by(Treatment) %>%
  summarise(num.networks=n(),M.Nodes=mean(Nodes),Sd.Nodes=sd(Nodes),M.Edges=mean(Edges), Sd.Edges=sd(Edges),M.MDegree=mean(MeanDegree),Sd.MDegree=sd(MeanDegree), M.Conn=mean(Connectance), Sd.Conn=sd(Connectance))
print(trt.table)


#############################################################################
# Examine if treatment affected number of edges, mean degree or connectance
#############################################################################

stepAIC(glm.nb(Edges~Treatment*MNA,data=Exp.data)) 
mod.edges<- glm.nb(Edges~Treatment + MNA,data=Exp.data) # best model 
summary(mod.edges) #  food and MNA significant 

stepAIC(lm(MeanDegree~Treatment*MNA,data=Exp.data))
meandegree.lm <- lm(MeanDegree~Treatment*MNA,data=Exp.data)
summary(meandegree.lm) 

stepAIC(lm(Connectance~Treatment*MNA,data=Exp.data))
connectance.lm <- lm(Connectance~Treatment*MNA,data=Exp.data)
summary(connectance.lm)

############################################################################
# 3- panel fig

pdf(file="regressionplots3panel.pdf",width=4,height=9,onefile=TRUE)
par(mfrow=c(3,1))

########### number of readings
betas <- readings.lm1$coefficients
plot(Exp.data$MNA,log(Exp.data$num.readings/Exp.data$num.days), type="n",ylab="log(readings per day)",xlab="Density of mice on trapping grid")
points(Exp.data$MNA[which(Exp.data$Treatment=="control")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="control")] ,col=control.color,pch=15)
points(Exp.data$MNA[which(Exp.data$Treatment=="bedding")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="bedding")] ,col=bedding.color,pch=17)
points(Exp.data$MNA[which(Exp.data$Treatment=="food")],log(Exp.data$num.readings/Exp.data$num.days)[which(Exp.data$Treatment=="food")] ,col=food.color,pch=19)
curve( betas[1]+ betas[2] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=control.color,add=TRUE,lwd=2) 
curve( betas[1] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=bedding.color,add=TRUE,lwd=2) 
curve( betas[1]+ betas[3] + betas[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=food.color,add=TRUE,lwd=2) 

########### MEan degree

betas.md <- meandegree.lm$coefficients
plot(Exp.data$MNA,Exp.data$MeanDegree, type="n",ylab="Mean degree",xlab="Density of mice on trapping grid")
points(Exp.data$MNA[which(Exp.data$Treatment=="control")],Exp.data$MeanDegree[which(Exp.data$Treatment=="control")] ,col=control.color,pch=15)
points(Exp.data$MNA[which(Exp.data$Treatment=="bedding")],Exp.data$MeanDegree[which(Exp.data$Treatment=="bedding")] ,col=bedding.color,pch=17)
points(Exp.data$MNA[which(Exp.data$Treatment=="food")],Exp.data$MeanDegree[which(Exp.data$Treatment=="food")] ,col=food.color,pch=19)
curve( betas.md[1]+ betas.md[2] + betas.md[4]*x + betas.md[5]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=control.color,add=TRUE,lwd=2) 
curve( betas.md[1] + betas.md[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=bedding.color,add=TRUE,lwd=2) 
curve( betas.md[1]+ betas.md[3] + betas.md[4]*x + betas.md[6]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=food.color,add=TRUE,lwd=2) 

### connectance
betas.c <- connectance.lm$coefficients
plot(Exp.data$MNA,Exp.data$Connectance, type="n",ylab="Connectance",xlab="Density of mice on trapping grid")
points(Exp.data$MNA[which(Exp.data$Treatment=="control")],Exp.data$Connectance[which(Exp.data$Treatment=="control")] ,col=control.color,pch=15)
points(Exp.data$MNA[which(Exp.data$Treatment=="bedding")],Exp.data$Connectance[which(Exp.data$Treatment=="bedding")] ,col=bedding.color,pch=17)
points(Exp.data$MNA[which(Exp.data$Treatment=="food")],Exp.data$Connectance[which(Exp.data$Treatment=="food")] ,col=food.color,pch=19)
curve( betas.c[1]+ betas.c[2] + betas.c[4]*x + betas.c[5]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=control.color,add=TRUE,lwd=2) 
curve( betas.c[1] + betas.c[4]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=bedding.color,add=TRUE,lwd=2) 
curve( betas.c[1]+ betas.c[3] + betas.c[4]*x + betas.c[6]*x ,from=min(Exp.data$MNA),to=max(Exp.data$MNA),col=food.color,add=TRUE,lwd=2)

dev.off()



#############################################################################
# now make one large network with all connections so can relate 
# degree to characteristics of individuals
#############################################################################

All.contact.list <- list()
all.matrices <- objects()[grep("corrected.contact.matrix.",objects())]
for(i in 1:length(all.matrices)){
  All.contact.list[[i]] <- get(all.matrices[i])
}

All.contact.matrix <- contact.all.fun(All.contact.list, weighted=FALSE)


all.graph <- graph.adjacency(All.contact.matrix, mode="undirected",diag=FALSE)
summary(all.graph)
sex <- character()
for(i in 1:length(V(all.graph)$name)){
  sex[i] <- unique(as.character(all.data$Sex[which(all.data$EarTag==V(all.graph)$name[i])]))
}
V(all.graph)$sex <- sex
V(all.graph)$color <- ifelse(V(all.graph)$sex=="M","steelblue","tomato")
plot(all.graph,vertex.label=NA,main="Full Network")

########################### Does it follow the 20-80 rule?
plot(0:37,degree.distribution(all.graph,cumulative=TRUE),ylab="Cumulative Degree Distribution",xlab="Degree")
# 81 individuals total. So 65th individual when ranked and higher would be in the top 20% by degree
sum(sort(individual.data$degree)[65:81]) / sum(individual.data$degree)
# top 20% have 51% of contacts


all.degree <- degree(all.graph)
all.betweenness <- betweenness(all.graph)
all.eigencentrality <- unlist(evcent(all.graph)$vector)

individual.data <- tibble(ID=V(all.graph)$name,degree=all.degree,betweenness=all.betweenness,EVCent=all.eigencentrality,sex=V(all.graph)$sex)

individual.data$num.readings <- numeric(dim(individual.data)[1])
for(i in 1:dim(individual.data)[1]){
  individual.data$num.readings[i] <- length(which(all.data$EarTag==individual.data$ID[i]))
}

boxplot(num.readings~sex,data=individual.data)
boxplot(log(num.readings)~sex,data=individual.data)

hist(individual.data$degree)
hist(log(individual.data$degree))
hist(individual.data$EVCent)
boxplot(degree~sex,data=individual.data)
boxplot(EVCent~sex,data=individual.data,ylab="eigenvector centrality")

mod.sex.degree <- lm(degree~sex, data=individual.data)
summary(mod.sex.degree) # sex is not significant predictor of degree (across treatments)
# but trend of males lower



individual.data %>%
  ggplot(aes(x=sex, y= degree,fill=sex)) +
  geom_violin()+
  scale_fill_viridis(discrete=TRUE)+
  theme_bw()+
  scale_y_continuous(name="Degree") +
  theme(
    legend.position="none", plot.title=element_text(size=11)) + 
  
  xlab("")





######################################################

pdf(file="TreatmentNetworks.pdf",width=6,height=3,onefile=TRUE)
par(mfrow=c(1,3))
plot(control.graph,vertex.label=NA,main="Control",vertex.color=control.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")
plot(bedding.graph,vertex.label=NA,main="Nesting material Treatment",vertex.color=bedding.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")
plot(food.graph,vertex.label=NA,main="Food Treatment",vertex.color=food.color2,vertex.frame.color="black",edge.width=2,edge.color="gray4")
dev.off()

####################################################

