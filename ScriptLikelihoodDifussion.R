############
#setwd("/Users/alessandrohernandez/Documents/Scripts")
library(deSolve)
library(rootSolve)

### Let's read all the stuff
#slim <- read.csv("slim.csv",header=TRUE)
#View(slim)
args = commandArgs(trailingOnly=TRUE)
print (args)
SNPFrequencyFileNumber = args[1]
SNPFrequencyFileTxt = paste("SimPro_",SNPFrequencyFileNumber, ".txt", sep = "")
slim = read.table(SNPFrequencyFileTxt,header=TRUE,sep="\t")
#head(SNPFrequencyFile)
#head(Medians)
#nrow(SNPFrequencyFile)
#nrow(Medians)
## Now let's reorder the data a little bit
#SNPFrequencyFile <- SNPFrequencyFile[SNPFrequencyFile$V3 != "no",]
#colnames(SNPFrequencyFile) <- c('Chr','Pos','Culture','A1', 'A2','AlleleFrequency','DerivedAlleles','AlleleCount')
#nrow(SNPFrequencyFile)
#MergedDataset <- merge(Medians, SNPFrequencyFile,by='Culture')
#MergedDataset<-MergedDataset[MergedDataset$Date_mean<5000,]
#slim <- MergedDataset[order(MergedDataset$Date_mean),]
### Get the tentative origin of the allele

LastOcurrenceData <- 0
i=0
for (i in 1:nrow(slim)){
    if (slim$DerivedAlleles[i] > 0){
        LastOcurrenceData <- i
    }
}

AlleleOriginLat <- slim$Lat[LastOcurrenceData]
AlleleOriginLong <- slim$Long[LastOcurrenceData]
AlleleOriginAge <- slim$Date_mean[LastOcurrenceData]
print(AlleleOriginAge)
print(AlleleOriginLat)
print(AlleleOriginLong)
######### Change coordinates

AlleleXAxis <-slim$Long #round(((slim$Lat - min(slim$Lat))/(max(slim$Lat)- min(slim$Lat))) * 99 + 1)
AlleleYAxis <-slim$Lat  #round(((slim$Long - min(slim$Long))/(max(slim$Long)- min(slim$Long))) * 99 + 1)

slim <- cbind(slim,AlleleXAxis)
slim <- cbind(slim,AlleleYAxis)

AlleleOriginLatX <-10 #round(((AlleleOriginLat - min(slim$Lat))/(max(slim$Lat)- min(slim$Lat))) * 99 + 1)
AlleleOriginLongY <-1 #round(((AlleleOriginLong - min(slim$Long))/(max(slim$Long)- min(slim$Long))) * 99 + 1)


######### Check this reference to see how the equation is solved https://cran.r-project.org/web/packages/rootSolve/vignettes/rootSolve.pdf

diffusion2D <- function(t, conc, par) {
Conc <- matrix(nrow = n, ncol = n, data = conc) # vector to 2-D matrix
dConc <- Conc*(1-Conc)*(Conc*d+s*(1-2*Conc)) # consumption
#dConc <- -r*Conc*Conc
BND <- rep(1, n) # boundary concentration

# constant production in certain cells
# dConc[ii]<- dConc[ii]+ p
#diffusion in X-direction; boundaries=imposed concentration

Flux <- -Dx * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]),
rep(0, n) )/dx
dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
#diffusion in Y-direction

Flux <- -Dy * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]),
rep(0, n))/dy
dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
return(list(as.vector(dConc)))
}

# parameters
dy <- dx <- 1 # grid size
Dy <- Dx <- 0.01 # diffusion coeff, X- and Y-direction
n <- 100
d <- 0.0
s <- 0.0
pars<-c(n,d,s)
# 10 random cells where substance is produced at rate p
# ii <- trunc(cbind(runif(10)*n+1, runif(10)*n+1))

### LL calculations

DifussionValuesToCheck <- c(0,0.01,0.02,0.03,0.04,0.05,0.06)
Counter <- 0

TimesToTest <- c()
Conc0 <- matrix(nrow = n, ncol = n, 0.)
for (i in 1:nrow(slim)){
    Time = round((AlleleOriginAge - slim$Date_mean[i]))
    if (Time >= 0){
    TimesToTest <- c(TimesToTest, Time)
    if (Time == 0){
        Conc0[slim$Lat[i], slim$Long[i]] = round(slim$AlleleFrequency[i],3)
    }
    }
}

###SelectionValuesToCheck <- c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5)

LL<-c()
for (j in 1:11){
  Dy <- Dx <- DifussionValuesToCheck[j]
  ##d = SelectionValuesToCheck[j]
  ##s = d / 2
  LL<-c(LL,0)
  N = 500 ### Population size per deme
    print (i)

    print(system.time(
    ST3 <- ode.2D(y = Conc0, times = c(0,TimesToTest), func = diffusion2D, parms = pars,
    dimens = c(n, n), method = rkMethod("rk45ck"))
    ))

    for (i in 1:length(TimesToTest)){
        if ( (ST3[i+1,(slim$AlleleXAxis[i]-1)*100 + slim$AlleleYAxis[i]] > 0) && (slim$AlleleCount[i] > 0) && (ST3[i+1,(slim$AlleleXAxis[i]-1)*100 + slim$AlleleYAxis[i]] <= 1) ){
        # print (i, log (dbinom(slim$DerivedAlleles[i],slim$AlleleCount[i], ST3[i+1,(slim$AlleleXAxis[i]-1)*100 + slim$AlleleYAxis[i]])))
    LL[j] <- LL[j] + log (dbinom(slim$DerivedAlleles[i],slim$AlleleCount[i], ST3[i+1,(slim$AlleleXAxis[i]-1)*100 + slim$AlleleYAxis[i]]))
    }}
}
LLToPrint <- t(LL)
print(LLToPrint)

LLFile <- paste("Scores_Likelihood/LLSelVentajosos_Slim", SNPFrequencyFileNumber, ".txt", sep = "")
write.table(LLToPrint,file=LLFile,append=TRUE,row.names=FALSE,col.names=FALSE)
