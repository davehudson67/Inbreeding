
#Code for Dave Hudson

require(adegenet)

micro<-read.csv("S:/Clare B/PhD/Data/June 2015 Data Files/Data/Inbreeding Chapter Sept16onwards/Feb17_MarkedIndividuals_MicroDrop_MicroSats.csv")
names(micro)




dat93<-micro
names(dat93)

dat93<- dat93[, -1]          #neccessary formatting 
names(dat93)

temp<-dat93
names(temp)
OUT = NULL                      #pasting diploid columns together to make a single column 
for(i in seq(1, ncol(temp), 2)) { 
  loc  <- paste(temp[, i], temp[, i+1], sep="") 
  loc2 <- as.numeric(as.character(loc)) 
  OUT  <- cbind(OUT, loc2) 
} 








dat<-data.frame(OUT)
rownames(dat)<-micro$X
col.names<-c("L01","L02","L03","L04","L05","L06","L07","L08","L09","L10","L11","L12","L13","L14","L15","L16","L17","L18","L19","L20","L21","L22")
colnames(dat)<-col.names




?df2genind
require(adegenet)

microloc<-df2genind(dat,ncode=6,missing=NA,ploidy=2,pop=NULL)
names(microloc)
microloc@ind.names
microloc@loc.names



str.names<-microloc@ind.names
micro.names<-as.character(micro$X)



library(inbreedR)

#Remove Badger Names
names(micro)
gen2<-micro[,-1]
rownames(gen2)<-micro[,1]
head(gen2)

#Convert to 1/0 format with 1 column per locus
gen.convert<-inbreedR::convert_raw(gen2)
head(gen.convert)
colnames(gen.convert)<-col.names

loci.names<-unique(microloc@loc.names)

######### Sims - this code was provided by Xav Harrison

#Set Number of Simulations - each run does many permutations for inbreeding so higher numbers take disproportionately longer
nsims<-1
fhet.corr<-numeric(nsims)

for (k in 1:nsims){   #for each permutation
  
  loci1<-sample(loci.names,11)    #randomly select half the loci of the 22
  loci2<-loci.names[which(!loci.names %in% loci1)]	#which loci are left 
  
  ##### Inbreeding Calcs
  
  #Subset the genind object to just the loci in 'loci1'
  gen.subset<-microloc[,loc=loci1]
  
  #Calculate Inbreeding
  gen.subset.inbreed<-inbreeding(gen.subset,N=100)
  
  
  #Calculate Mean Inbreeding Per Individual
  inbreed.means<-sapply(gen.subset.inbreed,mean)
  
  ######### MLH Calcs
  
  gen.convert.subset<-gen.convert[,loci2]
  hets<-inbreedR::MLH(gen.convert.subset)
  
  
  ######### Merge Datasets by ID and calculate correlation
  
  hetdata<-data.frame(id=names(hets),het=hets)
  inbreed.data<-data.frame(id=names(inbreed.means),f= inbreed.means)
  fhet<-merge(hetdata,inbreed.data)
  fhet2<-na.omit(fhet) #added this line in as NA's generated
  
  ########## Store Correlation
  fhet.corr[k]<-with(fhet2,cor(het,f),)
  
} #loop end
hist(fhet$het)

#Function to calculate summary stats
fhet.stats<-function(corrdata){
  return(c(mean=mean(corrdata),CI=quantile(corrdata,c(0.025,0.975))))
}



fhet.stats(fhet.corr)

adegenet::hw.test(microloc)
hwtest<-adegenet::HWE.test.genind(microloc,res="matrix",permut=TRUE)
hwtestdf<-as.data.frame(hwtest)
bon<-bon<-0.05/22

idx<-which(hwtestdf<bon,TRUE) 
toto<-summary(microloc)
names(toto)



bartlett.test(list(toto$Hexp,toto$Hobs))