###this is a bit crude but it does numerical integration from first principles to get area under median survival curve


dd<-readRDS("PostSurv_Sex3Inf_zAll.rds")
names(dd)
subdd<-dd[dd$Infection=="Infected Cub"&dd$Sex=="Male",c("Median","t")]

#calulcate area under whole curve using trapezoid rule
ddAUC<-0.5*1*(subdd[1,1]+2*sum(subdd[2:50,1])+subdd[51,1])
#this is an approximation whose accuracy will depend on the number of intervals
#comes from the trapezoid rule from https://www.lehigh.edu/~ineng2/clipper/notes/NumIntegration.htm
#0.5  h  [ f(x1) + 2  f(x1 + h) + 2  f(x1 + 2h) + f(x1 + 3h) ]

#then calculate area under curve conditioned on survival to age 30
#can change conditioned age by changing subscripts
ddAUC.age30<-(0.5*1*(subdd[31,1]+2*sum(subdd[32:50,1])+subdd[51,1]))/subdd[31,1]+30

###this is from model that ignores MLH. You can get same expected lifespan beyond 40 for each of the badger categories.
###challenge is what to do when MLH is in the model.
###is it possible to get a survival curve for a GIVEN value of MLH? If so, we can
###calculate ddAUC.age40 for each of a sequence of MLH values