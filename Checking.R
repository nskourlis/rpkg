
library("rpkg")
library("survival")
library("mstate")
library("tidyverse")
library("flexsurv")
library("msm")
library("mstate")

load("data/ebmt.rda", envir = parent.frame(), verbose = FALSE)
load("data/cav.rda", envir = parent.frame(), verbose = FALSE)



head(ebmt)
head(cav)


tmat <- transMat(x = list(c(2, 3),c(3), c() ), names = c("Transplant", "Platelet Recovery", "Relapse/Death" ) )


ebmt$age2=  recode(ebmt$age, ">40" =0, "20-40"=1,"<=20" =0 )
ebmt$age3=  recode(ebmt$age, ">40" =1, "20-40"=0,"<=20" =0 )

msebmt <- msprep(data = ebmt, trans = tmat, 
                 time = c(NA, "prtime", "rfstime"), status = c(NA, "prstat", "rfsstat"), keep=c("age2","age3"))

head(msebmt)



results3_days=rpkg::msboxes_R(data=msebmt,id= msebmt$id, yb=c(0.3,0.5,0.75),
                              xb=c(0.5,0.2,0.7),boxwidth=0.1,boxheight=0.1,
                              tmat.= tmat, tstop=msebmt$Tstop,scale=365.25,
                              jsonpath="data", name="msboxes_EBMT_R.json" ) 



results3_days



summary(cav)

cav$id=cav$PTNUM
cav$group=c(rbinom(nrow(cav),1,0.5))


tmat=matrix(NA,nrow=4,ncol=4)
tmat[1,2]=1
tmat[1,4]=2
tmat[2,1]=3
tmat[2,3]=4
tmat[2,4]=5
tmat[3,2]=6
tmat[3,4]=7


results3_days=msboxes_R(data=cav,id= cav$id, yb=c(0.3,0.5,0.6,0.75), msm=TRUE,
                        xb=c(0.5,0.2,0.7,0.3),boxwidth=0.1,boxheight=0.1,
                        tmat.= tmat, vartime=seq(0,10,by=1),scale=1,
                        jsonpath="C:/Users/niksko/Desktop/mstate3/datasets/json/msm/json_present_msm", name="msboxes_cav_R.json" ) 


#########################################################################################################################################





## Multi-state model analysis: Using flexsurv_json function together with flexsurv package

### Provide time vector

tgrid <- seq(1, 2, by = 1)   

### Provide transition matrix

tmat <- rbind(c(NA, 1, 2), c(NA, NA, 3), c(NA, NA, NA)) 


### Run transition specific hazard models: Clock forward approach and use of flexible parametric models


cfwei.list<-vector(3,mode="list")

for (i in 1:3) {
  
  cfwei.list[[i]]<-flexsurvreg(Surv(Tstart,Tstop,status)~age2+age3,subset=(trans==i),
                               dist="weibull",data=msebmt)
}





### Prediction for different covariate patterns (the 3 age categories)
wh1 <- which(msebmt$age2 == 0 & msebmt$age3 == 0)
pat1 <- msebmt[rep(wh1[1], 3), 9:10]
attr(pat1, "trans") <- tmat


wh2 <- which(msebmt$age2 == 1 & msebmt$age3 == 0)
pat2 <- msebmt[rep(wh2[1], 3), 9:10]
attr(pat2, "trans") <- tmat

wh3 <- which(msebmt$age2 == 0 & msebmt$age3 == 1)
pat3 <- msebmt[rep(wh3[1], 3), 9:10]
attr(pat3, "trans") <- tmat



results_cf <- rpkg::flexsurv_json( model=cfwei.list, vartime=seq(365.25,365.25,by=365.25), 
                                   qmat=tmat, process="Markov",
                                   totlos=TRUE, ci.json=FALSE, cl.json=0.95, B.json=10, tcovs=NULL,
                                   Mjson=20, variance=FALSE,
                                   covariates_list=list(pat1,pat2,pat3), 
                                   jsonpath="~",
                                   name="predictions_EBMT_flex.json" ) 
results_cf



pmatrix.fs(x=cfwei.list, trans=tmat, t =1, newdata=list(),
          B = 10, ci = "TRUE", cl = 0.95)

test <- function() {
  log("not a number")
  print("R does stop due to an error and never executes this line")
}

test()     # throws an error

#######################################################################################################


library("rpkg")
library("survival")
library("mstate")
library("tidyverse")
library("flexsurv")
library("msm")
library("mstate")

load("data/ebmt.rda", envir = parent.frame(), verbose = FALSE)
load("data/cav.rda", envir = parent.frame(), verbose = FALSE)




options(scipen = 999,"digits"=10)

head(cav)


### Renaming variable PTNUM to id

cav$id=cav$PTNUM

### Defining the transition matrix

tmat=matrix(NA,nrow=4,ncol=4)
tmat[1,2]=1; tmat[1,4]=2; tmat[2,1]=3; tmat[2,3]=4
tmat[2,4]=5; tmat[3,2]=6; tmat[3,4]=7

### Defining the transition matrix with initial values under an initial assumption
Q<- rbind(c(0,0.25,0,0.25),c(0.166,0,0.166,0.166),c(0,0.25,0,0.25),c(0,0,0,0))

### Getting initial Q matrix in a default way- Feed the hand made matrix 
q.crude<- crudeinits.msm(state~years, id,data=cav, qmatrix=Q)

### Apply the msm model
cavsex.msm<- msm(state~years, covariates=~1, id,data=cav,qmatrix=q.crude, deathexact = 4, control=list(trace=1,REPORT=1)) 
summary(cavsex.msm)

### Prediction for different covariate patterns (males and females)

results <- rpkg::msmjson(msm.model=cavsex.msm, vartime=seq(1,1,1), mat.init=q.crude,
                         totlos=TRUE, visit=TRUE, sojourn=TRUE, pnext=TRUE, efpt=TRUE, envisits=TRUE,
                         ci.json="normal", cl.json=0.95, B.json=10,
                         cores.json=NULL,piecewise.times.json=NULL, piecewise.covariates.json=NULL,num.integ.json=FALSE,
                         jsonpath="data",
                         name="predictions_cav_R.json" ) 

  


#################################################################################################
## Multi-state model analysis: Using semipar_mstate_json function together with mstate package

library("rpkg")
library("survival")
library("mstate")
library("tidyverse")
library("flexsurv")
library("msm")
library("mstate")

load("data/ebmt.rda", envir = parent.frame(), verbose = FALSE)
load("data/cav.rda", envir = parent.frame(), verbose = FALSE)



head(ebmt)
head(cav)


tmat <- transMat(x = list(c(2, 3),c(3), c() ), names = c("Transplant", "Platelet Recovery", "Relapse/Death" ) )


ebmt$age2=  recode(ebmt$age, ">40" =0, "20-40"=1,"<=20" =0 )
ebmt$age3=  recode(ebmt$age, ">40" =1, "20-40"=0,"<=20" =0 )

msebmt <- msprep(data = ebmt, trans = tmat, 
                 time = c(NA, "prtime", "rfstime"), status = c(NA, "prstat", "rfsstat"), keep=c("age2","age3"))


### Semi parametric analysis

#### Semi markov

crcox <- coxph(Surv(time, status) ~ strata(trans), data = msebmt)


#### Markov

cfcox <- coxph(Surv(Tstart, Tstop, status) ~strata(trans), data = msebmt)



wh1 <- which(msebmt$age2 == 0 & msebmt$age3 == 0)
pat1 <- msebmt[rep(wh1[1], 3), 9:10]
pat1$trans <- 1:3
attr(pat1, "trans") <- tmat
pat1$strata <- pat1$trans


wh2 <- which(msebmt$age2 == 1 & msebmt$age3 == 0)
pat2 <- msebmt[rep(wh2[1], 3), 9:10]
pat2$trans <- 1:3
attr(pat2, "trans") <- tmat
pat2$strata <- pat2$trans

wh3 <- which(msebmt$age2 == 0 & msebmt$age3 == 1)
pat3 <- msebmt[rep(wh3[1], 3), 9:10]
pat3$trans <- 1:3
attr(pat3, "trans") <- tmat
pat3$strata <- pat3$trans



results_semipar <- rpkg::mstatejson(x=cfcox,  qmat=tmat, process="Markov", 
                                    totlos=TRUE, ci.json=TRUE, cl.json=0.95, B.json=2,
                                    variance=FALSE, vartype="greenwood",
                                    covariates_list=list(pat1 ,pat2, pat3 ) , M=2,
                                    jsonpath="data",
                                    name="predictions_EBMT_mstate_fw.json")



results_semipar

results_semipar$timevar[[1]][1:10]
results_semipar$Nats
results_semipar$atlist
results_semipar$tmat









