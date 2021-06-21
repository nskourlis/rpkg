###pro-function 1 msm###

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data_msm PARAM_DESCRIPTION
#' @param id_msm PARAM_DESCRIPTION
#' @param timevar PARAM_DESCRIPTION
#' @param scale_msm PARAM_DESCRIPTION, Default: 1
#' @param Ntransitions PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname freq_func_total_msm
#' @export 
freq_func_total_msm<- function(data_msm, id_msm, timevar, scale_msm=1, Ntransitions) {
  
  freq=list()
  
  data_msm$id=id_msm
  
  
  if (is.null(data_msm$state)) return("provide a variable with the name state")
  
  timevar_used=timevar/scale_msm
  
  p=1
  
  for (j in 1:length(timevar_used) ) {
    
    datanew=data_msm[which(data_msm$years<= timevar_used[j]),]
    
    fr=as.matrix(statetable.msm(datanew$state,datanew$id,data=datanew))
    
    for (i in 1:ncol(fr)) {
      colnames(fr)[i]<-paste0("State"," ",   attributes(fr)$dimnames$to[i]) 
    }
    
    freq[[p]] =  as.data.frame(t(apply(fr, 2, function(x) sum(x))))
    colnames(freq[[p]])<- colnames(fr)
    
    p=p+1
  }
  
  frequence=  bind_rows(freq, .id = "column_label")
  frequence[is.na(frequence)]<-0
  
  vector=names(frequence)
  vectornew=vector[order(vector)]
  
  frequencies=frequence[ , vectornew]
  
  h= as.data.frame(matrix(nrow = length(timevar), ncol=Ntransitions, NA ))
  
  for (i in 1:Ntransitions) {
    colnames(h)[i]=paste0("h",i) 
  }
  
  frequencies=cbind(frequencies,h,timevar)
  frequencies
}

#####################################################################################


###pro function 2 regular freq_tot##########



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param msdata PARAM_DESCRIPTION
#' @param msid PARAM_DESCRIPTION
#' @param names_of_ststates PARAM_DESCRIPTION
#' @param values_ststates PARAM_DESCRIPTION
#' @param names_of_nastates PARAM_DESCRIPTION
#' @param values_nastates PARAM_DESCRIPTION
#' @param names_of_abstates PARAM_DESCRIPTION
#' @param values_abstates PARAM_DESCRIPTION
#' @param names_of_transitions PARAM_DESCRIPTION
#' @param values_of_transitions PARAM_DESCRIPTION
#' @param time PARAM_DESCRIPTION
#' @param timevar PARAM_DESCRIPTION
#' @param scale_inner PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname freq_func_total
#' @export 
freq_func_total <- function(msdata,msid,names_of_ststates, values_ststates,
                            names_of_nastates, values_nastates,
                            names_of_abstates,values_abstates,
                            names_of_transitions,values_of_transitions,
                            time, timevar,scale_inner=1) {
  

  library(viridis)
  library(reshape2)
  library(tidyverse)
  library(plyr)
  library(dplyr)
  
  fr=list()
  
  for (h in 1:length(timevar)) {
    #  from=dat$from, to=dat$to, trans=dat$trans,time=dat$rfstime, 
    
    # data<- read.csv("C:/Users/niksko/Desktop/mstate/jsonread/msset_ebmt.csv",header=TRUE, sep=";")
    options(scipen = 999)  
    data=msdata
    #attach(data)
    #names(data)
    
    data$time_tot=time       #######################################
    
    dis_id=vector()
    dis_id=unique(msid)   ########################
    
    listid=list()
    
    for (i in 1:length(dis_id)) {
      listid[[i]]=as.data.frame(data[data$id==dis_id[i],])
    }
    
    ############################################################################
    #Frequency of remaining in starting state until a certain time point##
    ############################################################################
    #Names of starting states
    names_ststates=names_of_ststates   ##############################################
    #Values of non absorbing states
    v_ststate=values_ststates         ##################################################  
    max(v_ststate)
    
    v_st=array(dim=c(length(dis_id),1,max(v_ststate)),"NA")
    for (k in v_ststate) {
      for (i in 1:length(dis_id)) {
        if (length(which(listid[[i]]$from==k & listid[[i]]$status==1 & listid[[i]]$time_tot<=timevar[h]*scale_inner))==0) {v_st[i,1,k]=TRUE} else {v_st[i,1,k]=FALSE}
      }
    }
    freq_stay_st=vector()
    for (k in v_ststate) {
      freq_stay_st[k]=length(which(v_st[,1,k]=="TRUE"))
    }
    freq_stay_st=freq_stay_st[!is.na(freq_stay_st)]
    freq_stay_st=matrix(nrow=1,ncol=length(freq_stay_st),freq_stay_st)
    freq_stay_st=freq_stay_st
    colnames(freq_stay_st)<-names_ststates
    freq_stay_st
    
    ############################################################################
    #Frequency of non absorbing intermediate states until a certain time point##
    ############################################################################
    #Names of non absorbing state
    names_nastates=names_of_nastates   ####################################
    #Values of non absorbing states
    v_nastate=values_nastates         ######################################
    
    max(v_nastate)
    
    v_na=array(dim=c(length(dis_id),1,max(v_nastate)),"NA")
    for (k in v_nastate) {
      for (i in 1:length(dis_id)) {
        if (length(which(listid[[i]]$from==k & listid[[i]]$status==1 & listid[[i]]$time_tot<=timevar[h]*scale_inner))==0 &
            length(which(listid[[i]]$from!=k & listid[[i]]$to==k & listid[[i]]$status==1 & listid[[i]]$time_tot<=timevar[h]*scale_inner))!=0) {
          v_na[i,1,k]=TRUE} else {v_na[i,1,k]=FALSE
          }
      }
    }
    
    freq_stay_na=vector()
    for (k in v_nastate) {
      freq_stay_na[k]=length(which(v_na[,1,k]=="TRUE"))
    }
    freq_stay_na=freq_stay_na[!is.na(freq_stay_na)]
    freq_stay_na=matrix(nrow=1,ncol=length(freq_stay_na),freq_stay_na)
    colnames(freq_stay_na)<-names_nastates
    freq_stay_na=freq_stay_na
    
    freq_stay_na
    
    ############################################################################
    #Frequency of ending up in absorbing states until a certain time point##
    ############################################################################
    #Names of  absorbing state
    names_abstates=names_of_abstates   #####################################
    #Values of non absorbing states 
    v_abstate=values_abstates           ###################################
    max(v_abstate)
    
    #Frequency of remaining in non absorbing states until a certain time point
    freq_stay_ab=vector()
    
    for (k in v_abstate) {
      
      freq_stay_ab[k]= length(which(data$to==k & data$status==1 & data$time_tot<=timevar[h]*scale_inner ))
    }
    
    freq_stay_ab=freq_stay_ab[!is.na(freq_stay_ab)]
    freq_stay_ab=matrix(nrow=1,ncol=length(freq_stay_ab),freq_stay_ab)
    freq_stay_ab=freq_stay_ab
    colnames(freq_stay_ab)<-names_abstates
    freq_stay_ab
    ######################################################
    frequencies_stay=cbind(freq_stay_st,freq_stay_na,freq_stay_ab)
    # frequencies_stay[ , order(names(frequencies_stay))]
    frequencies_stay
    #########################################################
    
    #Names of transitions
    names_transitions= names_of_transitions
    
    #Values of transitions
    v_trans=values_of_transitions
    
    trans=array(dim=c(length(dis_id),length(v_trans)),"NA")
    for (k in 1:length(v_trans)) {
      for (i in 1:length(dis_id)) {
        if (length(which(listid[[i]]$trans==k & listid[[i]]$status==1 & listid[[i]]$time_tot<=timevar[h]*scale_inner ))==1) 
        {trans[i,k]=TRUE} else {trans[i,k]=FALSE
        }
      }
    }
    
    freq_trans=vector()
    for (k in 1:length(v_trans)) {
      freq_trans[k]=length(which(trans[,v_trans[k]]=="TRUE"))
    }
    
    freq_trans=matrix(nrow=1,ncol=length(freq_trans),freq_trans)
    colnames(freq_trans)<-names_transitions
    freq_trans
    
    results=as.data.frame(cbind(frequencies_stay,freq_trans,timevar[h]))
    colnames(results)[ncol(results)]="timevar"
    
    fr[[h]]=results
    
  }
  
  
  fr_time=bind_rows(fr, .id = "time_label")
  
  
  
}

#######################################################################################


####msboxes###################


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param id PARAM_DESCRIPTION
#' @param yb PARAM_DESCRIPTION
#' @param xb PARAM_DESCRIPTION
#' @param boxwidth PARAM_DESCRIPTION
#' @param boxheight PARAM_DESCRIPTION
#' @param time_study PARAM_DESCRIPTION
#' @param vartime PARAM_DESCRIPTION
#' @param tmat. PARAM_DESCRIPTION
#' @param scale PARAM_DESCRIPTION
#' @param msm PARAM_DESCRIPTION, Default: FALSE
#' @param jsonpath PARAM_DESCRIPTION, Default: '~'
#' @param name PARAM_DESCRIPTION, Default: 'msboxes_R.json'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname msboxes_R
#' @export 
msboxes_R<- function(data,id,yb,xb,boxwidth,boxheight,time_study,vartime, tmat., scale, msm=FALSE,  jsonpath="~",name="msboxes_R.json" ) {
   
   
   library(viridis)
   library(reshape2)
   library(tidyverse)
   library(plyr)
   library(dplyr)
  
  ##########Read info from transition matrix #############################
  
  ####Number of transitions and number of states##############
  tmat_inner=tmat.
  tmat_inner[is.na(tmat_inner)] <- 0
  ntransitions=max(tmat_inner)
  nstates=ncol(tmat_inner)
  states=c(seq(1:nstates))
  #### Transitions values
  u=list()
  for (i in 1:nrow(tmat_inner)) {
    u[[i]]=unique(tmat_inner[i,][which(!is.na(tmat_inner[i,]))])
  }
  u_unlist=sort(unlist(u, recursive = TRUE, use.names = TRUE))
  transitions=u_unlist[u_unlist!=0]
  ######## Categorizing states###############
  state_kind=vector()
  for (i in 1:nstates) {
    if (length(which(tmat_inner[,i]==0))==nrow(tmat_inner)) {state_kind[i]="Starting"}
    else if (length(which(tmat_inner[i,]==0))==ncol(tmat_inner)) {state_kind[i]="Absorbing"}
    else {state_kind[i]="Intermediate"}
  }
  st_states=which(state_kind=="Starting")
  na_states=which(state_kind=="Intermediate")
  ab_states=which(state_kind=="Absorbing")
  
  #### preparing mat to save in json
  mat=tmat.
  mat[is.na(mat)] <- 0
  ntransitions=length(mat[which(mat>0)])
  nstates=ncol(mat)
  states=c(seq(1:nstates))
  p=1
  trmat=matrix(nrow=nstates, ncol=nstates, NA)
  
  for (i in 1:nstates) {
    for (k in 1:nstates) {
      
      if  (mat[i,k]>0)  {  
        
        trmat[i,k]=as.integer(p)
        p=p+1
      } 
      else if (mat[i,k]>0)  {  
        trmat[i,k]="NA"
      }
      
    }
  }

  ##############################################################################
  
  ##########Read info from data #############################
  
#  ntransitions=max(data$trans)
#  nstates= max(data$to)
#  transitions=order(unique(data$trans))  
#
#  d1=unique(from)
#  d2=unique(to)
#  states=unique(c(d1,d2))
#  states=order(states)
#  st_states=vector()
#  na_states=vector()
#  ab_states=vector()
#
#  for (i in 1:length(states)) {
#    if (length(which(data$from==i))==0) {ab_states[i]=states[i]}   }
#      ab_states=ab_states[!is.na(ab_states)]
#
#  for (i in 1:length(states)) {
#    if (length(which(data$to==i))==0) {st_states[i]=states[i]}   }
#    st_states=st_states[!is.na(st_states)]
#
#  na_states=states[which(states!=ab_states & states!=st_states)]
###################################################################################
tname=vector()
for (i in 1: ntransitions) {
  tname[i]=paste0("h",i) }

statename=vector()
for (i in 1: nstates) {
  statename[i]=paste0("State",i) }


tr_start_state=vector()
for (k in 1:ntransitions) {
  tr_start_state[k]=which(tmat_inner == k, arr.ind = TRUE)[1,1]
}

tr_end_state=vector()
for (k in 1:ntransitions) {
  tr_end_state[k]=which(tmat_inner == k, arr.ind = TRUE)[1,2]
}

#########Count initials in each state- Freq total does that anyway###########

#data=data[order(data$id,data$from,data$to),]
#
#dis_id=vector()
#dis_id=unique(data$id)   ########################
#
#listid=list()
#
#for (i in 1:length(dis_id)) {
#  listid[[i]]=as.data.frame(data[data$id==dis_id[i],])
#}
#
#init=array(dim=c(length(dis_id),1,nstates),"NA")
#for (k in 1:nstates) {
#  for (i in 1:length(dis_id)) {
#    if (length(which(listid[[i]]$from[1]==k))!=0) {init[i,1,k]=TRUE} else {init[i,1,k]=FALSE}
#  }
#}
#init_freq=vector()
#for (k in 1:nstates) {
#  init_freq[k]=length(which(init[,1,k]=="TRUE"))
#}
#init_freq

names_ststates=vector(); names_nastates=vector(); names_abstates=vector();names_transitions=vector();

for (i in st_states) {
  names_ststates[i]=paste0("State",i)
}
names_ststates=names_ststates[which(!is.na(names_ststates))]

for (i in na_states) {
  names_nastates[i]=paste0("State",i)
}
names_nastates=names_nastates[which(!is.na(names_nastates))]

for (i in ab_states) {
  names_abstates[i]=paste0("State",i)
}
names_abstates=names_abstates[which(!is.na(names_abstates))]

for (i in transitions) {
  names_transitions[i]=paste0("h",i)
}

#source(paste0(sourcefile,"/","freq_total.R"),local=T)$value
#source(paste0(sourcefile,"/","freq_total_msm.R"),local=T)$value

#source("C:/Users/niksko/Desktop/mstate3/R_flexsurv/functions/freq_total.R",local=T)$value
#source("C:/Users/niksko/Desktop/mstate3/R_flexsurv/functions/freq_total_msm.R",local=T)$value


if (msm==FALSE)  {
  
  p=freq_func_total(msdata=data,msid=id,
                    names_of_ststates=names_ststates, values_ststates=st_states,
                    names_of_nastates=names_nastates, values_nastates=na_states,
                    names_of_abstates=names_abstates, values_abstates=ab_states,
                    names_of_transitions=names_transitions,values_of_transitions=transitions,
                    time=time_study, timevar=vartime, scale_inner=scale)
  
  p$time_label= as.numeric(p$time_label)
  p=as.matrix(p)
  
  msm_inner=0
}

if (msm==TRUE)  {
 
  p= freq_func_total_msm(data_msm=data, id_msm= id, timevar=vartime, scale_msm=scale, Ntransitions=ntransitions) 
  names(p)[names(p) == "column_label"] <- "time_label"
  
  p$time_label= as.numeric(p$time_label)
  options(stringsAsFactors = TRUE)
  p=as.data.frame(p)
  p=as.matrix(p)
  msm_inner=1
}
 
 #####
 
 y=yb
 x=xb
 lefttoright=vector()
 toptobottom=vector()
 y1=vector()
 x1=vector()
 y2=vector()
 x2=vector()
 arrowstexty=vector()
 arrowstextx=vector()
 gradient=vector()
 cutoff=0.5
 textleft=vector()
 
 
 for (i in 1:ntransitions) {
    if (x[tr_start_state[i]]>x[tr_end_state[i]]) {textleft[i]=-1} 
    else {textleft[i]=1}
 }
 
 for (i in 1:ntransitions) {
    if (x[tr_start_state[i]]<x[tr_end_state[i]]) {lefttoright[i]=1} 
    else {lefttoright[i]=-1}
 }
 
 for (i in 1:ntransitions) {
    if (y[tr_start_state[i]]>y[tr_end_state[i]]) {toptobottom[i]=-1} 
    else {toptobottom[i]=1}
 }
 
 for (i in 1:ntransitions) {
    
    gradient[i]=abs(y[tr_start_state[i]]-y[tr_end_state[i]])/abs(x[tr_start_state[i]]-x[tr_end_state[i]])
 }
 
 
 for (i in 1:ntransitions) {
##Horizontal
if (y[tr_start_state[i]]==y[tr_end_state[i]]) {
    y1[i]=y[tr_start_state[i]]
    x1[i]=x[tr_start_state[i]]+lefttoright[i]*boxwidth/2
    y2[i]=y[tr_end_state[i]]
    x2[i]=x[tr_end_state[i]]  -lefttoright[i]*boxwidth/2
    arrowstexty[i]=y[tr_start_state[i]]+0.01
    arrowstextx[i]=(x[tr_start_state[i]]+x[tr_end_state[i]])/2
}
##Vertical
    
else if (x[tr_start_state[i]]==x[tr_end_state[i]]) {
   y1[i]=y[tr_start_state[i]]+toptobottom[i]*boxheight/2
   x1[i]=x[tr_start_state[i]]
   y2[i]=y[tr_end_state[i]]-toptobottom[i]*boxheight/2
   x2[i]=x[tr_end_state[i]]
   arrowstexty[i]=(y[tr_start_state[i]]+y[tr_end_state[i]])/2
   arrowstextx[i]=x[tr_start_state[i]]-0.01
}
##Dradient
else{
    
    if (gradient[i]<cutoff) {
    
      y1[i]=y[tr_start_state[i]]
      x1[i]=x[tr_start_state[i]]+lefttoright[i]*boxwidth/2
      y2[i]=y[tr_end_state[i]]
      x2[i]=x[tr_end_state[i]]  -lefttoright[i]*boxwidth/2
      arrowstexty[i]=(y[tr_start_state[i]]+y[tr_end_state[i]])/2 
      arrowstextx[i]=(x[tr_start_state[i]]+x[tr_end_state[i]])/2 +textleft[i]*0.02
    }
 
    else {
       y1[i]=y[tr_start_state[i]]+toptobottom[i]*boxheight/2
       x1[i]=x[tr_start_state[i]]
       y2[i]=y[tr_end_state[i]]-toptobottom[i]*boxheight/2
       x2[i]=x[tr_end_state[i]]
       arrowstexty[i]=(y[tr_start_state[i]]+y[tr_end_state[i]])/2 
       arrowstextx[i]=(x[tr_start_state[i]]+x[tr_end_state[i]])/2 +textleft[i]*0.02
    }
   
}
 }
 
# xb=xb-boxwidth/2
# yb=yb+boxheight/2
 
 #list(x1,y1,x2,y2,arrowstextx,arrowstexty)
 arrows=list(x1=x1,y1=y1,x2=x2,y2=y2)
 arrowstext=list(x=arrowstextx,y=arrowstexty)
 
# for (i in 1:max(p$time_label)) {
# p$timevar[i]=p$timevar[i]/scale
# }
 
 results=list(Nstates=nstates,Ntransitions=ntransitions, 
      xvalues=xb, yvalues=yb,
      boxwidth=boxwidth, boxheight=boxheight,
      statenames=statename, transnames=tname,
      arrows=arrows,arrowstext=arrowstext,tmat=trmat,msm=msm_inner,frequencies=p)
 
 library(RJSONIO)
 exportJson <- toJSON(results,  dataframe = "rows",matrix="rowmajor",force=TRUE, complex="list",flatten=TRUE)
 exportJson
 write(exportJson, paste0(jsonpath ,"/", name  ) )
 
 results
 
}


#results12=msboxes_R(data=msebmt,id= msebmt$id, yb=c(0.8,0.8,0.8,0.4,0.4,0.4),
#                    xb=c(0.2,0.5,0.8,0.2,0.5,0.8),boxwidth=0.1,boxheight=0.1,
#                    tmat.= tmat, time_study=msebmt$time,vartime=c(seq(0,2000,by=200)),scale=1) 
#


