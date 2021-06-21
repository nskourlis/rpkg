

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION, Default: cfcox
#' @param qmat PARAM_DESCRIPTION, Default: tmat
#' @param process PARAM_DESCRIPTION, Default: 'Markov'
#' @param totlos PARAM_DESCRIPTION, Default: TRUE
#' @param ci.json PARAM_DESCRIPTION, Default: TRUE
#' @param cl.json PARAM_DESCRIPTION, Default: 0.95
#' @param B.json PARAM_DESCRIPTION, Default: 100
#' @param variance PARAM_DESCRIPTION, Default: TRUE
#' @param vartype PARAM_DESCRIPTION, Default: 'aalen'
#' @param covariates_list PARAM_DESCRIPTION, Default: list()
#' @param M PARAM_DESCRIPTION, Default: 100
#' @param jsonpath PARAM_DESCRIPTION, Default: '~'
#' @param name PARAM_DESCRIPTION, Default: 'predictions_R_mstate.json'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stringi]{stri_sort}}
#' @rdname mstatejson
#' @export 
#' @importFrom stringi stri_sort
mstatejson <- function(x=cfcox, qmat=tmat, process="Markov", 
                       totlos=TRUE, ci.json=TRUE, cl.json=0.95, B.json=100,
                       variance=TRUE, vartype="aalen",
                       covariates_list=list(), M=100,  jsonpath="~",name="predictions_R_mstate.json"  )  {
  
  options(scipen = 999,"digits"=10)
  
  if (!require(survival)) install.packages("survival")
  if (!require(mstate)) install.packages("mstate")
  if (!require(tidyverse)) install.packages("tidyverse")
  library("survival")
  library("mstate")
  library("tidyverse")
  
  try( if (length(covariates_list)==0) stop("The user must provide covariate list patterns"))
  
  ################################################################  
  mat=qmat
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
  tmatlist=list()
  tmatlist[[1]]=trmat
  
  tr_start_state=vector()
  for (k in 1:ntransitions) {
    tr_start_state[k]=which(trmat == k, arr.ind = TRUE)[1,1]
  }
  
  tr_end_state=vector()
  for (k in 1:ntransitions) {
    tr_end_state[k]=which(trmat == k, arr.ind = TRUE)[1,2]
  }
  
  
  if (process=="Markov") clock="forward"
  if (process=="semiMarkov") clock="reset"
  
  
  ################################################################### 
  
  library("msm")
  library("stringi")
  library("RJSONIO")
  
  pred=list() 
  
  for (g in 1:length(covariates_list)) {
    
    #   g=1
    
    try( if (process!="Markov" & process!="semiMarkov") stop("Process must be either Markov or semiMarkov"))
    
    
    ### Calculate nstates, ntransitions, non absorbing states #### 
    nstates=ncol(qmat)
    states=c(seq(1:nstates))
    
    qmat2=qmat
    qmat2[is.na(qmat2)] <- 0
    Ntransitions=length(qmat2[which(qmat2>0)])
    
    
    modHaz <-  msfit(object=x, newdata=covariates_list[[g]],
                     variance=variance,vartype=vartype, trans=qmat)
    
    ###################################    
    #      ##Cumulative hazard #####
    ####################################      
    #semi parametric semi markov
    transr= reshape(modHaz$Haz, idvar = "time", timevar = "trans", direction = "wide")
    transr$timevar=transr$time
    transr=transr[,-1]
    namesh=vector()
    for (i in 1:Ntransitions) {namesh[i]=paste0("Haz_",tr_start_state[i],"_to_",tr_end_state[i])}
    names(transr)[1:Ntransitions]=namesh
    
    ############################################################################
    #### Timevar#####
    
    timevar=transr$time
    
    
    ########################################################################
    
    ########## probabilities ################
    
    if (process=="Markov") {
      
      probm=list()
      
      
      ### Prob markov
      probm <- probtrans(modHaz, predt = 0, direction = "forward")[[1]]
      probm= probm[,    which(!startsWith(names(probm), "se") )] 
      probm= probm[which(probm$pstate1>=0),] 
      probm$timevar=probm$time
      probm=probm[,-1]
      namesprob=vector()
      
      for (k in tr_start_state[1]) {
        for (j in states) {
          colnames(probm)[j]=paste0("P_",k,"_to_",j)
        }
        colnames(probm)[nstates+1]  ="timevar"  
      } 
      probm=probm[-1,]
      
      #############################################################
      #### Length of stay
      
      loslist=list()
      
      
      if (totlos==TRUE) {
        
        los_array= array(dim=c(length(timevar),nstates+1,nstates),NA)
        
        
        probm_all <- probtrans(modHaz, predt = 0, direction = "forward")
        
        for (j in 1:(nstates+1)) {
          
          for (r in 1:(nstates)) {
            
            for (i in 1:length(timevar)) {
              
              los.mstate <- ELOS(probm_all,tau=timevar[i])
              
              los_array[i,j,r]=rbind(t(los.mstate),timevar=timevar[i]  ) [j,r]
              
            }
          }
        }  
        
        
        
        for (k in states) {
          
          loslist[[k]]=los_array[,,k]
          loslist[[k]]=as.data.frame(loslist[[k]])
          
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist[[k]])[j]=paste0("Los_",k,"_to_",j)
            
          }
          
          colnames(loslist[[k]])[nstates+1]  ="timevar"  
          
        }
        
        loslist  
        
      }
      
    }   
    
    if (process=="semiMarkov") {
      
      ### Prob semi markov
      set.seed(1234)
      tv <- unique(modHaz$Haz$time)
      probm=list()
      p=1
      for (start in states) {
        probm[[p]]=  mssample(Haz=modHaz$Haz[which(modHaz[[1]][,2]>=0),], trans=qmat,
                              history=list(state=start,time=0,tstate=NULL),
                              beta.state=NULL, clock=clock, output="state",
                              M = 10, tvec = timevar, cens=NULL, do.trace=NULL)
        
        probm[[p]]$timevar=probm[[p]]$time
        probm[[p]]=probm[[p]][,-1]
        p=p+1   
      }
      
      
      for (k in states) {
        for (j in states) {
          colnames(probm[[k]])[j]=paste0("P_",k,"_to_",j)
        }
        colnames(probm[[k]])[nstates+1]  ="timevar"  
      } 
      probm
      #probm=probm[-1,]
      
      
      loslist=list()
      
      
      if (totlos==TRUE) {
        
        los_array= array(dim=c(length(timevar),nstates+1,nstates),NA)
        
        
        probm_all <- probtrans(modHaz, predt = 0, direction = "forward")
        
        for (j in 1:(nstates+1)) {
          
          for (r in 1:(nstates)) {
            
            for (i in 1:length(timevar)) {
              
              los.mstate <- ELOS(probm_all,tau=timevar[i])
              
              los_array[i,j,r]=rbind(t(los.mstate),timevar=timevar[i]  ) [j,r]
              
            }
          }
        }  
        
        
        
        for (k in states) {
          
          loslist[[k]]=los_array[,,k]
          loslist[[k]]=as.data.frame(loslist[[k]])
          
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist[[k]])[j]=paste0("Los_",k,"_to_",j)
            
          }
          
          colnames(loslist[[k]])[nstates+1]  ="timevar"  
          
        }
        loslist  
        
      }
      
    }
    is.cumhaz=1
    
    
    pred[[g]]=            list(probs=probm,
                               trans=transr, is.cumhaz=is.cumhaz,
                               los=loslist, timevar=timevar)   
  }
  
  rm(list=ls(pattern="^forMSTATE"))
  
  
  
  ########################################################################################################
  ######## Declare the list elements even if they stay empty ####
  
  pjson    =list()
  hjson=list()
  losjson=list()
  
  
  
  ########################################################################################################
  ####  Probabilities
  ########################################################################################################
  
  ##### Main estimates 
  if (process=="semiMarkov") {
    
    pjson=list()
    
    if (length(covariates_list)<=1) {
      p=1
      for (k in 1:nstates) {
        for (i in 1:nstates) {
          
          pjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$probs[[k]][,i])))
          assign(paste0("forMSTATE",names(pred[[1]]$probs[[k]])[i]),  pjson[[p]])
          p=p+1  
        }
      }
    }
    
    if (length(covariates_list)>1) {
      
      p=1
      for (k in 1:nstates) {
        for (i in 1:nstates) {
          
          tmplist=matrix(nrow =length(covariates_list) , ncol= length(as.vector(pred[[1]]$probs[[k]][,i])))
          tmplist[1,]=as.vector(pred[[1]]$probs[[k]][,i])
          
          for (j in 2:length(covariates_list)) {
            tmplist[j,]=as.vector(pred[[j]]$probs[[k]][,i])  
          }
          pjson[[p]]=tmplist
          
          assign(paste0("forMSTATE",names(pred[[1]]$probs[[k]])[i]),  pjson[[p]])
          
          p=p+1  
          
          
        }
      }
    }
    
    pvector=vector()
    pvector=ls(pattern = "forMSTATEP_._to_.$")  
    
    pvector=sub("forMSTATE","",ls(pattern = "^forMSTATEP_._to_.$") )
    
    pvector=stringi::stri_sort(pvector, numeric = TRUE)
    
    names(pjson)=pvector[1:(p-1)]
    
    pjson
    
  }
  
  if (process=="Markov") {
    
    pjson=list()
    
    if (length(covariates_list)<=1) {
      p=1
      for (k in 1:nstates) {
        for (i in 1:nstates) {
          
          pjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$probs[,i])))
          assign(paste0("forMSTATE",names(pred[[1]]$probs)[i]),  pjson[[p]])
          p=p+1  
        }
      }
    }
    
    if (length(covariates_list)>1) {
      
      p=1
      
      for (i in 1:nstates) {
        
        
        
        
        #          if (trans.specific=="Yes") { 
        #          
        #             tmplist=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
        #         
        #             tmplist[1,]= as.vector(pred[[1]]$probs[,i]) #[-length(as.vector(pred[[1]]$probs[,i]))]
        #      
        #           
        #           
        #           for (j in 2:length(covariates_list)) {
        #               tmplist[j,]=as.vector(pred[[j]]$probs[,i])#[-length(as.vector(pred[[j]]$probs[,i]))]
        #           }
        #       }    
        
        #     if (trans.specific=="No") { 
        
        tmplist=matrix(nrow =length(covariates_list) , ncol= length(timevar)) #+1
        
        tmplist[1,]=as.vector(pred[[1]]$probs[,i])
        
        
        
        for (j in 2:length(covariates_list)) {
          tmplist[j,]=as.vector(pred[[j]]$probs[,i])
        }
        #      } 
        
        pjson[[p]]=tmplist
        
        assign(paste0("forMSTATE",names(pred[[1]]$probs)[i]),  pjson[[p]])
        
        p=p+1  
        
        
      }
      
    }
    
    pvector=vector()
    pvector=ls(pattern = "forMSTATEP_._to_.$")  
    
    pvector=sub("forMSTATE","",ls(pattern = "^forMSTATEP_._to_.$") )
    
    pvector=stringi::stri_sort(pvector, numeric = TRUE)
    
    names(pjson)=pvector[1:(p-1)]
    
    pjson
    
  }
  
  
  #################################################################################################################
  
  ### Transitions main estimates ###
  hjson=list()
  
  #  if (trans.specific == "No") {
  
  if (length(covariates_list)<=1) {
    
    for (k in 1:ntransitions) {
      
      if (process=="Markov") {
        
        tmplisth=matrix(nrow =length(covariates_list) , ncol= length(timevar))
      }
      else if (process=="semiMarkov") {
        
        tmplisth=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
        
      }       
      
      tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
      hjson[[1]]=tmplisth
      assign(paste0("forMSTATE",names(pred[[1]]$trans)[i]),  hjson[[1]])
    }
    
  }
  
  if (length(covariates_list)>1) {
    
    for (k in 1:ntransitions) {
      
      if (process=="Markov") {
        
        tmplisth=matrix(nrow =length(covariates_list) , ncol= length(timevar))
      }
      else if (process=="semiMarkov") {
        
        tmplisth=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
        
      }
      
      tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
      
      
      
      for (j in 2:length(covariates_list)) {
        
        tmplisth[j,]=as.vector(pred[[j]]$trans[,k])  
      }
      
      
      hjson[[k]]=tmplisth
      assign(paste0("forMSTATE",names(pred[[1]]$trans)[k]),  hjson[[k]])
    }
  }
  
  hvector=vector()
  hvector=ls(pattern = "forMSTATEHaz_._to_.$")  
  
  hvector=sub("forMSTATE","",ls(pattern = "forMSTATEHaz_._to_.$") )
  
  hvector=stringi::stri_sort(hvector, numeric = TRUE)
  
  names(hjson)=hvector
  
  # }
  
  #     if (trans.specific == "Yes") {
  #       
  #       if (length(covariates_list)<=1) {
  #         
  #         for (k in 1:ntransitions) {
  #           
  #           if (process=="Markov") {
  #             
  #             tmplisth=matrix(nrow =length(covariates_list) , ncol= length(timevar))
  #           }
  #           else if (process=="semiMarkov") {
  #             
  #             tmplisth=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
  #             
  #           }   
  #           
  #           
  #           tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
  #           
  #           hjson[[1]]=tmplisth
  #           
  #           assign(paste0("forMSTATE",names(pred[[1]]$trans)[i]),  hjson[[1]])
  #         }
  #         
  #       }
  #       if (length(covariates_list)>1) {
  #         
  #         
  #         
  #         for (k in 1:ntransitions) {
  #           
  #           if (process=="Markov") {
  #             
  #             tmplisth=matrix(nrow =length(covariates_list) , ncol= length(timevar))
  #           }
  #           else if (process=="semiMarkov") {
  #             
  #             tmplisth=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
  #             
  #           }   
  #           
  #           tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
  #           
  #           for (j in 2:length(covariates_list)) {
  #             
  #             tmplisth[j,]=as.vector(pred[[j]]$trans[,k])  
  #           }
  #           
  #           
  #           hjson[[k]]=tmplisth
  #           assign(paste0("forMSTATE",names(pred[[1]]$trans)[k]),  hjson[[k]])
  #         }
  #       }
  #       
  #       hvector=vector()
  #       hvector=ls(pattern = "^forMSTATEh.$")  
  #       
  #       hvector=sub("forMSTATE","",ls(pattern = "^forMSTATEh.$") )
  #       
  #       hvector=stringi::stri_sort(hvector, numeric = TRUE)
  #       
  #       names(hjson)=hvector
  #       
  #     }
  #     
  
  ##############################################################################################################
  
  if (process=="Markov") { 
    
    
    if (totlos==TRUE) { 
      
      
      #   if (trans.specific == "No") {
      ### Los Main estimates #####
      
      losjson=list()
      
      if (length(covariates_list)<=1) {
        p=1
        for (k in 1:nstates) {
          for (i in 1:nstates) {
            
            losjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$los[[k]][,i])))
            assign(paste0("forMSTATE",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
            p=p+1  
          }
        }
      }
      
      if (length(covariates_list)>1) {
        
        p=1
        for (k in 1:nstates) {
          for (i in 1:nstates) {
            
            tmplistlos=matrix(nrow =length(covariates_list) , ncol= length(timevar))
            tmplistlos[1,]=as.vector(pred[[1]]$los[[k]][,i])
            
            for (j in 2:length(covariates_list)) {
              tmplistlos[j,]=as.vector(pred[[j]]$los[[k]][,i])  
            }
            losjson[[p]]=tmplistlos
            
            assign(paste0("forMSTATE",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
            
            p=p+1  
          }
        }
      }
      
      losvector=vector()
      losvector=ls(pattern = "forMSTATELos_._to_.$")  
      
      losvector=sub("forMSTATE","",ls(pattern = "^forMSTATELos_._to_.$") )
      
      losvector=stringi::stri_sort(losvector, numeric = TRUE)
      
      names(losjson)=losvector
      
      #   }
      
      #      if (trans.specific == "Yes") {
      #            ### Los Main estimates #####
      #            
      #            losjson=list()
      #            
      #            if (length(covariates_list)<=1) {
      #              p=1
      #              for (k in 1:nstates) {
      #                for (i in 1:nstates) {
      #                  
      #                  losjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$los[[k]][,i])))
      #                  assign(paste0("forMSTATE",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
      #                  p=p+1  
      #                }
      #              }
      #            }
      #            
      #            if (length(covariates_list)>1) {
      #              
      #              p=1
      #              for (k in 1:nstates) {
      #                for (i in 1:nstates) {
      #                  
      #                  tmplistlos=matrix(nrow =length(covariates_list) , ncol= length(transr$timevar))
      #                  tmplistlos[1,]=as.vector(pred[[1]]$los[[k]][,i])
      #                  
      #                  for (j in 2:length(covariates_list)) {
      #                    tmplistlos[j,]=as.vector(pred[[j]]$los[[k]][,i])  
      #                  }
      #                  losjson[[p]]=tmplistlos
      #                  
      #                  assign(paste0("forMSTATE",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
      #                  
      #                  p=p+1  
      #                }
      #              }
      #            }
      #            
      #            losvector=vector()
      #            losvector=ls(pattern = "forMSTATELos_._to_.$")  
      #            
      #            losvector=sub("forMSTATE","",ls(pattern = "^forMSTATELos_._to_.$") )
      #            
      #            losvector=stringi::stri_sort(losvector, numeric = TRUE)
      #            
      #            names(losjson)=losvector
      #            
      #            
      #            
      #            }
    }
  }
  
  
  
  ##########################################################################################
  
  ### Timevar###
  
  timejson=list()
  
  timejson[[1]]=timevar
  
  ######################################################################################################
  
  atlistjson=list()
  atlistmatrix=matrix(nrow=1, ncol=length(covariates_list), NA)
  
  
  for (j in 1:length(covariates_list) ) {
    
    name_paste=list()
    
    for (i in 1:ncol(covariates_list[[j]]) )  {
      
      name_paste[[i]]= paste0(names(covariates_list[[j]])[i]," ", as.character(covariates_list[[j]][1,i])  )
      
    }
    atlistmatrix[1,j]= paste(unlist(name_paste,recursive=FALSE), collapse = ' ')
  }
  
  atlistjson[[1]]=atlistmatrix
  
  
  
  ######################################################################################################
  is.cumhaz=1
  
  final_list=list(timevar=timejson, Nats=length(covariates_list) ,
                  atlist=atlistjson, tmat=tmatlist, is.cumhaz=is.cumhaz,
                  pjson,
                  hjson,
                  losjson
  )
  final_unlist= unlist(final_list, recursive=FALSE)
  
  exportJson <- toJSON(final_unlist,force = TRUE, flatten=TRUE)
  
  write(exportJson, paste0(jsonpath ,"/", name) )
  
  rm(list=ls(pattern="^forMSTATE"))
  
  final_list
  
  
}

