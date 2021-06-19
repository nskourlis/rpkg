
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param model PARAM_DESCRIPTION, Default: cffpm.list
#' @param vartime PARAM_DESCRIPTION, Default: seq(365.25, 730.5, by = 365.25)
#' @param qmat PARAM_DESCRIPTION, Default: tmat
#' @param process PARAM_DESCRIPTION, Default: 'Markov'
#' @param totlos PARAM_DESCRIPTION, Default: TRUE
#' @param ci.json PARAM_DESCRIPTION, Default: TRUE
#' @param cl.json PARAM_DESCRIPTION, Default: 0.95
#' @param B.json PARAM_DESCRIPTION, Default: 100
#' @param tcovs PARAM_DESCRIPTION, Default: NULL
#' @param Mjson PARAM_DESCRIPTION, Default: 100
#' @param variance PARAM_DESCRIPTION, Default: FALSE
#' @param covariates_list PARAM_DESCRIPTION, Default: list(pat1, pat2, pat3)
#' @param jsonpath PARAM_DESCRIPTION, Default: 'C:/Users/niksko/Desktop/mstate3/datasets/json/flexsurv/json_present_flexsurv'
#' @param name PARAM_DESCRIPTION, Default: 'predictions_EBMT_flex_fw.json'
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
#' @rdname flexsurv_json
#' @export 
#' @importFrom stringi stri_sort
flexsurv_json <- function( model=cffpm.list, vartime=seq(365.25,730.5,by=365.25), qmat=tmat, process="Markov",
                           totlos=TRUE, ci.json=TRUE, cl.json=0.95, B.json=100, tcovs=NULL,
                           Mjson=100, variance=FALSE,
                           covariates_list=list(pat1,pat2,pat3), 
                           jsonpath="C:/Users/niksko/Desktop/mstate3/datasets/json/flexsurv/json_present_flexsurv",
                           name="predictions_EBMT_flex_fw.json" )  {
  
  options(scipen = 999,"digits"=14)
  
  #if (!require(flexsurv)) install.packages("flexsurv")
  if (!require(survival)) install.packages("msm")
  if (!require(mstate)) install.packages("stringi")
  if (!require(tidyverse)) install.packages("RJSONIO")
  library("msm")
  library("stringi")
  library("RJSONIO")
  library("flexsurv")

  
  #  vartime=time/scale
  
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
  #Check <- newArgCheck()
  
  
  pred=list() 
  

  for (g in 1:length(covariates_list)) {
  #  g=1
    
    try( if (process!="Markov" & process!="semiMarkov") stop("Process must be either Markov or semiMarkov"))
    
    ### Calculate nstates, ntransitions #### 
     nstates=ncol(qmat)
      states=c(seq(1:nstates))
    
      tmat2=qmat
      tmat2[is.na(tmat2)] <- 0
     Ntransitions=length(tmat2[which(tmat2>0)])
    
    #################### Probabilities #######################################
    
    pm_array= array(dim=c(length(vartime),nstates+1,nstates),NA)
    pm_array_lci= array(dim=c(length(vartime),nstates+1,nstates),NA)
    pm_array_uci= array(dim=c(length(vartime),nstates+1,nstates),NA)
    
    #g=1
    if (process=="Markov") {
      pm_ci=list()
      for (k in 1:length(vartime)) {
            pm_ci[[k]]=pmatrix.fs(x=model, trans=qmat, t =vartime[k], newdata=covariates_list[[g]],
                            sing.inf = 1e+10, B = B.json, ci = ci.json, cl = cl.json)
  
      }
    }
    
    if (process=="semiMarkov") {
      pm_ci=list()
     #p=1
      for (k in 1:length(vartime)) {
        
        pm_ci[[k]]= pmatrix.simfs(x=model, trans=qmat, t = vartime[k],  newdata =covariates_list[[g]], ci = ci.json,
                                  tcovs = tcovs, B = B.json, cl = cl.json, M = Mjson)
     #   p=p+1  
      }
      
    }
    
    
    if (ci.json==FALSE) {
      for (i in 1: length(vartime)  ) {
        pm_array[i,,]=t(cbind(as.matrix(pm_ci[[i]]),vartime[i] ))
      }
    }
    
    if (ci.json==TRUE) {
      
      for (i in 1: length(vartime)  ) {
        pm_array[i,,]=    t(cbind(as.matrix(pm_ci[[i]]),vartime[i] ))
        pm_array_lci[i,,]=t(cbind(attributes(pm_ci[[i]])$lower,vartime[i]  ))
        pm_array_uci[i,,]=t(cbind(attributes(pm_ci[[i]])$upper,vartime[i]  ))
        
      }
    }
    
    
    pmlist=list()
    
    for (k in states) {
      pmlist[[k]]=pm_array[,,k]
      
      if (length(vartime)==1) {pmlist[[k]]=as.data.frame(t(as.matrix(pmlist[[k]])))}
      else if (length(vartime)>1) {pmlist[[k]]=as.data.frame(as.matrix(pmlist[[k]]))}
      
      # pmlist[[k]]=as.data.frame(as.matrix(pmlist[[k]]))
      #pmlist[[k]]=as.data.frame(t(as.matrix(pmlist[[k]])))
      #pmlist[[k]]=as.data.frame(pmlist[[k]])
    }
    
    for (k in states) {
      for (j in states) {
       
        colnames(pmlist[[k]])[j]=paste0("P_",k,"_to_",j)
      }
      colnames(pmlist[[k]])[nstates+1]  ="timevar"  
    } 
    pmlist
    
    pmlist_lci=list()
    pmlist_uci=list()
    
    
    if (ci.json==TRUE) {
      
      ##################### Probabilities lci #####################################################################
      
      
      for (k in states) {
        pmlist_lci[[k]]=pm_array_lci[,,k]
        
        if (length(vartime)==1) {pmlist_lci[[k]]=as.data.frame(t(as.matrix(pmlist_lci[[k]])))}
        else if (length(vartime)>1) {pmlist_lci[[k]]=as.data.frame(as.matrix(pmlist_lci[[k]]))}
        
        #pmlist_lci[[k]]=as.data.frame(as.matrix(pmlist_lci[[k]]))
        
      }
      
      for (k in states) {
        for (j in states) {
          colnames(pmlist_lci[[k]])[j]=paste0("P_",k,"_to_",j,"_lci")
        }
        colnames(pmlist_lci[[k]])[nstates+1]  ="timevar"  
      } 
      
      pmlist_lci
      
      ##################### Probabilities uci #####################################################################
      
      #pmlist_uci=list()
      for (k in states) {
        pmlist_uci[[k]]=pm_array_uci[,,k]
        
        if (length(vartime)==1) {pmlist_uci[[k]]=as.data.frame(t(as.matrix(pmlist_uci[[k]])))}
        else if (length(vartime)>1) {pmlist_uci[[k]]=as.data.frame(as.matrix(pmlist_uci[[k]]))}
        
      #  pmlist_uci[[k]]=as.data.frame(as.matrix(pmlist_uci[[k]]))
      }
      
      for (k in states) {
        for (j in states) {
          colnames(pmlist_uci[[k]])[j]=paste0("P_",k,"_to_",j,"_uci")
        }
        colnames(pmlist_uci[[k]])[nstates+1]  ="timevar"  
      } 
      
      pmlist_uci
      
      
    }
    
    
    
    
    #################### Cummulative transition intensities #######################################
    
    transm=msfit.flexsurvreg(object=model, t=vartime, newdata = covariates_list[[g]], variance = variance,
                             tvar = "trans", trans=qmat, B = B.json)$Haz
    
    
    
    trans= reshape(transm, idvar = "time", timevar = "trans", direction = "wide")
    trans$timevar=trans$time
    trans=trans[,-1]
    namesh=vector()
    for (i in 1:Ntransitions) {namesh[i]=paste0("Haz_",tr_start_state[i],"_to_",tr_end_state[i])}
    names(trans)[1:Ntransitions]=namesh
    
    
    is.cumhaz=1
    
    #################################################################################
    #################### Total length of stay #######################################
    #################################################################################
    
    #### Markov models###
    
    
    if (process=="Markov") {
      
      los_array= array(dim=c(length(vartime),nstates+1,nstates),NA)
      los_array_lci= array(dim=c(length(vartime),nstates+1,nstates),NA)
      los_array_uci= array(dim=c(length(vartime),nstates+1,nstates),NA)
      
      los_ci=list()
      
      for (k in 1:length(vartime)) {
      
      
          los_ci[[k]]= totlos.fs(x=model, trans=qmat, t = vartime[k], newdata =covariates_list[[g]] , ci = ci.json, tvar = "trans",
                        sing.inf = 1e+5, B = B.json, cl = cl.json)
      
      }
      
      if (ci.json==FALSE) {
        for (i in 1: length(vartime)  ) {
          los_array[i,,]=t(cbind(as.matrix(los_ci[[i]]),vartime[i] ))
        }
      }
      if (ci.json==TRUE) {
        
        for (i in 1: length(vartime)  ) {
          los_array[i,,]=    t(cbind(as.matrix(los_ci[[i]]),vartime[i] ))
          los_array_lci[i,,]=t(cbind(attributes(los_ci[[i]])$lower,vartime[i]  ))
          los_array_uci[i,,]=t(cbind(attributes(los_ci[[i]])$upper,vartime[i]  ))
          
        }
      }
      
      loslist=list()
      
      for (k in states) {
        loslist[[k]]=los_array[,,k]
        
        if (length(vartime)==1) {loslist[[k]]=as.data.frame(t(as.matrix(loslist[[k]])))}
        else if (length(vartime)>1) {loslist[[k]]=as.data.frame(as.matrix(loslist[[k]]))}
        
      }
      
      for (k in states) {
        for (j in states) {
          colnames(loslist[[k]])[j]=paste0("Los_",k,"_to_",j)
        }
        colnames(loslist[[k]])[nstates+1]  ="timevar"  
      } 
      loslist
      
      loslist_lci=list()
      loslist_uci=list()
      
      if (ci.json==TRUE) {
        
        ##################################################################################################
        ##################### Total length of stay lci ##################################################
        #################################################################################################
        
        for (k in states) {
          loslist_lci[[k]]=los_array_lci[,,k]
          
          if (length(vartime)==1) {loslist_lci[[k]]=as.data.frame(t(as.matrix(loslist_lci[[k]])))}
          else if (length(vartime)>1) {loslist_lci[[k]]=as.data.frame(as.matrix(loslist_lci[[k]]))}
          
         # loslist_lci[[k]]=as.data.frame(as.matrix(loslist_lci[[k]]))
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist_lci[[k]])[j]=paste0("Los_",k,"_to_",j,"_lci")
          }
          colnames(loslist_lci[[k]])[nstates+1]  ="timevar"  
        } 
        
        loslist_lci
        
        ################################################################################
        ##################### Total length of stay uci #################################
        ################################################################################# 
        #loslist_uci=list()
        for (k in states) {
          loslist_uci[[k]]=los_array_uci[,,k]
          
          if (length(vartime)==1) {loslist_uci[[k]]=as.data.frame(t(as.matrix(loslist_uci[[k]])))}
          else if (length(vartime)>1) {loslist_uci[[k]]=as.data.frame(as.matrix(loslist_uci[[k]]))}
          
         # loslist_uci[[k]]=as.data.frame(as.matrix(loslist_uci[[k]]))
          
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist_uci[[k]])[j]=paste0("Los_",k,"_to_",j,"_uci")
          }
          colnames(loslist_uci[[k]])[nstates+1]  ="timevar"  
        } 
        
        loslist_uci
      }
    }
    
    
    ###################################################################################################################
    
    ###########################################################################
    ############## Total length of stay #######################################
    
    ###### semi markov models ##############
    
    if (process=="semiMarkov") {
      
      
      los_array= array(dim=c(length(vartime),nstates+1,nstates),NA)
      los_array_lci= array(dim=c(length(vartime),nstates+1,nstates),NA)
      los_array_uci= array(dim=c(length(vartime),nstates+1,nstates),NA)
      
      
      
      los=list()
      los_ci=list()
      p=1
      for (j in states) {
        
        
        for (k in 1:length(vartime)) {
          
          los_ci[[k]]= totlos.simfs(x=model, trans=qmat, t = vartime[k], start = j, 
                                    newdata =covariates_list[[g]], ci = TRUE,
                                    tvar = "trans", tcovs = tcovs, group = NULL,
                                    M = Mjson, B = B.json, cl = cl.json)
          p=p+1  
          
        }
        los[[j]]=los_ci
      }
      
      
      
      if (ci.json==TRUE) {
        
        for (i in 1: length(vartime)  ) {
          for (j in 1: length(states)  ) {
            for (k in 1: length(states)  ) {
              los_array[i,k,j]=t(cbind(as.matrix(los[[j]][[i]][,1][k])))
              los_array_lci[i,k,j]=t(cbind(as.matrix(los[[j]][[i]][,2][k])))
              los_array_uci[i,k,j]=t(cbind(as.matrix(los[[j]][[i]][,3][k])))
            }
            los_array[i,nstates+1,j]=vartime[i]
            los_array_lci[i,nstates+1,j]=vartime[i]
            los_array_uci[i,nstates+1,j]=vartime[i]
            
          }
        }
        
      }    
      
      
      if (ci.json==FALSE) {
        for (i in 1: length(vartime)  ) {
          for (j in 1: length(states)  ) {
            for (k in 1: length(states)  ) {
              los_array[i,k,j]=t(cbind(as.matrix(los[[j]][[i]][k])))
            }
            los_array[i,nstates+1,j]=vartime[i]
          }
        }
        
      } 
      
      
      loslist=list()
      
      for (k in states) {
        loslist[[k]]=los_array[,,k]
        
        if (length(vartime)==1) {loslist[[k]]=as.data.frame(t(as.matrix(loslist[[k]])))}
        else if (length(vartime)>1) {loslist[[k]]=as.data.frame(as.matrix(loslist[[k]]))}
        
        
       # loslist[[k]]=as.data.frame(as.matrix(loslist[[k]]))
      }
      
      for (k in states) {
        for (j in states) {
          colnames(loslist[[k]])[j]=paste0("Los_",k,"_to_",j)
        }
        colnames(loslist[[k]])[nstates+1]  ="timevar"  
      } 
      loslist
      
      
      
      loslist_lci=list()
      loslist_uci=list()
      
      if (ci.json==TRUE) {
        
        ##############################################################################
        ##################### Total length of stay lci ###############################
        ##############################################################################
        
        for (k in states) {
          loslist_lci[[k]]=los_array_lci[,,k]
          
          
          if (length(vartime)==1) {loslist_lci[[k]]=as.data.frame(t(as.matrix(loslist_lci[[k]])))}
          else if (length(vartime)>1) {loslist_lci[[k]]=as.data.frame(as.matrix(loslist_lci[[k]]))}
          
          
         # loslist_lci[[k]]=as.data.frame(as.matrix(loslist_lci[[k]]))
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist_lci[[k]])[j]=paste0("Los_",k,"_to_",j,"_lci")
          }
          colnames(loslist_lci[[k]])[nstates+1]  ="timevar"  
        } 
        
        loslist_lci
        
        ##################################################################################
        ##################### Total length of stay uci ##################################
        ################################################################################
        
        loslist_uci=list()
        for (k in states) {
          loslist_uci[[k]]=los_array_uci[,,k]
          
          if (length(vartime)==1) {loslist_uci[[k]]=as.data.frame(t(as.matrix(loslist_uci[[k]])))}
          else if (length(vartime)>1) {loslist_uci[[k]]=as.data.frame(as.matrix(loslist_uci[[k]]))}
          
          
         # loslist_uci[[k]]=as.data.frame(as.matrix(loslist_uci[[k]]))
        }
        
        for (k in states) {
          for (j in states) {
            colnames(loslist_uci[[k]])[j]=paste0("Los_",k,"_to_",j,"_uci")
          }
          colnames(loslist_uci[[k]])[nstates+1]  ="timevar"  
        } 
        
        loslist_uci
      }
      
      
      
    }
    
    
    pred[[g]]=  list(probs=pmlist,probs_lci=pmlist_lci,probs_uci=pmlist_uci,
         trans=trans, is.cumhaz=is.cumhaz,
         los=loslist,los_lci=loslist_lci,los_uci=loslist_uci)
  }
  
  
  
  rm(list=ls(pattern="^forFlex"))

  ########################################################################################################
  ######## Declare the list elements even if they stay empty ####
  
  pjson    =list()
  pjson_lci=list()
  pjson_uci=list()
  hjson=list()
  losjson=list()
  losjson_lci=list()
  losjson_uci=list()
  
  
  ########################################################################################################
  ####  Probabilities
  ########################################################################################################
  
  ##### Main estimates 
  
  pjson=list()
  
  if (length(covariates_list)<=1) {
    p=1
    for (k in 1:nstates) {
      for (i in 1:nstates) {
        
        pjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$probs[[k]][,i])))
        pjson[[p]]= assign(paste0("forFlex",names(pred[[1]]$probs[[k]])[i]),  pjson[[p]])
        p=p+1  
      }
    }
  }
  
  if (length(covariates_list)>1) {
    
    p=1
    
     if (process=="Markov") {
       for (k in 1:nstates) {
         for (i in 1:nstates) {
           
           tmplist=matrix(nrow =length(covariates_list) , ncol= length(vartime) )
           tmplist[1,]=as.vector(pred[[1]]$probs[[k]][,i])
           
           for (j in 2:length(covariates_list)) {
             tmplist[j,]=as.vector(pred[[j]]$probs[[k]][,i])  
           }
           pjson[[p]]=tmplist
           
           assign(paste0("forFlex",names(pred[[1]]$probs[[k]])[i]),  pjson[[p]])
           
           p=p+1  
           
         }
       }
     }
    
     if (process=="semiMarkov") {
       for (k in states) {
         for (i in states) {
          
          tmplist=matrix(nrow =length(covariates_list) , ncol= length(vartime)  )
          tmplist[1,]=as.vector(pred[[1]]$probs[[k]][,i])
          
          for (j in 2:length(covariates_list)) {
            tmplist[j,]=as.vector(pred[[j]]$probs[[k]][,i])  
          }
          pjson[[p]]=tmplist
          
          assign(paste0("forFlex",names(pred[[1]]$probs[[k]])[i]),  pjson[[p]])
          
          p=p+1  
          
        }
      }
    }
  }
  
  pvector=vector()
  pvector=ls(pattern = "forFlexP_._to_.$")  
  
  pvector=sub("forFlex","",ls(pattern = "^forFlexP_._to_.$") )
  
  pvector=stringi::stri_sort(pvector, numeric = TRUE)
  
  names(pjson)=pvector[1:(p-1)]
  
  pjson_new=list()
  
  for (d in 1:length(pjson)) {
    
    pjson_new[[d]]=list()
    
    for (c in 1:length(covariates_list)) {
      
      pjson_new[[d]][[c]]=as.vector(pjson[[d]][c,])
    }
   
  }
  names(pjson_new)=names(pjson)
  pjson=pjson_new
  
 # for (d in 1:length(pjson)) {
#   pjson[[d]]=tapply(pjson[[d]],rep(1:nrow(pjson[[d]]),each=ncol(pjson[[d]])),function(i)i)
#    pjson[[d]]=as.list(as.data.frame( pjson[[d]]))
#  }
  
  
  ##### Lci possibilities estimates
  
  if (ci.json==TRUE ) {
    
    
    pjson_lci=list()
    
    if (length(covariates_list)<=1) {
      p=1
      for (k in states) {
        for (i in states) {
          
          pjson_lci[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$probs_lci[[k]][,i])))
          assign(paste0("forFlex",names(pred[[1]]$probs_lci[[k]])[i]),  pjson_lci[[p]])
          p=p+1  
        }
      }
    }
    
    if (length(covariates_list)>1) {
      
      p=1
      for (k in states) {
        for (i in states) {
          
          tmplist=matrix(nrow =length(covariates_list) , ncol= length(vartime))
          tmplist[1,]=as.vector(pred[[1]]$probs_lci[[k]][,i])
          
          for (j in 2:length(covariates_list)) {
            tmplist[j,]=as.vector(pred[[j]]$probs_lci[[k]][,i])  
          }
          pjson_lci[[p]]=tmplist
          
          assign(paste0("forFlex",names(pred[[1]]$probs_lci[[k]])[i]),  pjson_lci[[p]])
          
          p=p+1  
          
          
        }
      }
    }
    
    pvector_lci=vector()
    pvector_lci=ls(pattern = "^forFlexP.*lci$")  
    
    pvector_lci=sub("forFlex","",ls(pattern = "^forFlexP.*lci$") )
    
    pvector_lci=stringi::stri_sort(pvector_lci, numeric = TRUE)
    
    names(pjson_lci)=pvector_lci[1:(p-1)]
    

    pjson_lci_new=list()
    
    for (d in 1:length(pjson_lci)) {
      
      pjson_lci_new[[d]]=list()
      
      for (c in 1:length(covariates_list)) {
        
        pjson_lci_new[[d]][[c]]=as.vector(pjson_lci[[d]][c,])
      }
      
    }
    names(pjson_lci_new)=names(pjson_lci)
    pjson_lci=pjson_lci_new 
    
    
    ##### Uci possibilities estimates
    
    
    pjson_uci=list()
    
    if (length(covariates_list)<=1) {
      p=1
      for (k in states) {
        for (i in states) {
          
          pjson_uci[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$probs_uci[[k]][,i])))
          assign(paste0("forFlex",names(pred[[1]]$probs_uci[[k]])[i]),  pjson_uci[[p]])
          p=p+1  
        }
      }
    }
    
    if (length(covariates_list)>1) {
      
      p=1
      for (k in states) {
        for (i in states) {
          
          tmplist=matrix(nrow =length(covariates_list) , ncol= length(vartime))
          tmplist[1,]=as.vector(pred[[1]]$probs_uci[[k]][,i])
          
          for (j in 2:length(covariates_list)) {
            tmplist[j,]=as.vector(pred[[j]]$probs_uci[[k]][,i])  
          }
          pjson_uci[[p]]=tmplist
          
          assign(paste0("forFlex",names(pred[[1]]$probs_uci[[k]])[i]),  pjson_uci[[p]])
          
          p=p+1  
          
          
        }
      }
    }
    
    pvector_uci=vector()
    pvector_uci=ls(pattern = "^forFlexP.*uci$")  
    
    pvector_uci=sub("forFlex","",ls(pattern = "^forFlexP.*uci$") )
    
    pvector_uci=stringi::stri_sort(pvector_uci, numeric = TRUE)
    
    names(pjson_uci)=pvector_uci[1:(p-1)]
    
    
    pjson_uci_new=list()
    
    for (d in 1:length(pjson_uci)) {
      
      pjson_uci_new[[d]]=list()
      
      for (c in 1:length(covariates_list)) {
        
        pjson_uci_new[[d]][[c]]=as.vector(pjson_uci[[d]][c,])
      }
      
    }
    names(pjson_uci_new)=names(pjson_uci)
    pjson_uci=pjson_uci_new 
    
    
    
  } 
  
  #################################################################################################################
  
  ### Transitions main estimates ###
  hjson=list()
  
  if (length(covariates_list)<=1) {
    
    for (k in 1:ntransitions) {
      
      tmplisth=matrix(nrow =length(covariates_list) , ncol= length(vartime))
      tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
      hjson[[1]]=tmplisth
      assign(paste0("forFlex",names(pred[[1]]$trans)[i]),  hjson[[1]])
    }
    
  }
    
  if (length(covariates_list)>1) {
    
    
    
    for (k in 1:ntransitions) {
      
      tmplisth=matrix(nrow =length(covariates_list) , ncol= length(vartime))
      
      tmplisth[1,]=as.vector(pred[[1]]$trans[,k])
      
      for (j in 2:length(covariates_list)) {
        
        tmplisth[j,]=as.vector(pred[[j]]$trans[,k])  
      }
      
      
      hjson[[k]]=tmplisth
      assign(paste0("forFlex",names(pred[[1]]$trans)[k]),  hjson[[k]])
    }
  }
  
  hvector=vector()
  hvector=ls(pattern = "^forFlexHaz_._to_.$")  
  
  hvector=sub("forFlex","",ls(pattern = "^forFlexHaz_._to_.$" ) )
  
  hvector=stringi::stri_sort(hvector, numeric = TRUE)
  
  names(hjson)=hvector
  
 # for (d in 1:length(hjson)) {
    
#    hjson[[d]]=tapply(hjson[[d]],rep(1:nrow(hjson[[d]]),each=ncol(hjson[[d]])),function(i)i)
#    hjson[[d]]=as.list(as.data.frame( hjson[[d]]))
#  }
  
  hjson_new=list()
  
  for (d in 1:length(hjson)) {
    
    hjson_new[[d]]=list()
    
    for (c in 1:length(covariates_list)) {
      
      hjson_new[[d]][[c]]=as.vector(hjson[[d]][c,])
    }
    
  }
  
  names(hjson_new)=names(hjson)
  hjson=hjson_new  

  
  
  if (totlos==TRUE) {
    
    ### Los Main estimates #####
    losjson=list()
    
    if (length(covariates_list)<=1) {
      p=1
      for (k in 1:nstates) {
        for (i in 1:nstates) {
          
          losjson[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$los[[k]][,i])))
          assign(paste0("forFlex",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
          p=p+1  
        }
      }
    }
    
    if (length(covariates_list)>1) {
      
      p=1
      for (k in states) {
        for (i in states) {
          
          tmplistlos=matrix(nrow =length(covariates_list) , ncol= length(vartime))
          tmplistlos[1,]=as.vector(pred[[1]]$los[[k]][,i])
          
          for (j in 2:length(covariates_list)) {
            tmplistlos[j,]=as.vector(pred[[j]]$los[[k]][,i])  
          }
          losjson[[p]]=tmplistlos
          
          assign(paste0("forFlex",names(pred[[1]]$los[[k]])[i]),  losjson[[p]])
          
          p=p+1  
        }
      }
    }
    
    losvector=vector()
    losvector=ls(pattern = "^forFlexLos.*$")  
    
    losvector=sub("forFlex","",ls(pattern = "^forFlexLos.*$") )
    
    losvector=stringi::stri_sort(losvector, numeric = TRUE)
    
    names(losjson)=losvector
    
    #    for (d in 1:length(losjson)) {
    
    #     losjson[[d]]=tapply(losjson[[d]],rep(1:nrow(losjson[[d]]),each=ncol(losjson[[d]])),function(i)i)
    #     losjson[[d]]=as.list(as.data.frame( losjson[[d]]))
    #   }
    
    
    losjson_new=list()
    
    for (d in 1:length(losjson)) {
      
      losjson_new[[d]]=list()
      
      for (c in 1:length(covariates_list)) {
        
        losjson_new[[d]][[c]]=as.vector(losjson[[d]][c,])
      }
      
    }
    names(losjson_new)=names(losjson)
    losjson=losjson_new  
    
    
  

    
    if (ci.json==TRUE ) {
      
      
      ### Los lci estimates #####
      
      losjson_lci=list()
      
      if (length(covariates_list)<=1) {
        p=1
        for (k in 1:nstates) {
          for (i in 1:nstates) {
            
            losjson_lci[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$los_lci[[k]][,i])))
            assign(paste0("forFlex",names(pred[[1]]$los_lci[[k]])[i]),  losjson_lci[[p]])
            p=p+1  
          }
        }
      }
      
      if (length(covariates_list)>1) {
        
        p=1
        for (k in states) {
          for (i in states) {
            
            tmplistlos=matrix(nrow =length(covariates_list) , ncol= length(vartime))
            tmplistlos[1,]=as.vector(pred[[1]]$los_lci[[k]][,i])
            
            for (j in 2:length(covariates_list)) {
              tmplistlos[j,]=as.vector(pred[[j]]$los_lci[[k]][,i])  
            }
            losjson_lci[[p]]=tmplistlos
            
            assign(paste0("forFlex",names(pred[[1]]$los_lci[[k]])[i]),  losjson_lci[[p]])
            
            p=p+1  
          }
        }
      }
      
      losvector_lci=vector()
      losvector_lci=ls(pattern = "^forFlexLos.*lci$")  
      
      losvector_lci=sub("forFlex","",ls(pattern = "^forFlexLos.*lci$") )
      
      losvector_lci=stringi::stri_sort(losvector_lci, numeric = TRUE)
      
      names(losjson_lci)=losvector_lci
      
      #    for (d in 1:length(losjson_lci)) {
      
      #     losjson_lci[[d]]=tapply(losjson_lci[[d]],rep(1:nrow(losjson_lci[[d]]),each=ncol(losjson_lci[[d]])),function(i)i)
      #     losjson_lci[[d]]=as.list(as.data.frame( losjson_lci[[d]]))
      #   }
      
      
      losjson_lci_new=list()
      
      for (d in 1:length(losjson_lci)) {
        
        losjson_lci_new[[d]]=list()
        
        for (c in 1:length(covariates_list)) {
          
          losjson_lci_new[[d]][[c]]=as.vector(losjson_lci[[d]][c,])
        }
        
      }
      names(losjson_lci_new)=names(losjson_lci)
      losjson_lci=losjson_lci_new  
      

      
      
      ### Los uci estimates #####
      
      losjson_uci=list()
      
      if (length(covariates_list)<=1) {
        p=1
        for (k in 1:nstates) {
          for (i in 1:nstates) {
            
            losjson_uci[[p]]= as.data.frame(rbind(as.vector(pred[[1]]$los_uci[[k]][,i])))
            assign(paste0("forFlex",names(pred[[1]]$los_uci[[k]])[i]),  losjson_uci[[p]])
            p=p+1  
          }
        }
      }
      
      if (length(covariates_list)>1) {
        
        p=1
        for (k in states) {
          for (i in states) {
            
            tmplistlos=matrix(nrow =length(covariates_list) , ncol= length(vartime))
            tmplistlos[1,]=as.vector(pred[[1]]$los_uci[[k]][,i])
            
            for (j in 2:length(covariates_list)) {
              tmplistlos[j,]=as.vector(pred[[j]]$los_uci[[k]][,i])  
            }
            losjson_uci[[p]]=tmplistlos
            
            assign(paste0("forFlex",names(pred[[1]]$los_uci[[k]])[i]),  losjson_uci[[p]])
            
            p=p+1  
          }
        }
      }
      
      losvector_uci=vector()
      losvector_uci=ls(pattern = "^forFlexLos.*uci$")  
      
      losvector_uci=sub("forFlex","",ls(pattern = "^forFlexLos.*uci$") )
      
      losvector_uci=stringi::stri_sort(losvector_uci, numeric = TRUE)
      
      names(losjson_uci)=losvector_uci
      
  #    for (d in 1:length(losjson_uci)) {
        
   #     losjson_uci[[d]]=tapply(losjson_uci[[d]],rep(1:nrow(losjson_uci[[d]]),each=ncol(losjson_uci[[d]])),function(i)i)
   #     losjson_uci[[d]]=as.list(as.data.frame( losjson_uci[[d]]))
   #   }
      
      
      losjson_uci_new=list()
      
      for (d in 1:length(losjson_uci)) {
        
        losjson_uci_new[[d]]=list()
        
        for (c in 1:length(covariates_list)) {
          
          losjson_uci_new[[d]][[c]]=as.vector(losjson_uci[[d]][c,])
        }
        
      }
      names(losjson_uci_new)=names(losjson_uci)
      losjson_uci=losjson_uci_new  
      
      
 ###cijson
    }
  ###totlos  
  }     
  
  
  
  ##########################################################################################
  
  ### Timevar###
  
  timejson=list()
  
  timejson[[1]]=vartime
  
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
  
  atlistjson[[1]]= as.vector(atlistmatrix)
  
  
  
  ######################################################################################################
  is.cumhaz=1

  
  final_list=list(timevar=timejson, Nats=length(covariates_list) ,
                  atlist=atlistjson, tmat=tmatlist,  Ntransitions= Ntransitions,is.cumhaz=is.cumhaz, 
                  pjson,pjson_lci,pjson_uci,
                  hjson,
                  losjson,losjson_lci,losjson_uci )
  
 # final_list$timevar=as.vector(final_list$timevar)
  
  final_unlist= unlist(final_list, recursive=FALSE)
  
 # exportJson <- toJSON(final_unlist, pretty = TRUE,force = TRUE, flatten=TRUE, na='string')
  exportJson <- toJSON(final_unlist,force = TRUE, flatten=TRUE)
  
  
  write(exportJson, paste0(jsonpath,"/", name ) )
  
  rm(list=ls(pattern="^forFlex"))
  
  final_list
  
  }



