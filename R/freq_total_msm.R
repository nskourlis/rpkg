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

#freq_func_total_msm(data_msm=data, id_msm= id,  timevar=seq(1,10,1), scale_msm=scale, Ntransitions=ntransitions) 
