


#' Sensitivity analysis by factor
#' Calculates sobol indices along scanning one parameter
#' @param X1 data.frame of randomized model parameters for the sensitivity analysis
#' @param X2 data.frame of randomized model parameters for the sensitivity analysis
#' @param fixpar one row data.frame with required parameters not to be varied in the sensitivity analysis
#' @param scanvalues vector of values to be calculated for the gradient factor
#' @param scanpar character. Name of the Gradient factor
#' @param outpar name of the response variable (default "D18Olw")
#'
#' @return a list with two data frames (1) the estimations of the Sobol' first-order indices (2) the estimations of the Sobol' total sensitivity indices.
#' @export
#' @examples
#'
#' \dontrun{
#' parameters = get_parameter_definition ()
#' }

##********************* test importance of parameters to changes in air temperature

scan_sensitivity<-function(X1,X2,fixpar,scanvalues,scanpar,outpar="D18O_lw",element="O"){
  X1<-X1[,!names(X1) %in% c(names (fixpar),scanpar)]
  X2<-X2[,!names(X2) %in% c(names (fixpar),scanpar)]
  for (i in scanvalues){
    print (i)
    fixpar[scanpar]    <- i
    x <- sensitivity::sobol2007(model =plantIso_model, X1 = X1, X2 = X2, nboot = 1000, output = outpar,addpar=fixpar,element=element)
    if (i==scanvalues[1]){
      ds<-x$S$original
      dt<-x$T$original
    }else {
      ds<-cbind(ds,x$S$original)
      dt<-cbind(dt,x$T$original)
    }
  }
  rownames(ds)<-rownames(x$S)
  colnames(ds)<-scanvalues
  rownames(dt)<-rownames(x$T)
  colnames(dt)<-scanvalues
  output<-list(ds,dt)
  names(output)<-c('ds','dt')
  return (output)
}



#' Plot the relative sensitivity of each factor along a gradient.
#' This function plots the output of scan_sensitivity()
#' @param ds one of the objects in the list created by scan_sensitivity()
#' @param xlab data.frame of randomized model parameters for the sensitivity analysis
#' @param ylab data.frame of randomized model parameters for the sensitivity analysis
#' @param main data.frame of randomized model parameters for the sensitivity analysis
#' @param fcolors colors for the different factors (default colors as defined in #http://colorbrewer2.org/#type=diverging&scheme=Spectral&n=6)
#' @param parcol index of the gradient color, will be hold out from the plot (default=1)
#' @param legend boolean. plot legend (default=FALSE)
#' @importFrom  graphics plot polygon
#' @export
#' @examples
#' \dontrun{
#' plot_sensitivity(ds,scanpar,'normalized first-order indices','SA PMCG for D18O_lw',parcol = 1)
#' }


plot_sensitivity<-function(ds,xlab,ylab,main,fcolors=c('#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd'),parcol=1,legend=FALSE){
  ds[which(ds<0)]<-0
  scanvalues <- as.numeric(colnames(ds))
  csds<-colSums(ds)
  nds<-ds
  for (i in 1:nrow(ds)) {nds[i,]<-ds[i,]/csds}
  for (i in 2:nrow(ds)) {nds[i,]<-nds[i,]+nds[i-1,]}
  nds<-rbind(nds,(rep(0,length(scanvalues))))
  nds<-nds[c(nrow(nds),2:nrow(nds)-1),]
  fcolors<-fcolors[-parcol]
  plot(NA,xlim=range(scanvalues),ylim=c(0,1.05),bty='l',xaxs='i',yaxs='i',las=1,xlab=xlab,ylab=ylab,main=main)
  py <- rep(0,ncol(nds))
  for (i in 2:nrow(nds)-1){
    y<-c(nds[i,],rev(nds[i+1,]))
    polygon(c(scanvalues,rev(scanvalues)),y,border =NA,col = fcolors[i])
  }
  if(legend)legend(min(scanvalues),0.7,pch = NA,rev(rownames(ds)),xpd = TRUE,bty='n',fill=rev(fcolors[1:nrow(ds)]))
}

