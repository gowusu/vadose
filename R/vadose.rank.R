#' @title Ranking and  model performance assessment of vadose model object
#' 
#' @description This function ranks and boxplots \code{\link{BEST}},  \code{\link{OFEST}} or
#' \code{\link{lass3}} goodness of fit objects. Various models in  \code{\link{BEST}} 
#' and  \code{\link{OFEST}} are ranked by either median or mean or max or any of the fun objects. 
#' 
#'
#' @param data dataframe of goodness of fit indicators
#' @param var list of the indicators or variables that should be used to rank the models
#' @param model character. The column name indicating the models that are to be ranked. 
#' @param fun Either median or mean or max or min or any of the fun objects
#' @param xlab The x label of the boxplot
#' @param ylab The y label of the boxplot
#' @param main The title of the boxplot
#' @param ylim The y limits of the boxplot
#' @param las The orientation  of x labels
#' @param opar The par attributes
#' @param myrank a redundant variable for internal use
#' @author George Owusu
#' 
#'
#' @return rank and fun.summary indicating the aggregation of the data by fun
#' @export
#'
#' @examples
#' data(offinbest) #infiltration data
#' data(offinpsd) #Particle size Distribution data
#' data(htheta) #measured h and theta
#' PSD=offinpsd
#' Figure2=group.BEST(data = offinbest, group="TownName",PSD=PSD,layout=c(3,4),hlog=TRUE,plot=FALSE)
#' Figure3=group.OFEST(data = offinbest, hg="best",group="TownName",PSD=PSD,layout=c(3,4),hlog=TRUE,type = "nonlinear",plot=FALSE)
#' data=rbind(Figure2$statistics2,Figure3$statistics2)
#' ranking=vadose.rank(data)
#' ranking$rank
#' 
#' MBE=vadose.rank(data,var="MBE",ylim=c(-1,2),ylab=NULL,main=NULL)
#' MPE=ranking=vadose.rank(data,var="MPE",ylim=c(-20,30),ylab=NULL,main=NULL)
vadose.rank=function(data,var=c("r2","d","E","RMSE","ARE","MBE"),
  model="model",fun="median",xlab=NULL, ylab=c("R-square","Index of Agreement","Nash-Sutcliffe","RMSE","Average Relative Error","Mean Absolute Error")
  ,main=c("R-square","Index of Agreement","Nash-Sutcliffe","RMSE","Average Relative Error","Mean Bias Error"),
  ylim=list(c(0.99,1),c(0.98,1),c(0.9,1),c(0,2.5),c(0,30),c(-0.2,0.2)),
  las=2,opar=par(mfrow=c(2,3),mar=c(7,5,2,1)),myrank=NULL ){
  if(!is.null(data$statistics2)){
    data=data$statistics2
  }
  #myrank=NULL
  if(length(var)>1){
    fun.summary=NA
    opar
    val.length=length(var)
    for (i in 1:val.length) {
      #check input model
      if(length(model)==val.length){
        modeli=model[[i]] 
      }else{
        modeli=model 
      }
      
      if(length(fun)==val.length){
        funi=fun[[i]] 
      }else{
        funi=fun 
      }
      
      if(length(xlab)==val.length){
        xlabi=xlab[[i]] 
      }else{
        xlabi=xlab 
      }
      if(length(ylab)==val.length){
        ylabi=ylab[[i]] 
      }else{
        ylabi=ylab 
      }
      if(length(main)==val.length){
        maini=main[[i]] 
      }else{
        maini=main 
      }
      if(length(ylim)==val.length){
        ylimi=ylim[[i]] 
      }else{
        ylimi=ylim 
      }
      if(length(las)==val.length){
        lasi=las[[i]] 
      }else{
        lasi=las 
      }
      myrank=vadose.rank(data,var=var[i],myrank=myrank,model=modeli,
                         fun=funi,xlab=xlabi, ylab=ylabi,main=maini,ylim=ylimi,las=lasi) 
    }
    par(opar)
    rank<-sort(rowMeans(myrank))
    rank2<- data.frame(t(rank)[1,])
    row.names(rank2)<-names(rank)
    names(rank2)<-"rank"
    myrank=merge(myrank,rank2,by="row.names") 
    #row.names(myrank)=myrank[1]
    #aggdata=data[c(var,model)]
    #aggdata[[model]]=row.names(aggdata)
    agg<-vadose.tryCatch(aggregate(. ~ model, data = data[c(var,"model")], median))$value
    if(class(agg)[1]!="simpleError"){
      row.names(agg)<-agg[[1]]
      agg[[1]]=NULL
    }
    #agg<-agg[var]
    #myrank=merge(myrank,agg[var],by="row.names") 
    myrank=myrank[with(myrank, order(rank)), ]
    row.names(myrank)=myrank[[1]]
    myrank[[1]]=NULL
    list(rank=myrank,fun=fun,fun.summary=agg)
  }else{
    
    if(var=="MPE"||var=="MBE"){
      data[[model]] <-reorder(data[[model]], abs(data[[var]]), fun)
    }else{
      data[[model]] <-reorder(data[[model]], data[[var]], fun)
    }
    #if(!is.null(plot)){
    boxplot(data[[var]]~data[[model]],	xlab=xlab, ylab=ylab,main=main,ylim=ylim,las=las)
    #}
    this=levels(unique(data[[model]]))
    this2=data.frame(this,which(this>0))
    row.names(this2)=this2[[1]]
    this2[[1]]=NULL
    colnames(this2)=var
    if(var=="r2"||var=="r2_adj"||var=="Nash-Sutcliffe"||var=="d"||var=="E"||var=="F_value"||var=="d"){
      this2[[1]]=sort(this2[[1]],decreasing = TRUE) 
    }
    
    #print(this2)
    if(!is.null(myrank)){
      #print(myrank)
      this2=merge(myrank,this2,by="row.names") 
      row.names(this2)=this2[[1]]
      this2[[1]]=NULL
      #this2[name]=this2[row.names(this2)==row.names(myrank),]
    }else{
      this2
    }
    this2
  }
}