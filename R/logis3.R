#' @title R optmisation of a three-parameter logistic particle size distribution function
#' @description The function optimises 3-parameter logistic  particle size distribution function
#' of Skaggs et al(2001)
#' @author George Owusu
#' @inheritParams lass3
#' @inheritParams gompertz
#' @return 
#' \itemize{
#' \item{PD:} {  predicted fraction by mass of particles passing a particular diameter}
#' \item{sand:} { percentage of sand}
#' \item{clay:} { percentage of clay}
#' \item{silt:} { percentage of silt }
#' \item{a:} { optimised a parameter }
#' \item{b:} { optimised b parameter }
#' \item{e:} { optimised e parameter }
#' } 
#' @family Logistic functions
#' @seealso \code{\link{zhuang3}}, \code{\link{zhuang4}},\code{\link{lass3}},\code{\link{logarithmic}},
#' \code{\link{logistic2}}, \code{\link{logistic3}},\code{\link{logistic4}},\code{\link{fredlund3}},\code{\link{fredlund4}},
#'  \code{\link{jaky}},\code{\link{andersson}},\code{\link{gompertz}},
#' \code{\link{haverkamp}}
#' @export
#' @references 
#' Skaggs, T. H., Arya, L. M., Shouse, P. J., & Mohanty, B. P. (2001). 
#' Estimating Particle-Size Distribution from Limited Soil Texture Data. 
#' Soil Sci Soc Am J., 65, 1038-1044. 
#' @examples
#' data=read.csv(system.file("ext","sys","soil2.csv",package="vadose"))
#' single<- subset(data, ID=="30B20_1")
#' mod=logistic3 (data=single,p="sand",D="D",fr="Sand.")
#' plot(mod)
#' \dontrun{
#' #group simulation
#' mod2=logistic3 (data=data,D="D",fr="FractionSand",group="ID")
#' mod2$coef
#' }
#' @rdname logistic3

#generic function
logistic3<-function(data=NULL,D=NULL,fr=NULL,p=NULL,group=NULL,a=1,b=1,e=1) UseMethod ("logistic3")
#' @export
#' @rdname logistic3
#default function
logistic3.default<-function(data=NULL,D=NULL,fr=NULL,p=NULL,group=NULL,a=1,b=1,e=1)
{

  #stop warning from displaying
  options(warn=-1)
  PSDF=NULL
  #remove zeros
  #print("yes")
  if(!is.null(data)&&is.null(D)&&is.null(fr)){
    D=data$D
    fr=data$fr
    if(is.null(fr)){
      fr=data$F
    }
    if(is.null(fr)){
      fr=data$FR
    }
    if(is.null(fr)){
      fr=data$particle_fraction
    }
    if(is.null(fr)){
      fr=data$fraction
    }
    
    if(is.null(D)){
      D=data$d
    }
    if(is.null(D)){
      D=data$diameter
    }
    if(is.null(D)){
      D=data$Diameter
    }
    if(is.null(D)){
      D=data$particle_size
    }
    data=NULL
  }
  if(is.null(data)&&length(D)>1&&length(fr)>1){
    data=as.data.frame(cbind(D,fr))
    D=colnames(data)[1]
    fr=colnames(data)[2]
    #print(cbind(D,fr))
  }
  
  #predict fraction based on input
  if(is.null(data))
  {
    if(!is.null(fr))
    {
      return(print("please provide the diameter,D rather"))
    }
    fr=(e/(1+a*exp(-b*D)))
    return(print(list(fraction=fr)))
  }
  
  if(max(data[[D]])>10){
    data[[D]]=data[[D]]/1000
  }
  
  data[data==0]<-0.0001
  #group function
  if(!is.null(group)){
    aggdata =row.names(table(data[group]))
    #create group data frame###################################
    addoutput=data.frame(groupid=factor(),a=numeric(),b=numeric(),e=numeric(),porosity=numeric(),r2=numeric(),
                         RMSE=numeric(),NRMSE=numeric(),F_value=numeric(),df=numeric(),p_value=numeric(),AIC=numeric(),BIC=numeric(),r2_adj=numeric(),RMSE_adj=numeric(),NRMSE_adj=numeric(),nashSutcliffe=numeric(),texture=numeric())
    i=1
    while(i<=length(aggdata)){
      #print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      p2=p
      if(length(p)>1){
        p2=p[i]
      }
      mod1=logistic3(single,p=p2)
      gof1=gof(mod1)
      ######################################
      a=coef(mod1$logistic3)[[1]]
      b=coef(mod1$logistic3)[[2]]
      e=coef(mod1$logistic3)[[3]]
      
      
      txt1=texture(mod1)
      sand=mod1$sand
      silt=mod1$silt
      clay=mod1$clay
      #K=mod1$K
      #CU=mod1$CU
      
      
      #######################################################
      addoutput=rbind(addoutput,data.frame(groupid=aggdata[i],a=a,b=b,e=e,porosity=mod1$p,r2=gof1$r2,
                                           RMSE=gof1$RMSE,NRMSE=gof1$NRMSE,F_value=gof1$F_value,df=gof1$df,p_value=gof1$p_value,AIC=gof1$AIC,BIC=gof1$BIC,r2_adj=gof1$r2_adj,RMSE_adj=gof1$RMSE_adj,NRMSE_adj=gof1$NRMSE_adj,nashSutcliffe=gof1$nashSutcliffe,sand=sand,silt=silt,clay=clay,texture=txt1))
      i=i+1
    }
    factor=list(coef=addoutput)
  }
  else{
    
    if(is.null(p)){
      mod=logistic3 (data=data,p="loam",D=D,fr=fr)
      p=texture(mod)
    }
    
    #use published data
    if(!is.numeric(p)&&!is.null(p)){
      soil=ksat(para="soil",model=p)
      p=soil$para$p
      if(is.na(p)){
        p=soil$para$ths
      }
    }
    if(is.null(data)||is.null(D)||is.null(fr)){
      return(print("Data, diameter, D, and fraction, fr, of the particles are needed please"))
      break
    }
    if(max(data[[fr]])>2){
      data[[fr]]=data[[fr]]/100
    }
    ##############################################################
    PDf=paste(fr,"~(e/(1+a*exp(-b*",D,")))")
    
    ones=c(a=a,b=b,e=e)
    logistic3<- nlxb(PDf, start = ones, trace = FALSE, data = data)
    
    a=coef(logistic3)[1]
    b=coef(logistic3)[2]
    e=coef(logistic3)[3]
    
    
    
    PD=(e/(1+a*exp(-b*data[[D]])))
    x=seq(min(data[[D]]),max(data[[D]]),0.01)
    y=(e/(1+a*exp(-b*x)))
    
    
    i=min(data[[D]])
    sand=0
    while (i<=max(data[[D]])){
      
      PD2=(e/(1+a*exp(-b*i)))
      
      
      if(round(PD2[[1]],3)==0.1||round(PD2[[1]],2)==0.1||round(PD2[[1]],1)==0.1){
        D10=i
      }
      if(round(PD2[[1]],3)==0.6||round(PD2[[1]],2)==0.6||round(PD2[[1]],1)==0.6){
        D60=i
      }
      
      
      i=i+0.001
    }
    #CU=D60/D10
    # Hazen approximation
    #K=0.01*(D10^2)
    
    PSDF=TRUE
    
    
    
    #texture#########################################
    two=(e/(1+a*exp(-b*2)))*100
    silt=(e/(1+a*exp(-b*0.05)))*100
    clay=(e/(1+a*exp(-b*0.002)))*100
    sand=100-silt 
    silt=silt-clay
    clay=clay
   factor= list(PSDF=length(coef(logistic3)),x=x,y=y,logistic3=logistic3,predict=PD,D=data[[D]],a=a,b=b,fr=data[[fr]],p=p,sand=sand[[1]],silt=silt[[1]], clay=clay[[1]])
  }
factor$call<-match.call()

class(factor)<-"logistic3"
factor
}
#' @export
#' @rdname logistic3
#plot function
plot.logistic3<-function(x,main="Particle Size Distribution Function",xlab="Particle Diameter",ylab="Fraction",...)
{
  object=x
mod1=object
if(mod1$PSDF==TRUE||!is.null(mod1$PSDF==TRUE)){
plot(mod1$D,mod1$fr,xlab=xlab,ylab=ylab)
lines(mod1$x,mod1$y)
mtext(paste("R^2=",round(cor(mod1$predict,mod1$fr)^2,3)))
title(main)
}
}
#' @export
#' @rdname logistic3
coef.logistic3<-function(object,...)
{
x<-object$logistic3
coef=(coef(x))
return(coef)
}
#' @export
#' @rdname logistic3
#####################################################################
predict.logistic3<-function(object,D,...)
{
x<-object$logistic3

coef=(coef(x))
a=coef(x)[[1]]
b=coef(x)[[2]]
e=coef(x)[[3]]


PD2=(e/(1+a*exp(-b*D)))
#print(PD2)
return(PD2)
}


