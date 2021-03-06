#' @title R optmisation of Fredlund et al (2000) four-parameter particle size distribution function
#' @description The function optimises 4 parameters of Fredlund et al (2000) particle size distribution function
#' @author George Owusu
#' @inheritParams BEST
#' @param df an initial value of 'df' parameter to be optimised
#' @inheritParams lass3
#' @inheritParams fredlund3
#' @family Fredlund functions
#' @seealso \code{\link{zhuang3}}, \code{\link{zhuang4}},\code{\link{lass3}},\code{\link{logarithmic}},
#' \code{\link{logistic2}}, \code{\link{logistic3}},\code{\link{logistic4}},\code{\link{fredlund3}},
#'  \code{\link{jaky}},\code{\link{andersson}},\code{\link{gompertz}},
#' \code{\link{haverkamp}}
#' @return 
#' \itemize{
#' \item{PD:} {  predicted fraction by mass of particles passing a particular diameter}
#' \item{sand:} { percentage of sand}
#' \item{clay:} { percentage of clay}
#' \item{silt:} { percentage of silt }
#' \item{a:} { optimised a parameter }
#' \item{m:} { optimised m parameter }
#' \item{n:} { optimised n parameter }
#' \item{df:} { optimised df parameter }

#' } 
#' @export
#' @references 
#' Fredlund, M. D., Fredlund, D. G., & Wilson, G. W. (2000). 
#' An equation to represent grain-size distribution. Can Geotech J, 65, 638-648. 
#' @examples
#' data=read.csv(system.file("ext","sys","soil2.csv",package="vadose"))
#' single<- subset(data, ID=="30B20_1")
#' mod=fredlund4 (data=single,p="sand",D="D",fr="Sand.")
#' plot(mod)
#' 
#' #group simulation
#' mod2=fredlund4 (data=data,D="D",fr="FractionSand",group="ID")
#' mod2$coef
#' #generic function
#generic function
#' @rdname fredlund4
#' 
#generic function
fredlund4<-function(data=NULL,D=NULL,fr=NULL,p=NULL,group=NULL,a=1,m=1,n=1,df=1) UseMethod ("fredlund4")
#' @export
#' @rdname fredlund4
#default function
fredlund4.default<-function(data=NULL,D=NULL,fr=NULL,p=NULL,group=NULL,a=1,m=1,n=1,df=1)
{

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
    dmin=0.001
    fr=((1/(log(exp(1)+(a/D)^n))^m)*(1-(log(1+(df/D))/log(1+(df/dmin)))^7))
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
    addoutput=data.frame(groupid=factor(),a=numeric(),m=numeric(),n=numeric(),df1=numeric(),porosity=numeric(),r2=numeric(),
                         RMSE=numeric(),NRMSE=numeric(),F_value=numeric(),df=numeric(),p_value=numeric(),AIC=numeric(),BIC=numeric(),r2_adj=numeric(),RMSE_adj=numeric(),NRMSE_adj=numeric(),nashSutcliffe=numeric(),texture=numeric())
    i=1
    while(i<=length(aggdata)){
      #print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      p2=p
      if(length(p)>1){
        p2=p[i]
      }
      mod1=fredlund4(single,p=p2)
      gof1=gof(mod1)
      ######################################
      a=coef(mod1$fredlund4)[[1]]
      m=coef(mod1$fredlund4)[[2]]
      n=coef(mod1$fredlund4)[[3]]
      df1=coef(mod1$fredlund4)[[4]]
      
      
      txt1=texture(mod1)
      sand=mod1$sand
      silt=mod1$silt
      clay=mod1$clay
      #K=mod1$K
      #CU=mod1$CU
      
      
      #######################################################
      addoutput=rbind(addoutput,data.frame(groupid=aggdata[i],a=a,m=m,n=n,df1=df1,porosity=mod1$p,r2=gof1$r2,
                                           RMSE=gof1$RMSE,NRMSE=gof1$NRMSE,F_value=gof1$F_value,df=gof1$df,p_value=gof1$p_value,AIC=gof1$AIC,BIC=gof1$BIC,r2_adj=gof1$r2_adj,RMSE_adj=gof1$RMSE_adj,NRMSE_adj=gof1$NRMSE_adj,nashSutcliffe=gof1$nashSutcliffe,sand=sand,silt=silt,clay=clay,texture=txt1))
      i=i+1
    }
    factor=list(coef=addoutput)
  }
  else{
    
    if(is.null(p)){
      mod=fredlund4 (data=data,p="loam",D=D,fr=fr)
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
    dmin=min(data[[D]])
    PDf=paste(fr,"~((1/(log(exp(1)+(a/",D,")^n))^m)*(1-(log(1+(df/",D,"))/log(1+(df/",dmin,")))^7))")
    ones=c(a=a,m=m,n=n,df=df)
    fredlund4<- nlxb(PDf, start = ones, trace = FALSE, data = data)
    #print(coef(fredlund4))
    
    a=coef(fredlund4)[1]
    m=coef(fredlund4)[2]
    n=coef(fredlund4)[3]
    df=coef(fredlund4)[4]
    
    
    PD=((1/(log(exp(1)+(a/data[[D]])^n))^m)*(1-(log(1+(df/data[[D]]))/log(1+(df/dmin)))^7))
    x=seq(min(data[[D]]),max(data[[D]]),0.01)
    y=((1/(log(exp(1)+(a/x)^n))^m)*(1-(log(1+(df/x))/log(1+(df/dmin)))^7))
    
    
    
    i<-min(data[[D]])
    sand=0
    while (i<=max(data[[D]])){
      
      PD2=((1/(log(exp(1)+(a/i)^n))^m)*(1-(log(1+(df/i))/log(1+(df/dmin)))^7))
      
      
      
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
    two=((1/(log(exp(1)+(a/2)^n))^m)*(1-(log(1+(df/2))/log(1+(df/dmin)))^7))*100
    
    silt=((1/(log(exp(1)+(a/0.05)^n))^m)*(1-(log(1+(df/0.05))/log(1+(df/dmin)))^7))*100
    
    clay=((1/(log(exp(1)+(a/0.002)^n))^m)*(1-(log(1+(df/0.002))/log(1+(df/dmin)))^7))*100
    
    sand=100-silt 
    silt=silt-clay
    clay=clay
    #####################################
    factor<-list(PSDF=length(coef(fredlund4)),x=x,y=y,fredlund4=fredlund4,predict=PD,D=data[[D]],a=a,m=m,n=n,df=df,dmin=dmin,fr=data[[fr]],p=p,sand=sand[[1]],silt=silt[[1]], clay=clay[[1]])
  }

factor$call<-match.call()

class(factor)<-"fredlund4"
factor
}
#' @export
#' @rdname fredlund4

#plot function
plot.fredlund4<-function(x,main="Particle Size Distribution Function",xlab="Particle Diameter",ylab="Fraction",...)
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
#' @rdname fredlund4
coef.fredlund4<-function(object,...)
{
x<-object$fredlund4
coef=(coef(x))
return(coef)
}
#' @export
#' @rdname fredlund4
#####################################################################
predict.fredlund4<-function(object,D,...)
{
x<-object$fredlund4

a=coef(x)[[1]]
m=coef(x)[[2]]
n=coef(x)[[3]]
df=coef(x)[[4]]

dmin=object$dmin
PD2=((1/(log(exp(1)+(a/D)^n))^m)*(1-(log(1+(df/D))/log(1+(df/dmin)))^7))

#print(PD2)
return(PD2)
}

