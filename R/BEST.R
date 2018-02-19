#' @title R estimation of Beerkan estimation of soil transfer (BEST) 
#' @description These functions compute the various versions of BEST. 
#' \code{BEST} computes all the 3 versions of BEST:Lassabatere etal(2006), Yilmaz et al(2010) and Bagarello et al (2014).
#' \code{BESTlass} implements Lassabatere algorithm
#' \code{BESTyb} implements Lassabatere algorithms of Yilmaz and Bagarello et al. et al.
#' 
#' @param data dataframe. It can contain data with column names of "time" and "I"
#' @param legend A legend of the plot
#' @param time character or numeric. The name of time variable in the dataframe.
#' If the "data" parameter contains "time", this will be ignored.
#' The unit must be in seconds.
#' @param I character or numeric. The name of cumulative infiltration variable 
#' in the dataframe. If the "data" parameter contains "I", this will be ignored.
#' The unit must be in millimetres [mm].
#' @param S Initial value of soil sorptivity [L T^(-1/2)]. The default is 0.1.  
#' This will be optimise
#' @param n numeric. A shaping parameter for water retention curve. This can be calibrated.
#' If the parameter "PSD" is not set NULL \code{\link{lass3}} will be automatically used to estimate it.
#' @param m character. The water retention curve condition. 
#' It takes either "b" for Burdine condition or "m" for Mualem condition. 
#' @param y numeric. coefficient of BEST A parameter, commonly set at 0.75[-]
#' This can be calibrated.
#' @param b numeric. coefficient of BEST B parameter, commonly set at 0.6[-].
#' This can be calibrated.
#' @param r numeric. Radius of infiltrometer ring [mm]. The default is 75 mm
#' @param pb numeric. Bulk density[g cm^-3]. This can be calibrated.
#' @param tho numeric. initial volumetric soil water content. This can be calibrated.
#' @param thr numeric. Residual volumetric soil water content. This can be calibrated.
#' @param theta The  measured soil water content[m3m-3] at the modelled pressure levels.
#' It is used to check the fitness of the model.
#' @param PSD dataframe. Particle Size Distribution. A dataframe of two columns.
#' The first column should be labeled "D" in the range of 0.001 to 2 mm.
#' The second column should be named as "fr" i.e.' the fraction of Diameter. 
#' The range should be 0-1.
#' @param h list. The range of metric head desired or measured [mm].
#' @param steady The number of data points that should be used for steady state.
#' This can be observed during infiltration measurement
#' @param init numeric. The number of data points that should be initially used for optimisation of
#' Best equation. 
#' @param ths numeric. The saturated soil water content [m3/m3]. 
#' @param x a return object of the function.
#' @param main Title of the plot
#' @param xlab x label of the plot
#' @param ylab y label of the plot
#' @param ylab2 y label of the second plot (K)
#' @param hlog TRUE or FALSE. Whether metric potential should be log transformed
#' @param klog TRUE or FALSE. Whether hydraulic conductivity should be log transformed
#' @param col Color of the plot
#' @param units Units of the plot
#' @param mfrow The graphical layout of the plots
#' @param type The type of plot. It takes  "all" or "h"or "hk" or "k"
#' @param obs Observed data
#' @param est Estimated or modelled data
#' @param model A calibration model parameter. It takes "all" for \code{{BEST}} "Lassabatere" for 
#' \code{{BESTlass}} and "Yilmaz" for \code{{BESTyb}}, " Bagarello" \code{{BESTyb}} 
#' @param P Number of Parameters
#' @param plot whether a plot should be performed.
#' @param ... Any other graphical parameter
#' @inheritParams lass3
#' @inheritParams vg
#' @author George Owusu
#' @references 
#' \itemize{
#' \item{}{Lassabatere, L., Angulo-Jaramillo, R., Soria Ugalde, J. M., Cuenca, R., 
#' Braud, I., & Haverkamp, R. (2006). Beerkan Estimation of Soil Transfer 
#' Parameters through Infiltration Experiments-BEST. 
#' Soil Sci. Soc. Am. J., 70(2), 521-532. doi: 10.2136/sssaj2005.0026}
#' \item{}{Yilmaz, D., Lassabatere, L., Angulo-Jaramillo, R., Deneele, D., & Legret, M. (2010). 
#' Hydrodynamic Characterization of Basic Oxygen Furnace Slag through an 
#' Adapted BEST Method. Vadose Zone J., 9(1), 107-116. doi: 10.2136/vzj2009.0039}
#'  \item{}{Aiello, R., Bagarello, V., Barbagallo, S., Consoli, S., Di Prima, S., Giordano, 
#'  G., & Iovino, M. (2014). An assessment of the Beerkan method for determining the 
#'  hydraulic properties of a sandy loam soil. Geoderma, 235-236(0), 300-307. 
#'  doi: http://dx.doi.org/10.1016/j.geoderma.2014.07.024}
#' \item{}{Bagarello, V., Di Prima, S., & Iovino, M. (2014). 
#' Comparing Alternative Algorithms to Analyze the Beerkan Infiltration Experiment. 
#' Soil Sci. Soc. Am. J., 78(3), 724-736. doi: 10.2136/sssaj2013.06.0231}
#' }
#' @details 
#' The function can perform normal BEST computation with \code{BEST}.
#' 
#' @import nlmrt
#' @return The return parameters can be assessed with \code{coef.BEST}. 
#' The predicted soil water can be assessed with \code{predict.BEST}. 
#' In all the function returns BEST parameters A, B, C, Ks, S, slope, intercept, mod_theta, hg.
#' Because \code{BEST} combines \code{BESTlass}  and \code{BESTyb} 'l', 'y', 
#' and 'b' can be used to access Lassabatere etal(2006), Yilmaz et al(2010) 
#' and Bagarello et al (2014) parameters, respectively. 
#' For  the definitions of the output of PSD see \code{\link{lass3}}.
#' For the definitions of the output of goodness of fit tests see \code{\link{gof}}
#' \itemize{
#' \item{A:} {  constant[L^-1]}
#' \item{B:} { constant[L^-1]}
#' \item{C: } { constant[L^-1]}
#' \item{Ks,Ksl,Ksy,Ksb:} { saturated hydraulic conductivity for various algorithms [LT^-1] }
#' \item{S,Sl,Sy,Sb:} { Sorptivity for various algorithms [LT^-0.5]}
#' \item{slope,slopel,slopey:} { slope of steady state data for variousalgorithms}
#' \item{intercept,interceptl,intercepty:} { intercept of steady state data for variousalgorithms}
#' \item{mod_theta,mod_thetal,mod_thetay,mod_thetab:} { The estimated soil water content at 
#' different metric potentials for variousalgorithms }
#' \item{theta:} { The observed soil water content at different metric potentials  }
#' \item{K,Kl,Ky,Kb:} { The unsaturated hydraulic conductivity [LT^-1] }
#' \item{Kr,Krl,Kry,Krb:} { The ratio of K and Ks }
#' \item{hg,hgl,hgy,hgb:} { bubbling capillary pressure for various algorithms [L] }
#' \item{init:} { the number of data points that should be initially used for 
#' optimisation of Best equation [T]}
#' \item{q:} { infiltration rates [LT^-1] }
#' } 
#' @export
#'
#' @examples
#' \dontrun{
#' psd=read.csv(system.file("ext","sys","psd.csv",package="vadose"))
#' hf=read.csv(system.file("ext","sys","h.csv",package="vadose"))
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' mod1<-BEST(data=data,time="time",I="I",h=hf$h,PSD=psd)
#' 
#' pred=predict(mod1)
#' coefficients=coef(mod1)
#' 
#' # with measured theta and h
#' mod<-BEST(data=data,time="time",I="I",h=hf$h,theta=mod1$mod_thetab,PSD=psd,y=0.7,b=0.6)
#' print(gof.BEST(mod))
#' plot(mod,hlog="yes",klog=TRUE)
#' par(mfrow=c(1,2))
#' plot(mod,type="h")
#' plot(mod,type="psd")
#' gof1=gof.BEST(mod)
#' print(gof1)
#' 
#' #calibration
#' cal.mod<-cal.BEST(data=data,time="time",I="I",h=hf$h,theta=mod1$mod_thetab,PSD=psd,steady=3,
#' init=5,thr=0,ths=NULL,tho=0.162)
#' #calibrate all, set the model to 'all' Lassabatere etal(2006), 
#' Yilmaz et al(2010) and Bagarello et al (2014)
#' cal.modall<-cal.BEST(data=data,time="time",I="I",h=hf$h,theta=mod1$mod_thetab,PSD=psd,steady=3,
#' init=5,thr=0,ths=NULL,tho=0.162,model="all")
#' #get the parameters
#' parameters=coef.BEST(cal.modall)
#' S=parameters$S
#' A=parameters$A
#' 
#' }
#generic function
BEST<-function(data=NULL,time,I,S=0.1,n=NULL,m="b",y=0.75,b=0.6,r=75,pb=1.2,tho=0.169,
               thr=0,PSD=NULL,h=NULL,theta=NULL,steady=3,init=5,ths=NULL) UseMethod ("BEST")

#' @rdname BEST
#' @export
#default function
BEST.default<-function(data=NULL,time,I,S=0.1,n=NULL,m="b",y=0.75,b=0.6,r=75,pb=1.2,tho=0.169,
                       thr=0,PSD=NULL,h=NULL,theta=NULL,steady=3,init=5,ths=NULL)
{
  if(is.null(h)){
    x=1:42
    h= 1.1163*exp(0.3699*x) 
  }
  options(warn=-1)
  #decalare the group data incase group variable is not null
  addoutput=NULL
  #set parameters of optimsation functions######################
  if(is.null(data)){
    data=data.frame(cbind(I,time)) 
    names(data)=c("I","time") 
  }
  if(!is.null(data)){
    if(!is.null(data$I)){
      I="I"
    }
    if(!is.null(data$time)){
      time="time"
    }
    
    if(!is.null(data$m)){
      m=data$m
    }
    
    if(!is.null(data$y)){
      y=mean(data$y)
    }
    
    if(!is.null(data$b)){
      b=mean(data$b)
    }
    if(!is.null(data$pb)){
      pb=mean(data$pb)
    }
    
    
    
    if(!is.null(data$r)){
      r=mean(data$r)
    }
    
    if(!is.null(data$thr)){
      thr=mean(data$thr)
    }
    
    if(!is.null(data$tho)){
      tho=mean(data$tho)
    }
    
    if(!is.null(data$init)){
      init=mean(data$init)
    }
    
    if(!is.null(data$steady)){
      steady=mean(data$steady)
    }
    
    if(!is.null(data$ths)){
      ths=mean(data$ths)
    }
    
  }
  
  if(!is.null(names(h))){
    h3=h
    h=h3$h
    theta=h3$theta
  }
  
  if(!is.null(names(theta))){
    theta3=theta
    theta=theta3$theta
    h2=theta3$h
    if(is.null(h2)){
      h=h
    }else{
      h=h2
    }
  }
  #n from PSD
  
  p <- 1-(pb/2.65)
  psd2=NULL
  texture <- NULL
  sand<-NULL
  silt<-NULL
  clay<-NULL
  
  
  if(!is.null(PSD)&&is.null(n)){
    psd2<-(lass3(PSD,p=p))
    n<-psd2$n
    texture <- texture(psd2)
    sand<-psd2$sand
    silt<-psd2$silt
    clay<-psd2$clay
  }
  
  
  l<-BESTlass (data=data,time=time,I=I,S=S,n=n,m=m,y=y,
               b=b,r=r,pb=pb,tho=tho,thr=thr,PSD=PSD,h=h,steady=steady,
               init=init,ths=ths)
  yb<-BESTyb (data=data,time=time,I=I,S=S,n=n,
              m=m,y=y,b=b,r=r,pb=pb,tho=tho,thr=thr,PSD=PSD,h=h,
              steady=steady,ths=ths,init=init)
  
  factor=list(BESTy=yb$BEST,BESTl=l$BEST,data=data,time=time,tmaxl=l$tmax,
              ktotl=l$ktot,I=I,Sl=l$S,A=l$A,B=l$B,C=yb$C,Ksl=l$Ks,kstepl=l$kstep,hgl=l$hg,
              cpl=l$cp,slopel=l$slope,interceptl=l$intercept,n2=l$n2,m=l$m,n=n,ths=l$ths,thr=l$thr,
              mod_thetal=l$theta_mod,Krl=l$Kr,Kl=l$K,formularl=l$BESTF,
              tmaxy=yb$tmax,ktoty=yb$ktot,Sy=yb$S,Ksy=yb$Ks,kstepy=yb$kstep,hgy=yb$hg,
              cpy=yb$cp,slopey=yb$slope,intercepty=yb$intercept,tho=tho,
              mod_thetay=yb$mod_theta,Kry=yb$Kr,Ky=yb$K,formulary=yb$BESTF,
              Ksb=yb$KsB,Sb=yb$SB,hgb=yb$hgB,Krb=yb$KrB,Kb=yb$KB,mod_thetab=yb$mod_thetaB,
              psd=psd2, I_modl=l$I_mod,I_mody=yb$I_mod,I_modb=yb$I_modb,texture=texture,sand=sand,silt=silt,
              clay=clay,h=h,theta=theta,q=l$q,init=init,pb=pb,PSDF=NULL)
  
  factor$call<-match.call()
  
  class(factor)<-"BEST"
  factor
}
#' @rdname BEST
#' @export
BESTlass <- function(data=NULL,time,I,S=0.1,n=NULL,m="b",y=0.75,b=0.6,r=75,pb=1.2,tho=0.169,
                     thr=0,PSD=NULL,h=NULL, theta = NULL,steady=3,init=5,ths=NULL)
{
  if(is.null(h)){
    x=1:42
    h= 1.1163*exp(0.3699*x) 
  }
  
  if(is.null(data)){
    data=data.frame(cbind(I,time)) 
    names(data)=c("I","time")
  }
  if(!is.null(data)){
    if(!is.null(data$I)){
      I="I"
    }
    if(!is.null(data$time)){
      time="time"
    }
    
    if(!is.null(data$m)){
      m=data$m
    }
    
    if(!is.null(data$y)){
      y=mean(data$y)
    }
    
    if(!is.null(data$b)){
      b=mean(data$b)
    }
    if(!is.null(data$pb)){
      pb=mean(data$pb)
    }
    
    
    
    if(!is.null(data$r)){
      r=mean(data$r)
    }
    
    if(!is.null(data$thr)){
      thr=mean(data$thr)
    }
    
    if(!is.null(data$tho)){
      tho=mean(data$tho)
    }
    if(!is.null(data$init)){
      init=mean(data$init)
    }
    
    if(!is.null(data$steady)){
      steady=mean(data$steady)
    }
    
    if(!is.null(data$ths)){
      ths=mean(data$ths)
    }
    
  }
  
  if(!is.null(names(h))){
    h3=h
    h=h3$h
    theta=h3$theta
  }
  
  if(!is.null(names(theta))){
    theta3=theta
    theta=theta3$theta
    h2=theta3$h
    if(is.null(h2)){
      h=h
    }else{
      h=h2
    }
  }
  
  #reference
  #stop warnings from displaying
  options(warn=-1)
  #decalare the group data incase group variable is not null
  addoutput=NULL
  p <- 1-(pb/2.65)
  psd2=NULL
  texture <- NULL
  sand<-NULL
  silt<-NULL
  clay<-NULL
  if(!is.null(PSD)&&is.null(n)){
    psd2<-(lass3(PSD,p=p))
    n<-psd2$n
    texture <- texture(psd2)
    sand<-psd2$sand
    silt<-psd2$silt
    clay<-psd2$clay
  }
  #mualem condition
  if(m=="m"||m=="Mualem"||m=="mualem"){
    m=1-(1/n)
  }
  #burdine condition
  if(m=="b"||m=="Burdine"||m=="burdine"){
    m=1-(2/n)
  }
  
  n2=(2/(n*m))+3
  if(is.null(ths)){
    ths=1-(pb/2.65)
  }
  
  A=y/(r*(ths-tho))
  B=((2-b)/3)*(1-(tho/ths)^n2)+(tho/ths)^n2
  #C=(1/(2*(1-((tho/ths)^n2))*(1-b)))*log(1/b)
  
  #set parameters of optimsation functions######################
  ones <- c(S=S) # all ones start
  #cumulative function ###############################
  i=init
  Kso=NULL
  while(i<=nrow(data)){
    data2=data[1:i,]
    steady1=(nrow(data)-steady)+1
    steadydata=data[steady1:nrow(data),]
    newdata=as.data.frame(cbind(steadydata[time],steadydata[I]))
    colnames(newdata)=c("time","I")
    mod=lm(newdata$I~newdata$time)
    intercept=coef(mod)[[1]]
    slope=coef(mod)[[2]]
    #BESTF<-paste(I,"~((S*sqrt(time))+((",A,"*(S^2))+(",B,"*(i-(",A,"*(S^2)))))*time)")
    BESTF=paste(I,"~(S*sqrt(time)+((",A,"*S^2)+(",B,"*(",slope,"-(",A,"*S^2))))*time)")
    #BESTF=paste(I,"~((S*sqrt(time))+(",A,"*(1-",B,")*(S^2)+",slope,")*time)")
    
    BEST<- nlxb(BESTF, start = ones, trace = FALSE, data = data2)
    S=coef(BEST)[["S"]]
    Ks=slope-(A*(S^2))
    #Ks=is-(A*(S^2))
    
    tmax=(1/(4*((1-B)^2)))*(S/Ks)^2
    tmax=(1/(4*(1-B)^2))*(S/Ks)^2
    if(data[time][i,]<tmax){
      So=S
      Kso=Ks
      tmaxo=tmax
      ktot=data[time][i,]
      #print(cbind(i,S,tmax,Ks,data[time][i,]))
      kstep=i
    }else{
      #print(cbind(i,S,tmax,Ks,data[time][i,]))
      
      break
    }
    i=i+1
  }
  #steady=data[16:18,]
  if(is.null(Kso)){
    So=S
    Kso=Ks
    tmaxo=tmax
    ktot=data[time][i,]
    kstep=i
  }
  cp=gamma(1+(1/n))*(gamma((m*n2)-(1/n))/(gamma(m*n2))+(gamma((m*n2)+m-(1/n))/gamma((m*n2)+m)))
  hg=S^2/(cp*(ths-tho)*((1-(tho/ths)^n2)*Ks))
  
  theta_est=thr+((ths-thr)*((1+(h/hg)^n)^(-m)))
  Kr=(((theta_est-thr)/(ths-thr))^n2)
  K=Kr*Kso
  
  time2=data[[time]]
  I2=data[[I]]
  q=data.frame(t=numeric(),q=numeric())
  
  t2=1
  while(t2<=length(time2)){
    if(t2==1){
      time3=((time2[t2]^0.5+0^0.5)/2)^2
      I3=(I2[t2]-0)/(time2[t2]-0)
    }else{
      
      time3=((time2[t2]^0.5+time2[t2-1]^0.5)/2)^2
      I3=(I2[t2]-I2[t2-1])/(time2[t2]-time2[t2-1])
    }
    
    q=rbind(q,data.frame(t=time3,q=I3))
    t2=t2+1
    
  }
  
  I_mod=(So*sqrt(data[[time]])+((A*So^2)+(B*(slope-(A*So^2))))*data[[time]])
  #return varibales ########################################
  list(BEST=BEST,data=data,time=time,tmax=tmaxo,ktot=ktot,I=I,
       S=So,A=A,B=B,C=C,Ks=Kso,kstep=kstep,hg=hg,cp=cp,slope=slope,
       intercept=intercept,n2=n2,m=m,n=n,ths=ths,thr=thr,theta_mod=theta_est,I_mod=I_mod,
       Kr=Kr,K=K,q=q,formular=BESTF)
}

#' @rdname BEST
#' @export
BESTyb <- function(data=NULL,time,I,S=0.1,n=NULL,m="b",y=0.75,b=0.6,r=75,pb=1.2,tho=0.169,
                   thr=0,PSD=NULL,h=NULL,theta = NULL,steady=3,init=5,ths=NULL)
{
  #reference
  #stop warnings from displaying
  options(warn=-1)
  if(is.null(h)){
    x=1:42
    h= 1.1163*exp(0.3699*x) 
  }
  #decalare the group data incase group variable is not null
  if(is.null(data)){
    data=data.frame(cbind(I,time)) 
    names(data)=c("I","time")
  }
  if(!is.null(data)){
    if(!is.null(data$I)){
      I="I"
    }
    if(!is.null(data$time)){
      time="time"
    }
    
    if(!is.null(data$m)){
      m=data$m
    }
    
    if(!is.null(data$y)){
      y=mean(data$y)
    }
    
    if(!is.null(data$b)){
      b=mean(data$b)
    }
    if(!is.null(data$pb)){
      pb=mean(data$pb)
    }
    
    
    
    if(!is.null(data$r)){
      r=mean(data$r)
    }
    
    if(!is.null(data$thr)){
      thr=mean(data$thr)
    }
    if(!is.null(data$tho)){
      tho=mean(data$tho)
    }
    if(!is.null(data$init)){
      init=mean(data$init)
    }
    
    if(!is.null(data$steady)){
      steady=mean(data$steady)
    }
    
    if(!is.null(data$ths)){
      ths=mean(data$ths)
    }
    
  }
  
  if(!is.null(names(h))){
    h3=h
    h=h3$h
    theta=h3$theta
  }
  
  if(!is.null(names(theta))){
    theta3=theta
    theta=theta3$theta
    h2=theta3$h
    if(is.null(h2)){
      h=h
    }else{
      h=h2
    }
  }
  addoutput=NULL
  #set parameters of optimsation functions######################
  p <- 1-(pb/2.65)
  psd2=NULL
  texture <- NULL
  sand<-NULL
  silt<-NULL
  clay<-NULL
  if(!is.null(PSD)&&is.null(n)){
    psd2<-(lass3(PSD,p=p))
    n<-psd2$n
    texture <- texture(psd2)
    sand<-psd2$sand
    silt<-psd2$silt
    clay<-psd2$clay
  }
  #mualem condition
  if(m=="m"||m=="Mualem"||m=="mualem"){
    m=1-(1/n)
  }
  #burdine condition
  if(m=="b"||m=="Burdine"||m=="burdine"){
    m=1-(2/n)
  }
  
  n2=(2/(n*m))+3
  if(is.null(ths)){
    ths=1-(pb/2.65)
  }
  S
  A=y/(r*(ths-tho))
  B=((2-b)/3)*(1-(tho/ths)^n2)+(tho/ths)^n2
  C=(1/(2*(1-((tho/ths)^n2))*(1-b)))*log(1/b)
  
  
  
  ones <- c(S=S) # all ones start
  
  #cumulative function ###############################
  Kso=NULL
  i=init
  while(i<=nrow(data)){
    data2=data[1:i,]
    steady1=(nrow(data)-steady)+1
    steadydata=data[steady1:nrow(data),]
    newdata=as.data.frame(cbind(steadydata[time],steadydata[I]))
    colnames(newdata)=c("time","I")
    mod=lm(newdata$I~newdata$time)
    intercept=coef(mod)[[1]]
    slope=coef(mod)[[2]]
    BESTF=paste(I,"~((",A,"*S^2)+(",slope,"-",A,"*S^2)*time+(",C,"*(S^2/(",slope,"-",A,"*S^2))))")
    BESTF=paste(I,"~((",A,"*S^2)+(",C,"*((S^2)/",intercept,"))*time+(",C,"*(S^2/(",C,"*((S^2)/",intercept,")))))")
    BESTF=paste(I,"~S*sqrt(time)+(A*(S^2)+(B*C*(S^2)/intercept)*time)")
    BESTF=paste(I,"~((S*sqrt(time))+((",A,"*(S^2)+(",B,"*",C,"*","(S^2))/",intercept,")*time))")
    BEST<- nlxb(BESTF, start = ones, trace = FALSE, data = data2)
    S=coef(BEST)[["S"]]
    Ks=(C*S^2)/intercept
    #Ks=is-(A*(S^2))
    
    tmax=(1/(4*((1-B)^2)))*(S/Ks)^2
    tmax=(1/(4*(1-B)^2))*(S/Ks)^2
    if( data[time][i,]<tmax){
      #print(cbind(i,S,tmax,Ks,data[time][i,]))
      So=S
      Kso=Ks
      tmaxo=tmax
      ktot=data[time][i,]
      kstep=i
      #slopeo=slope
      #intercepto=intercept
    }else{
      #print(cbind(i,tmax,Ks,data[time][i,]))
      break
    }
    i=i+1
  }
  if(is.null(Kso)){
    So=S
    Kso=Ks
    tmaxo=tmax
    ktot=data[time][i,]
    kstep=i
  }
  cp=gamma(1+(1/n))*(gamma((m*n2)-(1/n))/(gamma(m*n2))+(gamma((m*n2)+m-(1/n))/gamma((m*n2)+m)))
  hg=S^2/(cp*(ths-tho)*((1-(tho/ths)^n2)*Ks))
  #Bagarello et al. (2013)
  SB=sqrt(slope/(A+(C/intercept)))
  KsB=(C*SB^2)/intercept
  hgB=SB^2/(cp*(ths-tho)*((1-(tho/ths)^n2)*KsB))
  
  theta_est=thr+((ths-thr)*((1+(h/hg)^n)^(-m)))
  Kr=(((theta_est-thr)/(ths-thr))^n2)
  K=Kr*Kso
  theta_estB=thr+((ths-thr)*((1+(h/hgB)^n)^(-m)))
  KrB=(((theta_est-thr)/(ths-thr))^n2)
  KB=Kr*KsB
  I_mod=((So*sqrt(data[[time]]))+((A*(So^2)+(B*C*(So^2))/intercept)*data[[time]]))
  I_modb=((SB*sqrt(data[[time]]))+((A*(SB^2)+(B*C*(SB^2))/intercept)*data[[time]]))
  
  #return varibales ########################################
  list(BEST=BEST,data=data,time=time,tmax=tmaxo,ktot=ktot,kstep=kstep,KsB=KsB,
       hgB=hgB,I=I,SB=SB,S=So,A=A,B=B,C=C,Ks=Kso,kstep=kstep,hg=hg,cp=cp,slope=slope,
       intercept=intercept,n2=n2,m=m,n=n,ths=ths,formular=BESTF,thr=thr,
       ths=ths,phi=m*n,mod_theta=theta_est,I_mod=I_mod,I_modb=I_modb,h=h,
       Kr=Kr,K=K,mod_thetaB=theta_estB,KrB=KrB,KB=KB)
}
#' @export
#' @rdname BEST
group.BEST=function(data = NULL, time, I, S = 0.1, n = NULL, m = "b",
                    y = 0.75, b = 0.6, r = 75, pb = 1.2, tho = 0.169, thr = 0,
                    PSD = NULL, h = NULL, theta = NULL, steady = 3,
                    init = 5, ths = NULL,group,plot=TRUE,layout=c(2,2)
                    ,hlog=NULL,klog=NULL,opar=par(mar=c(2,2,1.8,2))){
  
  if(is.null(h)){
    x=1:42
    h= 1.1163*exp(0.3699*x) 
  }
  if(plot==TRUE){
   # dev.new()
    opar1 <- par()
    par(oma=c(3,1.5,1,1.5))
    par(mfrow=layout) 
    legend=FALSE
    }
  h2=NULL
  theta2=NULL
  parameters=data.frame(group=factor(),Sl=numeric(),Sy=numeric(),Sb=numeric(),
                        Ksl=numeric(),Ksy=numeric(),Ksb=numeric(),
                        hgl=numeric(),hgy=numeric(),hgb=numeric(),
                        n=numeric(),m=numeric(),ths=numeric(),
                        tho=numeric(),thr=numeric(),A=numeric(),B=numeric()
                        ,C=numeric())
  
  
  gof1=data.frame()
  gof1box=data.frame()
  
  predict=data.frame(group=factor(),mod_thetal=numeric(),mod_thetay=numeric(),mod_thetab=numeric(),
                     Kl=numeric(),Ky=numeric(),Kb=numeric())
  aggdata =row.names(table(data[group]))
  #if(!is.null(data$n)){
  #  n=aggregate(data$n,data[group],FUN=mean)
  #}
  for(i in 1:length(aggdata)) {
    
    single=data[data[group]==aggdata[i],]
    if(!is.null(PSD)){
      PSD2=PSD[PSD[group]==aggdata[i],]
    }
    if(is.null(n)){
      if(!is.null(single$n)){
        n2=single$n[1]
        # print(n)
      }
    }else{
      if(length(n)==length(aggdata)){
        n2=n[i]
      }else{
        n2=n
      }
    }
    
    if(is.null(m)){
      if(!is.null(single$m)){
        m2=single$m[1]
      }
    }else{
      if(length(m)==length(aggdata)){
        m2=m[i]
      }else{
        m2=m
      }
    }  
    
    
      if(!is.null(single$y)){
        y2=single$y[1]
      
    }else{
      if(length(y)==length(aggdata)){
        y2=y[i]
      }else{
        y2=y
      }
    }
    
    
      if(!is.null(single$b)){
        b2=single$b[1]
      
    }else{
      if(length(b)==length(aggdata)){
        b2=b[i]
      }else{
        b2=b
      }
    }
    
    
      if(!is.null(single$r)){
        r2=single$r[1]
      
    }else{
      if(length(r)==length(aggdata)){
        r2=r[i]
      }else{
        r2=r
      }
    }
    
    
      if(!is.null(single$pb)){
        pb2=single$pb[1]
      
    }else{
      if(length(pb)==length(aggdata)){
        pb2=pb[i]
      }else{
        pb2=pb
      }
    }
    
    
      if(!is.null(single$tho)){
        tho2=single$tho[1]
      
    }else{
      if(length(tho)==length(aggdata)){
        tho2=tho[i]
      }else{
        tho2=tho
      }
    }
    
   
      if(!is.null(single$thr)){
        thr2=single$thr[1]
      
    }else{
      if(length(thr)==length(aggdata)){
        thr2=thr[i]
      }else{
        thr2=thr
      }
    }
    
    if(is.null(ths)){
      ths2=ths
      if(!is.null(single$ths)){
        ths2=single$ths[1]
      }
    }else{
      if(length(ths)==length(aggdata)){
        ths2=ths[i]
      }else{
        ths2=ths
      }
    }
    
    
      if(!is.null(single$steady)){
        steady2=single$steady[1]
      
    }else{
      if(length(steady)==length(aggdata)){
        steady2=steady[i]
      }else{
        steady2=steady
      }
    }
    
    
      if(!is.null(single$S)){
        S2=single$S[1]
      
    }else{
      if(length(S)==length(aggdata)){
        S2=S[i]
      }else{
        S2=S
      }
    }
    
      if(!is.null(single$init)){
        init2=single$init[1]
      
    }else{
      if(length(init)==length(aggdata)){
        init2=init[i]
      }else{
        init2=init
      }
    }
    if(!is.null(names(h))){
      h4=h[h[group]==aggdata[i],]
      h2=h4$h
      theta2=h4$theta
    }
    
    if(!is.null(names(theta))&&!is.null(theta)){
      theta4=theta[theta[group]==aggdata[i],]
      
      if(is.null(h2)){
        h2=theta4$h
      }
      if(is.null(theta2)){
        theta2=theta4$theta
      }
    }
    
    if(is.null(h2)){
      h2=h
    }
    
    if(is.null(theta2)){
      theta2=theta
    }
    
    if(is.null(theta2)&&is.null(h2)){
      h2=seq(0, 15000, by = 10)
      theta2=NULL
    }
    
      mod=BEST(data = single,  S = S2, n = n2, m = m2, y = y2,
             b = b2, r = r2, pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
             h = h2, theta = theta2, steady = steady2, init = init2,
             ths = ths2)
    
    gof2=vadose.tryCatch(gof.BEST(mod))$value
    #print(row.names(gof2))
    
    #print(gof2)
    if(class(gof2)[1]=="simpleError"){
      gof3=cbind(aggdata[i],NA,NA,NA)
      gof4=rbind(NA,NA,NA)
    }else{
      gof3=cbind(aggdata[i],gof2[1,],gof2[2,],gof2[3,])
      gof4=rbind(gof2[1,],gof2[2,],gof2[3,])
      
    }
    if(!is.null(theta2)){
      if(class(gof2)[1]=="simpleError"){
        gof3=cbind(aggdata[i],NA,NA,NA)
        gof4=rbind(NA,NA,NA)
        
      }else{
      gof3=cbind(aggdata[i],gof2[1,],gof2[2,],gof2[3,])
      gof4=rbind(gof2[1,],gof2[2,],gof2[3,])
      }
    }
    if(class(gof2)[1]!="simpleError"){
      
    gof1=rbind(gof1,data.frame(gof3))
    gof4$ID=aggdata[i]
    gof4$model=row.names(gof4)
    gof1box=rbind(gof1box,data.frame(gof4))
    #gof1box=cbind(gof1box,aggdata[i])
    }
    parameters=rbind(parameters,data.frame(group=aggdata[i],Sl=mod$Sl,Sy=mod$Sy,
                                           Sb=mod$Sb,Ksl=mod$Ksl,Ksy=mod$Ksy,Ksb=mod$Ksb,
                                           hgl=mod$hgl,hgy=mod$hgy,hgb=mod$hgb,
                                           n=mod$n,m=mod$m,ths=mod$ths,
                                           tho=mod$tho,thr=mod$thr,A=mod$A,B=mod$B
                                           ,C=mod$C))
    
    predict=rbind(predict,data.frame(group=aggdata[i],mod_thetal=mod$mod_thetal,
                                     mod_thetay=mod$mod_thetay,mod_thetab=mod$mod_thetab,
                                     Kl=mod$Kl,Ky=mod$Ky,Kb=mod$Kb,q=mean(mod$q$q)))
    
    if(plot==TRUE){
      plot(mod,type="hk",main=aggdata[i],legend=legend,hlog=hlog,klog=klog)
      }
    
    
    #coef3=coef(mod)
    #print(aggdata[i])
    #print(single$ID[[i]])
    #print(data.frame(coef3$Ks))
    
  }
  if(plot==TRUE){
    
    col=c("blue", "red", "darkgreen","white")
    plot_colors=col
    par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    #legend("bottom", c("IM", "IBD", "1R", "2R"), xpd = TRUE, horiz = TRUE, 
    #     inset = c(0,0), bty = "n", pch = c(4, 2, 15, 19), col = 1:4, cex = 2)
    if(is.null(theta2)){
      legend("bottom", (c("Bagarello","Lassabatere","Yilmaz","Water Content")), 
             cex=0.9,col=plot_colors,lty=1:3,lwd=2, bty="n",xpd = TRUE, horiz = TRUE, 
             inset = c(0,0))
    }else{
      legend("bottom", (c("Measured","Bagarello","Lassabatere","Yilmaz","Water Content")), 
             cex=0.9, pch=c("o","","",""),
             col=c("black",plot_colors),lty=0:3,lwd=2, bty="n",xpd = TRUE, horiz = TRUE, 
             inset = c(0,0))
    }
    par(opar1)
  }
  row.names(gof1)=NULL
  if(!is.null(theta2)){
    gof1name=paste(names(gof2[1,]),".Lassabatere_I",sep="")
    gof2name=paste(names(gof2[2,]),".Yilmaz_I",sep="")
    gof3name=paste(names(gof2[3,]),".Bagarello_I",sep="")
    gof4name=paste(names(gof2[4,]),".Lassabatere_THETA",sep="")
    gof5name=paste(names(gof2[5,]),".Bagarello_THETA",sep="")
    gof6name=paste(names(gof2[6,]),".Yilmaz_THETA",sep="")
    gofnames=cbind(gof4name,gof5name,gof6name)
  }else{
  gof1name=paste(names(gof2[1,]),"Lassabatere",sep="")
  gof2name=paste(names(gof2[2,]),"Yilmaz",sep="")
  gof3name=paste(names(gof2[3,]),"Bagarello",sep="")
  gofnames=cbind(gof1name,gof2name,gof3name)
  }
  #rownames(gof1)=vadose.tryCatch(gof1[[1]])
  #gof1[1]=NULL
  names(gof1)=vadose.tryCatch(gofnames)
  tabname=names(gof1[1])
  TEXT2=""
  tabname2=sub("c", TEXT2, tabname, ignore.case =TRUE, fixed=FALSE) 
 
  
  tabnamesplit=strsplit(tabname2, ",")
  group="ID"
  tabnamesplit2=c(group,tabnamesplit[[1]])
  tabname2=gsub("[[:punct:]]", TEXT2, tabnamesplit2, ignore.case =TRUE, fixed=FALSE)
  tabname2=gsub("[[:blank:]]", TEXT2, tabname2, ignore.case =TRUE, fixed=FALSE)
  
  tabname2=gsub("THETA", ".THETA", tabname2, ignore.case =TRUE, fixed=FALSE)
  
  names(gof1)=tabname2
  
  if(plot==TRUE){
    par(opar1)
  }
  
  #row.names(predict)=predict$group
  #predict$group=NULL
  #boxplot
  #boxgof=gof1
  #boxgof[1]=NULL
  #if(plot==TRUE){
    
  #}
  #print(gof1box)
  #print(names(boxgof))
  list(parameters=parameters,predict=predict,statistics=gof1,statistics2=gof1box)
}


#' @export
#' @rdname BEST
cal.BEST=function(data = NULL, time, I, S = 0.1, n = NULL, m = "b", 
                  y = c(0.6,0.8,0.05),b = c(0.1,0.9,0.1), r = 75, pb = c(1.0, 1.8,0.1),
                  tho = c(0.01,0.8,0.05), thr = c(0,0.2,0.05), PSD = NULL,
                  h = NULL, theta = NULL, steady = c(2,10,1), init = c(5,7,1),
                  ths = NULL,model="all",group=NULL){
  if(is.null(h)){
    x=1:42
    h= 1.1163*exp(0.3699*x) 
  }
  if(!is.null(group)){
    parameters=data.frame()
    parameters2=data.frame()
    
    fit=data.frame()
    #group="TownName"
    aggdata =row.names(table(data[group]))
    h2=NULL
    theta2=NULL
    for(i in 1:length(aggdata)) {
      
      single=data[data[group]==aggdata[i],]
      if(!is.null(PSD)){
        PSD2=PSD[PSD[group]==aggdata[i],]
      }
      if(is.null(n)){
        if(!is.null(single$n)){
          n2=single$n[1]
          # print(n)
        }
      }else{
        if(length(n)==length(aggdata)){
          n2=n[i]
        }else{
          n2=n
        }
      }
      
      if(is.null(m)){
        if(!is.null(single$m)){
          m2=single$m[1]
        }
      }else{
        if(length(m)==length(aggdata)){
          m2=m[i]
        }else{
          m2=m
        }
      }  
      
      
      if(!is.null(single$y)){
        y2=single$y[1]
        
      }else{
        if(length(y)==length(aggdata)){
          y2=y[i]
        }else{
          y2=y
        }
      }
      
      
      if(!is.null(single$b)){
        b2=single$b[1]
        
      }else{
        if(length(b)==length(aggdata)){
          b2=b[i]
        }else{
          b2=b
        }
      }
      
      
      if(!is.null(single$r)){
        r2=single$r[1]
        
      }else{
        if(length(r)==length(aggdata)){
          r2=r[i]
        }else{
          r2=r
        }
      }
      
      
      if(!is.null(single$pb)){
        pb2=single$pb[1]
        
      }else{
        if(length(pb)==length(aggdata)){
          pb2=pb[i]
        }else{
          pb2=pb
        }
      }
      
      
      if(!is.null(single$tho)){
        tho2=single$tho[1]
        
      }else{
        if(length(tho)==length(aggdata)){
          tho2=tho[i]
        }else{
          tho2=tho
        }
      }
      
      
      if(!is.null(single$thr)){
        thr2=single$thr[1]
        
      }else{
        if(length(thr)==length(aggdata)){
          thr2=thr[i]
        }else{
          thr2=thr
        }
      }
      
      if(is.null(ths)){
        ths2=ths
        if(!is.null(single$ths)){
          ths2=single$ths[1]
        }
      }else{
        if(length(ths)==length(aggdata)){
          ths2=ths[i]
        }else{
          ths2=ths
        }
      }
      
      
      if(!is.null(single$steady)){
        steady2=single$steady[1]
        
      }else{
        if(length(steady)==length(aggdata)){
          steady2=steady[i]
        }else{
          steady2=steady
        }
      }
      
      
      if(!is.null(single$S)){
        S2=single$S[1]
        
      }else{
        if(length(S)==length(aggdata)){
          S2=S[i]
        }else{
          S2=S
        }
      }
      
      if(!is.null(single$init)){
        init2=single$init[1]
        
      }else{
        if(length(init)==length(aggdata)){
          init2=init[i]
        }else{
          init2=init
        }
      }
      if(!is.null(names(h))){
        h4=h[h[group]==aggdata[i],]
        h2=h4$h
        theta2=h4$theta
      }
      
      if(!is.null(names(theta))&&!is.null(theta)){
        theta4=theta[theta[group]==aggdata[i],]
        
        if(is.null(h2)){
          h2=theta4$h
        }
        if(is.null(theta2)){
          theta2=theta4$theta
        }
      }
      
      if(is.null(h2)){
        h2=h
      }
      
      if(is.null(theta2)){
        theta2=theta
      }
      mod=cal.BEST(data = single,  S = S2, n = n2, m = m2, y = y2,
               b = b2, r = r2, pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
               h = h2, theta = theta2, steady = steady2, init = init2,
               ths = ths2)
       # print(mod$parameters)
        
        
        
      if(!is.null(mod$parameters)){
        thisparameters=mod$parameters
        thisparameters$group=aggdata[i]
        #print(mod$parameters)
      parameters= rbind(parameters,data.frame(thisparameters))
      print(aggdata[i])
      }
      
      if(!is.null(mod$parameters2)){
        
        thisparameters2=mod$parameters2
        thisparameters2$group=aggdata[i]
        parameters2=rbind(parameters2,data.frame(thisparameters2))
        # print(aggdata[i])
        
      }
      
      if(!is.null(mod$fit)){
        thisfit=mod$fit
        thisfit$model=row.names(thisfit)
        
        fit=rbind(fit,data.frame(cbind(aggdata[i],thisfit)))
              }
      
    }
    #row.names(parameters)=parameters$group
    #parameters$group=NULL
    #row.names(fit)=fit[[1]]
    #fit[1]=NULL
   
    fit$group=fit$aggdata.i
    fit[1]=NULL
    list(parameters=parameters,parameters2=parameters2,statistics=fit)
    
  }else{
    
    if(is.null(theta)){
    return  (print("Observation or measured data 'theta' is needed"))
  }
  
  
    para=NULL
    count=0 
    
    if(length(y)==1){
      y=c(y,y,y)
    }else{
      para=cbind(para,y=y)
      count=count+1
    }
    
    if(length(b)==1){
      b=c(b,b,b)
    }else{
      para=cbind(para,b=b)
      count=count+1
    }
    
    if(length(tho)==1){
      tho=c(tho,tho,tho)
    }else{
      para=cbind(para,tho=tho)
      count=count+1
    }
    
    n4=TRUE
    if(is.null(n)){
      n4=NULL
      n=c(1,1,1)
    }else{
      if(length(n)==1){
        n=c(n,n,n)
      }else{
        para=cbind(para,n=n)
        count=count+1
      }
    }
    if(length(thr)==1){
      thr=c(thr,thr,thr)
    }else{
      para=cbind(para,thr=thr)
      count=count+1
    }
    
    if(length(init)==1){
      init=c(init,init,init)
    }else{
      para=cbind(para,init=init)
      count=count+1
    }
    
    
    if(length(pb)==1){
      pb=c(pb,pb,pb)
      
    }else{
      para=cbind(para,pb=pb)
      count=count+1
    }
    
    
    
    if(length(steady)==1){
      steady=c(steady,steady,steady)
    }else{
      para=cbind(para,steady=steady)
      count=count+1
    }
    print("Iteration Started. Kindly wait.................")
    pnames=names(data.frame(para))
    count=length(pnames)
    #return(list(para=para))
    #print(pnames)
    if(count==0){
      count=1
    }
    rmse=100
    r2=0
    NRMSE=100
    d=-1
    init2=seq(init[1],init[2],init[3])
    steady2=seq(steady[1],steady[2],steady[3])
    b2=seq(b[1],b[2],b[3])
    y2=seq(y[1],y[2],y[3])
    pb2=seq(pb[1],pb[2],pb[3])
    tho2=seq(tho[1],tho[2],tho[3])
    thr2=seq(thr[1],thr[2],thr[3])
    n2=seq(n[1],n[2],n[3])
  fit2=NULL
  para2=NULL
  mod=NULL
  mod2=NULL
  parameters=NULL
  for(init1 in init2){
    for(steady1 in steady2) {
      for(tho1 in tho2) {
        for(thr1 in thr2){
          for(pb1 in pb2){
            for(y1 in y2){
              for(b1 in b2){
                for(n1 in n2){
                  if(is.null(n4)){
                    #print("YES")
                    n3=NULL
                  }else{
                    n3=n1
                  }
                  mod=BEST(data = data, time=time, I=I, S =S, n = n3, m = m, y = y1,b = b1, 
                           r = r, pb = pb1, tho = tho1, thr = thr1, PSD = PSD,
                           h = h, theta = theta, steady = steady1, init = init1,ths = ths)
                  #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
                  
                  fit=vadose.tryCatch(gof.BEST(mod,P=count))$value
                  if(class(fit)[1]!="simpleError"){
                    
                  
                  if(model=="all"||model=="All"){
                    fitrmse1=fit
                    fitrmse=mean(fitrmse1$RMSE[1:3]) 
                    fitr2=mean(fitrmse1$r2[1:3]) 
                    fitd=mean(fitrmse1$d[1:3]) 
                    fitNRMSE=mean(fitrmse1$NRMSE[1:3]) 
                  }else{
                    fitrmse1=fit[which(row.names(fit)==model),]
                    fitrmse=fitrmse1$RMSE 
                    fitr2=fitrmse1$r2 
                    fitd=fitrmse1$d 
                    fitNRMSE=fitrmse1$NRMSE 
                  }
                  #fitr2=fit$r2
                   # print(fitrmse)
                    #print(rmse)
                  if(fitrmse<rmse){
                    #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
                    parameters=(list(y=y1,b=b1,n=n3,tho=tho1,thr=thr1,pb=pb1,
                                     steady=steady1,init=init1)) 
                    #print(pnames)
                    para2=parameters[pnames]
                    print(data.frame(para2))
                    print(fitrmse1)
                    fit2=fitrmse1
                    r2=fit$r2
                    rmse=fitrmse
                    r2=fitr2
                    NRMSE=fitNRMSE
                    d=fitd
                    mod2=mod
                  }
                  
                  if(fitrmse==0){
                    break
                  }
                  }
                  #parameters=(list(y=y1,b=b1,tho=tho1,thr=thr1,pb=pb1,
                  #                steady=steady1,init=init1)) 
                  #para2=parameters[pnames]
                  #print(data.frame(para2))
                  #print(cbind(init1,steady1,b1))
                  #print(fitrmse1)
                  
                  
                  
                }
              }
              
            }
            
            
          }
          
          
        }
        
        
      }                            
    }
    
    
  }
  
  
  list(fit=fit2,statistics=fit2,parameters=parameters,parameters2=para2,mod=mod2)
}
}

#' @export
#' @rdname BEST
gof.BEST<-function(x=NULL,obs=NULL,est=NULL,P=NULL)
{
  object=x
  if(is.null(obs)){
    h=object$h
    mod_thetab=object$mod_thetab
    mod_thetal=object$mod_thetal
    mod_thetay=object$mod_thetay
    theta=object$theta
    obs=object$data[object$I]
    estb=object$I_modb
    estl=object$I_modl
    esty=object$I_mody
    time=object$data[object$time]
    gofb=(gof(object,c=obs,r=estb,P=P))
    gofl=(gof(object,c=obs,r=estl,P=P))
    gofy=(gof(object,c=obs,r=esty,P=P))
  }
  if(is.null(obs)){
    obs=object$fr
    est=object$predict
  }
  if(is.null(obs)){
    obs=object$theta
    est=object$predict2
  }
  if(is.null(obs)){
    obs=object$I
    est=object$predict
  }
  
  if(!is.null(theta)){  
    theta=data.frame(object$theta)
    #print(data.frame(theta))
    gofhb=(gof(object,c=theta,r=mod_thetab,P=P))
    gofhl=(gof(object,c=theta,r=mod_thetal,P=P))
    gofhy=(gof(object,c=theta,r=mod_thetay,P=P))
    data=rbind(gofhl,gofhy,gofhb)
    #row.names(data)=c("Lassabatere_Infiltration","Yilmaz_Infiltration","Bagarello_Infiltration","Lassabatere_Theta","Yilmaz_Theta","Bagarello_Theta")
    row.names(data)=c("Lassabatere","Yilmaz","Bagarello")
    
    }else{
    data=rbind(gofl,gofy,gofb)
    row.names(data)=c("Lassabatere","Yilmaz","Bagarello")
  }
  
  data
  #list(gof=data)
}

#' @export
#' @rdname BEST
coef.BEST<-function(object,...){
  
  if(is.null(object$A)){
    object=object$mod
  }
A=object$A
B=object$B
C=object$C
n=object$n
m=object$m
pb=object$pb
thr=object$thr
ths=object$ths
tho=object$tho
n2=object$n2

S=list(Sl=object$Sl[1],Sy=object$Sy[1],Sb=object$Sb[1])
Ks=list(Ksl=object$Ksl[1],Ksy=object$Ksy[1],Ksb=object$Ksb[1])
slope=list(slopel=object$slopel[1],slopey=object$slopey[1])
intercept=list(interceptl=object$interceptl[1],intercepty=object$intercepty[1])
hg=list(hgl=object$hgl[1],hgy=object$hgy[1],hgb=object$hgb[1])

list(S=S,Ks=Ks,hg=hg,slope=slope,intercept=intercept,A=A,B=B,C=C,n2=n2,n=n,m=m,pb=pb,thr=thr,ths=ths,tho=tho)

  
}

#' @export
#' @rdname BEST
predict.BEST<-function(object,h=NULL,...)
{
  if(is.null(object$A)){
    object=object$mod
  }
  if(is.null(h)){
    theta=data.frame(cbind(object$mod_thetal,object$mod_thetay,object$mod_thetab))
    names(theta)=c("theta.Lassabatere","theta.Yilmaz","theta.Bagarello")
    K=data.frame(cbind(object$Kl,object$Ky,object$Kb))
    names(K)=c("K.Lassabatere","K.Yilmaz","K.Bagarello")
    #return(predict)
    #predict=predict[[1]]
  }else{
    if(!is.null(object$Ks)){
      theta=object$thr+((object$ths-object$thr)*((1+(h/object$hg)^object$n)^(-object$m)))
      Kr=(((theta-object$thr)/(object$ths-object$thr))^object$n2)
      K=Kr*object$Ks
    }else{
      thetal=object$thr+((object$ths-object$thr)*((1+(h/object$hgl)^object$n)^(-object$m)))
      Krl=(((thetal-object$thr)/(object$ths-object$thr))^object$n2)
      Kl=Krl*object$Ksl
      
      thetay=object$thr+((object$ths-object$thr)*((1+(h/object$hgy)^object$n)^(-object$m)))
      Kry=(((thetay-object$thr)/(object$ths-object$thr))^object$n2)
      Ky=Kry*object$Ksy
      
      thetab=object$thr+((object$ths-object$thr)*((1+(h/object$hgb)^object$n)^(-object$m)))
      Krb=(((thetab-object$thr)/(object$ths-object$thr))^object$n2)
      Kb=Krb*object$Ksb
      theta=data.frame(cbind(thetal=thetal,thetay=thetay,thetab=thetab))
      names(theta)=c("theta.Lassabatere","theta.Yilmaz","theta.Bagarello")
      K=data.frame(cbind(Kl=Kl,Ky=Ky,Kb=Kb))
      names(K)=c("K.Lassabatere","K.Yilmaz","K.Bagarello")
    }
   
  }
  
  list(theta=theta,K=K)
  

}



#' @export
#' @rdname BEST
#plot function
plot.BEST<-function(x,main=NULL,xlab="Water Content",
                    ylab="",ylab2="",
                    hlog=NULL,klog=NULL,col=c("blue", "red", "darkgreen"),
                    units=c("s","mm"),mfrow=c(2,2),type="all",legend=TRUE,
                    opar=par(mar=c(2,2,1.5,2)),...)
{
  object=x
  #par(mfrow=c(2,2),mar=c(4, 4, 4, 4))de
  Kb=object$Kb
  Kl=object$Kl
  Ky=object$Ky
  h=object$h
  mod_thetab=object$mod_thetab
  mod_thetal=object$mod_thetal
  mod_thetay=object$mod_thetay
  theta=object$theta
  obs=object$data[[object$I]]
  estb=object$I_modb
  estl=object$I_modl
  esty=object$I_mody
  time=object$data[[object$time]]
  
  ###########################
  plot_colors=col
  #par(mar = c(5, 4, 4, 4) + 0.3)
  opar
  #par(mar=c(4.2, 3.8, 0.2, 0.2))
  #plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen")
  
  #ylab=paste(ylab,"(",units[2],")",sep="")
  #ylab2=paste(ylab2,"(",units[2],"/",units[1],")",sep="")
  
  ymax <- max(Kb, na.rm=TRUE)
  #dev.new(width=15, height=10)
  if(type=="all"){
    par(mfrow=mfrow)
  }
  
  if(type=="all"||type=="h"||type=="hk"||type=="k"){
    if(!is.null(theta)){
      if(!is.null(hlog)){
        plot(theta, h,xlab=xlab,ylab=ylab,log="y",main=main,...) # first plot
      }else{
        plot(theta, h,xlab=xlab,ylab=ylab,main=main,...) # first plot
        
      }
    }
    if(!is.null(hlog)){
      #ylab=paste(ylab,"[log]")
      if(!is.null(theta)){
        lines(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],log="y",...)
      }else{
        plot(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],main=main,log="y",...)
      }
      lines(mod_thetal,h,xlab=xlab,ylab=ylab,type="l", lty=2, lwd=2,col=plot_colors[2],log="y",...)
      lines(mod_thetay,h,xlab=xlab,ylab=ylab, type="l", lty=3, lwd=2,col=plot_colors[3],log="y",...)
      
      if(is.null(theta)){
        if(legend!=FALSE){
        legend("bottom", (c("Bagarello","Lassabatere","Yilmaz")), cex=0.8, col=plot_colors,lty=1:3, lwd=2, bty="n")
        }
          }else{
            if(legend!=FALSE){
              legend("bottom", (c("Measured","Bagarello","Lassabatere","Yilmaz")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
            }
      }
    }else{
      if(!is.null(theta)){
        lines(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],...)
        
      }else{
        plot(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],main=main,...)
      }
      
      lines(mod_thetal,h,type="l",lty=2, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[2],...)
      lines(mod_thetay,h,type="l",lty=3, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[3],...)
      if(is.null(theta)){
        if(legend!=FALSE){
                  legend("top", (c("Bagarello","Lassabatere","Yilmaz")), cex=0.8, col=plot_colors,lty=1:3, lwd=2, bty="n")
        }
          }else{
            if(legend!=FALSE){
              
        legend("top", (c("Measured","Bagarello","Lassabatere","Yilmaz")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
            }
      }
    }
    par(new = TRUE)
    if(!is.null(klog)){
      plot(mod_thetab, Kb, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[1],log="y",...)
      lines(mod_thetal, Kl, type="l",lty=2, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[2],log="y",...)
      lines(mod_thetay, Ky, type="l",lty=3, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[3],log="y",...)
      #ylab2=paste(ylab2,"[log]")
    }else{
      plot(mod_thetab, Kb, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[1],...)
      lines(mod_thetal, Kl,  type="l",lty=2, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[2],...)
      lines(mod_thetay, Ky, type="l",lty=3, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[3],...)
    }
    
    axis(side=4, ylim=c(0, ymax))
    mtext(ylab2, side=4, line=2,cex=0.8)
    index1=length(mod_thetab)
    index2=length(mod_thetab)-1
    leb=(mod_thetab[index2]-mod_thetab[index1])
    text(max(mod_thetab)-(0*6),max(Kb),"K")
    text(min(mod_thetab)+(min(mod_thetab)*0.30),max(Kb),"h")
    
  }
  
  if(type=="all"||type=="cum"||type=="I"||type=="I"){
    plot(time,obs,ylab="",main=paste("Cumulative infiltration(",units[2],")",sep=""),xlab=paste("Time(",units[1],")",sep=""))
    par(new = TRUE)
    #lines(time,estb,type="l",lty=1, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[1])
    lines(time,estl,type="l",lty=2, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[2])
    lines(time,esty,type="l",lty=3, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[3])
    #legend("bottomright", (c("Measured","Bagarello","Lassabatere","Yilmaz")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
    if(legend!=FALSE){
    legend("bottomright", (c("Measured","Lassabatere","Yilmaz","Bagarello")), cex=0.8, pch=c("o","",""),col=c("black",plot_colors[2],plot_colors[3]),lty=c(0,1,2,3),lwd=2, bty="n")
    }
  }
  if(type=="all"||type=="i"||type=="f"){
    plot(object$q,main="Infiltration rate",xlab=paste(names(object$q)[1],"(",units[1],")",sep=""),ylab=paste(names(object$q)[2],"(",units[2],"/",units[1],")",sep=""))
  }
  
  if(type=="all"||type=="psd"||type=="PSD"){
    plot(object$psd,xlab=2)
  }
  if(type=="all"){
    par(mfrow=c(1,1))
  }
}
