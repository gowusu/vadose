#input function ############################################
#default function
#' @title R Estimation of the Water Retention Shape Parameters from particle size distribution
#' @description The function implements Lassabatere et al (2006) particle size distribution function
#' for estimation of Water Retention curve. It also provides the parameters based on 
#' Minasny & McBratney (2007) studies.
#' @param data a dataframe of Particle Size Distribution. 
#' It can have column names of diameter of the particles "D" [mm], 
#' and fraction of the particles "fr" that ranges from 0-1.
#' @param D character. The corresponding diameter in the data. 
#' @param fr character. The corresponding fraction of the particles. 
#' @param N numeric. The initial value of the shape parameter.
#' @param M numeric. The initial value of the shape parameter.
#' @param Dg numeric. A scale parameter
#' @param p porosity[g/m^3]
#' @param S The percentage of sand in the soil. This is needed if there is no data on PSD.
#' @param C The percentage of clay in the soil. This is needed if there is no data on PSD
#' @param ksat Saturated hydraulic conductivity. If NULL, this will be internally  estimated with
#' \code{\link{ksat}}
#' @param group character. The name of the group variables if the data is from different areas.
#' @param object Model output object
#' @author George Owusu
#' @inheritParams BEST
#' @return 
#' \itemize{
#' \item{PD:} {  predicted fraction by mass of particles passing a particular diameter}
#' \item{poresizeindex:} { pore size distribution index[-]}
#' \item{D10: } { The diameter at which 10 per cent of the sample mass comprised of smaller particles }
#' \item{D60: } { The diameter at which 60 per cent of the sample mass comprised of smaller particles }
#' \item{sand:} { percentage of sand}
#' \item{clay:} { percentage of clay}
#' \item{silt:} { percentage of silt }
#' \item{n:} { shape parameter of van Genuchten water retention alpha parameter }
#' } 
#' @seealso \code{\link{zhuang3}}, \code{\link{zhuang4}},\code{\link{logarithmic}},
#' \code{\link{logistic2}}, \code{\link{logistic3}},\code{\link{fredlund3}},\code{\link{fredlund4}},
#'  \code{\link{jaky}},\code{\link{andersson}},\code{\link{gompertz}},
#' \code{\link{haverkamp}}
#' @export
#' @references 
#' \itemize{
#' \item{}{Lassabatere, L., Angulo-Jaramillo, R., Soria Ugalde, J. M., Cuenca, R., Braud, I., & Haverkamp, R. (2006). 
#' Beerkan Estimation of Soil Transfer Parameters through Infiltration 
#' Experiments-BEST. Soil Sci. Soc. Am. J., 70(2), 521-532. 
#' doi: 10.2136/sssaj2005.0026}
#' \item{}{Minasny, B., & McBratney, A. B. (2007). Estimating the Water Retention Shape Parameter from 
#' Sand and Clay Content Soil Sci. Soc. Am. J., 71(4), 1105-1110. doi: 10.2136/sssaj2006.0298N}
#' }
#' 
#' @examples
#' data=read.csv(system.file("ext","sys","soil2.csv",package="vadose"))
#' single<- subset(data, ID=="30B20_1")
#' mod=lass3 (data=single,p="sand",D="D",fr="Sand.")
#' plot(mod)
#' 
#' #group simulation
#' mod2=lass3 (data=data,D="D",fr="FractionSand",group="ID",ksat=TRUE)
#' mod2$coef
#' #generic function
lass3<-function(data=NULL,D=NULL,fr=NULL,N=0.1,M=0.1,Dg=0.1,p=NULL,S=NULL,C=NULL,ksat=NULL,group=NULL) UseMethod ("lass3")
#' @export
#' @rdname lass3
lass3.default<-function(data=NULL,D=NULL,fr=NULL,N=0.1,M=0.1,Dg=0.1,p=NULL,S=NULL,C=NULL,ksat=NULL,group=NULL)
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
      fr=data$I
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
      return(print("please provide the diameter,D"))
    }
    fr=(1+(Dg/D)^N)^-M
    return(print(list(fraction=fr)))
  }
  
  #convert from micrometers to mm
  if(max(data[[D]])>10){
    data[[D]]=data[[D]]/1000
  }
  data[data==0]<-0.0001
  #group function
  if(!is.null(group)){
    aggdata =row.names(table(data[group]))
    #create group data frame###################################
    addoutput=data.frame(groupid=factor(),N=numeric(),M=numeric(),Dg=numeric(),n=numeric(),m=numeric(),n2=numeric(),psi=numeric(),b=numeric(),porosity=numeric(),D10=numeric(),D60=numeric(),CU=numeric(),K=numeric(),KB=numeric(),r2=numeric(),
                         RMSE=numeric(),NRMSE=numeric(),F_value=numeric(),df=numeric(),p_value=numeric(),AIC=numeric(),BIC=numeric(),r2_adj=numeric(),RMSE_adj=numeric(),NRMSE_adj=numeric(),nashSutcliffe=numeric(),texture=numeric())
    i=1
    while(i<=length(aggdata)){
      #print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      p2=p
      if(length(p)>1){
        p2=p[i]
      }
      mod1=lass3(single,p=p2,ksat=ksat)
      gof1=gof(mod1)
      N=coef(mod1$lass3)[[1]]
      M=coef(mod1$lass3)[[2]]
      Dg=coef(mod1$lass3)[[3]]
      txt1=texture(mod1)
      sand=mod1$sand
      silt=mod1$silt
      clay=mod1$clay
      addoutput=rbind(addoutput,data.frame(groupid=aggdata[i],N=N,M=M,Dg=Dg,n=mod1$n,m=mod1$m,n2=mod1$n2,psi=mod1$poresizeindex,b=mod1$b,porosity=mod1$p,D10=mod1$D10,D60=mod1$D60,CU=mod1$CU,K=mod1$K,KB=mod1$KB,r2=gof1$r2,RMSE=gof1$RMSE,NRMSE=gof1$NRMSE,F_value=gof1$F_value,df=gof1$df,p_value=gof1$p_value,AIC=gof1$AIC,BIC=gof1$BIC,r2_adj=gof1$r2_adj,RMSE_adj=gof1$RMSE_adj,NRMSE_adj=gof1$NRMSE_adj,nashSutcliffe=gof1$nashSutcliffe,sand=sand,silt=silt,clay=clay,texture=txt1))
      
      i=i+1
    }
   factor= list(coef=addoutput)
  }
  else{
    
    if(!is.null(C)&!is.null(S)){
      x1=24.547-(0.238*S)-(0.082*C)
      x2=-3.569+(0.081*S)
      x3=0.694-(0.024*S)+(0.048*C)
      Sx1=1/(1+exp(-x1))
      Sx2=1/(1+exp(-x2))
      Sx3=1/(1+exp(-x3))
      n=2.18+0.11*(48.087-(44.954*Sx1)-(1.023*Sx2)-(3.896*Sx3))
      
      m=1-(2/n)
      poreindex=n*m
      psi=poreindex
      n2=(2/(n*m))+3
      b=1/(m*n)
      factor=list(m=m,n=n,n2=n2,poresizeindex=poreindex,b=b,psi=psi)
    }
    else{
      if(is.null(p)){
        mod=lass3 (data=data,p="loam",D=D,fr=fr)
        p=texture(mod)
      }
      
      #use published data
      if(!is.numeric(p)&&!is.null(p)){
        soil=ksat(para="soil",model=p)
        p=soil$para$p
        if(is.na(p)){
          p=soil$para$ths
        }
        #print(p)
      }
      
      if(is.null(data)||is.null(D)||is.null(fr)){
        return(print("Data, Porosity, p, diameter, D, and fraction, fr, of the particles are needed please"))
        break
      }
      if(max(data[[fr]])>2){
        data[[fr]]=data[[fr]]/100
      }
      S=NULL
      p=p
      s=0.001
      while(s<=1.01){
        thiss=((1-p)^s)+(p^(2*s))
        #print(c(round(thiss,3),s))
        if(round(thiss,2)==1){
          S=s
          break
        }
        s=s+0.001
      }
      s=S
      
      if(is.null(S)){
        s=s
      }
      k=((2*s)-1)/(2*s*(1-s))
      
      PDf=paste(fr,"~((1+(Dg/",D,")^N)^(-M))")
      ones=c(N=N,M=M,Dg=Dg)
      lass3<- nlxb(PDf, start = ones, trace = FALSE, data = data)
      
      N=coef(lass3)[1]
      M=coef(lass3)[2]
      Dg=coef(lass3)[3]
      PD=(1+(Dg/data[[D]])^N)^-M
      x=seq(min(data[[D]]),max(data[[D]]),0.01)
      y=(1+(Dg/x)^N)^-M
      
      #print(cor(PD,data["fraction"]))
      PM=(M*N)/(1+M)
      pm=PM*(1+k)^(-1)
      m=(1/pm)*(sqrt(1+(pm^2))-1)
      m=m[[1]]
      n=2/(1-m)
      n2=(2/(m*n))+3
      poreindex=n*m
      psi=poreindex
      
      b=1/(m*n)
      PSDF=TRUE
      D10=NULL
      D60=NULL
      CU=NULL
      K=NULL
      KM=NULL
      KB=NULL
      k1970=NULL
      D10=(((Dg^N)/((1/((0.1)^(1/M)))-1)))^(1/N)
      D60=(((Dg^N)/((1/((0.6)^(1/M)))-1)))^(1/N)
      D16<-(((Dg^N)/((1/((0.16)^(1/M)))-1)))^(1/N)
      D84<-(((Dg^N)/((1/((0.84)^(1/M)))-1)))^(1/N)
      D50<-(((Dg^N)/((1/((0.50)^(1/M)))-1)))^(1/N)
      D5<-(((Dg^N)/((1/((0.05)^(1/M)))-1)))^(1/N)
      D95<-(((Dg^N)/((1/((0.95)^(1/M)))-1)))^(1/N)
      D90<-(((Dg^N)/((1/((0.90)^(1/M)))-1)))^(1/N)
      
      D10=D10[[1]]
      D60=D60[[1]]
      
      phi16<- (-log2(D16))[[1]]
      phi84<- (-log2(D84))[[1]]
      phi50<- (-log2(D50))[[1]]
      phi5<- (-log2(D5))[[1]]
      phi95<- (-log2(D95))[[1]]
      phi90<- (-log2(D90))[[1]]
      phi10<- (-log2(D10))[[1]]
      
      CU=D60/D10
      # Hazen approximation
      K=0.01*(D10^2)
      km=0
      #The Krumbein and Monk equation
      #Geometric (graphic) Mean
      
      
      
      
      #print(cbind(phi5,phi16,phi50,phi84,phi95))
      GMe=((phi16+phi50+phi84)/3)
      #standard deviation
      phisd=((phi84-phi16)/4)+((phi95-phi5)/6.6)
      KM=(760*((GMe)^2)*exp(-1.3*phisd))*(9.87*10^-9)
      
      #Bergs model
      DM=median(data[[D]])
      KB=(5.1*10^-6)*(p^1.5)*DM*exp(-1.38*(phi90-phi10))
      
      r2=cor(PD,data[[fr]])^2
      
      #texture
      two=((1+(Dg/2)^N)^-M)*100
      silt=((1+(Dg/0.05)^N)^-M)*100 
      clay=((1+(Dg/0.002)^N)^-M)*100
      sand=100-silt 
      silt=silt-clay
      clay=clay
     factor<- list(PSDF=length(coef(lass3)),k=k,x=x,y=y,M=M,PM=PM,pm=pm,Dg=Dg,N=N,lass3=lass3,m=m,n=n,predict=PD,PD=PD,s=s,D=data[[D]],
           n2=n2,poresizeindex=poreindex,b=b,fr=data[[fr]],p=p,S=S,C=C,psi=psi,D10=D10,
           D60=D60,CU=CU,K=K,KB=KB,KM=KM,r2=r2,sand=sand[[1]],silt=silt[[1]], clay=clay[[1]])
    }
  }


factor$call<-match.call()

class(factor)<-"lass3"
factor
}

#' @rdname lass3
#' @export
#plot function
plot.lass3<-function(x,main="Particle Size Distribution Function",xlab="Particle Diameter",ylab="Fraction",...)
{
  object=x
mod1=object
if(mod1$PSDF==TRUE||!is.null(mod1$PSDF==TRUE)){
plot(mod1$D,mod1$fr,xlab=xlab,ylab=ylab)
lines(mod1$x,mod1$y)
if(xlab!="2"){
mtext(paste("R^2=",round(cor(mod1$predict,mod1$fr)^2,3)))
}
title(main)
}
}
#' @rdname lass3
#' @export
coef.lass3<-function(object,...)
{
x<-object$lass3
n<-object$n
n2<-object$n2
poresizeindex<-object$poresizeindex
b=object$b
m=object$m
coef=(coef(x))
print(coef)
print(list(n=n,n2=n2,b=b,m=m,poresizeindex=poresizeindex,coef=coef))
}

#' @rdname lass3
#' @export
predict.lass3<-function(object,D,...)
{
x<-object$lass3

coef=(coef(x))
N=coef(x)[[1]]
M=coef(x)[[2]]
Dg=coef(x)[[3]]
PD2=(1+(Dg/D)^N)^-M

return(PD2)
}


#' Prediction of USDA Soil texture classification from PSD object
#'
#' @param object A PSD object from \code{\link{lass3}}
#' @param sand In case object is NULL, percentage of Sand can be supplied
#' @param silt In case object is NULL, percentage of Silt can be supplied
#' @param clay Incase object is NULL, percentage of Clay can be supplied
#' @author George Owusu
#' 
#'
#' @return texture
#' @export
#'
#' @examples
#' file="C:/Users/GeoKings/Documents/George Owusu/UG PhD/Data/soil2.csv"
#' data=read.csv(file)
#' single<- subset(data, ID=="30B20_1")
#' mod=lass3 (data=single,p="sand",D="D",fr="Sand.")
#' texture(mod)
texture<-function(object=NULL,sand=NULL,silt=NULL,clay=NULL){

if(!is.null(object)){
sand=object$sand
silt=object$silt
clay=object$clay
}
texture=NULL


if(((silt + 1.5*clay) < 15)){
texture="Sand"
}

if(((silt + 1.5*clay >= 15) && (silt + 2*clay < 30))){
texture="Loamy sand"
}

if(((clay >= 7 && clay < 20) && (sand > 52) && ((silt + 2*clay) >= 30) || (clay < 7 && silt < 50 && (silt+2*clay)>=30))){
texture="Sandy Loam"
}


##################
if((clay >= 7 && clay < 27) && (silt >= 28 && silt < 50) && (sand <= 52))
{
texture="Loam"
}

if((silt >= 50 && (clay >= 12 && clay < 27)) || ((silt >= 50 && silt < 80) && clay < 12))
{
texture="Silty loam"
}

if((silt >= 80 && clay < 12)){
texture="Silt"
}

if(sand>=20&&sand<=45&&silt>=15&&silt<=52&clay>=27&&clay<=40){
texture="Clay loam"
}

if(((clay >= 20 && clay < 35) && (silt < 28) && (sand > 45)) ){
texture="Sandy clay loam"
}
if(((clay >= 27 && clay < 40) && (sand > 20 && sand <= 45))){
texture="Silty clay loam"
}

if((clay >= 35 && sand > 45)){
texture="Sandy clay"
}

if((clay >= 40 && silt >= 40)){
texture="Silty clay"
}
if((clay >= 40 && sand <= 45 && silt < 40)){
texture="Clay"
}
return(texture)

}
#gof function
#' Estimation of goodness of fit tests with model objects or comparing two models 
#'
#' @param object Function return object
#' @param c comparing model, optional
#' @param r reference model
#' @param P The number of parameters
#' @author George Owusu
#' @return 
#' \itemize{
#' \item{r2:} {  R-square [-]}
#'  \item{r2_adj:} { An adjusted R-square [-]}
#' \item{d:} { Index of agreement[-]}
#' \item{p_value:} { p value}
#' \item{sse:} { Sum of square error}
#' \item{sst:} { Sum of square total}
#' \item{AIC:} { AIC}
#' \item{BIC:} { BIC}
#' \item{RMSE:} { Root Mean Square Error}
#' \item{RMSE_adj:} { adjusted Root Mean Square Error}
#' \item{NRMSE} { Normalised Root Mean Square Error}
#' \item{NRMSE_adj:} { adjusted Normalised Root Mean Square Error}
#' \item{E:} { Nash Sutcliffe Coefficient}
#' \item{Cp:} { statistic (Mallows 1973)}
#' \item{F1:} { F-statistic}
#' \item{AAE:} { average absolute error}
#' \item{ARE:} { average relative error}
#' \item{MBE:} { Mean Bias  error}
#' \item{MPE:} { Mean Percentage error}
#' }
#' @inheritParams philip
#' @references 
#' Mallows, C. L. (1973). Some comments on Cp. Technometrics, 15, 661-675. 
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","soil2.csv",package="vadose"))
#' single<- subset(data, ID=="30B20_1")
#' mod=lass3 (data=single,p="sand",D="D",fr="Sand.")
#' fit=gof(mod)
#' r2=fit$r2
#' AIC=fit$AIC
#' print (fit)
#' 
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' philip1<-philip(data=data,time="time",I="I")
#' sw1<-sw(data=data,time="time",I="I")
#' thisgof=gof(c=gof(sw1),r=gof(philip1))


gof<-function(object=NULL,c=NULL,r=NULL,P=NULL)
{
  
if(!is.null(c)&&is.null(r)&&!is.null(object)){
c=object
r=c
}
if(!is.null(c$sse)&&!is.null(r$sse)&&!is.null(c$sse)){
ssec=c$sse
sser=r$sse
Cp=(ssec/(sser/(c$N-r$parameters)))-(c$N-(2*c$parameters))
cp1=(ssec/(sser/(c$N-r$parameters)))
cp2=(c$N-(2*c$parameters))
Cp=cp1-cp2
f1=(ssec-sser)/sser
f2=(r$N-r$parameters)/((c$N-c$parameters)-(r$N-r$parameters))
f=f1*f2
F1=((ssec-sser)/sser)*((r$N-r$parameters)/((c$N-c$parameters)-(r$N-r$parameters)))
list(F1=F1,Cp=Cp)
}else{

obs=c
est=r
if(is.null(obs)){
obs=object$fr
est=object$predict
}
if(is.null(obs)||is.null(est)){
obs=object$I
est=object$predict
}

if(is.null(obs)||is.null(est)){
  obs=object$theta
  est=object$mod_theta
}

if(class(object)=="BEST"||class(object)=="OFEST"){
  obs=obs[[1]]
}
if(class(object)=="valiantzas"){
  #est=est[1]
  #print(cbind(obs,est))
}
mod=lm(obs~est)
#print(mod)

#print(summary(mod))
#print("p-value")
x=mod
F_value=summary(mod)$fstatistic[[1]]
df=summary(mod)$fstatistic[[3]]
r2=summary(mod)[[8]]
f <- summary(mod)$fstatistic
p_value <- pf(f[1],f[2],f[3],lower.tail=F)
#print(p_value)
attributes(p_value) <- NULL
#Nash-Sutcliffe efficiency rating for model estiamtes https://github.com/USGS-R/rloadest/blob/master/R/nashSutcliffe.R
#nashSutcliffe <- function(obs, est, na.rm=TRUE) {
na.rm=TRUE

  E <- 1 - sum((obs - est)^2, na.rm=na.rm)/sum((obs-mean(obs,na.rm=na.rm))^2, na.rm=na.rm)
#root mean square error
#RMSE=sqrt(mean(mod$residuals^2))
N=nrow(object$data)

AAE=sum(abs(obs-est))/N
ARE=sum(abs((obs-est)/obs))*(100/N)
MBE=sum((est-obs))/N
#MPE=sum(((obs-est)/obs))*(100/N)
MPE=sum(((est-obs)/obs)*100)/N

sse=sum((obs-est)^2)
RMSE=sqrt(sse/N)
NRMSE=RMSE/(max(obs)-min(obs))
#AIC
AIC=AIC(mod)

sse=sum((obs-est)^2)
if(is.null(P)){
P=object$PSDF
if(is.null(P)){
  if(class(object)=="valiantzas"||class(object)=="philip"){
    P=length((object[[1]]))
  }else{
P=length(coef(object[[1]]))
}
}
}
N=length(obs)
rmse_adj=sqrt((sse/(N-P)))
NRMSE_adj=rmse_adj/(max(obs)-min(obs))

adj=N-P
sst=sum((obs-mean(obs))^2)
r2_adj=r2-(1-r2)*(P/(N-P-1))

#AIC=N*log(sse)+(2*P)

AIC=(N*(log(2*pi)+log(sse/adj)+1))+P
#Priestley, M.B. (1981), Spectral Analysis and Time Series, Academic Press. ISBN 0-12-564922-3 (p. 375)
BIC=(N*log(sse/N))+P*log(N)
#parameters=as.data.frame(list(r2=r2,r2_adj=r2_adj,N=N, p=P,F_value=F_value, df=df, p_value=p,AIC=AIC,BIC=BIC,RMSE=RMSE,RMSE_adj=rmse_adj,NRMSE=NRMSE,NRMSE_adj=NRMSE_adj,AAE=AAE,ARE=ARE,MBE=MBE,MPE=MPE,nashSutcliffe=E))
Psquare=sum((est-obs)^2)
obsmean=(abs(obs-mean(obs)))
estmean=(abs(est-mean(obs)))
obs.est=(obsmean+estmean)^2
#d=1-((sum(est-obs)^2)/sum(((abs(est-mean(obs)))+(abs(obs-mean(obs))))^2))
d=1-(Psquare/sum(obs.est))
parameters=list(r2=r2,r2_adj=r2_adj,N=N, parameters=P,d=d,NSE=E,F_value=F_value, df=df, p_value=p_value,sse=sse,sst=sst,AIC=AIC,BIC=BIC,RMSE=RMSE,AAE=AAE,ARE=ARE,MBE=MBE,MPE=MPE,RMSE_adj=rmse_adj,NRMSE=NRMSE,NRMSE_adj=NRMSE_adj,nashSutcliffe=E,E=E)

if(is.null(object$PSDF)){
  parameters=as.data.frame(parameters)
}
parameters
#parameters=as.data.frame(list(r2=r2,r2_adj=r2_adj,d=d, p_value=p,sse=sse,sst=sst,AIC=AIC,BIC=BIC,RMSE=RMSE,RMSE_adj=rmse_adj,NRMSE=NRMSE,NRMSE_adj=NRMSE_adj,nashSutcliffe=E,
#                             b=coef(mod)[[2]],AAE=AAE,ARE=ARE,MBE=MBE,MPE=MPE,r=cor(obs,est), intercept=coef(mod)[[1]],N=N, parameters=P,F_value=F_value, df=df))
#print(parameters)-
}
}

