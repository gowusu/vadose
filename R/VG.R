#input function ############################################
#' @title R estimation of van Genuchten water retention curve and hydraulic conductivity with internal Ks
#' @description This function estimate and optimise van Genuchten water retention and hydraulic conductivity.
#' @param data dataframe that contains pressure head "h" and the corresponding water content "theta"
#' @param h character. The name of the pressure head as in the data parameter above.
#' @param theta character. The name of the water content as in the data parameter above
#' @inheritParams BEST
#' @param Ks saturated hydraulic conductivity
#' @param alp van Genuchten water retention alpha parameter
#' @param para character. This is useful if the function is used for 
#' estimation instead of optimisation. It can be combined with Ks parameter to retrieve
#' hydraulic parameters from Carsel&Parrish (1988). If para is set to "soil" and Ks is set to
#' "Clay", the parameters of clay soil from Carsel&Parrish (1988) will be used for estimation.
#' see example below.
#' @param layout plot layout
#' @param kcol colour of hydraulic conductivity
#' @param hcol colour of pressure head
#' @author George Owusu
#'  
#' @inheritParams lass3
#' @references 
#' \itemize{
#' \item{}{van Genuchten, M. T. (1980). A Closed-form Equation for Predicting the
#'  Hydraulic Conductivity of Unsaturated Soils1. Soil Sci. Soc. Am. J., 44(5), 892-898. 
#'  doi: 10.2136/sssaj1980.03615995004400050002x}
#' \item{}{Carsel, R. F., & Parrish, R. S. (1988). Developing joint probability distributions of 
#' soil water retention characteristics. Water Resour. Res., 24, 755-769.}
#' }
#'
#' @return optimised parameters
#' \itemize{
#' \item{Ks} { saturated hydraulic conductivity [LT^-1] }
#' \item{mod_theta,mod_thetal,mod_thetay,mod_thetab:} { The estimated soil water content at 
#' different metric potentials}
#' \item{theta:} { The observed soil water content at different metric potentials  }
#' \item{K:} { The unsaturated hydraulic conductivity [LT^-1] }
#' \item{se:} { effective saturation}
#' \item{thr} { optimised residual water content[LTLT^-1] }
#' \item{ths} { optimised saturated water content[LTLT^-1] }
#' \item{alp} { optimised alpha value}
#' \item{n} {optimised n shaping parameter }
#' \item{Kr} { The ratio of K and Ks }
#' } 
#' @export
#'
#' @examples
#' \dontrun{
#'  #optimisation with single data
#' datartc=read.csv(system.file("ext","sys","retentionVG.csv",package="vadose"))
#' modrtc<-vg(data=datartc,h="h",theta="theta",thr=0.1, ths=0.1, alp=0.1, n=1)
#' plot(modrtc)
#' 
#' #optimisation with group data isric
#' 
#' data=read.csv(system.file("ext","sys","isric2.csv",package="vadose"))
#' #optimisation with group data
#' #used public initials and multiple Ks
#' modisrc<-vg(data=data,h="x",theta="y",m="b",thr=0.1, ths=0.1, alp=0.1, n=1,group="Sample",
#' Ks=c("Sand","Clay","Silt","silty clay loam"),para="soil")
#' modisrc<-vg(data=data,h="x",theta="y",m="b",thr=0.1, ths=0.3, alp=0.01, n=2,group="Sample")
#' plot(modisrc)
#' mod=vg(h=200)
#' }
#generic function
vg<-function(data=NULL,h,theta=NULL,thr=0.1, ths=0.1, alp=0.1, n=1,m="m",Ks,para=NULL,group=NULL) UseMethod ("vg")

#default function
#' @rdname vg
#' @export
vg.default<-function(data=NULL,h,theta=NULL,thr=0.1, ths=0.1, alp=0.1, n=1,m="m",Ks="Nasta",para=NULL,group=NULL)
{
  #stop warning from displaying
  options(warn=-1)
  #compute retention whene data is null but other parameters are provided
  if(data==""||is.null(data))
  {
    #use published parameters
    if(!is.null(para)){
      if(para=="soil"||para=="Soil"){
        thismod=ksat(para=para,model=Ks)
        thr=thismod$para$thr
        ths=thismod$para$ths
        alp=thismod$para$alp
        n=thismod$para$n
        
        if(!is.numeric(Ks)){
          Ks=thismod$para$Ks
        }
      }
    }
    if(is.null(Ks)){
      return(print("Please privide the correct Ks value. Check the soil texture name"))
    }
    thism=m
    PSDF=4
    
    if(is.numeric(m)){
      m=m
      PSDF=5
    }
    #mualem condition
    if(m=="m"||m=="Mualem"||m=="mualem"){
      m=1-(1/n)
    }
    #burdine condition
    if(m=="b"||m=="Burdine"||m=="burdine"){
      m=1-(2/n)
    }
    
    #predict water content
    predict2=thr+(ths-thr)/((1+(alp*h)^n)^(-m))
    
    #h2=seq(min(h), max(h), by = 0.05)
    h2=h
    predict=thr+(ths-thr)/((1+(alp*h2)^n)^(-m))
    
    se=(predict-thr)/(ths-thr)
    #print(se)
    Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)^2
    if(thism=="b"||thism=="Burdine"||thism=="burdine"){
      Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)
    }
    
    if(is.numeric(Ks)){
      Ks=Ks
    }
    else
    {
      th33=thr+(ths-thr)/((1+(alp*33)^n)^(m))###########
      th1500=thr+(ths-thr)/((1+(alp*1500)^n)^(m))
      if(length(Ks>1)){
        Ks=ksat(ths=ths,thr=thr,alpha=alp,hb=1/alp,th33=th33,th1500=th1500,C=NULL,f=ths,para=para,n=n,model=Ks[1])
      }
      else{
        Ks=ksat(ths=ths,thr=thr,alpha=alp,hb=1/alp,th33=th33,th1500=th1500,C=NULL,f=ths,para=para,n=n,model=Ks)
      }
      Ks=Ks$Ks
    }
    K=Ks*Kr
    return(list(PSDF=PSDF,theta=theta,predict=predict,h2=h2,predict2=predict2,thr=thr,ths=ths,alp=alp,n=n,h=h,h2=h2,K=K,Kr=Kr,Ks=Ks,se=se,water=predict))
  }
  
  
  #remove zeros
  data[data==0]<-0.0001
  
  #decalare the group data incase group variable is not null
  addoutput=NULL
  
  #use initial values of Carsel and Parrish [1988]
  if(!is.null(para)){
    if(para=="soil"||para=="Soil"){
      thismod=ksat(para=para,model=Ks)
      thr=thismod$para$thr
      ths=thismod$para$ths
      alp=thismod$para$alp
      n=thismod$para$n
    }
  }
  #set parameters of optimsation functions######################
  if(m==""){
    return(print("please indicate m; m can take a number or b or m"))
  }
  ones <- c(thr=thr, ths=ths, alp=alp, n=n) # all ones start
  thism=m
  if(is.numeric(m)){
    ones <- c(thr=thr, ths=ths, alp=alp, n=n,m=m) # all ones start
    m="m"
  }
  #mualem condition
  if(m=="m"||m=="Mualem"||m=="mualem"){
    m="1-(1/n)"
  }
  #burdine condition
  if(m=="b"||m=="Burdine"||m=="burdine"){
    m="1-(2/n)"
  }
  
  
  
  rate="yes"
  
  vgh=paste(theta,"~(thr+(ths-thr)/((1+(alp*",h,")^n)^(",m,")))")
  
  rate="no"
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    vg<- nlxb(vgh, start = ones, trace = FALSE, data = data)
    print(vg)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame###################################
    addoutput=data.frame(groupid=factor(),theta=numeric(),thetaEST=numeric(),h=numeric(),thr=numeric(),
                         ths=numeric(),alpha=numeric(),n=numeric(),m=numeric(),se=numeric(),K=numeric(),Kr=numeric(),Ks=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      #use initial values of Carsel and Parrish [1988]
      if(length(Ks)>1){
        if(!is.null(para)){
          if(para=="soil"||para=="Soil"){
            thismod=ksat(para=para,model=Ks[i])
            thr=thismod$para$thr
            ths=thismod$para$ths
            alp=thismod$para$alp
            n=thismod$para$n
            ones <- c(thr=thr, ths=ths, alp=alp, n=n) # all ones start
          }
        }
      }
      vg<- nlxb(vgh, start = ones, trace = FALSE, data = single)
      print(vg)
      print("....................................................................................")
      #vg paramters #################################
      groupdata=single[[group]]
      #prediction cumulative equation ##############################
      thr=coef(vg)[1]
      ths=coef(vg)[2]
      alp=coef(vg)[3]
      n=coef(vg)[4]
      theta3=single[[theta]]
      
      h3=single[[h]]
      
      #Ks=350.2
      if(is.numeric(thism)){
        m=coef(vg)[5]
      }
      
      if(thism=="m"||thism=="Mualem"||thism=="mualem"){
        m=1-(1/n)
      }
      
      if(thism=="b"||thism=="Burdine"||thism=="burdine"){
        m=1-(2/n)
      }
      
      ##########################################################
      h2=seq(min(h3), max(h3), by = 0.05)
      predict=thr+(ths-thr)/((1+(alp*h2)^n)^(m))###########
      predict2=thr+(ths-thr)/((1+(alp*h3)^n)^(m))###########
      
      se=(predict2-thr)/(ths-thr)
      se1=(predict-thr)/(ths-thr)
      
      Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)^2
      Kr1=(se1^0.5)*(1-(1-(se1^(1/m)))^m)^2
      
      if(thism=="b"||thism=="Burdine"||thism=="burdine"){
        Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)
        #De=((1-m)*se^(0.5-(m+1)/(2*m))/(2*alp*m(ths-thr)))*((1-se^(1/m))^-((m+1)/2)-((1-se^(1/m))^((m-1)/2)))
        Kr1=(se1^0.5)*(1-(1-(se1^(1/m)))^m)
      }
      
      
      
      if(is.numeric(Ks)){
        Ks=Ks
      }
      else
      {
        th33=thr+(ths-thr)/((1+(alp*33)^n)^(m))###########
        th1500=thr+(ths-thr)/((1+(alp*1500)^n)^(m))
        Ks=ksat(ths=ths,thr=thr,alpha=alp,hb=1/alp,th33=th33,th1500=th1500,C=NULL,f=ths,para=para,n=n,model=Ks)
        Ks=Ks$Ks
      }
      K=Ks*Kr
      K1=Ks*Kr1
      
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,theta=theta3,thetaEST=predict2,h=h3,thr=thr,ths=ths,alpha=alp,n=n,m=m,se=se,K=K,Ks=Ks,Kr=Kr))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    thr=coef(vg)[1]
    ths=coef(vg)[2]
    alp=coef(vg)[3]
    n=coef(vg)[4]
    theta=data[[theta]]
    h=data[[h]]
    
    #Ks=350.2
    if(is.numeric(thism)){
      m=coef(vg)[5]
    }
    
    if(thism=="m"||thism=="Mualem"||thism=="mualem"){
      m=1-(1/n)
    }
    
    if(thism=="b"||thism=="Burdine"||thism=="burdine"){
      m=1-(2/n)
    }
    ##########################################################
    h2=seq(min(h), max(h), by = 0.05)
    predict=thr+(ths-thr)/((1+(alp*h2)^n)^(m))###########
    predict2=thr+(ths-thr)/((1+(alp*h)^n)^(m))###########
    
    se=(predict-thr)/(ths-thr)
    Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)^2
    
    if(thism=="b"||thism=="Burdine"||thism=="burdine"){
      Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)
      #De=((1-m)*se^(0.5-(m+1)/(2*m))/(2*alp*m(ths-thr)))*((1-se^(1/m))^-((m+1)/2)-((1-se^(1/m))^((m-1)/2)))
    }
    
    if(is.numeric(Ks)){
      Ks=Ks
    }
    else
    {
      th33=thr+(ths-thr)/((1+(alp*33)^n)^(m))###########
      th1500=thr+(ths-thr)/((1+(alp*1500)^n)^(m))
      Ks=ksat(ths=ths,thr=thr,alpha=alp,hb=1/alp,th33=th33,th1500=th1500,C=NULL,f=ths,para=para,n=n,model=Ks)
      Ks=Ks$Ks
    }
    K=Ks*Kr
    ###############################
    
    
    
    #par(mar = c(5, 4, 4, 4) + 0.3)
    #plot(theta, h,xlab="Water Content") # first plot
    #lines(predict,h2)
    #par(new = TRUE)
    #plot(predict, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
    #axis(side=4, at = pretty(range(K)))
    #mtext("K", side=4, line=3)
  }
  else
  {
    ################################################
    thr=addoutput$thr
    ths=addoutput$ths
    alp=addoutput$alpha
    n=addoutput$n
    predict=addoutput$thetaEST
    h=addoutput$h
    Ks=addoutput$Ks
    ##############################################
  }
  if(!is.null(group)){
    #print(addoutput)
  }
  #return varibales ########################################
  factor<-list(PSDF=length(coef(vg)),vg=vg,data=data,theta=theta,thr=thr,ths=ths,alp=alp,n=n,h=h,K=K,Kr=Kr,Ks=Ks,se=se,water=predict,h2=h2,group=group,
       predict=predict,formular=vgh,addoutput=addoutput,mod_theta=predict,rate="no",predict2=predict2)

factor$call<-match.call()

class(factor)<-"vg"
factor
}

#predict function#########################
#' @rdname vg
#' @export
predict.vg<-function(object,h=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$h,x$predict2))
names(predict)=c("h","watercontent")

if(!is.null(h)&object$rate=="yes"){
predict=0.5*object$S*(h^-0.5)+object$A
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(h)&object$rate!="yes"){
predict=object$thr+(object$ths-object$thr)/((1+(object$alp*h)^object$n)^(1-(1/object$n)))
predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#cumulative###################################
if(!is.null(h)&object$rate!="yes"){
#predict=object$addoutput$S*(h^0.5)+object$addoutput$A*h
predict=object$addoutput$thr+(object$addoutput$ths-object$addoutput$thr)/((1+(object$addoutput$alp*h)^object$addoutput$n)^(1-(1/object$addoutput$n)))

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}

print(predict)
}
}
#' @rdname vg
#' @export
#plot function
plot.vg<-function(x,main=NULL,xlab="Water Content",
ylab="Tension Head,h(cm)",ylab2="Hydraulic Conductivity, K (cm/day)",
layout=NULL,kcol="green",hcol="black",hlog=NULL,klog=NULL,...)
{
object <- x
if(is.null(object$group)){
#par(mfrow=c(2,2),mar=c(4, 4, 4, 4))

K=object$K
h=object$h
h2=object$h2
predict=object$predict
theta=object$theta
#plot(theta,h,xlab="Water Content")
#lines(predict,h2,xlab="Water Content")
#plot(predict,K, type = "l",xlab="Water Content")


###########################
if(!is.null(theta)){
plot(theta, h,xlab=xlab,ylab=ylab,...) # first plot
lines(predict,h2)

par(new = TRUE)
plot(predict, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",...)

axis(side=4, at = pretty(range(K)))
mtext(ylab2, side=4, line=3)
index1=length(theta)
index2=length(theta)-1
leb=(theta[index2]-theta[index1])
text(max(theta)-leb,max(K),"K")
text(min(theta)+leb,max(K),"h")
}
if(is.null(theta)){
ymax <- max(K, na.rm=TRUE)
if(!is.null(hlog)){
plot(predict,h2,type="l",xlab=xlab,ylab=ylab,col=hcol,log="y",...)
}else{
plot(predict,h2,type="l",xlab=xlab,ylab=ylab,col=hcol,...)
}
par(new = TRUE)
if(!is.null(klog)){
plot(predict, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=kcol,log="y",...)
}else{
plot(predict, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=kcol,...)
}
axis(side=4, ylim=c(0, ymax),col.axis=kcol )
mtext(ylab2, side=4, line=3,col=kcol)
}
}

else
{
aggdata =row.names(table(object$addoutput$groupid))
data=object$addoutput

if(is.null(layout)){
lengthD=length(aggdata)
 #dev.new(width=15, height=10)
#dev.new(width=12, height=11)
}
op=NULL
if(lengthD==2){
op=par(mfrow=c(1,2),mar=c(2, 2, 2, 2))
}

if(lengthD==3){
op=par(mfrow=c(1,3),mar=c(2, 2, 2, 2))
}

if(lengthD==4){
op=par(mfrow=c(2,2),mar=c(4, 4, 2, 2))
}

if(lengthD==5){
op=par(mfrow=c(2,3),mar=c(2, 2, 2, 2))
}

if(lengthD==6){
op=par(mfrow=c(2,3),mar=c(2, 2, 2, 2))
}
if(lengthD==9){
op=par(mfrow=c(3,3),mar=c(2, 2, 2, 2))
}


if(lengthD>6){
op=par(mfrow=c(round(lengthD/2),round(lengthD/2)),mar=c(2, 2, 2, 2))
}

#print(length(aggdata))

#matrix plot
i=1
while(i<=length(aggdata)){
#label=aggdata[i]
#print (label)
single=data[data["groupid"]==aggdata[i],]
h=single$h
thetaEST=single$thetaEST
theta=single$theta
K=single$K
r2=cor(thetaEST,theta)^2
title=NULL
if(is.null(main)){
title=paste(aggdata[i],"(R2=",round(r2,3),")")
}
############################################################
thr=single$thr
ths=single$ths
alp=single$alp
n=single$n
Ks=single$Ks
se=single$se
Kr=single$Kr
K=single$K
m=single$m

h2=seq(min(h), max(h), by = 0.05)
predict2=thr+(ths-thr)/((1+(alp*h2)^n)^(m))###########
predict=thr+(ths-thr)/((1+(alp*h)^n)^(m))###########
se=(predict-thr)/(ths-thr)
se2=(predict2-thr)/(ths-thr)
Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)^2
Kr2=(se2^0.5)*(1-(1-(se2^(1/m)))^m)^2
K2=Ks*Kr2
log2="jj"
par(mar = c(5, 4, 4, 4) + 0.3)
if(log2=="h"){
log2="y"
}
if(log2=="K"||log2=="hK"||log2=="Kh"||log2=="kh"||log2=="hk"){
log2="y"
}
if(log2=="theta"||log2=="water"){
log2="x"
}

plot(theta, h,xlab="Water Content",ylab=ylab,main=main,cex.lab=0.7,...) # first plot
lines(predict2,h2)
par(new = TRUE)
plot(predict2, K2, type = "l",axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(K2)))
mtext(ylab2, side=4, line=3,cex=0.7)
title(title)
index1=length(theta)
index2=length(theta)-1
leb=(theta[index2]-theta[index1])*2
#print(leb)
text(max(theta)-leb,max(K),"K")
text(min(theta)+leb,max(K),"h")

i=i+1
}
par(op)
}
}

#summary function
#' @rdname vg
#' @export
summary.vg<-function(object,...)
{


x<-object$vg

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{
coef=aggregate(cbind(thr,ths,alpha, n,m,Ks) ~ groupid, data = object$addoutput, mean)
print(coef)
}

}

#print function
#' @rdname vg
#' @export
print.vg<-function(x,...)
{
object<-x
if(is.null(object$group)){
x<-object$vg
print((x))
}
else
{
coef=aggregate(cbind(thr,ths,alpha, n,m,Ks) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}

#coef function
#' @rdname vg
coef.vg<-function(object,...)
{
x<-object$vg
if(is.null(object$group)){
coef=(coef(x))

}
else
{
coef=aggregate(cbind(thr,ths,alpha, n,m,Ks) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

#statistics of error


#optimisation with group data
#modrtc<-vg(data=datartc,h="h",theta="theta",thr=0.1, ths=0.1, alp=0.1, n=1,group="group")
#optimisation with group data
#used public initials and multiple Ks
#modisrc<-vg(data=data,h="x",theta="y",m="b",thr=0.1, ths=0.1, alp=0.1, n=1,group="Sample",Ks=c("Sand","Clay","Silt","silty clay loam"),para="soil")
#modisrc<-vg(data=data,h="x",theta="y",m="b",thr=0.1, ths=0.3, alp=0.01, n=2,group="Sample")
#plot(modisrc)

#optimisation with single data
#modrtc<-vg(data=datartc,h="h",theta="theta",thr=0.1, ths=0.1, alp=0.1, n=1)
#plot(modrtc,type=c("K","h2"))
#x <-c(0,log(10),log(31),log(100),log(200),log(500),log(2500),log(15000))# get log (x)
#x <-c(0,10,31,100,200,500,2500,15000)# get log (x)
#y <-c(0.362,0.353,0.294,0.230,0.206,0.182,0.148,0.136)
#Brookdata <-data.frame(x,y)# make data frame
#optimisation with single with Carsel and Parrish [1988] initial values
#mod=vg(data=Brookdata,h="x",theta="y",m="b",thr=0.1, ths=0.5, alp=0.05, n=2,Ks="clay",para="soil")
#predicting h with Carsel and Parrish [1988] initial values
#mod=vg(data=NULL,h=h,theta="y",thr=0.1, ths=0.5, alp=0.05, n=2,Ks="sand",para="soil")
#h=c(100,650,1000,5000,10000,15000)
#predicting h with Carsel and Parrish [1988] initial values
#mod=vg(h=h,Ks="SILTY CLAY",para="soil")
#mod=vg(h=assin_breko$x,Ks="silty clay loam",para="soil")
#plot(mod)
#plot(mod)
#mod<-vg(data=data,h="x",theta="y",thr=0.1, ths=0.1, alp=0.1, n=1,m="b",group="Sample",Ks="sr")
#plot(mod)
#mod=vg(h=200)