#input function ############################################
#' @title Swartzendruber infiltration parameter optimisation in R
#' 
#' @description This function optimises Swartzendruber (1987)
#' cumulative infiltration (I) or infiltration rate (i) parameters: As, Ks and S. 
#'
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param As an optimised parameter 
#'
#' @return 
#' \itemize{
#' \item{As:} { optimised As parameter}
#' \item{Ks:} {  optimised saturated hydraulic conductivity[L/T]}
#' \item{S:} { optimised sorptivity [LT^-0.5]}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Swartzendruber, D. (1987). A quasi-solution of Richards equation for 
#' the downward infiltration of water into soil. Water Resour Res, 23, 809-817. 
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' sw1<-sw(data=data,time="time",I="I")
#' print(gof(sw1))
#' plot(sw1)
#' predict(sw1)
#' coef(sw1)
#' gof(sw1)
#' 
#' #infiltration rate
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' swr<-sw(data=assin_breko,time="minutes",I="cm.hr")
#' 
#' #group simulation
#' swg<-sw(data=data,time="minutes",I="cm.hr",group="ID")
#' coef(swg)
#generic function
sw<-function(data,time,I,As=0.1,S=0.1,Ks=0.1,group=NULL) UseMethod ("sw")############
#' @rdname sw
#' @export
sw.default <- function(data,time,I,As=0.1,S=0.1,Ks=0.1,group=NULL)
{
  
#https://books.google.com.gh/books?id=wpHF2M66QVEC&pg=PR13&lpg=PR13&dq=Jury++Soil+physics&source=bl&ots=Rvck5EDCjn&sig=gcCvVu8SMcSzhnDRQi8Cd-r1_5s&hl=en&sa=X&ei=vX42VdW4HY_haODsgNAK&ved=0CDQQ6AEwAw#v=onepage&q=Jury%20%20Soil%20physics&f=false
#reference
#stop warning from displaying
  options(warn=-1)
  
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
  }
#decalare the group data incase group variable is not null
addoutput=NULL
#set parameters of optimsation functions######################
ones <- c(As=As, S = S,Ks=Ks) # all ones start
#determine whether it is rate or cumulative data
rate="yes"

if(data[[I]][length(data[[I]])]>data[[I]][1])
{
#cumulative function ###############################
#swF<-paste(I,"~1-exp(-As*sqrt(",time,"))*(S/As)+(Ks*",time,")")
swF<-paste(I,"~(1-exp(-As*sqrt(",time,")))*(S/As)+(Ks*",time,")")

rate="no"
}
else
{
#rate function #######################################
swF<-paste(I,"~((Ks*",time,")+(S/As)*(1-exp(-As*(",time,"^0.5))))")
}
#check for grouped data and execute the optimisation function
if(is.null(group)){
sw<- nlxb(swF, start = ones, trace = FALSE, data = data)
#print(sw)

}
else
{
# write the group function
aggdata =row.names(table(data[group]))
#create group data frame###############################################
addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),As=numeric(),S=numeric(),Ks=numeric())
i=1
while(i<=length(aggdata)){
print(paste("Group Number:",aggdata[i]))
single=data[data[group]==aggdata[i],]
#group function
sw<- nlxb(swF, start = ones, trace = FALSE, data = single)
print(sw)
print("....................................................................................")
#sw paramters #################################
As=coef(sw)[1]
S=coef(sw)[2]
Ks=coef(sw)[3]
groupdata=single[[group]]
time2=single[[time]]
#prediction cumulative equation ##############################
predict2=(1-exp(-As*sqrt(time2)))*(S/As)+(Ks*time2)
if(rate=="yes"){
#predict rate###############################################
predict2=((Ks*time2)+(S/As)*(1-exp(-As*(time2^0.5))))

}
addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,As=As,S=S,Ks=Ks))

i=i+1
}
}

#equations for ungrouped data ############################
if(is.null(group)){
#prediction################################################
As=coef(sw)[1]
time2=data[[time]]
I=data[[I]]
time=data[[time]]
###########################################################
S=coef(sw)[2]
Ks=coef(sw)[3]##########################################
######################
predict=(1-exp(-As*sqrt(time2)))*(S/As)+(Ks*time2)



if(rate=="yes"){
#####################################################
predict=((Ks*time2)+(S/As)*(1-exp(-As*(time2^0.5))))
}
}
else
{
################################################
As=addoutput$As
I=data[[I]]
time=data[[time]]
##############################################
S=addoutput$S
Ks=addoutput$Ks
predict=addoutput$predict
time=addoutput$time
}
#return varibales ########################################
factor <- list(sw=sw,data=data,time=time,I=I,As=As,S=S,Ks=Ks,group=group,
predict=predict,rate=rate,formular=swF,addoutput=addoutput,output=addoutput)
factor$call<-match.call()

class(factor)<-"sw"
factor
}


#' @rdname sw
#' @export
#predict function#########################
predict.sw<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")
#rate
if(!is.null(time)&object$rate=="yes"){####################
predict=(1-exp(-object$As*sqrt(time)))*(object$S/object$As)+(object$Ks*time)
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=(1-exp(-object$As*sqrt(time)))*(object$S/object$As)+(object$Ks*time)

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=(1-exp(-object$As*sqrt(time)))*(object$S/object$As)+(object$Ks*time)

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=(1-exp(-object$As*sqrt(time)))*(object$S/object$As)+(object$Ks*time)
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}

print(predict)
}
}
#' @rdname sw
#' @export
#plot function
plot.sw<-function(x,xlab="Time(Minutes)",ylab="Cumulative (mm)",main=NULL,layout=NULL,...)
{

object <- x
S=x$S
time=x$time
I=x$I
rate=x$rate
op=par()
if(rate=="yes"&ylab=="Cumulative (cm)"){
ylab="rate(cm/mins)"
}

if(is.null(x$group)){
r2=cor(I,x$predict)^2
if(is.null(main)){

main=paste("R2=",round(r2,4))
}
plot(time,I,xlab=xlab,ylab=ylab,main=main,...)
#plot(time,I)

par(new=T)
predict=x$predict
lines(time,predict,col="red")
}
else
{
aggdata =row.names(table(object$addoutput$groupid))
data=object$addoutput

if(is.null(layout)){
lengthD=length(aggdata)
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
op=par(mfrow=c(3,3),mar=c(2, 2, 2, 2))
}

if(lengthD>6){
op=par(mfrow=c(round(lengthD/2),round(lengthD/2)),mar=c(2, 2, 2, 2))
}



}
#print(length(aggdata))

#matrix plot
i=1
while(i<=length(aggdata)){
#label=aggdata[i]
#print (label)
single=data[data["groupid"]==aggdata[i],]
I=single$observed
predict=single$predict
time=single$time
r2=cor(I,predict)^2
title=NULL
if(is.null(main)){
title=paste(aggdata[i],"(R2=",round(r2,3),")")
}
plot(time,I,main=main,xlab=xlab,ylab=ylab,...)
title(title)
par(new=T)
lines(time,predict,col="red")

i=i+1
}
par(op)

}
}

#summary function
summary.sw<-function(object)
{


x<-object$sw

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{########################
coef=aggregate(cbind(As, S,Ks) ~ groupid, data = object$addoutput, mean)##############
print(coef)
}

}
#' @rdname sw
#' @export
#print function
print.sw<-function(x,...)
{
  object <- x
  
if(is.null(object$group)){
x<-object$sw
print((x))
}
else
{######################################################
coef=aggregate(cbind(As, S,Ks) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @rdname sw
#' @export
#coef function
coef.sw<-function(object,...)
{
x<-object$sw
if(is.null(object$group)){
coef=(coef(x))

}
else################################################
{
coef=aggregate(cbind(As, S,Ks) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

#statistics of error

