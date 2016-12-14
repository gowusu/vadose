#input function ############################################
#' @title Brutsaert infiltration parameter optmisation in R
#' @description This function optimises Brutsaert, W. (1977)
#' cumulative infiltration (I): B, Ks and S. 

#' 
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param B brutsaert parameter
#' @param para numeric. Whether 2 or 3 parameter brutsaert model should be implemented.
#'
#' @return 
#' \itemize{
#' \item{Ks:} {  optimised saturated hydraulic conductivity [L/T]}
#' \item{S:} { optimised sorptivity [LT^-0.5]}
#' \item{B:} { optmised B parameter}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Brutsaert, W. (1977). Vertical infiltration in dry soil. Water Resour. Res., 13, 363-368. 
#' @export
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' brutsaert1<-brutsaert(data=data,time="time",I="I")
#' #print(gof(brutsaert1))
#' plot(brutsaert1)
#' predict(brutsaert1)
#' coef(brutsaert1)
#' #gof(brutsaert1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' brutsaertg<-brutsaert(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(brutsaertg)
#' @rdname brutsaert
brutsaert <- function(data=NULL,time,I,S=0.1,Ks=0.1,B=1,para=2,group=NULL) UseMethod ("brutsaert")
#' @export
#' @rdname brutsaert
brutsaert.default <- function(data=NULL,time,I,S=0.1,Ks=0.1,B=1,para=2,group=NULL)
{
#stop warning from displaying
options(warn=-1)

if(data==""||is.null(data))
{
predict=((Ks*time)+((S^2)/(B*Ks))*(1-(1/((1+(B*Ks*sqrt(time)))/S))))
return(predict)
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
  }
#decalare the group data incase group variable is not null
addoutput=NULL
#set parameters of optimsation functions######################
#determine whether it is rate or cumulative data
rate="yes"
if(data[[I]][length(data[[I]])]>data[[I]][1])
{
#cumulative function ###############################
brutsaertF<-paste(I,"~((Ks*time)+((S^2)/(B*Ks))*(1-(1/((1+(B*Ks*sqrt(time)))/S))))")
ones <- c(S = S, Ks = Ks) # all ones start
brutsaertF<-paste(I,"~((Ks*",time,")+((S^2)/(",B,"*Ks))*(1-(1/((1+(",B,"*Ks*sqrt(",time,")))/S))))")
if(para==3){
ones <- c(S = S, Ks = Ks,B=B) # all ones start
brutsaertF<-paste(I,"~((Ks*",time,")+((S^2)/(B*Ks))*(1-(1/((1+(B*Ks*sqrt(",time,")))/S))))")
}
rate="no"
}
else
{
#rate function #######################################

#brutsaertF<-paste(I,"~(0.5*S","*",time,"^(-0.5))+ Ks")

}
#check for grouped data and execute the optimisation function
if(is.null(group)){
brutsaert<- nlxb(brutsaertF, start = ones, trace = FALSE, data = data)
print(brutsaert)

}
else
{
# write the group function
aggdata =row.names(table(data[group]))
#create group data frame
addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),Ks=numeric(),S=numeric())
i=1
while(i<=length(aggdata)){
print(paste("Group Number:",aggdata[i]))
single=data[data[group]==aggdata[i],]
#group function
brutsaert<- nlxb(brutsaertF, start = ones, trace = FALSE, data = single)
print(brutsaert)
print("....................................................................................")
#brutsaert paramters #################################
S=coef(brutsaert)[1]
Ks=coef(brutsaert)[2]
groupdata=single[[group]]
time2=single[[time]]
#prediction cumulative equation ##############################
if(para==3){
B=coef(brutsaert)[3]##########################################
}
predict2=((Ks*time2)+((S^2)/(B*Ks))*(1-(1/((1+(B*Ks*sqrt(time2)))/S))))


if(rate=="yes"){
#predict rate###############################################
predict2=0.5*S*(time2^-0.5)+Ks
}
addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,Ks=Ks,S=S))

i=i+1
}
}

#equations for ungrouped data ############################
if(is.null(group)){
#prediction################################################
S=coef(brutsaert)[1]
time2=data[[time]]
I=data[[I]]
time=data[[time]]
###########################################################
Ks=coef(brutsaert)[2]##########################################
######################
if(para==3){
B=coef(brutsaert)[3]##########################################
}
predict=((Ks*time)+((S^2)/(B*Ks))*(1-(1/((1+(B*Ks*sqrt(time)))/S))))


if(rate=="yes"){
#####################################################
predict=0.5*S*(time2^-0.5)+Ks
}
}
else
{
################################################
S=addoutput$S
I=data[[I]]
time=data[[time]]
##############################################
Ks=addoutput$Ks
predict=addoutput$predict
time=addoutput$time

}
#return varibales ########################################
factor <- list(brutsaert=brutsaert,data=data,time=time,I=I,S=S,Ks=Ks,B=B,para=para,group=group,
predict=predict,rate=rate,formular=brutsaertF,addoutput=addoutput,output=addoutput)
factor$call<-match.call()

class(factor)<-"brutsaert"
factor
}

#generic function
 
#' @rdname brutsaert
#' @export
#predict function#########################
predict.brutsaert<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$Ks
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$Ks*time)+((object$S^2)/(object$B*object$Ks))*(1-(1/((1+(object$B*object$Ks*sqrt(time)))/object$S))))
predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$Ks

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$addoutput$Ks*time)+((object$addoutput$S^2)/(object$addoutput$B*object$addoutput$Ks))*(1-(1/((1+(object$addoutput$B*object$addoutput$Ks*sqrt(time)))/object$addoutput$S))))
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}

print(predict)
}
}
#' @export 
#' @rdname brutsaert
#plot function
plot.brutsaert<-function(x,xlab="Time(Minutes)",ylab="Cumulative (mm)",main=NULL,layout=NULL,...)
{

object<-x
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
#' @rdname brutsaert

#summary function
summary.brutsaert<-function(object)
{


x<-object$brutsaert

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{
coef=aggregate(cbind(Ks, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}

}
#' @rdname brutsaert

#print function
print.brutsaert<-function(x,...)
{
object=x
if(is.null(object$group)){
x<-object$brutsaert
print((x))
}
else
{
coef=aggregate(cbind(Ks, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}

#' @export 
#' @rdname brutsaert
#coef function
coef.brutsaert<-function(object,...)
{
x<-object$brutsaert
if(is.null(object$group)){
coef=(coef(x))

}
else
{
coef=aggregate(cbind(Ks, S) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

#statistics of error

