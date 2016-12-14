#input function ############################################
#' @title Valiantzas infiltration parameter optimization in R
#' 
#' @description This function optimises Valiantzas (2010)
#' cumulative infiltration (I)  parameters: Ks, S and ratio between philip A and Ks. 
#'
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param ratio a parameter expressing relationship between philip A and Ks
#' @return 
#' \itemize{
#' \item{Ks:} {  optimised saturated hydraulic conductivity}
#' \item{S:} { optimised sorptivity [LT^-0.5]}
#' \item{ratio:} { optimised ratio parameter}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Valiantzas, J. D. (2010). New linearized two-parameter infiltration equation 
#' for direct determination of conductivity and sorptivity. 
#' Journal of Hydrology, 387. doi: 10.1016/j.jhydrol.2009.12.049
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' Valiantzas1<-valiantzas(data=data,time="time",I="I")
#' plot(Valiantzas1)
#' predict(Valiantzas1)
#' coef(Valiantzas1)
#' 
#' 
#' #infiltration rate
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' Valiantzasr<-valiantzas(data=assin_breko,time="minutes",I="cm.hr")
#' 
#' #group simulation
#' Valiantzasg<-valiantzas(data=data,time="minutes",I="cm.hr",group="ID")
#' coef(Valiantzasg)
#' #generic function
valiantzas<-function(data=NULL,time,I,S=0.1,Ks=0.3,type="linear",ratio=0.5,group=NULL) UseMethod ("valiantzas")
#' @rdname valiantzas
#' @export
valiantzas.default <- function(data=NULL,time,I,S=0.1,Ks=0.3,type="linear",ratio=0.5,group=NULL)
{
#stop warning from displaying
options(warn=-1)

if(data==""||is.null(data))
{
#predict=0.5*S*(time^-0.5)+Ks
#predict=((Ks*time2)+(S*sqrt(time2)*(1-((Ks*time2)/single[[I]]))^(1-ratio)))

#return(print (list(ft=predict)))
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
ones <- c(S = S, Ks = Ks,ratio=ratio) # all ones start
#determine whether it is rate or cumulative data
rate="yes"
if(data[[I]][length(data[[I]])]>data[[I]][1])
{
#cumulative function ###############################
valiantzasF<-paste(I,"~((Ks*time)+(S*sqrt(time)*(1-((Ks*time)/I))^(1-ratio)))")
valiantzasF<-paste(I,"~((Ks*",time,")+(S*sqrt(",time,")*(1-((Ks*",time,")/",I,"))^(1-ratio)))")


rate="no"
}
else
{
#rate function #######################################
valiantzasF<-paste(I,"~(0.5*S","*",time,"^(-0.5))+ Ks")

}
#check for grouped data and execute the optimisation function
if(is.null(group)){
if(rate=="no"&&type=="linear"){
dependent=(data[[I]]^2)/(data[[time]])
independent=(data[[I]])
linear=lm(dependent~ independent)
S=sqrt(coef(linear)[[1]])
Ks=coef(linear)[[2]]
valiantzas=c(S,Ks)
}else{

valiantzas<- nlxb(valiantzasF, start = ones, trace = FALSE, data = data)
}
#print(valiantzas)

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
if(rate=="no"&&(type=="linear"||type=="lr")){
dependent=(single[[I]]^2)/(single[[time]])
independent=(single[[I]])
linear=lm(dependent~ independent)
S=sqrt(coef(linear)[[1]])
Ks=coef(linear)[[2]]
valiantzas=c(S,Ks)
}else{
valiantzas<- nlxb(valiantzasF, start = ones, trace = FALSE, data = single)
S=coef(valiantzas)[1]
Ks=coef(valiantzas)[2]
ratio=coef(valiantzas)[2]

}
print(valiantzas)
print("....................................................................................")
#valiantzas paramters #################################
groupdata=single[[group]]
time2=single[[time]]
#prediction cumulative equation ##############################
#predict2=S*(time2^0.5)+Ks*time2
predict2=((Ks*time2)+(S*sqrt(time2)*(1-((Ks*time2)/single[[I]]))^(1-ratio)))

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
if(rate=="yes"||type=="nonlinear"){
print(valiantzas)
S=coef(valiantzas)[1]
Ks=coef(valiantzas)[2]##########################################
ratio=coef(valiantzas)[3]
}


time2=data[[time]]
I=data[[I]]
time=data[[time]]
###########################################################
######################
if(type=="nonlinear"){
predict=((Ks*time)+(S*sqrt(time)*(1-((Ks*time)/I))^(1-ratio)))
}

if(type=="linear"){
predict=((0.5*Ks*time)+(S*sqrt(time)*(1+((0.5*Ks/S)^2)*time)^0.5))
}


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
factor <- list(valiantzas=valiantzas,data=data,time=time,I=I,S=S,Ks=Ks,group=group,
predict=predict,rate=rate,formular=valiantzasF,addoutput=addoutput,output=addoutput,type=type)

factor$call<-match.call()

class(factor)<-"valiantzas"
factor
}
#' @rdname valiantzas
#' @export
#predict function#########################
predict.valiantzas<-function(object,time=NULL,...)
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
if(object$type=="nonlinear"){
predict=((object$Ks*time)+(object$S*sqrt(time)*(1-((object$Ks*time)/object$I))^(1-object$ratio)))
}

if(object$type=="linear"){
predict=((0.5*object$Ks*time)+(object$S*sqrt(time)*(1+((0.5*object$Ks/object$S)^2)*time)^0.5))
}
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
predict=object$addoutput$S*(time^0.5)+object$addoutput$Ks*time
if(object$type=="nonlinear"){
predict=((object$addoutput$Ks*time)+(object$addoutput$S*sqrt(time)*(1-((object$addoutput$Ks*time)/object$I))^(1-object$addoutput$ratio)))
}

if(object$type=="linear"){
predict=((0.5*object$addoutput$Ks*time)+(object$addoutput$S*sqrt(time)*(1+((0.5*object$addoutput$Ks/object$addoutput$S)^2)*time)^0.5))
}
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}

print(predict)
}
}
#' @rdname valiantzas
#' @export
#plot function
plot.valiantzas<-function(x,xlab="Time(Minutes)",ylab="Cumulative (mm)",main=NULL,layout=NULL,...)
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
#' @rdname valiantzas
#' @export
#summary function
summary.valiantzas<-function(object,...)
{


x<-object$valiantzas

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
#' @rdname valiantzas
#' @export
#print function
print.valiantzas<-function(x,...)
{
  object <- x
if(is.null(object$group)){
x<-object$valiantzas
print((x))
}
else
{
coef=aggregate(cbind(Ks, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @rdname valiantzas
#' @export
#coef function
coef.valiantzas<-function(object,...)
{
x<-object$valiantzas
if(is.null(object$group)){
if(object$type=="nonlinear"||object$type=="nlr"||object$type=="nls"){
coef=(coef(x))
}
if(object$type=="linear"||object$type=="lr"){
coef=list(S=object$S,Ks=object$Ks)
}
}
else
{
coef=aggregate(cbind(Ks, S) ~ groupid, data = object$addoutput, mean)
}
print((coef=coef))
}

#statistics of error

