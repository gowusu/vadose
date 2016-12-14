#input function ############################################
#' @title R optimisation of philip infiltration parameters
#' @description This function can linearly or nonlinearly optimise philip
#' cumulative infiltration (I) or infiltration rate (i) parameters: A and S. 
#'
#' @inheritParams BEST
#' @inheritParams lass3
#' @param I cumulative infiltration (I) or infiltration rate (i) [mm]
#' @param S numeric. sorptivity parameter
#' @param A numeric.It is related Hydraulic Conductivity parameter
#' @param type character. It takes "linear" or "nonlinear" Philip equation
#' @inheritParams OFEST
#' @return 
#' \itemize{
#' \item{A:} { opitmised  A parameter}
#' \item{Ks:} {  opitmised  saturated hydraulic conductivity=0.5*A [L/T]}
#' \item{S:} { optimised sorptivity [LT^-0.5]}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @export
#' @references 
#' Philip, J. R. (1957). The theory of infiltration:Sorptivity and algebraic 
#' infiltration equations. Soil Science, 84, 257-264. 
#'
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' philip1<-philip(data=data,time="time",I="I")
#' print(gof(philip1))
#' plot(philip1)
#' predict(philip1)
#' coef(philip1)
#' gof(philip1)
#' 
#' #infiltration rate
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' philipr<-philip(data=assin_breko,time="minutes",I="cm.hr")
#' 
#' #group simulation
#' philipg<-philip(data=data,time="minutes",I="cm.hr",group="ID")
#' coef(philipg)
#' modprediction<-philip(time=1,S=0.1,A=0.5)
philip<-function(data=NULL,time,I,S=0.1,A=0.1,type="nonlinear",group=NULL) UseMethod ("philip")
#' @export
#' @rdname philip
philip.default <- function(data=NULL,time,I,S=0.1,A=0.1,type="nonlinear",group=NULL)
{
  
#stop warning from displaying
options(warn=-1)

if(data==""||is.null(data))
{
predict=0.5*S*(time^-0.5)+A
return(predict)
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
ones <- c(S = S, A = A) # all ones start
#determine whether it is rate or cumulative data
rate="yes"
if(data[[I]][length(data[[I]])]>data[[I]][1])
{
#cumulative function ###############################
philipF<-paste(I,"~(S","*",time,"^0.5)+(A","*",time,")")

rate="no"
}
else
{
#rate function #######################################
philipF<-paste(I,"~(0.5*S","*",time,"^(-0.5))+ A")

}
#check for grouped data and execute the optimisation function
if(is.null(group)){
if(type=="linear"){
dependent=data[[I]]/sqrt(data[[time]])
independent=sqrt(data[[time]])
if(rate=="yes"){
dependent=data[[I]]
independent=data[[time]]^-0.5
}
linear=lm(dependent~ independent)
S=coef(linear)[[1]]
A=coef(linear)[[2]]
if(rate=="yes"){
S=coef(linear)[[1]]*2
}
philip=c(S,A)
}else{

philip<- nlxb(philipF, start = ones, trace = FALSE, data = data)
}
#print(philip)

}
else
{
# write the group function
aggdata =row.names(table(data[group]))
#create group data frame
addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),A=numeric(),S=numeric())
i=1
while(i<=length(aggdata)){
print(paste("Group Number:",aggdata[i]))
single=data[data[group]==aggdata[i],]
#group function
if((type=="linear"||type=="lr")){
dependent=single[[I]]/sqrt(single[[time]])
independent=sqrt(single[[time]])
if(rate=="yes"){
dependent=single[[I]]
independent=single[[time]]^-0.5
}
linear=lm(dependent~ independent)
S=coef(linear)[[1]]
if(rate=="yes"){
S=coef(linear)[[1]]*2
}

A=coef(linear)[[2]]
philip=c(S,A)
}else{
philip<- nlxb(philipF, start = ones, trace = FALSE, data = single)
S=coef(philip)[1]
A=coef(philip)[2]
}
#print(philip)
#print("....................................................................................")
#philip paramters #################################
groupdata=single[[group]]
time2=single[[time]]
#prediction cumulative equation ##############################
predict2=S*(time2^0.5)+A*time2
if(rate=="yes"){
#predict rate###############################################
predict2=0.5*S*(time2^-0.5)+A
}
addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,A=A,S=S))

i=i+1
}
}

#equations for ungrouped data ############################
if(is.null(group)){
#prediction################################################
if(type=="nonlinear"){
S=coef(philip)[1]
A=coef(philip)[2]##########################################
}
time2=data[[time]]
I=data[[I]]
time=data[[time]]
###########################################################
######################
predict=S*(time2^0.5)+A*time2

if(rate=="yes"){
#####################################################
predict=0.5*S*(time2^-0.5)+A
}
}
else
{
################################################
S=addoutput$S
I=data[[I]]
time=data[[time]]
##############################################
A=addoutput$A
predict=addoutput$predict
time=addoutput$time

}
#return varibales ########################################
factor <- list(philip=philip,data=data,time=time,I=I,S=S,A=A,Ks=0.5*A,group=group,
predict=predict,rate=rate,formular=philipF,addoutput=addoutput,output=addoutput,type=type)
factor$call<-match.call()

class(factor)<-"philip"
factor
}

#generic function
#predict function#########################
#' @rdname philip
#' @export
predict.philip<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$A
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$S*(time^0.5)+object$A*time
predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$A

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$addoutput$S*(time^0.5)+object$addoutput$A*time
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
#' @rdname philip
#plot function
plot.philip<-function(x,xlab="Time(Minutes)",ylab="Cumulative (mm)",main=NULL,layout=NULL,...)
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
#' @rdname philip
#' @export
#summary function
summary.philip<-function(object,...)
{


x<-object$philip

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{
coef=aggregate(cbind(A, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}

}
#' @rdname philip
#' @export
#print function
print.philip<-function(x,...)
{
object <- x
if(is.null(object$group)){
x<-object$philip
print((x))
}
else
{
coef=aggregate(cbind(A, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @rdname philip
#' @export
#coef function
coef.philip<-function(object,...)
{
x<-object$philip
if(is.null(object$group)){
if(object$type=="nonlinear"||object$type=="nlr"||object$type=="nls"){
coef=(coef(x))
}
if(object$type=="linear"||object$type=="lr"){
coef=list(S=object$S,A=object$A)
}
}
else
{
coef=aggregate(cbind(A, S) ~ groupid, data = object$addoutput, mean)
}
print(list(coef=coef))
}

#statistics of error

