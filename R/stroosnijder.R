#input function ############################################
#' @title stroosnijder infiltration parameter optimisation in R
#' @description This function optimises stroosnijder (1976)
#' cumulative infiltration  (I): Ks and S. 
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @inheritParams brutsaert
#' @return 
#' \itemize{
#' \item{Ks:} {  optimised saturated hydraulic conductivity [L/T]}
#' \item{S:} { optimised sorptivity [LT^-0.5]}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' stroosnijder, L., 1976. Infiltratie en Herverdeling van Water in Gronde. 
#' Versl.Landbouwkd, Onderz. 847 
#' @export
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' stroosnijder1<-stroosnijder(data=data,time="time",I="I")
#' #print(gof(stroosnijder1))
#' plot(stroosnijder1)
#' predict(stroosnijder1)
#' coef(stroosnijder1)
#' #gof(stroosnijder1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' stroosnijderg<-stroosnijder(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(stroosnijderg)
#' @rdname stroosnijder
#generic function
stroosnijder<-function(data=NULL,time,I,S=0.1,K=0.1,para=2,group=NULL) UseMethod ("stroosnijder")
#' @export
#' @rdname stroosnijder
#default function
stroosnijder.default<-function(data=NULL,time,I,S=0.1,K=0.1,para=2,group=NULL)
{
  options(warn=-1)
  
  if(data==""||is.null(data))
  {
    predict=((K*time)+((3*S^2)/(4*K))*(1-(exp(-(4*K*(time^0.5))/(3*S)))))
    
    return(print (list(ft=predict)))
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
  ones <- c(S = S, K = K) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    #stroosnijderF<-paste(I,"~(K*",time,")+((3*S^2)/(4*K))*(1-exp(-(4*K*(",time,"^0.5))/(3*S)))")
    stroosnijderF<-paste(I,"~((K*",time,")+((3*S^2)/(4*K))*(1-(exp(-(4*K*(",time,"^0.5))/(3*S)))))")
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    
    stroosnijderF<-paste(I,"~(0.5*S","*",time,"^(-0.5))+ K")
    
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    stroosnijder<- nlxb(stroosnijderF, start = ones, trace = FALSE, data = data)
    print(stroosnijder)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),K=numeric(),S=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      stroosnijder<- nlxb(stroosnijderF, start = ones, trace = FALSE, data = single)
      print(stroosnijder)
      print("....................................................................................")
      #stroosnijder paramters #################################
      S=coef(stroosnijder)[1]
      K=coef(stroosnijder)[2]
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction cumulative equation ##############################
      predict2=((K*time2)+((3*S^2)/(4*K))*(1-(exp(-(4*K*(time2^0.5))/(3*S)))))
      
      if(rate=="yes"){
        #predict rate###############################################
        predict2=0.5*S*(time2^-0.5)+K
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,K=K,S=S))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    S=coef(stroosnijder)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    K=coef(stroosnijder)[2]##########################################
    ######################
    predict=((K*time2)+((3*S^2)/(4*K))*(1-(exp(-(4*K*(time2^0.5))/(3*S)))))
    
    
    if(rate=="yes"){
      #####################################################
      predict=0.5*S*(time2^-0.5)+K
    }
  }
  else
  {
    ################################################
    S=addoutput$S
    I=data[[I]]
    time=data[[time]]
    ##############################################
    K=addoutput$K
    predict=addoutput$predict
    time=addoutput$time
    
  }
  #return varibales ########################################
  factor<-list(stroosnijder=stroosnijder,data=data,time=time,I=I,S=S,K=K,group=group,
       predict=predict,rate=rate,formular=stroosnijderF,addoutput=addoutput,output=addoutput)

factor$call<-match.call()

class(factor)<-"stroosnijder"
factor
}
#' @export
#' @rdname stroosnijder
#predict function#########################
predict.stroosnijder<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$K
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$K*time)+((3*object$S^2)/(4*object$K))*(1-(exp(-(4*object$K*(time^0.5))/(3*object$S)))))

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=0.5*object$S*(time^-0.5)+object$K

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$addoutput$K*time)+((3*object$addoutput$S^2)/(4*object$addoutput$K))*(1-(exp(-(4*object$addoutput$K*(time^0.5))/(3*object$addoutput$S)))))

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
#' @rdname stroosnijder
#plot function
plot.stroosnijder<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
{
object=x
x<-object
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
#' @export
#' @rdname stroosnijder
#summary function
summary.stroosnijder<-function(object,...)
{


x<-object$stroosnijder

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{
coef=aggregate(cbind(K, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}

}
#' @export
#' @rdname stroosnijder
#print function
print.stroosnijder<-function(x,...)
{
object=x
if(is.null(object$group)){
x<-object$stroosnijder
print((x))
}
else
{
coef=aggregate(cbind(K, S) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @export
#' @rdname stroosnijder
#coef function
coef.stroosnijder<-function(object,...)
{
x<-object$stroosnijder
if(is.null(object$group)){
coef=(coef(x))

}
else
{
coef=aggregate(cbind(K, S) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

#statistics of error

