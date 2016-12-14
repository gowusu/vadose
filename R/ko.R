#input function ############################################
#' @title Kostiakov infiltration parameter optmisation in R
#' @description This function optimises Kostiakov (1932)
#' cumulative infiltration (I) parameters.
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @inheritParams brutsaert
#' @param a constant
#' @param k constant
#' @return 
#' \itemize{
#' \item{a:} {  optimised constant}
#' \item{k:} { optimised constant}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Kostiakov, A. N. (1932). On the dynamics of the coefficient of water-percolation in soils and on the 
#' necessity for studying it from a dynamic point of view for purposes of amelioration. 6.  
#' @export
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' kostiakov1<-kostiakov(data=data,time="time",I="I")
#' #print(gof(kostiakov1))
#' plot(kostiakov1)
#' predict(kostiakov1)
#' coef(kostiakov1)
#' #gof(kostiakov1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' kostiakovg<-kostiakov(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(kostiakovg)
#' @rdname kostiakov
#generic function

kostiakov<-function(data=NULL,time,I,a=0.1,k=0.1,group=NULL) UseMethod ("kostiakov")############
#' @export
#' @rdname kostiakov
#default function
kostiakov.default<-function(data=NULL,time,I,a=0.1,k=0.1,group=NULL)###################
{
  options(warn=-1)
  
  if(data==""||is.null(data))
  {
    predict=a*k*(time^(a-1))
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
  ones <- c(a=a, k = k) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    koF<-paste(I,"~k*",time,"^a")
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    koF<-paste(I,"~a*k*",time,"^(a-1)")
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    kostiakov<- nlxb(koF, start = ones, trace = FALSE, data = data)
    print(kostiakov)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame###############################################
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),k=numeric(),a=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      kostiakov<- nlxb(koF, start = ones, trace = FALSE, data = single)
      print(kostiakov)
      print("....................................................................................")
      #kostiakov paramters #################################
      a=coef(kostiakov)[1]
      k=coef(kostiakov)[2]
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction cumulative equation ##############################
      predict2=k*(time2^a)
      
      if(rate=="yes"){
        #predict rate###############################################
        predict2=a*k*(time2^(a-1))
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,k=k,a=a))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    a=coef(kostiakov)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    k=coef(kostiakov)[2]##########################################
    ######################
    predict=k*(time2^a)
    
    if(rate=="yes"){
      #####################################################
      predict=a*k*(time2^(a-1))
    }
  }
  else
  {
    ################################################
    a=addoutput$a
    I=data[[I]]
    time=data[[time]]
    ##############################################
    k=addoutput$k
    predict=addoutput$predict
    time=addoutput$time
  }
  #return varibales ########################################
  factor<-list(kostiakov=kostiakov,data=data,time=time,I=I,a=a,k=k,group=group,
       predict=predict,rate=rate,formular=koF,addoutput=addoutput,output=addoutput)


factor$call<-match.call()

class(factor)<-"kostiakov"
factor
}
#' @export
#' @rdname kostiakov
#predict function#########################
predict.kostiakov<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){####################
predict=object$a*object$k*(time^(object$a-1))
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$k*(time^object$a)
predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=object$a*object$k*(time^(object$a-1))
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$k*(time^object$a)
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
#' @rdname kostiakov
#plot function
plot.kostiakov<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
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
#' @rdname kostiakov
#summary function
summary.kostiakov<-function(object,...)
{


x<-object$kostiakov

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{########################
coef=aggregate(cbind(k, a) ~ groupid, data = object$addoutput, mean)##############
print(coef)
}

}

#print function
print.kostiakov<-function(x,...)
{
object=x
if(is.null(object$group)){
x<-object$kostiakov
print((x))
}
else
{######################################################
coef=aggregate(cbind(k, a) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}

#coef function
coef.kostiakov<-function(object,...)
{
x<-object$kostiakov
if(is.null(object$group)){
coef=(coef(x))

}
else################################################
{
coef=aggregate(cbind(k, a) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

