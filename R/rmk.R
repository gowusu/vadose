#input function ############################################
#' @title Parhi Revised modified kostiakov infiltration parameter optimisation in R
#' @description This function optimises Parhi et al (2007) revised modified kostiakov
#' cumulative infiltration (I) 
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @inheritParams brutsaert
#' @param a1 constant
#' @param a2 constant
#' @param b1 constant
#' @param b2 constant
#' @return 
#' \itemize{
#' \item{a1:} {  optimised constant}
#' \item{a2:} { optimised constant}
#' \item{b1}{optimised constant}
#' \item{b2}{optimised constant}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Parhi, P. K., Mishra, S. K., & Singh, R. (2007). A modification to Kostiakov and 
#' modified Kostiakov  infiltration modeles. Water Resour. Manage., 21, 1973-1989.   
#' @export
#' @examples
#' 
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' revised.modified.kostiakov1<-revised.modified.kostiakov(data=data,time="time",I="I")
#' #print(gof(revised.modified.kostiakov1))
#' plot(revised.modified.kostiakov1)
#' predict(revised.modified.kostiakov1)
#' coef(revised.modified.kostiakov1)
#' #gof(revised.modified.kostiakov1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' revised.modified.kostiakov1vg<-revised.modified.kostiakov(data=data,time="minutes",
#' I="CumInfil",group="ID")
#' coef(revised.modified.kostiakov1vg)
#' 
#' @rdname revised.modified.kostiakov
#' 
revised.modified.kostiakov<-function(data=NULL,time,I,a1=0.1,a2=0.1,b1=0.1,b2=0.1,group=NULL) UseMethod ("revised.modified.kostiakov")############
#' @export
#' @rdname revised.modified.kostiakov
#default function
revised.modified.kostiakov.default<-function(data=NULL,time,I,a1=0.1,a2=0.1,b1=0.1,b2=0.1,group=NULL)###################
{
  
  options(warn=-1)
  if(data==""||is.null(data))
  {
    predict=(a1*(time^b1)+a2*(time^(-b2)))
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
  ones <- c(a1=a1,a2=a2,b1=b1,b2=b2) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    #rmkF<-paste(I,"~((a1/(b+1))*(",time,"^(b1+1))+(a2/(1-b2))*",time,"^(1-b2))")
    rmkF<-paste(I,"~((a1/(b1+1))*(",time,"^(b1+1))+(a2/(1-b2))*(",time,"^(1-b2)))")
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    rmkF<-paste(I,"~(a1*(",time,"^b1)+a2*(",time,"^(-b2)))")
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    revised.modified.kostiakov<- nlxb(rmkF, start = ones, trace = FALSE, data = data)
    print(revised.modified.kostiakov)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame###############################################
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),a1=numeric(),a2=numeric(),b1=numeric(),b2=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      revised.modified.kostiakov<- nlxb(rmkF, start = ones, trace = FALSE, data = single)
      print(revised.modified.kostiakov)
      print("....................................................................................")
      #revised.modified.kostiakov paramters #################################
      a1=coef(revised.modified.kostiakov)[1]
      a2=coef(revised.modified.kostiakov)[2]
      b1=coef(revised.modified.kostiakov)[3]
      b2=coef(revised.modified.kostiakov)[4]
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction cumulative equation ##############################
      predict2=((a1/(b1+1))*(time2^(b1+1))+(a2/(1-b2))*(time2^(1-b2)))
      
      if(rate=="yes"){
        #predict rate###############################################
        predict2=(a1*(time2^b1)+a2*(time2^(-b2)))
        
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,a1=a1,a2=a2,b1=b1,b2=b2))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    a1=coef(revised.modified.kostiakov)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    a2=coef(revised.modified.kostiakov)[2]
    b1=coef(revised.modified.kostiakov)[3]
    b2=coef(revised.modified.kostiakov)[4]
    
    ##########################################
    ######################
    predict=((a1/(b1+1))*(time2^(b1+1))+(a2/(1-b2))*(time2^(1-b2)))
    
    if(rate=="yes"){
      #####################################################
      predict=(a1*(time2^b1)+a2*(time2^(-b2)))
    }
  }
  else
  {
    ################################################
    a1=addoutput$a1
    I=data[[I]]
    time=data[[time]]
    ##############################################
    a2=addoutput$a2
    b1=addoutput$b1
    b2=addoutput$b2
    
    predict=addoutput$predict
    time=addoutput$time
  }
  #return varibales ########################################
  factor=list(revised.modified.kostiakov=revised.modified.kostiakov,data=data,time=time,I=I,a1=a1,a2=a2,b1=b1,b2=b2,group=group,
       predict=predict,rate=rate,formular=rmkF,addoutput=addoutput,output=addoutput)


factor$call<-match.call()

class(factor)<-"revised.modified.kostiakov"
factor
}
#' @export
#' @rdname revised.modified.kostiakov
#predict function#########################
predict.revised.modified.kostiakov<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")
#rate
if(!is.null(time)&object$rate=="yes"){####################
predict=(object$a1*(time^object$b1)+object$a2*(time^(-object$b2)))

predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$a1/(object$b1+1))*(time^(object$b1+1))+(object$a2/(1-object$b2))*(time^(1-object$b2)))

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=(object$a1*(time^object$b1)+object$a2*(time^(-object$b2)))

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=((object$a1/(object$b1+1))*(time^(object$b1+1))+(object$a2/(1-object$b2))*(time^(1-object$b2)))
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
#' @rdname revised.modified.kostiakov
#plot function
plot.revised.modified.kostiakov<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
{
object=x
#x<-object
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
#' @rdname revised.modified.kostiakov
#summary function
summary.revised.modified.kostiakov<-function(object,...)
{


x<-object$revised.modified.kostiakov

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{########################
coef=aggregate(cbind(a1,a2,b1,b2) ~ groupid, data = object$addoutput, mean)##############
print(coef)
}

}

#print function
print.revised.modified.kostiakov<-function(object)
{

if(is.null(object$group)){
x<-object$revised.modified.kostiakov
print((x))
}
else
{######################################################
coef=aggregate(cbind(a1,a2,b1,b2) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @export
#' @rdname revised.modified.kostiakov
#coef function
coef.revised.modified.kostiakov<-function(object,...)
{
x<-object$revised.modified.kostiakov
if(is.null(object$group)){
coef=(coef(x))

}
else################################################
{
coef=aggregate(cbind(a1,a2,b1,b2) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}




