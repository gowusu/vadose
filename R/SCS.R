#input function ############################################
#input function ############################################
#' @title SCS infiltration parameter optmisation in R
#' @description This function optimises Jury et al (1991)
#' cumulative infiltration (I) parameters
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @inheritParams brutsaert
#' @param a constant
#' @param b constant
#' @return 
#' \itemize{
#' \item{a:} {  optimised constant}
#' \item{b:} { optimised constant}
#' \item{output:} { output of the group simulation}
#' }
#' @author George Owusu
#' @references 
#' Jury, W. A., Gardner, W. R., & H., G. W. (Eds.). (1991). 
#' Soil physics (5th ed.). New York (NY): John Wiley & Sons.  
#' @export
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' SCS1<-scs(data=data,time="time",I="I")
#' #print(gof(SCS1))
#' plot(SCS1)
#' predict(SCS1)
#' coef(SCS1)
#' #gof(SCS1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' assin_breko<- subset(data, ID=="41A20_1")
#' SCSg<-scs(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(SCSg)
#' @rdname SCS
#generic function
scs<-function(data,time,I,a=0.1,b=0.1,group=NULL) UseMethod ("scs")############
#' @export
#' @rdname SCS
#default function
scs.default<-function(data,time,I,a=0.1,b=0.1,group=NULL)###################
{
  options(warn=-1)
  #decalare the group data incase group variable is not null
  addoutput=NULL
  #set parameters of optimsation functions######################
  ones <- c(a=a, b= b) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    scsF<-paste(I,"~a*(",time,"^b)+0.6985")
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    scsF<-paste(I,"~(a*b*(",time,"^(b-1)))")
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    scs<- nlxb(scsF, start = ones, trace = FALSE, data = data)
    print(scs)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame###############################################
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),a=numeric(),b=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      scs<- nlxb(scsF, start = ones, trace = FALSE, data = single)
      print(scs)
      print("....................................................................................")
      #scs paramters #################################
      a=coef(scs)[1]
      b=coef(scs)[2]
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction cumulative equation ##############################
      predict2=a*(time2^b)+0.6985
      
      if(rate=="yes"){
        #predict rate###############################################
        predict2=(a*b*(time2^(b-1)))
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,a=a,b=b))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    a=coef(scs)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    b=coef(scs)[2]
    ##########################################
    ######################
    predict=a*(time2^b)+0.6985
    
    if(rate=="yes"){
      #####################################################
      predict=(a*b*(time2^(b-1)))
    }
  }
  else
  {
    ################################################
    a=addoutput$a
    I=data[[I]]
    time=data[[time]]
    ##############################################
    b=addoutput$b
    predict=addoutput$predict
    time=addoutput$time
  }
  #return varibales ########################################
  factor<-list(scs=scs,data=data,time=time,I=I,a=a,b=b,group=group,
       predict=predict,rate=rate,formular=scsF,addoutput=addoutput,output=addoutput)


factor$call<-match.call()

class(factor)<-"scs"
factor
}
#' @export
#' @rdname SCS
#predict function#########################
predict.scs<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")
#rate
if(!is.null(time)&object$rate=="yes"){####################
predict=(1-exp(-object$As*sqrt(time)))*(object$S/object$As)+(object$K*time)

predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$a*(time^object$b)+0.6985

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=object$a*object$k*(time^(object$a-1))###########
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$a*(time^object$b)+0.6985################
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
#' @rdname SCS
#plot function
plot.scs<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
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

#summary function
#' @export
#' @rdname SCS
summary.scs<-function(object,...)
{


x<-object$scs

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{########################
coef=aggregate(cbind(a,b) ~ groupid, data = object$addoutput, mean)##############
print(coef)
}

}

#print function
print.scs<-function(object)
{

if(is.null(object$group)){
x<-object$scs
print((x))
}
else
{######################################################
coef=aggregate(cbind(a,b) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}

#coef function
coef.scs<-function(object,...)
{
x<-object$scs
if(is.null(object$group)){
coef=(coef(x))

}
else################################################
{
coef=aggregate(cbind(a,b) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

