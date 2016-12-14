#' Horton infiltration parameter optmisation in R
#' 
#' @description This function optimises Horton (1940)
#'  infiltration parameters: fo,fc, and k. It also predicts infiltration.
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param fc final constant rate of infiltration at saturation  
#' @param fo parameter for initial infiltration capacity
#' @param k a constant depending primarily upon soil and vegetation
#' @author George Owusu
#' 
#' @references 
#' Horton, R. E. (1940). An approach toward a physical 
#' interpretation of infiltration capacity Proc. Sci. Soc. Amer., 5, 399-417. 
#'
#' @return
#' \itemize{
#' \item{fc:} {  final constant rate of infiltration at saturation}
#' \item{fo:} { parameter for initial infiltration capacity}
#' \item{k:} { a constant depending primarily upon soil and vegetation}
#' \item{predict:}{predicted infiltration}
#' \item{output:} { output of the group simulation}
#' }
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' horton1<-horton(data=data,time="time",I="I")
#' #print(gof(horton1))
#' plot(horton1)
#' predict(horton1)
#' coef(horton1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' horton1g<-horton(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(horton1g)
#' #generic function ######################################
#' @rdname horton
horton<-function(data=NULL,time,I,fc=0.1,fo=0.1,k=0.1,group=NULL) UseMethod ("horton")
#' @export
#' @rdname horton
#default function######################################

horton.default<-function(data=NULL,time,I,fc=0.1,fo=0.1,k=0.1,group=NULL)
{
  
  #stop warning from displaying
  options(warn=-1)
  
  if(data==""||is.null(data))
  {
    predict=fc+(fo-fc)*exp(-k*time)
    #print (predict)
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
  ones <- c(fc=fc, fo=fo, k=k) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    hortonF<-paste(I,"~fc*",time,"+((fo-fc)/k)+(1-exp(-k*",time,"))")
    
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    hortonF<-paste(I,"~fc+(fo-fc)*exp(-k*",time,")")
    
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    horton<- nlxb(hortonF, start = ones, trace = FALSE, data = data)
    print(horton)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),fc=numeric(),fo=numeric(),k=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      horton<- nlxb(hortonF, start = ones, trace = FALSE, data = single)
      print(horton)
      print("....................................................................................")
      #horton paramters #################################
      fc=coef(horton)[1]
      fo=coef(horton)[2]
      k=coef(horton)[3]
      
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction of cumulative equation ##############################
      predict2=fc*time2+((fo-fc)/k)+(1-exp(-k*time2))
      if(rate=="yes"){
        #predict rate###############################################
        predict2=fc+(fo-fc)*exp(-k*time2)
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,fc=fc,fo=fo,k=k))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    fc=coef(horton)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    fo=coef(horton)[2]##########################################
    k=coef(horton)[3]
    ###################### cumulative
    predict=fc*time2+((fo-fc)/k)+(1-exp(-k*time2))
    if(rate=="yes"){
      #####################################################
      predict=fc+(fo-fc)*exp(-k*time2)
      
    }
  }
  else
  {
    ################################################
    fc=addoutput$fc
    I=data[[I]]
    time=data[[time]]
    ##############################################
    fo=addoutput$fo
    k=addoutput$k
    predict=addoutput$predict
    time=addoutput$time
    
  }
  #return varibales ########################################
  factor<-list(horton=horton,data=data,time=time,I=I,fc=fc,fo=fo,k=k,group=group,
       predict=predict,rate=rate,formular=hortonF,addoutput=addoutput,output=addoutput)


factor$call<-match.call()

class(factor)<-"horton"
factor
}

#' @rdname horton
#' @export
#predict function
predict.horton<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){###########################
predict=object$fc+(object$fo-object$fc)*exp(-object$k*time)

predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$fc*time+((object$fo-object$fc)/object$k)+(1-exp(-object$k*time))

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=object$fc+(object$fo-object$fc)*exp(-object$k*time)

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$addoutput$fc*time+((object$addoutput$fo-object$addoutput$fc)/object$addoutput$k)+(1-exp(-object$addoutput$k*time))

predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))


}

print(predict)
}
}
#' @rdname horton
#' @export
#plot function
plot.horton<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
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
#matrix plot
i=1
while(i<=length(aggdata)){
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
#' @rdname horton
#summary function
summary.horton<-function(object)
{


x<-object$horton

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else######################################
{
coef=aggregate(cbind(fc, fo,k) ~ groupid, data = object$addoutput, mean)
print(coef)
}

}
#' @rdname horton
#print function
print.horton<-function(x,...)
{
  object=x
if(is.null(object$group)){
x<-object$horton
print((x))
}
else######################################
{
coef=aggregate(cbind(fc, fo,k) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}

#coef function
#' @export 
#' @rdname horton
coef.horton<-function(object,...)
{
x<-object$horton
if(is.null(object$group)){
coef=(coef(x))

}
else######################################
{
coef=aggregate(cbind(fc, fo,k) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}



