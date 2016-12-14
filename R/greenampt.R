#' Green and Ampt infiltration parameter optmisation in R
#' 
#' @description This function optimises Green and Ampt (1911)
#'  infiltration parameters: Ks and G. It also predicts infiltration.
#'
#' @inheritParams philip
#' @inheritParams OFEST
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param Ks Hydraulic conductivity  
#' @param G Green and Ampt parameter that is equivalent to Sorptivity
#' @author George Owusu
#' 
#' @references 
#' Green, W. A., & Ampt, G. A., 1911..4,1-24. (1911). Studies on soil  physics:1.
#' The flow of air and water through soils. Journal of Agricultural Science, 4(1-24).  
#' @return
#' \itemize{
#' \item{Ks:} {  Hydraulic conductivity [L]}
#' \item{G:} { Green and Ampt parameter that is equivalent to Sorptivity [LT^0.5]}
#' \item{predict:}{predicted infiltration}
#' \item{output:} { output of the group simulation}
#' }
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' greenampt1<-greenampt(data=data,time="time",I="I")
#' #print(gof(greenampt1))
#' #plot(greenampt1)
#' predict(greenampt1)
#' coef(greenampt1)
#' 
#' #group simulation
#' data=read.csv(system.file("ext","sys","infiltration2.csv",package="vadose"))
#' greenampt1g<-greenampt(data=data,time="minutes",I="CumInfil",group="ID")
#' coef(greenampt1g)
#' #generic function ######################################
#' @rdname greenampt
#generic function
greenampt<-function(data,time,I,Ks=0.1,G=0.1,group=NULL) UseMethod ("greenampt")############
#' @export
#' @rdname greenampt
#default function
greenampt.default<-function(data,time,I,Ks=0.1,G=0.1,group=NULL)###################
{
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
  ones <- c(Ks=Ks, G = G) # all ones start
  #determine whether it is rate or cumulative data
  rate="yes"
  f=NULL
  if(data[[I]][length(data[[I]])]<data[[I]][1])
  {
    return (print("the function accepts only cumulative infiltration in cm"))
    data$f=data[I]
    #data[I]=cumsum(data[I])
    f="yes"
  }
  if(data[[I]][length(data[[I]])]>data[[I]][1])
  {
    #cumulative function ###############################
    greenamptF<-paste(I,"~Ks*",time,"+G*log(1+(",I,"/G))")
    
    rate="no"
  }
  else
  {
    #rate function #######################################
    greenamptF<-paste(I,"~G*Ks*",time,"^(G-1)")
  }
  #check for grouped data and execute the optimisation function
  if(is.null(group)){
    greenampt<- nlxb(greenamptF, start = ones, trace = FALSE, data = data)
    print(greenampt)
    
  }
  else
  {
    # write the group function
    aggdata =row.names(table(data[group]))
    #create group data frame###############################################
    addoutput=data.frame(groupid=factor(),time=numeric(),observed=numeric(),predict=numeric(),Ks=numeric(),G=numeric())
    i=1
    while(i<=length(aggdata)){
      print(paste("Group Number:",aggdata[i]))
      single=data[data[group]==aggdata[i],]
      #group function
      greenampt<- nlxb(greenamptF, start = ones, trace = FALSE, data = single)
      print(greenampt)
      print("....................................................................................")
      #greenampt paramters #################################
      Ks=coef(greenampt)[1]
      G=coef(greenampt)[2]
      groupdata=single[[group]]
      time2=single[[time]]
      #prediction cumulative equation ##############################
      predict2=Ks*(time2)+(G*log(1+(single[[I]]/G)))
      
      
      if(rate=="yes"){
        #predict rate###############################################
        predict2=G*Ks*(time2^(G-1))
      }
      addoutput=rbind(addoutput,data.frame(groupid=groupdata,time=time2,observed=single[[I]],predict=predict2,Ks=Ks,G=G))
      
      i=i+1
    }
  }
  
  #equations for ungrouped data ############################
  if(is.null(group)){
    #prediction################################################
    Ks=coef(greenampt)[1]
    time2=data[[time]]
    I=data[[I]]
    time=data[[time]]
    ###########################################################
    G=coef(greenampt)[2]
    I=coef(greenampt)[3]##########################################
    ######################
    predict=Ks*(time2)+(G*log(1+(I/G)))
    predict3=Ks *((G/I)+1)
    
    if(!is.null(f)){
      predict=predict3
      I=data$f[[1]]
    }
    
    if(rate=="yes"){
      #####################################################
      predict=G*Ks*(time2^(G-1))
    }
  }
  else
  {
    ################################################
    G=addoutput$G
    I=data[[I]]
    time=data[[time]]
    ##############################################
    Ks=addoutput$Ks
    predict=addoutput$predict
    time=addoutput$time
  }
  #return varibales ########################################
  
factor<-list(greenampt=greenampt,data=data,time=time,I=I,Ks=Ks,G=G,group=group,
             predict=predict,rate=rate,formular=greenamptF,addoutput=addoutput,output=addoutput)

factor$call<-match.call()

class(factor)<-"greenampt"
factor
}
#' @export
#' @rdname greenampt
#predict function#########################
predict.greenampt<-function(object,time=NULL,...)
{
x<-object

if(is.null(object$group)){
predict=as.data.frame(cbind(x$time,x$predict))
names(predict)=c("time","predict")

if(!is.null(time)&object$rate=="yes"){####################
predict=object$G*object$Ks*(time^(object$G-1))
predict=predict[[1]]
}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$Ks*(time)+(object$G*log(1+(object$I/object$G)))

predict=predict[[1]]
}

print((predict))
}
else######################################
{
predict=object$addoutput

#rate######################################
if(!is.null(time)&object$rate=="yes"){
predict=object$G*object$Ks*(time^(object$G-1))
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}
#cumulative###################################
if(!is.null(time)&object$rate!="yes"){
predict=object$Ks*(time)+(object$G*log(1+(object$I/object$G)))
predict2= (data.frame(cbind(object$addoutput$groupid,predict)))
names(predict2)=c("groupid","predict")
predict=aggregate(predict2$predict,by=list(predict2$groupid),FUN=mean)
colnames(predict)=c("Group","Predict")
predict$Group=row.names(table(object$addoutput$groupid))

}

print(predict)
}
}

#plot function
#' @export
#' @rdname greenampt
plot.greenampt<-function(x,xlab="Time(Minutes)",ylab="Cumulative (cm)",main=NULL,layout=NULL,...)
{
object<-x
x<-object
G=x$G
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
title=paste(aggdata[i],"(R2=",round(r2,4),")")
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
#' @rdname greenampt
#summary function
summary.greenampt<-function(object,...)
{


x<-object$greenampt

if(is.null(object$group)){
summary1=summary(x)
print(summary1)
}
else
{########################
coef=aggregate(cbind(Ks, G) ~ groupid, data = object$addoutput, mean)##############
print(coef)
}

}

#print function
#' @export
#' @rdname greenampt
print.greenampt<-function(x,...)
{
object=x
if(is.null(object$group)){
x<-object$greenampt
print((x))
}
else
{######################################################
coef=aggregate(cbind(Ks, G) ~ groupid, data = object$addoutput, mean)
print(coef)
}


}
#' @export
#' @rdname greenampt
#coef function
coef.greenampt<-function(object,...)
{
x<-object$greenampt
if(is.null(object$group)){
coef=(coef(x))

}
else################################################
{
coef=aggregate(cbind(Ks, G) ~ groupid, data = object$addoutput, mean)
}
print(coef)
}

