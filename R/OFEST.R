#' @title The Offin Estimation of Soil Transfer parameters
#' @description This function uses \code{\link{philip}} or \code{\link{sw}} 
#' with \code{\link{lass3}}, \code{\link{ksat}} and \code{\link{vg}} to estimate 
#' water retention and hydraulic conductivity curves. See datails below.
#' 
#' @inheritParams BEST
#' @inheritParams lass3
#' @inheritParams ksat
#' @inheritParams vg
#' @param hg character "rawls" or "BEST" bubbling capillary pressure  algorithm. 
#' It takes either Rawls (1993) as in Dingman(2002) or BEST
#' @param model character. Either "philip" or "sw" (Swartzendruber) or "valiantzas" or "brutsaert".
#' The model for estimation of sorptivity (S) and hydraulic conductivity (Ks).
#' @param type character. Whether "linear" or "nonlinear" 
#' philip based equation to be applied.
#' @param K character. The type of hydraulic conductivity model. 
#' It takes "BC" for Brooks and Corey (1964) and VG for van Genuchten (1980).
#'
#' @return For the output of PSD see \code{\link{lass3}}.
#' For the output of goodness of fit tests see \code{\link{gof.OFEST}}. 
#' In the case of \code{cal.OFEST} the output can be assesed with mod object.
#' The return parameters can be assessed with \code{coef.OFEST}. 
#' The predicted soil water can be assessed  with \code{predict.OFEST} 
#' \itemize{
#' \item{Ks} { saturated hydraulic conductivity [LT^-1] }
#' \item{S} { Sorptivity  [LT^-0.5]}
#' \item{mod_theta,mod_thetal,mod_thetay,mod_thetab:} { The estimated soil water content at 
#' different metric potentials}
#' \item{theta:} { The observed soil water content at different metric potentials  }
#' \item{K:} { The unsaturated hydraulic conductivity [LT^-1] }
#' \item{hg:} { bubbling capillary pressure  [L] }
#' \item{q:} { infiltration rates [LT^-1] }
#' } 
#' @seealso \code{\link{BEST}}
#' @details The function first estimates sorptivity and Ks from 3D or early part of 
#' vertical cumulative infiltration data with \code{\link{philip}} or \code{\link{sw}}
#'  or \code{\link{valiantzas}} or \code{\link{brutsaert}}.
#' Second, \code{\link{lass3}} is used to estimate pore-size parameters, e.g. n, b. Mualem or
#' Burdine condition is applied to estimate m parameter of van Genuchten hydraulic conductivity 
#' curve. Third, a bubbling capillary pressure hb is estimated with Ks and S. 
#' Fourth, van Genuchten (1980) water retention curve and hydraulic conductivity 
#' are estimated. Finally, the predicted infiltration rate (q) is estimated.
#' @export
#' @author George Owusu
#' @references 
#' \itemize{
#' \item{}{Swartzendruber, D. (1987). A quasi-solution of Richards equation 
#' for the downward infiltration of water into soil. Water Resour Res, 23, 809-817.}
#' \item{}{Philip, J. R. (1957). The theory of infiltration:Sorptivity and 
#' algebraic infiltration equations. Soil Science, 84, 257-264.}
#' \item{}{Brutsaert, W. (1977). Vertical infiltration in dry soil. 
#' Water Resour. Res., 13, 363-368.}
#' \item{}{Valiantzas, J. D. (2010). New linearized two-parameter infiltration equation for 
#' direct determination of conductivity and sorptivity. Journal of Hydrology, 387. 
#' doi: 10.1016/j.jhydrol.2009.12.049}
#' \item{}{Rawls, W. J., Ahuja, L. R., Brakensiek, D. L., & Shirmohammadi, A. (1993). 
#' Infiltration and Soil Water Movement. 
#' In D. Maidment (Ed.), Handbook of Hydrology: McGraw-Hill Education.}
#' \item{}{Dingman, S. L. (2002). Physical Hydrology: Prentice-Hall.}
#' \item{}{van Genuchten, M. T. (1980). A Closed-form Equation for Predicting the
#'  Hydraulic Conductivity of Unsaturated Soils1. Soil Sci. Soc. Am. J., 44(5), 892-898. 
#'  doi: 10.2136/sssaj1980.03615995004400050002x}
#'  \item{}{Brooks, R. H., & Corey, A. T. (1964). Hydraulic properties of porous 
#'  medium Hydrology Paper (Vol. 3): Colorado State University.}
#' }
#' @examples
#' psd=read.csv(system.file("ext","sys","psd.csv",package="vadose"))
#' hf=read.csv(system.file("ext","sys","h.csv",package="vadose"))
#' data=read.csv(system.file("ext","sys","exampleBEST.csv",package="vadose"))
#' ofin1<-OFEST(data=data,time="time",I="I",h=hf$h,PSD=psd)
#' pred=predict(ofin1)
#' coefficients=coef(ofin1)
#' plot(ofin1,type="h")
#' plot(ofin1,type="all")
#' plot(ofin1,type="all",hlog=TRUE,klog=TRUE,kcol="green")
#' 
#' #with measured theta
#' ofin2<-OFEST(data=data,time="time",I="I",h=hf$h,PSD=psd,model="sw",theta=ofin1$mod_theta)
#' print(gof.OFEST(ofin2))
#' ofin1$K
#' \dontrun{
#' #calibration
#' ofin2.cal<-cal.OFEST(data=data,time="time",I="I",h=hf$h,PSD=psd,model="sw",theta=ofin1$mod_theta)
#' coef2=coef(ofin2.cal)
#' predict=predict(ofin2.cal,h=500)
#'  }
#generic function
OFEST<-function(data=NULL,time,I,n=NULL,m="b",pb=1.2,tho=0.169,
               thr=0,ths=NULL,theta=NULL,PSD=NULL,Ks=NULL,h=seq(0, 1500, by = 5)
               ,model="philip",type="nonlinear",hg="BEST",K="BC") UseMethod ("OFEST")

#' @rdname OFEST
#' @export
 OFEST.default<-function(data=NULL,time,I,n=NULL,m="b",pb=1.2,tho=0.169,
thr=0,ths=NULL,theta=NULL,PSD=NULL,Ks=NULL,h=NULL,
model="philip",type="nonlinear",hg="BEST",K="BC")
{
   if(is.null(h)){
     x=1:42
     h= 1.1163*exp(0.3699*x) 
   }
  thisK=K
   if(model=="all"||model=="All"||model=="ALL"){
     philipm=OFEST(data = data,   n = n, m = m,
                   pb = pb, tho = tho, thr = thr, PSD = PSD,
                   h = h, theta = theta, model = "philip", type = type,
                   ths = ths,hg=hg)
     #names(philipm)=paste("philip",names(philipm),sep="")
     sw=OFEST(data = data,   n = n, m = m,
                   pb = pb, tho = tho, thr = thr, PSD = PSD,
                   h = h, theta = theta, model = "sw", type = type,
                   ths = ths,hg=hg)
     #names(sw)=paste("sw",names(sw),sep="")
     
     valiantzasm=OFEST(data = data,   n = n, m = m,
              pb = pb, tho = tho, thr = thr, PSD = PSD,
              h = h, theta = theta, model = "valiantzas", type = type,
              ths = ths,hg=hg)
     #names(valiantzasm)=paste("valiantzas",names(valiantzasm),sep="")
     brutsaertm=OFEST(data = data,   n = n, m = m,
                       pb = pb, tho = tho, thr = thr, PSD = PSD,
                       h = h, theta = theta, model = "brutsaert", type = type,
                       ths = ths,hg=hg)
     #names(brutsaertm)=paste("brutsaert",names(brutsaertm),sep="")
     
     list(philip=philipm,sw=sw,valiantzas=valiantzasm,brutsaert=brutsaertm)
   }else{
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
     
     if(!is.null(data$n)){
       n=mean(data$n)
     }
     
     if(!is.null(data$pb)){
       pb=mean(data$pb)
     }
     
     if(!is.null(data$tho)){
       tho=mean(data$tho)
     }
     if(!is.null(data$thr)){
       thr=mean(data$thr)
     }
     if(!is.null(data$m)){
       m=data$m
     }
     if(!is.null(data$Ks)){
       Ks=data$Ks
     }
     if(!is.null(data$ths)){
       ths=data$ths
     }
     if(!is.null(data$model)){
       model=data$model
     }
     if(!is.null(data$type)){
       type=data$type
     }
     if(!is.null(data$hg)){
       hg=data$hg
     }
   }
     
     if(!is.null(names(h))){
       h3=h
       h=h3$h
       theta=h3$theta
     }
     
     if(!is.null(names(theta))){
       theta3=theta
       theta=theta3$theta
       h2=theta3$h
       if(is.null(h2)){
         h=h
       }else{
         h=h2
       }
     }
#lass3data=read.csv(file)
if(model=="philip"){
mod<-philip(data=data,time=time,I=I,type=type)
}

if(model=="sw"){
mod<-sw(data=data,time=time,I=I)
}
   
   if(model=="valiantzas"){
     mod<-valiantzas(data=data,time=time,I=I)
   }
   
   if(model=="brutsaert"){
     mod<-brutsaert(data=data,time=time,I=I)
   }
#gofthis=gof(mod)
#n from PSD
   S=mod$S
p=1-(pb/2.65)
if(is.null(Ks)){
Ks=mod$Ks
}
psd2=NULL
texture <- NULL
sand<-NULL
silt<-NULL
clay<-NULL
if(!is.null(n)){
  m=1-(2/n)
  b=1/(m*n)
  n2=(2/(m*n))+3
  
}

if(!is.null(PSD)&&is.null(n)){
psd2<-(lass3(PSD,p=p))
n<-psd2$n
b=psd2$b
n2<-psd2$n2
#m=mod1$m
texture <- texture(psd2)
sand<-psd2$sand
silt<-psd2$silt
clay<-psd2$clay
}
#mualem condition
thism=m
if(m=="m"||m=="Mualem"||m=="mualem"){
m=1-(1/n)
}
#burdine condition
if(m=="b"||m=="Burdine"||m=="burdine"){
m=1-(2/n)
}

if(is.null(ths)){
ths=1-(pb/2.65)
}
#print(hg)

if(hg=="rawls"){
hb=mod$S^2/((p-tho)*Ks*(((2*b)+3)/(b+3)))
hb=hb[[1]]
}else{
cp=gamma(1+(1/n))*(gamma((m*n2)-(1/n))/(gamma(m*n2))+(gamma((m*n2)+m-(1/n))/gamma((m*n2)+m)))
hb=S^2/(cp*(ths-tho)*((1-(tho/ths)^n2)*Ks))
hb=hb[[1]]
}
alp=hb^-1
#van Genuchten (1980)
#predict=thr+(ths-thr)/((1+(alp*h)^n)^(m))###########
predict=thr+((ths-thr)*((1+(h/hb)^n)^(-m)))


if(thisK=="vg"||thisK=="VG"||thisK=="van Genuchten"){
se=(predict-thr)/(ths-thr)
#print(se)
Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)^2
if(thism=="b"||thism=="Burdine"||thism=="burdine"){
  Kr=(se^0.5)*(1-(1-(se^(1/m)))^m)
}
K=Ks*Kr
}else{
  K=Ks*((predict-thr)/(ths-thr))^n2
  
}
#actual data
time2=data[[time]]
I2=data[[I]]
q=data.frame(t=numeric(),q=numeric())

t2=1
while(t2<=length(time2)){
if(t2==1){
time3=((time2[t2]^0.5+0^0.5)/2)^2
I3=(I2[t2]-0)/(time2[t2]-0)
}else{

time3=((time2[t2]^0.5+time2[t2-1]^0.5)/2)^2
I3=(I2[t2]-I2[t2-1])/(time2[t2]-time2[t2-1])
}

q=rbind(q,data.frame(t=time3,q=I3))
t2=t2+1

}


factor<-list(OFEST=mod,data=data,time=time,I=I,K=K,mod_theta=predict,Ks=Ks[[1]],
             S=S[[1]],theta=theta,h=h,q=q,psd=psd2,model=model,hg=hb,type=type,
             m=m,n=n,n2=n2,b=b,thr=thr,ths=ths,tho=tho,pb=pb)
factor$call<-match.call()

class(factor)<-"OFEST"
factor
}
}
 
 #' @export
 #' @rdname OFEST
 group.OFEST=function(data = NULL, time, I, n = NULL, m = "b", pb = 1.2, tho = 0.169,
                      thr = 0, ths = NULL, theta = NULL, PSD = NULL, Ks = NULL,
                      h = NULL, model = "all", type = "nonlinear",hg="BEST",
                      group,plot=TRUE,hlog=NULL,klog=NULL,layout=c(2,2),opar=par(mar=c(2,2,1.8,2))){
   
   if(is.null(h)){
     x=1:42
     h= 1.1163*exp(0.3699*x) 
   }
   
   if(plot==TRUE){
     #dev.new()
     opar1<-par()
     par(oma=c(3,1.5,1,1.5))
     par(mfrow=layout) 
     legend=FALSE
   }
   h2=NULL
   theta2=NULL
   parameters=data.frame(group=factor(), philipS=numeric(),philipKs=numeric(),
                         philiphg=numeric(),swS=numeric(),swKs=numeric(),
                         swhg=numeric(),valiantzasS=numeric(),valiantzasKs=numeric(),
                         valiantzashg=numeric(),brutsaertS=numeric(),brutsaertKs=numeric(),
                         brutsaerthg=numeric(), n=numeric(),m=numeric()
                         ,thr=numeric(),ths=numeric(),tho=numeric())
   parameters=data.frame()
   
   gof1=data.frame()
   gof1box=data.frame()
   predict=data.frame(group=factor(),mod_theta_philip=numeric(),mod_theta_sw=numeric(),
                      mod_theta_valiantzas=numeric(),mod_theta_brutsaert=numeric(),
                      K_philip=numeric(),K_sw=numeric(),K_valiantzas=numeric(),
                      K_brutsaert=numeric())
   
   predict=data.frame()
   predict=data.frame()
   aggdata =row.names(table(data[group]))
   coefall=NULL
   thisgroup="nothing145lhsgfsfs54ua"
   #if(!is.null(data$n)){
   #  n=aggregate(data$n,data[group],FUN=mean)
   #}
   for(i in 1:length(aggdata)) {
     
     single=data[data[group]==aggdata[i],]
     if(!is.null(PSD)){
       PSD2=PSD[PSD[group]==aggdata[i],]
     }
     if(is.null(n)){
       n2=n
       if(!is.null(single$n)){
         n2=single$n[1]
         # print(n)
       }
     }else{
       if(length(n)==length(aggdata)){
         n2=n[i]
       }else{
         n2=n
       }
     }
     
     
       if(!is.null(single$m)){
         m2=single$m[1]
       
     }else{
       if(length(m)==length(aggdata)){
         m2=m[i]
       }else{
         m2=m
       }
     }  
     
     
     
     
       if(!is.null(single$pb)){
         pb2=single$pb[1]
     }else{
       if(length(pb)==length(aggdata)){
         pb2=pb[i]
       }else{
         pb2=pb
       }
     }
     
     
       if(!is.null(single$tho)){
         tho2=single$tho[1]
       
     }else{
       if(length(tho)==length(aggdata)){
         tho2=tho[i]
       }else{
         tho2=tho
       }
     }
     
     
       if(!is.null(single$thr)){
         thr2=single$thr[1]
       
     }else{
       if(length(thr)==length(aggdata)){
         thr2=thr[i]
       }else{
         thr2=thr
       }
     }
     
     if(is.null(ths)){
       ths2=ths
       if(!is.null(single$ths)){
         ths2=single$ths[1]
       }
     }else{
       if(length(ths)==length(aggdata)){
         ths2=ths[i]
       }else{
         ths2=ths
       }
     }
     
     
     
     
     if(!is.null(names(h))){
       h4=h[h[group]==aggdata[i],]
       h2=h4$h
       theta2=h4$theta
     }
     
     if(!is.null(names(theta))&&!is.null(theta)){
       theta4=theta[theta[group]==aggdata[i],]
       
       if(is.null(h2)){
         h2=theta4$h
       }
       if(is.null(theta2)){
         theta2=theta4$theta
       }
     }
     
     if(is.null(h2)){
       h2=h
     }
     
     if(is.null(theta2)){
       theta2=theta
     }
    if(is.null(theta2)&&is.null(h2)){
      h2=seq(0, 15000, by = 10)
      theta2=NULL
    }
     mod=OFEST(data = single,   n = n2, m = m2,
               pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
               h = h2, theta = theta2, model = model, type = type,
               ths = ths2,hg=hg)
     
     
     if(!is.null(mod$philip)){
       
       coefPH=coef.OFEST(mod$philip)
       coefSW=coef.OFEST(mod$sw)
       coefva=coef.OFEST(mod$valiantzas)
       coefbr=coef.OFEST(mod$brutsaert)
       #print((mod$philip))
       gofPH=vadose.tryCatch(gof.OFEST(mod$philip))$value
       gofSW=vadose.tryCatch(gof.OFEST(mod$sw))$value
       gofva=vadose.tryCatch(gof.OFEST(mod$valiantzas))$value
       gofbr=vadose.tryCatch(gof.OFEST(mod$gofbr))$value
       #print(gofSW)
       if(class(gofPH)[1]=="simpleError"){
         gofPH=NA 
       }
       if(class(gofSW)[1]=="simpleError"){
         gofSW=NA 
       }
       if(class(gofva)[1]=="simpleError"){
         gofva=NA 
       }
       if(class(gofbr)[1]=="simpleError"){
         gofbr=NA 
       }
       
       
       #gofSW=gof.OFEST(mod$sw)
       #gofva=gof.OFEST(mod$valiantzas)
       #gofbr=gof.OFEST(mod$brutsaert)
       gofall=cbind(gofPH,gofSW,gofva,gofbr)
       #print(gofPH)
       coefall=cbind(coefPH,coefSW,coefva,coefbr)
       
       coefall2=cbind(aggdata[i],coefall[1,],coefall[2,],coefall[3,],coefall[4,])
       if(type!="linear"){
       gof4=rbind(gofPH,gofSW,gofva) 
       gof4$model=toupper(row.names(gof4))

       }else{
         gof4=rbind(gofPH,gofva)
         gof4$model=paste(toupper(row.names(gof4)),".L",sep="")
       }
       gof4$ID=aggdata[i]
       gof1box=rbind(gof1box,data.frame(gof4))
       #print(gof4)
       gof1=rbind(gof1,data.frame(group=aggdata[i],philip=gofPH,sw=gofSW,valiantzas=gofva,
                                  brutsaert= gofbr))
      # }
       parameters=rbind(parameters,data.frame(group=aggdata[i], 
                                              philipS=coefPH$S,philipKs=coefPH$Ks,
                                              philiphg=coefPH$hg,swS=coefSW$S,swKs=coefSW$Ks,
                                              swhg=coefSW$hg,valiantzasS=coefva$S,valiantzasKs=coefva$Ks,
                                              valiantzashg=coefva$hg,brutsaertS=coefbr$S,brutsaertKs=coefbr$Ks,
                                              brutsaerthg=coefbr$hg, n=mod$philip$n,m=mod$philip$m
                                              ,thr=mod$philip$thr,ths=mod$philip$ths,tho=mod$philip$tho))
       
       
       predict=rbind(predict,data.frame(group=aggdata[i],mod_theta_philip=mod$philip$mod_theta,
                                        mod_theta_sw=mod$sw$mod_theta,mod_theta_valiantzas=mod$valiantzas$mod_theta,
                                        mod_theta_brutsaert=mod$brutsaert$mod_theta,
                                        K_philip=mod$philip$K,K_sw=mod$sw$K,K_valiantzas=mod$valiantzas$K,
                                        K_brutsaert=mod$brutsaert$K,q=mean(mod$philip$q$q)))
     }else{
       
       predict=rbind(predict,data.frame(group=aggdata[i],mod_theta=mod$mod_theta,K=mod$K))
       coefPH=mod
       thisgof=vadose.tryCatch(gof.OFEST(mod))$value
       if(class(thisgof)[1]=="simpleError"){
         thisgof=NA 
       }
       gof1=rbind(gof1,data.frame(group=aggdata[i],thisgof))                                                                        
       parameters=rbind(parameters,data.frame(group=aggdata[i], 
                                              S=coefPH$S,Ks=coefPH$Ks,
                                              hg=coefPH$hg, n=mod$n,m=mod$m
                                              ,thr=mod$thr,ths=mod$ths,tho=mod$tho))
     }
     
     
     if(plot==TRUE){
       if(model=="all"|model=="All"||model=="ALL"){
         main1=aggdata[i]
       }else{
       main1=aggdata[i]
       }
      thisplot=vadose.tryCatch( plot.OFEST(mod,type="hk",main=main1,legend=legend,hlog=hlog,klog=klog,opar=opar))
     }
     
     
     #coef3=coef(mod)
     #print(aggdata[i])
     #print(single$ID[[i]])
     #print(data.frame(coef3$Ks))
     
   }
   
   if(plot==TRUE){
   col=c("blue", "red", "darkgreen","white")
   plot_colors=col
   par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
   plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
   #legend("bottom", c("IM", "IBD", "1R", "2R"), xpd = TRUE, horiz = TRUE, 
     #     inset = c(0,0), bty = "n", pch = c(4, 2, 15, 19), col = 1:4, cex = 2)
   if(is.null(theta2)){
     if(!is.null(mod$philip)){
     legend("bottom", (c("Philip","Swartzendruber","Valiantzas","Water Content")), 
            cex=0.9,col=plot_colors,lty=1:3,lwd=2, bty="n",xpd = TRUE, horiz = TRUE, 
            inset = c(0,0))
     }else{
       legend("bottom", (c("h","K")), 
              cex=0.9, 
              col=c("black","green"),lty=c(1,1),lwd=2, bty="n",xpd = TRUE, horiz = TRUE, 
              inset = c(0,0))
     }
   }else{
     if(is.null(mod$philip)){
       
   legend("bottom", (c("Measured","h","K")), 
          cex=0.9, pch=c("o","",""),
          col=c("black","black","green"),lty=c(0,1,1),lwd=2, bty="n",xpd = TRUE, horiz = TRUE, 
          inset = c(0,0))
     }else{
       legend("bottom", (c("Measured","Philip","Swartzendruber","Valiantzas","Water   Content")),
              cex=0.9, pch=c("o","","","",""),col=c("black",plot_colors),lty=c(0,1,2,3,4),lwd=2, bty="n"
              ,xpd = TRUE, horiz = TRUE, inset = c(0,0))  
       #text("top","Hooooooooooooooooooooo")
     }
   }
   par(opar1)
   }
   row.names(gof1)=NULL
   tabname=names(gof1)
   tabname2=gsub(".1", "_THETA", tabname, ignore.case =TRUE, fixed=FALSE)
   names(gof1)=tabname2
   if(plot==TRUE){
     par(opar1)
   }
   #row.names(predict)=predict$group
   #predict$group=NULL
   #print(gof1box)
   list(parameters=parameters,predict=predict,statistics=gof1,statistics2=gof1box)
 }
 
 #' @export
 #' @rdname OFEST
 cal.OFEST=function(data = NULL, time, I, n = NULL, m = "b", pb = c(1, 1.8,0.1), 
                    tho = c(0.01,0.8,0.05), thr =c(0,0.2,0.05), ths = NULL, 
                    theta = NULL, PSD = NULL, Ks = NULL,hg="BEST",
                    h = NULL, model = "all", type = "nonlinear",
                    group=NULL){
   #print(group[1])
   if(is.null(h)){
     x=1:42
     h= 1.1163*exp(0.3699*x) 
   }
   if(!is.null(group)){
     parameters=data.frame()
     parameters2=data.frame()
     fit=data.frame()
     fit4=data.frame()
     #group="TownName"
     aggdata =row.names(table(data[group]))
     h2=NULL
     theta2=NULL
     for(i in 1:length(aggdata)) {
       
       single=data[data[group]==aggdata[i],]
       if(!is.null(PSD)){
         PSD2=PSD[PSD[group]==aggdata[i],]
       }
       if(is.null(n)){
         n2=n
         if(!is.null(single$n)){
           n2=single$n[1]
           # print(n)
         }
       }else{
         if(length(n)==length(aggdata)){
           n2=n[i]
         }else{
           n2=n
         }
       }
       
       
       if(!is.null(single$m)){
         m2=single$m[1]
         
       }else{
         if(length(m)==length(aggdata)){
           m2=m[i]
         }else{
           m2=m
         }
       }  
       
       
       
       
       if(!is.null(single$pb)){
         pb2=single$pb[1]
       }else{
         if(length(pb)==length(aggdata)){
           pb2=pb[i]
         }else{
           pb2=pb
         }
       }
       
       
       if(!is.null(single$tho)){
         tho2=single$tho[1]
         
       }else{
         if(length(tho)==length(aggdata)){
           tho2=tho[i]
         }else{
           tho2=tho
         }
       }
       
       
       if(!is.null(single$thr)){
         thr2=single$thr[1]
         
       }else{
         if(length(thr)==length(aggdata)){
           thr2=thr[i]
         }else{
           thr2=thr
         }
       }
       
       if(is.null(ths)){
         ths2=ths
         if(!is.null(single$ths)){
           ths2=single$ths[1]
         }
       }else{
         if(length(ths)==length(aggdata)){
           ths2=ths[i]
         }else{
           ths2=ths
         }
       }
       
       if(!is.null(names(h))){
         h4=h[h[group]==aggdata[i],]
         h2=h4$h
         theta2=h4$theta
       }
       
       if(!is.null(names(theta))&&!is.null(theta)){
         theta4=theta[theta[group]==aggdata[i],]
         
         if(is.null(h2)){
           h2=theta4$h
         }
         if(is.null(theta2)){
           theta2=theta4$theta
         }
       }
       
       if(is.null(h2)){
         h2=h
       }
       
       if(is.null(theta2)){
         theta2=theta
       }
       
       
       #print(c(n2,tho2,thr2))
       if(model=="all"||model=="All"||model=="ALL"){
         print(aggdata[i])
         
         modph=cal.OFEST(data = single,   n = n2, m = m2,
                       pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
                       h = h2, theta = theta2, model = "philip", type = type,
                       ths = ths2,hg=hg)
         modsw=cal.OFEST(data = single,   n = n2, m = m2,
                         pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
                         h = h2, theta = theta2, model = "sw", type = type,
                         ths = ths2,hg=hg)
         modva=cal.OFEST(data = single,   n = n2, m = m2,
                         pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
                         h = h2, theta = theta2, model = "valiantzas", type = type,
                         ths = ths2,hg=hg)
         modbru=cal.OFEST(data = single,   n = n2, m = m2,
                         pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
                         h = h2, theta = theta2, model = "brutsaert", type = type,
                         ths = ths2,hg=hg)
        # if(i==1){
           #names(modph$parameters)=paste("philip",names(modph$parameters),sep=".") 
           #names(modsw$parameters)=paste("sw",names(modsw$parameters),sep=".") 
           #names(modva$parameters)=paste("valiantzas",names(modva$parameters),sep=".") 
           #names(modbru$parameters)=paste("brutsaert",names(modbru$parameters),sep=".") 
         #}
         
         thisparameters=data.frame(c(modph$parameters,modsw$parameters,modva$parameters))
         #thisparameters=cbind(philip=modph$parameters,sw=modsw$parameters,valiantzas=modva$parameters)
         thisparameters$group=aggdata[i]
         #print(thisparameters)
         #print(parameters)
         parameters=rbind(parameters,data.frame(thisparameters))
         
         thisparameters2=data.frame(c(modph$parameters2,modsw$parameters2,modva$parameters2))
         #thisparameters=cbind(philip=modph$parameters,sw=modsw$parameters,valiantzas=modva$parameters)
         thisparameters2$group=aggdata[i]
         #print(thisparameters)
         #print(parameters)
         parameters2=rbind(parameters2,data.frame(thisparameters2))
         
         thisfit=data.frame(c(modph$fit,modsw$fit,modva$fit))
         thisfit$group=aggdata[i]
         fit=rbind(fit,data.frame(thisfit))
         
         thisfit4=data.frame(rbind(modph$fit,modsw$fit,modva$fit))
         thisfit4$group=aggdata[i]
         thisfit4$model=row.names(thisfit4)
         fit4=rbind(fit4,data.frame(thisfit4))
         #fit4$model=row.names(fit4)
       }else
         {
       mod=cal.OFEST(data = single,   n = n2, m = m2,
                 pb = pb2, tho = tho2, thr = thr2, PSD = PSD2,
                 h = h2, theta = theta2, model = model, type = type,
                 ths = ths2,hg=hg)
       #print(mod$parameters)
       if(!is.null(mod$parameters)){
         
      thisparameters=mod$parameters
      thisparameters$group=aggdata[i]
       parameters=rbind(parameters,data.frame(thisparameters))
       print(aggdata[i])
       
       }
       if(!is.null(mod$parameters2)){
         
         thisparameters2=mod$parameters2
         thisparameters2$group=aggdata[i]
         parameters2=rbind(parameters2,data.frame(thisparameters2))
        # print(aggdata[i])
         
       }
       if(!is.null(mod$fit)){
        # print(mod$fit)
       fit=rbind(fit,data.frame(cbind(aggdata[i],mod$fit)))
       #fit4=rbind(fit4,data.frame(rbind(aggdata[i],mod$fit)))
       }
     }
       #plot(mod$mod)
   }
   if(model=="all"||model=="All"||model=="ALL"){
     row.names(parameters)=parameters$group
     parameters$group=NULL
     row.names(parameters2)=parameters2$group
     parameters2$group=NULL
     #print(fit)
     row.names(fit)=fit$group
     
     fit$group=NULL
     #column names
     tabname=names(fit)
     tabname2=gsub(".1", ".Swartzendruber.Infiltration", tabname, ignore.case =TRUE, fixed=TRUE)
     names(fit)=tabname2
     
     tabname=names(fit)
     tabname2=grep("Swartzendruber", tabname, ignore.case =FALSE, fixed=FALSE)
     thisname3=paste(tabname[1:length(tabname2)],".Philip.Infiltration",sep=".")
     tabname[1:length(tabname2)]=thisname3
     names(fit)=tabname
     
     tabname=names(fit)
     tabname2=gsub(".2", ".Valiantzas.Infiltration", tabname, ignore.case =TRUE, fixed=TRUE)
     names(fit)=tabname2
     
    
     
     tabname=names(fit)
     tabname2=gsub(".3", ".Philip.THETA", tabname, ignore.case =TRUE, fixed=TRUE)
     names(fit)=tabname2
     
     tabname=names(fit)
     tabname2=gsub(".4", ".Swartzendruber.THETA", tabname, ignore.case =TRUE, fixed=TRUE)
     names(fit)=tabname2
     
     tabname=names(fit)
     tabname2=gsub(".5", ".Valiantzas.THETA", tabname, ignore.case =TRUE, fixed=TRUE)
     names(fit)=tabname2
     
     #for paramters names
     tabname=names(parameters)
     tabname2=gsub(".1", ".Swartzendruber", tabname, ignore.case =TRUE, fixed=TRUE)
    names(parameters)=tabname2
     
     tabname=names(parameters)
     tabname2=gsub(".2", ".valiantzas", tabname, ignore.case =TRUE, fixed=TRUE)
     names(parameters)=tabname2
    
     tabname=names(parameters)
     tabname2=grep("Swartzendruber", tabname, ignore.case =FALSE, fixed=FALSE)
     thisname3=paste(tabname[1:length(tabname2)],"Philip",sep=".")
     tabname[1:length(tabname2)]=thisname3
     names(parameters)=tabname
     
     #for parameters 2 names
     tabname=names(parameters2)
     tabname2=gsub(".1", ".Swartzendruber", tabname, ignore.case =TRUE, fixed=TRUE)
     names(parameters2)=tabname2
     
     tabname=names(parameters2)
     tabname2=gsub(".2", ".valiantzas", tabname, ignore.case =TRUE, fixed=TRUE)
     names(parameters2)=tabname2
     
     tabname=names(parameters2)
     tabname2=grep("Swartzendruber", tabname, ignore.case =FALSE, fixed=FALSE)
     thisname3=paste(tabname[1:length(tabname2)],"Philip",sep=".")
     tabname[1:length(tabname2)]=thisname3
     names(parameters2)=tabname
   }else{
     row.names(parameters)=parameters$group
     parameters$group=NULL
     row.names(parameters2)=parameters2$group
     parameters2$group=NULL
     row.names(fit)=fit[[1]]
     
     tabname=names(fit)
     tabname2=gsub(".1", "_THETA", tabname, ignore.case =TRUE, fixed=FALSE)
     names(fit)=tabname2
     fit[1]=NULL
   }
     if(type=="linear"){
       fit4=vadose.tryCatch(fit4[which(fit4$model!="sw"),])$value
       fit4$model=paste(fit4$model,".L",sep="")
     }
     
     list(parameters=parameters,parameters2=parameters2,statistics=fit,fit=fit,statistics2=fit4)
   }else{
   
   if(is.null(theta)){
     return  (print("Observation or measured data 'theta' is needed"))
   }
   
   para=NULL
   count=0 
   
   
   if(length(tho)==1){
     tho=c(tho,tho,tho)
   }else{
     para=cbind(para,tho=tho)
     count=count+1
   }
   
   
   if(length(thr)==1){
     thr=c(thr,thr,thr)
   }else{
     para=cbind(para,thr=thr)
     count=count+1
   }
   n4=TRUE
   if(is.null(n)){
     n4=NULL
     n=c(1,1,1)
   }else{
     if(length(n)==1){
       n=c(n,n,n)
     }else{
       para=cbind(para,n=n)
       count=count+1
     }
   }
   
   
   if(length(pb)==1){
     pb=c(pb,pb,pb)
     
   }else{
     para=cbind(para,pb=pb)
     count=count+1
   }
   
   print("Iteration Started. Kindly wait.................")
   pnames=names(data.frame(para))
   count=length(pnames)
   #return(list(para=para))
   #print(count)
   rmse=100
   r2=0
   
   pb2=seq(pb[1],pb[2],pb[3])
   tho2=seq(tho[1],tho[2],tho[3])
   thr2=seq(thr[1],thr[2],thr[3])
   n2=seq(n[1],n[2],n[3])
   fit2=NULL
   para2=NULL
   mod=NULL
   mod2=NULL
   parameters=NULL
   for(tho1 in tho2) {
     for(thr1 in thr2){
       for(pb1 in pb2){
         for(n1 in n2){
           if(is.null(n4)){
             #print("YES")
             n3=NULL
           }else{
             n3=n1
           }
           mod=OFEST(data = data, time=time, I=I, n = n3, m = m, pb= pb1, tho=tho1,
                     thr = thr1, ths = ths, theta = theta, PSD = PSD, Ks = Ks,
                     h = h, model = model, type = type)
           #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
           #print(mod)
           fit=vadose.tryCatch(gof.OFEST(mod,P=count))$value
           if(class(fit)[1]!="simpleError"){
           #fit=gof.OFEST(mod,P=count)
           
             fitrmse=fit$RMSE
           
           fitr2=fit$r2
           if(fitrmse<rmse){
             #print(cbind(WP1,FC1,kcini1,kcmid1,kcend1,p1,Ze1,REW1,CN1,FC1,initgw1,a11,b11,a21,b21,a31,b31,a41,b41,rc1))   
             parameters=(list(tho=tho1,n=n3,thr=thr1,pb=pb1)) 
             para2=parameters[pnames]
             print(data.frame(para2))
             print(fit)
             r2=fitr2
             fit2=fit
             rmse=fitrmse
             mod2=mod
              }
           
           if(fitrmse==0){
             break
           }
           #parameters=(list(tho=tho1,n=n3,thr=thr1,pb=pb1)) 
           #para2=parameters[pnames]
           #print(data.frame(para2))
           #print(cbind(tho1,thr1,pb1))
           #print(fitrmse)
           #print(rmse)
           #print(fit$d)
           }
         }
       }
     }
   }
   
   list(fit=fit2,parameters=parameters,parameters2=para2,mod=mod2)
   }
 }
 
 #' @export
 #' @rdname OFEST
 gof.OFEST<-function(x=NULL,obs=NULL,est=NULL,P=NULL)
 {
   if(!is.null(x$philip)){
     
      ph=gof.OFEST(x$philip) 
      sw=gof.OFEST(x$sw) 
      va=gof.OFEST(x$valiantzas) 
      br=gof.OFEST(x$brutsaert) 
      list(philip=ph,sw=sw,valiantzas=va,brutsaert=br)
     #print("")
   }else{
   gofphilip=NULL
   theta=  x$theta
   
   
   #print(gofphilip)
   gofhb=NULL
   
   if(is.null(theta)){
   if(x$type=="nonlinear"){
    data=gof(x$OFEST)
  row.names(data)=x$model
   }else{
     data=gof(x$OFEST)
     row.names(data)=x$model
     #print(row.names(data))
   }
   }
   
   if(!is.null(theta)){  
     theta=data.frame(x$theta)
     #print(data.frame(theta))
     data=(gof(x,c=theta,r=x$mod_theta,P=P))
     row.names(data)=x$model
     
   }
  
   #print(data)
   data
   #list(gof=data)
 }
 }  
 #' @export
 #' @rdname OFEST
 coef.OFEST<-function(object,...){
   
   if(is.null(object$Ks)){
     object=object$mod
   }
   
   n=object$n
   m=object$m
   pb=object$pb
   thr=object$thr
   ths=object$ths
   tho=object$tho
   Ks=object$Ks
   S=object$S
   hg=object$hg
   n2=object$n2
   
   list(S=S,Ks=Ks,hg=hg,n=n,m=m,n2=n2,pb=pb,thr=thr,ths=ths,tho=tho)
   
   
 }
 
 #' @export
 #' @rdname OFEST
 predict.OFEST<-function(object,h=NULL, ...)
 {
   if(is.null(object$Ks)){
     object=object$mod
   }
   if(is.null(h)){
     theta=object$mod_theta
     K=object$K
     
   }else{
     theta=object$thr+((object$ths-object$thr)*((1+(h/object$hg)^object$n)^(-object$m)))
     Kr=(((theta-object$thr)/(object$ths-object$thr))^object$n2)
     K=Kr*object$Ks
   }
   
   list(theta=theta,K=K)
   
   
 }
 
#default function
#plot function
 #' @rdname OFEST
 #' @export
plot.OFEST<-function(x,main=NULL,xlab="Water Content",
                     ylab="",ylab2="",
                     layout=NULL,hlog=NULL,klog=NULL,kcol="darkgreen",
                     col=c("blue", "red", "darkgreen","gold"),
                     units=c("s","mm"),mfrow=c(2,2),type="all",legend=TRUE,opar=par(mar=c(2,2,1.5,2)),...)
{
  if(!is.null(x$philip)){
    object=x
    #par(mfrow=c(2,2),mar=c(4, 4, 4, 4))de
    Kb=object$philip$K
    Kl=object$sw$K
    Ky=object$valiantzas$K
    Kz=object$brutsaert$K
    h=object$philip$h
    mod_thetab=object$philip$mod_theta
    mod_thetal=object$sw$mod_theta
    mod_thetay=object$valiantzas$mod_theta
    mod_thetaz=object$brutsaert$mod_theta
    
    theta=object$philip$theta
    obs=object$philip$data[[object$philip$I]]
    estb=object$philip$OFEST$predict
    estl=object$sw$OFEST$predict
    esty=object$valiantzas$OFEST$predict
    estz=object$brutsaert$OFEST$predict
    
    time=object$philip$data[[object$philip$time]]
    
    ###########################
    plot_colors=col
    #par(mar = c(5, 4, 4, 4) + 0.3)
    opar
    
    #par(mar=c(4.2, 3.8, 0.2, 0.2))
    #plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen")
    
    #ylab=paste(ylab,"(",units[2],")",sep="")
    #ylab2=paste(ylab2,"(",units[2],"/",units[1],")",sep="")
    
    ymax <- max(Kb, na.rm=TRUE)
    #dev.new(width=15, height=10)
    if(type=="all"){
      par(mfrow=mfrow)
    }
    
    if(type=="all"||type=="h"||type=="hk"||type=="k"){
      if(!is.null(theta)){
        if(!is.null(hlog)){
          plot(theta, h,xlab=xlab,ylab=ylab,log="y",main=main,...) # first plot
        }else{
          plot(theta, h,xlab=xlab,ylab=ylab,main=main,...) # first plot
          
        }
      }
      if(!is.null(hlog)){
        #ylab=paste(ylab,"[log]")
        if(!is.null(theta)){
          lines(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],log="y",...)
        }else{
          plot(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],main=main,log="y",...)
        }
        lines(mod_thetal,h,xlab=xlab,ylab=ylab,type="l", lty=2, lwd=2,col=plot_colors[2],log="y",...)
        lines(mod_thetay,h,xlab=xlab,ylab=ylab, type="l", lty=3, lwd=2,col=plot_colors[3],log="y",...)
        if(legend!=FALSE){
        lines(mod_thetaz,h,xlab=xlab,ylab=ylab, type="l", lty=4, lwd=2,col=plot_colors[4],log="y",...)
        }
        if(is.null(theta)){
          if(legend!=FALSE){
          legend("bottom", (c("Philip","Swartzendruber","Valiantzas","Brutsaert")), cex=0.8, col=plot_colors,lty=1:3, lwd=2, bty="n")
          }
            }else{
              if(legend!=FALSE){
          legend("bottom", (c("Measured","Philip","Swartzendruber","Valiantzas","Brutsaert")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
              }
                }
      }else{
        if(!is.null(theta)){
          lines(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],...)
          
        }else{
          plot(mod_thetab,h,type="l",xlab=xlab,ylab=ylab,col=plot_colors[1],main=main,...)
        }
        
        lines(mod_thetal,h,type="l",lty=2, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[2],...)
        lines(mod_thetay,h,type="l",lty=3, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[3],...)
        if(legend!=FALSE){
          lines(mod_thetaz,h,type="l",lty=4, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[4],...)
        }
          if(is.null(theta)){
          if(legend!=FALSE){
          legend("top", (c("Philip","Swartzendruber","Valiantzas","Brutsaert")), cex=0.8, col=plot_colors,lty=1:3, lwd=2, bty="n")
          }
            }else{
          if(legend!=FALSE){
          legend("top", (c("Measured","Philip","Swartzendruber","Valiantzas","Brutsaert")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
          } 
        }
      }
      par(new = TRUE)
      if(!is.null(klog)){
        plot(mod_thetab, Kb, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[1],log="y",...)
        lines(mod_thetal, Kl, type="l",lty=2, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[2],log="y",...)
        lines(mod_thetay, Ky, type="l",lty=3, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[3],log="y",...)
        lines(mod_thetaz, Kz, type="l",lty=4, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[4],log="y",...)
        #ylab2=paste(ylab2,"[log]")
      }else{
        plot(mod_thetab, Kb, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[1],...)
        lines(mod_thetal, Kl,  type="l",lty=2, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[2],...)
        lines(mod_thetay, Ky, type="l",lty=3, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[3],...)
        lines(mod_thetaz, Kz, type="l",lty=4, lwd=2, axes = FALSE, bty = "n", xlab = "", ylab = "",col=plot_colors[3],...)
      }
      
      
      
      axis(side=4, ylim=c(0, ymax))
      mtext(ylab2, side=4, line=2,cex=0.8)
      index1=length(mod_thetay)
      index2=length(mod_thetay)-1
      leb=(mod_thetay[index2]-mod_thetay[index1])
      text(max(mod_thetay)-(0*6),max(Kb),"K")
      text(min(mod_thetay)+(0*6),max(Kb),"h")
      
    }
    
    if(type=="all"||type=="cum"||type=="I"||type=="I"){
      plot(time,obs,ylab="",main=paste("Cumulative infiltration(",units[2],")",sep=""),xlab=xlab)
      par(new = TRUE)
      #lines(time,estb,type="l",lty=1, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[1])
      lines(time,estl,type="l",lty=2, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[2])
      lines(time,esty,type="l",lty=3, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[3])
      lines(time,estz,type="l",lty=4, lwd=2,xlab=xlab,ylab=ylab,col=plot_colors[4])
      
      #legend("bottomright", (c("Measured","Bagarello","Lassabatere","Yilmaz")), cex=0.8, pch=c("o","","",""),col=c("black",plot_colors),lty=c(0,1,2,3),lwd=2, bty="n")
      if(legend!=FALSE){
      legend("bottomright", (c("Measured","Philip","Swartzendruber","Valiantzas","Brutsaert")), cex=0.8, pch=c("o","",""),col=c("black",plot_colors[2],plot_colors[3]),lty=c(0,1,2,3),lwd=2, bty="n")
      }
        }
    if(type=="all"||type=="i"||type=="f"){
      plot(object$philip$q,main="Infiltration rate",xlab=xlab,ylab=ylab)
    }
    
    if(type=="all"||type=="psd"||type=="PSD"){
      plot(object$philip$psd)
    }
    if(type=="all"){
      par(mfrow=c(1,1))
    }  
  }else{
  object<-x
K=object$K
h=object$h
mod_theta=object$mod_theta
theta=object$theta
obs=object$data[[object$I]]
time=object$data[[object$time]]

if(!is.null(units)){
ylab=paste(ylab,"(",units[2],")",sep="")
ylab2=paste(ylab2,"(",units[2],"/",units[1],")",sep="")
}else{
  ylab=ylab
  ylab2=ylab2
}
#print(ylab2)
par(mar = c(5, 4, 4, 4) + 0.3)
col="black"
ymax <- max(K, na.rm=TRUE)

if(type=="all"){
par(mfrow=mfrow)
}

if(type=="all"||type=="h"||type=="hk"||type=="k"){

if(!is.null(theta)){
if(!is.null(hlog)){
plot(theta, h,xlab=xlab,ylab=ylab,log="y",main=main,...) # first plot
}else{
plot(theta, h,xlab=xlab,ylab=ylab,main=main,...) # first plot

}
}
if(!is.null(hlog)){
#ylab=paste(ylab,"[log]")
if(!is.null(theta)){
lines(mod_theta,h,type="l",xlab=xlab,ylab=ylab,col=col,log="y",...)
}else{
plot(mod_theta,h,type="l",xlab=xlab,ylab=ylab,main=main,log="y",cex=0.8,...)
}

}else{
if(!is.null(theta)){
lines(mod_theta,h,type="l",xlab=xlab,ylab=ylab,col=col,...)

}else{
#plot(mod_theta,h,type="l",xlab=xlab,ylab=ylab,col=col,...)
}
if(!is.null(theta)){
lines(mod_theta,h,type="l",xlab=xlab,ylab=ylab,col=col,...)

}else{
plot(mod_theta,h,type="l",xlab=xlab,ylab=ylab,col=col,main=main,...)
}
}
par(new = TRUE)
if(!is.null(klog)){
plot(mod_theta, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",log="y",col=kcol,...)
#ylab2=paste(ylab2,"[log]")
}else{
plot(mod_theta, K, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",col=kcol,...)
}

axis(side=4, ylim=c(0, ymax))
mtext(ylab2, side=4, line=2,cex=0.8)
index1=length(mod_thetay)
index2=length(mod_thetay)-1
leb=(mod_thetay[index2]-mod_thetay[index1])
text(max(mod_thetay)-(0*6),max(Kb),"K")
text(min(mod_thetay)+(0*6),max(Kb),"h")

}

if(type=="all"||type=="psd"||type=="PSD"){
  if(!is.null(object$psd)){
plot(object$psd)
}
}  
if(type=="all"||type=="cum"||type=="I"||type=="I"){
#plot(time,obs,ylab="",main=paste("Cumulative infiltration(",units[2],")",sep=""),xlab=paste("Time(",units[1],")",sep=""))
}

if(type=="all"||type=="i"||type=="f"){
plot(object$q,main="Infiltration rate",xlab=xlab,ylab=ylab)
}

if(type=="all"||type=="model"){
plot(object$OFEST,main=paste("Cumulative infiltration(",units[2],")",sep=""),xlab=xlab)
}

if(type=="all"){
par(mfrow=c(1,1))
}
}
}
#ofin1<-OFEST(data=data,time="time",I="I",h=hf$h,PSD=psd)
#ofin2<-OFEST(data=data,time="time",I="I",h=hf$h,PSD=psd,model="sw",theta=ofin1$mod_theta)
#print(gof(ofin2))
#print(gof(ofin1)

#plot(ofin1,type="h")
#plot(ofin1,type="all")
#plot(ofin1,type="all",hlog=T,klog=T,kcol="green")

#ofin1$K
