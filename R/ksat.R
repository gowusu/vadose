#' @title Estimating saturated hydraulic conductivity with water retention parameters
#' @description This function estimates Ks based on 5 models that include (see references) 
#' Guarracino (2007), Nasta et al (2013), Mishra et al(1990), Rawls et al(1982) and 
#' Saxton & Rawls (2006). The model can also return tabulated Ks of Carsel and Parrish (1988)
#' when soil texture is given as a "model" parameters.
#' @inheritParams BEST
#' @param alpha van Genuchten's  water retention parameter. 
#' It is equal to 1/hb
#' @param hb the bubbling capillary pressure parameter of Brooks and Corey. It is equal to 1/apha.
#' @param th33 The water content at field capacity [m3/m3] .i.e metric potential of 33 kPa
#' @param th1500 The water content at wilting point [m3/m3] i.e metric potential of 1500 kPa
#' @param C clay percentage
#' @param f porosity
#' @param para a parameter for each model
#' @param n van Genuchten's n that is realted to pore-size distribution.
#' @param model model type. It can take  "Nasta" (default), "Mishra", "RAWLS", "rawls2006", 
#' "Guarracino". 
#' @references 
#' \itemize{
#' \item{}{Guarracino, L. (2007). Estimation of saturated hydraulic conductivity 
#' Ks from the van Genuchten shape parameter. Water Resour. Res., 43. 
#' doi: doi:10.1029/2006WR005766}
#' \item{}{Paolo Nasta, Jasper A. Vrugt, & Romano, N. (2013). Prediction of the 
#' saturated hydraulic conductivity from Brooks and Corey's 
#' water retention parameters. WATER RESOURCES RESEARCH, 49. doi: doi:10.1002/wrcr.20269}
#' \item{}{Mishra, S., & Parker, J. C. (1990). On the relation between 
#' saturated conductivity and capillary retention characteristics. 
#' Ground Water, 28, 775-777.}
#' \item{}{Rawls, W. J., Brakensiek, D. L., & Saxon, K. E. (1982). 
#' Estimation of soil water properties. Transactions of the ASAE, 25, 1316-1320.}
#' \item{}{Saxton, K. E., & Rawls, W. J. (2006). Soil Water Characteristic Estimates by 
#' Texture and Organic Matter for Hydrologic Solutions. 
#' SOIL SCI. SOC. AM. J., 70. doi: 10.2136/sssaj2005.0117}
#' }
#' @return 
#' \itemize{
#' \item{Ks:} {  saturated hydraulic conductivity}
#' \item{para:} { soil parameters}
#' }
#' @export
#'
#' @examples
#' data=read.csv(system.file("ext","sys","retentionVG.csv",package="vadose"))
#' mod<-vg(data=data,h="h",theta="theta",thr=0.1, ths=0.1, alp=0.1, n=1,Ks="nvr")
#' nvr=ksat(ths=mod$ths,hb=1/mod$alp,n=mod$n-1,model="Nasta")
#' print(nvr)
#' g=ksat(ths=mod$ths,alpha=mod$alp,model="Guarracino")
#' print(g$Ks) 

#' soil=ksat(ths=mod$ths,alpha=mod$alp,model="Clay")
#' soil=ksat(para="soil",model="Clay")

#' mp=ksat(ths=mod$ths,alpha=mod$alp,model="Mishra")
#' print(mp)

#' r=ksat(ths=mod$ths,hb=1/mod$alp,n=mod$n-1,model="RAWLS",thr=mod$thr,f=mod$ths,para=86)
#' print(r)

#' sr=ksat(ths=mod$ths,th33=predict(mod,33),th1500=predict(mod,1500),model="rawls2006")
#' print(sr)
#generic function
ksat<-function(ths=NULL,thr=0,alpha=NULL,hb=NULL,th33=NULL,th1500=NULL,C=NULL,f=NULL,para=NULL,n=NULL,model="nvr") UseMethod ("ksat")

#' @rdname ksat
#' @export
ksat.default<-function(ths=NULL,thr=0,alpha=NULL,hb=NULL,th33=NULL,th1500=NULL,C=NULL,f=NULL,para=NULL,n=NULL,model="nvr"){

#if(!is.null(para)&&para!="soil"){
#Nasta model
if((model=="NVR"||model=="nvr"||model=="nasta"||model=="NASTA"||model=="Nasta")){
if(is.null(ths)||is.null(hb)||is.null(n)){
return(print("Please provide thetaS, air entry pressure hb and 
shaping parameter, n. See Paolo Nasta, Jasper A. Vrugt, and Romano N. 2013. Prediction of the saturated hydraulic conductivity from Brooks and Corey's water retention parameters. Water Resources Research 49."))
}
if(is.null(para)){
para=0.0138
}
p2=n
Ks=9.579*10^5*para*(p2/(p2+2))*(ths/hb^2)

#return(Ks[[1]])
}

##Guarracino [2007] model
if(model=="G"||model=="g"||model=="Guarracino"||model=="guarracino"||model=="GUARRACINO"){
if(is.null(ths)||is.null(alpha)){
return(print("Please provide thetaS,alpha and parameter D 
see Guarracino L. 2007. Estimation of saturated hydraulic conductivity Ks from the van Genuchten shape parameter. Water Resour. Res. 43."))
}
D=para
if(is.null(para)){
D=1.9943
para=D
}

Ks=1.74*10^5*((2-D)/(4-D))*ths*(alpha^2)
#return(Ks[[1]])
}

##Mishra, 1990 model
if(model=="mp"||model=="MP"||model=="Mishra"||model=="mishra"||model=="MISHRA"){
if(is.null(ths)||is.null(alpha)){
return(print("Please provide thetaS,alpha and parameter mp see Mishra S, Parker J C. 1990. On the relation between saturated conductivity and capillary retention characteristics. Ground Water 28:775-777."))
}
mp=para
if(is.null(para)){
mp=39.09
para=mp
}

Ks=((3.89*10^5)/mp)*((ths^2.5)*alpha^2)
#return(Ks[[1]])
}

##Rawls, 1982 model
if(model=="r"||model=="R"||model=="rawls82"||model=="RAWLS"||model=="Rawls"){
if(is.null(f)||is.null(n)||is.null(thr)||is.null(hb)){
return(print("Please provide n, porosity, f, thr and parameter x1; see Rawls W J, Brakenseik D L, and Saxton K E. 1982. Estimation of soil water properties. Trans. ASAE 25:1316-1320."))
}
x1=para
if(is.null(para)){
x1=86
para=x1
}

the=f-thr
Ks=(x1*((the/hb)^2))*(n^2/((n+1)*(n+2)))
#return(Ks[[1]])
}

##Rawls, 1982 model
if(model=="sr"||model=="SR"||model=="rawls2006"){
if(is.null(ths)||is.null(th33)||is.null(th1500)){
return(print("Please provide ths, th33, thr and th1500; see Saxton K E, Rawls W J. 2006. Soil Water Characteristic Estimates by Texture and Organic Matter for Hydrologic Solutions. Soil Sci. Soc. Am. J. 70."))
}

B=(log(1500)-log(33))/(log(th33)-log(th1500))
slope=1/B
para=slope
Ks=1930*(ths-th33)^(3-slope)
#return(Ks[[1]])
}

#Carsel and Parrish 1988
  file=system.file("ext","sys","CarselandParrish.csv",package="vadose")
  
#file="C:/Users/GeoKings/Documents/George Owusu/UG PhD/Data/Carsel and Parrish.csv"
mydata=read.csv(file)
model2=toupper(model)
texture=toupper( as.character(mydata[[1]]))
pattern=match(model2, texture)
if(is.null(para)){
para="NULL"
}
if( !is.na(pattern)||para=="soil"||para=="Soil"){
#texture=c("Sand","Loamy Sand","Sandy Loam","Loam","Silt","Silt Loam","Sandy Clay Loam","Clay Loam","Silty Clay Loam","Sandy Clay","Silty Clay","Clay")
Ks=model
if( !is.na(pattern)){
Ks <- mydata[ which(mydata$Texture==model2), ]$Ks
}

if(!is.null(para)){
if(para=="soil"||para=="Soil"){
if( is.na(pattern)){
return(print("Please state the model value as soil texture"))
break
}
para <- mydata[ which(mydata$Texture==model2), ]
}
}
}

if( !is.na(pattern)){
  factor=list(Ks=Ks,para=para )

}
else
{
  factor=list(Ks=Ks[[1]],para=para )

}

factor$call<-match.call()

class(factor)<-"ksat"
factor
}



#mod<-vg(data=datartc,h="h",theta="theta",thr=0.1, ths=0.1, alp=0.1, n=1,Ks="nvr")
#nvr=ksat(ths=mod$ths,hb=1/mod$alp,n=mod$n-1,model="NVR")
#print(nvr)
#g=ksat(ths=mod$ths,alpha=mod$alp,model="g")
#print(g)

#soil=ksat(ths=mod$ths,alpha=mod$alp,model="Clay")
#soil=ksat(para="soil",model="Clay")

#mp=ksat(ths=mod$ths,alpha=mod$alp,model="mp")
#print(mp)

#r=ksat(ths=mod$ths,hb=1/mod$alp,n=mod$n-1,model="r",thr=mod$thr,f=mod$ths,para=86)
#print(r)

#sr=ksat(ths=mod$ths,th33=predict(mod,33),th1500=predict(mod,1500),model="sr")
#print(sr)
