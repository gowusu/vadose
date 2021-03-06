# R Vadose Package: estimation of soil water characteristics

The vadose package can derive soil water retention  and hydraulic conductivity curves using 
models from the Offin Estimation of Soil Transfer parameters (OFEST) and Beerkan Estimation of Soil Transfer parameters (BEST). 

Specifically, the package can perform the following functions:

+   2 Soil water characteristics Curves functions: OFEST and BEST
+   10 Many infiltration functions: green and Ampt,Philip, SCS (US-Soil Conservation Service), Kostiakov, Horton,Swartzendruber,Modified Kostiakov, Revised Modified Kostiakov, brutsaert and valiantzas, etc
+   15 particle size distribution functions: logistic models, van Genuchten models, Fredlund models, ogarithmic model, Andersson models, and Gompterz, etc
+   van Genuchten water retention and unsaturated hydraulic conductivity model with internal saturated hydraulic conductivity functions 

## Installation
You can install the latest development version from github with,
 <pre><code> if (packageVersion("devtools") < 1.6) {
 
  if (!require(devtools)) {
  install.packages("devtools")
  }
  install.packages("devtools")
  }
  if(!require(nlmrt)){
  install.packages("nlmrt")
  }
if(!require(httr)){
  install.packages("httr")
  }
if(!require(vadose)){
devtools::install_github("gowusu/vadose")
  }

</code></pre>


**Author**: George Owusu <owusugeorge@ug.edu.gh> This is part of the author's PhD program.
