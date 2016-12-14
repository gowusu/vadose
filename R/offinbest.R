#' Horizontal Infiltration data of Offin Basin
#'
#' A dataset containing horizontal infiltration data of Offin Basin.
#' 
#' @name offinbest
#' @docType data
#' @author George Owusu
#' @format A data frame with 302 rows and 9 variables:
#' \describe{
#'   \item{ID}{Unique ID of the data}
#'   \item{TownName}{Town Name}
#'   \item{time}{Recorded cumulative time [seconds]}
#'   \item{I}{Cumulative Infiltration [mm]}
#'   \item{pb}{Bulk density [g/cm3]}
#'   \item{r}{Radius of the cylinder ring}
#'   \item{tho}{initial soil water content [-]}
#'   \item{thr}{residual soil water content [-]}
#'   \item{n}{van Genuchten n parameter}
#' }
#' 
NULL

#' Measured water retention data 
#'
#' A dataset hydraulic heads and water content.
#' 
#' @name htheta
#' @docType data
#' @author George Owusu
#' @format A data frame with 90 rows and 4 variables:
#' \describe{
#'   \item{ID}{Unique ID of the data}
#'   \item{TownName}{Town Name}
#'   \item{hbar}{Metric potentials in bar}
#'   \item{h}{Metric potential in mm}
#'   \item{theta}{gravimetric soil water content [-]}
#' }
#' 
NULL

#' Particle Size Distribution Size of Offin Basin
#'
#' A dataset containing diameter and fractional proportion of soil particle in OFFin Basin
#' 
#' @name offinpsd
#' @docType data
#' @author George Owusu
#' @format A data frame with 311 rows and 4 variables:
#' \describe{
#'   \item{ID}{Uniquet ID of the data}
#'   \item{TownName}{Town Name}
#'   \item{D}{Particle diameter [mm]}
#'   \item{FractionSand}{proportion of particle sizes [0-1]}
#'   \item{fr}{proportion of particle sizes [0-1]}
#' }
#' 
NULL

#' Metadata of Infiltration data of Offin Basin
#'
#' A metadata  of  infiltration data of Offin Basin.
#' 
#' @name offinmeta
#' @docType data
#' @author George Owusu
#' @format A data frame with 15 rows and 18 variables:
#' \describe{
#'   \item{ID}{Unique ID of the data}
#'   \item{TownName}{Town Name}
#'   \item{longitude}{longitudes in degrees}
#'   \item{latitude}{latitudes in degrees}
#'   \item{Time}{Total time of infiltration in minutes }
#'   \item{N}{the Number of water volumes (cylinders)}
#'   \item{q}{the mean infiltration rate in mm/s}
#'   \item{Texture}{USDA soil texture}
#'   \item{tho}{Initial soil water content (cm3/cm3)}
#'   \item{thr}{residual water content(cm3/cm3)}
#'   \item{pb}{bulk density (g/cm3)}
#'   \item{n}{van Genuchten n parameter}
#'   \item{m}{van Genuchten m parameter}
#'   \item{n2}{Particle Size parameter}
#'   \item{b}{Particle Size distribution parameter}

#' }
#' 
NULL

