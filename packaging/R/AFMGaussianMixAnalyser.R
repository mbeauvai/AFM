require("data.table")
require("mixtools")

# normality tests
require(gridExtra)
require(ggplot2)


#if(getRversion() >= "3.1.0") utils::suppressForeignCheck(c("r", "roughness","x","predict.gstat"))


#' @title AFM image Gaussian Mix analysis class
#' 
#' @description \code{AFMImageGaussianMixAnalysis} handles an \code{\link{AFMImage}} Gaussian mix of heights analysis 
#' 
#' @slot minGaussianMix the minimum number of components to calculate
#' @slot maxGaussianMix the maximum number of components to calculate
#' @slot epsilonGaussianMix the convergence criterion
#' @slot gaussianMix a data.table to store the calculated Gaussian mixes
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImageGaussianMixAnalysis-class
#' @rdname AFMImageGaussianMixAnalysis-class
#' @author M.Beauvais
AFMImageGaussianMixAnalysis<-setClass("AFMImageGaussianMixAnalysis",
                              slots = c(
                                minGaussianMix="numeric",
                                maxGaussianMix="numeric",
                                epsilonGaussianMix="numeric",
                                gaussianMix="array",
                                updateProgress="function"),
                              validity = function(object) { 
                                return(TRUE)
                              }
)

#' Constructor method of AFMImageGaussianMixAnalysis Class.
#' 
#' @param .Object an AFMImageGaussianMixAnalysis object
#' @rdname AFMImageGaussianMixAnalysis-class
#' @export
setMethod("initialize",
          "AFMImageGaussianMixAnalysis",
          function(.Object) {
            .Object@minGaussianMix<-2
            .Object@maxGaussianMix<-2
            .Object@epsilonGaussianMix<-1e-4
            .Object@gaussianMix<-array()
            validObject(.Object) ## valide l'objet
            return(.Object)
          })

#' Wrapper function AFMImageGaussianMixAnalysis
#'
#' @rdname AFMImageGaussianMixAnalysis-class
#' @export
AFMImageGaussianMixAnalysis <- function() {
  return(new("AFMImageGaussianMixAnalysis"))
}

#' Method \code{GaussianMix} returns a data.table of Gaussian mixes
#' @name AFMImageGaussianMixAnalysis-class
#' @rdname AFMImageGaussianMixAnalysis-class
setGeneric("gaussianMix",function(object){standardGeneric("gaussianMix")})
setGeneric(name= "gaussianMix<-", 
           def= function(AFMImageGaussianMixAnalysis, value) {
             return(standardGeneric("gaussianMix<-"))
           })


#' @rdname AFMImageGaussianMixAnalysis-class
#' @aliases gaussianMix
#' @param object a \code{\link{AFMImageGaussianMixAnalysis}}
setMethod("gaussianMix",signature=signature(object='AFMImageGaussianMixAnalysis'),
          function(object) {
            return(object@gaussianMix)
          }
)
setReplaceMethod(f="gaussianMix",
                 signature(AFMImageGaussianMixAnalysis = "AFMImageGaussianMixAnalysis", value = "array"),
                 definition= function(AFMImageGaussianMixAnalysis, value) {
                   AFMImageGaussianMixAnalysis@gaussianMix <- value
                   return(AFMImageGaussianMixAnalysis)
                 })

#' Method \code{minGaussianMix} returns a data.table of Gaussian mixes
#' @name AFMImageGaussianMixAnalysis-class
#' @rdname AFMImageGaussianMixAnalysis-class
setGeneric("minGaussianMix",function(object){standardGeneric("minGaussianMix")})
setGeneric(name= "minGaussianMix<-", 
           def= function(AFMImageGaussianMixAnalysis, value) {
             return(standardGeneric("minGaussianMix<-"))
           })


#' @rdname AFMImageGaussianMixAnalysis-class
#' @aliases minGaussianMix
setMethod("minGaussianMix",signature=signature(object='AFMImageGaussianMixAnalysis'),
          function(object) {
            return(object@minGaussianMix)
          }
)
setReplaceMethod(f="minGaussianMix",
                 signature(AFMImageGaussianMixAnalysis = "AFMImageGaussianMixAnalysis", value = "numeric"),
                 definition= function(AFMImageGaussianMixAnalysis, value) {
                   AFMImageGaussianMixAnalysis@minGaussianMix <- value
                   return(AFMImageGaussianMixAnalysis)
                 })



#' Method \code{maxGaussianMix} returns a data.table of Gaussian mixes
#' @name AFMImageGaussianMixAnalysis-class
#' @rdname AFMImageGaussianMixAnalysis-class
setGeneric("maxGaussianMix",function(object){standardGeneric("maxGaussianMix")})
setGeneric(name= "maxGaussianMix<-", 
           def= function(AFMImageGaussianMixAnalysis, value) {
             return(standardGeneric("maxGaussianMix<-"))
           })


#' @rdname AFMImageGaussianMixAnalysis-class
#' @aliases maxGaussianMix
setMethod("maxGaussianMix",signature=signature(object='AFMImageGaussianMixAnalysis'),
          function(object) {
            return(object@maxGaussianMix)
          }
)
setReplaceMethod(f="maxGaussianMix",
                 signature(AFMImageGaussianMixAnalysis = "AFMImageGaussianMixAnalysis", value = "numeric"),
                 definition= function(AFMImageGaussianMixAnalysis, value) {
                   AFMImageGaussianMixAnalysis@maxGaussianMix <- value
                   return(AFMImageGaussianMixAnalysis)
                 })



#' Method \code{epsilonGaussianMix} returns a data.table of Gaussian mixes
#' @name AFMImageGaussianMixAnalysis-class
#' @rdname AFMImageGaussianMixAnalysis-class
setGeneric("epsilonGaussianMix",function(object){standardGeneric("epsilonGaussianMix")})
setGeneric(name= "epsilonGaussianMix<-", 
           def= function(AFMImageGaussianMixAnalysis, value) {
             return(standardGeneric("epsilonGaussianMix<-"))
           })


#' @rdname AFMImageGaussianMixAnalysis-class
#' @aliases epsilonGaussianMix
setMethod("epsilonGaussianMix",signature=signature(object='AFMImageGaussianMixAnalysis'),
          function(object) {
            return(object@epsilonGaussianMix)
          }
)
setReplaceMethod(f="epsilonGaussianMix",
                 signature(AFMImageGaussianMixAnalysis = "AFMImageGaussianMixAnalysis", value = "numeric"),
                 definition= function(AFMImageGaussianMixAnalysis, value) {
                   AFMImageGaussianMixAnalysis@epsilonGaussianMix <- value
                   return(AFMImageGaussianMixAnalysis)
                 })





#' Perform  the calculation for the Gaussian mixes
#' 
#' \code{\link{performGaussianMixCalculation}} perform all the calculation for PSD exploitation
#' @param AFMImageGaussianMixAnalysis an \code{\link{AFMImageGaussianMixAnalysis}} to manage and store the results of PSD analysis
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
#' @examples
#' library(AFM)
#' 
#' data(AFMImageCollagenNetwork)
#' 
#' AFMImage<-AFMImageCollagenNetwork
#' AFMImage@fullfilename<-"/Users/one/AFMImageCollagenNetwork.txt"
#' gMixAnalysis<-AFMImageGaussianMixAnalysis()
#' # Create a closure to update progress
#' gMixAnalysis@updateProgress<- function(value = NULL, detail = NULL, message = NULL) {
#'   if (exists("progressGaussianMix")){
#'     if (!is.null(message)) {
#'       progressGaussianMix$set(message = message, value = 0)
#'     }else{
#'       progressGaussianMix$set(value = value, detail = detail)
#'     }
#'   }
#' }
#' gMixAnalysis<-performGaussianMixCalculation(AFMImageGaussianMixAnalysis= gMixAnalysis, AFMImage)
#' print("done performGaussianMixCalculation")
performGaussianMixCalculation<-function(AFMImageGaussianMixAnalysis, AFMImage) {
  if (is.function(AFMImageGaussianMixAnalysis@updateProgress)&&
      !is.null(AFMImageGaussianMixAnalysis@updateProgress())) {
    AFMImageGaussianMixAnalysis@updateProgress(message="Calculating Gaussian Mixes")
  }
  
  #data(AFMImageCollagenNetwork)
  #AFMImage<-AFMImageCollagenNetwork
  
  # parameters
  min<-AFMImageGaussianMixAnalysis@minGaussianMix
  max<-AFMImageGaussianMixAnalysis@maxGaussianMix
  mepsilon<-AFMImageGaussianMixAnalysis@epsilonGaussianMix
  
  gaussianMixList = array(list(), max)

  min_height<- 0
  max_height<- 3000
  heights<-AFMImage@data$h
  heights<-heights[heights<(max_height/10)]

  # allH<-data.table(h=heights)
  # g<-ggplot(allH, aes(h)) + geom_histogram(binwidth = 0.1)
  # print(g)

  mixtureCounter<-0
  mixtureNumberOfComponents<-min
  for(mixtureNumberOfComponents in seq(min,max)){
    if (is.function(AFMImageGaussianMixAnalysis@updateProgress)&&
        !is.null(AFMImageGaussianMixAnalysis@updateProgress())) {
      mixtureCounter<-mixtureCounter+1
      AFMImageGaussianMixAnalysis@updateProgress(message=paste("Calculating Gaussian Mixes", mixtureCounter ,"/",(as.numeric(max)-as.numeric(min)+1)) , detail = paste0((as.numeric(max)-as.numeric(min)+1),"/",(as.numeric(max))))
      AFMImageGaussianMixAnalysis@updateProgress(value= 6)
      
    }
    

      heights.k<- mixtools::normalmixEM(heights,
                                        k=mixtureNumberOfComponents,
                                        arbmean = TRUE,
                                        ECM=TRUE,
                                        verb=TRUE,
                                        maxit=10000,
                                        epsilon=mepsilon)
      #heights.k
      gaussianMixList[[mixtureNumberOfComponents]]<-heights.k
  }
  
  AFMImageGaussianMixAnalysis@gaussianMix<-gaussianMixList
  
  #print(gaussianMixList)
  return(AFMImageGaussianMixAnalysis)
}

getGaussianMix<-function(exportDirectory, sampleName) {
  exportCsvFilename<-paste(sampleName,"-gaussian-mix.png", sep="")
  exportCsvFullFilename<-paste(exportDirectory, exportCsvFilename, sep="/")
  return(exportCsvFullFilename)
}

#' pnormmix distribution of a mixture of normals
#' 
#' @param q a vector of quantiles
#' @param mixture a gaussian mixture
#' @export
pnormmix <- function(q,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(q,component) {
    lambda[component]*pnorm(q,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,q=q)
  return(rowSums(pnorms))
}

#' dnormalmix density of a mixture of normals
#' 
#' @param x a vector of quantiles
#' @param mixture a gaussian mixture
#' @param log perform a log transsformation of the result
#' @export
dnormalmix <- function(x,mixture,log=FALSE) {
  lambda <- mixture$lambda
  k <- length(lambda)
  like.component <- function(x,component) {
    lambda[component]*dnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  likes <- sapply(1:k,like.component,x=x)
  d <- rowSums(likes)
  if (log) {
    d <- log(d)
  }
  return(d)
}

#' loglike sum of density of a mixture of normals
#' 
#' @param x a vector of quantiles
#' @param mixture a gaussian mixture
#' @export
loglike.normalmix <- function(x,mixture) {
  loglike <- dnormalmix(x,mixture,log=TRUE)
  return(sum(loglike))
}



