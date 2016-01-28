library(igraph)



#' AFM network edge class
#' 
#' @name AFMImageNetworkEdge-class
#' @rdname AFMImageNetworkEdge-class
#' @exportClass AFMImageNetworkEdge
AFMImageNetworkEdge<-setClass("AFMImageNetworkEdge",
                   slots = c(from="character",
                             to="character",
                             weight="numeric"
                   ))


#' Constructor method of AFMImageNetworkEdge Class.
#'
#' @rdname AFMImageNetworkEdge-class
#' @export
setMethod(f="initialize",
          signature="AFMImageNetworkEdge", 
          definition= function(.Object,
                               from,
                               to, 
                               weight
                              )  
          { 
            if (!missing(from)) .Object@from<-from
            if (!missing(to)) .Object@to<-to
            if (!missing(weight)) .Object@weight<-weight
            return(.Object)
          })

#' Wrapper function AFMImageNetworkEdge
#'
#' @rdname AFMImageNetworkEdge-class
#' @export
AFMImageNetworkEdge <- function(from,
                                to, 
                                weight) {
  return(new("AFMImage", from, to,  weight))
}


HASHSIZE<-512*512

#' @author M.Beauvais
#' @export
getVertexId<-function(AFMImage,x,y) {
  if ((x<1)||(x>AFMImage@samplesperline)||
      (y<1)||(y>AFMImage@lines)) return(-1)
  #print(paste("getVertexId",x,y,as.numeric(x+HASHSIZE*y)))
  #return(as.numeric(x+AFMImage@samplesperline*y))
  return(as.numeric(x+HASHSIZE*y))
  
}

#' @author M.Beauvais
#' @export
getCoordinatesFromVertexId<-function(AFMImage, vId) {
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(c(x,y))
}

#' @author M.Beauvais
#' @export
getNetworkGridLayout<-function(AFMImage, vId) {
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(data.table(x=x,y=y))
}

#' @author M.Beauvais
#' @export
existsEdge<-function(AFMImage, vertexId) {
  # print(vertexId)
  if ((vertexId<1)||(vertexId>(AFMImage@samplesperline+HASHSIZE*(AFMImage@lines-1)))) {
    # print("return FALSE")
    return(FALSE)
  }
  # print(vertexId)
  
  
  coordinates<-getCoordinatesFromVertexId(AFMImage, vertexId)
  # print(coordinnates)
  id<-coordinates[1]+AFMImage@samplesperline*coordinates[2]
  # print(id)
  if (AFMImage@data$h[id]>0) {
    # print("return TRUE")
    return(TRUE)
  }
  # print("return FALSE")
  return(FALSE)
}

#' @author M.Beauvais
#' @export
getSurroundingVertexesList<-function(AFMImage,x,y) {
#   print(x)
#   print(y)
  horizontalWeight<-AFMImage@hscansize/AFMImage@samplesperline
  verticalWeight<-AFMImage@vscansize/AFMImage@lines
  diagWeight<-sqrt((AFMImage@vscansize/AFMImage@lines)^2+(AFMImage@hscansize/AFMImage@samplesperline)^2)
  
  currentVertexId<-getVertexId(AFMImage,x,y)
  vList=data.table()
  #x+1 y
  nearVertexId<-getVertexId(AFMImage,x+1,y) 
  # print(nearVertexId)
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(horizontalWeight)))
  #x+1 y+1
  nearVertexId<-getVertexId(AFMImage,x+1,y+1) 
  #print(existsEdge(AFMImage, nearVertexId))
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x y+1
  nearVertexId<-getVertexId(AFMImage,x,y+1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(verticalWeight)))
  
  #x-1 y+1
  nearVertexId<-getVertexId(AFMImage,x-1,y+1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x-1 y
  nearVertexId<-getVertexId(AFMImage,x-1,y) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(horizontalWeight)))
  
  #x-1 y-1
  nearVertexId<-getVertexId(AFMImage,x-1,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  
  #x y-1
  nearVertexId<-getVertexId(AFMImage,x,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(verticalWeight)))
  
  #x+1 y-1
  nearVertexId<-getVertexId(AFMImage,x+1,y-1) 
  if (existsEdge(AFMImage, nearVertexId)) vList<-rbind(vList, data.table(from=as.character(currentVertexId), to=as.character(nearVertexId), weight=as.numeric(diagWeight)))
  return(vList)
}

#' @author M.Beauvais
#' @export
isAdjacentToBetterVertex<-function(AFMImage,x,y) {
  #   print(x)
  #   print(y)

  currentVertexId<-getVertexId(AFMImage,x,y) 
  currentH<-AFMImage@data$h[currentVertexId]
  
  if(currentH<=0) return(FALSE)
  
  #x+1 y
  nearVertexId<-getVertexId(AFMImage,x+1,y) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x+1 y+1
  nearVertexId<-getVertexId(AFMImage,x+1,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x y+1
  nearVertexId<-getVertexId(AFMImage,x,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y+1
  nearVertexId<-getVertexId(AFMImage,x-1,y+1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y
  nearVertexId<-getVertexId(AFMImage,x-1,y) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x-1 y-1
  nearVertexId<-getVertexId(AFMImage,x-1,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x y-1
  nearVertexId<-getVertexId(AFMImage,x,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  #x+1 y-1
  nearVertexId<-getVertexId(AFMImage,x+1,y-1) 
  if ((nearVertexId>0)&(currentH<=AFMImage@data$h[nearVertexId])) return(TRUE)
  
  return(FALSE)
}
