require(igraph)

HASHSIZE<-512*512


setOldClass("igraph")

#' AFM image networks analysis class
#' 
#' A S4 class to handle the networks calculation 
#'
#' @slot binaryAFMImage the AFMImage after transformation before analysis
#' @slot binaryAFMImageWithCircles the AFMImage after transformation with the spotted circles
#' @slot circlesTable a data.table of identified circles
#' @slot edgesTable  a data.table of edges
#' @slot fusionedNodesCorrespondance  a data.table of corresponsdance between intial node and fusioned node
#' @slot fusionedNodesEdgesTable a data.table of nodes fusioned because of intersecting
#' @slot isolatedNodesTable a data.table of isolated nodes
#' @slot heightNetworksslider used multiplier of heights to facilitate analysis
#' @slot filterNetworkssliderMin used filter minimum value to facilitate analysis
#' @slot filterNetworkssliderMax used filter maximum value to facilitate analysis
#' @slot originalGraph a list of \code{\link{igraph}}
#' @slot skeletonGraph a list of \code{\link{igraph}}
#' @slot shortestPaths a data.table of shortest paths
#' @slot libVersion version of the AFM library used to perform the analysis
#' @slot updateProgress a function to update a graphical user interface
#' @name AFMImageNetworksAnalysis-class
#' @rdname AFMImageNetworksAnalysis-class
#' @exportClass AFMImageNetworksAnalysis
#' @author M.Beauvais
AFMImageNetworksAnalysis<-setClass("AFMImageNetworksAnalysis",
                                   slots = c(
                                     binaryAFMImage="AFMImage",
                                     binaryAFMImageWithCircles="AFMImage",
                                     circlesTable="data.table",
                                     edgesTable="data.table",
                                     fusionedNodesCorrespondance="data.table",
                                     fusionedNodesEdgesTable="data.table",
                                     isolatedNodesList="numeric",
                                     heightNetworksslider="numeric",
                                     filterNetworkssliderMin="numeric",
                                     filterNetworkssliderMax="numeric",
                                     originalGraph="igraph", 
                                     skeletonGraph="igraph",
                                     shortestPaths="data.table",
                                     libVersion="character",
                                     updateProgress="function"))

#' Constructor method of AFMImageNetworksAnalysis Class.
#'
#' @param .Object an AFMImageNetworksAnalysis Class
#' @param binaryAFMImage the AFMImage after transformation before analysis
#' @param binaryAFMImageWithCircles the AFMImage after transformation with the spotted circles
#' @param circlesTable a data.table of identified circles
#' @param edgesTable  a data.table of edges
#' @param fusionedNodesCorrespondance  a data.table of correspon
#' @param fusionedNodesEdgesTable a data.table of corresponsdance between intial node and fusioned node
#' @param isolatedNodesList a data.table of isolated nodes
#' @param heightNetworksslider used multiplier of heights to facilitate analysis
#' @param filterNetworkssliderMin used filter minimum value to facilitate analysis
#' @param filterNetworkssliderMax used filter maximum value to facilitate analysis
#' @param originalGraph a list of \code{\link{igraph}}
#' @param skeletonGraph a list of \code{\link{igraph}}
#' @param shortestPaths a data.table of shortest path
#' @param libVersion version of the AFM library used to perform the analysis
#' @rdname AFMImageNetworksAnalysis-class
#' @export
setMethod("initialize", "AFMImageNetworksAnalysis", function(.Object, 
                                                             binaryAFMImage,
                                                             binaryAFMImageWithCircles,
                                                             circlesTable,
                                                             edgesTable,
                                                             fusionedNodesCorrespondance,
                                                             fusionedNodesEdgesTable,
                                                             isolatedNodesList,
                                                             heightNetworksslider,
                                                             filterNetworkssliderMin,
                                                             filterNetworkssliderMax,                                                             
                                                             originalGraph, 
                                                             skeletonGraph,
                                                             shortestPaths,
                                                             libVersion)  
{
  if(!missing(binaryAFMImage)) .Object@binaryAFMImage<-binaryAFMImage
  if(!missing(binaryAFMImageWithCircles)) .Object@binaryAFMImageWithCircles<-binaryAFMImageWithCircles
  if(!missing(circlesTable)) .Object@circlesTable<-circlesTable
  if(!missing(edgesTable)) .Object@edgesTable<-edgesTable
  if(!missing(fusionedNodesCorrespondance)) .Object@fusionedNodesCorrespondance<-fusionedNodesCorrespondance
  if(!missing(fusionedNodesEdgesTable)) .Object@fusionedNodesEdgesTable<-fusionedNodesEdgesTable
  if(!missing(isolatedNodesList)) .Object@isolatedNodesList<-isolatedNodesList
  if(!missing(originalGraph)) .Object@originalGraph<-originalGraph
  if(!missing(skeletonGraph)) .Object@skeletonGraph<-skeletonGraph
  if(!missing(shortestPaths)) .Object@shortestPaths<-shortestPaths
  if(!missing(heightNetworksslider)) .Object@heightNetworksslider<-heightNetworksslider
  if(!missing(filterNetworkssliderMin)) .Object@filterNetworkssliderMin<-filterNetworkssliderMin
  if(!missing(filterNetworkssliderMax)) .Object@filterNetworkssliderMax<-filterNetworkssliderMax
  if(!missing(libVersion)) .Object@libVersion<-libVersion
  #validObject(.Object)      
  return(.Object)
})


#' Wrapper function AFMImageNetworksAnalysis
#'
#' @rdname AFMImageNetworksAnalysis-class
#' @export
AFMImageNetworksAnalysis <- function() {
  return(new("AFMImageNetworksAnalysis"))
}

#' Multiply, filter the heights and make a binary AFMImage from the transformed AFMImage
#'
#' \code{transformAFMImageForNetworkAnalysis} update  \code{\link{AFMImageNetworksAnalysis}} making a binary AFMImage
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis n \code{\link{AFMImageNetworksAnalysis}} to store the results of the transformation
#' 
#' @name transformAFMImageForNetworkAnalysis
#' @rdname transformAFMImageForNetworkAnalysis-methods
#' @exportMethod transformAFMImageForNetworkAnalysis
#' @author M.Beauvais
setGeneric(name= "transformAFMImageForNetworkAnalysis", 
           def= function(AFMImageNetworksAnalysis, AFMImage) {
             return(standardGeneric("transformAFMImageForNetworkAnalysis"))
           })

#' @rdname transformAFMImageForNetworkAnalysis-methods
#' @aliases transformAFMImageForNetworkAnalysis,AFMImage-method
setMethod(f="transformAFMImageForNetworkAnalysis", "AFMImageNetworksAnalysis",
          definition= function(AFMImageNetworksAnalysis, AFMImage) {
            newAFMImage<-multiplyHeightsAFMImage(AFMImage, multiplier=AFMImageNetworksAnalysis@heightNetworksslider)
            newAFMImage<-filterAFMImage(newAFMImage,
                                        Min=AFMImageNetworksAnalysis@filterNetworkssliderMin,
                                        Max=AFMImageNetworksAnalysis@filterNetworkssliderMax)
            newAFMImage<-makeBinaryAFMImage(newAFMImage)
            AFMImageNetworksAnalysis@binaryAFMImage<-copy(newAFMImage)
            return(AFMImageNetworksAnalysis)
          })

#' Calculate networks on the surface
#'
#' \code{calculateNetworks} update  \code{\link{AFMImageNetworksAnalysis}}
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis n \code{\link{AFMImageNetworksAnalysis}} to store the results of networks analysis
#' 
#' @name calculateNetworks
#' @rdname calculateNetworks-methods
#' @exportMethod calculateNetworks
#' @author M.Beauvais
setGeneric(name= "calculateNetworks", 
           def= function(AFMImageNetworksAnalysis, AFMImage) {
             return(standardGeneric("calculateNetworks"))
           })

#' @rdname calculateNetworks-methods
#' @aliases calculateNetworks,AFMImage-method
setMethod(f="calculateNetworks", "AFMImageNetworksAnalysis",
          definition= function(AFMImageNetworksAnalysis, AFMImage) {
            
            counter<-0
            totalLength<-2
            if (!is.null(AFMImageNetworksAnalysis@updateProgress)&&
                is.function(AFMImageNetworksAnalysis@updateProgress)&&
                !is.null(AFMImageNetworksAnalysis@updateProgress())) {
              text <- paste0("Creating networks")
              AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
              
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0("Creating networks", round(counter, 2),"/",totalLength)
              AFMImageNetworksAnalysis@updateProgress(value= value, detail = text)
              print("update")
            }
            
            AFMImageNetworksAnalysis@originalGraph<-calculateIgraph(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis, AFMImage = AFMImage)
            
            if (!is.null(AFMImageNetworksAnalysis@updateProgress)&&
                is.function(AFMImageNetworksAnalysis@updateProgress)&&
                !is.null(AFMImageNetworksAnalysis@updateProgress())) {
              text <- paste0("Creating networks skeleton")
              AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
              
              counter<-counter+1
              value<-counter / totalLength
              text <- paste0("Creating networks", round(counter, 2),"/",totalLength)
              AFMImageNetworksAnalysis@updateProgress(value= value, detail = text)
              print("update")
            }
            
            AFMImageNetworksAnalysis<-calculateNetworkSkeleton(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis, AFMImage = AFMImage)
            
            return(AFMImageNetworksAnalysis)
          })


#' Get vertex id from x,y coordinates
#'
#' \code{getVertexId} return the vertexId
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
#' @author M.Beauvais
#' @export
getVertexId<-function(AFMImage,x,y) {
  if ((x<0)||(x>AFMImage@samplesperline)||
      (y<0)||(y>AFMImage@lines)) return(-1)
  #print(paste("getVertexId",x,y,as.numeric(x+HASHSIZE*y)))
  #return(as.numeric(x+AFMImage@samplesperline*y))
  return(as.numeric(x+HASHSIZE*y))
  
}

#' Get x,y coordinates from vertex id
#'
#' \code{getCoordinatesFromVertexId} return a list x,y coordinates
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vId the vertex id
#' @author M.Beauvais
#' @export
getCoordinatesFromVertexId<-function(AFMImage, vId) {
  # vertexId<-as.numeric(vId)
  # y<-floor(vertexId/HASHSIZE)
  # x<-vertexId-y*HASHSIZE
  # return(c(x,y))
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(data.table(vId=vId, coords.x1=x,coords.x2=y))
}

#' #' @export
#' getCoordinatesFromVertexId2<-function(AFMImage, vId) {
#'   vertexId<-as.numeric(vId)
#'   y<-floor(vertexId/HASHSIZE)
#'   x<-vertexId-y*HASHSIZE
#'   return(data.table(vId=vId, coords.x1=x,coords.x2=y))
#' }

#' Get getNetworkGridLayout
#'
#' \code{getNetworkGridLayout} return a list x,y coordinates
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vId the vertex id
#' @author M.Beauvais
#' @export
getNetworkGridLayout<-function(AFMImage, vId) {
  vertexId<-as.numeric(vId)
  y<-floor(vertexId/HASHSIZE)
  x<-vertexId-y*HASHSIZE
  return(data.table(x=x,y=y))
}

#' Does an edge exist ?
#'
#' \code{existsEdge} return TRUE if an edge exists for this vertex id
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param vertexId the vertex id
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

#' Get surrounding vertices from x,y coordinates
#'
#' \code{getSurroundingVertexesList} return the vertexId
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
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

#' isAdjacentToBetterVertex
#'
#' \code{isAdjacentToBetterVertex} return TRUE if vertex is adjacent to a better vertex
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param x coordinates in x axis
#' @param y coordinates in y axis
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

#' gridIgraphPlot
#'
#' \code{gridIgraphPlot} return TRUE if vertex is adjacent to a better vertex
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param g the networks
#' @author M.Beauvais
#' @export
gridIgraphPlot<-function(AFMImage, g){
  # define the layout matrix
  coordinatesVector<-getNetworkGridLayout(AFMImage, V(g)$name)
  #coordinatesVector
  
  l<-matrix(coordinatesVector$x ,byrow = TRUE)
  l<-cbind(l, coordinatesVector$y)
  #l
  
  # plot(all, layout=All_layout, vertex.size=2, vertex.label=V(All)$name,
  #      vertex.color="green", vertex.frame.color="red", edge.color="grey",  
  #      edge.arrow.size=0.01, rescale=TRUE,vertex.label=NA, vertex.label.dist=0.0,
  #      vertex.label.cex=0.5, add=FALSE,   vertex.label.font=.001)
  plot(g, layout=l, 
       vertex.shape="circle", vertex.size=2, vertex.label=NA, vertex.color="red", vertex.frame.color="red",
       edge.color="grey"
  )
  
}

#' Calculate iGraph from AFMImage
#'
#' \code{calculateIgraph} return 
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
calculateIgraph<-function(AFMImage, AFMImageNetworksAnalysis) {
  if (missing(AFMImageNetworksAnalysis)) {
    AFMImageNetworksAnalysis<-NULL
  }
  graphicalUpdate<-FALSE
  graphicalCounter<-0
  
  if (!is.null(AFMImageNetworksAnalysis)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress)&&
      is.function(AFMImageNetworksAnalysis@updateProgress)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress())) {
    graphicalUpdate<-TRUE
    totalLength<-AFMImage@samplesperline*(AFMImage@lines-1)
  }
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="1/2 - Generating edges list", value=0)
  }
  print(paste("Generating edge list"))
  
  counter<-1
  #edgeList=data.table()  
  edgeList <- vector("list", AFMImage@samplesperline*AFMImage@lines+1)
  
  for (x in seq(1: AFMImage@samplesperline)) {
    for (y in seq(1: (AFMImage@lines-1))) {
      currentVertexId<-getVertexId(AFMImage,x,y)
      if (existsEdge(AFMImage, currentVertexId)) {
        #edgeList<-rbind(edgeList, getSurroundingVertexesList(AFMImage,x,y))
        edgeList[[counter]] <- getSurroundingVertexesList(AFMImage,x,y)
        counter<-counter+1
      }
      if (graphicalUpdate) {
        graphicalCounter<-graphicalCounter+1
        if (graphicalCounter/100==floor(graphicalCounter/100)) {
          value<-graphicalCounter / totalLength
          text <- paste0(round(graphicalCounter, 2),"/",totalLength)
          AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
        }
      }
    }
  }
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="2/2 - Generating network", value=0)
  }
  
  newEdgeList<-rbindlist(edgeList)
  el=as.matrix(newEdgeList)
  print(paste("Creating graph"))
  g<-graph_from_edgelist(el[,1:2], directed=FALSE)
  print(paste("Created",counter,"vertices"))
  AFMImageNetworksAnalysis@originalGraph<-g
  return(g)
}

#' getListOfDiameters
#'
#' \code{getListOfDiameters} return 
#' 
#' @param g list of igraph networks
#' @author M.Beauvais
#' @export
getListOfDiameters<-function(g) {
  LIST_OF_DIAMETERS = c()
  listOfGraph=decompose(g)
  for(g in listOfGraph){
    LIST_OF_DIAMETERS=c(LIST_OF_DIAMETERS, diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL))
  }
  return(LIST_OF_DIAMETERS)  
}

#' canBeRemoved
#'
#' \code{canBeRemoved} return 
#' 
#' @param vertexId a vertex id
#' @param g a igraph
#' @param allVertices list of all vertices
#' @param DEGREE_LIMIT_FOR_CANDIDATE_VERTICE degree
#' 
#' @author M.Beauvais
#' @export
canBeRemoved<-function(vertexId, g, allVertices, DEGREE_LIMIT_FOR_CANDIDATE_VERTICE) {
  avList<-adjacent_vertices(g, v=c(vertexId), mode = c("all"))
  avListNew<-unique(avList[[vertexId]]$name)
  found<-NULL
  if (nrow(allVertices[, c("found"):=vertexId %in% avListNew & degree<(DEGREE_LIMIT_FOR_CANDIDATE_VERTICE+1)][found==TRUE])>0) {
    return(FALSE)
  }else{
    return(TRUE)
  }
  
}

#' calculateNetworkSkeleton
#'
#' \code{calculateNetworkSkeleton} return 
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
calculateNetworkSkeleton<-function(AFMImage, AFMImageNetworksAnalysis) {
  if (missing(AFMImageNetworksAnalysis)) {
    AFMImageNetworksAnalysis<-NULL
    return(new("list"))
  }
  
  g<-AFMImageNetworksAnalysis@originalGraph
  
  graphicalUpdate<-FALSE
  graphicalCounter<-0
  
  if (!is.null(AFMImageNetworksAnalysis)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress)&&
      is.function(AFMImageNetworksAnalysis@updateProgress)&&
      !is.null(AFMImageNetworksAnalysis@updateProgress())) {
    graphicalUpdate<-TRUE
    totalLength<-length(V(g))
    
  }
  
  
  DEGREE_LIMIT_FOR_CANDIDATE_VERTICE=4
  NUMBER_OF_NETWORKS = length(decompose(g))
  LIST_OF_DIAMETERS<-getListOfDiameters(g)
  print(LIST_OF_DIAMETERS)  
  
  
  #   distance_table(g, directed = FALSE)
  #   coreness(g)
  
  
  verticesThatCantBeRemovedList=c()
  print(paste("starting with ", length(V(g)), " vertices"))
  
  if (graphicalUpdate) {
    AFMImageNetworksAnalysis@updateProgress(message="1/1 - removing vertices and edges", value=0)
  }
  
  continueExploration<-TRUE
  while(continueExploration) {
    
    edgeList<-V(g)$name
    
    uniqueVerticesList<-unique(edgeList)
    uniqueVerticesList
    # degree de chaque noeud
    edgeDegreeList<-degree(g, v=uniqueVerticesList, mode = c("all"), loops = FALSE, normalized = FALSE)
    edgeDegreeList
    
    # liste ordonn'e9e croissante des noeuds en fonction du degree
    
    allVertices<-data.table(vertexId=uniqueVerticesList, degree=edgeDegreeList)
    # get-list of adjacent vertices with degree > 2 (can't remove if degree < 2)
    allVertices<-allVertices[order(degree)]
    listOfCandidateVertices<-allVertices[degree>DEGREE_LIMIT_FOR_CANDIDATE_VERTICE]
    
    listOfCandidateVertices<-listOfCandidateVertices[!listOfCandidateVertices$vertexId %in% verticesThatCantBeRemovedList]
    
    continueExploration<-FALSE
    if (nrow(listOfCandidateVertices)>0) {
      
      #             res<-sapply(listOfCandidateVertices$vertexId, canBeRemoved, g=g, allVertices=allVertices, simplify=F)
      #             vMatrix<-as.matrix(res, ncol=2)
      #             
      #             verticesToBeRemoved<-data.table(vertexId= rownames(vMatrix), toBeRemoved= vMatrix[,1])[toBeRemoved==TRUE]$vertexId
      #             print(paste("to be removed",verticesToBeRemoved))
      #             
      #             if (length(verticesToBeRemoved)>0) {
      #               g<-delete_vertices(g, c(verticesToBeRemoved))
      #               #continueExploration<-TRUE
      #               continueExploration<-continueExploration+1
      #             }
      #       
      for (vi in seq(1:nrow(listOfCandidateVertices))){
        onevertexId=listOfCandidateVertices$vertexId[vi]
        if (canBeRemoved(onevertexId, g=g, allVertices=allVertices, DEGREE_LIMIT_FOR_CANDIDATE_VERTICE=DEGREE_LIMIT_FOR_CANDIDATE_VERTICE)) {
          vId<-listOfCandidateVertices$vertexId[vi]
          
          # store the list of adjacent vertices of the node before deleting it
          avList<-unique(adjacent_vertices(g, v=c(vId), mode = c("all"))[[vId]]$name)
          
          
          g<-delete_vertices(g, listOfCandidateVertices$vertexId[vi])
          continueExploration<-TRUE
          
          NEW_LIST_OF_DIAMETERS=getListOfDiameters(g)
          #print(NEW_LIST_OF_DIAMETERS)  
          
          # did the vertex removal split the network or diminish the diameter
          if ((length(decompose(g))>NUMBER_OF_NETWORKS)||(!identical(LIST_OF_DIAMETERS,NEW_LIST_OF_DIAMETERS))) {
            print (paste("should not have removed", vId))
            verticesThatCantBeRemovedList=c(verticesThatCantBeRemovedList, listOfCandidateVertices$vertexId[vi])
            
            g<-g+vertices(as.numeric(vId))
            
            listOfEdges=c()
            for(j in seq(1,length(avList))) {
              listOfEdges=c(listOfEdges, vId, avList[j], avList[j],vId)
            }
            g<-g+edges(listOfEdges)
          }else{
            print("61")
            NEW_LIST_OF_DIAMETERS=getListOfDiameters(g)
            if ((!identical(LIST_OF_DIAMETERS,NEW_LIST_OF_DIAMETERS))) {
              print (paste("should not have removed", vId))
              verticesThatCantBeRemovedList=c(verticesThatCantBeRemovedList, listOfCandidateVertices$vertexId[vi])
              
              g<-g+vertices(as.numeric(vId))
              
              listOfEdges=c()
              for(j in seq(1,length(avList))) {
                listOfEdges=c(listOfEdges, vId, avList[j], avList[j],vId)
              }
              
              g<-g+edges(listOfEdges)
            }
            break
          }
        }
        
      }
      if (graphicalUpdate) {
        graphicalCounter<-graphicalCounter+1
        value<-graphicalCounter / totalLength
        text <- paste0(round(graphicalCounter, 2),"/",totalLength)
        AFMImageNetworksAnalysis@updateProgress(value= 0, detail = text)
      }
      
    }else{
      continueExploration<-FALSE
    }
  }
  print(paste("ending with ", length(V(g)), " vertices"))
  
  AFMImageNetworksAnalysis@skeletonGraph<-g
  
  return(AFMImageNetworksAnalysis)
}

#' Calculate topology image (TBC)
#'
#' \code{getTopologyAFMImage} return the global topological distance
#' 
#' @param BinaryAFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy in a binary format 0 or 1 values for heigths
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}} from Atomic Force Microscopy
#' @author M.Beauvais
#' @export
getTopologyAFMImage<-function(BinaryAFMImage, AFMImageNetworksAnalysis){
  
  filterVector<-unlist(BinaryAFMImage@data$h)
  
  topology<-c()
  
  
  for (x in 1:BinaryAFMImage@samplesperline) {
    for (y in 1:BinaryAFMImage@lines) {
      if(x==1) {
        bX=seq(from=0, to=BinaryAFMImage@samplesperline-1, by=1)
      }else{
        if (x==BinaryAFMImage@samplesperline) {
          bX=seq(from=x-1, to=0, by=-1)
        }else{
          bX=seq(from=x-1, to=0, by=-1)
          bX=c(bX, seq(from=1, to=BinaryAFMImage@samplesperline-x, by=1))
        }
      }
      # bX
      
      if(y==1) {
        bY=seq(from=0, to=BinaryAFMImage@lines-1, by=1)
      }else{
        if (y==BinaryAFMImage@lines) {
          bY=seq(from=y-1, to=0, by=-1)
        }else{
          bY=seq(from=y-1, to=0, by=-1)
          bY=c(bY, seq(from=1, to=BinaryAFMImage@lines-y, by=1))
        }
      }
      # bY
      
      
      bX=BinaryAFMImage@hscansize*bX
      bY=BinaryAFMImage@vscansize*bY
      
      bX<-matrix(rep(bX,BinaryAFMImage@lines), ncol=BinaryAFMImage@lines, byrow=TRUE )
      bY<-matrix(rep(bY,BinaryAFMImage@samplesperline), ncol=BinaryAFMImage@samplesperline, byrow=FALSE )
      
      nm=as.numeric(1/sqrt(bX^2+bY^2))
      nm[is.infinite(nm)]<-0
      #nm*filterVector
      res<-sum(nm*filterVector)
      topology<-c(topology,res)
      #print(res)
      
    }
  }
  
  
  scanby<-BinaryAFMImage@scansize/BinaryAFMImage@samplesperline
  endScan<-BinaryAFMImage@scansize*(1-1/BinaryAFMImage@samplesperline)
  
  topologyAFMImage<-AFMImage(
    data = data.table(x = rep(seq(0,endScan, by= scanby), times = BinaryAFMImage@lines),
                      y = rep(seq(0,endScan, by= scanby), each = BinaryAFMImage@samplesperline),
                      h = topology),
    samplesperline = BinaryAFMImage@samplesperline, lines = BinaryAFMImage@lines,
    vscansize = BinaryAFMImage@vscansize, hscansize = BinaryAFMImage@hscansize, scansize = BinaryAFMImage@scansize,
    fullfilename = BinaryAFMImage@fullfilename )
  
  
  
  return(topologyAFMImage)
  
}

#' get a segment of points thanks to Bresenham line algorithm
#'
#' \code{getBresenham2DSegment} return the Bresenham segment in 2D from extremities coordinates
#' 
#' @param x1 abscissa coordinates of the first point
#' @param y1 ordinate coordinates of the first point
#' @param x2 abscissa coordinates of the second point
#' @param y2 ordinate coordinates of the second point
#' @return a data.table of points - data.table(x, y)
#' @author M.Beauvais
#' @export
getBresenham2DSegment<-function(x1, y1, x2, y2) {
  resX=c()
  resY=c()
  
  dx<-x2-x1
  dy<-y2-y1
  
  #print(paste("getBresenham2DSegment",dx,dy))
  
  if (dx !=0) {
    if (dx > 0) {
      if (dy !=0) {
        if (dy > 0) {
          if (dx >= dy) {
            e<-dx
            dx <- e  * 2 
            dy <- dy * 2  
            while(TRUE){
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 + 1
              if (x1 == x2) break
              e <- e - dy
              if (e < 0) {
                y1 <- y1 + 1
                e <- e + dx 
              }
            }
          } else {
            e <- dy
            dy <- e * 2
            dx <- dx * 2 
            while(TRUE){ 
              resX=c(resX,x1); resY=c(resY, y1)
              y1 <- y1 + 1
              if (y1 == y2) break
              e <- e - dx
              if (e < 0) {
                x1 <- x1 + 1 
                e <- e + dy
              }
            }
          }
        }else if (dy < 0){ # dy < 0 (et dx > 0)
          
          
          if (dx >= -dy) {
            e <- dx
            dx <- e * 2
            dy <- dy * 2
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 + 1
              if (x1 == x2) break
              e <- e + dy
              if (e < 0) {
                y1 <- y1 - 1 
                e <- e + dx
              }
            }
          } else{
            e <- dy
            dy <- e * 2 
            dx <- dx * 2
            #print(c(e,dy,dx))
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              #print(c(x1, y1))
              y1 <- y1 - 1
              if (y1 == y2) break
              e <- e - dx
              #print(paste(c("e",e)))
              if (e < 0) {
                x1 <- x1 + 1
                if(x1>x2) x1=x2 # MB !!!
                e <- e - dy
                #print(paste(c("e",e)))
              }
            }
          }
          
        }
      }  else if (dy == 0){ # dy = 0 (et dx > 0)
        while(x1 != x2) {
          resX=c(resX,x1); resY=c(resY, y1) 
          x1 <- x1 + 1
        }
      }
    }else if (dx<0) {  # dx < 0
      dy <- y2 - y1
      if (dy != 0) {
        if (dy > 0) {
          if (-dx >= dy) {
            e <- dx
            dx <- e * 2 
            dy <- dy * 2  
            while(TRUE){
              resX=c(resX,x1); resY=c(resY, y1) 
              x1 <- x1 - 1
              if (x1 == x2) break
              e <- e + dy
              if (e >= 0) {
                y1 <- y1 + 1 
                e <- e + dx 
              }
            }
          }else{
            e <- dy
            dy <- e * 2
            dx <- dx * 2 
            while(TRUE){ 
              resX=c(resX,x1); resY=c(resY, y1) 
              y1 <- y1 + 1
              if ( y1 == y2) break 
              e <- e + dx
              if (e <= 0) {
                x1 <- x1 - 1  
                e <- e + dy 
              }
            }
          }
        }else if(dy <0) {  # dy < 0 (et dx < 0)
          if (dx <= dy) {
            e <- dx
            dx <- e * 2 
            dy <- dy * 2  
            while(TRUE){  
              resX=c(resX,x1); resY=c(resY, y1)
              x1 <- x1 - 1
              if (x1 == x2) break
              e <- e - dy
              if (e >= 0) {
                y1 <- y1 - 1
                e <- e + dx 
              }
            }
          } else { 
            e <- dy
            dy <- e * 2 
            dx <- dx * 2 
            
            while(TRUE){
              resX=c(resX,x1); resY=c(resY, y1)
              y1 <- y1 - 1
              if ( y1 == y2 ) break
              e <- e - dx
              if (e >= 0) {
                x1 <- x1 - 1
                e <- e + dy
              }
            }
          }
        } 
      } else if (dy==0) {  # dy = 0 (et dx < 0)
        while(x1!=x2) {
          resX=c(resX,x1); resY=c(resY, y1)
          x1 <- x1 - 1
        }
      }
    }
  } else if (dx==0) {  # dx = 0
    dy <- y2 - y1
    if (dy != 0) {
      if (dy > 0) {
        while(y1 != y2) {
          resX=c(resX,x1); resY=c(resY, y1)
          y1 <- y1 + 1
        } 
        
      } else if (dy < 0) { # dy < 0 (et dx = 0)
        while(y1!=y2) {
          resX=c(resX,x1); resY=c(resY, y1)
          y1 <- y1 - 1
        }
        
      }
      
    }
    
  }
  resX=c(resX,x2); resY=c(resY, y2)
  pts = data.table(x=resX, y=resY)
  
  return(pts)
}

#' identify largest circles in binary image
#'
#' \code{identifyCirclesInBinaryImage} return TRUE if vertex is adjacent to a better vertex
#' 
#' @param AFMImageNetworksAnalysis a \code{\link{AFMImageNetworksAnalysis}}
#' @param smallBranchesTreatment TRUE if the algorithm should try to identify very small branches
#' @return AFMImageNetworksAnalysis the \code{\link{AFMImageNetworksAnalysis}} instance
#' @author M.Beauvais
#' @export
identifyNodesWithCircles<-function(AFMImageNetworksAnalysis, smallBranchesTreatment) {
  spDistsN1<-nbOfCircles<-maxArea<-h<-NULL
    
  newCircleAFMImage<-copy(AFMImageNetworksAnalysis@binaryAFMImage)
  newCircleAFMImage2<-copy(AFMImageNetworksAnalysis@binaryAFMImage)
  circleRadius<-15
  iteration<-0
  rm(avgDT)
  
  while(circleRadius>0) {
    
    iteration=iteration+1
    circleRadius=circleRadius-1
    if (circleRadius>0) {
      center<-c(circleRadius, circleRadius)
      blockSize<-circleRadius*2+1
      pts = SpatialPoints(cbind(rep(0:(blockSize-1),blockSize), rep(0:(blockSize-1),1,each= blockSize)))
      # pts
      nm <- spDistsN1(pts, center, longlat=FALSE)
      # nm
      # nm<circleRadius
      
      # find all blocks in image
      # and check if the circle with biggest radius inside the block exists in the image
      # if yes, set all the height of all the points inside circle to 10
      heights<-newCircleAFMImage@data$h
      
      binaryAFMImageMatrix<-matrix(heights, ncol=newCircleAFMImage@samplesperline)
      newBlockAFMImageMatrix<-matrix(heights, ncol=newCircleAFMImage@samplesperline)
      
      for (x in seq(1:(nrow(binaryAFMImageMatrix)-blockSize))) {
        for (y in seq(1:(ncol(binaryAFMImageMatrix)-blockSize))) {
          tempMatrix<-binaryAFMImageMatrix[x:(x+blockSize),y:(y+blockSize)]
          
          if (all(as.vector(tempMatrix)[nm<=circleRadius] == 1) == TRUE) {
            #print (paste(x,y))
            #newBlockAFMImageMatrix[x:(x+blockSize),y:(y+blockSize)]<-rep(5, (blockSize+1)*(blockSize+1) )
            #newBlockAFMImageMatrix[x+circleRadius+1, y+circleRadius+1]<-10
            newBlockAFMImageMatrix[x+circleRadius, y+circleRadius]<-10
            # as.vector(newBlockAFMImageMatrix[x:(x+blockSize),y:(y+blockSize)])
            # 
            # for(i in seq(x:blockSize)) {
            #   for(j in seq(y:blockSize)) {
            #     print (paste(i,j))
            #     newBlockAFMImageMatrix[i,j]<-5    
            #   }
            # }
          }
        }
      }
      
      
      newBlockAFMImage<-copy(newCircleAFMImage)
      newBlockAFMImage@data$h<-as.vector(newBlockAFMImageMatrix)
      #displayIn3D(newBlockAFMImage)
      
      newBlockAFMImage2<-copy(newBlockAFMImage)
      newBlockAFMImage2@data$h[newBlockAFMImage2@data$h<3]<-0
      #displayIn3D(newBlockAFMImage2)
      
      
      
      newBlockAFMImageMatrix2<-matrix(newBlockAFMImage2@data$h, ncol=newBlockAFMImage2@samplesperline)
      # get coordinates of non 0 elements
      nonZeroElements<-which(newBlockAFMImageMatrix2!=0,arr.ind = T)
      
      
      if(nrow(nonZeroElements)!=0) {
        print(paste("circleRadius",circleRadius))
        lat<-nonZeroElements[,1]
        lon<-nonZeroElements[,2]
        
        DBSCAN <- dbscan(cbind(lat, lon), eps = 1.5, MinPts = 1)
        plot(-lon, -lat, col = DBSCAN$cluster, pch = 20)
        
        #Sys.sleep(10)
        
        # which points is in which cluster ?
        
        #print(DBSCAN$cluster)
        # DBSCAN$eps
        # DBSCAN$minPts
        
        nodesDT=data.table(lon,lat,cluster=DBSCAN$cluster)
        setkeyv(nodesDT, "cluster")
        
        
        # number of cluster
        nbOfClusters<-length(unique(nodesDT$cluster))
        print(paste("nbOfClusters",nbOfClusters))
        if (nbOfClusters>0) {
          # number of points per cluster
          nbPointsPerCluster<-nodesDT[, nbOfCircles:=sum(lon!=0), by = c("cluster")]
          # size of cluster
          nodesDT[, maxArea:=(max(lon)-min(lon))*(max(lat)-min(lat)), by = c("cluster")]
          print(nodesDT)
          
          if (circleRadius>0) {
            nodesToBeRemoved<-nodesDT[maxArea<pi*(circleRadius^2),]
          }else{
            nodesToBeRemoved<-nodesDT
          }
          #nodesToBeRemoved<-nodesDT[nbOfCircles==1]
          #if(nrow(nodesToBeRemoved)!=0) {
          
          print(nodesToBeRemoved)        
          newCircleAFMImage@data$h[nodesToBeRemoved$lon+1+nodesToBeRemoved$lat*newCircleAFMImage@samplesperline]<-0
          
          for(oneCenter in seq(1, nrow(nodesToBeRemoved))) {
            center<-c(nodesToBeRemoved[oneCenter,]$lat, nodesToBeRemoved[oneCenter,]$lon)
            #pts = SpatialPoints(cbind(rep(1:blockSize2,blockSize2)+center[1], rep(1:blockSize2,1,each= blockSize2)+center[2]))
            
            
            # Use a bigger circle that will be removed from image
            # in order to exclude other nodes that could be very near
            circleRadius2=circleRadius+2
            blockSize2=circleRadius2*2+1
            pts = SpatialPoints(cbind(rep(0:(blockSize2-1),blockSize2)+center[1]-circleRadius2, rep(0:(blockSize2-1),1,each= blockSize2)+center[2]-circleRadius2))
            #pts = SpatialPoints(cbind(rep(1:blockSize2,blockSize2)+center[1]-circleRadius2, rep(1:blockSize2,1,each= blockSize2)+center[2]-circleRadius2))
            pts<-pts[pts$coords.x1>0&pts$coords.x1<newCircleAFMImage2@lines&pts$coords.x2>0&pts$coords.x2<newCircleAFMImage2@samplesperline]
            
            nm <- spDistsN1(pts, center, longlat=FALSE)
            # points that are inside the circle
            listOfPointsInsideCircle<-pts[nm<=circleRadius2]
            #listOfPointsInsideCircle$coords.x1
            #newCircleAFMImage@data$h[nodesToBeRemoved$lon+virSize+(nodesToBeRemoved$lat+virSize)*newCircleAFMImage@samplesperline]<-0
            #MB
            newCircleAFMImage@data$h[listOfPointsInsideCircle$coords.x1+1+(listOfPointsInsideCircle$coords.x2)*newCircleAFMImage@samplesperline]<-0
          }
          displayIn3D(newCircleAFMImage)
          # displayIn3D(newCircleAFMImage2)
          
          # for(virSize in seq(1, circleRadius)) {
          #   
          #   
          #   
          #   if (all(as.vector(tempMatrix)[nm<circleRadius] == 1) == TRUE)
          #   
          #   
          #   newCircleAFMImage@data$h[nodesToBeRemoved$lon+virSize+(nodesToBeRemoved$lat+virSize)*newCircleAFMImage@samplesperline]<-0
          #   newCircleAFMImage@data$h[nodesToBeRemoved$lon+virSize+(nodesToBeRemoved$lat-virSize)*newCircleAFMImage@samplesperline]<-0
          #   newCircleAFMImage@data$h[nodesToBeRemoved$lon-virSize+(nodesToBeRemoved$lat+virSize)*newCircleAFMImage@samplesperline]<-0
          #   newCircleAFMImage@data$h[nodesToBeRemoved$lon-virSize+(nodesToBeRemoved$lat-virSize)*newCircleAFMImage@samplesperline]<-0
          # }
          #displayIn3D(newCircleAFMImage)
          
          if (!exists("avgDT")) {  
            avgDT<-nodesToBeRemoved[, c(mean(lon)), by = c("cluster")]
            #avgDT<-data.table(cluster=avgDT$cluster, lon=round(avgDT$V1,0), lat=round(nodesToBeRemoved[, mean(lat), by = c("cluster")]$V1, 0), circleRadius=circleRadius)
            avgDT<-data.table(lon=round(avgDT$V1,0), lat=round(nodesToBeRemoved[, mean(lat), by = c("cluster")]$V1, 0), circleRadius=circleRadius)        
            print(avgDT)
            avgDT2<-copy(avgDT)
          }else{
            avgDT2<-nodesToBeRemoved[, c(mean(lon)), by = c("cluster")]
            #avgDT2<-data.table(cluster=avgDT2$cluster, lon=round(avgDT2$V1,0), lat=round(nodesToBeRemoved[, mean(lat), by = c("cluster")]$V1, 0), circleRadius=circleRadius)
            avgDT2<-data.table(lon=round(avgDT2$V1,0), lat=round(nodesToBeRemoved[, mean(lat), by = c("cluster")]$V1, 0), circleRadius=circleRadius)
            avgDT<-rbind(avgDT,avgDT2)
          }
          #print(avgDT2)
          print(paste("circleRadius=",circleRadius,"- nb of centers",nrow(avgDT2)))
          for(oneCenter in seq(1, nrow(avgDT2))) {
            center<-c(avgDT2[oneCenter,]$lat, avgDT2[oneCenter,]$lon)
            #center<-c(0,0)
            #pts = SpatialPoints(cbind(rep(1:blockSize,blockSize)+center[1]-circleRadius, rep(1:blockSize,1,each= blockSize)+center[2]-circleRadius))
            pts = SpatialPoints(cbind(rep(0:(blockSize2-1),blockSize2)+center[1]-circleRadius2, rep(0:(blockSize2-1),1,each= blockSize2)+center[2]-circleRadius2))
            pts<-pts[pts$coords.x1>0&pts$coords.x1<newCircleAFMImage2@lines&pts$coords.x2>0&pts$coords.x2<newCircleAFMImage2@samplesperline]
            nm <- spDistsN1(pts, center, longlat=FALSE)
            # points that are inside the circle
            listOfPointsInsideCircle<-pts[nm<=circleRadius]
            #print(listOfPointsInsideCircle$coords.x1)
            newCircleAFMImage2@data$h[listOfPointsInsideCircle$coords.x1+1+(listOfPointsInsideCircle$coords.x2)*newCircleAFMImage2@samplesperline]<-180+iteration*50
            
          }
          # }else{
          #   print(paste("all nodes are near each other for", circleRadius))
          # }
        }else{
          print(paste("no circle found for", circleRadius))
        }
      }
    }else{
      
      if (smallBranchesTreatment) {
        print("circleRadius=0 take the left points")
        
        untreatedPoints<-newCircleAFMImage@data[h!=0]
        
        while (nrow(untreatedPoints)>1) {
          
          
          print(paste("nrow(untreatedPoints)",nrow(untreatedPoints)))
          morePoints<-untreatedPoints[1,]
          nodesToBeRemoved<-data.table(lon=morePoints$x/(newCircleAFMImage@hscansize/newCircleAFMImage@samplesperline),
                                       lat=morePoints$y/(newCircleAFMImage@vscansize/newCircleAFMImage@lines),
                                       circleRadius=circleRadius)
          
          print(nodesToBeRemoved)
          #      Sys.sleep(5)
          
          # for(oneCenter in seq(1, nrow(nodesToBeRemoved))) {
          oneCenter<-1
          center<-c(nodesToBeRemoved[oneCenter,]$lat, nodesToBeRemoved[oneCenter,]$lon)
          #pts = SpatialPoints(cbind(rep(1:blockSize2,blockSize2)+center[1], rep(1:blockSize2,1,each= blockSize2)+center[2]))
          
          
          # Use a bigger circle that will be removed from image
          # in order to exclude other nodes that could be very near
          circleRadius2=circleRadius+4
          blockSize2=circleRadius2*2+1
          pts = SpatialPoints(cbind(rep(0:(blockSize2-1),blockSize2)+center[1]-circleRadius2, rep(0:(blockSize2-1),1,each= blockSize2)+center[2]-circleRadius2))
          
          # remove points that are outside and on the edge
          edgePixelThickness<-2
          pts<-pts[pts$coords.x1>edgePixelThickness&
                     pts$coords.x1<(newCircleAFMImage2@lines-edgePixelThickness)&
                     pts$coords.x2>edgePixelThickness&
                     pts$coords.x2<(newCircleAFMImage2@samplesperline-edgePixelThickness)]
          
          nm <- spDistsN1(pts, center, longlat=FALSE)
          # points that are inside the circle
          listOfPointsInsideCircle<-pts[nm<=circleRadius2,]
          # add center of the circle
          listOfPointsInsideCircle<-rbind(listOfPointsInsideCircle,
                                          SpatialPoints(cbind(center[1], center[2]))
          )
          
          
          
          #newCircleAFMImage@data$h[listOfPointsInsideCircle$coords.x1+(listOfPointsInsideCircle$coords.x2)*newCircleAFMImage@samplesperline]<-0
          # remove center too
          # newCircleAFMImage@data[x %in% morePoints$x& y %in% morePoints$y,]$h<-0
          
          diskPoints<-data.table(x=listOfPointsInsideCircle$coords.x2*newCircleAFMImage@vscansize/newCircleAFMImage@lines,
                                 y=listOfPointsInsideCircle$coords.x1*newCircleAFMImage@hscansize/newCircleAFMImage@samplesperline)
          #print(diskPoints)
          newCircleAFMImage@data[x %in% diskPoints$x & y %in% diskPoints$y,]$h<-0
          
          # }
          #displayIn3D(newCircleAFMImage)
          # displayIn3D(newCircleAFMImage2)
          
          nodesToBeRemoved<-data.table(lon=nodesToBeRemoved$lat, lat=nodesToBeRemoved$lon, circleRadius=nodesToBeRemoved$circleRadius)
          if (!exists("avgDT")) {
            avgDT<-copy(nodesToBeRemoved)
            #avgDT$h<-NULL
            #avgDT<-cbind(nodesToBeRemoved, data.table( circleRadius=circleRadius))
            print(avgDT)
            avgDT2<-copy(avgDT)
          }else{
            avgDT2<-copy(nodesToBeRemoved)
            #avgDT2$h<-NULL
            avgDT<-rbind(avgDT,avgDT2)
          }
          #print(avgDT2)
          
          newCircleAFMImage2@data$h[nodesToBeRemoved$lat+1+(nodesToBeRemoved$lon)*newCircleAFMImage2@samplesperline]<-180+iteration*50
          
          untreatedPoints<-newCircleAFMImage@data[h!=0]
        }
        print(paste("circleRadius=",circleRadius,"- nb of centers",nrow(avgDT2)))
        
      }
      
      
      
    }
  }
  
  AFMImageNetworksAnalysis@binaryAFMImageWithCircles<-copy(newCircleAFMImage2)
  AFMImageNetworksAnalysis@circlesTable<-copy(avgDT)
  return(AFMImageNetworksAnalysis)
}

####################################################################################################

#' getIntersectionPointWithBorder to be described
#'
#' \code{getIntersectionPointWithBorder} return a data.table
#' 
#' @param AFMImage a \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param center center
#' @param r radius
#' @param deg degree
#' @author M.Beauvais
#' @export
getIntersectionPointWithBorder<-function(AFMImage, center, r, deg) {
  theta <- (deg * pi) / (180)
  x = center$lon + r * cos(theta)
  y = center$lat + r * sin(theta)
  
  pt=data.table(lat=y, lon=x)
  return(pt)
}

#' get a triangle starting from center, two segments of length r with angles deg1 and deg2 
#'
#' \code{getTriangle} return a data.table points of a triangle
#' 
#' @param AFMImage a \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param center center
#' @param r length of segment
#' @param deg1 angle 1
#' @param deg2 angel 2
#' @author M.Beauvais
#' @export
getTriangle<-function(AFMImage, center, r, deg1, deg2) {
  pt1=getIntersectionPointWithBorder(AFMImage, center, r, deg1)
  pt2=getIntersectionPointWithBorder(AFMImage, center, r, deg2)
  
  trianglePts=data.table(lon=c(center$lon, pt1$lon, pt2$lon,center$lon), lat=c(center$lat, pt1$lat, pt2$lat,center$lat))
  return(trianglePts)
}

#' existsSegment checks if a segment exists in an AFMImage; check if all the heights at the segment coordinates are different to zero.
#'
#' \code{existsSegment} return a boolean
#' 
#' @param AFMImage a \code{\link{AFMImage}} from Atomic Force Microscopy or a binary \code{\link{AFMImage}}
#' @param segment a data.table coming from the getBresenham2Dsegment #TODO Segment class
#' @return TRUE if all the heights of the segment are different from zero
#' @author M.Beauvais
#' @export
existsSegment<-function(AFMImage, segment) {
  #print(segment)
  res<-!any(AFMImage@data$h[segment$y+1+segment$x*AFMImage@samplesperline]==0)
  #print(res)
  return(res)
}
#test existsSegment(binaryAFMImage,       segment= getBresenham2DSegment(10, 9,11, 9))
# existsSegment(binaryAFMImage, segment= getBresenham2DSegment(504,358,511,335))
# binaryAFMImage@samplesperline
# segment= getBresenham2DSegment(504,358,511,335)
# segment
# binaryAFMImage@data$h[segment$y+1+segment$x*binaryAFMImage@samplesperline]

#################################################################################
#
# center$lon and center$lat
#Test
#getCircleSpatialPoints(binaryAFMImage, data.table(lon=226, lat=344), 10)
#getCircleSpatialPoints(binaryAFMImage, center= data.table(lon=20, lat=10), circleRadius=5)
#getCircleSpatialPoints(binaryAFMImage, center= data.table(lon=20, lat=10), circleRadius=0)
#getCircleSpatialPoints(binaryAFMImage, center= data.table(lon=20, lat=10), circleRadius=1)
#center= data.table(lon=20, lat=10)

#' get the spatial points on the circle including the center of the circle
#' 
#' @param binaryAFMImage a binary \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param center the center of the circle with center$lon as the x coordinates and center$lat as the y coordinates
#' @param circleRadius the radius of the circle
#' @return a \code{\link{SpatialPoints}} object of all the points of the circle including the center of the circle  
#' @export
#' @author M.Beauvais
getCircleSpatialPoints<-function(binaryAFMImage, center, circleRadius) {
  spDistsN1<-NULL
  
  if (circleRadius<0) {
    stop("getCircleSpatialPoints - the radius is inferior to 0")
    return()
  }
  if (circleRadius>0) {
    blockSize<-circleRadius*2+1
    
    pts = SpatialPoints(cbind(rep(1:blockSize,blockSize)+center$lon-circleRadius-1, rep(1:blockSize,1,each= blockSize)+center$lat-circleRadius-1))
    #print(pts)
    pts<-pts[pts$coords.x1>0&pts$coords.x1<binaryAFMImage@lines&pts$coords.x2>0&pts$coords.x2<binaryAFMImage@samplesperline]
    #plot(pts)
    nm <- spDistsN1(matrix(c(pts$coords.x1, pts$coords.x2), ncol=2), c(center$lon, center$lat), longlat=FALSE)
    #print(nm)
    
    centerAllpoints<-pts[nm==circleRadius]
    #plot(centerAllpoints)
    
    centerAllpoints<-SpatialPoints(cbind(
      c(centerAllpoints$coords.x1, center$lon),
      c(centerAllpoints$coords.x2, center$lat)
    ))
  }else{
    centerAllpoints<-SpatialPoints(cbind(center$lon, center$lat))
  }
  return(centerAllpoints)
}

#################################################################################
# test
#Test
#AreNodesConnected(binaryAFMImage, data.table(lon=226, lat=344), 10, data.table(lon=25, lat=344), 5)
#AreNodesConnected(binaryAFMImage, data.table(lon=226, lat=344), 10, data.table(lon=25, lat=344), 0)
#AreNodesConnected(binaryAFMImage, center1, circleRadius1, data.table(lon=pt$coords.x1, lat=pt$coords.x2), pt$circleRadius)


# AreNodesConnected(binaryAFMImage, data.table(lon=76, lat=60), 1, data.table(lon=79, lat=65), 0)
# AreNodesConnected(binaryAFMImage, data.table(lon=76, lat=60), 0, data.table(lon=79, lat=65), 0)
# 
# circle1AllPoints<-getCircleSpatialPoints(binaryAFMImage, data.table(lon=76, lat=60), 1)
# circle1AllPoints<-circle1AllPoints[which(binaryAFMImage@data$h[circle1AllPoints$coords.x2+1+circle1AllPoints$coords.x1*binaryAFMImage@samplesperline]!=0)]
# circle1AllPoints
# 
# existsSegment(binaryAFMImage, segment= getBresenham2DSegment(76,60,79,65))
# binaryAFMImage@data$h[segment$y+1+segment$x*binaryAFMImage@samplesperline]
# which(binaryAFMImage@data$h[segment$y+1+segment$x*AFMImage@samplesperline]!=0)

#' check if nodes represented by circles are connected. The function defines all the possible segments between the circles and check if at least one segment exists.
#' 
#' @param binaryAFMImage a binary \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param center1 the center of the circle with center$lon as the x coordinates and center$lat as the y coordinates
#' @param radius1 the radius of the circle
#' @param center2 the center of the circle with center$lon as the x coordinates and center$lat as the y coordinates
#' @param radius2 the radius of the circle
#' @return TRUE if the nodes are connected
#' @export
#' @author M.Beauvais
AreNodesConnected<-function(binaryAFMImage, center1, radius1, center2, radius2) {
  # print(center1)
  # print(radius1)
  # print(center2)
  # print(radius2)
  
  if (radius1>0)  circle1AllPoints<-getCircleSpatialPoints(binaryAFMImage, center1, radius1)
  else{
    circle1AllPoints<-getCircleSpatialPoints(binaryAFMImage, center1, 1)
    circle1AllPoints<-circle1AllPoints[which(binaryAFMImage@data$h[circle1AllPoints$coords.x2+1+circle1AllPoints$coords.x1*binaryAFMImage@samplesperline]!=0)]
  }
  # print(circle1AllPoints)
  # print(length(circle1AllPoints))
  #plot(circle1AllPoints)
  
  if (radius1>0)  circle2AllPoints<-getCircleSpatialPoints(binaryAFMImage, center2, radius2)
  else{
    circle2AllPoints<-getCircleSpatialPoints(binaryAFMImage, center2, 1)
    circle2AllPoints<-circle2AllPoints[which(binaryAFMImage@data$h[circle2AllPoints$coords.x2+1+circle2AllPoints$coords.x1*binaryAFMImage@samplesperline]!=0)]
  }
  
  
  # print(circle2AllPoints)
  # print(length(circle2AllPoints))
  
  for (circlePt1Nb in seq(1, length(circle1AllPoints))) {
    circlePt1<-circle1AllPoints[circlePt1Nb,]
    
    for (circlePt2Nb in seq(1, length(circle2AllPoints))) {
      circlePt2<-circle2AllPoints[circlePt2Nb,]
      segment<-getBresenham2DSegment(circlePt1$coords.x1, circlePt1$coords.x2,
                                     circlePt2$coords.x1, circlePt2$coords.x2)
      #print(segment)
      if (existsSegment(binaryAFMImage, segment)) {
        print(paste("segment exists",center1$lon, center1$lat,":",center2$lon, center2$lat))
        return(TRUE)
      }
    }
  }
  return(FALSE)
}
##############################################################################################

#' calculate the angle between two vectors
#' 
#' @param x a vector
#' @param y a vector
#' @return the angle between the vectors
#' @export
#' @author M.Beauvais
getAngle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  # print(dot.prod)
  # print(norm.x)
  # print(norm.y)
  theta <- acos(dot.prod / (norm.x * norm.y))
  if (is.nan(theta)) theta=0
  return(as.numeric(theta))
}
# test
# getAngle(c(2,12), c(1,6))
# getAngle(c(2,12), c(4,24))

#' check if all the angles between one edge and a list of edges is superior to a specified value.
#' 
#' @param binaryAFMImage a binary \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param edge1 one edge
#' @param edges2 list of edges
#' @param minAngle the minimum angle value 
#' @return TRUE if all the angle are superior to the specified value
#' @export
#' @author M.Beauvais
isAngleBetweenEdgesAlwaysSuperiorToMinAngle<-function(binaryAFMImage, edge1, edges2, minAngle) {
  #print(edge1)
  #print(edges2)
  
  coordsFromEdge1=getCoordinatesFromVertexId(binaryAFMImage, as.numeric(edge1$from))
  coordsToEdge1=getCoordinatesFromVertexId(binaryAFMImage, as.numeric(edge1$to))
  
  coordsFromEdges2=getCoordinatesFromVertexId(binaryAFMImage, as.numeric(edges2$from))
  coordsToEdges2=getCoordinatesFromVertexId(binaryAFMImage, as.numeric(edges2$to))
  
  # allYCoordinates<-cbind(coordsFromEdges2, coordsToEdges2)
  # print(allYCoordinates)
  x=c(coordsToEdge1$coords.x1-coordsFromEdge1$coords.x1, 
      coordsToEdge1$coords.x2-coordsFromEdge1$coords.x2)
  
  for (y in seq(1,nrow(coordsToEdges2))){
    y=c(coordsToEdges2[y,]$coords.x1-coordsFromEdges2[y,]$coords.x1, coordsToEdges2[y,]$coords.x2-coordsFromEdges2[y,]$coords.x2)
    
    angle<-getAngle(x,y)
    if (angle>pi) angle<-angle-pi
    # print(angle)
    #print(paste("x=",x,"y=",y,"angle=",180*angle/pi, "degrees"))
    if (angle<minAngle) {
      #print(paste(c(x,"->",y)))
      return(FALSE)
    }
  }
  return(TRUE)
}
# isAngleBetweenEdgesAlwaysSuperiorToMinAngle(edge1=data.table(from=vid1, to=vid2), edges2=existingEdgesVid1,0.52)
# 
# existingEdges<-data.table(from = c("6553685"), to = c("2097229"),arrows = c("to"))
# isAngleBetweenEdgesAlwaysSuperiorToMinAngle(edge1=data.table(from=1835079, to=6553685), edges2=existingEdges,0.52)
# isAngleBetweenEdgesAlwaysSuperiorToMinAngle(edge1=data.table(from=6553685, to=1835079), edges2=existingEdges,0.52)

#' display the network of nodes and edges
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param edges list of edges
#' @param isolates list of isolated edges
#' @export
#' @author M.Beauvais
displaygridIgraphPlotFromEdges<-function(AFMImage, edges, isolates) {
  #print(edges)
  alledges2<-as.vector(t(matrix(c(edges$from,edges$to),ncol=2)))
  vnodes<-unique(c(edges$from, edges$to))
  vnodes<-data.frame(id=vnodes, label = vnodes)
  
  # coords=getCoordinatesFromVertexId(binaryAFMImage, as.numeric(levels(vnodes$id))[vnodes$id])
  # coords
  # coords[order(coords.x1),]
  
  g<-graph(edges= alledges2,isolates=isolates, directed=FALSE)
  gridIgraphPlot(AFMImage, g)
}


#' display the network of nodes and edges
#' 
#' @param AFMImageNetworksAnalysis an \code{\link{AFMImageNetworksAnalysis}}
#' @export
#' @author M.Beauvais
displaygridIgraphPlot<-function(AFMImageNetworksAnalysis) {
  displaygridIgraphPlotFromEdges(AFMImageNetworksAnalysis@binaryAFMImage,
                                 AFMImageNetworksAnalysis@fusionedNodesEdgesTable,
                                 AFMImageNetworksAnalysis@isolatedNodesList)
  
}



#' display the network of nodes and edges
#' 
#' @param AFMImageNetworksAnalysis a \code{\link{AFMImageNetworksAnalysis}}
#' @param MAX_DISTANCE the maximum distance between nodes to check if nodes are connected. Default value is 40.
#' @param ANGLE_MIN the minimum angle between edges to accept or reject an edge
#' @export
#' @author M.Beauvais
identifyEdgesFromCircles<-function(AFMImageNetworksAnalysis, MAX_DISTANCE=40, ANGLE_MIN=(pi/3.25)) {
  spDistsN1<-from<-to<-NULL
  
  TRACEID=6000000
  #ANGLE_MIN<-pi/5
  
  #binaryAFMImage, avgDT,
  
  
  vedges<-data.table(from = c(""), to = c(""),arrows = c("to"))
  edgesDiscarded=data.table(from = c(""), to = c(""),arrows = c("to"))
  
  alledges<-c()
  allvertices<-c()
  
  # for all the nodes of the future networks
  
  # for all the points in the circles in the plot
  # identify if a link is available with all the plot
  
  allNodesAsSpatialPoints = SpatialPoints(cbind(AFMImageNetworksAnalysis@circlesTable$lon, AFMImageNetworksAnalysis@circlesTable$lat))
  allNodesAsSpatialPoints$circleRadius<-AFMImageNetworksAnalysis@circlesTable$circleRadius
  for (centerId in seq(1, nrow(AFMImageNetworksAnalysis@circlesTable))) {
    #centerId=6
    #centerId=2
    print(paste0(centerId," / ", nrow(AFMImageNetworksAnalysis@circlesTable)))
    
    center1= AFMImageNetworksAnalysis@circlesTable[centerId,]
    circleRadius1=AFMImageNetworksAnalysis@circlesTable[centerId,]$circleRadius
    
    vid1<-getVertexId(AFMImageNetworksAnalysis@binaryAFMImage,center1$lon, center1$lat)
    allvertices<-c(allvertices, vid1)
    
    # Use a bigger circle that will be removed from image
    # in order to exclude other nodes that could be very near
    
    
    
    # no filter by taking only the closest points
    # calculate distance with all other points
    
    # take only the points that are no farther than xx
    #allNodesAsSpatialPoints$coords.x1
    #which(allNodesAsSpatialPoints$coords.x1!=center1$lon&allNodesAsSpatialPoints$coords.x2!=center1$lat)
    #otherNodes<-allNodesAsSpatialPoints[which(allNodesAsSpatialPoints$coords.x1!=center1$lon|allNodesAsSpatialPoints$coords.x2!=center1$lat),]
    #otherNodes<-allNodesAsSpatialPoints[which(allNodesAsSpatialPoints$coords.x1!=center1$lon&allNodesAsSpatialPoints$coords.x2!=center1$lat),]  
    otherNodes<-allNodesAsSpatialPoints[-centerId,]  
    #otherNodes<-allNodesAsSpatialPoints[!(allNodesAsSpatialPoints$coords.x1==center1$lon&allNodesAsSpatialPoints$coords.x2==center1$lat)]
    otherNodes$dist<-spDistsN1(matrix(c(otherNodes$coords.x1, otherNodes$coords.x2), ncol=2), c(center1$lon, center1$lat), longlat=FALSE)
    #otherNodes
    
    otherNodes<-otherNodes[with(otherNodes, order(otherNodes@data$dist)), ]
    otherNodes<-otherNodes[otherNodes@data$dist<MAX_DISTANCE,]
    #print(otherNodes)
    if (centerId== TRACEID) print(otherNodes)
    #○if (centerId> TRACEID) return()
    if (centerId>= TRACEID) Sys.sleep(6)
    
    if (vid1 %in% c("18350160")) print("hhhhh")
    #Sys.sleep(1.5)
    for (centerId2Nb in seq(1, nrow(otherNodes))) {
      #centerId2Nb<-1
      #centerId2Nb<-10
      pt<-otherNodes[centerId2Nb,]
      vid2<-getVertexId(AFMImageNetworksAnalysis@binaryAFMImage,pt$coords.x1, pt$coords.x2)
      
      
      if (nrow(edgesDiscarded[from %in% c(vid1, vid2)& to %in% c(vid1, vid2)])==0) {
        # print("edge not discarded")
        # print(c(vid1, vid2))
        #print("edgesDiscarded")
        #print(edgesDiscarded)
        
        if (AreNodesConnected(AFMImageNetworksAnalysis@binaryAFMImage, center1, circleRadius1, data.table(lon=pt$coords.x1, lat=pt$coords.x2), pt$circleRadius)) {
          
          #print(paste("segment exists",center1$lon, center1$lat,":",pt$coords.x1, pt$coords.x2))
          #listOfSegments=c(listOfSegments, segment)
          #ncoui=ncoui+1
          
          # find all the edges starting by vid1 or vid2
          #existingEdges<-vedges[from %in% c(vid1, vid2) | to %in% c(vid1, vid2)]
          existingEdges<-vedges[from %in% c(vid1) | to %in% c(vid2)]
          #existingEdges<-rbind(existingEdges, edgesDiscarded[from %in% c(vid1) | to %in% c(vid2)])
          existingEdges<-rbind(existingEdges, edgesDiscarded[from %in% c(vid1)])
          
          existingEdges<-existingEdges[!(from %in% vid1 & to %in% vid2)]
          existingEdges<-existingEdges[!(from %in% vid2 & to %in% vid1)]
          
          existingEdgesVid1<-existingEdges[from %in% vid1 | to %in% vid1]
          existingEdgesVid2<-existingEdges[from %in% vid2 | to %in% vid2]
          
          
          # print("existingEdgesVid1")
          # print(existingEdgesVid1)
          # print(nrow(existingEdgesVid1))
          # print("existingEdgesVid2")
          # print(existingEdgesVid2)
          # print(nrow(existingEdgesVid2))
          
          
          keepEdge<-TRUE
          # compare angles
          if ((is.data.table(existingEdgesVid1) && nrow(existingEdgesVid1)!=0)&&
              !isAngleBetweenEdgesAlwaysSuperiorToMinAngle(
                edge1=data.table(from=vid1, to=vid2), edges2=existingEdgesVid1, minAngle= ANGLE_MIN)
          ) {
            # print("do not keep vd1")
            keepEdge<-FALSE
          }
          if ((is.data.table(existingEdgesVid2) && nrow(existingEdgesVid2)!=0) &&
              !isAngleBetweenEdgesAlwaysSuperiorToMinAngle(
                edge1=data.table(from=vid1, to=vid2), edges2=existingEdgesVid2, minAngle= ANGLE_MIN)
          ) {
            # print("do not keep vd2")
            keepEdge<-FALSE
          }
          
          
          
          if (keepEdge){
            # else keep the new edge        
            vedges<-rbind(vedges, data.table(from = vid1, to = vid2,arrows = c("to")))
            edgesDiscarded<-rbind(edgesDiscarded, data.table(from = vid2, to = vid1,arrows = c("to")))
            
            alledges=c(alledges,vid1,vid2)
            #pts = SpatialPoints(cbind(pt$coords.x1, pt$coords.x2))        
            #result<-rbind(result, pts)
            #print("keep edge")
            displaygridIgraphPlotFromEdges(AFMImageNetworksAnalysis@binaryAFMImage, edges=vedges[-1,],  isolates = c())
            #cbind(getCoordinatesFromVertexId(AFMImageNetworksAnalysis@binaryAFMImage,vedges$from),getCoordinatesFromVertexId(AFMImageNetworksAnalysis@binaryAFMImage,vedges$to))
            #if (centerId2Nb> 6) Sys.sleep(30)
            if (centerId== TRACEID) Sys.sleep(2*1)
          }else{
            # angle is too low discard edge
            edgesDiscarded<-rbind(edgesDiscarded, data.table(from = vid1, to = vid2,arrows = c("to")))
            #print("discard edge")
            if (centerId== TRACEID) Sys.sleep(3*1)
            
          }
          
          #displaygridIgraphPlotFromEdges(AFMImageNetworksAnalysis@binaryAFMImage, vedges[-1,])
          #cbind(getCoordinatesFromVertexId(AFMImageNetworksAnalysis@binaryAFMImage,vedges$from),getCoordinatesFromVertexId(AFMImageNetworksAnalysis@binaryAFMImage,vedges$to))
          #Sys.sleep(3)
        }else{
          # print(paste("segment does not exist",center1$lon, center1$lat,":",pt$coords.x1, pt$coords.x2))
        }
      }else{
        # print(paste(c("edge already discarded",center1, "->",pt$coords.x1, pt$coords.x2,"-",vid1,"->",vid2)))
      }
    }
    #plot(result)
  }
  
  AFMImageNetworksAnalysis@edgesTable<-copy(vedges[-1,])
  return(AFMImageNetworksAnalysis)
}


######################
# manage the fusion of nodes which circles instersect
# keep all the circles, manage a fusion table
# node id / fusion id

#' fusion the nodes that are intersecting
#' 
#' @param AFMImageNetworksAnalysis the AFMImageNetworksAnalysis instance
#' @return a list of edges with fusioned nodes
#' @export
#' @author M.Beauvais
fusionCloseNodes<-function(AFMImageNetworksAnalysis) {
  
  spDistsN1<-group<-mean_lon<-lon<-mean_lat<-lat<-vertexId<-from<-to<-vedges<-NULL
  
  AFMImageNetworksAnalysis@circlesTable
  AFMImageNetworksAnalysis@circlesTable$group<-rep(0, nrow(AFMImageNetworksAnalysis@circlesTable))
  groupNumber<-0
  for (centerId in seq(1, nrow(AFMImageNetworksAnalysis@circlesTable))) {
    #centerId=6
    print(paste0(centerId," / ", nrow(AFMImageNetworksAnalysis@circlesTable)))
    
    center<- AFMImageNetworksAnalysis@circlesTable[centerId,]
    
    radiusVector<-center$circleRadius+AFMImageNetworksAnalysis@circlesTable$circleRadius
    
    distVector<-spDistsN1(matrix(c(AFMImageNetworksAnalysis@circlesTable$lon,AFMImageNetworksAnalysis@circlesTable$lat),ncol=2),
                          matrix(c(center$lon,center$lat),ncol=2),
                          longlat=FALSE)
    intersectVector<-distVector-radiusVector-2
    
    # print(radiusVector)
    # print(distVector)
    print(intersectVector)
    listOfIntersect<-which(intersectVector<0)
    #print(listOfIntersect)
    if (length(listOfIntersect)>1) {
      print("to be grouped")
      if (all(AFMImageNetworksAnalysis@circlesTable[listOfIntersect,]$group==0)) {
        groupNumber<-groupNumber+1
        AFMImageNetworksAnalysis@circlesTable[listOfIntersect,]$group<-groupNumber
      }else{
        print("special")
        #print(AFMImageNetworksAnalysis@circlesTable[listOfIntersect&group!=0])
        existingGroupNumber<-AFMImageNetworksAnalysis@circlesTable[listOfIntersect,][group!=0,][1]$group
        AFMImageNetworksAnalysis@circlesTable[listOfIntersect,]$group<-existingGroupNumber
      }
      #print(AFMImageNetworksAnalysis@circlesTable[listOfIntersect])
    }
    
  }
  AFMImageNetworksAnalysis@circlesTable
  
  
  nbOfNodesToFusion<-length(unique(AFMImageNetworksAnalysis@circlesTable[group!=0,]$group))
  nbOfNodesToFusion
  if (nbOfNodesToFusion>0) {
    # define new coordinates for all points
    # wh<-which(AFMImageNetworksAnalysis@circlesTable$group==0)
    # AFMImageNetworksAnalysis@circlesTable[wh]$new_lat<-AFMImageNetworksAnalysis@circlesTable[wh]$lat
    # AFMImageNetworksAnalysis@circlesTable[wh]$new_lon<-AFMImageNetworksAnalysis@circlesTable[wh]$lon
    # AFMImageNetworksAnalysis@circlesTable
    
    AFMImageNetworksAnalysis@circlesTable[, mean_lon:=floor(mean(lon)), by=group] 
    AFMImageNetworksAnalysis@circlesTable[, mean_lat:=floor(mean(lat)), by=group]
    AFMImageNetworksAnalysis@circlesTable[group==0, mean_lon:=lon] 
    AFMImageNetworksAnalysis@circlesTable[group==0, mean_lat:=lat] 
    AFMImageNetworksAnalysis@circlesTable
    
    
    # define edge correspondance table
    
    newvedges<-data.table(vertexId=getVertexId(AFMImageNetworksAnalysis@binaryAFMImage, AFMImageNetworksAnalysis@circlesTable[group!=0,]$lon, AFMImageNetworksAnalysis@circlesTable[group!=0,]$lat),
                          new_vertexId=getVertexId(AFMImageNetworksAnalysis@binaryAFMImage, AFMImageNetworksAnalysis@circlesTable[group!=0,]$mean_lon, AFMImageNetworksAnalysis@circlesTable[group!=0,]$mean_lat))
    newvedges
    setkey(newvedges, vertexId)
    
    # tranform the isolated nodes
    isolates<-AFMImageNetworksAnalysis@isolatedNodesList
    isolates %in% newvedges$vertexId
    newvedges
    onewh<-which(isolates %in% newvedges$vertexId)
    for(index in onewh) {
      print(index)
      oldvertexId<-isolates[index]
      print(oldvertexId)
      newVertexId<-newvedges[vertexId %in% oldvertexId]$new_vertexId
      print(newVertexId)
      isolates<-replace(isolates, isolates==oldvertexId, as.character(newVertexId))
    }
    isolates<-unique(isolates)
    isolates
    
    # tranform the edges with the fusioned edge
    newvedges2<-copy(AFMImageNetworksAnalysis@edgesTable)
    newvedges2
    
    onewh<-which(newvedges2$from %in% newvedges$vertexId)
    for(index in onewh) {
      print(index)
      oldvertexId<-newvedges2[index,]$from
      print(oldvertexId)
      newVertexId<-newvedges[vertexId %in% oldvertexId,]$new_vertexId
      print(newVertexId)
      newvedges2[index, from:=as.character(newVertexId)]
    }
    newvedges2
    
    onewh<-which(newvedges2$to %in% newvedges$vertexId)
    for(index in onewh) {
      print(index)
      oldvertexId<-newvedges2[index,]$to
      print(oldvertexId)
      newVertexId<-newvedges[vertexId %in% oldvertexId,]$new_vertexId
      print(newVertexId)
      newvedges2[index, to:=as.character(newVertexId)]
    }
    newvedges2
  }else{
    newvedges2<-vedges
  }
  #print(newvedges2)
  
  
  AFMImageNetworksAnalysis@fusionedNodesCorrespondance<-copy(newvedges)
  if (typeof(newvedges2) %in% c("data.table")) {
    AFMImageNetworksAnalysis@fusionedNodesEdgesTable<-copy(newvedges2)
  }else{
    AFMImageNetworksAnalysis@fusionedNodesEdgesTable<-copy(AFMImageNetworksAnalysis@edgesTable) 
  }
  return(AFMImageNetworksAnalysis)
}

#' identify isolated nodes comparing the list of edges and the list of nodes
#' 
#' @param AFMImageNetworksAnalysis the AFMImageNetworksAnalysis instance
#' @return the updated instance of AFMImageNetworksAnalysis
#' @export
#' @author M.Beauvais
identifyIsolatedNodes<-function(AFMImageNetworksAnalysis) {
  isolates<-getVertexId(AFMImageNetworksAnalysis@binaryAFMImage, AFMImageNetworksAnalysis@circlesTable$lon, AFMImageNetworksAnalysis@circlesTable$lat)
  print(isolates)
  vedges<-AFMImageNetworksAnalysis@edgesTable
  AFMImageNetworksAnalysis@isolatedNodesList<-isolates[!isolates %in% vedges$from & !isolates %in% vedges$to]
  return(AFMImageNetworksAnalysis)
}
# AFMImageNetworksAnalysis<-identifyNodesWithCircles(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis,
#                                                    smallBranchesTreatment = TRUE)
# AFMImageNetworksAnalysis<-identifyEdgesFromCircles(AFMImageNetworksAnalysis= AFMImageNetworksAnalysis)
# AFMImageNetworksAnalysis<-identifyIsolatedNodes(AFMImageNetworksAnalysis)
# AFMImageNetworksAnalysis<-getEdgesAfterNodesFusion(AFMImageNetworksAnalysis)


#' calculate the physical distances between nodes
#' 
#' @param AFMImage an \code{\link{AFMImage}} from Atomic Force Microscopy
#' @param pathVidVector a network path
#' @return the physical distance the extrmities of the path
#' @export
#' @author M.Beauvais
calculatePhysicalDistanceFromPath<-function(AFMImage, pathVidVector) {
  physicalDistance<-0
  hscale<-AFMImage@hscansize/AFMImage@samplesperline
  vscale<-AFMImage@vscansize/AFMImage@lines
  #print(c(hscale, vscale))
  vid1<-pathVidVector[1]
  for (pathInd in seq(2, length(pathVidVector))) {
    vid2<-pathVidVector[pathInd]
    vid1Coords<-getCoordinatesFromVertexId(AFMImage, vid1)
    vid2Coords<-getCoordinatesFromVertexId(AFMImage, vid2)
    
    # print(c(vid1Coords, vid2Coords))  
    # print(sqrt((hscale*(vid1Coords$coords.x1-vid2Coords$coords.x1))^2+(vscale*(vid1Coords$coords.x2-vid2Coords$coords.x2))^2))  
    physicalDistance<-physicalDistance+sqrt((hscale*(vid1Coords$coords.x1-vid2Coords$coords.x1))^2+(vscale*(vid1Coords$coords.x2-vid2Coords$coords.x2))^2)
    vid1<-pathVidVector[pathInd]
  }
  #print(physicalDistance)
  return(physicalDistance)
}
# TODO check if strsplit return results
#path<-strsplit(directedConnectedNodesDT[1,]$shortest_path,"-")[[1]]
#calculatePhysicalDistanceFromPath(newAFMImage, path)


#' create the igraph graph from the fusionned nodes and the edges
#' 
#' @param AFMImageNetworksAnalysis a \code{\link{AFMImageNetworksAnalysis}}
#' @param fusioned TRUE if graph is created from fusioned nodes
#' @export
#' @author M.Beauvais
createGraph<-function(AFMImageNetworksAnalysis, fusioned=TRUE) {
  isolatedNodesList<-AFMImageNetworksAnalysis@isolatedNodesList
  from<-to<-NULL
  ultimateNetwork<-copy(AFMImageNetworksAnalysis@fusionedNodesEdgesTable[from!=to,])
  
  isolates<-isolatedNodesList[!isolatedNodesList %in% ultimateNetwork$from & !isolatedNodesList %in% ultimateNetwork$to]
  
  
  alledges2<-as.vector(t(matrix(c(ultimateNetwork$from,ultimateNetwork$to),ncol=2)))
  # vnodes<-unique(c(ultimateNetwork$from, ultimateNetwork$to))
  # vnodes<-data.frame(id=ultimateNetwork, label = ultimateNetwork)
  
  # TODO
  g<-graph(edges=alledges2, directed=FALSE, isolates=isolates)
  if (fusioned) AFMImageNetworksAnalysis@skeletonGraph<-g
  return(AFMImageNetworksAnalysis)
}

#' calculate the  shortest path between adjacent nodes
#' Calculateshrotest path between all nodes of degree superior to 3
#' identify nodes that are connected by nodes which degrees is inferior to 2
#' and calculate the distance between these nodes
#' 
#' @param AFMImageNetworksAnalysis a \code{\link{AFMImageNetworksAnalysis}}
#' @param fusioned TRUE if graph is created from fusioned nodes
#' @export
#' @author M.Beauvais
calculateShortestPaths<-function(AFMImageNetworksAnalysis, fusioned=TRUE) {
  requireNamespace("parallel")
  
  
  workerFunc <- function(vid1index, binaryAFMImage, g, nodesAnalysisDT) {
    #vid1index<-1
    requireNamespace("igraph")
    requireNamespace("AFM")
    requireNamespace("data.table")
    
    directedConnectedNodesDT<-data.table(vid1="", vid2="", shortest_path="",physicalDistance=0)
    nb<-0
    print(paste0(vid1index," / ", nrow(nodesAnalysisDT)))
    #TODO calculate distance matrix
    # vid2index by distance
    # when nbOfShortestPath for vid1 reach degree of node then break
    nbOfShortestPath<-0
    vid1Degree<-nodesAnalysisDT[vid1index,]$node_degree
    
    for (vid2index in seq(1, nrow(nodesAnalysisDT))) {
      print(paste0(vid1index," / ", nrow(nodesAnalysisDT),"-",vid2index))
      vid1<- nodesAnalysisDT[vid1index,]$vid
      vid2<- nodesAnalysisDT[vid2index,]$vid
      if (!vid1 %in% vid2) {
        allPath<-all_shortest_paths(g, vid1, vid2)
        print(allPath)
        print(is.null(allPath$res))
        if (length(allPath$res)>0) {
          for(pathIndex in seq(1,length(allPath$res))) {
            #if (length(allPath$res)>0) {
            path<-allPath$res[[pathIndex]]$name
            #TODO is it working if two points are in separate graph ?
            numberOfNodesInShortestPath<-length(which(path %in% nodesAnalysisDT$vid == TRUE))
            #print(allPath$res[[1]]$name %in% nodesAnalysisDT$vid)
            if (numberOfNodesInShortestPath==2) {
              #keep
              print(c("interresting", vid1, vid2))
              nbOfShortestPath<-nbOfShortestPath+1
              print(path)
              physicalDistance<-calculatePhysicalDistanceFromPath(binaryAFMImage, path)
              print(physicalDistance)
              #totalPhysicalDistance<-totalPhysicalDistance+physicalDistance
              
              # TODO fill with reverse path
              directedConnectedNodesDT<-rbindlist(list(directedConnectedNodesDT, 
                                                       data.table(vid1=vid1, vid2=vid2, shortest_path=paste0(path, collapse = "-"),physicalDistance=physicalDistance)))
              if (nbOfShortestPath==vid1Degree) break;
            }
            #}
          }
        }
      }
    }
    return(directedConnectedNodesDT[-1,]) 
  }
  
  # start   
  #☺fusioned<-TRUE
  if (fusioned) g<-AFMImageNetworksAnalysis@skeletonGraph
  else g<-AFMImageNetworksAnalysis@originalGraph
  
  node_degree<-NULL
  verticesAnalysisDT<-data.table(vid=V(g)$name, node_degree=unname(degree(g)))
  nodesAnalysisDT<-copy(verticesAnalysisDT[node_degree>2])
  values <- seq(1, nrow(nodesAnalysisDT))
  #values <- seq(1, 28)
  
  ## Number of workers (R processes) to use:
  # Calculate the number of cores
  numWorkers <- parallel::detectCores() - 1
  ## Set up the ’cluster’
  cl <- parallel::makeCluster(numWorkers, type = "PSOCK")
  
  start.time <- Sys.time()
  print(start.time)
  ## Parallel calculation (parLapply):
  res <- parallel::parLapply(cl, values, workerFunc, 
                   AFMImageNetworksAnalysis@binaryAFMImage, 
                   g, 
                   nodesAnalysisDT)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  ## Shut down cluster
  parallel::stopCluster(cl)
  directedConnectedNodesDT<-rbindlist(res)
  
  AFMImageNetworksAnalysis@shortestPaths<-directedConnectedNodesDT
  
  return(AFMImageNetworksAnalysis)
}

