#' @title \code{plot_bipartite_graph}
#'
#' @description \code{plot_bipartite_graph} allows you to plot to networks 
#'              as a bipartite graph 
#'
#' @param g: data.frame with 2 columns, 'From' = starting vertex; 'To' = ending vertex
#' @param circleLayout: boolean indicating whether to use a circle or bipartite layout
#' @param directed: boolean indicating whether to 
#' @param tname: title name

#'
#' @export
plot_bipartite_graph <- function(g, circleLayout=T, directed=F, tname=NULL){
  
  g <- igraph::graph.data.frame(g, directed=directed)
  
  ## motify colors 
  V(g)$type <- igraph::bipartite_mapping(g)$type
  V(g)$color <- ifelse(igraph::V(g)$type, "lightblue", "salmon")
  V(g)$shape <- ifelse(igraph::V(g)$type, "circle", "square")
  E(g)$color <- "lightgray"
  
  ## plot in 2 different formats
  if(circleLayout){
    igraph::plot.igraph(g, layout = igraph::layout.circle, main=tname)
  } else {
    igraph::plot.igraph(g, layout = igraph::layout.bipartite, main=tname)  
  }
  
}

#' @title \code{plot_tripartite_graph}
#'
#' @description \code{plot_bipartite_graph} allows you to plot to networks 
#'              as a bipartite graph 
#'
#' @param g: data.frame with 2 columns, 'From' = starting vertex; 'To' = ending vertex
#' @param directed: boolean indicating whether to 
#' @param tname: title name
#' @export


plot_tripartite_graph <- function(df, directed=T, tname='Ligand-Motif-Gene'){
  
  ## create graphs from data frame
  g <- igraph::graph.data.frame(df, directed=directed)
  
  ## Specify Order of Layers
  layer = rep(3, length(V(g)$name))
  layer[grep("L_",V(g)$name)]=1
  layer[grep("M_",V(g)$name)]=2
  
  names = V(g)$name
  names = sub("G_|L_|M_","", names)
  V(g)$name = names
  
  ## Create Layout 
  layout = layout_with_sugiyama(g, layers=layer)
    ## plot graph
  plot(g,
       layout=cbind(layer, layout$layout[,1]),
       vertex.shape=c("circle","circle","circle")[layer],
       vertex.size=c(0,0,0)[layer],
       vertex.label.dist=c(-3,0,3)[layer],
       vertex.label.degree=0, 
       main=tname)
}