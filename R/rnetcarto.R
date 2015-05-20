#' @useDynLib rnetcarto
NULL
#' Compute modularity and modularity roles for graphs using simulated
#' annealing
#'
#' @docType package
#' @name rnetcarto


#' @title Computes modularity and modularity roles from a network.
#'
#' @param interactions a list of named species 
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export 
netcarto <- function(interactions)
{
    # Number of edges
    E = length(interactions[[1]])

    # Read the weigth if they are supplied
    if (length(interactions)==3){
        weight = interactions[[3]]
    }
    else{
        weight = numeric(E) + 1
    }

    # Convert the species names to integer
    fct = factor(c(interactions[[1]], interactions[[2]]))
    idx = as.integer(fct) - 1L

    # Call rgraphlib
    ans <- .Call("netcarto", idx[1:E], idx[(E+1):(2*E)], weight)

    # Build the dataframe 
    df = data.frame(levels(fct), ans[[1]], ans[[2]], ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
}


#' @title Computes modularity and modularity roles from a bipartite network.
#'
#' @param web a matrix representing the interactions observed between
#' two groups of species (respectively rows and columns).
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export 
bipartmod <- function(web) 
{
    # Removing empty columns and lines. 
    web = web[rowSums(web==0)!=0, colSums(web==0)!=0]

    # Get non zero positions.
    non_zero <- which(!web == 0)
    
    # Run rgraph 
    ans = .Call("bipartmod",
        row(web)[non_zero]-1L, col(web)[non_zero]-1L,
        web[cbind(row(web)[non_zero], col(web)[non_zero])])

    # Build the dataframe 
    df = data.frame(rownames(web),ans[[1]],ans[[2]],ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
} 
