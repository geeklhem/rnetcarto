#' @useDynLib rnetcarto, CBipartmod, CNetcarto
NULL
#' Compute modularity and modularity roles for graphs using simulated
#' annealing
#'
#' @docType package
#' @name rnetcarto


#' @title Computes modularity and modularity roles from a network.
#'
#' @param interactions a (2,S)-shaped list of named species
#' @param seed Seed for the random number generator: Must be a
#' positive integer.
#' @param iterfac At each temperature of the simulated annealing
#' (SA), the program performs fN^2 individual-node updates (involving the
#' movement of a single node from one module to another) and fN
#' collective updates (involving the merging of two modules and the split
#' of a module). The number "f" is the iteration factor.
#' @param symmetric If !=0 all edges a->b are copied b->a.
#' @param coolingfac Temperature cooling factor. 
#' @param auto_link If !=0 allows self looping edges a->a
#' @param add_weight If !=0 weights are summed if the edge already exist.
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export 
netcarto <- function(interactions,
                     seed=as.integer(floor(runif(1, 1,100000001))),
                     iterfac=1.0,
                     symmetric=TRUE,
                     coolingfac=0.995,
                     auto_link=FALSE,
                     add_weight=FALSE)
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
    ans <- .Call("netcarto", idx[1:E], idx[(E+1):(2*E)], weight,
                 coolingfac, seed,
                 iterfac, as.integer(symmetric), as.integer(auto_link), as.integer(add_weight))

    # Build the dataframe 
    df = data.frame(levels(fct), ans[[1]], ans[[2]], ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
}


#' @title Computes modularity and modularity roles from a bipartite network.
#'
#' @param web a matrix representing the interactions observed between
#' two groups of species (respectively rows and columns).
#' @param seed Seed for the random number generator: Must be a
#' positive integer.
#' @param iterfac At each temperature of the simulated annealing
#' (SA), the program performs fN^2 individual-node updates (involving the
#' movement of a single node from one module to another) and fN
#' collective updates (involving the merging of two modules and the split
#' of a module). The number "f" is the iteration factor.
#' @param coolingfac Temperature cooling factor.
#' @param degree_based If TRUE, use the degree based roles metrics otherwise
#' use the strength based metric. 
#' @param weighted If TRUE, use the weighted modularity definition. 
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export 
bipartmod <- function(web,
                      seed=as.integer(floor(runif(1, 1,100000001))),
                      iterfac=1.0,
                      coolingfac=.995,
                      degree_based=FALSE, weighted=TRUE)
{
    # Removing empty columns and lines. 
    web = web[rowSums(web==0)!=0, colSums(web==0)!=0]

    # Get non zero positions.
    non_zero <- which(!web == 0)
    
    # Run rgraph 
    ans = .Call("CBipartmod",
        row(web)[non_zero]-1L, col(web)[non_zero]-1L,
        web[cbind(row(web)[non_zero], col(web)[non_zero])],
        seed, iterfac, coolingfac, as.integer(degree_based), as.integer(weighted), 0)

    # Build the dataframe 
    df = data.frame(rownames(web),ans[[1]],ans[[2]],ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
} 




