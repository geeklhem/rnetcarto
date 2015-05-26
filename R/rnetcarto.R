#' @useDynLib rnetcarto, CBipartmod, CNetcarto
NULL
#' Compute modularity and modularity roles for graphs using simulated
#' annealing
#'
#' @docType package
#' @name rnetcarto


#' @title Computes modularity and modularity roles from a network.
#'
#' @param web network either as a square adjacency matrix or a list
#' describing E interactions a->b: the first (resp. second) element is
#' the vector of the labels of a (resp. b), the third (optional) is
#' the vector of interaction weights.
#' @param seed Seed for the random number generator: Must be a
#' positive integer.
#' @param iterfac At each temperature of the simulated annealing
#' (SA), the program performs fN^2 individual-node updates (involving the
#' movement of a single node from one module to another) and fN
#' collective updates (involving the merging of two modules and the split
#' of a module). The number "f" is the iteration factor.
#' @param symmetric If TRUE all edges a->b are copied b->a.
#' @param coolingfac Temperature cooling factor.
#' @param auto_link If TRUE allows self looping edges a->a
#' @param add_weight If TRUE weights are summed if the edge already exist.
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export
netcarto <- function(web,
                     seed=as.integer(floor(runif(1, 1,100000001))),
                     iterfac=1.0,
                     symmetric=TRUE,
                     coolingfac=0.995,
                     auto_link=FALSE,
                     add_weight=FALSE)
{

    # Read matrix of list input..
    if(is.matrix(web)){
        # Sanity checks...
        if(ncol(web) != nrow(web)){
            stop("Input matrix web must be square")
        }
        if (all(rownames(web) != colnames(web))){
            warning("Columns and row names are not matching, are you sure this is an adjacency matrix ?")
        }
        if (is.null(rownames(web))){
            rownames(web) = 1:nrow(web)
        }
        # Removing empty columns and lines.
        web = web[rowSums(web==0)!=ncol(web), colSums(web==0)!=nrow(web)]

        # Get non zero positions.
        non_zero <- which(!web == 0)

        # Parameters
        nodes1 = row(web)[non_zero]-1L
        nodes2 = col(web)[non_zero]-1L
        weights = web[cbind(row(web)[non_zero], col(web)[non_zero])]
        names = rownames(web)

    } else if (is.list(web)){

        E = length(web[[1]]) # Number of edges

        # Read the weigth if they are supplied
        if (length(web)==3){
            weights = web[[3]]
        } else if (length(web)==2){
            weights = numeric(E) + 1
        } else{
            stop("Input list should be of length 2 (unweighted edges) or 3 (weighted edges)");
        }

        # Convert the species names to integer

        fct = factor(c(web[[1]], web[[2]]))
        idx = as.integer(fct) - 1L
        names = levels(fct)
        nodes1 = idx[1:E]
        nodes2 = idx[(E+1):(2*E)]
    }


    # Call rgraphlib
    ans <- .Call("CNetcarto",as.integer(nodes1), as.integer(nodes2) , as.numeric(weights),
                 as.numeric(coolingfac), as.integer(seed),
                 as.numeric(iterfac), as.integer(symmetric), as.integer(auto_link), as.integer(add_weight))

    # Build the dataframe
    df = data.frame(names, ans[[1]], ans[[2]], ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
}


#' @title Computes modularity and modularity roles from a bipartite network.
#'
#' @param web interactions between two set of species either as an
#' adjacency matrix or a list describing E interactions, the first
#' element is the vector of species in the first set, the second the
#' vector of the corresponding species in the second set and the third
#' (optional) is the vector of interaction weights.
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
#' @param add_weight If TRUE weights are summed if the edge already
#' exist (only meaningful if the list input format is used).
#' @return A list. The first element is a dataframe with the name,
#' module, z-score, and participation coefficient for each row of the
#' input matrix. The second element is the modularity of this
#' partition.
#' @export
bipartmod <- function(web,
                      seed=as.integer(floor(runif(1, 1,100000001))),
                      iterfac=1.0,
                      coolingfac=.950,
                      degree_based=FALSE,
                      weighted=TRUE,
                      add_weight=TRUE)
{
    # Read the input matrix or list.
    if(is.matrix(web)){
        if (is.null(rownames(web))){
            rownames(web) = 1:nrow(web)
        }

        # Removing empty columns and lines.
        web = web[rowSums(web==0)!=ncol(web), colSums(web==0)!=nrow(web)]

        # Get non zero positions.
        non_zero <- which(!web == 0)

        # Parameters
        nodes1 = row(web)[non_zero]-1L
        nodes2 = col(web)[non_zero]-1L
        weights = web[cbind(row(web)[non_zero], col(web)[non_zero])]
        names = rownames(web)
    } else if (is.list(web)){
        if (length(web[[1]]) != length(web[[2]])){
            stop("Bad labels number: all elements of the input list should have the same length.")}

        # Read the weigth if they are supplied
        if (length(web)==3){
            if (length(web[[3]])==length(web[[1]])){
                weights = web[[3]]
            }else{
                stop("Bad weight number: all elements of the input list should have the same length.")
            }
        } else if (length(web)==2){
            weights = numeric(length(web[[1]])) + 1
        } else{
            stop("Input list web should be of length 2 (unweighted edges) or 3 (weighted edges).");
        }

        # Convert the species names to integer
        fct_par1= factor(web[[1]])
        fct_par2= factor(web[[2]])
        nodes1 = as.integer(fct_par1) - 1L
        nodes2 = as.integer(fct_par2) - 1L
        names = levels(fct_par1)
    } else{
        stop("web should be a list or a matrix.")
    }

                                        # Run rgraph
    ans = .Call("CBipartmod",
        as.integer(nodes1), as.integer(nodes2) ,as.numeric(weights),
        as.integer(seed), as.numeric(iterfac), as.numeric(coolingfac),
        as.integer(degree_based), as.integer(weighted),
        as.integer(add_weight))

    
    # Build the dataframe
    df = data.frame(names,ans[[1]],ans[[2]],ans[[3]])
    names(df) <- c("name","module","z-score","participation")
    return(list(df,ans[[4]]))
}
