#' Cohesion in Finite Metric Space
#' 
#' Modularity is included..
#' 
#' @references 
#' \insertRef{chang_mathematical_2016}{DAS}
#' 
#' @seealso \code{\link{reldist}}
#' 
#' @export
cohesion <- function(x, label=NULL){
  ##################################################
  # Case Branching
  if ((length(label)==0)&&(is.null(label))){ # just need point-level cohesion
    # use reldist
    relout = DAS::reldist(x, label=NULL)
    RDy  = relout$average
    RDxy = relout$points
    n    = length(RDy)
    
    # compute point-level cohesion
    output.mat = array(0,c(n,n))
    for (i in 1:n){
      for (j in 1:n){
        output.mat[i,j] = RDy[j] - RDxy[i,j]
      }
    }
    
    # return the output
    output = list()
    output$points = output.mat
    return(output)
  } else {
    # use reldist
    relout = DAS::reldist(x, label)
    RDy  = relout$average
    RDxy = relout$points
    RDSS = relout$sets
    n    = length(RDy)
    
    # compute point-level cohesion
    output.mat = array(0,c(n,n))
    for (i in 1:n){
      for (j in 1:n){
        output.mat[i,j] = RDy[j] - RDxy[i,j]
      }
    }
    
    # compute set-level cohesion
    flabel = as.factor(label)
    llevel = base::levels(flabel)
    label  = as.integer(flabel)
    ulabel = unique(label)
    nlabel = length(ulabel)
    output.set = array(0,c(nlabel,nlabel))
    for (i in 1:nlabel){
      idSi = which(label==ulabel[i])
      for (j in 1:nlabel){
        idSj = which(label==ulabel[j])
        output.set[i,j] = sum(output.mat[idSi,idSj])
      }
    }
    
    # compute modularity of a clustering
    Q = base::sum(base::diag(output.set))
    
    # compute normalized modularity
    nQ = 0
    for (i in 1:nlabel){
      idSi = which(label==ulabel[i])
      nQ = nQ + sum(output.mat[idSi,idSi])/(length(idSi))
    }
    
    # return the output
    output = list()
    output$points = output.mat
    output$sets   = output.set
    output$modularity = Q
    output$normalized = nQ  # normalized modularity
    return(output)
  }
}


# # a 2-dimensional example
# x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
# colnames(x) <- c("x", "y")
# (cl <- kmeans(x, 2))
# clcluster = cl$cluster
# cohoutput = cohesion(dist(x), label=clcluster)
