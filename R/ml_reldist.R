#' Relative Distance in Finite Metric Space
#' 
#' @references 
#' \insertRef{chang_mathematical_2016}{DAS}
#' 
#' @export
reldist <- function(x, label=NULL){
  ##################################################
  # Check Input and Transform
  x = check_input(x)
  n = nrow(x)
  
  ##################################################
  # Compute : relative distance between points
  output.mat <- array(0,c(n,n))
  for (i in 1:n){
    for (j in 1:n){
      output.mat[i,j] = x[i,j] - base::mean(as.vector(x[i,]))
    }
  }
  
  ##################################################
  # Compute : average relative distance
  output.vec = as.vector(base::colMeans(output.mat))
  
  ##################################################
  # Compute : label, if provided
  if ((length(label)==0)&&(is.null(label))){
    # return
    output = list()
    output$points  = output.mat
    output$average = output.vec
    return(output)
  } else {
    # check the input
    if (!(is.vector(label)&&(length(label)==n))){
      stop("* reldist : 'label' should be a vector of class membership.")
    }
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
        output.set[i,j] = sum(output.mat[idSi,idSj])/(length(idSi)*length(idSj))
      }
    }
    rownames(output.set) = llevel
    colnames(output.set) = llevel
    
    # return
    output = list()
    output$points  = output.mat
    output$average = output.vec
    output$sets    = output.set
    return(output)
  }
}

# # a 2-dimensional example
# x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
# colnames(x) <- c("x", "y")
# (cl <- kmeans(x, 2))
# clcluster = cl$cluster
# reloutput = reldist(dist(x), label=clcluster)
  