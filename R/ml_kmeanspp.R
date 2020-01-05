#' k-means++ seeding
#' 
#' 
#' @export
kmeanspp <- function(x, k=2){
  ##################################################3
  # Check Input and Transform
  x  = check_input(x)
  n  = nrow(x)
  K  = round(k)
  if (K >= n){
    stop("* kmeanspp : 'k' should be smaller than the number of observations.")
  }
  if (K < 2){
    stop("* kmeanspp : 'k' should be larger than 1.")
  }
  id.now = 1:n
  
  ##################################################3
  # Computation
  #   initialize
  id.center = base::sample(id.now, 1)
  id.now    = base::setdiff(id.now, id.center)
  #   iterate
  for (i in 1:(K-1)){
    # compute distance to the nearest
    tmpdmat = x[id.now, id.center]
    if (i==1){
      d2vec = as.vector(tmpdmat)^2
      d2vec = d2vec/base::sum(d2vec)
    } else {
      d2vec = as.vector(base::apply(tmpdmat, 1, base::min))^2
      d2vec = d2vec/base::sum(d2vec)
    }
    # sample one
    id.tmp = base::sample(id.now, 1, prob=d2vec)
    # update
    id.center = c(id.center, id.tmp)
    id.now    = base::setdiff(id.now, id.tmp)
  }
  #   let's compute label
  dmat    = x[,id.center]
  cluster = base::apply(dmat, 1, base::which.min)
  
  ##################################################
  # Return
  output = list()
  output$cluster = cluster
  return(output)
}

# # personal experiment
# library(mlbench)
# cassini = mlbench::mlbench.smiley(200)
# mydata  = cassini$x
# mydist  = stats::dist(mydata)
# mylabel = cassini$classes
# 
# kpp3 = kmeanspp(mydist, k=3)$cluster
# kpp4 = kmeanspp(mydist, k=4)$cluster
# kpp5 = kmeanspp(mydist, k=5)$cluster
# kpp6 = kmeanspp(mydist, k=6)$cluster
# kpp7 = kmeanspp(mydist, k=7)$cluster
# 
# x11()
# par(mfrow=c(2,3), pty="s")
# plot(mydata[,1], mydata[,2], type="p", col=mylabel, pch=19, main="true", xlab="", ylab="")
# plot(mydata[,1], mydata[,2], type="p", col=kpp3, pch=19, main="kmeanspp 3", xlab="", ylab="")
# plot(mydata[,1], mydata[,2], type="p", col=kpp4, pch=19, main="kmeanspp 4", xlab="", ylab="")
# plot(mydata[,1], mydata[,2], type="p", col=kpp5, pch=19, main="kmeanspp 5", xlab="", ylab="")
# plot(mydata[,1], mydata[,2], type="p", col=kpp6, pch=19, main="kmeanspp 6", xlab="", ylab="")
# plot(mydata[,1], mydata[,2], type="p", col=kpp7, pch=19, main="kmeanspp 7", xlab="", ylab="")
