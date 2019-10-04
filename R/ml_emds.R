#' Euclified Multidimensional Scaling
#' 
#' strategy 1 : transitive closure of the triangle inequality (labdsv)
#' strategy 2 : Non-Euclidean or Non-metric Measures Can Be Informative; adding positive numbers to all off-diagonal entries
#' 
#' 
#' @export
emds <- function(x, ndim=2, method=c("closure","gram")){
  ##################################################3
  # Check Input and Transform
  ndim = round(ndim)
  k  = as.integer(ndim)
  x  = check_input(x)
  n  = nrow(x)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* DAS::emds - 'ndim' should be an integer in [1,nrow(x)). ")
  }

  method = match.arg(method) 
  mydim  = round(ndim)
  
  ##################################################3
  # Branching
  if (nef(x) < 100*.Machine$double.eps){ # if Euclidean, okay
    output = cmds(x, ndim=mydim)
  } else { # if not Euclidean
    if (method=="closure"){ # strategy 1 : transitive closure of the triangle inequality
      xnew   = as.matrix(labdsv::euclidify(stats::as.dist(x))) # well it seems to work well..
      output = cmds(xnew, ndim = mydim)
    } else {                # strategy 2 : add positive numbers to all off-diagonal entries
      gamma0 = emds_gamma0(x)
      ggrid  = seq(from=min(0.001, gamma0/1000), to=(gamma0*0.999), length.out=20) # just try 20 cases
      vgrid  = rep(0,20)
      for (i in 1:20){
        xtmp = x + ggrid[i]
        diag(xtmp) = rep(0,nrow(xtmp))
        vgrid[i] = nef(xtmp)
      }
      idopts = which.min(vgrid)
      if (length(idopts)>1){ # if multiple, use the first one.
        idopts = idopts[1]
      }
      optgamma   = ggrid[idopts]
      xnew       = x + optgamma
      diag(xnew) = rep(0,nrow(xnew))
      output     = cmds(xnew, ndim = mydim)
    } 
  }
  
  ##################################################3
  # Report 
  return(output)
}


# example -----------------------------------------------------------------
# library(labdsv)
# data(bryceveg) # returns a vegetation data.frame
# dis.bc <- as.matrix(dsvdis(bryceveg,'bray/curtis')) # calculate a Bray/Curtis
# 
# out.cmds <- cmds(dis.bc, ndim=2)$embed
# out.emds1 <- emds(dis.bc, ndim=2, method="closure")$embed
# out.emds2 <- emds(dis.bc, ndim=2, method="gram")$embed
# 
# par(mfrow=c(1,3),pty="s")
# plot(out.cmds, main="cmds")
# plot(out.emds1, main="emds::closure")
# plot(out.emds2, main="emds::gram")

#' 
#' dis.bc2 = dis.bc + 2
#' diag(dis.bc2) = 0
#' 
#' n  = nrow(dis.bc2)
#' D2 = dis.bc2^2
#' H  = diag(n)- (1/n)*outer(rep(1,n),rep(1,n))
#' J  = -0.5*(H%*%D2%*%H)
#' min(eigen(J)$values)
