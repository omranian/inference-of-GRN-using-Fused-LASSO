
cubic.spline.interpol <- function(data,domain,x)
{
  res <- matrix(1,nrow(data),length(x))
  for (i in 1:nrow(data))
  {
    res[i,] <- splint(domain,as.numeric(data[i, ]),x)
  }
  return(res)
}