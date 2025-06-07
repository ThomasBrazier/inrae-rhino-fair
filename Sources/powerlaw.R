dpowerlaw <- function (x, xmin, alpha, log = FALSE) 
{
  if (log) {
    pdf = log(alpha - 1) - log(xmin) - alpha * (log(x/xmin))
    pdf[x < xmin] = -Inf
  }
  else {
    pdf = (alpha - 1)/xmin * (x/xmin)^(-alpha)
    pdf[x < xmin] = 0
  }
  pdf
}

ppowerlaw <- function (q, xmin, alpha, lower.tail = FALSE) 
{
  cdf = 1 - (q/xmin)^(-alpha + 1)
  if (!lower.tail) 
    cdf = 1 - cdf
  cdf[q < round(xmin)] = 0
  cdf
}

qpowerlaw <- function(p,xmin,alpha){
  alpha*p^(1/(1-xmin))
}

