#Digital computing integrals
quad_Integral <- function(f,a,b)
{
  l<-(b-a)/2
  mid<-(b+a)/2
  return(l*( f(l/sqrt(3) +mid) +
               f(-l/sqrt(3) +mid) ))
}

#returns ei for n elements on given interval [a,b]
base_Function <- function(a,b, n,i, derivative)
{
  interval=(b-a)/n
  
  
  if (!derivative)
  {
    if (i==1)
    {
      e <- function(x)
      {
        if (x>=a && x<=a+interval)
          return((a+interval-x)/interval)
        return(0)
      }
    }
    else if (i>1 && i<n)
    {
      i=i-1
      e <- function(x)
      {
        if (x>=a + (i-1)*interval && x<=a + i*interval)
          return((x- a-(i-1)*interval)/interval)
        else if (x>=a + i*interval && x<=a + (i+1)*interval)
          return((a+(i+1)*interval -x)/interval)
        return(0)
      }
    }
    else if (i==n)
    {
      e <- function(x) 0
    }
  }
  else
  {
    if (i==1)
    {
      e <- function(x)
      {
        if (x>=a && x<=a+interval)
          return(-1/interval)
        return(0)
      }
    }
    else if (i>1 && i<n)
    {
      i=i-1
      e <- function(x)
      {
        if (x>=a + (i-1)*interval && x<=a + i*interval)
          return(1/interval)
        else if (x>=a + i*interval && x<=a + (i+1)*interval)
          return(-1/interval)
        return(0)
      }
    }
    else if (i==n)
    {
      e<- function(x) 0
    }
  }
  return(e)
}

E <- function(x)
{
  if (x>=0 && x<=1) 
    return(3)
  if (x>1&&x<=2)
    return(5)
  return(0)
}

construct_Matrix <- function(n,a,b)
{
  m<-matrix(nrow=n-1,ncol=n-1)
  for (i in 1:(n-1))
  {
    func=function(x) E(x)*base_Function(a,b,n,i,TRUE)(x)*base_Function(a,b,n,1,TRUE)(x)
    B=quad_Integral(func,a,a+(b-a)/n)
    m[1,i]=B-E(0)*base_Function(a,b,n,1,FALSE)(0)*base_Function(a,b,n,i,FALSE)(0)
    m[i,1]=m[1,i]
  }

  for (i in 2:(n-1))
  {
    for (j in 2:(n-1))
    {
      func=function(x) E(x)*base_Function(a,b,n,j,TRUE)(x)*base_Function(a,b,n,i,TRUE)(x)
      B=quad_Integral(func,a+(i-2)*(b-a)/n,a+(i)*(b-a)/n)
      m[i,j] <- B-E(0)*base_Function(a,b,n,j,FALSE)(0)*base_Function(a,b,n,i,FALSE)(0)
    }
  }
  return(m)
}

constructVector <- function(n,a,b)
{
  v=matrix(nrow=n-1,ncol=1)
  for(i in 1:(n-1))
    v[i,1]=-10*E(0)*base_Function(a,b,n,i,FALSE)(0)
  return(v)
}

constructFunction <- function(n,a,b)
{
  B=construct_Matrix(n,a,b)
  v=constructVector(n,a,b)
  solution=solve(B,v)
  print(solution)


  answer=function(x)
  {
    s=0
    for (i in 1:(n-1))
      s=s+base_Function(a,b,n,i,FALSE)(x)*solution[i,1]
    return(s)
  }
  return(answer)
}

f=constructFunction(100,0,2)




