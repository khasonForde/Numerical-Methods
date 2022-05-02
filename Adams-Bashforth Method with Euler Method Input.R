#Adams-Bashforth Method with Euler Method Input-Multistep Methods
#Function in Differential Equation
f<-function(x,y){
  return(1-x+4*y)
}

#Initialization of parameters for Euler Method
x<-0
y<-1
h<-0.1
n<-1

#Euler Method
Euler<-function(x,y,h,n){
  yValues<-vector(mode="numeric",(n))
  yValues[1]=y
  df<-NULL
  for(i in 1:(n)){
    k1<-f(x,y)
    y<-y+h*k1
    x<-x+h
    yValues[i]<-y
  }
  return(yValues[1])
}

#Parameters for Adams-Bashforth Method
ABy0<-y
ABy1<-ImprovedEuler(x,y,h,n)

#Adams-Bashforth Method
AdamsBashforthSecondOrder<-function(x0,y0,y1,h,n){
  df<-NULL
  k1<-f(x0,y0)
  k2<-f(x0+h,y1)
  y2<-y1+((3/2)*h*k2)-((1/2)*h*k1)
  k3<-f(x+(n+1)*h,y)
  x1<-x0+h
  df<-rbind(df, data.frame(X=x1+h,Y=y2))
  AdamsBashforthMethod<-df
  cat("The value of f at x-1 and y-1 is", k1, "\n")
  cat("The value of f at x and y", k2, "\n")
  cat("The value of f at x+1 and y+1", k3, "\n")
  cat("Adams-Bashforth Second Order Method Approximation of y+1 if x-1=", x0," with y-1=", ABy0, "and x=", x1,"with y=", ABy1, ":\n")
  return(AdamsBashforthMethod)
}

AdamsBashforthSecondOrder(x,ABy0,ABy1,h,1)