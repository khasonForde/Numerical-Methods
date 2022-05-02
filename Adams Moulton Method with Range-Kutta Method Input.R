#Adams-Moulton Method with Range-Kutta Method Input-Multistep Methods
f<-function(x,y){
  return(1-x+4*y)
}

#Initialization of parameters for Range-Kutta Method
x<-0
y<-1
h<-0.1
n<-1

#Range-Kutta Method
RangeKutta<-function(x,y,h,n){
  yValues<-vector(mode="numeric",(n))
  yValues[1]=y
  df<-NULL
  for(i in 1:(n)){
    k1<-f(x,y)
    k2<-f(x+(0.5*h),y+(0.5*h*k1))
    k3<-f(x+(0.5*h),y+(0.5*h*k2))
    k4<-f(x+h,y+h*k3)
    y<-y+(h/6)*(k1+2*k2+2*k3+k4)
    x<-x+h
    yValues[i]<-y
  }
  return(yValues[1])
}

#Parameters for Adams-Moulton Method
AMy0<-y
AMy1<-RangeKutta(x,y,h,n)

#Adams-Moulton Second Order Method
AdamsMoultonSecondOrder<-function(x0,y0,y1,h,n){
  df<-NULL
  k1<-f(x0,y0)
  k2<-f(x0+h,y1)
  x1<-x0+h
  NewApproxy1<-y0+((1/2)*h*k1)-((1/2)*h*k2)
  df<-rbind(df, data.frame(X=x1,Y=NewApproxy1))
  AdamsMoultonMethod<-df
  cat("The value of f at x and y is", k1, "\n")
  cat("The value of f at x+1 and y+1 is", k2, "\n")
  cat("Adams-Moutlton Second Order Method Approximation for y+1 at x=", x0, "with y=", y0, "and x+1=", x1, "with y+1=", y1, ":\n" )
  return(AdamsMoultonMethod)
}
AdamsMoultonSecondOrder(x,AMy0,AMy1,h,1)