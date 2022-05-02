#Adams-Moulton Fourth Order Method with Range-Kutta Method Input
sink(file="")
f<-function(x,y){
  return()
}

#Initialization of parameters for Range-Kutta Method
x<-
y<-
h<-
n<-

#Range-Kutta Method
RangeKutta<-function(x,y,h,n){
  yValues<-vector(mode="numeric",n)
  df<-NULL
  for(i in 1:n){
    k1<-f(x,y)
    k2<-f(x+(0.5*h),y+(0.5*h*k1))
    k3<-f(x+(0.5*h),y+(0.5*h*k2))
    k4<-f(x+h,y+h*k3)
    y<-y+(h/6)*(k1+2*k2+2*k3+k4)
    x<-x+h
    yValues[i]<-y
  }
  return(yValues)
}

AMValues<-RangeKutta(x,y,h,n)
AMValues

#Adams-Moulton Fourth Order Method
AdamsMoultonFourthOrder<-function(x0,yNMinus2,yNMinus1,yN,yNPlus1,h,n){
  x1<-x0+(n-3)*h
  x2<-x0+(n-2)*h
  x3<-x0+(n-1)*h
  x4<-x0+n*h
  k1<-f(x1,yNMinus2)
  k2<-f(x2,yNMinus1)
  k3<-f(x3,yN)
  k4=f(x4,yNPlus1)
  yNPlus1<-yN +(h/24)*((9*k4)+(19*k3)-(5*k2)+k1)
  cat("The value of Y+1 at x=", x4, "is:")
  return(yNPlus1)
}
AdamsMoultonFourthOrder(x,AMValues[n-3],AMValues[n-2],AMValues[n-1],AMValues[n],h,n)
sink()
