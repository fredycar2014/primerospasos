primerospasos
=============

Código base de aprendizaje de R

    ##########################################
    #######----PRIMEROS PASOS EN R-----#######
    ####--FREDY GONZALO CARVAJAL SERRANO--####
    ##########################################

    #####   GENERALIDADES   #####
    
# CITANDO R
citation()

# SOLICITAR AYUDA
?read.table
help.search("data input")

# ENCONTRAR INFORMACIÓN
find("abline")
apropos("lm")
    
# INSTALANDO PAQUETES
install.packages("graphics")

# EJEMPLOS DE AYUDA
example(lm)    
  
# DEMOS DE AYUDA
demo(graphics)

# CONTENIDO DE PAQUETES EN R
library(help=spatial)

# EDITOR DE DATOS
library(MASS)
attach(bacteria)
fix(bacteria) 

# VARIABLES QUE SE HAN CREADO EN LA SESIÓN
objects()  

# VARIABLES ACTIVAS (ATTACH)
search()
   
# REMOVER VARIABLES CREADAS    
rm(x,y,z)
detach(bacteria)
    
# PARA DESHACERSE DE TODO!!
rm(list=ls())   
    
# MATEMÁTICAS
    
x<-10
    log(x) #log to base e of x
    exp(x) #antilog of x (ex)
    log(x,2) #log to base n of x
    log10(x) #log to base 10 of x
    sqrt(x) #square root of x
    factorial(x) #x! = x × (x - 1) × (x - 2)×· · ·×3 × 2
    choose(25,x) #binomial coefficients n!/(x! (n  x)!)
    gamma(x) #(x), for real x (x1)!, for integer x
    lgamma(x) #natural log of (x)
    floor(x) #greatest integer less than x
    ceiling(x) #smallest integer greater than x
    trunc(x) #closest integer to x between x and 0, e.g. trunc(1.5) = 1, trunc(1.5) = 1;
    round(x, digits=2) #round the value of x to an integer
    signif(x, digits=6) #give x to 6 digits in scientific notation
    runif(x) #generates n random numbers between 0 and 1 from a uniform distribution
    cos(x) #cosine of x in radians
    sin(x) #sine of x in radians
    tan(x) #tangent of x in radians
    acos(1/x) 
    asin(x^(-1))
    atan(x)
    acosh(x); asinh(1/x); atanh(1/x)
    abs(-x) #the absolute value of x, ignoring the minus sign if there is one
    sin(tan(x))
    x^sin(x)
    exp(-sin(x))   
pi
cos(pi/2)
9 %/% 4
9 %% 4
15421 %% 7 == 0
15421/7
    
# VARIABLES
    data <- read.table("c://data/datas/daphnia.txt", header=T)
    attach(data)
    head(data)
    is.factor(Water)
    is.numeric(Growth.rate)
    is.numeric(Detergent)
    levels(Detergent)
    nlevels(Detergent)
    length(levels(Detergent))
    as.vector(unclass(Daphnia))
    as.integer(unclass(Growth.rate))
y <- 10
x <- sqrt(y)
x*x == 10
x*x-10   

# SECUENCIAS
    
    0:10
    20:2
    seq(20,30,0.5)
    seq(30,20,-5)
    N <- c(55,76,92,103,84,88,121,91,65,77,99)
    seq(from=0.04,by=0.01,length=11)
    seq(0.04,by=0.01,along=N)
    rep(9,5)
    rep(1:4,3)
    rep(1:4, each=3)
    rep(1:4, each=3, times=3)
    rep(1:4,1:4)
    
# GENERAR NIVELES DE FACTORES
    
    gl(4,3)
    gl(4,3,24)
    
    Temp <- gl(2, 2, 24, labels = c("Low", "High"))
    Soft <- gl(3, 8, 24, labels = c("Hard","Medium","Soft"))
    M.user <- gl(2, 4, 24, labels = c("N", "Y"))
    Brand <- gl(2, 1, 24, labels = c("X", "M"))
    data.frame(Temp,Soft,M.user,Brand)
    
    is.finite(10)
    is.infinite(10)
    is.infinite(Inf)
   
# TRABAJANDO CON SUBÍNDICES LÓGICOS
    
    x <- 0:10
    sum(x)
    sum(x<5)
    sum(x[x<5])
    
    y <- c(8,3,5,7,6,6,8,9,2,3,9,4,10,4,11)
    sort(y)
    rev(sort(y))
    
    x <- c(2,3,4,1,5,8,2,3,7,5,7)
    which.max(x)
    which.min(x)
    
# FUNCIONES VECTORIALES
    
    x <- c(8,3,5,7,6,6,8,9,2,3,9,4,10,4,11)
    mean(x)
    max(x)
    min(x)
    sum(x)
    mean(x)
    median(x)
    range(x)
    var(x)
    cor(x,y)
    sort(x)
    rank(x)
    order(x)
    quantile(x)
    cumsum(x)
    cumprod(x)
    cummax(x)
    cummin(x)
    pmax(x,y,z)
    pmin(x,y,z)
    cumprod(1:5)
    cumsum(1:5)
    
    counts <- rnbinom(10000,mu=0.92,size=1.1)
    counts[1:30]
    table(counts)
    
    data<-read.table("c:/data/temperatures.txt",header=T)
    attach(data)
    names(data)
    tapply(temperature,month,mean)
    tapply(temperature,month,var)
    tapply(temperature,month,min)
    tapply(temperature,month,function(x) sqrt(var(x)/length(x)))
    tapply(temperature,list(yr,month),mean)[,1:6]
    tapply(temperature,yr,mean,na.rm=TRUE)
    tapply(temperature,yr,mean,trim=0.2)
    
    data<-read.table("c:/data/pHDaphnia.txt",header=T)
    names(data)
    aggregate(Growth.rate~Water,data,mean)
    aggregate(Growth.rate~Water+Detergent,data,mean)
    aggregate(cbind(pH,Growth.rate)~Water+Detergent,data,mean)

# DATAFRAMES
    
    pacienteID <- c(1, 2, 3, 4)
    EDAD <- c(25, 34, 28, 52)
    DIABETES <- c("T-1", "T-2", "T-1", "T-1")
    ESTADO <- c("Pobre", "Mejorando", "Excelente", "Pobre")
    pacientedata <- data.frame(pacienteID, EDAD, DIABETES, ESTADO)
    pacientedata

    worms <- read.csv("c:/data/worms.csv", sep=";", header=TRUE)
    worms
    names(worms)
    attach(worms)
    worms[7]
    worms[,1:3]
    worms[5:15,]
    worms[Area>3 & Slope <3,]
    worms[Area>3 & Slope <3 & Worm.density >4,]
    worms[order(Area),]
    worms[order(Area),c(2,3,5,7)]
    worms[rev(order(worms[,5])),c(5,7)]
    summary(worms)
    with(worms,tapply(Worm.density,Vegetation,mean))
    aggregate(worms[,c(2,3,5,7)],list(Vegetation),mean)
    aggregate(worms[,c(2,3,5,7)],list(Comunidad=Vegetation),mean)
    
    dataz <- read.table("c:/data/das.txt", sep=";", header=TRUE)
    attach(dataz)
    head(dataz)
    plot(dataz)
    which(dataz > 10)
    
    data4 <- read.csv("c:/data/weather..csv", sep=";", header = TRUE)
    attach(data4)
    head(data4)
    plot(factor(month),upper)
    
    
    
# LISTAS
    
    g <- "Mi Primera Lista"
    h <- c(25, 26, 18, 39)
    j <- matrix(1:10, nrow=5)
    k <- c("uno", "dos", "tres")
    mylist <- list(title=g, ages=h, j, k)
    mylist

# ENTRADA DE DATOS
    
    midata <- data.frame(EDAD=numeric(0), GENERO=character(0), PESO=numeric(0))
    midata <- edit(midata)
    
# LA FUNCIÓN SAMPLE
    
    y <- c(8,3,5,7,6,6,8,9,2,3,9,4,10,4,11)
    sample(y)
    sample(y,5)

# LOOPS Y REPAEATS
    
    for (i in 1:5) print(i^2)
    fibonacci <- function(n) {
      +       a <- 1
      +       b <- 0
      +       while(n>0)
        +       {swap <- a
                 +        a <- a+b
                 +        b <- swap
                 +        n <- n-1 }
      +       b }
    sapply(1:10,fibonacci)

# FECHAS Y TIEMPO
    
    Sys.time()
    as.numeric(Sys.time())
    y <- strptime("01/06/1965",format="%d/%m/%Y")
    weekdays(y)
    
    #####   ENTRADA DE DATOS   #####
    
    x <-scan()
    
    getwd()
    
    data <- read.table("c://data/phDaphnia.txt", header=T)
    data
    data <- read.delim("c:/data/phDaphnia.txt")
    data

    install.packages("RODBC")
    library(RODBC)
    channel <- odbcConnect("northwind")
    
# GRÁFICOS
    
    mtcars
    attach(mtcars)
    summary(mpg)
    summary(mtcars)
    plot(wt, mpg)
    abline(lm(mpg~wt))
    title("Regresión de MPG sobre Weight")
    detach(mtcars)
    
    jpeg("mygraph.jpg")
    attach(mtcars)
    plot(wt, mpg)
    abline(lm(mpg~wt))
    title("Regresión de MPG sobre Peso")
    detach(mtcars)
    dev.off()
    
    dosis <- c(20, 30, 40, 45, 60)
    drogaA <- c(16, 20, 27, 40, 60)
    drogaB <- c(15, 18, 25, 31, 40)
    plot(dosis, drogaA, type="b")
    
    opar <- par(no.readonly=TRUE)
    par(lty=1, pch=20)
    plot(dosis, drogaA, lwd=3, cex=2, type="b")
    title("Mi Gráfico")
    par(opar)
    
    dose <- c(20, 30, 40, 45, 60)
    drugA <- c(16, 20, 27, 40, 60)
    drugB <- c(15, 18, 25, 31, 40)
    opar <- par(no.readonly=TRUE)
    par(pin=c(2, 3))
    par(lwd=2, cex=1.5)
    par(cex.axis=.75, font.axis=3)
    plot(dose, drugA, type="b", pch=19, lty=2, col="red")
    plot(dose, drugB, type="b", pch=23, lty=6, col="blue", bg="green")
    par(opar)
    
    plot(dosis, drogaA, type="b",
         col="red", lty=1, pch=2, lwd=1,
         main="Pruebas Clínicas para Droga A",
         sub="Este es un Dato Hipotético",
         xlab="Dosis", ylab="Drga Respuesta",
         xlim=c(0, 60), ylim=c(0, 70))
  
    x <- c(1:10)
    y <- x
    z <- 10/x
    opar <- par(no.readonly=TRUE)
    par(mar=c(5, 4, 4, 8) + 0.1)
    plot(x, y, type="b",
         pch=21, col="red",
         yaxt="n", lty=3, ann=FALSE)
    lines(x, z, type="b", pch=22, col="blue", lty=2)
    axis(2, at=x, labels=x, col.axis="red", las=2)
    axis(4, at=z, labels=round(z, digits=2),
         col.axis="blue", las=2, cex.axis=0.7, tck=-.01)
    mtext("y=1/x", side=4, line=3, cex.lab=1, las=2, col="blue")
    title("An Example of Creative Axes",
          xlab="X values",
          ylab="Y=X")
    par(opar)
    
    
    
    data1 <- read.table("c:/data/scatter1.txt",header=T)
    attach(data1)
    names(data1)
    plot(x1,y1,col="red")
    plot(x1,y1,col="red",xlab="Variable Explicatoria",ylab="Variable Respuesta")
    abline(lm(y1~x1))
    data2 <- read.table("c:/data/scatter2.txt",header=T)
    attach(data2)
    names(data2)
    abline(lm(y2~x2))
    plot(c(x1,x2),c(y1,y2),xlab="Variable Explicatoria",
         ylab="Variable Respuesta",type="n")
    points(x1,y1,col="red")
    points(x2,y2,col="blue",pch=16)
    abline(lm(y1~x1))
    abline(lm(y2~x2))
    range(c(x1,x2))
    range(c(y1,y2))
    plot(c(x1,x2),c(y1,y2),xlim=c(0,100),ylim=c(0,70),
      xlab="Variable Explicatoria",ylab="Variable Respuesta",type="n")
    points(x1,y1,col="red")
    points(x2,y2,col="blue",pch=16)
    abline(lm(y1~x1))
    abline(lm(y2~x2))
    
    # COLOR DE SÍMBOLOS
    
    plot(0:9,0:9,pch=16,type="n",xaxt="n",yaxt="n",ylab="col",xlab="bg")
    axis(1,at=1:8)
    axis(2,at=1:8)
    for (i in 1:8) points(1:8,rep(i,8),pch=c(21,22,24),bg=1:8,col=i)
    
    # ADICIONANDO TEXTO A GRÁFICOS
    
    map.places <- read.csv("c:/data/map.places.csv",header=T)
    attach(map.places)
    names(map.places)
    
    map.data <- read.csv("c:/data/bowens.csv",header=T)
    attach(map.data)
    names(map.data)
    
    nn <- ifelse(north<60,north+100,north)
    
    windows(9,7)
    plot(c(20,100),c(60,110),type="n",xlab="",ylab="",xaxt="n", yaxt="n")
    for (i in 1:length(wanted)){
      ii <- which(place == as.character(wanted[i]))
      text(east[ii], nn[ii], as.character(place[ii]), cex = 0.6) }
    
    
    library(qcc)
    #make 2 plots in 1 figure
    par(mfrow=c(1,2))
    
    #points have base value of 10 w/ normally distributed error
    lugnuts <- rep(10, 100) + rnorm(100, mean=0, sd=0.5)
    qcc(lugnuts, type="xbar.one", center=10, add.stats=FALSE,
        title="1st Batch", xlab="i-th lugnut produced")
    
    #first 90 points have base value of 10 w/ normally distributed error,
    #last 10 points have base value of 11 w/ normally distributed error
    lugnuts <- c(rep(10, 90), rep(11, 10)) + rnorm(100, mean=0, sd=0.5)
    qcc(lugnuts, type="xbar.one", center=10, add.stats=FALSE,
        title="2nd Batch", xlab="i-th lugnut produced")
    
    #example using holdout/test sets
    lugnuts <- rep(10, 100) + rnorm(100, mean=0, sd=0.5)
    qcc(lugnuts, newdata=rep(11, 10) + rnorm(10, mean=0, sd=0.5),
        type="xbar.one", center=10, add.stats=FALSE, title="2nd Batch", xlab="i-th lugnut produced")
    
    data(package="MSQC")
    example(MSQC)
    
    library("qcc")
    set.seed(20)
    x <- round(rnorm(120,20,2),2)
    length <- matrix(x, ncol = 4, byrow = TRUE)
    par(mfrow = c(1,2))
    qcc(length, type = "xbar", std.dev = "RMSDF"); qcc(length, type = "R")
    qcc(length, type = "R")

    library(ggplot2)
    qplot(mtcars$wt, mtcars$mpg)
    qplot(wt, mpg, data=mtcars)
    # This is equivalent to:
    ggplot(mtcars, aes(x=wt, y=mpg)) + geom_point()
    plot(pressure$temperature, pressure$pressure, type="l")
    plot(pressure$temperature, pressure$pressure, type="l")
    points(pressure$temperature, pressure$pressure)
    lines(pressure$temperature, pressure$pressure/2, col="red")
    points(pressure$temperature, pressure$pressure/2, col="red")
    
    install.packages("ggplot2")
    install.packages("gcookbook")
