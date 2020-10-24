###############################
###############################
###############################
#####                     #####
##### Week 1:             #####
#####   13 Jul 2020       #####
#####                     #####
###############################
###############################
###############################

# Getting the package "maps"
install.packages("maps")

map("nz")
map("world")

# about to evaluate 2+2
2+2

3+3 # this is 3+3

# Continuation line

1*2*3*
  4*5*6*7*8*9*10

(1*2*3 
  *4*5*6*7*8*9*10)

factorial(10)

prod(1:10)
sum(1:10)

x <- prod(1:10)


# Division
5/7



print(sqrt(7), digits=22)

x <- sqrt(7)

x^2 - 7



25/7 # floating point division
25 %% 7 # the remainder
25 %/% 7 # the dividend

# Rounding
pi

round(pi,5)
round(pi,1)
round(pi,0)
round(pi,-1)
round(pi,-2)

floor(pi)
ceiling(pi)

is.double(pi)
is.logical(pi)
is.logical(is.logical(pi))

x <- 5
is.double(x)

x <- "Hello"
is.double(x)
is.character(x)

x <- 5
y <- 10
x + y

x/(x+y)

1/0
sqrt(-1)

hello <- "goodbye"
hello

hello <- "k"
hello

cc <- "what is going on?!876986798!!!!!"
cc
cc <- 'what is going on?!\n"876986798!!!!!'
cc

# Vectors and arrays

x <- c(1, 6, 3, 2)
x

x <- 1:100
x

x <- sqrt(1:100)
x

x <- c(4, 6, 8, 4, 33)

y <- c(1,1,11,1)

z <- c(x,y)
z

z <- c(x,x,x,y,x,y)
z

length(z)
min(z)
max(z)
sort(z)
range(z)

c(min(z), max(z))


2 + 2
x <- 2
y <- 2
x + y

x <- c(1,2,4,8,16)
y <- c(0,1,1,0,1)

x + y

x + 3
x + c(3,3,3,3,3)

x <- c(4,2,3,1)
y <- c(1,5)

x + y
x + c(1,5,1,5)


x <- c(4,2,3,1,7)
y <- c(1,5)
x + y

x = y   #this one is BAD
x
x <- y # much better

# comparison
3==6
3>6
3>=6

x <- 1:4
y <- sqrt(1:4)
x
y

x==y
x>y

sqrt(0.25)

x
x[1]
x[3]
x <- runif(10)
x

x[1]
x[2]
x[1:3]
x>0.4

###x[1,4]  no!
x[c(1,4,1,1,1,1)]




###############################
###############################
###############################
#####                     #####
##### Week 2:             #####
#####   20 Jul 2020       #####
#####                     #####
###############################
###############################
###############################

x <- c(1, 4, 2, 7, 5, 4)

x[1]
x[5]
x[10]
length(x)
x[length(x)]

rev(x)[1]
rev(x)[1][1]


x
x[2:4][3]
x[2:4][c(2,1)]

# sequences
1:10
1:length(x)
1:3
x[1:3]

x[3:1]

x
x[1:length(x)]
rev(x)
x[length(x):1]

i <- c(3,2,5)
x[i]

x
which.min(x)
which.max(x)

x[which.max(x)]
max(x)

nn <- c("Mary","Whetu","Kien","Mike")
scores <- c(90, 76, 55, 99)

max(scores)
which.max(scores)
nn[which.max(scores)]

nn[which.min(scores)]
min(scores)

x <- rnorm(10, 5, 1)
x

range(x)
idx <- c(which.min(x), which.max(x))
idx
x[idx]

# negative indices
x[4]
x[-4]  # all elements except 4th element

length(x)
length(x[-4])

max(x)
x <- x[-which.max(x)]
max(x)

x
x[1:3]
x[-(1:3)]

x[-1:3]
-1:3

-(1:3)
x[-(1:3)]


y <- c(1,54,3,5,55)
y

y[2] <- 0
y

y[1:3] <- 66
y

y 
y<60

y[y<60]

y[y<60] <- 66
y

y <- c(-6, 4, -4, 66, 7, 4, 3)
# logical vector
idx <- y<0
idx
y[idx]
y[idx] <- 0
y

# named elements
nn <- c("Mary","Whetu","Kien","Mike")
scores <- c(90, 76, 55, 99)

scores <- c("Mary"=90, "Whetu"=76, "Kien"=55, "Mike"=99)
scores

scores[4]
scores["Mike"]

names(scores)
scores[c("Mike","Whetu","Mike")]

postgrads <- c("Mary", "Mike")
scores[postgrads]
max(scores[postgrads])
min(scores[postgrads])

scores
names(scores)
scores[names(scores)]

sort(names(scores))
scores[sort(names(scores))]
sort(scores)
?sort

sort(scores)
sort(scores, decreasing=TRUE)

rev(sort(scores))

# Sequences

100:120
98:76

seq(from=100, to=120)
?seq

seq(from=100, to=120, by=7)

# repetition


rep(5, 60)
rep("A", 7)
rep(nn, 2)
rep(1:10,6)

rep(1:10, each=6)

rep(1:10, 10:1)

nn
rep(nn, c(5,4,4,3))

# Arrays

x <- (1:12)^2
x

length(x)

a <- array(x, dim=c(3,4))
a

# a matrix is a 2D array
m <- matrix(x, nrow=3)
m
?matrix
b <- array(x, dim=c(3,2,2))
b

# subsets of an array
a

# individual element
a[2,4]
# full row
a[3,]
# FUll column
a[,3]

a[,-3]
a[,-c(1,4)]
a[,c(1,4)]

a[3,c(1,4)]

a[,c(1,4)]
a[,c(1,4)][3,]

array(1:12, dim=c(3,8))

array(1:12, dim=c(3,7))

a
t(a) # transpose

x
x[c(1,4)]

a
dimnames(a) <- list(c("Mike","Julie","Susan"), 
                    c("MATH132","STAT193","DATA101","BIOL113"))
a

a[1,]
a["Mike",]
a[,c("DATA101","STAT193")]


a <- array(round(runif(12, 45,99)), dim=c(3,4))
dimnames(a) <- list(c("Mike","Julie","Susan"), 
                    c("MATH132","STAT193","DATA101","BIOL113"))
a

max(a)
max(a[,"STAT193"])
max(a[,"MATH132"])

min(a["Julie",])

min(a)

apply(a, 1, min)
apply(a, 2, min)

?apply


# Matrix operations

a <- array(c(1,2,3,2,4,3), dim=c(3,2))
a

b <- array(c(2, 3, 5, 22, 3, 2, 3), dim=c(3,2))
b

a*b

a
b
# element wise multiplication
a*b

# matrix multiplication
c <- t(b)
a
c
a %*% c

# matrix inverse
dd <- a%*%c
solve(dd)

a
b

# binding columnwise
cbind(a,b)
# binding rowwise
rbind(a,b)

cbind(a, cbind( 1:3, 4:6 ))
cbind(a, 1:3, 4:6)

rbind(1:2, 1:4)

# More on logical vectors
a <- c(1, 3, 4, -1)

a<0

a[a<0] <- 7
a

a <- c(1, 3, 4, -1, -4, -5)
b <- a[a<0]
b

a[a<0] <- c(55,66,77)
a

a[a<0]

a <- c(1, 3, 4, -1, -4, -5)
a

# two ways to count the number of negative elements
length(a[a<0])
sum(a<0)

any(a<0)
all(a<0)

all(a<10)

a
a==4  # compares each element to the value 4
a=4   # assigns a the value 4
a <- 4

a <- c(1, 3, 4, -1, -4, -5)
b <- c(10, 20, 40, 30, 22, 0)

a>0
b
b[a>0]



smat <- array(round(runif(12, 45,99)), dim=c(3,4))
dimnames(smat) <- list(c("Mike","Julie","Susan"), 
                       c("MATH132","STAT193","DATA101","BIOL113"))
smat

s132 <- smat[,"MATH132"]
s101 <- smat[,"DATA101"]

s132
s132>57
s101[s132>57]

rownames(smat)
colnames(smat)

rownames(smat)[smat[,"MATH132"]>57]

s132 <- smat[,"MATH132"]
high <- s132>57
studentnames <- rownames(smat)

s132
high
studentnames
studentnames[high]

high
!high  # ! is NOT : flips TRUE and FALSE

TRUE
!TRUE

!high
studentnames[!high]

ss <- round(runif(12, 30,80))
ss

ss[ss>60]
ss[!(ss>60)]
ss[ss==60]
ss[ss==31]

ss[ss!=31]
ss[ss==31]
ss[ss>31]
ss[ss>=31]
ss[!(ss>=41)]
# same as 
ss[ss<41]

ss>=40
ss<=60
ss>=40 & ss<=60  # & is AND
# can't write 40<=ss<=60
ss<40 | ss>60  # | is OR

ss
ss[ss<40]
ss[ss>60]
ss[ss<40 | ss>60]
ss[!(ss>=40 & ss<=60)]

?xor

# categorical variables
a <- rownames(smat)
a

a=="Mike"
a=="mike"

a
a %in% c("Mike","Julie","Terence","Brad")

c("Mike","Julie","Terence","Brad") %in% a

b <- c("Mike","Julie","Terence","Brad")
a[a %in% c("Mike","Julie","Terence","Brad")]

a[a%in%b]


b[b%in%a]
b[!(b%in%a)]

a
sort(a)
sort(a, decreasing=TRUE)

# Special Values in R

# NA= not available
x <- c(1, 3, 2, NA, 5)
x

pops <- c(Auckland=1495000, Wellington=405000, Christchurch=389700,
          Dunedin=NA)
pops

sum(pops)
# can't do this:
pops==NA
# instead must do this:
is.na(pops)

pops[is.na(pops)]
pops[!is.na(pops)]
sum(pops[!is.na(pops)])

sum(pops)
sum(pops, na.rm=TRUE)

min(pops, na.rm=TRUE)

# NULL

# NaN = Not a Number
sqrt(-1)

b <- sqrt(c(1,2,3,-1))
b

b[!is.nan(b)]


sort(pops)

NA
NULL
NaN
Inf

x <- c(1,4,6,0,1)
1/x

30*(Inf+9)
Inf-90

Inf-Inf

-Inf

-1/0

1/Inf

#########################
# Reading data
# usually store in a data.frame

cities <- data.frame(city=c("Wellington","Auckland","Dunedin"),
                     pop=c(405000, 1495000, 118500),
                     area=c(444, 1086, 255))
cities

is.data.frame(cities)

dim(cities)
names(cities) # gives the column names

colnames(cities)
rownames(cities)
dimnames(cities) # a list object with two entries

dimnames(cities)
dimnames(cities)[[1]]
dimnames(cities)[[2]]

# accessing a column
cities[,3]
cities[,"area"]
cities$area


cities[2,]
cities["2",]


# Existing data sets inside R (or in an R package)
data()

trees # is an existing dataset
trees <- trees # makes a local copy of the existing data set trees

library(MASS)
data(Animals)
dim(Animals)

Animals
rownames(Animals)
Animals["Goat",]
Animals[4,]

Animals["Goat",2]
Animals["Goat","brain"]
Animals$brain
Animals$brain[4]
Animals$brain["Goat"] # this one fails

# Read an external file
list.files()
getwd()

# a csv file:
titanic <- read.csv("titanic.csv", stringsAsFactors=FALSE)
titanic

nrow(titanic)
head(titanic)


titanic <- read.table("titanic.csv", stringsAsFactors=FALSE,
                      sep=",", header=TRUE)

barplot(table(titanic$pclass))

# Webscraping


#### Writing out data sets
atab

write.csv(atab, file="airport.csv")

atab
barplot(table(atab$Airline), col='green', cex.names=0.5)


## Summarising data
surf <- read.csv("surf.csv", stringsAsFactors=FALSE)
dim(surf)

# top rows of a data frame
head(surf)
# structure of a data frame
str(surf)

surf[1,]
surf$Gender
surf[,"Gender"]
surf[,2]

range(surf$Age)
min(surf$Age)

surf[which.min(surf$Age),]
surf[which.max(surf$Income),]

table(surf$Age)
surf$Age==min(surf$Age)
mean(surf$Income[surf$Age==min(surf$Age)])

surf$Age[surf$Age==min(surf$Age)]
mean(surf$Age[surf$Age==min(surf$Age)])

surf[surf$Age==min(surf$Age),]

# deleting a column
surf$X
surf$X <- NULL
surf$X
surf[1:2,]

surf$Hello

# a brief summary of the columns in a data frame
summary(surf)

summary(surf$Age)
summary(surf$Gender)

mean(surf$Age)
mean(surf$Gender)

summary(c(1,2,6,NA))
mean(c(1,2,6,NA))
mean(c(1,2,6,NA), na.rm=TRUE)

quantile(surf$Age)
?quantile

quantile(surf$Income, probs=c(0.05, 0.05, 0.20, 0.50 ,0.95))

boxplot(c(surf$Age,65,70), col="red", ylab="Age")








###############################
###############################
###############################
#####                     #####
##### Week 3:             #####
#####   27 Jul 2020       #####
#####                     #####
###############################
###############################
###############################

surf <- read.csv("surf.csv", stringsAsFactors=FALSE)

surf <- read.table("surf.csv", sep=",", header=TRUE, 
                   stringsAsFactors=FALSE, fill=TRUE)
surf$X <- NULL
surf

boxplot(surf$Age, ylab="Age (yr)", main="Age Distribution",
        col="purple", ylim=c(0,60))

boxplot(surf$Age, ylab="Age (yr)", main="Age Distribution",
        col="purple", ylim=c(0,max(surf$Age)))

boxplot(c(surf$Age,100), ylab="Age (yr)", main="Age Distribution",
        col="purple")
?boxplot

hist(surf$Age, xlab="Age (yr)", main="Age Distribution", col="grey")

hist(surf$Age, xlab="Age (yr)", main="Age\nDistribution", col="grey")

mean(surf$Age)
summary(surf$Age)

# Categorical variables
surf$Gender

surf$Qualification
is.character(surf$Qualification)

mean(surf$Qualification) # doesn't make sense

tt <- table(surf$Qualification) # counts within categories
tt


tt["degree"]
tt[c("degree","school")]

sum(tt[c("degree","school")])

tt
barplot(table(surf$Qualification), 
        cex.names=0.8, cex.axis=0.8,
        ylab="Frequency", main="Qualification\nDistribution")
# cex = character expansion factor

surf <- read.csv("surf.csv", stringsAsFactors=FALSE)
is.character(surf$Qualification) # TRUE
is.numeric(surf$Qualification)   # FALSE
as.numeric(surf$Qualification)   # lots of NA's


surf$Qualification <- factor(surf$Qualification)
is.character(surf$Qualification) # FALSE
is.numeric(surf$Qualification)   # FALSE
as.numeric(surf$Qualification)   # lots of NA's
is.factor(surf$Qualification)    # TRUE
levels(surf$Qualification)
# default order of factors is alphabetical

# remove the current factor levels
surf$Qualification <- as.character(surf$Qualification)
is.factor(surf$Qualification)

surf$Qualification <- factor(surf$Qualification,
                             levels=c("none","school",
                                      "vocational","degree"))

barplot(table(surf$Qualification), 
        cex.names=0.8, cex.axis=0.8,
        ylab="Frequency", main="Qualification\nDistribution")

surf$Qualification

# beware that we have to be exact in specifying the levels
surf$Qualification <- as.character(surf$Qualification)
surf$Qualification <- factor(surf$Qualification,
                             levels=c("None","School",
                                      "Voc","Uni"))
# degree was misspelled here - caused a problem
barplot(table(surf$Qualification), 
        cex.names=0.8, cex.axis=0.8,
        ylab="Frequency", main="Qualification\nDistribution")
surf$Qualification

surf <- read.csv("surf.csv", stringsAsFactors=FALSE)

levs <- unique(surf$Qualification)
levs
levs <- levs[c(3,2,1,4)]
levs

surf$Qualification <- factor(surf$Qualification, levels=levs)
barplot(table(surf$Qualification), 
        cex.names=0.8, cex.axis=0.8,
        ylab="Frequency", main="Qualification\nDistribution")

levels(surf$Qualification) <- c("None","School","Voc.","Uni")
barplot(table(surf$Qualification), 
        cex.names=0.8, cex.axis=0.8,
        ylab="Frequency", main="Qualification\nDistribution")

table(surf$Qualification)
any(is.na(surf$Qualification))

surf$Qualification[1] <- NA
any(is.na(surf$Qualification))

table(surf$Qualification)
table(surf$Qualification, exclude=NULL)

# Two variables at once
# Two numerical

plot(surf$Age, surf$Income, pch=16, col="red", 
     cex=surf$Hours/70, xlab="Age", ylab="Weekly Income")

split(surf$Income, surf$Qualification)
# creates a list object - one element of the list for each 
# level of Qualification

boxplot(split(surf$Income, surf$Gender), col="blue")

table(surf$Qualification, surf$Gender, surf$Marital)

barplot(table(surf$Qualification, surf$Gender))
barplot(table(surf$Qualification, surf$Gender), beside=TRUE)
barplot(table(surf$Qualification, surf$Gender), beside=TRUE,
        legend=TRUE, col=c("red","yellow","brown","green"))

levels(surf$Qualification)

as.numeric(surf$Qualification)

# these two statements do exactly the same thing...
as.character(surf$Qualification)
levels(surf$Qualification)[as.numeric(surf$Qualification)]

###############################################
# Functions

sqrt(5)

square <- function(x) {
  return(x^2)
}

h <- square(5)
h


square <- function(x) {
  # Square the value x
  return(x^2)
}

square(3)
square

h

?square
?sqrt

square(3)
square(x=3)

args(square)


v <- 1:10
square(v)

a <- array(1:12, dim=c(3,4))
a

square(a)

# arguments
args(square)

square(5)
square(c(1,4,-2))
square() # error

square <- function(x=0) {
  # Square the value x
  return(x^2)
}

square()

args(mean)
?mean

x <- rnorm(10, mean=0, sd=3)
x <- c(x, NA)
x

mean(x, na.rm=TRUE)
# x is identified by its position
# na.rm is identifed by its name

mean(x=x, na.rm=TRUE) # fully name matching
mean(na.rm=TRUE, x=x) # fully name matching

y <- c(3,4,2)
mean(x=y, na.rm=TRUE)
mean(y, na.rm=TRUE)
x

mean(na.rm=TRUE, y)
mean(y)

mean(y, na.rm=TRUE)
mean(y, na.r=TRUE)
mean(y, na.=TRUE)
mean(y, na=TRUE)
mean(y, n=TRUE)

# Reserved words
?reserved

sqrt <- 4
sqrt

x <- c(3,4,1,2,3)
x

a <- 1
b <- 20
c <- 12
d <- 1

# these break
for <- 3
TRUE <- 3

Inf <- 4 # breaks
inf <- 4 # OK


# Scope
sqrt
sqrt(9)

library(maps)
ls(pos=10)

# refer explicitl to sqrt() in the base package
base::sqrt

library(dplyr)
# masks the filter() function in stats
stats::filter
dplyr::filter
filter

filter <- function(x) return(x)

filter("a")
filter(27)
filter(filter)
rm(filter)
filter() # error

trees

# make a new object in the global environ with the value of 
# the first thing you can find called "trees"
trees <- trees


trees <- "h"
trees <- trees
rm(trees)
trees <- trees

trees <- "h"
trees
datasets::trees

# scoping inside functions
square(5)




square <- function(x) {
  # Square the value x
  print(x)
  print(jjk)
  return(x^2)
}

square(34)

jj <- square(34)

square(c(6,6,4))



badsquare <- function() {
  # Square the value x
  return(x^2)
}

x <- 5
badsquare()

# More programming

# add up the elements of a vector

x <- c(2,4,3,1)
n <- length(x)
n
# sum(x)
1:n

sumx <- function(x) {
  # Add up the elements of a vector x
  n <- length(x)
  result <- 0
  
  for(i in 1:n) {
    print(i)
    result <- result + x[i]
  }
  
  return(result)
}

s <- sumx(1:10)

sum(c(1,4,2,1))


for(i in 1:4) {
  print(i)
}

digits <- 0:9
letters
LETTERS
digits

for(c in LETTERS) {
  print(c)
}

surf <- read.csv("surf.csv", stringsAsFactors = FALSE)
names(surf)
surf[1:3,]

for(nn in names(surf)) {
  print(nn)
  print(length(nn))
  print(nchar(nn))
}

names(surf)
nn <- names(surf)[1]
nn

nchar(nn)
length(nn)
length(names(surf))
nchar(names(surf))



sumx <- function(x) {
  # Add up the elements of a vector x
  if(!is.numeric(x)) {
    stop("x needs to be numeric!!!! (idiot)")
  }
  n <- length(x)
  result <- 0
  if(n>0) {
    for(i in 1:n) {
      result <- result + x[i] 
    }
  }
  return(result)
}

sumx(1:100)
sumx(names(surf))

x <- numeric()
length(x)

x <- 1:3
x
x <- x[-2]
x
x <- x[-c(1,2)]
x

sumx(x)

sum(c(1,NA,3), na.rm=TRUE)
sumx(c(1,NA,3))



sumx <- function(x, na.rm=FALSE, verbose=FALSE) {
  # Add up the elements of a vector x
  if(!is.numeric(x)) {
    stop("x needs to be numeric!!!! (idiot)")
  }
  n <- length(x)
  result <- 0
  if(n>0) {
    for(i in 1:n) {
      if(verbose) print(i) 
      if(is.na(x[i]) & na.rm==TRUE) {
        next # go to the next iteration of the loop
      } else {
        result <- result + x[i] 
      }
      
    }
  }
  return(result)
}


sumx(1:10)
sumx(c(1:10,NA))
sumx(c(1:10,NA), na.rm=TRUE)
sumx(c(1:10,NA), na.rm=TRUE, verbose=TRUE)
sumx(c(1:10,NA), na.rm=TRUE, verbose=FALSE)


x <- 7
y <- sum(runif(10))

x <- 7; y <- sum(runif(10))

a <- 3; b <- 4
c <- sqrt(a^2+b^2)
c
a <- 3
c <- 1
if(c==5) {
  print("Hello")
} else if(c>5) {
  print("Big")
} else if(c==1) {
  print("1")
} else if(a==3) {
  print("oops")
} else if(c>5) {
  print("Big!!!!")
} else {
  print("Small")
}

if(c==5) 
  print("Hello")

ss <- function(x) return(x+2)
ss(8)

# while()

n <- 5
while(n>0) {
  print(n)
}

while(TRUE) {
  print("forever.")
}

while(FALSE) {
  print("forever.")
}


#finished <- FALSE
#while(!finished) {
#  
#  if() {
#    finished <- TRUE
#  }
#  
#}

# Applying functions

amat <- array(rnorm(100),dim=c(25,4))
amat
hist(amat)

mean(amat)

# mean in every row
apply(amat, 1, mean)  # dimension 1 = row
# mean in every column
apply(amat, 2, mean)  # dimension 2 = column
# max in each column
apply(amat, 2, max)  # dimension 2 = column

apply(amat, 2, 
      function(x) {
        return(sort(x)[1:8])
      })


# family of apply-type functions
surf[1:3,]

mean(surf$Income)

tapply(surf$Income, surf$Marital, mean)

tapply(surf$Income, surf$Marital, length)
tapply(surf$Hours, surf$Marital, length)
table(surf$Marital)

summary(surf$Income)
tapply(surf$Income, surf$Marital, summary)









###############################
###############################
###############################
#####                     #####
##### Week 4:             #####
#####   3 Aug 2020        #####
#####                     #####
###############################
###############################
###############################

trees <- trees


write.csv(trees, file="trees.csv")

write.table(trees, file="trees.csv",
            row.names=FALSE,
            sep=",", 
            na="")

tt <- table(surf$Marital, surf$Ethnicity)
tt

write.table(tt, file="ME.csv",
            row.names=TRUE,
            col.names=NA,
            sep=",", 
            na="")

# tab separated 
write.table(tt, file="ME.txt",
            row.names=TRUE,
            col.names=NA,
            sep="\t", 
            na="")

# Writing to text files

sink("notes.txt")
print("The surf dataset has")
print(nrow(surf))
print("rows")
sink()

print(nrow(surf))

sink("notes.txt") # deletes the existing file
print("hello")
sink()

sink("notes.txt", append=TRUE) # adds to the end of an existing file
print("there")
sink()

sink()

# graphics files
#   jpeg
#   png
#   svg
#   tiff
#   postscript

stdwidth <- 1000; stdheight <- 600
jpeg("hoursincome.jpg", width=stdwidth, height=stdheight)
plot(surf$Hours, surf$Income, xlab="Hours", ylab="Income",
     main="Hours and Income", pch=16, col="blue",
     las=2)
dev.off()

# par() - customising graphics
par()


par(mar=c(7,7,7,7))

par(mar=c(5.1,4.1,4.1,2.1))
plot(surf$Hours, surf$Income, xlab="Hours", ylab="Income",
     main="Hours and Income", pch=16, col="blue",
     las=2)


par(mfrow=c(2,1))
plot(surf$Hours, surf$Income, xlab="Hours", ylab="Income",
     main="Hours and Income", pch=16, col="blue",
     las=2)
plot(surf$Hours, surf$Age, xlab="Hours", ylab="Age",
     main="Hours and Age", pch=16, col="blue",
     las=2)

stdwidth <- 1000; stdheight <- 1000
jpeg("fourgraphs.jpg", width=stdwidth, height=stdheight)
par(mfrow=c(2,2))
plot(surf$Hours, surf$Income, xlab="Hours", ylab="Income",
     main="Hours and Income", pch=16, col="blue",
     las=2)
plot(surf$Hours, surf$Age, xlab="Hours", ylab="Age",
     main="Hours and Age", pch=16, col="blue",
     las=2)
plot(surf$Hours, surf$Income, xlab="Hours", ylab="Income",
     main="Hours and Income", pch=16, col="blue",
     las=2)
plot(surf$Hours, surf$Age, xlab="Hours", ylab="Age",
     main="Hours and Age", pch=16, col="blue",
     las=2)
dev.off()


surf <- read.csv("surf.csv", stringsAsFactors = FALSE)

par(mfcol=c(1,1))
plot(surf$Hours, surf$Age, xlab="Hours", ylab="Age",
     main="Hours and Age", pch=16, col="blue",
     las=2)

surf <- read.csv("H:/_links/data202/2020T2-prvt/Richard/notes/week3/surf.csv", stringsAsFactors = FALSE)

list.files()


write.csv(Orange, "orangetrees.csv")
file.exists("orangetrees.csv")
file.copy("orangetrees.csv", "a.csv")
file.remove("a.csv")
list.files()

fname <- "orangetrees.csv"
file.copy(fname, "a.csv")

system("ls")
system("date")
system('whoami')

aa <- Sys.getenv("PATH")

.Machine
.Machine$double.xmax

.Platform

R.Version()

s <- read.csv("surf.csv", stringsAsFactors=FALSE)
list.files("data")


### Strings ####

print("The number of rows of surf is")
print(nrow(surf))


airlines <- read.csv("airlines.csv", stringsAsFactors=FALSE,
                     encoding="latin1")
nrow(airlines) # 6107
airlines[1:2,]

# Exact matching of a string
airlines[airlines$Country=="New Zealand",]

# substrings
nchar("New Zealand") # 11
nz <- "New Zealand"
nz

substring(nz,1,4)
substring(nz,1)
substring(nz,6)
substring(nz,6,10)

countries <- sort(unique(airlines$Country))
countries

nchar(countries)
substring(countries,1,4)

thai <- airlines[airlines$Country=="Thailand",]
nrow(thai) # 54
THAI <- airlines[airlines$Country=="THAILAND",]
nrow(THAI) # 1

nz
tolower(nz)
toupper(nz)

allthai <- airlines[tolower(airlines$Country)=="thailand",]
nrow(allthai)

airlines[tolower(substring(airlines$Airline,1,3))=="roy",]

?grep
grep("royal", tolower(airlines$Airline))  # location numbers
grepl("royal", tolower(airlines$Airline)) # vector of TRUE and FALSE values

airlines[grep("royal", tolower(airlines$Airline)),]

# regular expressions
grep("Ro?a", airlines$Airline)
grep("Ro??a", airlines$Airline)
grep("Ro*a", airlines$Airline)
grep("^Ro*a", airlines$Airline)

# paste, paste0

a <- c("Fred","Wilma","Barney","Betty")
b <- c("Flintstone","Flintstone","Rubble","Rubble")

paste(a,b,sep=" ")
paste(a,b,sep="")
paste(a,b,sep="::::")
paste(a,b,b,a,sep="::::")

paste("Hello my name is",a,sep=" ")

paste("The number of rows in surf is",nrow(surf))

# paste0 has no separator
paste0(a," ",b)

print(paste(a,b))

cat(paste(a,b,"\n"))

i <- 1
i
print(i)

for(i in 1:10) {
  i
}

for(i in 1:10) {
  print(i)
}

for(i in 1:10) {
  cat(i)
}

for(i in 1:10) {
  cat("[1] ")
  cat(i)
  cat("\n")
}

r <- 4
area <- pi*r^2

cat("I have a circle with radius "); cat(r); cat(", and its area is "); cat(area); cat("\n")

paste("I have a circle with radius ", r, ", and its area is ",
      area, "\n", sep="")

cat(paste("I have a circle with radius ", r, ", and its area is ",
          area, "\n", sep=""))

sprintf("I have a circle with radius %.3f, and its area is %.3f\n",
        r, area)

# %s for strings
# %f for real numbers (decimals)
# %d for integers
# %% for a single percent sign

a 
b
ages <- c(44, 26, 32, 34)
scores <- c(23, 56, 21, 5)

cat(sprintf("%s %s is %d years old and scored %d%% on that test\n", 
            a, b, ages, scores))

# left justified strings
cat(sprintf("%-8s %-10s is %.1f years old and scored %d%% on that test\n", 
            a, b, ages, scores))

# right justified strings
cat(sprintf("%8s %10s is %.1f years old and scored %2d%% on that test\n", 
            a, b, ages, scores))

cat(sprintf("%8s %10s is %.1f years old and scored %02d%% on that test\n", 
            a, b, ages, scores))

times <- c("01:22:23.4", "01:20:16.3", "00:58:22.1", "1:12:33.0")
times

# to get the hours:
substring(times,1,2) # doesn't work

a <- t(array(as.numeric(unlist(strsplit(times, split=":"))),dim=c(3,4)))
colnames(a) <- c("hours","mins","secs")
a

airlines[1:3,]
unique(unlist(strsplit(airlines$Airline[1:10],split=" ")))

a <- c("Fred","Wilma","Barney","Betty")
b <- c("Flintstone","Flintstone","Rubble","Rubble")
c <- paste(a,b)
c

strsplit(c, split=" ")

d <- as.data.frame(strsplit(c, split=" "))
names(d) <- 1:4
d[1,]
d[2,]
install.packages('wordcloud')
library(wordcloud)
words <- unlist(strsplit(airlines$Airline,split=" "))
d <- as.data.frame(table(words))
d[1:10,]

odx <- order(d$Freq,decreasing=TRUE)
d[odx,][1:10,]

wordcloud(d$words, freq=d$Freq, min.freq=1, max.words=100)

wordcloud(d$words, freq=d$Freq, min.freq=1, max.words=100, 
          colors=brewer.pal(8,"Dark2"))
title(sprintf("the time is %s",date()))

#######  Reshaping datasets ######## 
Orange
dim(Orange)
Orange[1:3,]

par()$mar
par(mar=c(5.1,4.1,4.1,2.1))
plot(Orange$age, Orange$circumference, pch=as.character(Orange$Tree), 
     xlab="Age", ylab="Circumference")


Orange

plot(Orange$age, Orange$circumference, pch=16,
     xlab="Age", ylab="Circumference")

by(Orange, 
   Orange$Tree,
   function(hello) {
     print(hello)
     lines(hello$age, hello$circumference)
   })


wideOrange <- reshape(Orange,
                      idvar="Tree",
                      timevar="age",
                      v.name="circumference",
                      direction="wide")
wideOrange

names(wideOrange)[2:8]

longOrange <- reshape(wideOrange,
                      idvar="Tree",
                      varying=names(wideOrange)[2:8],
                      timevar="age",
                      direction="long")
longOrange



###############################
###############################
###############################
#####                     #####
##### Week 5:             #####
#####   10 Aug 2020       #####
#####                     #####
###############################
###############################
###############################
Sys.Date()
Sys.time()
system("date")

birthdays <- c("23/11/1963","1/1/1970","3/12/1977","7/7/1990","1/1/2000","5/7/2019")
birthdays
# this are character strings

bd <- as.Date(birthdays, format="%d/%m/%Y")
typeof(bd)
class(bd)

as.numeric(bd) # gives number of days since 1/1/1970

format(bd, format="%d-%m-%Y")
format(bd, format="%m-%d-%Y")
format(bd, format="%m-%d-%y")
format(bd, format="%a, %m-%b-%y")
format(bd, format="%A, %m %B %Y")

# Date times
birthdays <- c("23/11/1963 12:30:00",
               "1/1/1970 00:00:00","3/12/1977 22:02:00",
               "7/7/1990 00:00:00","1/1/2000 12:12:12","5/7/2019 25:10:10")
birthdays

bdt <- as.POSIXct(birthdays, format="%d/%m/%Y %H:%M:%S")
bdt

bdtl <- as.POSIXlt(birthdays, format="%d/%m/%Y %H:%M:%S")
bdtl

as.numeric(bdt)
as.numeric(bdtl)
bdtl


sat <- read.csv("active-satellites.csv", stringsAsFactors=FALSE)

sat$Date.of.Launch[1:3]

sat$launchdate <- as.Date(sat$Date.of.Launch, format="%m/%d/%Y")
typeof(sat$launchdate)
class(sat$launchdate)

hist(sat$launchdate, breaks="years")

# can treat dates like numbers
max(sat$launchdate, na.rm=TRUE)
min(sat$launchdate, na.rm=TRUE)

max(sat$launchdate, na.rm=TRUE)-min(sat$launchdate, na.rm=TRUE)

any(is.na(sat$launchdate))
sat[is.na(sat$launchdate),]

file.remove("test_db.sqlite")

library(RSQLite)

# open a connection to a SQLite database
test_conn <- dbConnect(RSQLite::SQLite(), "test_db.sqlite")

# close the connection
dbDisconnect(test_conn)


# reconnect
test_conn <- dbConnect(RSQLite::SQLite(), "test_db.sqlite")


dbListTables(test_conn)
system("ls -l")

trees <- trees

# write this into the database
dbWriteTable(test_conn, "treedata", trees)
dbListTables(test_conn)
system("ls -l")

dbExistsTable(test_conn, "treedata")
dbExistsTable(test_conn, "red")

dbListFields(test_conn, "treedata")
names(trees)

hh <- dbReadTable(test_conn, "treedata")
hh

dbWriteTable(test_conn, "treedata", hh, overwrite=TRUE)
dbWriteTable(test_conn, "treedata", hh)

dbRemoveTable(test_conn, "red")
dbRemoveTable(test_conn, "treedata")

surf <- read.csv("surf.csv", stringsAsFactors=FALSE)
dbWriteTable(test_conn, "surfshort", surf[1:10,1:8], overwrite=TRUE)

dbListTables(test_conn)
dbReadTable(test_conn, "surfshort")

# * = everything (all rows, all columns)
# SELECT * FROM ...
dbGetQuery(test_conn, "SELECT * FROM surfshort")

ss <- dbGetQuery(test_conn, "SELECT * FROM surfshort")
ss









###############################
###############################
###############################
#####                     #####
##### Week 9:             #####
#####   21 September 2020 #####
#####                     #####
###############################
###############################
###############################


############################
## Demonstrate the use of ##
## the sample() function. ##
############################

# Set seed for replicability of results.
set.seed(0)

# Randomly sample from the outcomes of 1, 3, and 7.
sample(x = c(1, 3, 7), size = 8, replace = TRUE, prob = c(0.1, 0.3, 0.6))

# Randomly sample from the outcomes of "a", "Cat"", and "e9f".
sample(x = c("a", "Cat", "e9f"), size = 4, replace = TRUE, prob = c(0.1, 0.3, 0.6))

####################################
## Simulate flipping a fair coin. ##
## Use the sample() function.     ##
####################################

N <- 10000

set.seed(0)

outcomes <- sample(x = c("0", "1","2","3"), size = N, replace = TRUE, prob = c(0.503, 0.325,0.166,0.006))

coin.flips <- data.frame(table(outcomes))

coin.flips$Rel.Freq <- coin.flips$Freq / N


####################################
## Simulate flipping a fair coin. ##
## Use the runif() function.      ##
####################################

# Specify the number of coin flips.
N <- 1000000

# Set the seed for replicability of results.
set.seed(0)

# Generate N random numbers uniformly between 0 and 1.
r <- runif(N)

# Construct an empty vector in which to store the outcomes of the flips.
outcomes <- rep(NA, N)

# Use a "for" loop to assign outcomes corresponding to each random number between 0 and 1.
for(i in 1 : N)
{
  if(r[i] <= 0.5)
  {
    outcomes[i] <- "Heads"
  } else
  {
    outcomes[i] <- "Tails"
  }
}

# Alternative way to assign outcomes corresponding to each random number between 0 and 1 using the ifelse() function.
outcomes <- ifelse(test = (r <= 0.5), yes = "Heads", no = "Tails")

# Construct a data frame of the frequency of heads and tails.
coin.flips <- data.frame(table(outcomes))
# Calculate the relative frequency of heads and tails and add as a variable to the data frame.
coin.flips$Rel.Freq <- coin.flips$Freq / N
coin.flips

# Construct a barplot of the relative frequency of heads and tails.
barplot(coin.flips$Rel.Freq)

####################################
# Simulate flipping a biased coin. #
# Use the sample() function.       #
####################################


outcomes<- sample(x = c(0, 1, 2,3), size = 10000, replace = TRUE, prob = c(0.503, 0.325, 0.166,0.006))

var(outcomes)

coin.flips <- data.frame(table(outcomes))
coin.flips$Rel.Freq <- coin.flips$Freq / 10000
coin.flips

EX <- sum(coin.flips$Rel.Freq * coin.flips$Freq)
EX
new.x <- (coin.flips$Freq - EX) ^ 2
rbind(new.x, coin.flips$Rel.Freq)

VX <- sum(coin.flips$Rel.Freq * new.x)
VX
# Construct a barplot of the relative frequency of heads and tails.
barplot(coin.flips$Rel.Freq)

####################################
# Simulate flipping a biased coin. #
# Use the runif() function.        #
####################################

# Specify the number of coin flips.
N <- 1000000

# Set the seed for replicability of results.
set.seed(0)

# Generate N random numbers uniformly between 0 and 1.
r <- runif(N)

# Construct an empty vector in which to store the outcomes of the flips.
outcomes <- rep(NA, N)

# Use a "for" loop to assign outcomes corresponding to each random number between 0 and 1.
for(i in 1 : N)
{
  if(r[i] <= 0.3)
  {
    outcomes[i] <- "Heads"
  } else
  {
    outcomes[i] <- "Tails"
  }
}

# Alternative way to assign outcomes corresponding to each random number between 0 and 1 using the ifelse() function.
outcomes <- ifelse(test = (r <= 0.3), yes = "Heads", no = "Tails")

# Construct a data frame of the frequency of heads and tails.
coin.flips <- data.frame(table(outcomes))
# Calculate the relative frequency of heads and tails and add as a variable to the data frame.
coin.flips$Rel.Freq <- coin.flips$Freq / N
coin.flips

# Construct a barplot of the relative frequency of heads and tails.
barplot(coin.flips$Rel.Freq)

############################################
# Simulate four outcomes with different    #
# probabilities Use the sample() function. #
############################################

# Consider outcomes of 1, 2, 3, and 4.
# Probabilities for these outcomes are 0.2, 0.3, 0.1, and 0.4.

# Specify number of times to sample.
N <- 1000000

# Set the seed for replicability of results.
set.seed(1)

outcomes <- sample(x = 1 : 4, size = N, replace = TRUE, prob = c(0.2, 0.3, 0.1, 0.4))

# Construct a data frame of the frequency of the outcomes of 1 to 4.
random.draws <- data.frame(table(outcomes))
# Calculate the relative frequency of the outcomes of 1 to 4 and add as a variable to the data frame.
random.draws$Rel.Freq <- random.draws$Freq / N
random.draws

# Construct a barplot of the relative frequency of the outcomes of 1 to 4.
barplot(random.draws$Rel.Freq)

############################################
# Simulate four outcomes with different.   #
# probabilities Use the runif() function.  #
############################################

# Consider outcomes of 1, 2, 3, and 4.
# Probabilities for these outcomes are 0.2, 0.3, 0.1, and 0.4.

# Specify the number of times to sample.
N <- 1000000

# Set the seed for replicability of results.
set.seed(0)

# Generate N random numbers uniformly between 0 and 1.
r <- runif(N)

# Construct an empty vector in which to store the outcomes of the random draws.
outcomes <- rep(NA, N)

# Use a "for" loop to assign outcomes corresponding to each random number between 0 and 1.
for(i in 1 : N)
{
  if(r[i] <= 0.2)
  {
    outcomes[i] <- 1
  } else if(r[i] <= 0.5)
  {
    outcomes[i] <- 2
  } else if(r[i] <= 0.6)
  {
    outcomes[i] <- 3
  } else
  {
    outcomes[i] <- 4		
  }
}

# Construct a data frame of the frequency of the outcomes of 1 to 4.
random.draws <- data.frame(table(outcomes))
# Calculate the relative frequency of the outcomes of 1 to 4 and add as a variable to the data frame.
random.draws$Rel.Freq <- random.draws$Freq / N
random.draws

# Construct a barplot of the relative frequency of the outcomes of 1 to 4.
barplot(random.draws$Rel.Freq)





#********************#********************#********************#********************#********************#********************
#*#********************#********************#********************#********************#********************#********************
#*#********************#********************#********************#********************#********************#********************
  ###############################
###############################
###############################
#####                     #####
##### Week 10:            #####
#####   28 September 2020 #####
#####                     #####
###############################
###############################
###############################


# Construct the probability mass function.
x <- 0 : 3
px <- c(0.25 ,0.35, 0.1, 0.3)

# Combine x and px to visualise the probability mass function.
rbind(x, px)

####################
## Expected value ##
####################

# Calculate the expected value.
EX <- sum(px * x)
EX

# Approximate the expected value using simulation.
N <- 1000000

# Set the seed.
set.seed(0)

outcomes <- sample(x = 0 : 3, size = N, replace = TRUE, prob = px)
mean(outcomes)

##############
## Variance ##
##############

# Calculate the variance.
new.x <- (x - EX) ^ 2
rbind(new.x, px)

VX <- sum(px * new.x)
VX

# Calculate the standard deviation.
SD <- sqrt(VX)
SD

# Approximate the variance and standard deviation using simulation.
N <- 1000000

# Set the seed.
set.seed(0)

outcomes <- sample(x = 0 : 3, size = N, replace = TRUE, prob = px)
var(outcomes) # Variance
sd(outcomes) # Standard deviation

##############################
## R's built-in probability ##
## distribution functions.  ##
##############################

# Generate N = 1,000,000 random outcomes from a Normal(0, 1) distribution.
outcomes <- rnorm(N)

# Visualise the distribution of these N = 1,000,000 outcomes.
hist(outcomes) # Histogram
plot(density(outcomes)) # Density plot

# Generate N = 1,000,000 random Uniform(0, 1) outcomes, representing cumulative probabilities.
r <- runif(N)
# Use the qnorm() function to inverse transform these cumulative probabilities to corresponding quantile outcomes for a Normal(0, 1) distribution.
outcomes2 <- qnorm(r)

# Visualise the distribution of these N = 1,000,000 outcomes in red on the previous density plot.
plot(density(outcomes))
points(density(outcomes2), type = "l", col = 2)

##########################################
## Discrete uniform [a, b] distribution ##
##########################################

# Approximate the expected value using simulation.
N <- 1000000

# Set the seed.
set.seed(0)

outcomes <- sample(x = 1 : 6, size = N, replace = TRUE)
mean(outcomes)

# Function to calculate probabilities for a discrete uniform distribution
dunifdis <- function(x, a, b)
{
  if(x >= a & x <= b & x == round(x))
  {
    prob <- 1 / ((b - a) + 1)
  } else
  {
    prob <- 0
  }
  
  return(prob)
}

dunifdis(x = 1, a = 1, b = 6)
dunifdis(x = 1.12, a = 1, b = 6)

############################
## Bernoulli distribution ##
############################

# A Bernoulli distribution is an instance of a binomial distribution (discussed shortly) 
# where the sample size is 1 (specifically, X∼Ber(p) is the same as X∼Bin(1,p)), 
# so functions for probabilities, quantiles, and random generation of outcomes 
# for a Bernoulli distribution use those of the binomial distribution.

# Calculate probabilities for individual outcomes.
dbinom(x = 0, size = 1, prob = 0.7) # P(X = 0)
dbinom(x = 1, size = 1, prob = 0.7) # P(X = 1)
dbinom(x = 0.5, size = 1, prob = 0.7) # P(X = 0.5)

pbinom(q = 0.5, size = 1, prob = 0.7) # P(X <= 0.5)

# Convert cumulative probabilities to corresponding outcomes.
qbinom(0.2, size = 1, prob = 0.7)
qbinom(0.5, size = 1, prob = 0.7)

# Generate 10 random outcomes for a Ber(0.7) distribution.
rbinom(n = 10, size = 1, prob = 0.7)

# Approximate the expected value and variance.
outcomes <- rbinom(n = 1000000, size = 1, prob = 0.7)
mean(outcomes)
var(outcomes)

###########################
## Binomial distribution ##
###########################

# Store possible outcomes and probabilities for a Bin(10, 0.7) distribution.
x <- 0 : 10
px <- dbinom(0 : 10, size = 10, prob = 0.7)

# Visualise the probability mass function.
rbind(x, px)

# Calculate the expected value.
EX <- sum(px * x)
EX

# Approximate the expected value and variance.
# Generate 1,000,000 outcomes for a Bin(10, 0.7) distribution.
outcomes <- rbinom(1000000, size = 10, prob = 0.7)
mean(outcomes)
var(outcomes)

# Generate 1 outcome for a Bin(10, 0.7) distribution.
rbinom(1, size = 10, prob = 0.7)
# Generate 10 Ber(0.7) outcomes and sum to produce 1 outcome for a Bin(10, 0.7) distribution.
sum(rbinom(10, size = 1, prob = 0.7))
dbinom(10,size = 1, prob = 0.7)
##########################
## Poisson distribution ##
##########################

# Calculate Poisson probabilities
dpois(0 : 5, lambda = 4) # P(X = 0), P(X = 1), ..., P(X = 5)
ppois(5, lambda = 4) # P(X <= 5)
p <- sum(dpois(0 : 5, lambda = 4)) # P(X <= 5)
p

qpois(p, lambda = 4) # Inverse transformation.

##############################
## Continuous distributions ##
##############################

# Demonstrate how R generates normal random variables using the rnorm() function.
r <- runif(1000000) # Generate numbers uniformly between 0 and 1.
outcomes <- qnorm(r) # Use qnorm() for inverse transformation.

# Visualise the resulting outcomes using a histogram and density plot.
hist(outcomes)
plot(density(outcomes))


###############################
###############################
###############################
#####                     #####
##### Week 10:            #####
#####   28 September 2020 #####
#####                     #####
###############################
###############################
###############################


########################
## Coin-Flipping Game ##
########################

########################
## Monty Hall Problem ##
########################

# Specify the number of times to play the game.
N <- 1000000

# Create vectors of length N in which to store game results.
win.stick <- rep(NA, N)
win.change <- rep(NA, N)

# Set the seed for replicability of results.
set.seed(0)

# Use a "for" loop to play the game multiple times.
for(i in 1 : N)
{
  # Prize gets placed behind a door.
  prize.door <- sample(1 : 3, size = 1)
  
  # Contestant selects a door.
  contestant.door <- sample(1 : 3, size = 1)
  
  # Monty Hall reveals a door without the prize.
  monty.options <- 1 : 3
  
  # Approach 1: Have Monty systematically choose the first (or only) door available to him.
  #monty.door <- (monty.options[monty.options != prize.door & monty.options != contestant.door])[1]
  
  # Approach 2: Have Monty randomly sample one of the doors available to him.
  monty.options <- monty.options[monty.options != prize.door & monty.options != contestant.door]
  
  if(length(monty.options) > 1)
  {
    monty.door <- sample(monty.options, size = 1)
    
  } else
  {
    monty.door <- monty.options
  }
  
  # Option 1: Stick with initial choice.
  win.stick[i] <- contestant.door == prize.door
  
  # Option 2: Change doors.
  contestant.options <- 1 : 3
  contestant.door <- contestant.options[contestant.options != monty.door & contestant.options != contestant.door]
  
  win.change[i] <- contestant.door == prize.door
}

# Observe relative frequencies of winning for the two strategies. 
table(win.stick) / N
table(win.change) / N

#############################
## Newsboy Freddie Problem ##
#############################

# Specify values for how much Freddy pays for newspapers, how much he sells them for, and how much he is reimbursed for unsold newspapers.
freddy.cost <- 1.50
selling.price <- 2.50
reimbursement <- 0.50

# Set the seed for replicability of results.
set.seed(0)

# Specify number of simulations
N <- 10000

# Construct a vector in which to store profits.
profit <- rep(NA, N)

# Cycle through all possible numbers of newspapers that Freddy can buy to see the average profitability associated with each number of newspapers.
for(buy in 40 : 70)
{
  # Carry out a simulation of N random days of sales.
  for(i in 1 : N)
  {
    # Determine demand for newspapers for a given day.
    demand <- sample(40 : 70, size = 1)
    
    # Freddy sells his newspapers.
    sell <- min(buy, demand)
    
    # Determine the number of unsold newspapers.
    unsold <- buy - sell
    
    # Determine the profit.
    profit[i] <- selling.price * sell - freddy.cost * buy  + reimbursement * unsold
  }
  
  # Print out the number of newspapers Freddy bought and what his average profitability was.
  print(c(buy, mean(profit)))
}



###############################
###############################
###############################
#####                     #####
##### Week 11:            #####
#####   6 October 2020    #####
#####                     #####
###############################
###############################
###############################


#####################
## Social networks ##
#####################

# Adjacency matrix of reported friendships between individuals.
friendship.aj <- matrix(c(NA, 0, 1, 1, 1, 
                          0, NA, 0, 0, 0, 
                          1, 0, NA, 0, 1, 
                          1, 1, 0, NA, 1, 
                          1, 0, 0, 0, NA), 
                        nrow = 5, byrow = TRUE)
friendship.aj

# Edge list of reported friendships between individuals.
friendship.el <- matrix(c(1, 3, 
                          1, 4, 
                          1, 5, 
                          3, 1, 
                          3, 5, 
                          4, 1, 
                          4, 2, 
                          4, 5, 
                          5, 1), 
                        ncol = 2, byrow = TRUE)
friendship.el

# Load the "network" package.
library(network)

# Save the network to a "network" object.
friendship.net <- network(friendship.aj)
friendship.net
friendship.net <- network(friendship.el)
friendship.net

# Extract the adjacency matrix.
as.matrix(friendship.net)
as.sociomatrix(friendship.net)

# Extract the edge list.
as.edgelist(friendship.net)
as.matrix(friendship.net, matrix.type = "edgelist")

########################
## Bipartite Networks ##
########################

# Affiliation matrix of students and rooms in which they have classes.
affiliation.mat <- matrix(c(1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 
                            0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 
                            0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 
                            0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                            1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0), 
                          nrow = 5, byrow = TRUE, 
                          dimnames = list(c("Julio", "Kimiko", "Larry", "Monoa", "Nhung"), 
                                          c("AM101", "ELT206", "GBLT1", "GBLT2", "HMLT205", "HMLT206", "MCLT101", "MCLT102", "MLT220", "NKLT303", "RHLT2", "SUMT228", "TTRLT1")))
affiliation.mat

# Edge list of students and rooms in which they have classes.
affiliation.el <- matrix(c(1, 6, 
                           1, 10,
                           1, 11, 
                           1, 17, 
                           2, 7, 
                           2, 12, 
                           2, 13, 
                           2, 17, 
                           2, 18, 
                           3, 8, 
                           3, 9, 
                           3, 14, 
                           3, 15, 
                           3, 16, 
                           4, 8, 
                           4, 9, 
                           5, 6, 
                           5, 11, 
                           5, 17), 
                         ncol = 2, byrow = TRUE, 
                         dimnames = list(NULL, c("Actor", "Event")))
affiliation.el

# Construct a bipartite network object from an affiliation matrix.
bipartite.net <- network(affiliation.mat)
bipartite.net

# Example where the network() function confuses an affiliation matrix for an adjacency matrix, so the "bipartite" argument is required.
bipartite.net <- network(affiliation.mat[, 1 : 5], directed = FALSE, bipartite = 5)
bipartite.net
as.matrix(bipartite.net)
# Construct a bipartite network object from an edge list.
bipartite.net <- network(affiliation.el, directed = FALSE, bipartite = 5)
bipartite.net

# Extract the affiliation matrix and edge list.
as.matrix(bipartite.net)
as.edgelist(bipartite.net)

# Extract the expanded sociomatrix/adjacency matrix.
as.matrix(bipartite.net, expand.bipartite = TRUE)

##############################
## Node and edge attributes ##
##############################

# Node attributes.
friendship.net %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")
friendship.net %v% "vertex.names"
friendship.net %v% "age" <- c(19, 18, 20, 22, 21)
friendship.net %v% "age"
# Edge attribute.
friendship.net %e% "type" <- c("Best friend", "Dating", "Friend", 
                               "Best friend", "Friend", "Dating", 
                               "Friend", "Dating", "Friend")
# 在有1的位置放入vertex type
friendship.net %e% "type"

# Extract adjacency matrix and edgelist with edge attributes.
as.matrix(friendship.net, attrname = "type")
as.matrix(friendship.net, matrix.type = "edgelist", attrname = "type")


# Distance sociomatrix
dist.sociomatrix <- matrix(c(0, 764.71, 1430.13, 113.89, 211.43, 643.19, 
                             764.71, 0, 361.56, 676.47, 714.66, 440.23, 
                             1430.13, 361.56, 0, 982.42, 1023.96, 615.57, 
                             113.89, 676.47, 982.42, 0, 78.66, 391.49, 
                             211.43, 714.66, 1023.96, 78.66, 0, 417.5, 
                             643.19, 440.23, 615.57, 391.49, 417.5, 0), 
                           nrow = 6, byrow = TRUE, dimname = list(c("Auckland", "Christchurch", "Dunedin", "Hamilton", "Tauranga", "Wellington"), c("Auckland", "Christchurch", "Dunedin", "Hamilton", "Tauranga", "Wellington")))

# Construct a network object with edge weights saved to an edge attribute called "distance".
dist.net <- network(dist.sociomatrix, directed = FALSE, ignore.eval = FALSE, names.eval = "distance")
dist.net

# Extract the sociomatrix.
as.matrix(dist.net)
as.matrix(dist.net, attrname = "distance")

###########################
## Network Visualisation ##
###########################

# Edge list of reported friendships between individuals.
friendship.el <- matrix(c(1, 3, 
                          1, 4, 
                          1, 5, 
                          3, 1, 
                          3, 5, 
                          4, 1, 
                          4, 2, 
                          4, 5, 
                          5, 1), 
                        ncol = 2, byrow = TRUE) # byrow !!!
friendship.el

# Load the "network" package.
library(network)

# Save the network to a "network" object.
friendship.net <- network(friendship.el)

# Node attributes.
friendship.net %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")
friendship.net %v% "age" <- c(19, 18, 20, 22, 21)
# Edge attribute.
friendship.net %e% "type" <- c("Best friend", "Dating", "Friend", "Best friend", "Friend", "Dating", "Friend", "Dating", "Friend")

# Visualise the network.
plot(friendship.net, displaylabels = TRUE, vertex.cex = 2, edge.lwd = 2, vertex.col = friendship.net %v% "age", edge.col = 2)

########################
## Add Health dataset ##
########################

# Load the "rda" package to make use of datasets.
library(rda)
# Read in the "addhealth9" dataset from the "rda" package.
data(addhealth9)


# Construct a "network" object using the edge list.
addhealth.net <- network(addhealth9$E, ignore.eval = FALSE, names.eval = "intensity")
# Add node attributes for sex, race, and grade.
addhealth.net %v% "sex" <- addhealth9$X[, 1]
addhealth.net %v% "race" <- addhealth9$X[, 2]
addhealth.net %v% "grade" <- addhealth9$X[, 3]

# Visualise the network.
plot(addhealth.net)

# Decorate the graph to highlight friendships by sex and race.
plot(addhealth.net, 
     vertex.sides = ifelse(addhealth.net %v% "sex" == 1, 50, 3),
     vertex.col = addhealth.net %v% "race",
     edge.col = gray(0.5))

# Decorate the graph to highlight friendships by sex and grade.
plot(addhealth.net, 
     vertex.sides = ifelse(addhealth.net %v% "sex" == 1, 50, 3),
     vertex.col = addhealth.net %v% "grade",
     edge.col = gray(0.5))

########################
## Network statistics ##
########################

# Adjacency matrix of reported friendships between individuals.
friendship.aj <- matrix(c(NA, 0, 1, 1, 1, 
                          0, NA, 0, 0, 0, 
                          1, 0, NA, 0, 1, 
                          1, 1, 0, NA, 1, 
                          1, 0, 0, 0, NA), 
                        nrow = 5, byrow = TRUE)

# Construct both directed and undirected versions of the friendship network.
friendship.net.dir <- network(friendship.aj)
friendship.net.undir <- network(friendship.aj, directed = FALSE)

friendship.net.dir %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")
friendship.net.undir %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")

# Adjacency
is.adjacent(friendship.net.undir, 2, 4)
is.adjacent(friendship.net.undir, 4, 2)

is.adjacent(friendship.net.dir, 2, 4)
is.adjacent(friendship.net.dir, 4, 2)

# Incidence
get.inducedSubgraph(friendship.net.undir, eid = 1) %v% "vertex.names"


###############################
###############################
###############################
#####                     #####
##### Week 12:            #####
#####   13 October 2020    #####
#####                     #####
###############################
###############################
###############################


########################
## Network statistics ##
########################

# Load the "network" package.
library(network)

# Adjacency matrix of reported friendships between individuals.
friendship.aj <- matrix(c(NA, 0, 1, 1, 1, 
                          0, NA, 0, 0, 0, 
                          1, 0, NA, 0, 1, 
                          1, 1, 0, NA, 1, 
                          1, 0, 0, 0, NA), 
                        nrow = 5, byrow = TRUE)

# Construct both directed and undirected versions of the friendship network.
friendship.net.dir <- network(friendship.aj)
friendship.net.undir <- network(friendship.aj, directed = FALSE)

# Store names corresponding to nodes.
friendship.net.dir %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")
friendship.net.undir %v% "vertex.names" <- c("Alice", "Bob", "Carol", "Dennis", "Eunice")

# Visualise the directed and undirected networks.
plot(friendship.net.dir, displaylabel = TRUE, vertex.cex = 3, edge.lwd = 2, arrowhead.cex = 1.5)

plot(friendship.net.undir, displaylabel = TRUE, vertex.cex = 3, edge.lwd = 2)

#####################
## Social dynamics ##
#####################

# Load the "ergm" and "sna" packages.
library(ergm)
library(sna)

# "sna" functions: grecip, dyad.census, gtrans
grecip(friendship.net.dir)
dyad.census(friendship.net.dir)
gtrans(friendship.net.dir)
gtrans(friendship.net.dir, measure = "weakcensus")

grecip(friendship.net.undir)
dyad.census(friendship.net.undir)
gtrans(friendship.net.undir, measure = "weakcensus")

# "ergm" terms: asymmetric, mutual, transitive triangle
summary(friendship.net.dir ~ asymmetric + mutual + transitive)
summary(friendship.net.undir ~ triangles)

###########
## Flows ##
###########

# "sna" functions: degree, gden
degree(friendship.net.dir, cmode = "indegree")
degree(friendship.net.dir, cmode = "outdegree")
degree(friendship.net.dir) # Total = indegree + outdegree
degree(friendship.net.undir, gmode = "graph")
gden(friendship.net.undir)
gden(friendship.net.dir)

network.size(friendship.net.dir)

# "ergm" terms: degree/idegree/odegree, density
summary(friendship.net.undir ~ degree(0 : 4) + density)
summary(friendship.net.dir ~ idegree(0 : 4) + odegree(0 : 4) + density)

##################
## Connectivity ##
##################

# "sna" functions: reachability, geodist, is.connected, components, cutpoints
reachability(friendship.net.dir)
reachability(friendship.net.undir)
geodist(friendship.net.undir)
geodist(friendship.net.dir)
is.connected(friendship.net.dir)
is.connected(friendship.net.undir)
cutpoints(friendship.net.undir, return.indicator = TRUE) #Dennis is the cutpoint

# "ergm" function: ergm.geodistdist
ergm.geodistdist(friendship.net.dir)

# Confirm that the directed network without Bob is connected.
reduced.network <- get.inducedSubgraph(friendship.net.dir, v = c(1, 3, 4, 5)) 
# get rid of Bob then it is connected
is.connected(reduced.network)

################
## Centrality ##
################

# "sna" functions: degree, centralization
centralization(friendship.net.undir, FUN = "degree", mode = "graph")
