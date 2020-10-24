library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggthemes)
library(zoo)
library(network)
library(ergm)
library(sna)
library(lubridate)


dim(HairEyeColor)
HEC <- HairEyeColor

HairEyeColor

count(HairEyeColor)

c <- as.data.frame(HairEyeColor)

h <- length(unique(c$Hair[c$Sex == "Female"]))
e <- length(unique(c$Eye[c$Sex == "Female"]))

f <- c(h,e)
f
z <- matrix(f)
z
m <- matrix(x, nrow=3)

c$Hair[c$Sex == "Female"]

barplot(table(c$Hair,c$Eye),beside = TRUE,ylab="Frequency", main="Color of eye and hair")


  pp <- function(x,y) {
    if (length(x)!=length(y)){
      stop("length of two elements should be same")
    }
    nx <- length(x)
    s <- 0
    for(i in 1:nx) {
      s <- s + x[i]/y[i]
      
    }
    return(s)
  }

 pp(c(2,2),c(3,1,3))

 weekdaynames <- c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")
 
 weekdaynames[grep("s",weekdaynames)]
 
 
 weekdaynames[grep("s",weekdaynames,fixed=TRUE)]
 
 
 weekdaynames[grep("s",tolower(weekdaynames))]
 
 
 weekdaynames[grep("s",weekdaynames,ignore.case=TRUE)]
 
 
 flights <- read.csv("statsNZ_tourism_flights_arrivals.csv", header=TRUE, stringsAsFactors = FALSE)
 flights$Date <- as.Date(as.yearqtr(flights$YearQuarter, format="%YQ%q"))
 

 flights
 flights.long <- pivot_longer(flights, cols=Adelaide:Vancouver, names_to="Origin", values_to="NumFlights")
 
 flights.median. <- median(flights.long$Date)
 flights.median.
 
 tapply(flights.long$NumFlights, flights.long$Date, median)
 

 
   new2 <- group_by(flights.long, YearQuarter,Origin)# group_by convert the data frame to a tibble 
   flights.max<- summarise(new2, Max=max(NumFlights),.groups = 'drop')
   
   flights.max<- arrange(flights.max, desc(Max)) 
 head(flights.max,5)
 
 flights.long4 <- filter(flights.long, flights.long$Origin %in% c("Sydney", "Melbourne", "Brisbane" , "Coolangatta"))
 flights.long4$Origin <- factor(flights.long4$Origin, levels= c("2", "1", "3", "4")) 
 levels(flights.long4$Origin)
 
 
 
 ggplot(earnings_long) + 
   geom_line(aes(x=Date, color=Industries, y=Average.Weekly.Earnings),stat="identity",size = 2)+ scale_color_colorblind() + labs(x = "year", y = "avg earning", title = "avg earning by year")
 
 ggplot(flights.long4) + 
   geom_line(aes(x=Date,linetype=Origin,y=NumFlights),stat="identity")
 
 
 flights.long4.2015 <- read.csv("flightslong2015.csv",header=TRUE)
 flights.long4.2015
 
 ggplot(flights.long4.2015) + 
   geom_bar(aes(x=Date,y=NumFlights,fill=Origin),stat="identity",position="dodge")
 +scale_y_continous(labels=percent)
 
 library(network)
                      m <- matrix(c(NA, 1, 0, 0, 0,  
                                     0, NA, 1, 0, 0,  
                                     1, 1, NA, 1, 0, 
                                     1, 0, 0, NA, 0,
                                     0, 0, 0, 1, NA),
                                   nrow = 5, byrow = TRUE,
                                   dimnames = list(c("1", "2", "3", "4", "5"),
                                                   c("1", "2", "3", "4", "5")))
 shipping_net <- network(m)
 plot(shipping_net,displaylabels=TRUE)
 dyad.census(shipping_net)
 degree(shipping_net, cmode = "indegree")
 is.adjacent(shipping_net, 4, 5)
 
 
 
 
 
 
