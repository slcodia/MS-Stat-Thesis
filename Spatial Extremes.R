# SPATIAL EXTREMES
library(SpatialExtremes)
library(tidyverse)
# WIND: Spatial GEV - the rough way ===================
require(maps)
rm(list=ls())
data(wind)

par(mar = rep(0, 4))
maps::map(xlim = c(0, 9), ylim = c(47.5, 57.5), 
          fill = TRUE, col = "grey")
points(coord[,1:2], pch = 15)

#Exploratory plot

symbolplot(wind, coord, which = "mean")
symbolplot(wind, coord)
border <- function(add = FALSE){
    maps::map(xlim = c(0,9), ylim = c(47.5, 57.5), add = add)}
symbolplot(wind,coord, plot.border = border)


# Defining the trend surface function using R base
loc.form <- y~lon*lat
scale.form <- y~1
shape.form <- y~1

# convert coordinates from degrees to km (rough estimate)
coord[,1:2] <- 111*coord[,1:2]

# fit the model
M0 <- fitspatgev(wind, scale(coord, scale = FALSE), 
                 loc.form, scale.form, shape.form)
M0

View(coord)

# prediction
x <- seq(min(coord[,1]), max(coord[,1]), length = 100)
y <- seq(min(coord[,2]), max(coord[,2]), length = 100)

grid <- expand.grid(x,y); colnames(grid) <- c("lon","lat")

# since we scaled the data, we need to use the same transformation
grid[,1] <- grid[,1] - mean(coord[,1])
grid[,2] <- grid[,2] - mean(coord[,2])

ans <- predict(M0, newdata = grid, ret.per = 20)$Q20

maps::map(xlim =range(x)/111, ylim = range(y)/111)

image(x/111, y/111, matrix(ans,100), add = TRUE, col = cm.colors(64))
contour(x/111,y/111, matrix(ans,100), add = TRUE); maps::map(add = TRUE)


# RAIN: Bayesian Hierarchical Modeling ============
rm(list=ls())
View(rain)
data(rainfall) # rainfall
data(swissalt) # elevation in Switzerland

par(mar=rep(0,4))
image(lon.vec, lat.vec, alt.mat, col = terrain.colors(64), asp = 1,
      bty = "n", xlab = "n", axes = FALSE)
swiss(add = TRUE, city = TRUE)
points(coord, pch = 15)

# setup - running the Gibbs Sampler
hyper <- list()

#Multivariate normal for the Beta
hyper$betaMeans <-  list(loc   = rep(0,3), 
                         scale = rep(0,3), 
                         shape = 0)

hyper$betaIcov  <-  list(loc   = diag(rep(1/1000,3)), 
                         scale = diag(rep(1/1000,3)), 
                         shape = 1/10)
# Inverse gamma distributions sills
hyper$sills     <-  list(loc   = c(1,12), 
                         scale = c(1,1), 
                         shape = c(1,0.04))

# Gamma distributions for ranges and smoothing parameters
hyper$ranges    <-  list(loc   = c(5,3),
                         scale = c(1,1),
                         shape = c(5,3))
hyper$smooths   <-  list(loc   = c(1,1),
                         scale = c(1,1),
                         shape = c(1,1))
#proposal distribution for the GEV, range, and smooth parameters
prop <- list(gev = c(3, 0.1, 0.3),
             ranges = c(1,0.8, 1.2),
             smooths = rep(0,3))

# initial state of for the chain
start <- list(sills = c(10, 10, 0.5), 
              ranges = c(20, 10, 10),
              smooths = c(1, 1, 1),
              beta = list(loc = c(25, 0, 0),
                          scale = c(33, 0, 0),
                          shape = 0.001))

loc.form   <- y~lon+lat
scale.form <- y~lon+lat
shape.form <- y~1
View(rain)
chain <- latent(rain, coord[,1:2], "powexp",
                loc.form, scale.form, shape.form, 
                hyper = hyper, prop = prop,
                start = start, n = 100, burn.in = 500, thin = 5)
View(coord[,1:2])
chain
View(rain)
View(coord)


# prediction of the pointwise 20-year return level
map.latent(chain, ret.per = 20, plot.contour = FALSE)
View(chain)
# once in a 20 year maximum rainfall


# using PH data ----

# reading data
rain.ph <- read.csv("maxdataspread.csv",fileEncoding = "UTF-8-BOM")
drop <- c("BASCO.RADAR","BALER.RADAR","ITBAYAT",
          "SANGLEY.POINT","MEDILLIN.CEBU","TANAY", "COTOBATO")
rain.ph <- rain.ph[-c(1:21), !(names(rain.ph) %in% drop)]

coord.ph <- read.csv("rainfalldata/PAGASA Weather Stations.csv",fileEncoding = "UTF-8-BOM")
coord.ph$STATION.NAME <- gsub(" ",".",coord.ph$STATION.NAME)

# cleaning data
rain.ph.clean <- rain.ph[-(1:3),] %>% 
    '['(,colSums(is.na(.)) == 0) %>%
    '['(,-c(1,2)) %>%
    as.matrix()

keep.stations <- colnames(rain.ph.clean)
coord.ph.clean <- subset(coord.ph, STATION.NAME %in% keep.stations)%>%
    `[`(,-c(1)) %>%
    as.matrix()


# running the MCMC algorithm
chain.ph <- latent(rain.ph.clean, coord.ph.clean[,1:2], "powexp",
                   loc.form, scale.form, shape.form, 
                   hyper = hyper, prop = prop,
                   start = start, n = 100, burn.in = 500, thin = 5)

map.latent.ph <- map.latent(chain.ph, x=seq(114,128, length.out = 100), y=seq(2,21,length.out=100),ret.per = 20, plot.contour = F)
View(map.latent.ph$ci.up)

par(mar = rep(0,4))
ph.map <- maps::map(xlim=c(114,128), ylim = c(2,21))
undebug(map.latent)

# COPULA MODEL
n.site <- 30
n.obs <- 50

coord <- matrix(runif(2*n.site, -10, 10), ncol = 2)
colnames(coord) <- c("lon", "lat")




# TEMPERATURE: max-stable process
rm(list = ls())

data(USHCNTemp) # get the USHCN temperature data

par(mar = rep(0, 4))
maps::map("usa")
points(metadata[,c("lon","lat")], pch = 20)



# Nonparametric

