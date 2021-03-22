####################################################
## Coded by Darwin Pucha Cofrep
##  Plot rings  V2.3
##  Natural shape of rings
## plotRings Make a plot and/or animation of a cross section based on up to four ring-width series. 
## Gives basic summary statistics (e.g. annual basal area, mean ring width) of an approximated stem disc.
## Date: may.2018
## Updates: 24.may.2018, 
## Last update: 21.mar.2021
###################################################



rm(list = ls())     # clear objects 
graphics.off()      # close/clean graphics windows  



#############################################################################################

plotRings2 <- function(trwN= NA_real_, trwS = NA_real_,
                      trwE = NA_real_, trwW = NA_real_,
                      year,
                      ring.ticks = FALSE,
                      col.ring.ticks = "grey40",
                      tick.freq = 10,
                      axis.plot = F,
                      length.unit = "mm",
                      animation = FALSE, 
                      sys.sleep = 0.1, 
                      year.labels = TRUE, 
                      d2pith = NA,
                      col.inrings = "grey", 
                      lwd = 1,
                     #col.outring = "black",      # This argument is needed to adapt to this function
                     # x.rings = "none", col.x.rings = "red",      # This argument is needed to adapt to this function
                      xy.lim = "auto.lim",
                      saveGIF=FALSE, fname="GIF_plotRings.gif",   # This argument is needed to adapt to this function
                      title = "Annual tree growth",
                      text2 = NA,
                      species.name = NA
                      ) {

  op <- par(no.readonly=TRUE)  # save default graphing parameters to reset
  
  ## Creating a data.frame
  'TRW <- data.frame(row.names = year, trwN = trwN, 
                    trwS = if (exists("trwS") == TRUE) 
                      trwS
                    else NA, trwE = if (exists("trwE") == TRUE) 
                      trwE
                    else NA, trwW = if (exists("trwW") == TRUE) 
                      trwW
                    else NA)'
  
  TRW <- data.frame(row.names = year, 
                    trwN = if (exists("trwN") == TRUE)
                      trwN, 
                    trwS = if (exists("trwS") == TRUE) 
                      trwS
                    else NA, trwE = if (exists("trwE") == TRUE) 
                      trwE
                    else NA, trwW = if (exists("trwW") == TRUE) 
                      trwW
                    else NA)
 
  #Copy TRW columns
   TRW[, c("N", "S", "E", "W")] <- TRW[,1:4]
  
  # Fill the empty columns with the mean
    # detect the empty columns
      emptyCol <- apply(TRW[, 5:8], 2, function(myCol) {sum(!is.na(myCol)) == 0 } )
  # fill with the mean from the columns with values, T = empty cols, F = cols with data
      TRW[, names(emptyCol)[emptyCol == T]] <- rowMeans(TRW[, names(emptyCol)[emptyCol == F] ])
  
  ## Setting the length unit of ring measurement
  if(length.unit == "mm") 
    TRW[, 5:8] <- TRW[, 5:8]    
  else if(length.unit == "1/10 mm") 
    TRW[, 5:8] <- TRW[, 5:8]/10 
  else if(length.unit == "1/100 mm") 
    TRW[, 5:8] <- TRW[, 5:8]/100
  else if(length.unit == "1/1000 mm") 
    TRW[, 5:8] <- TRW[, 5:8]/1000 
  
  ## Setting the unit of d2pith
  if(length.unit == "mm") 
    d2pith <- d2pith
  else if(length.unit == "1/10 mm") 
    d2pith <- d2pith/10 
  else if(length.unit == "1/100 mm") 
    d2pith <- d2pith/100
  else if(length.unit == "1/1000 mm") 
    d2pith <- d2pith/1000
  
  
TRW <- TRW[as.logical((rowSums(is.na(TRW))-length(TRW))),] # It is to remove rows with NAs across all rows
  

# Distance to pith (d2pith)
# Add a new row with the d2pith value at the first position if this argument is assigned
if (!is.na(d2pith)) {  
  TRW.d2pith <- TRW[1, 5:8]
  TRW.d2pith[1,] <- NA
  rownames(TRW.d2pith)[1] <- as.numeric(rownames(TRW))[1]-1 
  TRW.d2pith[1,] <- d2pith/2  # Total width length of the pith
  TRW <- rbind(TRW.d2pith, TRW)
}


# Data NA filling
library(sinkr)

TRW.NAfill <- eof(TRW[,5:8], recursive=TRUE) #  Convert to RSEOF - "Recursively Subtracted Empirical Orthogonal Functions" 
TRW.NAfill_rcn <- eofRecon(TRW.NAfill) # NA data Reconstruction
TRW.NAfill_rcn <- as.data.frame(TRW.NAfill_rcn) # as.data.frame

# matplot(TRW, type ="l")  # Plot visualization of original data
# matplot(TRW.NAfill_rcn, type ="l") # Plot visualization wiht NA reconstructed data values

TRW[,5:8] <- TRW.NAfill_rcn[, 1:4]  # Replacing data to the original data.frame



#data arrangement for the four cardinal axes

data <- data.frame(N = cumsum(TRW[,8])*-1,
                   S = cumsum(TRW[,7]), 
                   E = cumsum(TRW[,6])*-1,
                   W = cumsum(TRW[,5]) )


################
# Code based on the Answer of cuttlefish44 in
# http://stackoverflow.com/questions/40196151/how-to-draw-a-shape-ellipse-or-oval-following-some-points-and-calculate-its-ar/40201111#40201111


# change all data into xy coordinates and make ring-factor
library(reshape2); library(dplyr)

data <- t(data)
colnames(data) <-  rownames(TRW[, 5:8])          # ring-factor
df <- melt(data, value.name = "x")        # change into long-form
df$y <- df$x                              # make xy coordinates
df[df$Var1=="S"|df$Var1=="N", "y"] <- 0
df[df$Var1=="E"|df$Var1=="W", "x"] <- 0


#calculation of center coordinates, ox & oy
center <- df %>% group_by(Var2) %>% summarize(sum(x)/2, sum(y)/2) %>% as.data.frame()

#calculation of parameters of ellipse; semi-major and -minor axis, ra & rb
opt.f <- function(par, subset, center) {     # target function
  ox <- center[[1]]                          # par[1] and par[2] are ra and rb
  oy <- center[[2]]
  x <- subset$x
  y <- subset$y
  sum(abs((x - ox)^2/par[1]^2 + (y - oy)^2/par[2]^2 - 1))   # from ellipse equation
}

lev <- levels(as.factor(df$Var2))

## search parameters
res <- sapply(1:length(lev), function(a) 
  optim(c(1,1), opt.f, subset = subset(df, Var2 == lev[a]), 
        center = center[a, 2:3], control = list(reltol = 1.0e-12)))

res  # result. you can get detail by res[,1etc]. values are not 0 but much nearly 0

#function to plot (Probably some packages have similar one)
radian <- function(degree) degree/180*pi
plot.ellipse <- function(ox, oy, ra, rb, phi=0, start=0, end=360, length=100, func=lines, ...) {
  theta <- c(seq(radian(start), radian(end), length=length), radian(end))
  if (phi == 0) {
    func(ra*cos(theta)+ox, rb*sin(theta)+oy, ...)
  } else {
    x <- ra*cos(theta)
    y <- rb*sin(theta)
    phi <- radian(phi)
    cosine <- cos(phi)
    sine <- sin(phi)
    func(cosine*x-sine*y+ox, sine*x+cosine*y+oy, ...)
  }
}

# Auto limit
if(xy.lim == "auto.lim") {
max.xy <-  max(abs(df[3:4]))*1.1
tick.lim <- ceiling(max.xy/100)*100  # rounding number for ticks
ticks <- seq(-tick.lim, tick.lim, tick.freq) 
ticks[ticks == 0] <- NA   # remove the value 0
}
else
{ max.xy <- xy.lim
  ticks <- seq(-xy.lim, xy.lim, tick.freq) 
ticks[ticks == 0] <- NA   # remove the value 0
}

### Plot rings

if(animation == TRUE) {
  for(a in 1:length(lev)) {
  par(mar=c(1,1,1,1)+0.1, mai=c(1,1,1,1), cex=1, xaxs="i",yaxs="i")
  
    plot(0, type="n", xlim=c(-max.xy, max.xy), ylim =c(-max.xy, max.xy), asp=1, xlab="", 
     ylab="", axes = FALSE,  bg = 'transparent',
     main = mtext(bquote(~bold(.(title))),line=3.5,adj=0.5, side=3, cex=1.5), 
     sub=if(!is.na(species.name)) mtext(bquote(~plain(.("(")) ~italic(.(species.name)) ~plain(.(")"))),
     line=2.5,adj=0.5, side=3, cex=1) )

  mtext("South", side = 1, line = 0.5); mtext("West", side = 2, line = 0.1)
  mtext("North", side = 3, line = 0.5); mtext("East", side = 4, line = 0.1)
  mtext(text2, side = 1, line = 0, at = 500, cex= 0.85, col = "grey")

  # Add year labels
  leg.year <- as.numeric(row.names(TRW))
  
  if(year.labels == TRUE) {legend('topright', legend=leg.year[a], inset = 0.01, cex=1.2, bg = "white", bty="o", box.col = "white")}
 
   # Axes
  if(axis.plot == TRUE) {axis(1, pos=0, labels = T, at = ticks);  axis(2, pos=0, las=2, labels = T, at = ticks) } else
    {axis(1, pos=0, labels = F, at = c(-max.xy, max.xy));  axis(2, pos=0, las=2, labels = F, at = c(-max.xy, max.xy)) }
 # if(ring.ticks == TRUE) {points(df$x, df$y, pch=as.character(rep(0:20, each = 4)), col=col.ring.ticks, cex=0.8)}  #add ticks marks
  if(ring.ticks == TRUE) {points(df$x, df$y, pch=3, col=col.ring.ticks, cex=0.8)}  #add ticks marks
 
  plot.ellipse(ox = center[a, 2], oy = center[a, 3], col = c(col.inrings), 
               ra = res[,a]$par[1], rb = res[,a]$par[2], length=300, lwd = 4)

  Sys.sleep(sys.sleep) 
  par(new=TRUE)
   }  
  par(op) # restore original graphic settings
  } 

  else
  {for(a in 1:length(lev)) {
      par(mar=c(1,1,1,1)+0.1, mai=c(1,1,1,1), cex=1, xaxs="i",yaxs="i")
      
      plot(0, type="n", xlim=c(-max.xy, max.xy), ylim =c(-max.xy, max.xy), asp=1, xlab="", 
           ylab="", axes = FALSE, 
           main = mtext(bquote(~bold(.(title))),line=3.5,adj=0.5, side=3, cex=1.5), 
           sub=if(!is.na(species.name)) mtext(bquote(~plain(.("(")) ~italic(.(species.name)) ~plain(.(")"))),
                                              line=2.5,adj=0.5, side=3, cex=1), lwd = lwd )
      
      mtext("South", side = 1, line = 0.5); mtext("West", side = 2, line = 0.1)
      mtext("North", side = 3, line = 0.5); mtext("East", side = 4, line = 0.1)
      mtext(text2, side = 1, line = 0, at = 500, cex= 0.85, col = "grey")
      
      #testing for a secondary x axis
      #if(axis.plot == TRUE) {axis(1, pos=0, labels = T, at = ticks);  axis(2, pos=0, las=2, labels = T, at = ticks) } 
      if(axis.plot == TRUE) {
        axis.labels.interval <- 100
        axis(side=1, line = 2, at = seq(-max.xy, max.xy, axis.labels.interval), 
             labels = seq(-max.xy, max.xy, axis.labels.interval), col = "black", 
             col.ticks="black", col.axis="black", cex.axis = 0.7)
        mtext("Widht [mm]",1,line=4,at=0.2,col="black")
        abline(h = 0, v = 0)
        }      
        else
        {axis(1, pos=0, labels = F, at = c(-max.xy, max.xy));  axis(2, pos=0, las=2, labels = F, at = c(-max.xy, max.xy)) }
      # if(ring.ticks == TRUE) {points(df$x, df$y, pch=as.character(rep(0:20, each = 4)), col=col.ring.ticks, cex=0.8)}  #add ticks marks
      if(ring.ticks == TRUE) {points(df$x, df$y, pch=3, col=col.ring.ticks, cex=0.8)}  #add ticks marks
      
      
      
      plot.ellipse(ox = center[a, 2], oy = center[a, 3], col = c(col.inrings), 
                   ra = res[,a]$par[1], rb = res[,a]$par[2], length=300, lwd = lwd)
      # Add year labels
      leg.year <- as.numeric(row.names(TRW))
      if(year.labels == TRUE) {legend('topright', legend=leg.year[a], inset = 0.01, cex=1.2, bg = "white", bty="o", box.col = "white")}
      
      par(new=TRUE)
  }
par(op) # restore original graphic settings    
  }    
  

###################################################
###################################################
#  CALCULATIONS

# trw means
TRW$trw.means <- rowMeans(TRW[,5:8], na.rm = TRUE)

# Accumulative trw.means
TRW$trw.mean.acc <- cumsum(TRW$trw.means)

# Accumulative mean Diameter (trw.means x2)
TRW$diam.acc.mean <- TRW$trw.mean.acc*2 

# Accumulative Diameter N-S

# Accumulative Diameter E-W

# Accumulative trw.means - MAI
#TRW$year.seq  <- seq(1, length(year), 1)
#TRW$trw.acc.MAI <- TRW$trw.acc/TRW$year.seq

TRW$year.seq  <- seq(1, nrow(TRW), 1)
TRW$trw.acc.MAI <- TRW$trw.mean.acc/TRW$year.seq


# Accumulative trw.means - CAI
TRW$trw.CAI <- ffcsaps(TRW$trw.means, x = as.numeric(row.names(TRW)), nyrs = 12)



# Basal Area calculation
# BA cummulative
TRW$Cum.B.Area <- sapply(res[1,], function(a) pi * a[1] * a[2])
# BA individual (inside out)
TRW$B.Area.ind  <-c(TRW$Cum.B.Area[1], TRW$Cum.B.Area[2:nrow(TRW)] - TRW$Cum.B.Area[1:nrow(TRW)-1])  



## Print Report:  
print("REPORT")
TRW


}




