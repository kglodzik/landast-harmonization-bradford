# -------------------------------------------------------------------------------------
# Landsat Harmonization Analysis - Bradford Forest
# This script builds on processing steps in Google Earth Engine (GEE) and ArcGIS Pro
# It harmonizes surface reflectance values across Landsat 5, 7, and 8 to enable
# consistent LAI estimation across time.

# GEE steps:
#  - Extracted cloud-free Landsat imagery over Bradford Forest
#  - Identified time-matched image pairs within ~8 days of each other
#  - Exporting R and NIR bands from each sensor
# GEE script: https://code.earthengine.google.com/ccea0823362ae6d5b4dc6e8b59d60317

# ArcGIS steps (can be replicated in any spatial software):
#  - Generated random points at a minimize separation to avoid spatial autocorrelation
#    (we used 350 m based on screening using Incremental Spatial Autocorrelation tool)
#  - Extracted Red and NIR band values from each Landsat image export to those points

# Reproducing this workflow:
#  - Use included sample_data/ for Bradford Forest for a script demonstration
#  - Or create your own datasets of date-matched NIR and Red reflectance values

# NOTE:
# Variable names like 'regNIR' and 'regRed' are reused and overwritten throughout this script
# Similarly, intermediate datasets like 'L5andL7' and 'L8andL7' are filtered and refined
#  iteratively while maintaining their name
# Therefore be mindful when backtracking or re-running earlier parts non-sequentially

library(reshape2)
library(scales)

setwd('C:\\Users\\kglod\\UFL Dropbox\\Katie Glodzik\\Postdoc\\OtherProjects\\LandsatContinuity\\ArcGISoutputsForR\\sample_data')
L5_NIR     = read.csv('L5_NIR_BandValues_Bradford.csv'      ,header=TRUE); L5_NIR    =L5_NIR[,-2]
L5_Red     = read.csv('L5_Red_BandValues_Bradford.csv'      ,header=TRUE); L5_Red    =L5_Red[,-2]
L7for5_NIR = read.csv('L7withL5_NIR_BandValues_Bradford.csv',header=TRUE); L7for5_NIR=L7for5_NIR[,-2]
L7for5_Red = read.csv('L7withL5_Red_BandValues_Bradford.csv',header=TRUE); L7for5_Red=L7for5_Red[,-2]
L7for8_NIR = read.csv('L7withL8_NIR_BandValues_Bradford.csv',header=TRUE); L7for8_NIR=L7for8_NIR[,-2]
L7for8_Red = read.csv('L7withL8_Red_BandValues_Bradford.csv',header=TRUE); L7for8_Red=L7for8_Red[,-2]
L8_NIR     = read.csv('L8_NIR_BandValues_Bradford.csv'      ,header=TRUE); L8_NIR    =L8_NIR[,-2]
L8_Red     = read.csv('L8_Red_BandValues_Bradford.csv'      ,header=TRUE); L8_Red    =L8_Red[,-2]

names(L5_NIR)[1] = names(L7for5_NIR)[1] = names(L7for8_NIR)[1] = names(L8_NIR)[1] =
names(L5_Red)[1] = names(L7for5_Red)[1] = names(L7for8_Red)[1] = names(L8_Red)[1] = 'pixelID'

# Remove rows where all values were nan or 0 (most needing removal are nan's). Mostly for Landsat 7
# instances where points fell in the black NoData striping that occurs in Landsat 7 images.
L5_NIR    =L5_NIR    [(rowSums(is.na(L5_NIR))    <(ncol(L5_NIR)-1))     & (rowSums(L5_NIR[,-1],na.rm=TRUE)    !=0),] 
L5_Red    =L5_Red    [(rowSums(is.na(L5_Red))    <(ncol(L5_Red)-1))     & (rowSums(L5_Red[,-1],na.rm=TRUE)    !=0),]
L7for5_NIR=L7for5_NIR[(rowSums(is.na(L7for5_NIR))<(ncol(L7for5_NIR)-1)) & (rowSums(L7for5_NIR[,-1],na.rm=TRUE)!=0),]
L7for5_Red=L7for5_Red[(rowSums(is.na(L7for5_Red))<(ncol(L7for5_Red)-1)) & (rowSums(L7for5_Red[,-1],na.rm=TRUE)!=0),]
L7for8_NIR=L7for8_NIR[(rowSums(is.na(L7for8_NIR))<(ncol(L7for8_NIR)-1)) & (rowSums(L7for8_NIR[,-1],na.rm=TRUE)!=0),]
L7for8_Red=L7for8_Red[(rowSums(is.na(L7for8_Red))<(ncol(L7for8_Red)-1)) & (rowSums(L7for8_Red[,-1],na.rm=TRUE)!=0),]
L8_NIR    =L8_NIR    [(rowSums(is.na(L8_NIR))    <(ncol(L8_NIR)-1))     & (rowSums(L8_NIR[,-1],na.rm=TRUE)    !=0),]
L8_Red    =L8_Red    [(rowSums(is.na(L8_Red))    <(ncol(L8_Red)-1))     & (rowSums(L8_Red[,-1],na.rm=TRUE)    !=0),]

# It was good to have some redundancy to field names at first, as an accuracy check, but lets simplify a little now
names(L5_NIR)    =gsub('_NIR','',names(L5_NIR))
names(L5_Red)    =gsub('_R'  ,'',names(L5_Red))
names(L7for5_NIR)=gsub('_NIR','',names(L7for5_NIR))
names(L7for5_Red)=gsub('_R'  ,'',names(L7for5_Red))
names(L7for8_NIR)=gsub('_NIR','',names(L7for8_NIR))
names(L7for8_Red)=gsub('_R'  ,'',names(L7for8_Red))
names(L8_NIR)    =gsub('_NIR','',names(L8_NIR))
names(L8_Red)    =gsub('_R'  ,'',names(L8_Red))

#############################################################################
## Landsat 7 and 8
#############################################################################

L7L8DateMatches = read.csv('L7toL8MatchesTable_Bradford.csv', header=TRUE)

## Get each of the L5 and L7 datasets into long data frame form
L8_NIRlong    = melt(L8_NIR,     na.rm=FALSE, value.name='L8NIR',id="pixelID"); names(L8_NIRlong)[2]    ="L8date"
L7for8_NIRlong= melt(L7for8_NIR, na.rm=FALSE, value.name="L7NIR",id="pixelID"); names(L7for8_NIRlong)[2]="L7date"
L8_Redlong    = melt(L8_Red,     na.rm=FALSE, value.name='L8Red',id="pixelID"); names(L8_Redlong)[2]    ="L8date"
L7for8_Redlong= melt(L7for8_Red, na.rm=FALSE, value.name='L7Red',id="pixelID"); names(L7for8_Redlong)[2]="L7date"

# Add a column that gives the L7 date to be matched to each L8 date
L8andL7NIR = merge(L8_NIRlong, L7L8DateMatches,by='L8date',all.x=TRUE)

# Add a column that gives the L7 NIR values (found by linking to the L7 dates just added)
L8andL7NIR = merge(L8andL7NIR, L7for8_NIRlong ,by=c('pixelID','L7date'),all.x=TRUE)

# You now have a table with L8 & L7 NIR values. Now we add L8 & L7 Red values
L8andL7pre    = merge(L8andL7NIR, L8_Redlong    ,by=c('pixelID','L8date'),all.x=TRUE)
L8andL7pre    = merge(L8andL7pre, L7for8_Redlong,by=c('pixelID','L7date'),all.x=TRUE)

# Remove L7_20230324 records. The L7 striping must have been bad for that one because there are only 3 records.
L8andL7pre = L8andL7pre[L8andL7pre$L7date!="L7_20230324",]

# Copy the date frame and begin noise removal on the new data frame
L8andL7 = L8andL7pre[rowSums(!is.na(L8andL7pre))==ncol(L8andL7pre), ]

# For comparison, check regression now, BEFORE eliminating data after 2020 or other bad dates
regNIR = lm(L8NIR~L7NIR,data=L8andL7pre); summary(regNIR)$adj.r.squared
regRed = lm(L8Red~L7Red,data=L8andL7pre); summary(regRed)$adj.r.squared

## Check  whether there are any date pairs that should be excluded because of major differences
## in values, which would suggest a large-scale change happened at the study area in between
## the matched images

# First calculate point-level differences between matched pairs
L8andL7$NIRdiff = abs(L8andL7$L8NIR - L8andL7$L7NIR)
L8andL7$Reddiff = abs(L8andL7$L8Red - L8andL7$L7Red)

# Aggregate point differences by date
L7dates = unique(L8andL7$L7dateunique) ; L7dates

L8andL7NIRdiff_agg = aggregate(L8andL7$NIRdiff, list(L8andL7$L7dateunique),FUN=mean,na.rm=TRUE); 
L8andL7NIRdiff_agg$x = round(L8andL7NIRdiff_agg$x,4)
L8andL7NIRdiff_agg[order(-L8andL7NIRdiff_agg$x),] ## take a look -- any pairs with unusually high differences?

L8andL7Reddiff_agg = aggregate(L8andL7$Reddiff, list(L8andL7$L7dateunique),FUN=mean,na.rm=TRUE) 
L8andL7Reddiff_agg$x = round(L8andL7Reddiff_agg$x,4)
L8andL7Reddiff_agg[order(-L8andL7Reddiff_agg$x),] ## take a look -- any pairs with unusually high differences?

# Also calculate R-squared values within each date pair and compile them
L7L8NIR_Rsqs = L7L8Red_Rsqs = c()
for (x in L7dates){
  rsq = with(L8andL7[L8andL7$L7dateunique==x,], round(summary(lm(L8NIR~L7NIR))$adj.r.squared,3))
  L7L8NIR_Rsqs[which(L7dates==x)]=rsq
}
for (x in L7dates){
  rsq = with(L8andL7[L8andL7$L7dateunique==x,], round(summary(lm(L8Red~L7Red))$adj.r.squared,3))
  L7L8Red_Rsqs[which(L7dates==x)]=rsq
}

# Get the R-squared values and average differences for Red and NIR together in a data frame
L7L8Rsqs = data.frame(cbind(L7dates,L7L8NIR_Rsqs,L7L8Red_Rsqs))

L7L8Rsqs = merge(L7L8Rsqs, L8andL7NIRdiff_agg, by.x = 'L7dates', by.y='Group.1')
names(L7L8Rsqs)[length(L7L8Rsqs)]='L7L8NIR_diff'

L7L8Rsqs = merge(L7L8Rsqs, L8andL7Reddiff_agg, by.x = 'L7dates', by.y='Group.1')
names(L7L8Rsqs)[length(L7L8Rsqs)]='L7L8Red_diff'

# Graphing R-squared values and average differences for individual dates
xseq = seq(1,nrow(L7L8Rsqs))
par(mfrow=c(2,1),mar=c(5.2,2.9,0.7,3.1),mgp=c(1.8,0.5,0))
plot(L7L8Rsqs$L7L8NIR_Rsqs ~ xseq,xaxt='n',pch=16,type='b',tck=-0.02,ylab='NIR R-sq (black)',xlab='',ylim=c(0.33,0.96))
axis(1,at=xseq,labels=L7L8Rsqs$L7dates,cex.axis=0.8,las=2,tck=-0.015)
par(new=TRUE);
plot(L7L8Rsqs$L7L8NIR_diff ~ xseq,xaxt='n',yaxt='n',pch=16,type='b',tck=-0.02,col='blue',ylab='',xlab='',ylim=c(0.0032,0.0400))
axis(4) ; mtext('NIR mean |difference|', side=4,line=1.8)
#  L7_20181205, L7_20220404b looks bad
# L7_20221202 little bad but leave and see if later pixel-by-pixel filtering fixes it

plot(L7L8Rsqs$L7L8Red_Rsqs ~ xseq,xaxt='n',pch=16,type='b',tck=-0.02,ylab='Red R-sq (black)',xlab='',ylim=c(0.33,0.96))
axis(1,at=xseq,labels=L7L8Rsqs$L7dates,cex.axis=0.8,las=2,tck=-0.015)
par(new=TRUE);
plot(L7L8Rsqs$L7L8Red_diff ~ xseq,xaxt='n',yaxt='n',pch=16,type='b',tck=-0.02,col='blue',ylab='',xlab='',ylim=c(0.0032,0.0143))
axis(4) ; mtext('Red mean |difference|', side=4,line=1.8)
#  L7_20170116 looks bad

# Eliminate data after 2020. Spectral drift means Landsat 7 data is less reliable starting in 2021
# https://pubs.usgs.gov/publication/70224264
L8andL7$year = as.numeric(substr(L8andL7$L7date,4,7))
L8andL7=L8andL7[L8andL7$year<=2020,]

# Regression AFTER eliminating data after 2020
regNIR = lm(L8NIR~L7NIR,data=L8andL7); summary(regNIR)$adj.r.squared
regRed = lm(L8Red~L7Red,data=L8andL7); summary(regRed)$adj.r.squared

# Eliminate the other bad pairs
L8andL7=subset(L8andL7, subset = L7date != 'L7_20181205' & L7date != 'L7_20170116')

# Regression AFTER eliminating those other bad pairs
regNIR = lm(L8NIR~L7NIR,data=L8andL7); summary(regNIR)$adj.r.squared
regRed = lm(L8Red~L7Red,data=L8andL7); summary(regRed)$adj.r.squared

L8andL7=L8andL7[rowSums(!is.na(L8andL7))==ncol(L8andL7), ]

## More noise filtering at the pixel level, based on cases where L8 and L7 found very different values 
## for pixels, which again suggest real changes on the ground in between matched dates, though
## at a smaller scale now. First we use a screening method based on Roy 2016, which only removes 
## extremely misalligned data
L8andL7$NIRnoisecalc = abs(L8andL7$L8NIR-L8andL7$L7NIR) / (0.5*abs(L8andL7$L8NIR+L8andL7$L7NIR))
L8andL7$NIRnoisy     = ifelse(L8andL7$NIRnoisecalc > 1, 1, 0)

L8andL7$Rednoisecalc = abs(L8andL7$L8Red-L8andL7$L7Red) / (0.5*abs(L8andL7$L8Red+L8andL7$L7Red))
L8andL7$Rednoisy     = ifelse(L8andL7$Rednoisecalc > 1, 1, 0)

# set "noisy" field to 1 in cases where either band was flagged
L8andL7$noisy = ifelse((L8andL7$NIRnoisy | L8andL7$Rednoisy) == 1, 1, 0)

mean(L8andL7$noisy,na.rm=TRUE) # the portion removed via Roy 2016 adapted method. It will likely be tiny

# Environment clean up
L8andL7 = within(L8andL7, rm(NIRnoisy,Rednoisy,NIRnoisecalc,Rednoisecalc))
rm(L8_NIRlong,L7for8_NIRlong,L8_Redlong,L7for8_Redlong,L8andL7NIR)

## We need to do more here because there are still major disparities that suggest some type
#  of change on the ground between the paired Landsat dates.
#  Additional screening will be based simply on absolute differences.

## To set an unacceptable amount of disparity, calculate mean and st. dev. of the L8-to-L7 differences
#  and we'll say the mean + 3*sd is the cut-off. Anything higher is eliminated.
#  That is, we define an acceptable level of difference between date-matched pixels,
#  and if the difference exceeds that, we assume real change happened on the ground and eliminate it
L8andL7$NIRnoisy2 = ifelse(L8andL7$NIRdiff > 
                             mean(L8andL7$NIRdiff,na.rm=TRUE) + 3*sd(L8andL7$NIRdiff,na.rm=TRUE), 1, 0)
L8andL7$Rednoisy2 = ifelse(L8andL7$Reddiff > 
                             mean(L8andL7$Reddiff,na.rm=TRUE) + 3*sd(L8andL7$Reddiff,na.rm=TRUE), 1, 0)

# set "noisy2" field to 1 in cases where either band was flagged
L8andL7$noisy2 = ifelse((L8andL7$NIRnoisy2 | L8andL7$Rednoisy2) == 1, 1, 0) # Set to 1 if either is high
mean(L8andL7$noisy2,na.rm=TRUE) # the portion removed by new noise screening

# Environment Clean up
L8andL7 = within(L8andL7, rm(NIRdiff,Reddiff,NIRnoisy2,Rednoisy2))

# Regression BEFORE eliminating noisy data
regNIR = lm(L8NIR~L7NIR,data=L8andL7); summary(regNIR)$adj.r.squared
regRed = lm(L8Red~L7Red,data=L8andL7); summary(regRed)$adj.r.squared

# Remove the noisy data
L8andL7final=L8andL7[L8andL7$noisy==0 & L8andL7$noisy2==0,]
L8andL7final=L8andL7final[rowSums(!is.na(L8andL7final))==ncol(L8andL7final), ]
nrow(L8andL7final)

# Regression AFTER eliminating noisy data
regNIR_L7asL8 = lm(L8NIR~L7NIR,data=L8andL7final); summary(regNIR_L7asL8)$adj.r.squared
regRed_L7asL8 = lm(L8Red~L7Red,data=L8andL7final); summary(regRed_L7asL8)$adj.r.squared


# GRAPHING
# Graphing scatterplots at individual dates
L7dates = unique(L8andL7pre$L7dateunique) ; L7dates
par(mfrow=c(6,6),mar=c(2.2,2.1,0.7,0.6),mgp=c(1.0,0.1,0))
L7L8NIR_Rsqs = L7L8Red_Rsqs = c()

for (x in L7dates){
  with(L8andL7pre[L8andL7pre$L7dateunique==x,], 
      plot(L8NIR~L7NIR, pch=3,col=alpha('navy',0.2),xlim=c(0.08,0.395), ylim=c(0.08,0.395),tck=-0.03))
      title(main=x,line=-0.85)
}

for (x in L7dates){
  with(L8andL7pre[L8andL7pre$L7dateunique==x,], 
       plot(L8Red~L7Red, pch=3,col=alpha('navy',0.2),xlim=c(0.009,0.165),ylim=c(0.009,0.165),tck=-0.01))
       title(main=x,line=-0.85)
}


# Uncomment to graph giant scatterplots (TAKES A FEW MOMENTS TO RUN)
#par(mfrow=c(3,2),mar=c(2.2,2.1,0.7,0.6),mgp=c(1.0,0.1,0))
#with(L8andL7final,  plot(L8NIR~L7NIR, pch=3,col=alpha('navy',0.2),xlim=c(0.10,0.395),ylim=c(0.10,0.395)))
#  abline(a=0,b=1,lwd=2) ; abline(regNIR_L7asL8,col='tomato3',lwd=2)
#with(L8andL7final,  plot(L8Red~L7Red, pch=3,col=alpha('navy',0.2),xlim=c(0.009,0.165),ylim=c(0.009,0.165)))
#  abline(a=0,b=1,lwd=2) ; abline(regRed_L7asL8,col='tomato3',lwd=2)


#############################################################################
## Landsat 5 and 7
#############################################################################

L5L7DateMatches = read.csv('L5toL7MatchesTable_Bradford.csv', header=TRUE)

## Get each of the L5 and L7 datasets into long data frame form
L5_NIRlong    = melt(L5_NIR,     na.rm=FALSE, value.name='L5NIR',id="pixelID"); names(L5_NIRlong)[2]    ="L5date"
L7for5_NIRlong= melt(L7for5_NIR, na.rm=FALSE, value.name="L7NIR",id="pixelID"); names(L7for5_NIRlong)[2]="L7date"
L5_Redlong    = melt(L5_Red,     na.rm=FALSE, value.name='L5Red',id="pixelID"); names(L5_Redlong)[2]    ="L5date"
L7for5_Redlong= melt(L7for5_Red, na.rm=FALSE, value.name='L7Red',id="pixelID"); names(L7for5_Redlong)[2]="L7date"

# Add a column that gives the L7 date to be matched to each L5 date
L5andL7NIR = merge(L5_NIRlong, L5L7DateMatches,by='L5date',all.x=TRUE)

# Add a column that gives the L7 NIR values (found by linking to the L7 dates just added)
L5andL7NIR = merge(L5andL7NIR, L7for5_NIRlong ,by=c('pixelID','L7date'),all.x=TRUE)

# You now have a table with L5 & L7 NIR values. Now we add L5 & L7 Red values
L5andL7pre    = merge(L5andL7NIR, L5_Redlong    ,by=c('pixelID','L5date'),all.x=TRUE)

L5andL7pre    = merge(L5andL7pre, L7for5_Redlong,by=c('pixelID','L7date'),all.x=TRUE)

# Copy the date frame and begin noise removal on the new data frame
L5andL7 = L5andL7pre[rowSums(!is.na(L5andL7pre))==ncol(L5andL7pre), ]

# For comparison, check regression now, BEFORE eliminating bad dates
regNIR = lm(L7NIR~L5NIR,data=L5andL7); summary(regNIR)$adj.r.squared
regRed = lm(L7Red~L5Red,data=L5andL7); summary(regRed)$adj.r.squared

## Check  whether there are any date pairs that should be excluded because of major differences
## in values, which would suggest a large-scale change happened at the study area in between
## the matched images

# First calculate point-level differences between matched pairs
L5andL7$NIRdiff = abs(L5andL7$L5NIR - L5andL7$L7NIR)
L5andL7$Reddiff = abs(L5andL7$L5Red - L5andL7$L7Red)

# Aggregate point differences by date

L7dates = unique(L5andL7$L7dateunique) ; L7dates

L5andL7NIRdiff_agg = aggregate(L5andL7$NIRdiff, list(L5andL7$L7dateunique), FUN=mean,na.rm=TRUE) 
L5andL7NIRdiff_agg$x = round(L5andL7NIRdiff_agg$x,4)
L5andL7NIRdiff_agg[order(-L5andL7NIRdiff_agg$x),] ## take a look -- any pairs with unusually high differences?

L5andL7Reddiff_agg = aggregate(L5andL7$Reddiff, list(L5andL7$L7dateunique), FUN=mean,na.rm=TRUE) 
L5andL7Reddiff_agg$x = round(L5andL7Reddiff_agg$x,4)
L5andL7Reddiff_agg[order(-L5andL7Reddiff_agg$x),] ## take a look -- any pairs with unusually high differences?

# Also calculate R-squared values within each date pair and compile them
L5L7NIR_Rsqs = L5L7Red_Rsqs = c()
for (x in L7dates){
  rsq = with(L5andL7pre[L5andL7pre$L7dateunique==x,], round(summary(lm(L5NIR~L7NIR))$adj.r.squared,3))
  L5L7NIR_Rsqs[which(L7dates==x)]=rsq
}
for (x in L7dates){
  rsq = with(L5andL7pre[L5andL7pre$L7dateunique==x,], round(summary(lm(L5Red~L7Red))$adj.r.squared,3))
  L5L7Red_Rsqs[which(L7dates==x)]=rsq
}

# Get the R-squared values and average differences for Red and NIR together in a data frame
L5L7Rsqs = data.frame(cbind(L7dates,L5L7NIR_Rsqs,L5L7Red_Rsqs))

L5L7Rsqs = merge(L5L7Rsqs, L5andL7NIRdiff_agg, by.x = 'L7dates', by.y='Group.1')
names(L5L7Rsqs)[length(L5L7Rsqs)]='L5L7NIR_diff'

L5L7Rsqs = merge(L5L7Rsqs, L5andL7Reddiff_agg, by.x = 'L7dates', by.y='Group.1')
names(L5L7Rsqs)[length(L5L7Rsqs)]='L5L7Red_diff'

# Graphing R-squared values and average differences for individual dates
xseq = seq(1,nrow(L5L7Rsqs))
par(mfrow=c(2,1),mar=c(5.2,2.8,0.7,2.9),mgp=c(1.8,0.5,0))
plot(L5L7Rsqs$L5L7NIR_Rsqs ~ xseq,xaxt='n',pch=16,type='b',tck=-0.02,ylab='NIR R-sq',xlab='',ylim=c(0.33,0.96))
axis(1,at=xseq,labels=L5L7Rsqs$L7dates,cex.axis=0.8,las=2,tck=-0.015)
par(new=TRUE);
plot(L5L7Rsqs$L5L7NIR_diff ~ xseq,xaxt='n',yaxt='n',pch=16,type='b',tck=-0.02,col='blue',ylab='',xlab='',ylim=c(0.0032,0.400))
axis(4) ; mtext('NIR diff', side=4,line=1.8)
# L7_20090227 and L7_20101231 look bad

plot(L5L7Rsqs$L5L7Red_Rsqs ~ xseq,xaxt='n',pch=16,type='b',tck=-0.02,ylab='Red R-sq',xlab='',ylim=c(0.33,0.96))
axis(1,at=xseq,labels=L5L7Rsqs$L7dates,cex.axis=0.8,las=2,tck=-0.015)
par(new=TRUE);
plot(L5L7Rsqs$L5L7Red_diff ~ xseq,xaxt='n',yaxt='n',pch=16,type='b',tck=-0.02,col='blue',ylab='',xlab='',ylim=c(0.0032,0.0143))
axis(4) ; mtext('Red diff', side=4,line=1.8)
# bit of judgement call for first couple dates but we'll keep all becausee R-squares are very good
# and the magnitude of overall Red Difference is small (note right axis values)

## Eliminate the bad pairs
L5andL7=subset(L5andL7, subset = L7date != 'L7_20090227' & L7date != 'L7_20101231')
L5andL7=L5andL7[rowSums(!is.na(L5andL7))==ncol(L5andL7), ]

#Regression AFTER eliminating the bad pairs
regNIR = lm(L7NIR~L5NIR,data=L5andL7); summary(regNIR)$adj.r.squared
regRed = lm(L7Red~L5Red,data=L5andL7); summary(regRed)$adj.r.squared

## More noise filtering at the pixel level, based on cases where L8 and L7 found very different values 
## for pixels, which again suggest real changes on the ground in between matched dates, though
## at a smaller scale now. First we use a screening method based on Roy 2016, which only removes 
## extremely misalligned data
L5andL7$NIRnoisecalc = abs(L5andL7$L5NIR-L5andL7$L7NIR) / (0.5*abs(L5andL7$L5NIR+L5andL7$L7NIR))
L5andL7$NIRnoisy     = ifelse(L5andL7$NIRnoisecalc > 1, 1, 0) 

L5andL7$Rednoisecalc = abs(L5andL7$L5Red-L5andL7$L7Red) / (0.5*abs(L5andL7$L5Red+L5andL7$L7Red))
L5andL7$Rednoisy     = ifelse(L5andL7$Rednoisecalc > 1, 1, 0) 

# set "noisy" field to 1 in cases where either band was flagged
L5andL7$noisy = ifelse((L5andL7$NIRnoisy | L5andL7$Rednoisy) == 1, 1, 0)

sum(L5andL7$noisy,na.rm=TRUE) # the portion removed via Roy 2016 adapted method. It will likely be tiny

# Environment clean up
L5andL7 = within(L5andL7, rm(NIRnoisy,Rednoisy,NIRnoisecalc,Rednoisecalc))
rm(L5_NIRlong,L7for5_NIRlong,L5_Redlong,L7for5_Redlong,L5andL7NIR)

## We need to do more here because there are still major disparities that suggest some type
#  of change on the ground between the paired Landsat dates.
#  Additional screening will be based simply on absolute differences.

## To set an unacceptable amount of disparity, calculate mean and st. dev. of the L8-to-L7 differences
#  and we'll say the mean + 3*sd is the cut-off. Anything higher is eliminated.
#  That is, we define an acceptable level of difference between date-matched pixels,
#  and if the difference exceeds that, we assume real change happened on the ground and eliminate it
L5andL7$NIRnoisy2 = ifelse(L5andL7$NIRdiff > 
                             mean(L5andL7$NIRdiff,na.rm=TRUE) + 3*sd(L5andL7$NIRdiff,na.rm=TRUE), 1, 0)
L5andL7$Rednoisy2 = ifelse(L5andL7$Reddiff > 
                             mean(L5andL7$Reddiff,na.rm=TRUE) + 3*sd(L5andL7$Reddiff,na.rm=TRUE), 1, 0)

L5andL7$noisy2 = ifelse((L5andL7$NIRnoisy2 | L5andL7$Rednoisy2) == 1, 1, 0) # Set to 1 if either is high
L8andL7=L8andL7[rowSums(!is.na(L8andL7))==ncol(L8andL7), ]

mean(L5andL7$noisy2,na.rm=TRUE) # the portion removed

# Environment clean up
#L5andL7 = within(L5andL7, rm(NIRdiff,Reddiff,NIRnoisy2,Rednoisy2))

# Regression BEFORE eliminating noisy data
regNIR = lm(L7NIR~L5NIR,data=L5andL7); summary(regNIR)$adj.r.squared
regRed = lm(L7Red~L5Red,data=L5andL7); summary(regRed)$adj.r.squared

# Remove the noisy data
L5andL7final= L5andL7[L5andL7$noisy2!=1,]
L5andL7final=L5andL7final[rowSums(!is.na(L5andL7final))==ncol(L5andL7final), ]
nrow(L5andL7final)

# Regression AFTER eliminating noisy data
regNIR_L5asL7 = lm(L7NIR~L5NIR,data=L5andL7final); summary(regNIR_L5asL7)$adj.r.squared
regRed_L5asL7 = lm(L7Red~L5Red,data=L5andL7final); summary(regRed_L5asL7)$adj.r.squared

# GRAPHING
# Graphing scatterplots for individual dates
L7dates = unique(L5andL7$L7dateunique) ; L7dates
par(mfrow=c(5,7),mar=c(2.2,2.1,0.7,0.6),mgp=c(1.0,0.1,0))
L5L7NIR_Rsqs = L5L7Red_Rsqs = c()

for (x in L7dates){
  with(L5andL7[L5andL7$L7dateunique==x,], 
      plot(L5NIR~L7NIR, pch=3,col=alpha('navy',0.2),xlim=c(0.08,0.395), ylim=c(0.08,0.395),tck=-0.03))
    title(main=x,line=-0.85)
}

for (x in L7dates){
  with(L5andL7[L5andL7$L7dateunique==x,], 
       plot(L5Red~L7Red, pch=3,col=alpha('navy',0.2),xlim=c(0.009,0.165),ylim=c(0.009,0.165),tck=-0.01))
  title(main=x,line=-0.85)
}

# Uncomment to graph giant scatterplots (TAKES A FEW MOMENTS TO RUN)
#par(mfrow=c(1,2),mar=c(2.2,2.1,0.7,0.6),mgp=c(1.0,0.1,0))
#with(L5andL7final,plot(L7NIR~L5NIR, pch=3,col=alpha('navy',0.2),xlim=c(0.10,0.395),ylim=c(0.10,0.395)))
#  abline(a=0,b=1,lwd=2) ; abline(regNIR,col='tomato3',lwd=2)
#with(L5andL7final,plot(L7Red~L5Red, pch=3,col=alpha('navy',0.2),xlim=c(0.008,0.157),ylim=c(0.008,0.157)))
#  abline(a=0,b=1,lwd=2) ; abline(regRed,col='tomato3',lwd=2)

#################################################################################################

## Good so far. Now we have regressions that predicts L8-equivalent values from L7 values.
## We have regressions that predicts L7-equivalent values from L5 values.

## Now we need to apply the L7-as-L8 regression to L7 to derive them as L8-equivalent.
## Then, regress these L8-equivalent L7 values against L5 values, to create a regression equation
## that predicts L8-equivalent values from L5

head(L5andL7final ,2)

# Apply the L8-equivalent-from-L7 equation to L7, deriving a dataset of L8 equivalent values
L5andL7final$L8equivL7NIR = predict(regNIR_L7asL8,  L5andL7final)
L5andL7final$L8equivL7Red = predict(regRed_L7asL8,  L5andL7final)

# Regress the new L8-equivalent values against the L5 values
regNIR_L5asL8 = lm(L8equivL7NIR~L5NIR,data=L5andL7final) ; summary(regNIR)
regRed_L5asL8 = lm(L8equivL7Red~L5Red,data=L5andL7final) ; summary(regRed)

head(L5andL7final)

#Uncomment to run (these take a few moments)
#par(mfrow=c(2,1),mar=c(1.8,2.1,0.5,0.6),mgp=c(1.0,0.1,0))
#plot(L5andL7final$L8equivL7NIR~L5andL7final$L5NIR,pch=3,cex=0.8,col=alpha('dodgerblue3',0.4),xlim=c(0.10,0.395),ylim=c(0.10,0.395),xlab='L5 NIR',ylab='L8 equivalent NIR')
#abline(a=0,b=1,lwd=2,lty='longdash')
#abline(regNIR,col='red',lwd=2)

#plot(L5andL7final$L8equivL7Red~L5andL7final$L5Red,pch=3,cex=0.8,col=alpha('dodgerblue3',0.1),xlim=c(0.008,0.157),ylim=c(0.008,0.157),xlab='L5 Red',ylab='L8 equivalent Red') 
#abline(a=0,b=1,lwd=2,lty='longdash') 
#abline(regRed,col='red',lwd=2)


############################
## Final tegression equations which can be applied to Landsat 5 and 7 spectral reflectance
#  to harmonize them with Landsat 8 spectral reflectance

summary(regNIR_L7asL8)
summary(regRed_L7asL8)
summary(regNIR_L5asL8)
summary(regRed_L5asL8)
