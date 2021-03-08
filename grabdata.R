pacman::p_load(devtools)
#install_github("macvt/BSE4304/BSEHyrdoModels")


# Installing the packages we will play with today
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,raster,soilDB,rgdal)
pacman::p_load(EcoHydRology,curl,httr,rnoaa)

#######Flowgage#############
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id, begin_date="2015-01-01",end_date="2021-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2020,]$id[2]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)


# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# Compare your precipitation to the flow out of your basin
mean(modeldata$Qmm)
mean(modeldata$P)
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0
summary(modeldata)


modeldata$HillslopeAboveExcess=0
summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

TopSlope = TMWB_Model(fnc_TMWB = TopSlope,func_DAWC =0.3, 
                      func_z=500,fnc_fcres=.3)
MidSlope$HillslopeAboveExcess=TopSlope$Excess
# Higher slope, medium ksat, fcres=0.5 
MidSlope = TMWB_Model(fnc_TMWB = MidSlope,func_DAWC =0.3, 
                      func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$HillslopeAboveExcess=MidSlope$Excess
BotSlope = TMWB_Model(fnc_TMWB = BotSlope,func_DAWC=.3,
                      func_z=1000,fnc_fcres=.2)

#AW Plots HW1
plot(BotSlope$date,BotSlope$AW,type="l",col=1,xlab="Date",ylab="AW (mm)")
lines(MidSlope$date,MidSlope$AW,type="l",col=2)
lines(TopSlope$date,TopSlope$AW,type="l",col=3)
# Excess Plots HW1
plot(BotSlope$date,BotSlope$Excess,type="l",col=1,xlab="Date",ylab="Excess (mm)")
lines(MidSlope$date,MidSlope$Excess,type="l",col=2)
lines(TopSlope$date,TopSlope$Excess,type="l",col=3)



plot(BotSlope$date,BotSlope$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlope$date,BotSlope$ET,type="l",col=2)
lines(MidSlope$date,MidSlope$ET,type="l",col=3)
lines(TopSlope$date,TopSlope$ET,type="l",col=4)

# Fun Cumulative Summary
plot(BotSlope$date,cumsum(BotSlope$ET),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlope$date,cumsum(MidSlope$ET),col="red")
lines(TopSlope$date,cumsum(TopSlope$ET),col="green")

plot(BotSlope$date,BotSlope$Qpred,type="l")
lines(BotSlope$date,BotSlope$Qmm, col=2)
NSeff(BotSlope$Qmm,BotSlope$Qpred)


####Graduate Student
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
download.file(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
mystream=subset(streams,GNIS_ID=="01478950")
mybbox=c(mystream@bbox)
mysoil = mapunit_geom_by_ll_bbox(mybbox)
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
# THE ONLY DIFFERENCE BETWEEN LAST WEEKS LAB AND THE HW3 IS:
# We modify these following lines from our Lab4 exercise to access 
# the "resdept" from the "corestrictions" table as described in the PDF.
# So, from last week we modify this line:
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
# To HW3 needs (here we are using "cr" for corestrictions where before 
# we used "ch" for "chorizon"... naming should reflect the data (though 
# many would use the full name so they can remember when read their code)
q_co2cr = paste("SELECT cokey,resdept_r FROM corestrictions WHERE cokey IN ", cokey_statement, sep="")
co2cr = SDA_query(q_co2cr)
View(co2cr)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2cr=merge(mu2co,co2cr)
mu2crmax=aggregate(mu2cr,list(mu2cr$mukey),max)
depth2restrict<-merge(mysoil,mu2crmax,by="mukey")
# And then a simple plot will do...
plot(mystream,col="red",lwd=4)
plot(depth2restrict,col=topo.colors(10),add=T)
plot(mystream,col="red",lwd=4,add=T)
# Though you want to explore nicer plots
pols1 <- list("sp.lines", as(mystream, 'SpatialLines'), col = "red", lwd = 4)
spplot(depth2restrict, sp.layout=list(pols1), zcol="resdept_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=colorRampPalette(c("grey87","royalblue4"))(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])
# or other colormap
spplot(depth2restrict, sp.layout=list(pols1), zcol="resdept_r", xlab="Longitude", ylab="Latitude", main="Depth to Restrictive Layer", colorkey=T, col.regions=topo.colors(100), checkEmptyRC=T, add=T,xlim=mystream@bbox[1,],ylim=mystream@bbox[2,])



modeldata$HillslopeAboveExcess=0
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata


TopSlopeCN=CNModel(TopSlopeCN, CNavg = 64)
TopSlopeCN = CNModel(fnc_CNModel = TopSlopeCN, CNavg = 64,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlopeCN$HillslopeAboveExcess=TopSlopeCN$Excess
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNModel(fnc_CNModel = MidSlopeCN, CNavg = 90,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlopeCN$HillslopeAboveExcess=MidSlopeCN$Excess
BotSlopeCN = CNModel(fnc_CNModel = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)


#AW Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$AW,type="l",col= "Black", xlab="Date",ylab="AW (mm)")
lines(MidSlopeCN$date,MidSlopeCN$AW,type="l",col= "Red")
lines(TopSlopeCN$date,TopSlopeCN$AW,type="l",col= "Green")
legend("bottomright", legend=c("Bottom Slope", "Medium Slope" , "Top Slope"),
       col=c("Black", "Red" , "Green"), lty=1:2, cex=0.6)
# Excess Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$Excess,type="l",col= "Black",xlab="Date",ylab="Excess (mm)")
lines(MidSlopeCN$date,MidSlopeCN$Excess,type="l",col="Red")
lines(TopSlopeCN$date,TopSlopeCN$Excess,type="l",col="Green")
legend("topleft", legend=c("Bottom Slope", "Medium Slope" , "Top Slope"),
       col=c("Black", "Red" , "Green"), lty=1:2, cex=0.6)
# PET and ET HW2
plot(BotSlopeCN$date,BotSlopeCN$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlopeCN$date,BotSlopeCN$ET,type="l",col="Red") 
lines(MidSlopeCN$date,MidSlopeCN$ET,type="l",col="Green")
lines(TopSlopeCN$date,TopSlopeCN$ET,type="l",col="Blue")
legend("topleft", legend=c("PET", "Bottom Slope", "Medium Slope" , "Top Slope"),
       col=c("Black","Red" , "Green","Blue"), lty=1:2, cex=0.6)
# or as cumulative summations
plot(TopSlopeCN$date,cumsum(BotSlopeCN$PET),type="l",
     xlab="Date",ylab="(P)ET")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$ET),col="red")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$ET),col="green")
lines(BotSlopeCN$date,cumsum(BotSlopeCN$ET),col="blue")


# Cumulative Summary of QPred is very informative
plot(BotSlopeCN$date,cumsum(BotSlopeCN$Excess),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$Excess),col="red")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$Excess),col="green")

# Model Performance 
plot(BotSlopeCN$date,BotSlopeCN$Qpred,type="l")
##lines(BotSlopeCN$)
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

##############HW2
# LickRun VA

# Plot the flow gages measured flow against each of the Top,Mid,Bot components

plot(myflowgage$flowdata$mdate,myflowgage$flowdata$Qmm,type="l",main="Predicted vs observed discharge, Lick Run VA",
     sub= "Black = Observed, Red = TMWB, Green = CN",
     xlab="Date",ylab="Flow (mm)")

lines(BotSlope$date,BotSlope$Qpred,col="red")

lines(BotSlopeCN$date,BotSlopeCN$Qpred,col="green")
NSE=function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs,
                                                         na.rm=TRUE))^2, na.rm=TRUE))}
print(NSE(myflowgage$flowdata$Qmm,BotSlope$Qpred))
print(NSE(myflowgage$flowdata$Qmm,BotSlopeCN$Qpred))
