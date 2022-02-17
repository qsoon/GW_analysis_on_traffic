library(tidyverse)
library(sf)
library(sp)
library(rgdal)
library(GWmodel)
library(RColorBrewer)
library(car)

# read shp file of Seoul
seoul_shp = st_read("/home/kyu9510/courseworks/spatialstat/
                    LARD_ADM_SECT_SGG_seoul/LARD_ADM_SECT_SGG_11.shp")
seoul_shp$SGG_NM = iconv(seoul_shp$SGG_NM, from = "CP949", to = "UTF-8",
                         sub = NA, mark = TRUE, toRaw = FALSE)

seoul_shp$SGG_NM[c(15,16,18)] = c("노원구", "도봉구", "성북구")

TAdata = as_Spatial(seoul_shp) # convert into SpatialPolygonsDataFrame

# read TA data
data = read.csv("/home/kyu9510/courseworks/spatialstat/TAdata.csv", 
                stringsAsFactors=FALSE) # 2020년 

# preprocessing
data = sapply(data[,2:ncol(data)], 
              function(x) as.numeric(gsub(",", "", x)))[nrow(data):2,]

TAdata@data = cbind(TAdata@data, data)

TAdata@data$TANUM = TAdata@data$TANUM/100
TAdata@data$POPDEN = TAdata@data$POPDEN/1000
TAdata@data$PLAREA = TAdata@data$PLAREA/10000
TAdata@data$CAPROP = 100*(TAdata@data$AGE0_9 + TAdata@data$AGE65_) / 
  (TAdata@data$AGE10_24 + TAdata@data$AGE25_44 + TAdata@data$AGE45_64)
TAdata@data$CARNUM = TAdata@data$CARNUM/10000

TAdata@data$X = unlist(lapply(TAdata@polygons, function(x) x@labpt[1]))
TAdata@data$Y = unlist(lapply(TAdata@polygons, function(x) x@labpt[2]))


### variable explanation ###
# TANUM : 교통사고발생건수 (100회)
# POPDEN : 인구밀도 (1000명 / km^2)
# PLAREA : 주차장 면적 (0.01km^2)
# CAPROP : 어린이노약자 비율 (%) 
# CROSSWALK2017 : 2017년 기준 도로 1km 당 횡단보도 수
# ROAD_RATIO : 도로율, 도로면적/시가화면적
# (행정구역면적에서 공원, 임야,하천 등을 제외한 면적) (%)
# CARNUM : 차량등록대수 (10000회)

# GW summary statistics
gw.ss.bs = gwss(TAdata, vars = c("TANUM","POPDEN","PLAREA","CAPROP",
                                 "CROSSWALK2017","ROAD_RATIO","CARNUM"),
                kernel = "bisquare",adaptive=TRUE,bw=5)

mypalette1 = brewer.pal(5, "YlOrRd")

spplot(gw.ss.bs$SDF, "TANUM_LM", key.space="right", col.regions=mypalette1, 
       cuts=4, main="GW mean for TANUM",
       colorkey = list(
         labels=list(
           at = c(seq(10,30,5),mean(TAdata@data[,"TANUM"])),
           labels = c(seq(10,30,5),"Global")
         )
       ))


mypalette2 <-brewer.pal(5, "YlGnBu")
spplot(gw.ss.bs$SDF, "POPDEN_LSD", key.space="right", col.regions=mypalette2, 
       cuts=4, main="GW standard deviation for POPDEN",
       colorkey = list(
         labels=list(
           at = c(seq(1,5,1),sd(TAdata@data[,"POPDEN"])),
           labels = c(seq(1,5,1),"Global")
         )
       ))


# Global regression
lm.global = lm(TANUM ~ POPDEN + PLAREA + CAPROP +
                 CROSSWALK2017 + ROAD_RATIO + CARNUM, data = TAdata)
summary(lm.global)
vif(lm.global) # identify multicollinearity

# remove variables with multicollinearity
lm.global = lm(TANUM ~ POPDEN + CAPROP +
                 CROSSWALK2017 + ROAD_RATIO + CARNUM, data = TAdata) 
summary(lm.global)
vif(lm.global)
lm.global.summary = summary(lm.global)

# GW regression
var.y = "TANUM"
var.x = c("POPDEN", "CAPROP", "CROSSWALK2017", "ROAD_RATIO", "CARNUM")

# select model
model.sel = model.selection.gwr(var.y, var.x, data=TAdata, 
                                kernel="bisquare", adaptive=TRUE, bw=8)

sorted.models = model.sort.gwr(model.sel, numVars = length(var.x),
                               ruler.vector = model.sel[[2]][,2])

model.list = sorted.models[[1]]
model.view.gwr(var.y, var.x, model.list = model.list)

plot(sorted.models[[2]][,2], col = "black", pch = 20, lty = 5,
     main = "Alternative view of GWR model selection procedure", ylab = "AICc",
     xlab = "Model number", type = "b")

bw.gwr.TA <- bw.gwr(TANUM ~ POPDEN + CAPROP + CROSSWALK2017 + 
                      ROAD_RATIO + CARNUM, data = TAdata, approach = "AICc",
                    kernel = "bisquare", adaptive = TRUE)

# implement GW regression
gwr.res = gwr.basic(TANUM ~ POPDEN + CAPROP + CROSSWALK2017 + 
                      ROAD_RATIO + CARNUM, data = TAdata,
                    bw = bw.gwr.TA, kernel = "bisquare", adaptive = TRUE)

print(gwr.res)
names(gwr.res$SDF)

# map GW estimates
mypalette3 <-brewer.pal(5, "Spectral")

spplot(gwr.res$SDF, "POPDEN", key.space = "right",
       col.regions = mypalette3, cuts=4,
       main = "GW regression coefficient estimates for POPDEN",
       sp.layout=map.layout)

spplot(gwr.res$SDF, "CAPROP", key.space = "right",
       col.regions = mypalette3, cuts=4,
       main = "GW regression coefficient estimates for CAPROP",
       sp.layout=map.layout)

spplot(gwr.res$SDF, "CROSSWALK2017", key.space = "right",
       col.regions = mypalette3, cuts=4,
       main = "GW regression coefficient estimates for CROSSWALK2017",
       sp.layout=map.layout)

spplot(gwr.res$SDF, "ROAD_RATIO", key.space = "right",
       col.regions = mypalette3, cuts=4,
       main = "GW regression coefficient estimates for ROAD_RATIO",
       sp.layout=map.layout)

spplot(gwr.res$SDF, "CARNUM", key.space = "right",
       col.regions = mypalette3, cuts=4,
       main = "GW regression coefficient estimates for CARNUM",
       sp.layout=map.layout)


# global PCA
TAdata.scaled = scale(TAdata@data[,c("POPDEN","CAPROP",
                                     "CROSSWALK2017",
                                     "ROAD_RATIO","CARNUM")])

TAdata.scaled = as.data.frame(TAdata.scaled)

gblpca = princomp(TAdata.scaled, cor = FALSE)

gblpca = princomp(TAdata@data[,c("POPDEN","CAPROP","CROSSWALK2017",
                                 "ROAD_RATIO","CARNUM")], cor = FALSE)

gblpca.result = rbind(gblpca$sdev^2, 
                      (gblpca$sdev^2 / sum(gblpca$sdev^2))*100, 
                      unclass(gblpca$loadings)) 

rownames(gblpca.result)[c(1,2)] = 
  c("Eigenvalues","Percentage of total variation")
colnames(gblpca.result) = paste("PC", c(1:5), sep="")

gblpca.result

# bandwidth selection
coords = as.matrix(cbind(TAdata$X, TAdata$Y)) # coordinate matrix
TA.scaled.spdf = SpatialPointsDataFrame(coords, TAdata.scaled)

bw.gwpca.basic = bw.gwpca(TA.scaled.spdf, 
                          vars = colnames(TAdata[,c("POPDEN","CAPROP",
                                                    "CROSSWALK2017",
                                                    "ROAD_RATIO",
                                                    "CARNUM")]@data), k = 5, 
                          robust = FALSE, adaptive = TRUE) 


bw.gwpca.basic = bw.gwpca(TAdata[,c("POPDEN","CAPROP","CROSSWALK2017",
                                    "ROAD_RATIO","CARNUM")], 
                          vars = colnames(TAdata[,c("POPDEN","CAPROP",
                                                    "CROSSWALK2017",
                                                    "ROAD_RATIO",
                                                    "CARNUM")]@data), k = 5, 
                          robust = FALSE, adaptive = TRUE) 

# adaptive=TRUE : include k nearest neighbors

bw.gwpca.basic # selected = 16, using default kernel = bisquare

# GWPCA
gwpca.basic = gwpca(TA.scaled.spdf, 
                    vars = colnames(TAdata[,c("POPDEN","CAPROP",
                                              "CROSSWALK2017",
                                              "ROAD_RATIO",
                                              "CARNUM")]@data),
                    bw = bw.gwpca.basic, k = 5, 
                    robust = FALSE, adaptive = TRUE) 

gwpca.basic = gwpca(TAdata[,c("POPDEN","CAPROP","CROSSWALK2017",
                              "ROAD_RATIO","CARNUM")], 
                    vars = colnames(TAdata[,c("POPDEN","CAPROP",
                                              "CROSSWALK2017",
                                              "ROAD_RATIO",
                                              "CARNUM")]@data),
                    bw = bw.gwpca.basic, k = 5, 
                    robust = FALSE, adaptive = TRUE) 
# k = the number of retained components

prop.var = function(gwpca.obj, n.components) { # calculate var proportion
  return((rowSums(gwpca.obj$var[, 1:n.components]) / 
            rowSums(gwpca.obj$var)) * 100)
}

# 25 location, first 3 PCs'var proportion
var.gwpca.basic = prop.var(gwpca.basic, 3) 

TAdata$var.gwpca.basic = var.gwpca.basic 


mypalette4 <-brewer.pal(6, "RdYlBu") # for gwpca PC 1~3 var.prop
# plot proportion of total variancce for first three local components
# add global
spplot(TAdata, "var.gwpca.basic", key.space = "right", 
       col.regions = mypalette4, cuts = 5, 
       main = "PTV for local components 1 to 3 (GWPCA)", 
       colorkey = list(
         labels=list(
           at = c(seq(88,97,2),91.4),
           labels = c(seq(88,96,2),"Global")
         )
       ))

# 25*5 matrix, 1st PC loadings on each location
loadings.pc1.gwpca.basic = gwpca.basic$loadings[,,1] 

# 25-vector, winning variable for pc1 on each location
win.item1.gwpca.basic = max.col(abs(loadings.pc1.gwpca.basic)) 

TAdata$win.item1.gwpca.basic = win.item1.gwpca.basic 

mypalette5 <-brewer.pal(5, "Paired")

# plot winning variable for local PC1 
# add global  
spplot(TAdata, "win.item1.gwpca.basic", key.space = "right", 
       key=list(x = 0.15, y = 0.93, corner = c(0.5, 0.5),title="",
                points = list(pch=16, col=c(mypalette5, mypalette5[5]), cex=1.5),
                text = list(c("POPDEN","CAPROP","CROSSWALK2017",
                              "ROAD_RATIO","CARNUM","Global"), cex=1.5)),
       col.regions = mypalette5, at = c(1, 2, 3, 4, 5, 6),
       main = "Winning variable: highest abs. loading on local Comp.1 (basic)",
       colorkey = FALSE)


# glyph plot
plot(TAdata, main = "Multivariate glyphs of loadings for local PC1")
coords = as.matrix(cbind(TAdata$X, TAdata$Y))
gwpca.glyph.plot(loadings.pc1.gwpca.basic, coords, add=TRUE, r1=20)
# add global glyph plot
gwpca.glyph.plot(rbind(unclass(gblpca$loadings)[,1],0), rbind(
  c(min(coords[,1])+5000,max(coords[,2])),0), add=TRUE, r1=200) 
text(min(coords[,1])+3000,max(coords[,2])-4000,"Global")


