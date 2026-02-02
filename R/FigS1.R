gc();rm(list=ls())

################################################################################
#required packages
list.of.packages <- c(
    "pdftools","plotly","metR",
    "dplyr", "tidyr","here","rlist","ggtext","gridExtra","grid",
    "lattice","viridis","performance","MuMIn","raster","RColorBrewer",
    "rgl","patchwork","cowplot","ggpubr","ggnewscale","plot.matrix",
    "DT","ggthemes","fields") 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

################################################################################
#define the data repository
dir.in=here("Output/full_model") #to change accordingly to the location of the data
dir.out=here("Figs") #to change. It's the repository where the results are saved

# Load data
mydataset <- read.csv2(here("Data","gen_data_final_fonseca.csv"),
                       sep=",",dec=".",h=T) 

#Data selection
## Latitude data
mydatatogo <- mydataset  %>%
    dplyr::filter(Type == "LAT",
                  shift_vel_sign == "pospos" | shift_vel_sign == "negneg", # Select only shifts in the same direction of velocity
                  SHIFT != 0, # remove non-significant shifts
                  # Nucleotide_diversity > 0 # select only GD values > 0
    ) %>% 
    mutate(
        # Climate velocity
        vel = as.numeric(velocity),
        vel_abs = abs(vel),
        vel_abs_log = log(vel_abs),
        vel_abs_log1p = log1p(vel_abs),
        
        # Shift
        SHIFT = SHIFT, # to deal with zero shift
        SHIFT_abs = abs(SHIFT),
        SHIFT_abs_log = log(SHIFT_abs),
        SHIFT_abs_log1p = log1p(SHIFT_abs),
        
        # Genetic diversity
        GD = Nucleotide_diversity, 
        GD_sqrt = sqrt(GD),
        GD_log = log(GD),
        GD_log1p = log1p(GD),
        
        # Methods
        Lat = abs(Latitude),
        Lat_band = round(Lat,0),
        ID.area = scale(ID.area),
        DUR = scale(DUR),
        LogExtent = log(Extent),
        START = scale(START),
        Param = factor(Param),
        Group = factor(Group),
        spp = factor(spp), 
        ExtentF = cut(Extent,
                      ordered_result = TRUE,
                      breaks=seq(min(Extent), max(Extent), length.out=10),
                      include.lowest=TRUE),
        NtempUnitsF = cut(NtempUnits,
                          ordered_result = TRUE,
                          breaks=seq(min(NtempUnits), max(NtempUnits), length.out=10),
                          include.lowest=TRUE)) %>%
    
    dplyr::select(
        # Genetic diversity
        GD, GD_log, GD_log1p, GD_sqrt, TajimasD,
        # Shift
        SHIFT, SHIFT_abs, SHIFT_abs_log, SHIFT_abs_log1p, 
        # SHIFT_cor, SHIFT_cor_abs, SHIFT_cor_abs_log, SHIFT_cor_raw, SHIFT_abs_log_scale,
        # Velocity
        vel, vel_abs, vel_abs_log, vel_abs_log1p, 
        trend.mean,
        # Methods + Taxonomy
        Article_ID, 
        Hemisphere,
        shift_vel_sign,
        Lat, # latitudinal position where GD was collected
        Lat_band,
        DUR, Nperiodes, LogNtempUnits, NtempUnits, Extent, LogExtent, ContinuousGrain, Quality, PrAb, ExtentF, NtempUnitsF,
        Param, Group, spp, Class, Order, Family, Genus, 
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data, Article_ID
    ) 

# transform continuous variables
cont_vars <- c(1:14,18, 20:26)

mydatatogo[,cont_vars] <- lapply(mydatatogo[,cont_vars], as.numeric)
mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

##eliminating species with uncertain obs
mydatatogo=subset(mydatatogo,spp!="Agrilus planipennis")
mydatatogo=subset(mydatatogo,spp!="Chrysodeixis eriosoma")
mydatatogo=droplevels(mydatatogo)

## Filter Classes with at least 5 species per Param
n_sps = 10

test <- mydatatogo %>%
    group_by(Class,Param) %>%
    summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    dplyr::filter(Class_Param %in% test$Class_Param) %>%
    dplyr::select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 


############################Figure S1###########################################

params <- c("ParamTE","ParamCE","ParamLE")
params2 <- c("TE","O","LE")
params3 <- c(":ParamTE","",":ParamLE")
param4 <- c("Trailing edge", "Centroid", "Leading edge")

png(paste0(dir.out,"/FigS1.png"),unit="cm",width=20,height=6.666,res=300)#,width=547,height=360

par(mfrow = c(1,3)) 

for(j in 1:length(unique(params))){
    
    c1=read.csv2(here(dir.in,"summary_coeff.csv"),
                 sep=";",dec=".",h=T)
    c1a=subset(c1,model=="allW")
    
    if(params2[j] == "O"){
        int=c1a$median[c1a$var=="(Intercept)"]
    } else {
        int=c1a$median[c1a$var=="(Intercept)"]+c1a$median[c1a$var==params[j]]
    }
    
    dsel=subset(mydatatogo,Param==params2[j])
    
    dsel$GD2=round(dsel$GD,4)
    dsel$vel_abs2=round(dsel$vel_abs,1)
    kd <- with(dsel, 
               MASS::kde2d(GD2,vel_abs2, 
                           n = 50,
                           lims=c(round(range(mydatatogo$GD),4),
                                  round(range(mydatatogo$vel_abs),1))))
    
    
    v1=seq(round(min(dsel$vel_abs),1),round(max(dsel$vel_abs),1),by=0.1)
    v2=seq(round(min(dsel$GD),4),round(max(dsel$GD),4),length.out=nrow(mydatatogo))  
    for(i in 1:length(kd$y)){
        pred1=exp(int+(c1a$median[c1a$var=="scale(GD)"]+c1a$median[c1a$var==paste0("scale(GD)",params3[j])])*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                      (c1a$median[c1a$var=="scale(vel_abs)"]+c1a$median[c1a$var==paste0("scale(vel_abs)",params3[j])])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))+
                      (c1a$median[c1a$var=="scale(vel_abs):scale(GD)"]+c1a$median[c1a$var==paste0("scale(vel_abs):scale(GD)",params3[j])])*((kd$y[i]-mean(mydatatogo$vel_abs))/sd(mydatatogo$vel_abs))*((kd$x-mean(mydatatogo$GD))/sd(mydatatogo$GD))+
                      c1a$median[c1a$var=="LogExtent"]*mean(mydatatogo$LogExtent)+c1a$median[c1a$var=="LogNtempUnits"]*mean(mydatatogo$LogNtempUnits)+c1a$median[c1a$var=="ContinuousGrain"]*2) 
        pred1=data.frame(GD=kd$x,pred1,
                         vel_abs=kd$y[i],freq=kd$z[,i])
        if(i==1){
            pred2=pred1
        }else{
            pred2=rbind(pred2,pred1)
        }
    }
    
    
    pred2=pred2[order(pred2$vel_abs),]
    pred2=pred2[order(pred2$GD),]
    z1=matrix(pred2$freq,ncol=length(unique(pred2$GD)),byrow=T)
    jet.colors2 <- colorRampPalette( c("red","orange","yellow","green") ) 
    pal2 <- jet.colors2(100)
    z=z1
    z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
    col.ind2 <- cut(z.facet.center, 100)
    col.ind2[z.facet.center<0.01]="#FFFFFF"
    s1=1:round(max(z.facet.center),0)
    col.ind3 <- cut(s1, 100)
    
    par(mar=c(0,2,2,4),
        mgp = c(8, 10, 0))
    
    with(pred2,
         persp(x=sort(unique(GD)),
               y=sort(unique(vel_abs)),
               z=matrix(pred1,ncol=length(unique(GD)),byrow=T),
               col=pal2[col.ind2],
               main = param4[j],
               ticktype = "detailed", phi=30, theta=-40,border="grey80",lwd=0.2,
               xlab='\n\nGenetic diversity',
               ylab='\n\nClimate change velocity',
               zlab='\n\nPredicted\nrange shift velocity',
               cex.lab=0.75,
               cex.axis=0.75))
    
    r1=rasterFromXYZ(pred2[,c("GD","vel_abs","pred1")])
    rX2=r1
    rX2[]=sample(1:100,size=ncell(r1),rep=T)
    
    plot(rX2, horizontal=F, legend.only=TRUE, col=pal2,
         smallplot=c(.85, .87, .4, .8), maxpixels=2000000,
         axis.args=list(tck=-0.5,
                        at=as.numeric(col.ind3[s1%in%c(1,seq(20,100,by=20))]),
                        labels=F),
         legend.args=list(text="Density of observations", 
                          side=4,font=2, 
                          line=1.2, 
                          cex=0.65))
    
    if(params2[j] == "TE"){
        plot(rX2, horizontal=F, legend.only=TRUE, col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(20,100,by=20))]),
                            labels=c(1,seq(20,100,by=20)),
                            cex.axis=0.5))
    }
    if(params2[j] == "O"){
        plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(20,60,by=20),max(s1))]),
                            labels=c(1,seq(20,60,by=20),max(s1)),
                            cex.axis=0.5))
    }
    if(params2[j] == "LE"){
        plot(rX2, horizontal=F, legend.only=TRUE,col=pal2,
             smallplot=c(.85, .87, .4, .8),maxpixels=2000000,
             axis.args=list(tck=F,
                            lwd=0,
                            line=-0.50,
                            at=as.numeric(col.ind3[s1%in%c(1,seq(25,175,by=25),max(s1))]),
                            labels=c(1,seq(25,175,by=25),max(s1)),
                            cex.axis=0.5))
    }
    
}



dev.off()



