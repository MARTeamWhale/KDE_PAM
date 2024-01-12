# creat KDE maps similar to ArcGIS...
options(scipen = 999)
library(stars)
library(mapsf)
library(stringr)

#NBW Map KDE --------
Bb_KDE = rast(here::here("output/Bb_sf_kernel_density.tif"))
plot(Bb_KDE)

#figure out quantiles excluding 0
Bb_KDE = st_as_stars(Bb_KDE)%>%rename(KDE = Bb_sf_kernel_density.tif)%>%mutate(KDE = ifelse(KDE == 0, NA, KDE))
breaks_quant = mf_get_breaks(x = c(min(Bb_KDE$KDE, na.rm = T),
                                   Bb_KDE$KDE),  breaks = "q6")
breaks_quant = c(breaks_quant[c(2,3,4,5,6,7)])

#do something funky to create quantile classes with NA values
Bb_KDE_quant <- mutate(Bb_KDE, breaks = cut(KDE, breaks_quant)) 
Bb_KDE_quant = rast(Bb_KDE_quant[2])

plot((Bb_KDE_quant))

# raster to polygon

Bb_KDEQ_ImHab = st_make_valid(st_as_sf(terra::as.polygons(Bb_KDE_quant), merge = T))

plot(st_geometry(Bb_KDEQ_ImHab))

#get rid of weird characters and make clean orderd factor for levels without hard coding
Bb_KDEQ_ImHab = Bb_KDEQ_ImHab%>%mutate(quants = str_replace_all(breaks, "[\\(\\]]", "" ))%>%
  mutate(quants = as.numeric(gsub(",.*", "", quants)))%>%mutate(quants = factor(quants))

Bb_KDEQ_ImHab$quants

Bb_KDEQ_ImHab = Bb_KDEQ_ImHab%>%mutate(Quantile = c(0.05, 0.25, 0.50, 0.75, 0.95))

write_sf(Bb_KDEQ_ImHab, "output/KDE_NBW_Quant.shp")

#only 50% and 95%
breaks_quant = mf_get_breaks(x = c(min(Bb_KDE$KDE, na.rm = T),
                                   Bb_KDE$KDE),  breaks = "q6")
breaks_5090 = c(breaks_quant[c(4,6,7)])

#do something funky to create quantile classes with NA values
Bb_KDE_5090 <- mutate(Bb_KDE, breaks = cut(KDE, breaks_5090)) 
Bb_KDE_5090 = rast(Bb_KDE_5090[2])

plot((Bb_KDE_5090))

# raster to polygon

Bb_KDE_5090_sf = st_make_valid(st_as_sf(terra::as.polygons(Bb_KDE_5090), merge = T))

plot(st_geometry(Bb_KDE_5090_sf))

#get rid of weird characters and make clean orderd factor for levels without hard coding
Bb_KDE_5090_sf = Bb_KDE_5090_sf%>%mutate(quants = str_replace_all(lyr.1, "[\\(\\]]", "" ))%>%
  mutate(quants = as.numeric(gsub(",.*", "", quants)))%>%mutate(quants = factor(quants))

Bb_KDE_5090_sf$quants

Bb_KDE_5090_sf = Bb_KDE_5090_sf%>%mutate(Quantile = c(0.50, 0.95))

write_sf(Bb_KDE_5090_sf, "Shapes/KDE/Bb_KDE_5090.shp")


#SBW Map KDE --------
SBW_KDE = rast(here::here("Shapes/KDE/KDE_SBW2.tif"))
plot(SBW_KDE)

#figure out quantiles excluding 0
SBW_KDE = st_as_stars(SBW_KDE)%>%rename(KDE = KDE_SBW2.tif)%>%mutate(KDE = ifelse(KDE == 0, NA, KDE))
breaks_quant = mf_get_breaks(x = c(min(SBW_KDE$KDE, na.rm = T),
                                   SBW_KDE$KDE),  breaks = "q6")
breaks_quant = c(breaks_quant[c(2,3,4,5,6,7)])

#do something funky to create quantile classes with NA values
SBW_KDE_quant <- mutate(SBW_KDE, breaks = cut(KDE, breaks_quant)) 
SBW_KDE_quant = rast(SBW_KDE_quant[2])

plot((SBW_KDE_quant))

# raster to polygon

SBW_KDEQ_ImHab = st_make_valid(st_as_sf(terra::as.polygons(SBW_KDE_quant), merge = T))

plot(st_geometry(SBW_KDEQ_ImHab))

#get rid of weird characters and make clean orderd factor for levels without hard coding
SBW_KDEQ_ImHab = SBW_KDEQ_ImHab%>%mutate(quants = str_replace_all(lyr.1, "[\\(\\]]", "" ))%>%
  mutate(quants = as.numeric(gsub(",.*", "", quants)))%>%mutate(quants = factor(quants))

SBW_KDEQ_ImHab$quants

SBW_KDEQ_ImHab = SBW_KDEQ_ImHab%>%mutate(Quantile = c(0.05, 0.25, 0.50, 0.75, 0.95))

write_sf(SBW_KDEQ_ImHab, "Shapes/KDE/KDE_SBW_Quant.shp")

#only 50% and 95%
breaks_quant = mf_get_breaks(x = c(min(SBW_KDE$KDE, na.rm = T),
                                   SBW_KDE$KDE),  breaks = "q6")
breaks_5090 = c(breaks_quant[c(4,6,7)])

#do something funky to create quantile classes with NA values
SBW_KDE_5090 <- mutate(SBW_KDE, breaks = cut(KDE, breaks_5090)) 
SBW_KDE_5090 = rast(SBW_KDE_5090[2])

plot((SBW_KDE_5090))

# raster to polygon

SBW_KDE_5090_sf = st_make_valid(st_as_sf(terra::as.polygons(SBW_KDE_5090), merge = T))

plot(st_geometry(SBW_KDE_5090_sf))

#get rid of weird characters and make clean orderd factor for levels without hard coding
SBW_KDE_5090_sf = SBW_KDE_5090_sf%>%mutate(quants = str_replace_all(lyr.1, "[\\(\\]]", "" ))%>%
  mutate(quants = as.numeric(gsub(",.*", "", quants)))%>%mutate(quants = factor(quants))

SBW_KDE_5090_sf$quants

SBW_KDE_5090_sf = SBW_KDE_5090_sf%>%mutate(Quantile = c(0.50, 0.95))

write_sf(SBW_KDE_5090_sf, "Shapes/KDE/SBW_KDE_5090.shp")




#PLOT---------

#pallet
predpal="mako"
pal = rev(get(predpal)(7))
show_col(pal)
pal = c(pal[c(2,3,4,5,6)])

# Get the vertical and horizontal limits
ext <- extent( Bb_KDEQ_ImHab )
# Get x and y limits
lims <- list( x=c(ext@xmin-100, ext@xmax+200000), y=c(ext@ymin-100, ext@ymax+10000) )

#map it
m1 <-ggplot() +
  theme_bw()+
  
  # add land region--
  geom_sf(  data = landUTM, color=NA, fill="grey50") +
  # add contours--
  geom_sf(data = cont_UTM%>%dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), col = "grey50", linewidth  = 0.2) +
  # 
  # #KDE - with quantiles
  geom_sf(data  = SBW_KDEQ_ImHab, aes(fill = quants), col = NA,  alpha = .9, na.rm = T)+
  
  # SBW  PAM AMAR detects
  geom_sf(data = BW_detectsUTM%>%filter( source == "AMAR"), col = "black", shape = 24,
          fill = "yellow", size = 1.5, alpha = .5) +
  
  # set map limits
  coord_sf(lims_method = "orthogonal",
           xlim=lims$x, ylim=lims$y, expand = F
  )+
  
  
  # format axes
  ylab("") + 
  xlab("") +
  
  #theme KDE
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position= c(0.1, 0.9) , legend.background = element_rect(fill = NA), 
    legend.key = element_rect(fill = NA), plot.margin = margin(0, 0, 0, 0), 
    plot.title = element_text(hjust=0.5), 
    plot.subtitle  = element_text(hjust=0.5), axis.title = element_blank())+
  
  scale_fill_manual(values = (pal),
                    na.value = NA, na.translate = F,
                    labels=c("0.05","0.25", "0.50","0.75","0.95" ) ,
                    name = "",
                    guide_legend(title.position = "bottom") )+
  # add scale bar
  annotation_scale(
    location = "br",
    width_hint = 0.25,
    text_cex = 0.6,
    bar_cols = c("grey40", "white")
  ) 

m1


ggsave(here::here("FIGS/Fig6h_SBW_KDE.png"), m1, width = 11, height = 8.5, units = "in")