library(tidyverse)
library(dplyr)
library(ggplot2)
library(rstudioapi)
library(lubridate)

# set working directory ---------------------------------------------------
wd <- dirname(getSourceEditorContext()$path)
filedirectory <- paste0(dirname(wd),"/Files/Outputs")
setwd(filedirectory)
#read sim file
files_read = list.files(filedirectory, pattern = ".csv", full.names = T)
#use future
library(future.apply)
plan(multisession, gc = TRUE)
temp_data = future_lapply(files_read, read_csv, show_col_types = FALSE)
require(DescTools)
big_data = DoCall(rbind, temp_data)
# if (exists("big_data") == T){
#   rm(temp_data)
#   gc()
# }
big_data$Date <- as.Date(big_data$Date, format = "%m/%d/%Y")
big_data <- big_data %>% select(-NUTS2)
str(big_data)
#add dew point
require(weathermetrics)
big_data$DewPoint <- humidity.to.dewpoint(rh = big_data$RHmean, t = big_data$Tmean, temperature.metric = "celsius")
#add VPD 
require(plantecophys)
big_data$VPD <- RHtoVPD(big_data$RHmean, big_data$Tmean, Pa = 101)
#add Year
big_data$Year <- year(big_data$Date) 
#compute NUTS2
# library(tidygeocoder)
# revgeo <- big_data[,c("Lat","Lon")] %>%
#   reverse_geocode(lat = Lat, long = Lon, method = 'osm',return_coords = F,
#                   address = address_found, full_results = F) %>% 
#   as_tibble()
# write_csv(revgeo, "reverse_geo.csv")
# revgeo <- read_csv("C:/Users/loren/Documents/Tesi_master_AP/reverse_geo.csv")
# 
# library(stringr)
# NUTS <- as_tibble(str_split_fixed(revgeo$address_found, ", ", 10))

coord <- unique(big_data[,c("Lat","Lon")])
#map of regions
require(terra)
italy = raster::getData('GADM', country = 'Italy', level = 1, type = "sp")
region_names <-  as.data.frame(terra::extract(italy, coord[,c("Lon","Lat")]))
region_names <- toupper(region_names$NAME_1)
region_names <- cbind(coord,NUTS2=region_names)
#check missing NUTS2
region_names %>% filter(is.na(NUTS2))
big_data %>% filter(Lat == 44.25 & Lon == 7) %>% select(NUTS3) %>% mutate(NUTS3 =as.factor(NUTS3)) %>% sapply(levels)
#add missing region
region_names$NUTS2[is.na(region_names$NUTS2)] <- "PIEMONTE"
#abbreviate ragion names
region_names$NUTS2 <- as.factor(region_names$NUTS2)
levels(region_names$NUTS2) <- c("ABR","PUG","BAS","CAL","CAM","EMR","FVG","LAZ","LIG","LOM","MAR","MOL","PIE","SAR","SIC","TOS","TAD","UMB","VDA","VEN")

#add elevation
require(elevatr)
elev <- get_elev_point(region_names[c("Lon","Lat")], prj = "EPSG:4326", src = "aws",z = 5, override_size_check = TRUE)
region_names$Elev <- elev$elevation
#add NUTS2 and elevation to dataframe
big_data <- left_join(big_data,region_names, by = c("Lat", "Lon"))
#compute basic infection risk 

big_data <- big_data %>% rowwise() %>% mutate(InfRisk = (Rule310+Epi+Ipi+Dmcast+
                                                           Magarey+Rossi+misfits+
                                                           laore)/8) %>% ungroup()

# select data from the infective period -----------------------------------

data_season <- big_data %>% mutate(Month = month(Date)) %>% filter(Month %in% c(3:9))

# if (exists("data_season") == T){
#   rm(big_data)
#   gc()
# }

# compute distance metrics for the binary output of the models ------------

#multicore core
mod_bd <- data_season %>% dplyr::select(NUTS2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,misfits,laore) %>% group_by(NUTS2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "faith"))

hc <- hclust(mod_bd[[1]])
plot(hc)

#multicore
require(foreach)
require(doParallel)
temp2 <- data_season %>% dplyr::select(NUTS2,Year,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,misfits,laore) %>% 
  split(list(.$NUTS2))

cl<-makeCluster(10)
registerDoParallel(cl)

by.chunk <- foreach(df = temp2) %dopar% {
  parallelDist::parDist(as.matrix(t(df[,3:10])),method = "faith")
}

hc2 <- hclust(by.chunk[[1]])
plot(hc2)

# compute model summary variables ---------------------------------------

summary_models <- data_season %>% group_by(Lat,Lon,Year) %>% 
  mutate(
  Sum310 = cumsum(Rule310),
  SumEpi = cumsum(Epi),
  SumIpi = cumsum(Ipi),
  SumDmc = cumsum(Dmcast),
  SumMag = cumsum(Magarey),
  SumRos = cumsum(Rossi),
  SumMis = cumsum(misfits),
  SumLao = cumsum(laore),
  InfRisk = cumsum(InfRisk)
)

# compute season averages -----------------------------------------------

season_ave <- summary_models %>% 
  group_by(Lat,Lon,NUTS3,NUTS2,Year) %>%
  summarise(Tx_mean = mean(Tmax),
            Tx_max = max(Tmax),
            Tn_mean = mean(Tmin),
            Tn_min = min(Tmin),
            Lw_tot = sum(LW),
            RHx_mean = mean(RHmax),
            RHn_mean = mean(RHmin),
            DP_mean = mean(DewPoint),
            VPD_mean = mean(VPD),
            R310 = max(Sum310),
            Epi = max(SumEpi),
            Ipi = max(SumIpi),
            Dmc = max(SumDmc),
            Mag = max(SumMag),
            Ros = max(SumRos),
            Mis = max(SumMis),
            Lao = max(SumLao),.groups = 'drop') %>% 
  rowwise() %>% 
  mutate(Ave = mean(c(R310,Epi,Ipi,Dmc,Mag,Ros,Mis,Lao)),
         Std = sd(c(R310,Epi,Ipi,Dmc,Mag,Ros,Mis,Lao)),
         Med = median(c(R310,Epi,Ipi,Dmc,Mag,Ros,Mis,Lao)),
         NUTS3 = as.factor(NUTS3),
         NUTS2 = as.factor(NUTS2), Year = as.factor(Year)) %>% 
  ungroup()
#count rainy days
rain_summary <- summary_models %>% 
  group_by(Lat,Lon,NUTS3,NUTS2,Year) %>% 
  filter(Prec > 0.2) %>% 
  summarise(Prec_tot = sum(Prec),
            Rd_count = n(),
            R_ind = Prec_tot/Rd_count,
            .groups = 'drop') %>% 
  ungroup()
#append
season_ave <- add_column(season_ave, rain_summary[,c("Prec_tot","Rd_count","R_ind")], .after = "Year")


# Compute NUTS2 averages --------------------------------------------------

NUTS2_ave <- season_ave %>% group_by(NUTS2, Year) %>% 
  summarise(Tx_mean = mean(Tx_mean),
            Tx_max = mean(Tx_max),
            Tn_mean = mean(Tn_mean),
            Tn_min = mean(Tn_min),
            Lw_tot = mean(Lw_tot),
            RHx_mean = mean(RHx_mean),
            RHn_mean = mean(RHn_mean),
            DP_mean = mean(DP_mean),
            VPD_mean = mean(VPD_mean),
            Prec_tot = mean(Prec_tot),
            Rd_count = mean(Rd_count),
            R_ind = mean(R_ind),
            R310 = mean(R310),
            Epi = mean(Epi),
            Ipi = mean(Ipi),
            Dmc = mean(Dmc),
            Mag = mean(Mag),
            Ros = mean(Ros),
            Mis = mean(Mis),
            Lao = mean(Lao),
            Ave = mean(Ave),
            Std = mean(Std),
            Med = mean(Med),
            .groups = 'drop')
  

# PCA ---------------------------------------------------------------------
colnames(NUTS2_ave)
activevar <- NUTS2_ave %>%  select(c("Tx_mean":"R_ind","Ave","Med"))

#check for the presence of NA
NAval <- activevar %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_each(funs(sum(is.na(.)))) 

library(factoextra)
library(FactoMineR)
res.pca <- PCA(activevar, ncp = 3, graph = F, scale.unit = T)

# Supplementary variables -------------------------------------------------

res.pca2 <- PCA(NUTS2_ave %>% select(c("NUTS2":"R_ind","Ave","Med")), ncp = 3, graph = T, scale.unit = T, quali.sup = c("NUTS2","Year"))

# access pca results ------------------------------------------------------
# Eigenvalues
eig.val <- get_eigenvalue(res.pca2)
# Results for Variables
res.var <- get_pca_var(res.pca2)
res.var$coord          # Coordinates
res.var$cor            # Correlation between variables and pC
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca2)
head(res.ind$coord, 100)          # Coordinates
head(res.ind$contrib, 100)        # Contributions to the PCs
head(res.ind$cos2,100)           # Quality of representation


# PCA visualizations ------------------------------------------------------
#Correlation between variables
require(FactoMineR)
require(corrplot)
corrplot(cor(activevar), is.corr=T)
#scree plot
fviz_eig(res.pca2)
#quality of representation of the variables 
require(corrplot)
corrplot(res.var$cos2, is.corr=F)
# Total cos2 of variables on Dim.1 and Dim.2
fviz_cos2(res.pca2, choice = "var", axes = 1)
fviz_cos2(res.pca2, choice = "var", axes = 2)
#graph of variables
fviz_pca_var(res.pca2,max.overlaps = 5,
             col.var = "cos2", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
#contribution of variables to PC
corrplot(res.var$contrib, is.corr=FALSE)   
fviz_contrib(res.pca2, choice = "var", axes = 3)
#Quality and contribution of individuals and variables: biplot
fviz_pca_biplot(res.pca2,  col.ind = "cos2", geom = "point",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

# clustering of the features obtained by PCA ------------------------------

res.hcpc <- HCPC(res.pca2, nb.clust = -1, min =4, max = 10, graph = T, method = "ward", consol = T, iter.max = 1000, metric = "manhattan", graph.scale	
= "sqrt-inertia", description = T)

#dendrogram
fviz_dend(res.hcpc, 
          show_labels = F,                
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
)
#
fviz_cluster(res.hcpc,
             repel = F,
             show.clust.cent = F, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map",
)

# plot for printing -------------------------------------------------------
quali.var <- res.pca2$quali.sup$coord %>% as.data.frame() %>%  rownames_to_column ("Var")  

options(ggrepel.max.overlaps = 20)
fviz_pca_biplot(res.pca2, 
                # Individuals
                geom.ind = "point",
                repel = TRUE,
                col.ind = res.hcpc$data.clust$clust,
                pointshape = 19,
                pointsize = 1.5,
                palette = "lancet",
                addEllipses = TRUE,
                select.var = list(contrib = 6),
                ellipse.type ="convex",
                # Variables
                col.var = "gray48",
                title = NULL,
                max.overlaps = 20,
                legend.title = list(fill = "Clusters", color = "Clusters")
)+
  ggrepel::geom_label_repel(data =quali.var[1:20,], aes(x=Dim.1,y=Dim.2, label = Var), max.overlaps = 30, color = "black")

toplot <- data.frame(
x = res.pca2$ind$coord[,1],
y = res.pca2$ind$coord[,2],
z = res.pca2$ind$coord[,3],
col = as.factor(res.hcpc$data.clust$clust))

require(plotly)
plot_ly(toplot, x = ~x, y = ~y, z = ~z, color = ~col, alpha=0.8) %>%
  add_markers(size = 1) %>% 
  layout(scene = list(xaxis = list(title = 'PC 1'),
                      yaxis = list(title = 'PC 2'),
                      zaxis = list(title = 'PC 3')))
#make custom color palette
nuts <- NUTS2_ave$NUTS2
year <- NUTS2_ave$Year
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col=sample(color, nlevels(nuts))
require(plotly)
ggplotly(fviz_pca_biplot(res.pca2,  fill.ind = nuts,label = "none",map = "rowgab",addEllipses = TRUE, ellipse.type = "convex",geom.ind = "point",pointshape = 21, pointsize = 1, gradient.cols =  "RdYlBu"))


# Other clustering methods ------------------------------------------------

#res.pca.l <-prcomp(activevar, center=T, scale.=T,rank. = 3)#keep only the first two components


results <- as.data.frame(res.pca2$ind$coord)

require(fastcluster)
require(parallelDist)
#res_dist <- parDist(scale(as.matrix(results)),method = "euclidean") ###scale(as.matrix(results)),method = "euclidean" is the same as computing mahalanobis distances on PCA-transformed data:
#see: https://stats.stackexchange.com/questions/166525/is-mahalanobis-distance-equivalent-to-the-euclidean-one-on-the-pca-rotated-data
res_dist <- parDist(as.matrix(results),method = "manhattan")
res.clust <- fastcluster::hclust(res_dist, method = "ward.D2")
#number of clusters (using manhattan distances)
nb1 <- fviz_nbclust(results, cluster::clara, method = "wss", diss = res_dist) +
  geom_vline(xintercept = 4, linetype = 2)+labs(title = NULL) +theme_classic()+
  annotate("text", x = 4.2, y = 10000, label = "Optimal N° of clusters", angle = 90)

nb2 <- fviz_nbclust(results, cluster::clara, method = "silhouette", diss = res_dist) +  geom_vline(xintercept = 4, linetype = 2)+labs(title = NULL)+theme_classic() + theme(plot.margin = unit(c(5.5,5.5,5.5,5.5), "pt"))+
  annotate("text", x = 4.2, y = 0.11, label = "Optimal N° of clusters", angle = 90)

nb3 <- fviz_nbclust(results, cluster::clara, method = "gap_stat", diss = res_dist) +  geom_vline(xintercept = 4, linetype = 2)+labs(title = NULL) +theme_classic()+theme(plot.margin = unit(c(5.5,5.5,5.5,14), "pt"))+
  annotate("text", x = 4.2, y = 0.55, label = "Optimal N° of clusters", angle = 90)
#dendrograms
plot(res.clust)
fviz_dend(res.clust, k = 3,
          palette = "lancet",            
          rect = TRUE, rect_fill = T,
          rect_border = "lancet", 
          labels_track_height = 0,
          lower_rect = 0,
          show_labels = T,
          main = NULL,
          lwd = 0.3,
          rect_lty = 2,
          ggtheme = theme_bw())
#check silhouette of hierarchical clustering
res.clust2 <- hcut(res_dist, method = "ward.D2",  k = 4)
fviz_silhouette(res.clust2) #not that great

#dbscan
require(dbscan)
kNNdistplot(results, k = 10)
abline(h=1.3, col = "red", lty = 2)
res.db <- dbscan(results, eps = 1, minPts = 5)

fviz_pca_biplot(res.pca2, fill.ind = as.factor(res.db$cluster), label = "none",map = "rowgab",addEllipses = TRUE, ellipse.type = "convex",geom.ind = "point",col.ind = "black",palette = "jco",alpha.var ="contrib", col.var = "contrib",col.quali.sup = "red",
                gradient.cols = "RdYlBu",pointshape = 21, pointsize = 1)

#clara
require(cluster)
clara.res <- clara(results, 4,sampsize = 200,samples = 200, pamLike = T, metric = "manhattan")
fviz_pca_biplot(res.pca2, 
                # Individuals
                geom.ind = "point",
                repel = TRUE,
                col.ind = as.factor(clara.res$clustering),
                pointshape = 19,
                pointsize = 1.5,
                palette = "lancet",
                addEllipses = TRUE,
                select.var = list(contrib = 6),
                ellipse.type ="convex",
                # Variables
                col.var = "gray48",
                title = NULL,
                max.overlaps = 20,
                legend.title = list(fill = "Clusters", color = "Clusters"))+
  ggrepel::geom_label_repel(data = quali.var[1:20,], aes(x=Dim.1,y=Dim.2, label = Var), max.overlaps = 30, color = "black")
#check silhouette of clara clustering
fviz_silhouette(clara.res) #better

#hybrid hierarchical and kmeans clustering (the same as HCPC)
res.clust2
hcn <- res.clust2$cluster
# Compute cluster centers
clus.centers <- aggregate(results, list(hcn), median)
# Remove the first column
clus.centers <- clus.centers[, -1]
hck.res2 <- eclust(results, "kmeans", k = clus.centers, graph = FALSE, seed = 123,  iter.max = 1000)
#silhouette
scales::show_col(pal_lancet("lanonc")(9))
colors <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF")
names(colors) <- c(1,2,3,4)

fviz_silhouette(hck.res2)+labs(title = NULL, subtitle = NULL)+scale_color_manual(values = colors)+scale_fill_manual(values = colors)+labs(fill = "Cluster", color = "Cluster")

#agnes
res.agnes <- eclust(results,FUNcluster = "agnes", k = 4, graph = F, hc_metric = "manhattan", hc_method = "ward.D2")
fviz_silhouette(res.agnes)

#fanny
res.fanny <- eclust(results,FUNcluster = "fanny", k = 4, graph = F, hc_metric = "manhattan", hc_method = "ward.D2", maxit = 500)
fviz_silhouette(res.fanny)


# plot clustering results -------------------------------------------------

#plot PCA results with clustering of hybrid approach 
fviz_pca_biplot(res.pca2, 
                # Individuals
                geom.ind = "point",
                repel = TRUE,
                col.ind = as.factor(hck.res2$cluster),
                pointshape = 19, pointsize = 1.5,
                palette = "lancet",
                addEllipses = TRUE,
                select.var = list(contrib = 6),
                ellipse.type ="convex",
                # Variables
                col.var = "gray48",
                max.overlaps = 20,
                legend.title = list(fill = "Clusters", color = "Clusters")
)+
  ggrepel::geom_label_repel(data = quali.var[1:20,], aes(x=Dim.1,y=Dim.2, label = Var), max.overlaps = 30, color = "black")

#plot with plotly
require(rgl)
toplot <- data.frame(
  x = results$Dim.1,
  y = results$Dim.2,
  z = results$Dim.3,
  col = as.factor(clusters))

require(plotly)
plot_ly(toplot, x = ~x, y = ~y, z = ~z, color = ~col, alpha=0.8) %>%
  add_markers(size = 1) %>% 
  layout(scene = list(xaxis = list(title = 'PC 1'),
                      yaxis = list(title = 'PC 2'),
                      zaxis = list(title = 'PC 3')))
# clustering tendency -----------------------------------------------------

#visual data
require(hopkins)
hopkins(results)
#visual data
x11()
fviz_dist(res_dist, show_labels = FALSE, order = T,gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

# clustering on original data ---------------------------------------------

#scale data 
require(caret)
preProcValues <-preProcess(activevar, method = c("center", "scale", "YeoJohnson"))
actvar_proc <- predict(preProcValues,activevar)
res_dist <- parDist(as.matrix(scale(actvar_proc)),method = "euclidean")
res.clust <- fastcluster::hclust(res_dist, method = "ward.D2")
plot(res.clust)
clusters = cutree(res.clust, k = 3)
plot(actvar_proc, pch=20, col=rainbow(4, alpha=c(0.2,0.2,0.2,1))[clusters])

# clustering tendency -----------------------------------------------------

#visual data
require(hopkins)
hopkins(activevar)
#visual data
x11()
fviz_dist(res_dist, show_labels = FALSE)


# Uniform manifold approximation projection -------------------------------

library(umap)
res_umap <- umap(actvar_proc)
plot(res_umap$data)


# t-sne -------------------------------------------------------------------

library(tsne)
res_tsne <- tsne(as.matrix(actvar_proc), initial_dims = 4, max_iter = 300, perplexity = 50, whiten = T)
plot(res_tsne)



# Factor Analysis of Mixed Data -------------------------------------------

library(FactoMineR)
res.famd <- FAMD(NUTS2_ave[,c(1:14,23,25)], ncp = 10, graph = FALSE)
fviz_screeplot(res.famd)
fviz_famd_var(res.famd, repel = TRUE)
fviz_contrib(res.famd, "var", axes = 1)
fviz_contrib(res.famd, "var", axes = 2)
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
fviz_famd_ind(res.famd, col.ind = "cos2", label = F,
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)


res.hcpc2 <- HCPC(res.famd, nb.clust = -1, min =4, max = 10, graph = TRUE, method = "ward", consol = T, iter.max = 500, metric = "manhattan", graph.scale	
                 = "sqrt-inertia", description = T)


# compute distance metrics by climate type (cluster) ----------------------
NUTS2_ave$ClimType <- res.hcpc$data.clust$clust
NUTS2_ave$ClimType2 <- hck.res2$cluster
compare <- NUTS2_ave[,c("NUTS2","ClimType","ClimType2")]

cluster_regions <- compare %>% group_by(NUTS2) %>% summarise(MedClim1 = as.factor(median(as.numeric(ClimType))),
                                                              MedClim2 = as.factor(median(as.numeric(ClimType2))))
levels(cluster_regions$MedClim1)<- c("Mild-Dry","Mild-Wet","Hot-Dry","Cold-Wet")
levels(cluster_regions$MedClim2)<- c("Mild-Dry","Hot-Dry","Mild-Wet","Cold-Wet")
#for printing
clustersp <- as.factor(hck.res2$cluster)
levels(clustersp) <- c("Mild-Wet","Mild-Dry","Hot-Dry","Cold-Wet")

fviz_pca_biplot(res.pca2, 
                # Individuals
                geom.ind = "point",
                repel = TRUE,
                col.ind = clustersp,
                #pointshape = 19, 
                pointsize = 1.5,
                palette = "lancet",
                addEllipses = TRUE,
                select.var = list(contrib = 6),
                ellipse.type ="convex",
                title = NULL,
                labelsize = 4,
                mean.point = FALSE,
                # Variables
                col.var = "gray48",
                set.seed = 345,
                legend.title = list(color = "Clusters", fill = "Clusters", shape = "Clusters")
)+
ggrepel::geom_label_repel(data = quali.var[1:20,], aes(x=Dim.1,y=Dim.2, label = Var, color = as.factor(cluster_regions$MedClim2)), max.overlaps = 30, segment.colour = "black", show.legend = F)+
  xlab("Principal Component 1 (62.1%)")+
  ylab("Principal Component 2 (20.9%)")+
  theme_bw()
  
#save
ggsave(
  "Biplot.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

#silhouette
colors <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF")
names(colors) <- c(1,2,3,4)
nb4 <- fviz_silhouette(hck.res2, alpha = 0.1, print.summary = T)+labs(title = NULL, subtitle = NULL)+scale_color_manual(values = colors, labels = unique(clustersp))+scale_fill_manual(values = colors, labels = unique(clustersp))+labs(fill = "Clusters", color = "Clusters")+theme_classic()+theme(axis.text.x = element_text(color = "white"), axis.ticks.length.x = unit(0, "cm"), plot.margin = unit(c(5.5,5.5,5.5,2), "pt"), legend.position = c(0.5, 0.9)) +labs(x = "")+guides(colour = guide_legend(nrow = 1))+annotate("text", x = 100, y = 0.41, label = "Average Silhouette")
#grid 
require(ggpubr)
ggarrange(nb1,nb2,nb3,nb4,ncol =2,nrow = 2)
ggsave(
  "ClusterVal.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

data_season <- left_join(data_season,cluster_regions, by = "NUTS2") %>% 
  rename(Laore = laore) %>% 
  rename(Misfits = misfits)

#multicore
faithdist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "faith"))

jaccardist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "binary"))

bbdist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "braun-blanquet"))

dicedist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "dice"))

fagerdist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "fager"))

hammandist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "hamman"))

kuldist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "kulczynski2"))

micheldist <-data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "michael"))

mountdist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "mountford"))

mozleydist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "mozley"))

ochiaidist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "ochiai"))

phidist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "phi"))

russeldist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "russel"))

tanimotodist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "tanimoto"))

yuledist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "yule2"))

stilesdist <- data_season %>% dplyr::select(MedClim2,Rule310,Epi,Ipi,Dmcast,Magarey,Rossi,Misfits,Laore) %>% group_by(MedClim2) %>% group_map(~ parallelDist::parDist(as.matrix(t(.x)),method = "stiles"))

levels(data_season$MedClim2)
# Hierarchical clustering of binary distance estimates --------------------
#jaccard (co-occurrence coefficient that excludes negative matches)
hclust.jac <- lapply(jaccardist, hclust, method = "ward.D2")
p1 <- lapply(hclust.jac, function(x) fviz_dend(x, k = 2,
                                               palette = "jco",            
                                               rect = TRUE, rect_fill = F,
                                               rect_border = "jco", 
                                               labels_track_height = -0,
                                               lower_rect = -0.4,
                                               show_labels = T,
                                               main = NULL,
                                               cex = 0.6,
                                               lwd = 0.8,
                                               rect_lty = 3, set.seed = 134,
                                               ggtheme = theme_classic()))
require(ggpubr)
require(ggsci)
arr.jac <- ggarrange(p1[[1]],p1[[2]],p1[[3]],p1[[4]],ncol =2,nrow = 2,
                     labels = c("Mild-Dry","Hot-Dry","Mild-Wet","Cold-Wet"), font.label = c(size = 10), label.x = 0.75, label.y = 1.02)
annotate_figure(arr.jac, text_grob("Hierarchical Clustering on Jaccard Distances"))
#save
ggsave(
  "Jacdist.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)
#phi (correlation coefficient)
hclust.phi <- lapply(phidist, hclust, method = "ward.D2")
p2 <- lapply(hclust.phi, function(x) fviz_dend(x, k = 3,
                                              palette = "aaas",            
                                              rect = TRUE, rect_fill = F,
                                              rect_border = "aaas", 
                                              labels_track_height = -1,
                                              lower_rect = -0.4,
                                              show_labels = T,
                                              main = NULL,
                                              cex = 0.8,
                                              lwd = 0.8,horiz = T,
                                              rect_lty = 3, set.seed = 1345,
                                              ggtheme = theme_classic()))
  
require(ggpubr)
require(ggsci)
arr.phi <- ggarrange(p2[[2]],p2[[4]],p2[[1]],p2[[3]],ncol =2,nrow = 2,
                     labels = c("Hot-Dry","Cold-Wet","Mild-Dry","Mild-Wet"), font.label = c(size = 10), label.x = 0.75, label.y = 1.02)
annotate_figure(arr.phi, top = text_grob("Hierarchical Clustering on φ Distances",
                                         hjust = 0, x = 0.05, vjust = 0.1, face = "bold", size = 12))
#save
ggsave(
  "Phidist.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)
# Michael (co-occurrence coefficient that includes negative matches)
hclust.mic <- lapply(micheldist, hclust, method = "ward.D2")
p3 <- lapply(hclust.mic, function(x) fviz_dend(x, k = 2,
                                               palette = "jco",            
                                               rect = TRUE, rect_fill = F,
                                               rect_border = "jco", 
                                               labels_track_height = 0,
                                               lower_rect = -0.4,
                                               show_labels = T,
                                               main = NULL,
                                               cex = 0.6,
                                               lwd = 0.8,
                                               rect_lty = 3, set.seed = 134,
                                               ggtheme = theme_classic()))
require(ggpubr)
require(ggsci)
arr.mic <- ggarrange(p3[[1]],p3[[2]],p3[[3]],p3[[4]],ncol =2,nrow = 2,
                     labels = c("Mild-Dry","Hot-Dry","Mild-Wet","Cold-Wet"), font.label = c(size = 10), label.x = 0.75, label.y = 1.02)
annotate_figure(arr.mic, text_grob("Hierarchical Clustering on Michael Distances"))
#save
ggsave(
  "Micdist.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

# boxplot -----------------------------------------------------------------

season_ave <- left_join(season_ave,cluster_regions, by = "NUTS2")
colnames(season_ave)
seasave_piv <- season_ave %>% pivot_longer(cols = c("R310","Epi","Ipi","Dmc","Mag","Ros","Mis","Lao"),names_to = "Model",values_to = "YearSum")

seasave_piv$Model <- factor(seasave_piv$Model,     # Reorder factor levels
                            rev(c("Mag", "Mis", "R310","Epi","Ipi","Dmc","Ros","Lao")))

#cut data for Mis and Mag Models for best representation
seasave_MM <- seasave_piv %>% filter(Model %in% c("Mag", "Mis"))
seasave_MM <- seasave_MM %>% filter(YearSum < 30)
seasave_piv <- seasave_piv %>% filter(!Model %in% c("Mag", "Mis"))
seasave_piv <- rbind(seasave_piv,seasave_MM)
#cut data for DMC and Epi and Lao Models for best representation
seasave_DEL <- seasave_piv %>% filter(Model %in% c("Dmc", "Epi", "Lao"))
seasave_DEL <- seasave_DEL %>% filter(YearSum < 105)
seasave_piv <- seasave_piv %>% filter(!Model %in% c("Dmc", "Epi", "Lao"))
seasave_piv <- rbind(seasave_piv,seasave_DEL)
#cut data for Ipi and Ros and R310 Models for best representation
seasave_IRR <- seasave_piv %>% filter(Model %in% c("Ipi", "Ros", "R310"))
seasave_IRR <- seasave_IRR %>% filter(YearSum < 70)
seasave_piv <- seasave_piv %>% filter(!Model %in% c("Ipi", "Ros", "R310"))
seasave_piv <- rbind(seasave_piv,seasave_IRR)
#compute mean by group
OVE <- seasave_piv %>% group_by(Model) %>% summarise(MN = median(YearSum), Q25 = quantile(YearSum, 0.25), Q75 = quantile(YearSum, 0.75)) %>% arrange(-MN)


seasave_piv %>% ggplot()+
  geom_hline(data = OVE, aes(yintercept = MN), color = "darkgray", linetype = "solid")+
  geom_hline(data = OVE, aes(yintercept = Q25), color = "darkgray", linetype = "dashed")+
  geom_hline(data = OVE, aes(yintercept = Q75), color = "darkgray", linetype = "dashed")+
  geom_rect(data = OVE,aes(xmin = -Inf, xmax = Inf, ymin = Q25, ymax = Q75), alpha = .3, fill = "darkgray")+
  geom_violin(aes(y = YearSum, x = reorder(MedClim2,YearSum,FUN = median), color = MedClim2, fill = MedClim2),position = "dodge",scale = "width", width = 0.8, trim = T)+
  geom_boxplot(aes(y=YearSum, x= reorder(MedClim2,YearSum,FUN = median),alpha = 1), width = 0.2, show.legend = F,outlier.shape = NA)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  #scale_y_continuous(limits = c(0, 100)) +
  theme_minimal()+
  labs(color = "Model", fill = "Model")+
  ylab("Average season (Mar-Sept) infection days")+
  facet_wrap(vars(factor(Model, levels=c(OVE$Model))), scales = "free")+
  # nrow = 2, ncol = 2)+
  theme(legend.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
    axis.title.x = element_blank(),
    plot.background  = element_rect(color = "white"),
    strip.text.x = element_text(size = 10, face = "bold"))+ 
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1))
#save
ggsave(
  "BoxPlotModels.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

# Cumulated infections plot -----------------------------------------------
#prepare data

write_csv(left_join(summary_models,cluster_regions, by = "NUTS2")[,c(1:4,39:40,43:51,53)], "summary_models.csv")
test1 <- big_data  %>% 
  mutate(Doy = as.factor(data.table::yday(Date)), .after = Date)%>% filter(Doy == 188)
test2 <- big_data  %>% 
  mutate(Doy = as.factor(data.table::yday(Date)), .after = Date)%>% filter(Doy == 189)
levels(as.factor(test2$Year))
dim(test2)
dim(test1)
levels(as.factor(test1$Year))
diff <- anti_join(test2, test1)

ToPlotCum <- left_join(big_data,cluster_regions, by = "NUTS2") %>% 
  mutate(Doy = as.factor(data.table::yday(Date)), .after = Date) %>% 
  group_by(Doy,MedClim2) %>% 
  summarise(
  R310 = mean(Sum310),
  Sd310 = sd(Sum310)*0.25,
  Epi = mean(SumEpi),
  SdEpi = sd(SumEpi)*0.25,
  Ipi = mean(SumIpi),
  SdIpi = sd(SumIpi)*0.25,
  Dmc = mean(SumDmc),
  SdDmc = sd(SumDmc)*0.25,
  Mag = mean(SumMag),
  SdMag = sd(SumMag)*0.25,
  Ros = mean(SumRos),
  SdRos = sd(SumRos)*0.25,
  Mis = mean(SumMis),
  SdMis = sd(SumMis)*0.25,
  Lao = mean(SumLao),
  SdLao = sd(SumLao)*0.25,
)
ToPlotCumA <- ToPlotCum %>% 
  pivot_longer(cols = c("R310","Epi","Ipi","Dmc","Mag","Ros","Mis","Lao"),names_to = "Model",values_to = "YearCumSum") 
ToPlotCumB <- ToPlotCum %>% 
  pivot_longer(cols = c("Sd310","SdEpi","SdIpi","SdDmc","SdMag","SdRos","SdMis","SdLao"),names_to = "Model",values_to = "YearSd") 

ToPlotCum <- cbind(ToPlotCumA[,c(1:2,11:12)], YearSd = ToPlotCumB$YearSd)

ggplotly(ToPlotCum %>% ggplot()+geom_line(aes(x=as.numeric(Doy), y=YearCumSum, color = Model), size = 0.2)+
  scale_alpha(range = c(0, 0.4))+
  geom_ribbon(aes(x=as.numeric(Doy), ymin = YearCumSum-YearSd, ymax = YearCumSum+YearSd, fill = Model, alpha = 1),show.legend = F)+
  facet_wrap(vars(MedClim2), scales = "free")+
  theme_light())
# analyse yearly weather data ---------------------------------------------

year_ave <- season_ave %>% 
  group_by(Year) %>% 
  summarise(Tx_mean = mean(Tx_mean),
            Tx_max = mean(Tx_max),
            Tn_mean = mean(Tn_mean),
            Tn_min = mean(Tn_min),
            Tm_mean = mean(Tx_mean+Tn_mean)/2,
            Lw_tot = mean(Lw_tot),
            RHx_mean = mean(RHx_mean),
            RHn_mean = mean(RHn_mean),
            DP_mean = mean(DP_mean),
            VPD_mean = mean(VPD_mean),
            Prec_tot = mean(Prec_tot),
            Rd_count = mean(Rd_count),
            R_ind = mean(R_ind),
            .groups = 'drop') %>% 
  ungroup() %>% 
  mutate(T_mean = mean(Tm_mean),
         T_min = min(Tn_mean),
         T_max = max(Tn_mean),
         Prec_mean = mean(Prec_tot),
         Prec_max = max(Prec_tot),
         Prec_min = min(Prec_tot),
         Tdev = Tm_mean- T_mean ,
         Pdev = Prec_tot- Prec_mean)

require(RColorBrewer)
require(colorspace)
require(ggtext)
mycolors <- colorRampPalette(brewer.pal(9, "Spectral"), interpolate = "spline")(30)
mycolors <- darken(mycolors, 0.2)
scale_size <- seq(4,10, length.out = 30)
labels <- year_ave %>% filter(Year %in% c(1996,2012,2003,2002))
names(scale_size) <- as.factor(c(1991:2020))

ggplot(data = year_ave,aes(y=Tdev, x=Pdev, fill = Year, size = Year))+
  geom_hline(yintercept=0, linetype = "dashed")+
  geom_vline(xintercept=0, linetype = "dashed")+
  ggrepel::geom_text_repel(data = labels, aes(label = Year, fill = NULL, size = NULL), 
                           show.legend = F, box.padding = 1.5, seed = 101)+
  geom_point(pch = 21, color = "white")+
  scale_fill_manual(values = rev(mycolors)) +
  xlim(c(-180,180))+ ylim(c(-2,2))+
  scale_size_manual(values = scale_size)+
  theme_bw()+
  ylab("Season (Mar-Sept) T°C - Average T°C (1991-2020)")+
  xlab("Season (Mar-Sept) Prec (mm) - Average Prec (mm) (1991-2020)")+
  geom_richtext(aes(label = "Warm - Wet", x = 100, y = 1.8), size = 5, fill = NA,label.color = NA)+
  geom_richtext(aes(label = "Cold - Wet", x = 100, y = -1.5), size = 5, fill = NA,label.color = NA)+
  geom_richtext(aes(label = "Cold - Dry", x = -100, y = -1.5), size = 5, fill = NA,label.color = NA)+
  geom_richtext(aes(label = "Warm - Dry", x = -100, y = 1.8), size = 5, fill = NA,label.color = NA)+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.5,"line"))+ 
  guides(fill = guide_legend(nrow = 3))

#save
ggsave(
  "ClimatePlot.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 1.5,
  width = 16,
  height = 12,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)
#
year_ave %>% filter(Year == 2012) %>% select(Tm_mean, Tdev, Prec_tot, Pdev)

# visualize on map --------------------------------------------------------

library(giscoR)
library(sf)
library(ggsn)
library(ggspatial)

#get map of italy wth NUT levels (sf object)
gisco<-gisco_get_nuts(
  year = "2021",
  epsg = "4326",
  cache = TRUE,
  update_cache = FALSE,
  cache_dir = NULL,
  verbose = FALSE,
  resolution = "1", #resolution
  spatialtype = "RG",
  country = "IT", #Country
  nuts_id = NULL,
  nuts_level = "3") #%>%  #NUTS level  


#grid of points
grid <- gisco %>% 
  st_bbox() %>% 
  st_make_grid(square=T,cellsize = 0.3, crs = 4326) %>%
  st_sf() 

#prepare data 
modOut <-  season_ave %>%
  mutate(grid=paste0(Lat,'_',Lon)) %>% 
  group_by(Year, grid, Lat, Lon) %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) 

outMaps <- grid %>% st_join(modOut, join = st_contains) %>% 
  filter(!is.na(Year))

#filter the selected years
FoutMaps <- outMaps %>% filter(Year %in% c(1996,2012,2003,2002))
#map



ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Ave),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Average of the Models)",
       fill = "Seasonal\n Infections\n (n)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"),
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


ggsave(
  "MapAverageModel.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)


ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Std),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1),
                       breaks = c(0,5,10,15,20,25))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Deviation Among the Models)",
       fill = "Seasonal\n Infections\n (n)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggsave(
  "MapStdModel.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)


# maps of weather parameters ----------------------------------------------

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Prec_tot),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=0.1), 
                       breaks = c(50,100,200, 400, 700, 1300,  2200),
                       labels = c(50,100,200, 400, 700, 1300,  2200))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = NULL,
       fill = "Prec \n(mm)")+
  annotation_scale(location = "tr", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold",hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))
ggsave(
  "MapAvePrec.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Tn_mean),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", 
                       trans = scales::modulus_trans(p=1.5),
                       breaks = c(0,6,10, 14, 17, 20),
                       labels = c(0,6,10, 14, 17, 20))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = NULL,
       fill = "T°C")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggsave(
  "MapAveTemp.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)
# maps of separate models on selected years -------------------------------



ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Ros),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Rossi)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


ggsave(
  "MapRossiModel.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)


ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Mis),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Misfits)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggsave(
  "MapMisfitsModel.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Dmc),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Dimcast)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


ggsave(
  "MapDimcastModel.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)



ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=R310),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (R310)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Epi),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Epi)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Ipi),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Ipi)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Mag),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Magarey)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


ggplot(data = FoutMaps)+
  geom_raster(aes(x=Lon, y=Lat, fill=Lao),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1))+
  geom_sf(data=gisco, alpha =0, size = 0.2, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Year)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Laore)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 2,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

# Map of model comparison for year 2020 (very wet) ------------------------

outMaps_2002 <- outMaps %>% filter(Year %in% c(2002))
outMaps_2002 <- outMaps_2002 %>% pivot_longer(cols = c("R310","Epi","Ipi","Dmc","Mag","Ros","Mis","Lao"),names_to = "Model",values_to = "YearSum")
orderMaps2002 <- outMaps_2002 %>% group_by(Model) %>% summarise(Mean = mean(YearSum)) %>% arrange(-Mean)

ggplot(data = outMaps_2002)+
  geom_raster(aes(x=Lon, y=Lat, fill=YearSum),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1), limits= c(0,60), oob=scales::squish)+
  geom_sf(data=gisco, alpha =0, size = 0.3, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(vars(factor(Model, levels=c(orderMaps2002$Model))), nrow = 2)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Model Comparison - 2002)",
       fill = "Seasonal\n Infections\n (n)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))

ggsave(
  "MapModelComparison.png",
  device = "png",
  path = "C:/Users/loren/Documents/Tesi_master_AP",
  scale = 2.5,
  width = 16,
  height = 10,
  units = c("cm"),
  dpi = 600,
  limitsize = F
)

outMaps_2003 <- outMaps %>% filter(Year %in% c(2003))
outMaps_2003 <- outMaps_2003 %>% pivot_longer(cols = c("R310","Epi","Ipi","Dmc","Mag","Ros","Mis","Lao"),names_to = "Model",values_to = "YearSum")

ggplot(data = outMaps_2003)+
  geom_raster(aes(x=Lon, y=Lat, fill=YearSum),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1), limits= c(0,60), oob=scales::squish)+
  geom_sf(data=gisco, alpha =0, size = 0.3, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Model, nrow = 2)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Model Comparison - 2003)",
       fill = "Seasonal\n Infections\n (n)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


outMaps_1999 <- outMaps %>% filter(Year %in% c(1999))
outMaps_1999 <- outMaps_1999 %>% pivot_longer(cols = c("R310","Epi","Ipi","Dmc","Mag","Ros","Mis","Lao"),names_to = "Model",values_to = "YearSum")

ggplot(data = outMaps_1999)+
  geom_raster(aes(x=Lon, y=Lat, fill=YearSum),interpolate = F)+
  #scale_fill_distiller(palette = "Spectral", trans = scales::boxcox_trans(p=1.1))+
  scale_fill_viridis_c(option = "turbo", trans = scales::boxcox_trans(p=1), limits= c(0,60), oob=scales::squish)+
  geom_sf(data=gisco, alpha =0, size = 0.3, color = "black")+
  #geom_sf(data=grid,alpha=0.2)+
  #geom_sf_text(data = gisco,aes(label=NUTS_NAME))+
  theme_minimal()+
  facet_wrap(~Model, nrow = 2)+
  #xlab("Longitude (°)")+
  #ylab("Latitude (°)")+
  labs(title = "Primary Infection Events (Model Comparison - 1999)",
       fill = "Seasonal\n Infections\n (n)")+
  annotation_scale(location = "bl", style = "bar", width_hint = 0.3, pad_y = unit(0.76,"cm"),pad_x = unit(0.2, "cm")) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_y = unit(1.8,"cm"),pad_x = unit(1.5, "cm"), 
                         style = north_arrow_orienteering, width = unit(1.5,"cm"), height =  unit(1.5,"cm"))+
  theme(legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title.align = 0.5,
        title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = 0))


# Fit explorative model ---------------------------------------------------

#compute mean values among models (binary)
data_season <- data_season %>% mutate(MeanRisk = ifelse(InfRisk > 0.5, 1, 0))
#linear model
m <-glm(MeanRisk ~ Tmean, data_season, family = "binomial")
summary(m)

m1 <- glm(MeanRisk ~ Prec, data_season, family = "binomial")
summary(m1)

m2 <- glm(MeanRisk ~ LW, data_season, family = "binomial")
summary(m2)

m3 <- glm(MeanRisk ~ Tmin, data_season, family = "binomial")
summary(m3)

m4 <- glm(MeanRisk ~ RHmean, data_season, family = "binomial")
summary(m4)

m5 <-glm(MeanRisk ~ Prec + Tmin + LW , data_season, family = "binomial")
summary(m5)
#use the InfRisk var
m6 <-glm(InfRisk ~ Prec + Tmin + LW , data_season, family = quasibinomial('logit'))
summary(m6)

#betaregression
require(betareg)
data_season <- data_season %>% mutate(InfRisk =replace(InfRisk, InfRisk==1, 1-0.00000001)) %>% 
  mutate(InfRisk =replace(InfRisk, InfRisk==0, 0+0.00000001))
m7 <- betareg(InfRisk ~ Prec + Tmin + LW , data = data_season, link = "logit")

test <- (data_season$InfRisk*(length(data_season$InfRisk)-1)+0.5)/length(data_season$InfRisk)
test
range(test)
plot(test)
