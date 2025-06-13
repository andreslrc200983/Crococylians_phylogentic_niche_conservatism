######################################################################################################################################
######################################################################################################################################
##########################           PHYLOGENETIC NICHE CONSERVATISM IN AMERICAN CROCODYLIANS           ##############################
######################################################################################################################################
######################################################################################################################################

# The tables used in the analysis can be found in the Excel sheet within the GitHub repository:
# https://doi.org/10.15468/dd.vmdytt

######################################################################################################################################
# Install (if needed) and load packages required for the analysis 
######################################################################################################################################
library(sdm); library(scrubr) ; library(biomod2); library(rgbif); library(dismo); library(dplyr); library(mapview); library(plyr)
library(ecospat); library(ade4); library(terra); library(usdm); library(ape); library(phangorn); library(phytools); library(geiger)
library(ggtree); library(ggpmisc); library(ggplot2); library(vegan); library(reshape2); library(ggrepel); library(gridExtra); 
library(ENMTools); library(scales); library(grid)

######################################################################################################################################
# DOWNLOAD OCCURRENCE RECORDS FROM GBIF AND CLEANING
# For practical purposes, we showed the code for 1 species: Caiman yacare. This portion of the code was repeated for the other
# 10 American crocodylians. The complete uncleaned occurrences can be found in the file: "1_GBIF_original_records_downloaded.xlsx"
######################################################################################################################################
gbif("Caiman", "yacare", download = F) # Check the points available in GBIF
C_ya <- gbif("Caiman", "yacare", download = T) # Download the available points in GBIF
class(C_ya) # Check the class of the object. Should be "dataframe"
dim(C_ya) # Check the dimension of the the "C_ya" object
table(C_ya$basisOfRecord) # Check the type of data available in relation to the occurrences
C_ya_filter <- C_ya %>% dplyr::filter(basisOfRecord %in% c("HUMAN_OBSERVATION", "MATERIAL_CITATION", "PRESERVED_SPECIMEN"))
head(C_ya_filter)
C_ya_gbif <- C_ya_filter %>% dplyr::select(country, adm1, cloc,  references, year, eventDate, identifiedBy, recordedBy, locality,
                                           lat, lon, georeferenceRemarks, acceptedScientificName, basisOfRecord, 
                                           coordinateUncertaintyInMeters, individualCount, datasetKey) # Select columns of interest
head(C_ya_gbif) # Check the "C_ya_gbif" object, which should just contain all previous selected columns
nrow(C_ya_gbif) # Checking the number of records of C_ya_gbif with all data
unique(C_ya_gbif[c("year")]) # Check the unique values of years present in the data frame
C_ya_gbif_1 <- C_ya_gbif [complete.cases(C_ya_gbif[ , c(5,10)]),] # Delete "NA" records of "year" and "lat" (col 5,10)
C_ya_gbif_2 <- dplyr::filter(C_ya_gbif_1, year >= 1973) # Keep records with year >= 1973
C_ya_gbif_3 <- dplyr::filter(C_ya_gbif_2, coordinateUncertaintyInMeters <= 2000) # Keep records with uncertainty < 2000 m
unique(C_ya_gbif_3[c("country")]) # Check the unique values of "country"
C_ya_gbif_3 <- C_ya_gbif_3[!(C_ya_gbif_3$country=="Costa Rica"),] # Delete countries where species does not inhabit
C_ya_gbif_3 <- C_ya_gbif_3[!(C_ya_gbif_3$country=="Venezuela"),] # Delete countries where species does not inhabit
C_ya_gbif_3 <- C_ya_gbif_3[!(C_ya_gbif_3$country=="Ecuador"),] # Delete countries where species does not inhabit
C_ya_gbif_3 <- C_ya_gbif_3[!(C_ya_gbif_3$country=="Uruguay"),] # Delete countries where species does not inhabit
unique(C_ya_gbif_3[c("country")]) # Check the remaining values of "country" are correct 
C_ya_gbif_3 <- dframe(C_ya_gbif_3) %>% coord_impossible() # Delete occurrences which coordinates are impossible
C_ya_gbif_3 <- dframe(C_ya_gbif_3) %>% coord_unlikely() # Delete occurrences which coordinates are unlikely
nrow(C_ya_gbif_3) # Check the number of locations after previous filtering
duplicated(C_ya_gbif_3) # Identify duplicated records
C_ya_gbif_final <- unique(C_ya_gbif_3) # # Remove duplicate elements
nrow(C_ya_gbif_final) # Check the final number of locations 

# Visualize data in interactive map
head(C_ya_gbif_final)
C_ya_gbif_final_sp <- C_ya_gbif_final %>% dplyr::select(lon, lat) # Select "lat" and "lon" columns
head(C_ya_gbif_final_sp)
names(C_ya_gbif_final_sp)[1] <- "lon" # Assign the name "lon" to column lon (useful when original names are different)
names(C_ya_gbif_final_sp)[2] <- "lat" # Assign the name "lat" to column lat (useful when original names are different)
C_ya_gbif_final_sp$species <- 1 # Insert a column "Species" with values "1"
class(C_ya_gbif_final_sp) ## Check the class. Should be "dataframe"
coordinates(C_ya_gbif_final_sp) <- c("lon", "lat") # Convert dataframe into spatial points
class(C_ya_gbif_final_sp) # Check the class. Should be SpatialPointsDataFrame
proj4string(C_ya_gbif_final_sp) <- projection(raster()) # Add projection to points
mapview(C_ya_gbif_final_sp) # View the points in interactive map

# Export and save final dataframe as a .csv file
write.csv(C_ya_gbif_final, ".../C_ya_1_gbif.csv")

# Register the Derived Dataset for GBIF citation
C_ya_report <- table(C_ya_gbif_final$datasetKey) # Check the type of data available in relation to the occurrences
write.csv(C_ya_report, ".../C_ya_GBIF_credits.csv")

######################################################################################################################################
# After manual cleaning of GBIF records, we combined this dataset with records founded in published and gray literature. This dataset
# was then rarefied to obtain the final occurrences for each species (see details in main manuscript). Final occurrence datasets for
# each species can be found in the file: "2_Final_crocs_occurrences.xlsx".
######################################################################################################################################

######################################################################################################################################
# The final occurrences for each species was used within the "biomod2" package to delineate the calibration area for each species. The
# calibration areas for each species and the study area extent (see details in main manuscripts) can be found in the file: 
# "Crocs_shapefiles_geojson.zip". Unzip this folder and you will be able to read this shapefiles in ArcMap or Rstudio:
# ARCMAP:
A_mi <- st_read(".../Crocs_shapefiles_geojson/A_mi.geojson") # Load the file
# You can create a new file to save shapfiles. For example "Crocs_arcmap"
st_write(A_mi, ".../Crocs_shapefiles_geojson/Crocs_arcmap/A_mi.shp", delete_dsn = TRUE) # You can use the polygon in Arcmap.
# RSTUDIO:
# In R
A_mi <- st_read(".../Crocs_shapefiles_geojson/A_mi.geojson") # Load the file
st_geometry_type(A_mi) # Confirm it's a polygon
st_crs(A_mi) # Check CRS
plot(st_geometry(A_mi), col = "cyan", border = "black", main = "A. mississippiensis") # Plot only the geometry (single clean map)
######################################################################################################################################

######################################################################################################################################
# GEOGRAPHIC OVERLAP
# With the calibration areas of each species, we calculated the species pairwise geographic overlap (See details in manuscript), and
# saved the vvalues in a .csv file. This file can be found in the file: "Tables_for_R/00_Geographic_overlap_results.csv". We will use
# this file for the following code
######################################################################################################################################
data <- read.csv("...00_Geographic_overlap_results.csv", row.names = 1)

data_long <- data.frame(
  Species_1 = rep(rownames(data), times = ncol(data)),
  Species_2 = rep(colnames(data), each = nrow(data)),
  Overlap = as.vector(t(data))) # Convert the matrix into a long format

data_long$Species_2 <- gsub("\\.", " ", data_long$Species_2) # Standardize species names in Species_2 by replacing periods with spaces

species_order <- c("Crocodylus acutus", "Crocodylus intermedius", "Crocodylus rhombifer", "Crocodylus moreletii",
                   "Alligator mississippiensis", "Paleosuchus trigonatus", "Paleosuchus palpebrosus", "Caiman latirostris",
                   "Caiman crocodilus", "Caiman yacare", "Melanosuchus niger") # Define the custom order of species for Species_2

data_long$Species_2 <- factor(data_long$Species_2, levels = species_order)  # Correct order for Species_2
data_long$Species_1 <- factor(data_long$Species_1, levels = rev(species_order))  # Reverse order for Species_1

# Separate the lower (bottom) and upper (top) parts of the data
lower_data <- data_long[upper.tri(data, diag = TRUE), ]  # Lower triangle for lower
upper_data <- data_long[lower.tri(data, diag = TRUE), ]  # Upper triangle for upper

# Add a column indicating the type of data (lower or upper part of the matrix)
lower_data$Type <- "lower"
upper_data$Type <- "upper"

# Combine both datasets into a single data frame
combined_data <- bind_rows(lower_data, upper_data) 

# Filter the data for lower's D and upper I
lower_data <- subset(combined_data, Type == "lower")
upper_data <- subset(combined_data, Type == "upper")

# Adjust the data to set diagonal values of 1 as a separate group
lower_data <- lower_data %>%
  mutate(FillColor = ifelse(Species_1 == Species_2 & Overlap == 100, NA, Overlap))  # Exclude diagonal from the gradient scale

# FIGURE 3.3 This code was used for the manuscript figure.
Figure_3_3 <- ggplot(lower_data, aes(x = Species_2, y = Species_1)) +
  geom_tile(data = filter(lower_data, Species_1 == Species_2), fill = "black") +
  geom_tile(data = filter(lower_data, Species_1 != Species_2), aes(fill = FillColor)) +
  scale_fill_gradient2(midpoint = 20, low = "lightcyan", mid = "darkcyan", high = "darkcyan", na.value = "white",
                       limits = c(0, 23), oob = scales::squish) +
  geom_text(aes(label = sprintf("%.2f", Overlap)), color = "black", size = 5, fontface = "bold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 15, face = "bold.italic", color = "black", vjust = 1.5),
        axis.text.y = element_text(size = 15, face = "bold.italic", color = "black"), 
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(10, 100, 50, 10),
        legend.position = "none",  
        plot.title = element_blank()) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off")
Figure_3_3
ggsave(".../Figures/Fig_3_4.jpg", plot = Figure_3_3, width = 15, height = 10, dpi = 300)

######################################################################################################################################
# ENVIRONMENTAL OVERLAP
# With the calibration areas of each species and for the entire study extent, we extracted the 19 bioclimatic variables in ArcMap (see
# details in manuscript). This climatic variables were saved as .tif files.
######################################################################################################################################
par(mfrow=c(3,4)) # Define 3 rows and 4 columns to display the extent of predictors

###   Entire study extent for all 11 American crocodylians (calibration areas combined)
Var_Crocs <- list.files(".../Crocs/Variables", pattern = "tif$", full.names = TRUE) # Path of the folder that contains variables
Var_Crocs <- terra::rast(Var_Crocs)
plot(Var_Crocs[[1]])

Var_A_mi <- list.files(".../A_mi/Variables", pattern = "tif$", full.names = TRUE) # Alligator mississippiensis - A_mi
Var_A_mi <- terra::rast(Var_A_mi)
plot(Var_A_mi[[1]])
Var_C_ac <- list.files(".../C_ac/Variables", pattern = "tif$", full.names = TRUE) # Crocodylus acutus - C_ac  
Var_C_ac <- terra::rast(Var_C_ac)
plot(Var_C_ac[[1]])
Var_C_cr <- list.files(".../C_cr/Variables", pattern = "tif$", full.names = TRUE) # Caiman crocodilus - C_cr   
Var_C_cr <- terra::rast(Var_C_cr)
plot(Var_C_cr[[1]])
Var_C_in <- list.files(".../C_in/Variables", pattern = "tif$", full.names = TRUE) # Crocodylus intermedius - C_in
Var_C_in <- terra::rast(Var_C_in)
plot(Var_C_in[[1]])
Var_C_la <- list.files(".../C_la/Variables", pattern = "tif$", full.names = TRUE) # Caiman latirostris - C_la
Var_C_la <- terra::rast(Var_C_la)
plot(Var_C_la[[1]])
Var_C_mo <- list.files(".../C_mo/Variables", pattern = "tif$", full.names = TRUE) # Crocodylus moreletii - C_mo 
Var_C_mo <- terra::rast(Var_C_mo)
plot(Var_C_mo[[1]])
Var_C_rh <- list.files(".../C_rh/Variables", pattern = "tif$", full.names = TRUE) # Crocodylus rhombifer - C_rh
Var_C_rh <- terra::rast(Var_C_rh)
plot(Var_C_rh[[1]])
Var_C_ya <- list.files(".../C_ya/Variables", pattern = "tif$", full.names = TRUE) # Caiman yacare - C_ya
Var_C_ya <- terra::rast(Var_C_ya)
plot(Var_C_ya[[1]])
Var_M_ni <- list.files(".../M_ni/Variables", pattern = "tif$", full.names = TRUE) # Melanosuchus niger - M_ni
Var_M_ni <- terra::rast(Var_M_ni)
plot(Var_M_ni[[1]])
Var_P_pa <- list.files(".../P_pa/Variables", pattern = "tif$", full.names = TRUE) # Paleosuchus palpebrosus - P_pa
Var_P_pa <- terra::rast(Var_P_pa)
plot(Var_P_pa[[1]])
Var_P_tr <- list.files(".../P_tr/Variables", pattern = "tif$", full.names = TRUE) # Paleosuchus trigonatus - P_tr
Var_P_tr <- terra::rast(Var_P_tr)
plot(Var_P_tr[[1]])

dev.off()
######################################################################################################################################
# We will use the previous variables to extract bioclimatic values and convert them into matrices with no missing values
######################################################################################################################################
Crocs_matrix <- values(Var_Crocs) # Extract values to matrix
Crocs_matrix <- Crocs_matrix[complete.cases(Crocs_matrix), ] # Clean out missing values

A_mi_matrix <- values(Var_A_mi)
A_mi_matrix <- A_mi_matrix[complete.cases(A_mi_matrix), ]

C_ac_matrix <- values(Var_C_ac)
C_ac_matrix <- C_ac_matrix[complete.cases(C_ac_matrix), ]

C_cr_matrix <- values(Var_C_cr)
C_cr_matrix <- C_cr_matrix[complete.cases(C_cr_matrix), ]

C_in_matrix <- values(Var_C_in)
C_in_matrix <- C_in_matrix[complete.cases(C_in_matrix), ]

C_la_matrix <- values(Var_C_la)
C_la_matrix <- C_la_matrix[complete.cases(C_la_matrix), ]

C_mo_matrix <- values(Var_C_mo)
C_mo_matrix <- C_mo_matrix[complete.cases(C_mo_matrix), ]

C_rh_matrix <- values(Var_C_rh)
C_rh_matrix <- C_rh_matrix[complete.cases(C_rh_matrix), ]

C_ya_matrix <- values(Var_C_ya)
C_ya_matrix <- C_ya_matrix[complete.cases(C_ya_matrix), ]

M_ni_matrix <- values(Var_M_ni)
M_ni_matrix <- M_ni_matrix[complete.cases(M_ni_matrix), ]

P_pa_matrix <- values(Var_P_pa)
P_pa_matrix <- P_pa_matrix[complete.cases(P_pa_matrix), ]

P_tr_matrix <- values(Var_P_tr)
P_tr_matrix <- P_tr_matrix[complete.cases(P_tr_matrix), ]

######################################################################################################################################
# We will load the final occurrence datasets for each species and extract the values of the predictor variables in each location.
# NOTE 1: Clean final occurrences can be found in "2_Final_crocs_occurrences.xlsx". You need to download the excel file, and save each
# tab (e.g., A_mi_3, C_Cr_3, etc.) as an independent .csv file.
######################################################################################################################################
A_mi_occ <- read.csv(".../A_mi_3.csv") # Load the .csv file with occurrences of Alligator mississippiensis - A_mi
head(A_mi_occ) # Check the loaded file 
nrow(A_mi_occ) # Chech the number of occurrences
A_mi_occ_xy <- terra::extract(x = Var_A_mi, y = data.frame(A_mi_occ[,c('long','lat')])) # Extract predictor values in occurrences
head(A_mi_occ_xy) # Check the values extracted. Might be some NA values
A_mi_occ_var <- cbind(A_mi_occ, A_mi_occ_xy) # Merge occurrence localities and predictor values
head(A_mi_occ_var) # Check the object
nrow(A_mi_occ_var) # Check the number of occurrences
A_mi_occ_var <- A_mi_occ_var[complete.cases(A_mi_occ_var), ] # Eliminate any rows that has NA values
head(A_mi_occ_var) # Check the object after NA elimination
nrow(A_mi_occ_var) # Check the number of occurrences with predictor values. Should be lower tha previos one
write.csv(A_mi_occ_var, ".../A_mi_3_occ_var.csv") # Export and save final dataframe as a .csv file

C_ac_occ <- read.csv(".../C_ac_3.csv") # Crocodylus acutus - C_ac 
C_ac_occ_xy <- terra::extract(x = Var_C_ac, y = data.frame(C_ac_occ[,c('long','lat')]))
C_ac_occ_var <- cbind(C_ac_occ, C_ac_occ_xy) 
C_ac_occ_var <- C_ac_occ_var[complete.cases(C_ac_occ_var), ]
write.csv(C_ac_occ_var, ".../C_ac_3_occ_var.csv")

C_cr_occ <- read.csv(".../C_cr_3.csv") # Caiman crocodilus - C_cr 
C_cr_occ_xy <- terra::extract(x = Var_C_cr, y = data.frame(C_cr_occ[,c('long','lat')])) 
C_cr_occ_var <- cbind(C_cr_occ, C_cr_occ_xy) 
C_cr_occ_var <- C_cr_occ_var[complete.cases(C_cr_occ_var), ]
write.csv(C_cr_occ_var, ".../C_cr_3_occ_var.csv")

C_in_occ <- read.csv(".../C_in_3.csv") # Crocodylus intermedius - C_in
C_in_occ_xy <- terra::extract(x = Var_C_in, y = data.frame(C_in_occ[,c('long','lat')])) 
C_in_occ_var <- cbind(C_in_occ, C_in_occ_xy) 
C_in_occ_var <- C_in_occ_var[complete.cases(C_in_occ_var), ]
write.csv(C_in_occ_var, ".../C_in_3_occ_var.csv")

C_la_occ <- read.csv(".../C_la_3.csv") # Caiman latirostris - C_la   
C_la_occ_xy <- terra::extract(x = Var_C_la, y = data.frame(C_la_occ[,c('long','lat')])) 
C_la_occ_var <- cbind(C_la_occ, C_la_occ_xy) 
C_la_occ_var <- C_la_occ_var[complete.cases(C_la_occ_var), ]
write.csv(C_la_occ_var, ".../C_la_3_occ_var.csv")

C_mo_occ <- read.csv(".../C_mo_3.csv") # Crocodylus moreletii - C_mo
C_mo_occ_xy <- terra::extract(x = Var_C_mo, y = data.frame(C_mo_occ[,c('long','lat')])) 
C_mo_occ_var <- cbind(C_mo_occ, C_mo_occ_xy) 
C_mo_occ_var <- C_mo_occ_var[complete.cases(C_mo_occ_var), ]
write.csv(C_mo_occ_var, ".../C_mo_3_occ_var.csv")

C_rh_occ <- read.csv(".../C_rh_3.csv") # Crocodylus rhombifer - C_rh 
C_rh_occ_xy <- terra::extract(x = Var_C_rh, y = data.frame(C_rh_occ[,c('long','lat')])) 
C_rh_occ_var <- cbind(C_rh_occ, C_rh_occ_xy) 
C_rh_occ_var <- C_rh_occ_var[complete.cases(C_rh_occ_var), ]
write.csv(C_rh_occ_var, ".../C_rh_3_occ_var.csv")

C_ya_occ <- read.csv(".../C_ya_3.csv") # Caiman yacare - C_ya
C_ya_occ_xy <- terra::extract(x = Var_C_ya, y = data.frame(C_ya_occ[,c('long','lat')])) 
C_ya_occ_var <- cbind(C_ya_occ, C_ya_occ_xy) 
C_ya_occ_var <- C_ya_occ_var[complete.cases(C_ya_occ_var), ]
write.csv(C_ya_occ_var, ".../C_ya_3_occ_var.csv")

M_ni_occ <- read.csv(".../M_ni_3.csv") # Melanosuchus niger - M_ni 
M_ni_occ_xy <- terra::extract(x = Var_M_ni, y = data.frame(M_ni_occ[,c('long','lat')])) 
M_ni_occ_var <- cbind(M_ni_occ, M_ni_occ_xy) 
M_ni_occ_var <- M_ni_occ_var[complete.cases(M_ni_occ_var), ]
write.csv(M_ni_occ_var, ".../M_ni_3_occ_var.csv")

P_pa_occ <- read.csv(".../P_pa_3.csv") # Paleosuchus palpebrosus - P_pa
P_pa_occ_xy <- terra::extract(x = Var_P_pa, y = data.frame(P_pa_occ[,c('long','lat')])) 
P_pa_occ_var <- cbind(P_pa_occ, P_pa_occ_xy) 
P_pa_occ_var <- P_pa_occ_var[complete.cases(P_pa_occ_var), ]
write.csv(P_pa_occ_var, ".../P_pa_3_occ_var.csv")

P_tr_occ <- read.csv(".../P_tr_3.csv") # Paleosuchus trigonatus - P_tr
P_tr_occ_xy <- terra::extract(x = Var_P_tr, y = data.frame(P_tr_occ[,c('long','lat')])) 
P_tr_occ_var <- cbind(P_tr_occ, P_tr_occ_xy) 
P_tr_occ_var <- P_tr_occ_var[complete.cases(P_tr_occ_var), ]
write.csv(P_tr_occ_var, ".../P_tr_3_occ_var.csv")

######################################################################################################################################
# We will load the final occurrences with the values of the predictors. Files can be found in "3_Final_crocs_occ_with_variables.xlsx".
# You need to download the excel file, and save each tab (e.g., A_mi_occ_var, C_Cr_occ_var, etc.) as an independent .csv file.
######################################################################################################################################
A_mi_occ_var <- read.csv("D:/THESIS_CH3/Species/A_mi/Occ_tab/A_mi_3_occ_var.csv", row.names = 1) # Load the .csv file
C_ac_occ_var <- read.csv("D:/THESIS_CH3/Species/C_ac/Occ_tab/C_ac_3_occ_var.csv", row.names = 1) 
C_cr_occ_var <- read.csv("D:/THESIS_CH3/Species/C_cr/Occ_tab/C_cr_3_occ_var.csv", row.names = 1)
C_in_occ_var <- read.csv("D:/THESIS_CH3/Species/C_in/Occ_tab/C_in_3_occ_var.csv", row.names = 1)
C_la_occ_var <- read.csv("D:/THESIS_CH3/Species/C_la/Occ_tab/C_la_3_occ_var.csv", row.names = 1)
C_mo_occ_var <- read.csv("D:/THESIS_CH3/Species/C_mo/Occ_tab/C_mo_3_occ_var.csv", row.names = 1)
C_rh_occ_var <- read.csv("D:/THESIS_CH3/Species/C_rh/Occ_tab/C_rh_3_occ_var.csv", row.names = 1)
C_ya_occ_var <- read.csv("D:/THESIS_CH3/Species/C_ya/Occ_tab/C_ya_3_occ_var.csv", row.names = 1)
M_ni_occ_var <- read.csv("D:/THESIS_CH3/Species/M_ni/Occ_tab/M_ni_3_occ_var.csv", row.names = 1)
P_pa_occ_var <- read.csv("D:/THESIS_CH3/Species/P_pa/Occ_tab/P_pa_3_occ_var.csv", row.names = 1)
P_tr_occ_var <- read.csv("D:/THESIS_CH3/Species/P_tr/Occ_tab/P_tr_3_occ_var.csv", row.names = 1)

######################################################################################################################################
# NICHE QUANTIFICATION VIA PRINCIPAL COMPONENT ANALYSIS - PCA
# We will do a PCA of the environmental data of the entire study extent, each species range, including occurrence densities
######################################################################################################################################
# Entire study extent (All crocodylians environment)
Crocs_pca.env <- dudi.pca(Crocs_matrix, center = TRUE, scale = TRUE, scannf = F, nf = 2)
ecospat.plot.contrib(contrib = Crocs_pca.env$co, eigen = Crocs_pca.env$eig)
Crocs_scores <- Crocs_pca.env$li # PCA scores for the whole study area

# Plot the PCA axes as lines passing through the origin
plot(0, 0, type = "n", 
     xlab = paste("PC1 (", round(Crocs_pca.env$eig[1] / sum(Crocs_pca.env$eig) * 100, 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(Crocs_pca.env$eig[2] / sum(Crocs_pca.env$eig) * 100, 2), "%)", sep = ""),
     xlim = c(-1, 1), ylim = c(-1, 1), 
     main = NULL)  # Remove the title

abline(h = 0, v = 0, col = "darkgray", lwd = 1) # Draw the PCA axes (PC1 and PC2) as lines through the origin
arrows(0, 0, Crocs_pca.env$co[, 1], Crocs_pca.env$co[, 2], col = "black", length = 0.1, angle = 15, lwd = 1)  # Black thick arrows
text(Crocs_pca.env$co[, 1] * 1.1, Crocs_pca.env$co[, 2] * 1.1, labels = rownames(Crocs_pca.env$co), col = "bakc", cex = 0.9, font = 2)

# Species occurrences
scores_A_mi <- suprow(Crocs_pca.env, A_mi_occ_var[which(A_mi_occ_var[,3]==1),5:23])$li # (column 3 = 1 [presences])
head(scores_A_mi) # check that there is an object with 2 columns "Axis1" and "Axis2"
print(dim(A_mi_occ_var))  # Shows number of rows and columns
print(names(A_mi_occ_var))  # Shows column names

scores_C_ac <- suprow(Crocs_pca.env, C_ac_occ_var[which(C_ac_occ_var[,3]==1),5:23])$li
scores_C_cr <- suprow(Crocs_pca.env, C_cr_occ_var[which(C_cr_occ_var[,3]==1),5:23])$li
scores_C_in <- suprow(Crocs_pca.env, C_in_occ_var[which(C_in_occ_var[,3]==1),5:23])$li
scores_C_la <- suprow(Crocs_pca.env, C_la_occ_var[which(C_la_occ_var[,3]==1),5:23])$li
scores_C_mo <- suprow(Crocs_pca.env, C_mo_occ_var[which(C_mo_occ_var[,3]==1),5:23])$li
scores_C_rh <- suprow(Crocs_pca.env, C_rh_occ_var[which(C_rh_occ_var[,3]==1),5:23])$li
scores_C_ya <- suprow(Crocs_pca.env, C_ya_occ_var[which(C_ya_occ_var[,3]==1),5:23])$li
scores_M_ni <- suprow(Crocs_pca.env, M_ni_occ_var[which(M_ni_occ_var[,3]==1),5:23])$li
scores_P_pa <- suprow(Crocs_pca.env, P_pa_occ_var[which(P_pa_occ_var[,3]==1),5:23])$li
scores_P_tr <- suprow(Crocs_pca.env, P_tr_occ_var[which(P_tr_occ_var[,3]==1),5:23])$li

# Species ranges environment
scores_A_mi_dist <- suprow(Crocs_pca.env, A_mi_matrix)$li
head(scores_A_mi_dist) # check that there is an object with 2 columns "Axis1" and "Axis2"

scores_C_ac_dist <- suprow(Crocs_pca.env, C_ac_matrix)$li
scores_C_cr_dist <- suprow(Crocs_pca.env, C_cr_matrix)$li
scores_C_in_dist <- suprow(Crocs_pca.env, C_in_matrix)$li
scores_C_la_dist <- suprow(Crocs_pca.env, C_la_matrix)$li
scores_C_mo_dist <- suprow(Crocs_pca.env, C_mo_matrix)$li
scores_C_rh_dist <- suprow(Crocs_pca.env, C_rh_matrix)$li
scores_C_ya_dist <- suprow(Crocs_pca.env, C_ya_matrix)$li
scores_M_ni_dist <- suprow(Crocs_pca.env, M_ni_matrix)$li
scores_P_pa_dist <- suprow(Crocs_pca.env, P_pa_matrix)$li
scores_P_tr_dist <- suprow(Crocs_pca.env, P_tr_matrix)$li

######################################################################################################################################
# GRIDDING NICHES: PLOT THE NATIVE NICHES OF EACH SPECIES
######################################################################################################################################
grid_clim_A_mi <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_A_mi_dist, sp = scores_A_mi, R=300, th.sp = 0)
grid_clim_C_ac <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_ac_dist, sp = scores_C_ac, R=300, th.sp = 0)
grid_clim_C_cr <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_cr_dist, sp = scores_C_cr, R=300, th.sp = 0)
grid_clim_C_in <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_in_dist, sp = scores_C_in, R=300, th.sp = 0)
grid_clim_C_la <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_la_dist, sp = scores_C_la, R=300, th.sp = 0)
grid_clim_C_mo <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_mo_dist, sp = scores_C_mo, R=300, th.sp = 0)
grid_clim_C_rh <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_rh_dist, sp = scores_C_rh, R=300, th.sp = 0)
grid_clim_C_ya <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_C_ya_dist, sp = scores_C_ya, R=300, th.sp = 0)
grid_clim_M_ni <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_M_ni_dist, sp = scores_M_ni, R=300, th.sp = 0)
grid_clim_P_pa <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_P_pa_dist, sp = scores_P_pa, R=300, th.sp = 0)
grid_clim_P_tr <- ecospat.grid.clim.dyn(glob = Crocs_scores, glob1 = scores_P_tr_dist, sp = scores_P_tr, R=300, th.sp = 0)
# Run this codes to view the plots
ecospat.plot.niche(grid_clim_C_ac, title=expression(bolditalic("C. acutus")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_C_in, title=expression(bolditalic("C. intermedius")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_C_rh, title=expression(bolditalic("C. rhombifer")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_C_mo, title=expression(bolditalic("C. moreletii")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_A_mi, title=expression(bolditalic("A. mississippiensis")), name.axis1= "", name.axis2 = "PC2", cor = F)
ecospat.plot.niche(grid_clim_P_tr, title=expression(bolditalic("P. trigonatus")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_P_pa, title=expression(bolditalic("P. palpebrosus")), name.axis1= "", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_C_la, title=expression(bolditalic("C. latirostris")), name.axis1= "PC1", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_C_cr, title=expression(bolditalic("C. crocodilus")), name.axis1= "PC1", name.axis2 = "PC2", cor = F)
ecospat.plot.niche(grid_clim_C_ya, title=expression(bolditalic("C. yacare")), name.axis1= "PC1", name.axis2 = "", cor = F)
ecospat.plot.niche(grid_clim_M_ni, title=expression(bolditalic("M. niger")), name.axis1= "PC1", name.axis2 = "", cor = F)

######################################################################################################################################
# CALCULATE THE NICHE BREADTH USING ENMTools OF EACH OF THE 11 SPECIES USING THE LEVINS INVERSE CONCENTRATION NICHE BREADTH 
######################################################################################################################################
species_names <- c("A_mi", "C_ac", "C_cr", "C_in", "C_la", "C_mo", "C_rh", "C_ya", "M_ni", "P_pa", "P_tr") # Define species names

niche_breadth_values <- data.frame(Species = species_names, B1 = NA, B2 = NA) # Create an empty data frame to store the results

# Calculate niche breadth for each species and store the results
niche_breadth_values$B1[1] <- raster.breadth(grid_clim_A_mi$z)$B1
niche_breadth_values$B2[1] <- raster.breadth(grid_clim_A_mi$z)$B2

niche_breadth_values$B1[2] <- raster.breadth(grid_clim_C_ac$z)$B1
niche_breadth_values$B2[2] <- raster.breadth(grid_clim_C_ac$z)$B2

niche_breadth_values$B1[3] <- raster.breadth(grid_clim_C_cr$z)$B1
niche_breadth_values$B2[3] <- raster.breadth(grid_clim_C_cr$z)$B2

niche_breadth_values$B1[4] <- raster.breadth(grid_clim_C_in$z)$B1
niche_breadth_values$B2[4] <- raster.breadth(grid_clim_C_in$z)$B2

niche_breadth_values$B1[5] <- raster.breadth(grid_clim_C_la$z)$B1
niche_breadth_values$B2[5] <- raster.breadth(grid_clim_C_la$z)$B2

niche_breadth_values$B1[6] <- raster.breadth(grid_clim_C_mo$z)$B1
niche_breadth_values$B2[6] <- raster.breadth(grid_clim_C_mo$z)$B2

niche_breadth_values$B1[7] <- raster.breadth(grid_clim_C_rh$z)$B1
niche_breadth_values$B2[7] <- raster.breadth(grid_clim_C_rh$z)$B2

niche_breadth_values$B1[8] <- raster.breadth(grid_clim_C_ya$z)$B1
niche_breadth_values$B2[8] <- raster.breadth(grid_clim_C_ya$z)$B2

niche_breadth_values$B1[9] <- raster.breadth(grid_clim_M_ni$z)$B1
niche_breadth_values$B2[9] <- raster.breadth(grid_clim_M_ni$z)$B2

niche_breadth_values$B1[10] <- raster.breadth(grid_clim_P_pa$z)$B1
niche_breadth_values$B2[10] <- raster.breadth(grid_clim_P_pa$z)$B2

niche_breadth_values$B1[11] <- raster.breadth(grid_clim_P_tr$z)$B1
niche_breadth_values$B2[11] <- raster.breadth(grid_clim_P_tr$z)$B2

# Save the results as a CSV file
write.csv(niche_breadth_values, ".../01_Niche_breadth_results.csv", row.names = FALSE)

######################################################################################################################################
# Relationship between geographic range and niche breadth. We will use the file "Tables_for_R/Range_breadth.csv".
# This file has 2 columns: Go: Geographic occupied range; NB1: Levins' niche breadth.
########################################################################################################################
n_breadth <- read.csv("../Tables_for_R/Range_breadth.csv", row.names = 1)
n_breadth$log_Go <- log10(n_breadth$Go) # Log-transform Go (before correlation and regression)

# Including 11 species
cor_test_log <- cor.test(n_breadth$log_Go, n_breadth$NB1, method = "pearson") # Pearson correlation on log-transformed Go and NB1
cor_test_log$estimate

# Not including C. rhombifer
n_breadth_no_outlier <- n_breadth[n_breadth$species != "Crocodylus rhombifer", ]
n_breadth_no_outlier$log_Go <- log10(n_breadth_no_outlier$Go)
cor_test_no_outlier <- cor.test(n_breadth_no_outlier$log_Go, n_breadth_no_outlier$NB1, method = "pearson")
cor_test_no_outlier$estimate

# Pearson correlation coefficients
r_all_log <- cor_test_log$estimate  # r for all species
r_no_outlier_log <- cor_test_no_outlier$estimate  # r for excluding Crocodylus rhombifer

# Create equation labels for both regressions (including r)
equation_label_all_log <- paste("y = ", round(coef(lm(NB1 ~ log(Go), data = n_breadth))[[2]], 2), 
                                "x + ", round(coef(lm(NB1 ~ log(Go), data = n_breadth))[[1]], 2), sep="")
equation_label_no_outlier_log <- paste("y = ", round(coef(lm(NB1 ~ log(Go), data = n_breadth_no_outlier))[[2]], 2), 
                                       "x + ", round(coef(lm(NB1 ~ log(Go), data = n_breadth_no_outlier))[[1]], 2), sep="")

# Create R² and p-value labels
r2_label_all_log <- paste("R² = ", round(summary(lm(NB1 ~ log(Go), data = n_breadth))$r.squared, 3), sep="")
r2_label_no_outlier_log <- paste("R² = ", round(summary(lm(NB1 ~ log(Go), data = n_breadth_no_outlier))$r.squared, 3), sep="")

p_value_label_all_log <- paste("p = ",
                               format(round(summary(lm(NB1 ~ log(Go),data = n_breadth))$coefficients[2, 4], 4),scientific=F),sep = "")
p_value_label_no_outlier_log <- paste("p = ",
                                      round(summary(lm(NB1 ~ log(Go), data = n_breadth_no_outlier))$coefficients[2, 4],4),sep="")

# Pearson correlation coefficient labels (in italics)
r_label_all_log <- paste("r = ", round(r_all_log, 3), sep="")
r_label_no_outlier_log <- paste("r = ", round(r_no_outlier_log, 3), sep="")

# Figure 3.5. This code was used to create the manuscript figure.
Figure_3_5 <- ggplot(n_breadth, aes(x = Go, y = NB1)) +
  geom_smooth(data = n_breadth, aes(x = Go, y = NB1), method = "lm", se = FALSE, color = "darkcyan", 
              linetype = "dashed", size = 0.6, alpha = 0.5, inherit.aes = FALSE) +  # Regression line behind
  geom_smooth(data = n_breadth_no_outlier, aes(x = Go, y = NB1), method = "lm", se = FALSE, color = "gray50", 
              linetype = "dashed", size = 0.6, alpha = 0.5, inherit.aes = FALSE) +  # Same settings for the blue line
  geom_point(size = 4) +
  geom_text_repel(aes(label = species), size = 6, box.padding = 0.8, max.overlaps = Inf, segment.size = 0.2, 
                  family = "Times New Roman", fontface = "italic", fill = NA, color = "black") +
  annotate("text", x = 3000, y = 1, label = equation_label_all_log, size = 7, color = "darkcyan", fontface = "italic") +
  annotate("text", x = 3000, y = 0.96, label = r2_label_all_log, size = 7, color = "darkcyan", fontface = "italic") +
  annotate("text", x = 3000, y = 0.92, label = r_label_all_log, size = 7, color = "darkcyan", fontface = "italic") +  
  annotate("text", x = 3000, y = 0.88, label = p_value_label_all_log, size = 7, color = "darkcyan", fontface = "italic") +
  annotate("text", x = 3000, y = 0.78, label = equation_label_no_outlier_log, size = 7, color = "gray50", fontface = "italic") +
  annotate("text", x = 3000, y = 0.74, label = r2_label_no_outlier_log, size = 7, color = "gray50", fontface = "italic") +
  annotate("text", x = 3000, y = 0.7, label = r_label_no_outlier_log, size = 7, color = "gray50", fontface = "italic") +  
  annotate("text", x = 3000, y = 0.66, label = p_value_label_no_outlier_log, size = 7, color = "gray50", fontface = "italic") +
  labs(x = "Geographic range (km²)", 
       y = bquote(bold("Niche breadth - Standardized inverse concentration metric (Levins' " ~ B[1] * ")"))) +
  theme_bw() +
  theme(axis.line = element_line(color = "black"), panel.grid = element_blank(), plot.margin = margin(10, 10, 10, 10),
        axis.title = element_text(size = 17, face = "bold"), 
        axis.title.y = element_text(margin = margin(r = 10), vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 10), hjust = 0.5),
        panel.border = element_blank(), axis.text = element_text(size = 13)) +
  coord_cartesian(xlim = c(1000, 10000000), ylim = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "b")
plot(Figure_3_5)
ggsave(".../Figures/Fig_3_5.jpg", plot = plot_Go_NB1, width = 15, height = 10, dpi = 500)

######################################################################################################################################
# CALCULATE THE SCHOENER'S D - INDEX OF NICHE OVERLAP
######################################################################################################################################
# Alligator mississippiensis  vs  other 10 species
D.A_mi_C_ac <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_ac, cor=T)$D
D.A_mi_C_cr <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_cr, cor=T)$D
D.A_mi_C_in <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_in, cor=T)$D
D.A_mi_C_la <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_la, cor=T)$D
D.A_mi_C_mo <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_mo, cor=T)$D
D.A_mi_C_rh <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_rh, cor=T)$D
D.A_mi_C_ya <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_ya, cor=T)$D
D.A_mi_M_ni <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_M_ni, cor=T)$D
D.A_mi_P_pa <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_P_pa, cor=T)$D
D.A_mi_P_tr <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_P_tr, cor=T)$D
A_mi_list_D <- data.frame(A_mi_value = c(1, D.A_mi_C_ac, D.A_mi_C_cr, D.A_mi_C_in, D.A_mi_C_la, D.A_mi_C_mo,
                                       D.A_mi_C_rh, D.A_mi_C_ya, D.A_mi_M_ni, D.A_mi_P_pa, D.A_mi_P_tr))
# Crocodylus acutus  vs  other 9 species
D.C_ac_C_cr <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_cr, cor=T)$D
D.C_ac_C_in <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_in, cor=T)$D
D.C_ac_C_la <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_la, cor=T)$D
D.C_ac_C_mo <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_mo, cor=T)$D
D.C_ac_C_rh <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_rh, cor=T)$D
D.C_ac_C_ya <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_ya, cor=T)$D
D.C_ac_M_ni <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_M_ni, cor=T)$D
D.C_ac_P_pa <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_P_pa, cor=T)$D
D.C_ac_P_tr <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_P_tr, cor=T)$D
C_ac_list_D <- data.frame(C_ac_value = c(1,1,D.C_ac_C_cr, D.C_ac_C_in, D.C_ac_C_la, D.C_ac_C_mo, D.C_ac_C_rh,
                                       D.C_ac_C_ya, D.C_ac_M_ni, D.C_ac_P_pa, D.C_ac_P_tr))
# Caiman crocodilus  vs  other 8 species
D.C_cr_C_in <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_in, cor=T)$D
D.C_cr_C_la <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_la, cor=T)$D
D.C_cr_C_mo <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_mo, cor=T)$D
D.C_cr_C_rh <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_rh, cor=T)$D
D.C_cr_C_ya <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_ya, cor=T)$D
D.C_cr_M_ni <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_M_ni, cor=T)$D
D.C_cr_P_pa <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_P_pa, cor=T)$D
D.C_cr_P_tr <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_P_tr, cor=T)$D
C_cr_list_D <- data.frame(C_cr_value = c(1,1,1, D.C_cr_C_in, D.C_cr_C_la, D.C_cr_C_mo, D.C_cr_C_rh, D.C_cr_C_ya,
                                       D.C_cr_M_ni, D.C_cr_P_pa, D.C_cr_P_tr))
# Crocodylus intermedius  vs  other 7 species
D.C_in_C_la <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_la, cor=T)$D
D.C_in_C_mo <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_mo, cor=T)$D
D.C_in_C_rh <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_rh, cor=T)$D
D.C_in_C_ya <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_ya, cor=T)$D
D.C_in_M_ni <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_M_ni, cor=T)$D
D.C_in_P_pa <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_P_pa, cor=T)$D
D.C_in_P_tr <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_P_tr, cor=T)$D
C_in_list_D <- data.frame(C_in_values = c(1,1,1,1, D.C_in_C_la, D.C_in_C_mo, D.C_in_C_rh, D.C_in_C_ya, D.C_in_M_ni,
                                        D.C_in_P_pa, D.C_in_P_tr))
# Caiman latirostris  vs  other 6 species
D.C_la_C_mo <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_mo, cor=T)$D
D.C_la_C_rh <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_rh, cor=T)$D
D.C_la_C_ya <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_ya, cor=T)$D
D.C_la_M_ni <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_M_ni, cor=T)$D
D.C_la_P_pa <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_P_pa, cor=T)$D
D.C_la_P_tr <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_P_tr, cor=T)$D
C_la_list_D <- data.frame(C_la_values = c(1,1,1,1,1, D.C_la_C_mo, D.C_la_C_rh, D.C_la_C_ya, D.C_la_M_ni, D.C_la_P_pa,
                                        D.C_la_P_tr))
# Crocodylus moreletii  vs  other 5 species
D.C_mo_C_rh <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_C_rh, cor=T)$D
D.C_mo_C_ya <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_C_ya, cor=T)$D
D.C_mo_M_ni <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_M_ni, cor=T)$D
D.C_mo_P_pa <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_P_pa, cor=T)$D
D.C_mo_P_tr <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_P_tr, cor=T)$D
C_mo_list_D <- data.frame(C_mo_values = c(1,1,1,1,1,1, D.C_mo_C_rh, D.C_mo_C_ya, D.C_mo_M_ni, D.C_mo_P_pa, D.C_mo_P_tr))
# Crocodylus rhombifer  vs  other 4 species
D.C_rh_C_ya <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_C_ya, cor=T)$D
D.C_rh_M_ni <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_M_ni, cor=T)$D
D.C_rh_P_pa <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_P_pa, cor=T)$D
D.C_rh_P_tr <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_P_tr, cor=T)$D
C_rh_list_D <- data.frame(C_rh_values = c(1,1,1,1,1,1,1, D.C_rh_C_ya, D.C_rh_M_ni, D.C_rh_P_pa, D.C_rh_P_tr))
# Caiman yacare  vs  other 3 species
D.C_ya_M_ni <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_M_ni, cor=T)$D
D.C_ya_P_pa <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_P_pa, cor=T)$D
D.C_ya_P_tr <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_P_tr, cor=T)$D
C_ya_list_D <- data.frame(C_ya_values = c(1,1,1,1,1,1,1,1, D.C_ya_M_ni, D.C_ya_P_pa, D.C_ya_P_tr))
# Melanosuchus niger  vs  other 2 species
D.M_ni_P_pa <- ecospat.niche.overlap (grid_clim_M_ni, grid_clim_P_pa, cor=T)$D
D.M_ni_P_tr <- ecospat.niche.overlap (grid_clim_M_ni, grid_clim_P_tr, cor=T)$D
M_ni_list_D <- data.frame(M_ni_values = c(1,1,1,1,1,1,1,1,1, D.M_ni_P_pa, D.M_ni_P_tr))
# Paleosuchus palpebrosus  vs  other 1 species
D.P_pa_P_tr <- ecospat.niche.overlap (grid_clim_P_pa, grid_clim_P_tr, cor=T)$D
P_pa_list_D <- data.frame(P_pa_values = c(1,1,1,1,1,1,1,1,1,1, D.P_pa_P_tr))

Niche_overlap_crocs_D <- data.frame(A_mi_list_D, C_ac_list_D, C_cr_list_D, C_in_list_D, C_la_list_D, C_mo_list_D,
                                   C_rh_list_D, C_ya_list_D, M_ni_list_D, P_pa_list_D)
write.csv(Niche_overlap_crocs_D, ".../02_Niche_overlap_Schoener_D_corr.csv")

######################################################################################################################################
# CALCULATE THE HELLINGER'S I - INDEX OF NICHE OVERLAP
######################################################################################################################################
# Alligator mississippiensis  vs  other 10 species
I.A_mi_C_ac <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_ac, cor=T)$I
I.A_mi_C_cr <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_cr, cor=T)$I
I.A_mi_C_in <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_in, cor=T)$I
I.A_mi_C_la <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_la, cor=T)$I
I.A_mi_C_mo <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_mo, cor=T)$I
I.A_mi_C_rh <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_rh, cor=T)$I
I.A_mi_C_ya <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_C_ya, cor=T)$I
I.A_mi_M_ni <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_M_ni, cor=T)$I
I.A_mi_P_pa <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_P_pa, cor=T)$I
I.A_mi_P_tr <- ecospat.niche.overlap (grid_clim_A_mi, grid_clim_P_tr, cor=T)$I
A_mi_list_I <- data.frame(A_mi_value = c(1, I.A_mi_C_ac, I.A_mi_C_cr, I.A_mi_C_in, I.A_mi_C_la, I.A_mi_C_mo,
                                       I.A_mi_C_rh, I.A_mi_C_ya, I.A_mi_M_ni, I.A_mi_P_pa, I.A_mi_P_tr))
# Crocodylus acutus  vs  other 9 species
I.C_ac_C_cr <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_cr, cor=T)$I
I.C_ac_C_in <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_in, cor=T)$I
I.C_ac_C_la <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_la, cor=T)$I
I.C_ac_C_mo <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_mo, cor=T)$I
I.C_ac_C_rh <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_rh, cor=T)$I
I.C_ac_C_ya <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_C_ya, cor=T)$I
I.C_ac_M_ni <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_M_ni, cor=T)$I
I.C_ac_P_pa <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_P_pa, cor=T)$I
I.C_ac_P_tr <- ecospat.niche.overlap (grid_clim_C_ac, grid_clim_P_tr, cor=T)$I
C_ac_list_I <- data.frame(C_ac_value = c(1,1,I.C_ac_C_cr, I.C_ac_C_in, I.C_ac_C_la, I.C_ac_C_mo, I.C_ac_C_rh,
                                       I.C_ac_C_ya, I.C_ac_M_ni, I.C_ac_P_pa, I.C_ac_P_tr))
# Caiman crocodilus  vs  other 8 species
I.C_cr_C_in <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_in, cor=T)$I
I.C_cr_C_la <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_la, cor=T)$I
I.C_cr_C_mo <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_mo, cor=T)$I
I.C_cr_C_rh <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_rh, cor=T)$I
I.C_cr_C_ya <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_C_ya, cor=T)$I
I.C_cr_M_ni <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_M_ni, cor=T)$I
I.C_cr_P_pa <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_P_pa, cor=T)$I
I.C_cr_P_tr <- ecospat.niche.overlap (grid_clim_C_cr, grid_clim_P_tr, cor=T)$I
C_cr_list_I <- data.frame(C_cr_value = c(1,1,1, I.C_cr_C_in, I.C_cr_C_la, I.C_cr_C_mo, I.C_cr_C_rh, I.C_cr_C_ya,
                                       I.C_cr_M_ni, I.C_cr_P_pa, I.C_cr_P_tr))
# Crocodylus intermedius  vs  other 7 species
I.C_in_C_la <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_la, cor=T)$I
I.C_in_C_mo <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_mo, cor=T)$I
I.C_in_C_rh <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_rh, cor=T)$I
I.C_in_C_ya <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_C_ya, cor=T)$I
I.C_in_M_ni <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_M_ni, cor=T)$I
I.C_in_P_pa <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_P_pa, cor=T)$I
I.C_in_P_tr <- ecospat.niche.overlap (grid_clim_C_in, grid_clim_P_tr, cor=T)$I
C_in_list_I <- data.frame(C_in_values = c(1,1,1,1, I.C_in_C_la, I.C_in_C_mo, I.C_in_C_rh, I.C_in_C_ya, I.C_in_M_ni,
                                        I.C_in_P_pa, I.C_in_P_tr))
# Caiman latirostris  vs  other 6 species
I.C_la_C_mo <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_mo, cor=T)$I
I.C_la_C_rh <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_rh, cor=T)$I
I.C_la_C_ya <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_C_ya, cor=T)$I
I.C_la_M_ni <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_M_ni, cor=T)$I
I.C_la_P_pa <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_P_pa, cor=T)$I
I.C_la_P_tr <- ecospat.niche.overlap (grid_clim_C_la, grid_clim_P_tr, cor=T)$I
C_la_list_I <- data.frame(C_la_values = c(1,1,1,1,1, I.C_la_C_mo, I.C_la_C_rh, I.C_la_C_ya, I.C_la_M_ni, I.C_la_P_pa,
                                        I.C_la_P_tr))
# Crocodylus moreletii  vs  other 5 species
I.C_mo_C_rh <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_C_rh, cor=T)$I
I.C_mo_C_ya <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_C_ya, cor=T)$I
I.C_mo_M_ni <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_M_ni, cor=T)$I
I.C_mo_P_pa <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_P_pa, cor=T)$I
I.C_mo_P_tr <- ecospat.niche.overlap (grid_clim_C_mo, grid_clim_P_tr, cor=T)$I
C_mo_list_I <- data.frame(C_mo_values = c(1,1,1,1,1,1, I.C_mo_C_rh, I.C_mo_C_ya, I.C_mo_M_ni, I.C_mo_P_pa, I.C_mo_P_tr))
# Crocodylus rhombifer  vs  other 4 species
I.C_rh_C_ya <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_C_ya, cor=T)$I
I.C_rh_M_ni <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_M_ni, cor=T)$I
I.C_rh_P_pa <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_P_pa, cor=T)$I
I.C_rh_P_tr <- ecospat.niche.overlap (grid_clim_C_rh, grid_clim_P_tr, cor=T)$I
C_rh_list_I <- data.frame(C_rh_values = c(1,1,1,1,1,1,1, I.C_rh_C_ya, I.C_rh_M_ni, I.C_rh_P_pa, I.C_rh_P_tr))
# Caiman yacare  vs  other 3 species
I.C_ya_M_ni <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_M_ni, cor=T)$I
I.C_ya_P_pa <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_P_pa, cor=T)$I
I.C_ya_P_tr <- ecospat.niche.overlap (grid_clim_C_ya, grid_clim_P_tr, cor=T)$I
C_ya_list_I <- data.frame(C_ya_values = c(1,1,1,1,1,1,1,1, I.C_ya_M_ni, I.C_ya_P_pa, I.C_ya_P_tr))
# Melanosuchus niger  vs  other 2 species
I.M_ni_P_pa <- ecospat.niche.overlap (grid_clim_M_ni, grid_clim_P_pa, cor=T)$I
I.M_ni_P_tr <- ecospat.niche.overlap (grid_clim_M_ni, grid_clim_P_tr, cor=T)$I
M_ni_list_I <- data.frame(M_ni_values = c(1,1,1,1,1,1,1,1,1, I.M_ni_P_pa, I.M_ni_P_tr))
# Paleosuchus palpebrosus  vs  other 1 species
I.P_pa_P_tr <- ecospat.niche.overlap (grid_clim_P_pa, grid_clim_P_tr, cor=T)$I
P_pa_list_I <- data.frame(P_pa_values = c(1,1,1,1,1,1,1,1,1,1, I.P_pa_P_tr))

Niche_overlap_crocs_I <- data.frame(A_mi_list_I, C_ac_list_I, C_cr_list_I, C_in_list_I, C_la_list_I, C_mo_list_I,
                                  C_rh_list_I, C_ya_list_I, M_ni_list_I, P_pa_list_I)
write.csv(Niche_overlap_crocs_I, ".../03_Niche_overlap_Hellinger_I_corr.csv")

################################################################################################################################
# CORRELATION PLOTS FOR SCHOENER'S D AND HELLINGER'S I
# We have combined the data saved in files "02_Niche_overlap_Schoener_D_corr.csv" and "03_Niche_overlap_Hellinger_I_corr.csv".
# The data combined was saved in the file "04_Correlation_plot_D_and_I.csv", which will be used to prepare Figure 3.6.  
################################################################################################################################
data <- read.csv(".../Tables_for_R/04_Correlation_plot_D_and_I.csv", row.names = 1) # Load the data

data_long <- data.frame(Species_1 = rep(rownames(data), times = ncol(data)),
                        Species_2 = rep(colnames(data), each = nrow(data)),
                        Overlap = as.vector(t(data))) # Convert the matrix into a long format

data_long$Species_2 <- gsub("\\.", " ", data_long$Species_2) # Standardize species names in Species_2 by replacing periods with spaces

# Define the custom order of species for Species_2
species_order <- c("Crocodylus acutus", "Crocodylus intermedius", "Crocodylus rhombifer", "Crocodylus moreletii",
                   "Alligator mississippiensis", "Paleosuchus trigonatus", "Paleosuchus palpebrosus", "Caiman latirostris",
                   "Caiman crocodilus", "Caiman yacare", "Melanosuchus niger")

# Reorder the species in both Species_1 and Species_2
data_long$Species_2 <- factor(data_long$Species_2, levels = species_order)  # Correct order for Species_2
data_long$Species_1 <- factor(data_long$Species_1, levels = rev(species_order))  # Reverse order for Species_1

# Separate the Schoener's D (bottom) and Hellinger I (top) parts of the data
schoener_data <- data_long[upper.tri(data, diag = TRUE), ]  # Lower triangle for Schoener's D
hellinger_data <- data_long[lower.tri(data, diag = TRUE), ]  # Upper triangle for Hellinger I

# Add a column indicating the type of data (Schoener's D or Hellinger I)
schoener_data$Type <- "Schoener's D"
hellinger_data$Type <- "Hellinger I"

# Combine both datasets into a single data frame
combined_data <- bind_rows(schoener_data, hellinger_data)

# Filter the data for Schoener's D and Hellinger I
schoener_data <- subset(combined_data, Type == "Schoener's D")
hellinger_data <- subset(combined_data, Type == "Hellinger I")

# Adjust the data to set diagonal values of 1 as a separate group
schoener_data <- schoener_data %>%
  mutate(FillColor = ifelse(Species_1 == Species_2 & Overlap == 1, NA, Overlap))  # Exclude diagonal from the gradient scale
max(schoener_data$FillColor, na.rm = T)

# Figure 3.6 This code was used to prepare the manuscript plot for Schoener's D (lower part)
Figure_3_6_schoener_plot <- ggplot(schoener_data, aes(x = Species_2, y = Species_1)) +
  geom_tile(data = filter(schoener_data, Species_1 == Species_2), fill = "black") +
  geom_tile(data = filter(schoener_data, Species_1 != Species_2), aes(fill = FillColor)) +
  scale_fill_gradient2(midpoint = 0.15, low = "lightcyan", mid = "darkcyan", high = "black", na.value = "grey",
                       limits = c(0, 0.3), oob = scales::squish) +
  geom_text(aes(label = sprintf("%.3f", Overlap)), color = "black", size = 5, fontface = "bold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 15, face = "bold.italic", color = "black", vjust = 1.5),
        axis.text.y = element_text(size = 15, face = "bold.italic", color = "black"), 
        axis.title = element_blank(),
        panel.grid = element_blank(), plot.margin = margin(10, 100, 50, 10), legend.position = "none",
        plot.title = element_blank()) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off")
Figure_3_6_schoener_plot
ggsave(".../Figures/Figure_3_6_Schoener_plot.jpg", plot = Figure_3_6_schoener_plot, width = 15, height = 10, dpi = 500)


# Adjust the data to set diagonal values of 1 as a separate group
hellinger_data <- hellinger_data %>%
  mutate(FillColor = ifelse(Species_1 == Species_2 & Overlap == 1, NA, Overlap))  # Exclude diagonal from the gradient scale
max(hellinger_data$FillColor, na.rm = T)

# Figure 3.6 This code was used to prepare the manuscript plot for Hellinger's I (upper part)
Figure_3_6_hellinger_plot <- ggplot(hellinger_data, aes(x = Species_2, y = Species_1)) +
  geom_tile(data = filter(hellinger_data, Species_1 == Species_2), fill = "black") +
  geom_tile(data = filter(hellinger_data, Species_1 != Species_2), aes(fill = FillColor)) +
  scale_fill_gradient2(midpoint = 0.25, low = "gray96", mid = "darkgray", high = "black", na.value = "grey",
                       limits = c(0, 0.5), oob = scales::squish) +
  geom_text(aes(label = sprintf("%.3f", Overlap)), color = "black", size = 5, fontface = "bold") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 15, face = "bold.italic", color = "black", vjust = 1.5),
        axis.text.y = element_text(size = 15, face = "bold.italic", color = "black"),
        axis.title = element_blank(),
        panel.grid = element_blank(), plot.margin = margin(10, 100, 50, 10), legend.position = "none",
        plot.title = element_blank()) +
  scale_x_discrete(position = "top") +
  coord_cartesian(clip = "off")
Figure_3_6_hellinger_plot
ggsave(".../Figures/Figure_3_6_Hellinger_plot.jpg", plot = Figure_3_6_hellinger_plot, width = 15, height = 10, dpi = 500)

######################################################################################################################################
# NICHE EQUIVALENCY TEST.
# Following recommendations of Warren et al. (2008), we used 100 replications for the equivalency test.
# NOTE 2: This analysis is computationally intensive and may require significant processing time. For this research, we used a
# High Prerfomenace Computing Cluster (HPCC), and took approximately 13 hours.
######################################################################################################################################

# Alligator mississippiensis  vs  other 10 species ###################################################################################
eq.test_A_mi_C_ac <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_ac, rep=100, ncores = 50)
eq.test_A_mi_C_cr <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_cr, rep=100, ncores = 50)
eq.test_A_mi_C_in <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_in, rep=100, ncores = 50)
eq.test_A_mi_C_la <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_la, rep=100, ncores = 50)
eq.test_A_mi_C_mo <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_mo, rep=100, ncores = 50)
eq.test_A_mi_C_rh <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_A_mi_C_ya <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_A_mi_M_ni <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_A_mi_P_pa <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_A_mi_P_tr <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_A_mi_list <- data.frame (eq.test_A_mi_vs_sp = c(1, eq.test_A_mi_C_ac$obs$D, eq.test_A_mi_C_cr$obs$D, eq.test_A_mi_C_in$obs$D,
                                 eq.test_A_mi_C_la$obs$D, eq.test_A_mi_C_mo$obs$D, eq.test_A_mi_C_rh$obs$D, eq.test_A_mi_C_ya$obs$D, 
                                 eq.test_A_mi_M_ni$obs$D, eq.test_A_mi_P_pa$obs$D, eq.test_A_mi_P_tr$obs$D, 1, eq.test_A_mi_C_ac$obs$I, 
                                 eq.test_A_mi_C_cr$obs$I, eq.test_A_mi_C_in$obs$I, eq.test_A_mi_C_la$obs$I, eq.test_A_mi_C_mo$obs$I, 
                                 eq.test_A_mi_C_rh$obs$I, eq.test_A_mi_C_ya$obs$I, eq.test_A_mi_M_ni$obs$I, eq.test_A_mi_P_pa$obs$I, 
                                 eq.test_A_mi_P_tr$obs$I, 1, eq.test_A_mi_C_ac$p.D, eq.test_A_mi_C_cr$p.D, eq.test_A_mi_C_in$p.D, 
                                 eq.test_A_mi_C_la$p.D, eq.test_A_mi_C_mo$p.D, eq.test_A_mi_C_rh$p.D, eq.test_A_mi_C_ya$p.D, 
                                 eq.test_A_mi_M_ni$p.D, eq.test_A_mi_P_pa$p.D, eq.test_A_mi_P_tr$p.D, 1, eq.test_A_mi_C_ac$p.I, 
                                 eq.test_A_mi_C_cr$p.I, eq.test_A_mi_C_in$p.I, eq.test_A_mi_C_la$p.I, eq.test_A_mi_C_mo$p.I, 
                                 eq.test_A_mi_C_rh$p.I, eq.test_A_mi_C_ya$p.I, eq.test_A_mi_M_ni$p.I, eq.test_A_mi_P_pa$p.I,
                                 eq.test_A_mi_P_tr$p.I))
######## D ########
trace("ecospat.plot.overlap.test", edit = TRUE) # You can edit this function to adjust each plot (e.g., bar colors as darkcyan)

ecospat.plot.overlap.test(eq.test_A_mi_C_ac, "D", ""); ecospat.plot.overlap.test(eq.test_A_mi_C_cr, "D", "") 
ecospat.plot.overlap.test(eq.test_A_mi_C_in, "D", ""); ecospat.plot.overlap.test(eq.test_A_mi_C_la, "D", "") 
ecospat.plot.overlap.test(eq.test_A_mi_C_mo, "D", "")
eq.test_A_mi_C_rh_1 <- eq.test_A_mi_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_A_mi_C_rh_1$sim <- eq.test_A_mi_C_rh_1$sim[complete.cases(eq.test_A_mi_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_A_mi_C_rh_1, "D", ""); ecospat.plot.overlap.test(eq.test_A_mi_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_A_mi_M_ni, "D", ""); ecospat.plot.overlap.test(eq.test_A_mi_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_A_mi_P_tr, "D", "")

# NOTE: Ecospat provides the p-values (previosuly saved in "eq.test_A_mi_list"). You can calculate p-values manually with this code.
# It can be used to corroborate the values saved in the object "eq.test_A_mi_list" (optional).
null_distribution_A_mi_C_ac <- eq.test_A_mi_C_ac$sim$D
observed_D_A_mi_C_ac <- eq.test_A_mi_C_ac$obs$D
p_value_manual_A_mi_C_ac <- if (observed_D_A_mi_C_ac > mean(null_distribution_A_mi_C_ac)) {
  (sum(null_distribution_A_mi_C_ac >= observed_D_A_mi_C_ac) + 1) / (length(null_distribution_A_mi_C_ac) + 1)} else {
  (sum(null_distribution_A_mi_C_ac <= observed_D_A_mi_C_ac) + 1) / (length(null_distribution_A_mi_C_ac) + 1)}
p_value_manual_A_mi_C_ac

null_distribution_A_mi_C_cr <- eq.test_A_mi_C_cr$sim$D
observed_D_A_mi_C_cr <- eq.test_A_mi_C_cr$obs$D
p_value_manual_A_mi_C_cr <- if (observed_D_A_mi_C_cr > mean(null_distribution_A_mi_C_cr)) {
  (sum(null_distribution_A_mi_C_cr >= observed_D_A_mi_C_cr) + 1) / (length(null_distribution_A_mi_C_cr) + 1)} else {
  (sum(null_distribution_A_mi_C_cr <= observed_D_A_mi_C_cr) + 1) / (length(null_distribution_A_mi_C_cr) + 1)}
p_value_manual_A_mi_C_cr

null_distribution_A_mi_C_in <- eq.test_A_mi_C_in$sim$D
observed_D_A_mi_C_in <- eq.test_A_mi_C_in$obs$D
p_value_manual_A_mi_C_in <- if (observed_D_A_mi_C_in > mean(null_distribution_A_mi_C_in)) {
  (sum(null_distribution_A_mi_C_in >= observed_D_A_mi_C_in) + 1) / (length(null_distribution_A_mi_C_in) + 1)} else {
  (sum(null_distribution_A_mi_C_in <= observed_D_A_mi_C_in) + 1) / (length(null_distribution_A_mi_C_in) + 1)}
p_value_manual_A_mi_C_in

null_distribution_A_mi_C_la <- eq.test_A_mi_C_la$sim$D
observed_D_A_mi_C_la <- eq.test_A_mi_C_la$obs$D
p_value_manual_A_mi_C_la <- if (observed_D_A_mi_C_la > mean(null_distribution_A_mi_C_la)) {
  (sum(null_distribution_A_mi_C_la >= observed_D_A_mi_C_la) + 1) / (length(null_distribution_A_mi_C_la) + 1)} else {
  (sum(null_distribution_A_mi_C_la <= observed_D_A_mi_C_la) + 1) / (length(null_distribution_A_mi_C_la) + 1)}
p_value_manual_A_mi_C_la

null_distribution_A_mi_C_mo <- eq.test_A_mi_C_mo$sim$D
observed_D_A_mi_C_mo <- eq.test_A_mi_C_mo$obs$D
p_value_manual_A_mi_C_mo <- if (observed_D_A_mi_C_mo > mean(null_distribution_A_mi_C_mo)) {
  (sum(null_distribution_A_mi_C_mo >= observed_D_A_mi_C_mo) + 1) / (length(null_distribution_A_mi_C_mo) + 1)} else {
  (sum(null_distribution_A_mi_C_mo <= observed_D_A_mi_C_mo) + 1) / (length(null_distribution_A_mi_C_mo) + 1)}
p_value_manual_A_mi_C_mo

null_distribution_A_mi_C_rh_1 <- eq.test_A_mi_C_rh_1$sim$D
observed_D_A_mi_C_rh_1 <- eq.test_A_mi_C_rh_1$obs$D
p_value_manual_A_mi_C_rh_1 <- if (observed_D_A_mi_C_rh_1 > mean(null_distribution_A_mi_C_rh_1)) {
  (sum(null_distribution_A_mi_C_rh_1 >= observed_D_A_mi_C_rh_1) + 1) / (length(null_distribution_A_mi_C_rh_1) + 1)} else {
  (sum(null_distribution_A_mi_C_rh_1 <= observed_D_A_mi_C_rh_1) + 1) / (length(null_distribution_A_mi_C_rh_1) + 1)}
p_value_manual_A_mi_C_rh_1

null_distribution_A_mi_C_ya <- eq.test_A_mi_C_ya$sim$D
observed_D_A_mi_C_ya <- eq.test_A_mi_C_ya$obs$D
p_value_manual_A_mi_C_ya <- if (observed_D_A_mi_C_ya > mean(null_distribution_A_mi_C_ya)) {
  (sum(null_distribution_A_mi_C_ya >= observed_D_A_mi_C_ya) + 1) / (length(null_distribution_A_mi_C_ya) + 1)} else {
  (sum(null_distribution_A_mi_C_ya <= observed_D_A_mi_C_ya) + 1) / (length(null_distribution_A_mi_C_ya) + 1)}
p_value_manual_A_mi_C_ya

null_distribution_A_mi_M_ni <- eq.test_A_mi_M_ni$sim$D
observed_D_A_mi_M_ni <- eq.test_A_mi_M_ni$obs$D
p_value_manual_A_mi_M_ni <- if (observed_D_A_mi_M_ni > mean(null_distribution_A_mi_M_ni)) {
  (sum(null_distribution_A_mi_M_ni >= observed_D_A_mi_M_ni) + 1) / (length(null_distribution_A_mi_M_ni) + 1)} else {
  (sum(null_distribution_A_mi_M_ni <= observed_D_A_mi_M_ni) + 1) / (length(null_distribution_A_mi_M_ni) + 1)}
p_value_manual_A_mi_M_ni

null_distribution_A_mi_P_pa <- eq.test_A_mi_P_pa$sim$D
observed_D_A_mi_P_pa <- eq.test_A_mi_P_pa$obs$D
p_value_manual_A_mi_P_pa <- if (observed_D_A_mi_P_pa > mean(null_distribution_A_mi_P_pa)) {
  (sum(null_distribution_A_mi_P_pa >= observed_D_A_mi_P_pa) + 1) / (length(null_distribution_A_mi_P_pa) + 1)} else {
  (sum(null_distribution_A_mi_P_pa <= observed_D_A_mi_P_pa) + 1) / (length(null_distribution_A_mi_P_pa) + 1)}
p_value_manual_A_mi_P_pa

null_distribution_A_mi_P_tr <- eq.test_A_mi_P_tr$sim$D
observed_D_A_mi_P_tr <- eq.test_A_mi_P_tr$obs$D
p_value_manual_A_mi_P_tr <- if (observed_D_A_mi_P_tr > mean(null_distribution_A_mi_P_tr)) {
  (sum(null_distribution_A_mi_P_tr >= observed_D_A_mi_P_tr) + 1) / (length(null_distribution_A_mi_P_tr) + 1)} else {
  (sum(null_distribution_A_mi_P_tr <= observed_D_A_mi_P_tr) + 1) / (length(null_distribution_A_mi_P_tr) + 1)}
p_value_manual_A_mi_P_tr

######## I ########
trace("ecospat.plot.overlap.test", edit = TRUE) # You can edit this function to adjust each plot (e.g., bar colors as gray)

ecospat.plot.overlap.test(eq.test_A_mi_C_ac, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_C_cr, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_C_in, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_C_la, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_C_mo, "I", "")
eq.test_A_mi_C_rh_1 <- eq.test_A_mi_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_A_mi_C_rh_1$sim <- eq.test_A_mi_C_rh_1$sim[complete.cases(eq.test_A_mi_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_A_mi_C_rh_1, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_A_mi_P_tr, "I", "")

# Calculate p values manually
null_distribution_A_mi_C_ac <- eq.test_A_mi_C_ac$sim$I
observed_D_A_mi_C_ac <- eq.test_A_mi_C_ac$obs$I
p_value_manual_A_mi_C_ac <- if (observed_D_A_mi_C_ac > mean(null_distribution_A_mi_C_ac)) {
  (sum(null_distribution_A_mi_C_ac >= observed_D_A_mi_C_ac) + 1) / (length(null_distribution_A_mi_C_ac) + 1)} else {
    (sum(null_distribution_A_mi_C_ac <= observed_D_A_mi_C_ac) + 1) / (length(null_distribution_A_mi_C_ac) + 1)}
p_value_manual_A_mi_C_ac

null_distribution_A_mi_C_cr <- eq.test_A_mi_C_cr$sim$I
observed_D_A_mi_C_cr <- eq.test_A_mi_C_cr$obs$I
p_value_manual_A_mi_C_cr <- if (observed_D_A_mi_C_cr > mean(null_distribution_A_mi_C_cr)) {
  (sum(null_distribution_A_mi_C_cr >= observed_D_A_mi_C_cr) + 1) / (length(null_distribution_A_mi_C_cr) + 1)} else {
    (sum(null_distribution_A_mi_C_cr <= observed_D_A_mi_C_cr) + 1) / (length(null_distribution_A_mi_C_cr) + 1)}
p_value_manual_A_mi_C_cr

null_distribution_A_mi_C_in <- eq.test_A_mi_C_in$sim$I
observed_D_A_mi_C_in <- eq.test_A_mi_C_in$obs$I
p_value_manual_A_mi_C_in <- if (observed_D_A_mi_C_in > mean(null_distribution_A_mi_C_in)) {
  (sum(null_distribution_A_mi_C_in >= observed_D_A_mi_C_in) + 1) / (length(null_distribution_A_mi_C_in) + 1)} else {
    (sum(null_distribution_A_mi_C_in <= observed_D_A_mi_C_in) + 1) / (length(null_distribution_A_mi_C_in) + 1)}
p_value_manual_A_mi_C_in

null_distribution_A_mi_C_la <- eq.test_A_mi_C_la$sim$I
observed_D_A_mi_C_la <- eq.test_A_mi_C_la$obs$I
p_value_manual_A_mi_C_la <- if (observed_D_A_mi_C_la > mean(null_distribution_A_mi_C_la)) {
  (sum(null_distribution_A_mi_C_la >= observed_D_A_mi_C_la) + 1) / (length(null_distribution_A_mi_C_la) + 1)} else {
    (sum(null_distribution_A_mi_C_la <= observed_D_A_mi_C_la) + 1) / (length(null_distribution_A_mi_C_la) + 1)}
p_value_manual_A_mi_C_la

null_distribution_A_mi_C_mo <- eq.test_A_mi_C_mo$sim$I
observed_D_A_mi_C_mo <- eq.test_A_mi_C_mo$obs$I
p_value_manual_A_mi_C_mo <- if (observed_D_A_mi_C_mo > mean(null_distribution_A_mi_C_mo)) {
  (sum(null_distribution_A_mi_C_mo >= observed_D_A_mi_C_mo) + 1) / (length(null_distribution_A_mi_C_mo) + 1)} else {
    (sum(null_distribution_A_mi_C_mo <= observed_D_A_mi_C_mo) + 1) / (length(null_distribution_A_mi_C_mo) + 1)}
p_value_manual_A_mi_C_mo

null_distribution_A_mi_C_rh_1 <- eq.test_A_mi_C_rh_1$sim$I
observed_D_A_mi_C_rh_1 <- eq.test_A_mi_C_rh_1$obs$I
p_value_manual_A_mi_C_rh_1 <- if (observed_D_A_mi_C_rh_1 > mean(null_distribution_A_mi_C_rh_1)) {
  (sum(null_distribution_A_mi_C_rh_1 >= observed_D_A_mi_C_rh_1) + 1) / (length(null_distribution_A_mi_C_rh_1) + 1)} else {
    (sum(null_distribution_A_mi_C_rh_1 <= observed_D_A_mi_C_rh_1) + 1) / (length(null_distribution_A_mi_C_rh_1) + 1)}
p_value_manual_A_mi_C_rh_1

null_distribution_A_mi_C_ya <- eq.test_A_mi_C_ya$sim$I
observed_D_A_mi_C_ya <- eq.test_A_mi_C_ya$obs$I
p_value_manual_A_mi_C_ya <- if (observed_D_A_mi_C_ya > mean(null_distribution_A_mi_C_ya)) {
  (sum(null_distribution_A_mi_C_ya >= observed_D_A_mi_C_ya) + 1) / (length(null_distribution_A_mi_C_ya) + 1)} else {
    (sum(null_distribution_A_mi_C_ya <= observed_D_A_mi_C_ya) + 1) / (length(null_distribution_A_mi_C_ya) + 1)}
p_value_manual_A_mi_C_ya

null_distribution_A_mi_M_ni <- eq.test_A_mi_M_ni$sim$I
observed_D_A_mi_M_ni <- eq.test_A_mi_M_ni$obs$I
p_value_manual_A_mi_M_ni <- if (observed_D_A_mi_M_ni > mean(null_distribution_A_mi_M_ni)) {
  (sum(null_distribution_A_mi_M_ni >= observed_D_A_mi_M_ni) + 1) / (length(null_distribution_A_mi_M_ni) + 1)} else {
    (sum(null_distribution_A_mi_M_ni <= observed_D_A_mi_M_ni) + 1) / (length(null_distribution_A_mi_M_ni) + 1)}
p_value_manual_A_mi_M_ni

null_distribution_A_mi_P_pa <- eq.test_A_mi_P_pa$sim$I
observed_D_A_mi_P_pa <- eq.test_A_mi_P_pa$obs$I
p_value_manual_A_mi_P_pa <- if (observed_D_A_mi_P_pa > mean(null_distribution_A_mi_P_pa)) {
  (sum(null_distribution_A_mi_P_pa >= observed_D_A_mi_P_pa) + 1) / (length(null_distribution_A_mi_P_pa) + 1)} else {
    (sum(null_distribution_A_mi_P_pa <= observed_D_A_mi_P_pa) + 1) / (length(null_distribution_A_mi_P_pa) + 1)}
p_value_manual_A_mi_P_pa

null_distribution_A_mi_P_tr <- eq.test_A_mi_P_tr$sim$I
observed_D_A_mi_P_tr <- eq.test_A_mi_P_tr$obs$I
p_value_manual_A_mi_P_tr <- if (observed_D_A_mi_P_tr > mean(null_distribution_A_mi_P_tr)) {
  (sum(null_distribution_A_mi_P_tr >= observed_D_A_mi_P_tr) + 1) / (length(null_distribution_A_mi_P_tr) + 1)} else {
    (sum(null_distribution_A_mi_P_tr <= observed_D_A_mi_P_tr) + 1) / (length(null_distribution_A_mi_P_tr) + 1)}
p_value_manual_A_mi_P_tr

# Crocodylus acutus  vs  other 9 species #############################################################################################
eq.test_C_ac_C_cr <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_cr, rep=100, ncores = 50)
eq.test_C_ac_C_in <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_in, rep=100, ncores = 50)
eq.test_C_ac_C_la <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_la, rep=100, ncores = 50)
eq.test_C_ac_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_mo, rep=100, ncores = 50)
eq.test_C_ac_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_C_ac_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_ac_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_ac_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_ac_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_ac_list <- data.frame (eq.test_C_ac_vs_sp = c(1, 1, eq.test_C_ac_C_cr$obs$D, eq.test_C_ac_C_in$obs$D, eq.test_C_ac_C_la$obs$D,
                                 eq.test_C_ac_C_mo$obs$D, eq.test_C_ac_C_rh$obs$D, eq.test_C_ac_C_ya$obs$D, eq.test_C_ac_M_ni$obs$D,
                                 eq.test_C_ac_P_pa$obs$D, eq.test_C_ac_P_tr$obs$D, 1, 1, eq.test_C_ac_C_cr$obs$I, eq.test_C_ac_C_in$obs$I, 
                                 eq.test_C_ac_C_la$obs$I, eq.test_C_ac_C_mo$obs$I, eq.test_C_ac_C_rh$obs$I, eq.test_C_ac_C_ya$obs$I, 
                                 eq.test_C_ac_M_ni$obs$I, eq.test_C_ac_P_pa$obs$I, eq.test_C_ac_P_tr$obs$I, 1, 1, eq.test_C_ac_C_cr$p.D, 
                                 eq.test_C_ac_C_in$p.D, eq.test_C_ac_C_la$p.D, eq.test_C_ac_C_mo$p.D, eq.test_C_ac_C_rh$p.D, 
                                 eq.test_C_ac_C_ya$p.D, eq.test_C_ac_M_ni$p.D, eq.test_C_ac_P_pa$p.D, eq.test_C_ac_P_tr$p.D, 1, 1, 
                                 eq.test_C_ac_C_cr$p.I, eq.test_C_ac_C_in$p.I, eq.test_C_ac_C_la$p.I, eq.test_C_ac_C_mo$p.I, 
                                 eq.test_C_ac_C_rh$p.I, eq.test_C_ac_C_ya$p.I, eq.test_C_ac_M_ni$p.I, eq.test_C_ac_P_pa$p.I, 
                                 eq.test_C_ac_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_ac_C_cr, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_in, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_la, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_mo, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_rh, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_ac_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_ac_C_cr, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_in, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_la, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_mo, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_rh, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_C_ac_P_tr, "I", "")

# Caiman crocodilus  vs  other 8 species #############################################################################################
eq.test_C_cr_C_in <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_in, rep=100, ncores = 50)
eq.test_C_cr_C_la <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_la, rep=100, ncores = 50)
eq.test_C_cr_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_mo, rep=100, ncores = 50)
eq.test_C_cr_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_C_cr_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_cr_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_cr_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_cr_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_cr_list <- data.frame (eq.test_C_cr_vs_sp = c(1, 1, 1, eq.test_C_cr_C_in$obs$D, eq.test_C_cr_C_la$obs$D, eq.test_C_cr_C_mo$obs$D, 
                                 eq.test_C_cr_C_rh$obs$D, eq.test_C_cr_C_ya$obs$D, eq.test_C_cr_M_ni$obs$D, eq.test_C_cr_P_pa$obs$D, 
                                 eq.test_C_cr_P_tr$obs$D, 1, 1, 1, eq.test_C_cr_C_in$obs$I, eq.test_C_cr_C_la$obs$I, eq.test_C_cr_C_mo$obs$I, 
                                 eq.test_C_cr_C_rh$obs$I, eq.test_C_cr_C_ya$obs$I, eq.test_C_cr_M_ni$obs$I, eq.test_C_cr_P_pa$obs$I, 
                                 eq.test_C_cr_P_tr$obs$I, 1, 1, 1, eq.test_C_cr_C_in$p.D, eq.test_C_cr_C_la$p.D, eq.test_C_cr_C_mo$p.D, 
                                 eq.test_C_cr_C_rh$p.D, eq.test_C_cr_C_ya$p.D, eq.test_C_cr_M_ni$p.D, eq.test_C_cr_P_pa$p.D, 
                                 eq.test_C_cr_P_tr$p.D, 1, 1, 1, eq.test_C_cr_C_in$p.I, eq.test_C_cr_C_la$p.I, eq.test_C_cr_C_mo$p.I, 
                                 eq.test_C_cr_C_rh$p.I, eq.test_C_cr_C_ya$p.I, eq.test_C_cr_M_ni$p.I, eq.test_C_cr_P_pa$p.I, 
                                 eq.test_C_cr_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_cr_C_in, "D", "") 
ecospat.plot.overlap.test(eq.test_C_cr_C_la, "D", "") 
ecospat.plot.overlap.test(eq.test_C_cr_C_mo, "D", "")
eq.test_C_cr_C_rh_1 <- eq.test_C_cr_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_cr_C_rh_1$sim <- eq.test_C_cr_C_rh_1$sim[complete.cases(eq.test_C_cr_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_cr_C_rh_1, "D", "")
ecospat.plot.overlap.test(eq.test_C_cr_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_cr_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_cr_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_cr_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_cr_C_in, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_C_la, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_C_mo, "I", "")
eq.test_C_cr_C_rh_1 <- eq.test_C_cr_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_cr_C_rh_1$sim <- eq.test_C_cr_C_rh_1$sim[complete.cases(eq.test_C_cr_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_cr_C_rh_1, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_C_cr_P_tr, "I", "")

# Crocodylus intermedius  vs  other 7 species ########################################################################################
eq.test_C_in_C_la <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_la, rep=100, ncores = 50)
eq.test_C_in_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_mo, rep=100, ncores = 50)
eq.test_C_in_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_C_in_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_in_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_in_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_in_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_in_list <- data.frame (eq.test_C_in_vs_sp = c(1, 1, 1, 1, eq.test_C_in_C_la$obs$D, eq.test_C_in_C_mo$obs$D, eq.test_C_in_C_rh$obs$D,
                                 eq.test_C_in_C_ya$obs$D, eq.test_C_in_M_ni$obs$D, eq.test_C_in_P_pa$obs$D, eq.test_C_in_P_tr$obs$D,
                                 1, 1, 1, 1, eq.test_C_in_C_la$obs$I, eq.test_C_in_C_mo$obs$I, eq.test_C_in_C_rh$obs$I, eq.test_C_in_C_ya$obs$I, 
                                 eq.test_C_in_M_ni$obs$I, eq.test_C_in_P_pa$obs$I, eq.test_C_in_P_tr$obs$I, 1, 1, 1, 1, eq.test_C_in_C_la$p.D,
                                 eq.test_C_in_C_mo$p.D, eq.test_C_in_C_rh$p.D, eq.test_C_in_C_ya$p.D, eq.test_C_in_M_ni$p.D, eq.test_C_in_P_pa$p.D, 
                                 eq.test_C_in_P_tr$p.D, 1, 1, 1, 1, eq.test_C_in_C_la$p.I, eq.test_C_in_C_mo$p.I, eq.test_C_in_C_rh$p.I, 
                                 eq.test_C_in_C_ya$p.I, eq.test_C_in_M_ni$p.I, eq.test_C_in_P_pa$p.I, eq.test_C_in_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_in_C_la, "D", "") 
ecospat.plot.overlap.test(eq.test_C_in_C_mo, "D", "")
eq.test_C_in_C_rh_1 <- eq.test_C_in_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_in_C_rh_1$sim <- eq.test_C_in_C_rh_1$sim[complete.cases(eq.test_C_in_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_in_C_rh_1, "D", "")
ecospat.plot.overlap.test(eq.test_C_in_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_in_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_in_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_in_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_in_C_la, "I", "")
ecospat.plot.overlap.test(eq.test_C_in_C_mo, "I", "")
eq.test_C_in_C_rh_1 <- eq.test_C_in_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_in_C_rh_1$sim <- eq.test_C_in_C_rh_1$sim[complete.cases(eq.test_C_in_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_in_C_rh_1, "I", "")
ecospat.plot.overlap.test(eq.test_C_in_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_C_in_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_C_in_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_C_in_P_tr, "I", "")

# Caiman latirostris  vs  other 6 species ############################################################################################
eq.test_C_la_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_mo, rep=100, ncores = 50)
eq.test_C_la_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_C_la_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_la_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_la_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_la_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_la_list <- data.frame (eq.test_C_la_vs_sp = c(1, 1, 1, 1, 1, eq.test_C_la_C_mo$obs$D, eq.test_C_la_C_rh$obs$D, eq.test_C_la_C_ya$obs$D, 
                                 eq.test_C_la_M_ni$obs$D, eq.test_C_la_P_pa$obs$D, eq.test_C_la_P_tr$obs$D, 1, 1, 1, 1, 1, eq.test_C_la_C_mo$obs$I, 
                                 eq.test_C_la_C_rh$obs$I, eq.test_C_la_C_ya$obs$I, eq.test_C_la_M_ni$obs$I, eq.test_C_la_P_pa$obs$I, 
                                 eq.test_C_la_P_tr$obs$I, 1, 1, 1, 1, 1, eq.test_C_la_C_mo$p.D, eq.test_C_la_C_rh$p.D, eq.test_C_la_C_ya$p.D, 
                                 eq.test_C_la_M_ni$p.D, eq.test_C_la_P_pa$p.D, eq.test_C_la_P_tr$p.D, 1, 1, 1, 1, 1, eq.test_C_la_C_mo$p.I, 
                                 eq.test_C_la_C_rh$p.I, eq.test_C_la_C_ya$p.I, eq.test_C_la_M_ni$p.I, eq.test_C_la_P_pa$p.I, eq.test_C_la_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_la_C_mo, "D", "")
eq.test_C_la_C_rh_1 <- eq.test_C_la_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_la_C_rh_1$sim <- eq.test_C_la_C_rh_1$sim[complete.cases(eq.test_C_la_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_la_C_rh_1, "D", "")
ecospat.plot.overlap.test(eq.test_C_la_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_la_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_la_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_la_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_la_C_mo, "I", "")
eq.test_C_la_C_rh_1 <- eq.test_C_la_C_rh # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_la_C_rh_1$sim <- eq.test_C_la_C_rh_1$sim[complete.cases(eq.test_C_la_C_rh_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_la_C_rh_1, "I", "")
ecospat.plot.overlap.test(eq.test_C_la_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_C_la_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_C_la_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_C_la_P_tr, "I", "")

# Crocodylus moreletii  vs  other 5 species ##########################################################################################
eq.test_C_mo_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_C_rh, rep=100, ncores = 50)
eq.test_C_mo_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_mo_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_mo_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_mo_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_mo_list <- data.frame (eq.test_C_mo_vs_sp = c(1, 1, 1, 1, 1, 1, eq.test_C_mo_C_rh$obs$D, eq.test_C_mo_C_ya$obs$D, eq.test_C_mo_M_ni$obs$D,
                                 eq.test_C_mo_P_pa$obs$D, eq.test_C_mo_P_tr$obs$D, 1, 1, 1, 1, 1, 1, eq.test_C_mo_C_rh$obs$I, eq.test_C_mo_C_ya$obs$I, 
                                 eq.test_C_mo_M_ni$obs$I, eq.test_C_mo_P_pa$obs$I, eq.test_C_mo_P_tr$obs$I, 1, 1, 1, 1, 1, 1, eq.test_C_mo_C_rh$p.D, 
                                 eq.test_C_mo_C_ya$p.D, eq.test_C_mo_M_ni$p.D, eq.test_C_mo_P_pa$p.D, eq.test_C_mo_P_tr$p.D, 1, 1, 1, 1, 1, 1, 
                                 eq.test_C_mo_C_rh$p.I, eq.test_C_mo_C_ya$p.I, eq.test_C_mo_M_ni$p.I, eq.test_C_mo_P_pa$p.I, eq.test_C_mo_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_mo_C_rh, "D", "")
ecospat.plot.overlap.test(eq.test_C_mo_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_mo_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_mo_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_mo_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_mo_C_rh, "I", "")
ecospat.plot.overlap.test(eq.test_C_mo_C_ya, "I", "")
ecospat.plot.overlap.test(eq.test_C_mo_M_ni, "I", "")
ecospat.plot.overlap.test(eq.test_C_mo_P_pa, "I", "")
ecospat.plot.overlap.test(eq.test_C_mo_P_tr, "I", "")

# Crocodylus rhombifer  vs  other 4 species ##########################################################################################
eq.test_C_rh_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_C_ya, rep=100, ncores = 50)
eq.test_C_rh_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_rh_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_rh_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_rh_list <- data.frame (eq.test_C_rh_vs_sp = c(1, 1, 1, 1, 1, 1, 1, eq.test_C_rh_C_ya$obs$D, eq.test_C_rh_M_ni$obs$D, 
                                 eq.test_C_rh_P_pa$obs$D, eq.test_C_rh_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, eq.test_C_rh_C_ya$obs$I, 
                                 eq.test_C_rh_M_ni$obs$I, eq.test_C_rh_P_pa$obs$I, eq.test_C_rh_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 
                                 eq.test_C_rh_C_ya$p.D, eq.test_C_rh_M_ni$p.D, eq.test_C_rh_P_pa$p.D, eq.test_C_rh_P_tr$p.D, 1, 1, 1, 1, 
                                 1, 1, 1, eq.test_C_rh_C_ya$p.I, eq.test_C_rh_M_ni$p.I, eq.test_C_rh_P_pa$p.I, eq.test_C_rh_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_rh_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_rh_M_ni, "D", "")
eq.test_C_rh_P_pa_1 <- eq.test_C_rh_P_pa # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_rh_P_pa_1$sim <- eq.test_C_rh_P_pa_1$sim[complete.cases(eq.test_C_rh_P_pa_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_rh_P_pa_1, "D", "")
eq.test_C_rh_P_tr_1 <- eq.test_C_rh_P_tr # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_rh_P_tr_1$sim <- eq.test_C_rh_P_tr_1$sim[complete.cases(eq.test_C_rh_P_tr_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_rh_P_tr_1, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_rh_C_ya, "D", "")
ecospat.plot.overlap.test(eq.test_C_rh_M_ni, "D", "")
eq.test_C_rh_P_pa_1 <- eq.test_C_rh_P_pa # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_rh_P_pa_1$sim <- eq.test_C_rh_P_pa_1$sim[complete.cases(eq.test_C_rh_P_pa_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_rh_P_pa_1, "D", "")
eq.test_C_rh_P_tr_1 <- eq.test_C_rh_P_tr # Since there are NA values, we need to only use complete cases. So we create a new object
eq.test_C_rh_P_tr_1$sim <- eq.test_C_rh_P_tr_1$sim[complete.cases(eq.test_C_rh_P_tr_1$sim), ]
ecospat.plot.overlap.test(eq.test_C_rh_P_tr_1, "D", "")

# Caiman yacare  vs  other 3 species #################################################################################################
eq.test_C_ya_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_M_ni, rep=100, ncores = 50)
eq.test_C_ya_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_C_ya_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_C_ya_list <- data.frame (eq.test_C_ya_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, eq.test_C_ya_M_ni$obs$D, eq.test_C_ya_P_pa$obs$D, 
                                 eq.test_C_ya_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_C_ya_M_ni$obs$I, eq.test_C_ya_P_pa$obs$I, 
                                 eq.test_C_ya_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_C_ya_M_ni$p.D, eq.test_C_ya_P_pa$p.D, 
                                 eq.test_C_ya_P_tr$p.D, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_C_ya_M_ni$p.I, eq.test_C_ya_P_pa$p.I, 
                                 eq.test_C_ya_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_C_ya_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_ya_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_ya_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_C_ya_M_ni, "D", "")
ecospat.plot.overlap.test(eq.test_C_ya_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_C_ya_P_tr, "D", "")

# Melanosuchus niger  vs  other 2 species ############################################################################################
eq.test_M_ni_P_pa <- ecospat.niche.equivalency.test (grid_clim_M_ni, grid_clim_P_pa, rep=100, ncores = 50)
eq.test_M_ni_P_tr <- ecospat.niche.equivalency.test (grid_clim_M_ni, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_M_ni_list <- data.frame (eq.test_M_ni_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_M_ni_P_pa$obs$D, eq.test_M_ni_P_tr$obs$D,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_M_ni_P_pa$obs$I, eq.test_M_ni_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 eq.test_M_ni_P_pa$p.D, eq.test_M_ni_P_tr$p.D, 1, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_M_ni_P_pa$p.I, 
                                 eq.test_M_ni_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_M_ni_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_M_ni_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_M_ni_P_pa, "D", "")
ecospat.plot.overlap.test(eq.test_M_ni_P_tr, "D", "")

# Paleosuchus palpebrosus  vs  other 1 species #######################################################################################
eq.test_P_pa_P_tr <- ecospat.niche.equivalency.test (grid_clim_P_pa, grid_clim_P_tr, rep=100, ncores = 50)

eq.test_P_pa_list <- data.frame (eq.test_P_pa_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_P_pa_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, 1, 
                                 1, 1, eq.test_P_pa_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, eq.test_P_pa_P_tr$p.D, 1, 1, 1, 1, 1, 1, 
                                 1, 1, 1, 1, eq.test_P_pa_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(eq.test_P_pa_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(eq.test_P_pa_P_tr, "D", "")

Niche_overlap_equivalency_test <- data.frame(eq.test_A_mi_list, eq.test_C_ac_list, eq.test_C_cr_list, eq.test_C_in_list, eq.test_C_la_list,
                                             eq.test_C_mo_list, eq.test_C_rh_list, eq.test_C_ya_list, eq.test_M_ni_list, eq.test_P_pa_list)
write.csv(Niche_overlap_equivalency_test, ".../Tables/Niche_overlap_equivalency_test.csv")

######################################################################################################################################
# NICHE SIMILARITY TEST.
# Following recommendations of Warren et al. (2008), we used 100 replications for the similarity test.
# Similarity test varies depending of the position of the species compared within the code, thus we will do 2 similarity analyses:
# Similarity 1: a -> b
# Similarity 2: b -> a
# NOTE 3: This analysis is computationally intensive and may require significant processing time. For this research, we used a
# High Prerfomenace Computing Cluster (HPCC), and took approximately 8 hours.
######################################################################################################################################

######################################################################################################################################
# Similarity 1: a -> b
######################################################################################################################################

# Alligator mississippiensis  vs  other 10 species ###################################################################################
sim.test_A_mi_C_ac <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_ac, rep=100, ncores = 50)
sim.test_A_mi_C_cr <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_cr, rep=100, ncores = 50)
sim.test_A_mi_C_in <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_in, rep=100, ncores = 50)
sim.test_A_mi_C_la <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_la, rep=100, ncores = 50)
sim.test_A_mi_C_mo <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_mo, rep=100, ncores = 50)
sim.test_A_mi_C_rh <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_A_mi_C_ya <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_A_mi_M_ni <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_A_mi_P_pa <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_A_mi_P_tr <- ecospat.niche.equivalency.test (grid_clim_A_mi, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_A_mi_list <- data.frame (sim.test_A_mi_vs_sp = c(1, sim.test_A_mi_C_ac$obs$D, sim.test_A_mi_C_cr$obs$D, sim.test_A_mi_C_in$obs$D, 
                                  sim.test_A_mi_C_la$obs$D, sim.test_A_mi_C_mo$obs$D, sim.test_A_mi_C_rh$obs$D, sim.test_A_mi_C_ya$obs$D, 
                                  sim.test_A_mi_M_ni$obs$D, sim.test_A_mi_P_pa$obs$D, sim.test_A_mi_P_tr$obs$D, 1, sim.test_A_mi_C_ac$obs$I, 
                                  sim.test_A_mi_C_cr$obs$I, sim.test_A_mi_C_in$obs$I, sim.test_A_mi_C_la$obs$I, sim.test_A_mi_C_mo$obs$I, 
                                  sim.test_A_mi_C_rh$obs$I, sim.test_A_mi_C_ya$obs$I, sim.test_A_mi_M_ni$obs$I, sim.test_A_mi_P_pa$obs$I, 
                                  sim.test_A_mi_P_tr$obs$I, 1, sim.test_A_mi_C_ac$p.D, sim.test_A_mi_C_cr$p.D, sim.test_A_mi_C_in$p.D, 
                                  sim.test_A_mi_C_la$p.D, sim.test_A_mi_C_mo$p.D, sim.test_A_mi_C_rh$p.D, sim.test_A_mi_C_ya$p.D, 
                                  sim.test_A_mi_M_ni$p.D, sim.test_A_mi_P_pa$p.D, sim.test_A_mi_P_tr$p.D, 1, sim.test_A_mi_C_ac$p.I, 
                                  sim.test_A_mi_C_cr$p.I, sim.test_A_mi_C_in$p.I, sim.test_A_mi_C_la$p.I, sim.test_A_mi_C_mo$p.I, 
                                  sim.test_A_mi_C_rh$p.I, sim.test_A_mi_C_ya$p.I, sim.test_A_mi_M_ni$p.I, sim.test_A_mi_P_pa$p.I, 
                                  sim.test_A_mi_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_A_mi_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_cr, "D", "") 
ecospat.plot.overlap.test(sim.test_A_mi_C_in, "D", "") 
ecospat.plot.overlap.test(sim.test_A_mi_C_la, "D", "") 
ecospat.plot.overlap.test(sim.test_A_mi_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_A_mi_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_A_mi_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_A_mi_P_tr, "I", "")

# Crocodylus acutus  vs  other 9 species #############################################################################################
sim.test_C_ac_C_cr <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_cr, rep=100, ncores = 50)
sim.test_C_ac_C_in <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_in, rep=100, ncores = 50)
sim.test_C_ac_C_la <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_la, rep=100, ncores = 50)
sim.test_C_ac_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_mo, rep=100, ncores = 50)
sim.test_C_ac_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_C_ac_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_ac_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_ac_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_ac_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_ac, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_ac_list <- data.frame (sim.test_C_ac_vs_sp = c(1, 1, sim.test_C_ac_C_cr$obs$D, sim.test_C_ac_C_in$obs$D, sim.test_C_ac_C_la$obs$D,
                                  sim.test_C_ac_C_mo$obs$D, sim.test_C_ac_C_rh$obs$D, sim.test_C_ac_C_ya$obs$D, sim.test_C_ac_M_ni$obs$D,
                                  sim.test_C_ac_P_pa$obs$D, sim.test_C_ac_P_tr$obs$D, 1, 1, sim.test_C_ac_C_cr$obs$I, sim.test_C_ac_C_in$obs$I, 
                                  sim.test_C_ac_C_la$obs$I, sim.test_C_ac_C_mo$obs$I, sim.test_C_ac_C_rh$obs$I, sim.test_C_ac_C_ya$obs$I, 
                                  sim.test_C_ac_M_ni$obs$I, sim.test_C_ac_P_pa$obs$I, sim.test_C_ac_P_tr$obs$I, 1, 1, sim.test_C_ac_C_cr$p.D, 
                                  sim.test_C_ac_C_in$p.D, sim.test_C_ac_C_la$p.D, sim.test_C_ac_C_mo$p.D, sim.test_C_ac_C_rh$p.D, 
                                  sim.test_C_ac_C_ya$p.D, sim.test_C_ac_M_ni$p.D, sim.test_C_ac_P_pa$p.D, sim.test_C_ac_P_tr$p.D, 1, 1, 
                                  sim.test_C_ac_C_cr$p.I, sim.test_C_ac_C_in$p.I, sim.test_C_ac_C_la$p.I, sim.test_C_ac_C_mo$p.I, 
                                  sim.test_C_ac_C_rh$p.I, sim.test_C_ac_C_ya$p.I, sim.test_C_ac_M_ni$p.I, sim.test_C_ac_P_pa$p.I, 
                                  sim.test_C_ac_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_ac_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_ac_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_ac_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_C_ac_P_tr, "I", "")

# Caiman crocodilus  vs  other 8 species #############################################################################################
sim.test_C_cr_C_in <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_in, rep=100, ncores = 50)
sim.test_C_cr_C_la <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_la, rep=100, ncores = 50)
sim.test_C_cr_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_mo, rep=100, ncores = 50)
sim.test_C_cr_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_C_cr_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_cr_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_cr_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_cr_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_cr, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_cr_list <- data.frame (sim.test_C_cr_vs_sp = c(1, 1, 1, sim.test_C_cr_C_in$obs$D, sim.test_C_cr_C_la$obs$D, 
                                  sim.test_C_cr_C_mo$obs$D, sim.test_C_cr_C_rh$obs$D, sim.test_C_cr_C_ya$obs$D, sim.test_C_cr_M_ni$obs$D,
                                  sim.test_C_cr_P_pa$obs$D, sim.test_C_cr_P_tr$obs$D, 1, 1, 1, sim.test_C_cr_C_in$obs$I, sim.test_C_cr_C_la$obs$I,
                                  sim.test_C_cr_C_mo$obs$I, sim.test_C_cr_C_rh$obs$I, sim.test_C_cr_C_ya$obs$I,  sim.test_C_cr_M_ni$obs$I,
                                  sim.test_C_cr_P_pa$obs$I, sim.test_C_cr_P_tr$obs$I, 1, 1, 1, sim.test_C_cr_C_in$p.D, sim.test_C_cr_C_la$p.D,
                                  sim.test_C_cr_C_mo$p.D, sim.test_C_cr_C_rh$p.D, sim.test_C_cr_C_ya$p.D, sim.test_C_cr_M_ni$p.D,
                                  sim.test_C_cr_P_pa$p.D, sim.test_C_cr_P_tr$p.D, 1, 1, 1, sim.test_C_cr_C_in$p.I, sim.test_C_cr_C_la$p.I,
                                  sim.test_C_cr_C_mo$p.I, sim.test_C_cr_C_rh$p.I, sim.test_C_cr_C_ya$p.I, sim.test_C_cr_M_ni$p.I,
                                  sim.test_C_cr_P_pa$p.I, sim.test_C_cr_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_cr_C_in, "D", "") 
ecospat.plot.overlap.test(sim.test_C_cr_C_la, "D", "") 
ecospat.plot.overlap.test(sim.test_C_cr_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_cr_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_P_tr, "I", "")

# Crocodylus intermedius  vs  other 7 species ########################################################################################
sim.test_C_in_C_la <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_la, rep=100, ncores = 50)
sim.test_C_in_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_mo, rep=100, ncores = 50)
sim.test_C_in_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_C_in_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_in_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_in_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_in_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_in, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_in_list <- data.frame (sim.test_C_in_vs_sp = c(1, 1, 1, 1, sim.test_C_in_C_la$obs$D, sim.test_C_in_C_mo$obs$D, 
                                  sim.test_C_in_C_rh$obs$D, sim.test_C_in_C_ya$obs$D, sim.test_C_in_M_ni$obs$D, sim.test_C_in_P_pa$obs$D, 
                                  sim.test_C_in_P_tr$obs$D, 1, 1, 1, 1, sim.test_C_in_C_la$obs$I, sim.test_C_in_C_mo$obs$I, 
                                  sim.test_C_in_C_rh$obs$I, sim.test_C_in_C_ya$obs$I, sim.test_C_in_M_ni$obs$I, sim.test_C_in_P_pa$obs$I, 
                                  sim.test_C_in_P_tr$obs$I, 1, 1, 1, 1, sim.test_C_in_C_la$p.D, sim.test_C_in_C_mo$p.D, sim.test_C_in_C_rh$p.D, 
                                  sim.test_C_in_C_ya$p.D, sim.test_C_in_M_ni$p.D, sim.test_C_in_P_pa$p.D, sim.test_C_in_P_tr$p.D,
                                  1, 1, 1, 1, sim.test_C_in_C_la$p.I, sim.test_C_in_C_mo$p.I, sim.test_C_in_C_rh$p.I, sim.test_C_in_C_ya$p.I, 
                                  sim.test_C_in_M_ni$p.I, sim.test_C_in_P_pa$p.I, sim.test_C_in_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_in_C_la, "D", "") 
ecospat.plot.overlap.test(sim.test_C_in_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_C_in_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_C_in_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_in_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_in_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_in_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_in_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_P_tr, "I", "")

# Caiman latirostris  vs  other 6 species ############################################################################################
sim.test_C_la_C_mo <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_mo, rep=100, ncores = 50)
sim.test_C_la_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_C_la_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_la_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_la_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_la_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_la, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_la_list <- data.frame (sim.test_C_la_vs_sp = c(1, 1, 1, 1, 1, sim.test_C_la_C_mo$obs$D, sim.test_C_la_C_rh$obs$D, 
                                  sim.test_C_la_C_ya$obs$D, sim.test_C_la_M_ni$obs$D, sim.test_C_la_P_pa$obs$D, sim.test_C_la_P_tr$obs$D,
                                  1, 1, 1, 1, 1, sim.test_C_la_C_mo$obs$I, sim.test_C_la_C_rh$obs$I, sim.test_C_la_C_ya$obs$I, 
                                  sim.test_C_la_M_ni$obs$I, sim.test_C_la_P_pa$obs$I, sim.test_C_la_P_tr$obs$I, 1, 1, 1, 1, 1,
                                  sim.test_C_la_C_mo$p.D, sim.test_C_la_C_rh$p.D, sim.test_C_la_C_ya$p.D, sim.test_C_la_M_ni$p.D,
                                  sim.test_C_la_P_pa$p.D, sim.test_C_la_P_tr$p.D, 1, 1, 1, 1, 1, sim.test_C_la_C_mo$p.I, sim.test_C_la_C_rh$p.I, 
                                  sim.test_C_la_C_ya$p.I, sim.test_C_la_M_ni$p.I, sim.test_C_la_P_pa$p.I, sim.test_C_la_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_la_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_C_la_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_C_la_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_la_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_la_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_la_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_la_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_P_tr, "I", "")

# Crocodylus moreletii  vs  other 5 species ##########################################################################################
sim.test_C_mo_C_rh <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_C_rh, rep=100, ncores = 50)
sim.test_C_mo_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_mo_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_mo_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_mo_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_mo, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_mo_list <- data.frame (sim.test_C_mo_vs_sp = c(1, 1, 1, 1, 1, 1, sim.test_C_mo_C_rh$obs$D, sim.test_C_mo_C_ya$obs$D, 
                                  sim.test_C_mo_M_ni$obs$D, sim.test_C_mo_P_pa$obs$D, sim.test_C_mo_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 
                                  sim.test_C_mo_C_rh$obs$I, sim.test_C_mo_C_ya$obs$I, sim.test_C_mo_M_ni$obs$I, sim.test_C_mo_P_pa$obs$I, 
                                  sim.test_C_mo_P_tr$obs$I, 1, 1, 1, 1, 1, 1, sim.test_C_mo_C_rh$p.D, sim.test_C_mo_C_ya$p.D, 
                                  sim.test_C_mo_M_ni$p.D, sim.test_C_mo_P_pa$p.D, sim.test_C_mo_P_tr$p.D, 1, 1, 1, 1, 1, 1, 
                                  sim.test_C_mo_C_rh$p.I, sim.test_C_mo_C_ya$p.I, sim.test_C_mo_M_ni$p.I, sim.test_C_mo_P_pa$p.I, 
                                  sim.test_C_mo_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_mo_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_C_mo_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_mo_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_mo_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_mo_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_mo_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_P_pa, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_P_tr, "I", "")

# Crocodylus rhombifer  vs  other 4 species ##########################################################################################
sim.test_C_rh_C_ya <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_C_ya, rep=100, ncores = 50)
sim.test_C_rh_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_rh_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_rh_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_rh, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_rh_list <- data.frame (sim.test_C_rh_vs_sp = c(1, 1, 1, 1, 1, 1, 1, sim.test_C_rh_C_ya$obs$D, sim.test_C_rh_M_ni$obs$D,
                                  sim.test_C_rh_P_pa$obs$D, sim.test_C_rh_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, sim.test_C_rh_C_ya$obs$I, 
                                  sim.test_C_rh_M_ni$obs$I, sim.test_C_rh_P_pa$obs$I, sim.test_C_rh_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 
                                  sim.test_C_rh_C_ya$p.D, sim.test_C_rh_M_ni$p.D, sim.test_C_rh_P_pa$p.D, sim.test_C_rh_P_tr$p.D,
                                  1, 1, 1, 1, 1, 1, 1, sim.test_C_rh_C_ya$p.I, sim.test_C_rh_M_ni$p.I, sim.test_C_rh_P_pa$p.I, 
                                  sim.test_C_rh_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_rh_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_rh_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_P_tr, "D", "")

# Caiman yacare  vs  other 3 species #################################################################################################
sim.test_C_ya_M_ni <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_M_ni, rep=100, ncores = 50)
sim.test_C_ya_P_pa <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_C_ya_P_tr <- ecospat.niche.equivalency.test (grid_clim_C_ya, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_C_ya_list <- data.frame (sim.test_C_ya_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_M_ni$obs$D, sim.test_C_ya_P_pa$obs$D, 
                                  sim.test_C_ya_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_M_ni$obs$I, sim.test_C_ya_P_pa$obs$I, 
                                  sim.test_C_ya_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_M_ni$p.D, sim.test_C_ya_P_pa$p.D, 
                                  sim.test_C_ya_P_tr$p.D, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_M_ni$p.I, sim.test_C_ya_P_pa$p.I, 
                                  sim.test_C_ya_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_ya_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_ya_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_P_tr, "D", "")

# Melanosuchus niger  vs  other 2 species ############################################################################################
sim.test_M_ni_P_pa <- ecospat.niche.equivalency.test (grid_clim_M_ni, grid_clim_P_pa, rep=100, ncores = 50)
sim.test_M_ni_P_tr <- ecospat.niche.equivalency.test (grid_clim_M_ni, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_M_ni_list <- data.frame (sim.test_M_ni_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_P_pa$obs$D, sim.test_M_ni_P_tr$obs$D,
                                  1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_P_pa$obs$I, sim.test_M_ni_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  sim.test_M_ni_P_pa$p.D, sim.test_M_ni_P_tr$p.D, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_P_pa$p.I, 
                                  sim.test_M_ni_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_M_ni_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_M_ni_P_pa, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_P_tr, "D", "")

# Paleosuchus palpebrosus  vs  other 1 species #######################################################################################
sim.test_P_pa_P_tr <- ecospat.niche.equivalency.test (grid_clim_P_pa, grid_clim_P_tr, rep=100, ncores = 50)

sim.test_P_pa_list <- data.frame (sim.test_P_pa_vs_sp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_pa_P_tr$obs$D, 1, 1, 1, 1, 1, 1, 1, 
                                  1, 1, 1, sim.test_P_pa_P_tr$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_pa_P_tr$p.D, 1, 1, 1, 1, 
                                  1, 1, 1, 1, 1, 1, sim.test_P_pa_P_tr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_P_pa_P_tr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_P_pa_P_tr, "D", "")

Niche_overlap_similarity_test_1 <- data.frame(sim.test_A_mi_list, sim.test_C_ac_list, sim.test_C_cr_list, sim.test_C_in_list, sim.test_C_la_list,
                                            sim.test_C_mo_list, sim.test_C_rh_list, sim.test_C_ya_list, sim.test_M_ni_list, sim.test_P_pa_list)
write.csv(Niche_overlap_similarity_test_1, ".../Niche_overlap_similarity_test_1.csv")

######################################################################################################################################
# Similarity 2: b -> a
######################################################################################################################################

# Other 10 species vs Alligator mississippiensis #####################################################################################
sim.test_C_ac_A_mi <- ecospat.niche.similarity.test (grid_clim_C_ac, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_cr_A_mi <- ecospat.niche.similarity.test (grid_clim_C_cr, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_in_A_mi <- ecospat.niche.similarity.test (grid_clim_C_in, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_la_A_mi <- ecospat.niche.similarity.test (grid_clim_C_la, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_mo_A_mi <- ecospat.niche.similarity.test (grid_clim_C_mo, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_rh_A_mi <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_C_ya_A_mi <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_M_ni_A_mi <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_P_pa_A_mi <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_A_mi, rep=100,  ncores = 50)
sim.test_P_tr_A_mi <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_A_mi, rep=100,  ncores = 50)

sim.test_A_mi_2_list <- data.frame (sim.test_sp_vs_A_mi = c(1, sim.test_C_ac_A_mi$obs$D, sim.test_C_cr_A_mi$obs$D, sim.test_C_in_A_mi$obs$D, 
                                    sim.test_C_la_A_mi$obs$D, sim.test_C_mo_A_mi$obs$D, sim.test_C_rh_A_mi$obs$D, sim.test_C_ya_A_mi$obs$D, 
                                    sim.test_M_ni_A_mi$obs$D, sim.test_P_pa_A_mi$obs$D, sim.test_P_tr_A_mi$obs$D, 1, sim.test_C_ac_A_mi$obs$I, 
                                    sim.test_C_cr_A_mi$obs$I, sim.test_C_in_A_mi$obs$I, sim.test_C_la_A_mi$obs$I, sim.test_C_mo_A_mi$obs$I, 
                                    sim.test_C_rh_A_mi$obs$I, sim.test_C_ya_A_mi$obs$I, sim.test_M_ni_A_mi$obs$I, sim.test_P_pa_A_mi$obs$I, 
                                    sim.test_P_tr_A_mi$obs$I, 1, sim.test_C_ac_A_mi$p.D, sim.test_C_cr_A_mi$p.D, sim.test_C_in_A_mi$p.D, 
                                    sim.test_C_la_A_mi$p.D, sim.test_C_mo_A_mi$p.D, sim.test_C_rh_A_mi$p.D, sim.test_C_ya_A_mi$p.D, 
                                    sim.test_M_ni_A_mi$p.D, sim.test_P_pa_A_mi$p.D, sim.test_P_tr_A_mi$p.D, 1, sim.test_C_ac_A_mi$p.I, 
                                    sim.test_C_cr_A_mi$p.I, sim.test_C_in_A_mi$p.I, sim.test_C_la_A_mi$p.I, sim.test_C_mo_A_mi$p.I, 
                                    sim.test_C_rh_A_mi$p.I, sim.test_C_ya_A_mi$p.I, sim.test_M_ni_A_mi$p.I, sim.test_P_pa_A_mi$p.I, 
                                    sim.test_P_tr_A_mi$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_ac_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_C_cr_A_mi, "D", "") 
ecospat.plot.overlap.test(sim.test_C_in_A_mi, "D", "") 
ecospat.plot.overlap.test(sim.test_C_la_A_mi, "D", "") 
ecospat.plot.overlap.test(sim.test_C_mo_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_A_mi, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_A_mi, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_ac_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_cr_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_rh_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_A_mi, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_A_mi, "I", "")

# Other 9 species vs Crocodylus acutus ###############################################################################################
sim.test_C_cr_C_ac <- ecospat.niche.similarity.test (grid_clim_C_cr, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_C_in_C_ac <- ecospat.niche.similarity.test (grid_clim_C_in, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_C_la_C_ac <- ecospat.niche.similarity.test (grid_clim_C_la, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_C_mo_C_ac <- ecospat.niche.similarity.test (grid_clim_C_mo, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_C_rh_C_ac <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_C_ya_C_ac <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_M_ni_C_ac <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_P_pa_C_ac <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_ac, rep=100,  ncores = 50)
sim.test_P_tr_C_ac <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_ac, rep=100,  ncores = 50)

sim.test_C_ac_2_list <- data.frame (sim.test_sp_vs_C_ac = c(1, 1, sim.test_C_cr_C_ac$obs$D, sim.test_C_in_C_ac$obs$D, sim.test_C_la_C_ac$obs$D,
                                    sim.test_C_mo_C_ac$obs$D, sim.test_C_rh_C_ac$obs$D, sim.test_C_ya_C_ac$obs$D, sim.test_M_ni_C_ac$obs$D, 
                                    sim.test_P_pa_C_ac$obs$D, sim.test_P_tr_C_ac$obs$D, 1, 1, sim.test_C_cr_C_ac$obs$I, sim.test_C_in_C_ac$obs$I, 
                                    sim.test_C_la_C_ac$obs$I, sim.test_C_mo_C_ac$obs$I, sim.test_C_rh_C_ac$obs$I, sim.test_C_ya_C_ac$obs$I,
                                    sim.test_M_ni_C_ac$obs$I, sim.test_P_pa_C_ac$obs$I, sim.test_P_tr_C_ac$obs$I, 1, 1, sim.test_C_cr_C_ac$p.D, 
                                    sim.test_C_in_C_ac$p.D, sim.test_C_la_C_ac$p.D, sim.test_C_mo_C_ac$p.D, sim.test_C_rh_C_ac$p.D, 
                                    sim.test_C_ya_C_ac$p.D, sim.test_M_ni_C_ac$p.D, sim.test_P_pa_C_ac$p.D, sim.test_P_tr_C_ac$p.D,
                                    1, 1, sim.test_C_cr_C_ac$p.I, sim.test_C_in_C_ac$p.I, sim.test_C_la_C_ac$p.I, sim.test_C_mo_C_ac$p.I, 
                                    sim.test_C_rh_C_ac$p.I, sim.test_C_ya_C_ac$p.I, sim.test_M_ni_C_ac$p.I, sim.test_P_pa_C_ac$p.I,
                                    sim.test_P_tr_C_ac$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_cr_C_ac, "D", "") 
ecospat.plot.overlap.test(sim.test_C_in_C_ac, "D", "") 
ecospat.plot.overlap.test(sim.test_C_la_C_ac, "D", "") 
ecospat.plot.overlap.test(sim.test_C_mo_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_ac, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_ac, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_cr_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_C_in_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_ac, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_ac, "I", "")

# Other 8 species vs Caiman crocodilus ###############################################################################################
sim.test_C_in_C_cr <- ecospat.niche.similarity.test (grid_clim_C_in, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_C_la_C_cr <- ecospat.niche.similarity.test (grid_clim_C_la, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_C_mo_C_cr <- ecospat.niche.similarity.test (grid_clim_C_mo, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_C_rh_C_cr <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_C_ya_C_cr <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_M_ni_C_cr <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_P_pa_C_cr <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_cr, rep=100,  ncores = 50)
sim.test_P_tr_C_cr <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_cr, rep=100,  ncores = 50)

sim.test_C_cr_2_list <- data.frame (sim.test_sp_vs_C_cr = c(1, 1, 1, sim.test_C_in_C_cr$obs$D, sim.test_C_la_C_cr$obs$D, sim.test_C_mo_C_cr$obs$D,
                                    sim.test_C_rh_C_cr$obs$D, sim.test_C_ya_C_cr$obs$D, sim.test_M_ni_C_cr$obs$D, sim.test_P_pa_C_cr$obs$D, 
                                    sim.test_P_tr_C_cr$obs$D, 1, 1, 1, sim.test_C_in_C_cr$obs$I, sim.test_C_la_C_cr$obs$I, sim.test_C_mo_C_cr$obs$I,
                                    sim.test_C_rh_C_cr$obs$I, sim.test_C_ya_C_cr$obs$I, sim.test_M_ni_C_cr$obs$I, sim.test_P_pa_C_cr$obs$I, 
                                    sim.test_P_tr_C_cr$obs$I, 1, 1, 1, sim.test_C_in_C_cr$p.D, sim.test_C_la_C_cr$p.D, sim.test_C_mo_C_cr$p.D,
                                    sim.test_C_rh_C_cr$p.D, sim.test_C_ya_C_cr$p.D, sim.test_M_ni_C_cr$p.D, sim.test_P_pa_C_cr$p.D, 
                                    sim.test_P_tr_C_cr$p.D, 1, 1, 1, sim.test_C_in_C_cr$p.I, sim.test_C_la_C_cr$p.I, sim.test_C_mo_C_cr$p.I,
                                    sim.test_C_rh_C_cr$p.I, sim.test_C_ya_C_cr$p.I, sim.test_M_ni_C_cr$p.I, sim.test_P_pa_C_cr$p.I, 
                                    sim.test_P_tr_C_cr$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_in_C_cr, "D", "") 
ecospat.plot.overlap.test(sim.test_C_la_C_cr, "D", "") 
ecospat.plot.overlap.test(sim.test_C_mo_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_cr, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_cr, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_in_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_C_la_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_cr, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_cr, "I", "")

# Other 7 species vs Crocodylus intermedius ##########################################################################################
sim.test_C_la_C_in <- ecospat.niche.similarity.test (grid_clim_C_la, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_C_mo_C_in <- ecospat.niche.similarity.test (grid_clim_C_mo, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_C_rh_C_in <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_C_ya_C_in <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_M_ni_C_in <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_P_pa_C_in <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_in, rep=100,  ncores = 50)
sim.test_P_tr_C_in <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_in, rep=100,  ncores = 50)

sim.test_C_in_2_list <- data.frame (sim.test_sp_vs_C_in = c(1, 1, 1, 1, sim.test_C_la_C_in$obs$D, sim.test_C_mo_C_in$obs$D, 
                                    sim.test_C_rh_C_in$obs$D, sim.test_C_ya_C_in$obs$D, sim.test_M_ni_C_in$obs$D, sim.test_P_pa_C_in$obs$D, 
                                    sim.test_P_tr_C_in$obs$D, 1, 1, 1, 1, sim.test_C_la_C_in$obs$I, sim.test_C_mo_C_in$obs$I, 
                                    sim.test_C_rh_C_in$obs$I, sim.test_C_ya_C_in$obs$I, sim.test_M_ni_C_in$obs$I, sim.test_P_pa_C_in$obs$I, 
                                    sim.test_P_tr_C_in$obs$I, 1, 1, 1, 1, sim.test_C_la_C_in$p.D, sim.test_C_mo_C_in$p.D, sim.test_C_rh_C_in$p.D,
                                    sim.test_C_ya_C_in$p.D, sim.test_M_ni_C_in$p.D, sim.test_P_pa_C_in$p.D, sim.test_P_tr_C_in$p.D, 1, 1, 1, 1, 
                                    sim.test_C_la_C_in$p.I, sim.test_C_mo_C_in$p.I, sim.test_C_rh_C_in$p.I, sim.test_C_ya_C_in$p.I,
                                    sim.test_M_ni_C_in$p.I, sim.test_P_pa_C_in$p.I, sim.test_P_tr_C_in$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_la_C_in, "D", "") 
ecospat.plot.overlap.test(sim.test_C_mo_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_in, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_in, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_la_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_C_mo_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_in, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_in, "I", "")

# Other 6 species vs Caiman latirostris ##############################################################################################
sim.test_C_mo_C_la <- ecospat.niche.similarity.test (grid_clim_C_mo, grid_clim_C_la, rep=100,  ncores = 50)
sim.test_C_rh_C_la <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_C_la, rep=100,  ncores = 50)
sim.test_C_ya_C_la <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_la, rep=100,  ncores = 50)
sim.test_M_ni_C_la <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_la, rep=100,  ncores = 50)
sim.test_P_pa_C_la <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_la, rep=100,  ncores = 50)
sim.test_P_tr_C_la <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_la, rep=100,  ncores = 50)

sim.test_C_la_2_list <- data.frame (sim.test_sp_vs_C_la = c(1, 1, 1, 1, 1, sim.test_C_mo_C_la$obs$D, sim.test_C_rh_C_la$obs$D, 
                                    sim.test_C_ya_C_la$obs$D, sim.test_M_ni_C_la$obs$D, sim.test_P_pa_C_la$obs$D, sim.test_P_tr_C_la$obs$D,
                                    1, 1, 1, 1, 1, sim.test_C_mo_C_la$obs$I, sim.test_C_rh_C_la$obs$I, sim.test_C_ya_C_la$obs$I, 
                                    sim.test_M_ni_C_la$obs$I, sim.test_P_pa_C_la$obs$I, sim.test_P_tr_C_la$obs$I, 1, 1, 1, 1, 1,
                                    sim.test_C_mo_C_la$p.D, sim.test_C_rh_C_la$p.D, sim.test_C_ya_C_la$p.D, sim.test_M_ni_C_la$p.D,
                                    sim.test_P_pa_C_la$p.D, sim.test_P_tr_C_la$p.D, 1, 1, 1, 1, 1, sim.test_C_mo_C_la$p.I, 
                                    sim.test_C_rh_C_la$p.I, sim.test_C_ya_C_la$p.I, sim.test_M_ni_C_la$p.I, sim.test_P_pa_C_la$p.I, 
                                    sim.test_P_tr_C_la$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_mo_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_la, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_la, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_mo_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_C_rh_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_la, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_la, "I", "")

# Other 5 species vs Crocodylus moreletii ############################################################################################
sim.test_C_rh_C_mo <- ecospat.niche.similarity.test (grid_clim_C_rh, grid_clim_C_mo, rep=100,  ncores = 50)
sim.test_C_ya_C_mo <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_mo, rep=100,  ncores = 50)
sim.test_M_ni_C_mo <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_mo, rep=100,  ncores = 50)
sim.test_P_pa_C_mo <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_mo, rep=100,  ncores = 50)
sim.test_P_tr_C_mo <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_mo, rep=100,  ncores = 50)

sim.test_C_mo_2_list <- data.frame (sim.test_sp_vs_C_mo = c(1, 1, 1, 1, 1, 1, sim.test_C_rh_C_mo$obs$D, sim.test_C_ya_C_mo$obs$D, 
                                    sim.test_M_ni_C_mo$obs$D, sim.test_P_pa_C_mo$obs$D, sim.test_P_tr_C_mo$obs$D, 1, 1, 1, 1, 1, 1, 
                                    sim.test_C_rh_C_mo$obs$I, sim.test_C_ya_C_mo$obs$I, sim.test_M_ni_C_mo$obs$I, sim.test_P_pa_C_mo$obs$I, 
                                    sim.test_P_tr_C_mo$obs$I, 1, 1, 1, 1, 1, 1, sim.test_C_rh_C_mo$p.D, sim.test_C_ya_C_mo$p.D, 
                                    sim.test_M_ni_C_mo$p.D, sim.test_P_pa_C_mo$p.D, sim.test_P_tr_C_mo$p.D, 1, 1, 1, 1, 1, 1, 
                                    sim.test_C_rh_C_mo$p.I, sim.test_C_ya_C_mo$p.I, sim.test_M_ni_C_mo$p.I, sim.test_P_pa_C_mo$p.I, 
                                    sim.test_P_tr_C_mo$p.I))
######## D ########
ecospat.plot.overlap.test(sim.test_C_rh_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_mo, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_mo, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_rh_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_C_ya_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_mo, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_mo, "I", "")

# Other 4 species vs Crocodylus rhombifer ############################################################################################
sim.test_C_ya_C_rh <- ecospat.niche.similarity.test (grid_clim_C_ya, grid_clim_C_rh, rep=100,  ncores = 50)
sim.test_M_ni_C_rh <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_rh, rep=100,  ncores = 50)
sim.test_P_pa_C_rh <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_rh, rep=100,  ncores = 50)
sim.test_P_tr_C_rh <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_rh, rep=100,  ncores = 50)

sim.test_C_rh_2_list <- data.frame (sim.test_sp_vs_C_rh = c(1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_C_rh$obs$D, sim.test_M_ni_C_rh$obs$D,
                                    sim.test_P_pa_C_rh$obs$D, sim.test_P_tr_C_rh$obs$D, 1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_C_rh$obs$I, 
                                    sim.test_M_ni_C_rh$obs$I, sim.test_P_pa_C_rh$obs$I, sim.test_P_tr_C_rh$obs$I, 1, 1, 1, 1, 1, 1, 1, 
                                    sim.test_C_ya_C_rh$p.D, sim.test_M_ni_C_rh$p.D, sim.test_P_pa_C_rh$p.D, sim.test_P_tr_C_rh$p.D,
                                    1, 1, 1, 1, 1, 1, 1, sim.test_C_ya_C_rh$p.I, sim.test_M_ni_C_rh$p.I, sim.test_P_pa_C_rh$p.I, 
                                    sim.test_P_tr_C_rh$p.I))
####### D ########
ecospat.plot.overlap.test(sim.test_C_ya_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_rh, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_rh, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_C_ya_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_M_ni_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_rh, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_rh, "I", "")

# Other 3 species vs Caiman yacare ###################################################################################################
sim.test_M_ni_C_ya <- ecospat.niche.similarity.test (grid_clim_M_ni, grid_clim_C_ya, rep=100,  ncores = 50)
sim.test_P_pa_C_ya <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_C_ya, rep=100,  ncores = 50)
sim.test_P_tr_C_ya <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_C_ya, rep=100,  ncores = 50)

sim.test_C_ya_2_list <- data.frame (sim.test_sp_vs_C_ya = c(1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_C_ya$obs$D, sim.test_P_pa_C_ya$obs$D, 
                                    sim.test_P_tr_C_ya$obs$D, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_C_ya$obs$I, sim.test_P_pa_C_ya$obs$I, 
                                    sim.test_P_tr_C_ya$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_C_ya$p.D, sim.test_P_pa_C_ya$p.D, 
                                    sim.test_P_tr_C_ya$p.D, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_M_ni_C_ya$p.I, sim.test_P_pa_C_ya$p.I, 
                                    sim.test_P_tr_C_ya$p.I))
####### D ########
ecospat.plot.overlap.test(sim.test_M_ni_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_ya, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_ya, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_M_ni_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_P_pa_C_ya, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_C_ya, "I", "")

# Other 2 species vs Melanosuchus niger ##############################################################################################
sim.test_P_pa_M_ni <- ecospat.niche.similarity.test (grid_clim_P_pa, grid_clim_M_ni, rep=100,  ncores = 50)
sim.test_P_tr_M_ni <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_M_ni, rep=100,  ncores = 50)

sim.test_M_ni_2_list <- data.frame (sim.test_sp_vs_M_ni = c(1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_pa_M_ni$obs$D, sim.test_P_tr_M_ni$obs$D,
                                    1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_pa_M_ni$obs$I, sim.test_P_tr_M_ni$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                    sim.test_P_pa_M_ni$p.D, sim.test_P_tr_M_ni$p.D, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_pa_M_ni$p.I, 
                                    sim.test_P_tr_M_ni$p.I))
####### D ########
ecospat.plot.overlap.test(sim.test_P_pa_M_ni, "D", "")
ecospat.plot.overlap.test(sim.test_P_tr_M_ni, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_P_pa_M_ni, "I", "")
ecospat.plot.overlap.test(sim.test_P_tr_M_ni, "I", "")

# Other 1 species vs Paleosuchus palpebrosus #########################################################################################
sim.test_P_tr_P_pa <- ecospat.niche.similarity.test (grid_clim_P_tr, grid_clim_P_pa, rep=100,  ncores = 50)

sim.test_P_pa_2_list <- data.frame (sim.test_sp_vs_P_pa = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_tr_P_pa$obs$D, 1, 1, 1, 1, 1, 1, 
                                    1, 1, 1, 1, sim.test_P_tr_P_pa$obs$I, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_tr_P_pa$p.D, 1, 1, 
                                    1, 1, 1, 1, 1, 1, 1, 1, sim.test_P_tr_P_pa$p.I))
####### D ########
ecospat.plot.overlap.test(sim.test_P_tr_P_pa, "D", "")
######## I ########
ecospat.plot.overlap.test(sim.test_P_tr_P_pa, "I", "")

######################################################################################################################################
# PHYLOGENETIC INDEPENDENT CONTRAST TO ASSESS NICHE CONSERVATISM
# The phylogenetic tree we will use was adapted from Oaks (2011). This tree can be found inside the folder "Crocs_tree" saved as the 
# name of "Croc_tree_nexus.tre".
######################################################################################################################################
Croc_tree_nex <- read.nexus(".../Crocs_tree/Croc_tree_nexus.tre")
Croc_tree_nex
plotTree(Croc_tree_nex, ftype = "i", fsize = 1, lwd =1, mar = c(2, 1.1, 3.1,1.1))
axisPhylo()
nodelabels(cex = 0.5) #label the nodes in the phylogenetic tree
Ntip(Croc_tree_nex) 
edgelabels(round(Croc_tree_nex$edge.length, 2), pos = 3, frame = "none", cex = 0.6)

Croc_1_dist_phylo <- cophenetic.phylo(Croc_tree_nex) # Estimate the phylogenetic patristic distance
Croc_1_dist_phylo
class(Croc_1_dist_phylo)
write.csv(Croc_1_dist_phylo, ".../Tables_for_R/05_Phylo_distance.csv") # Matrix using species names

######################################################################################################################################
# MANTEL TEST
# Perform a correlation between Niche overlap (Scchoener's D, Hellinger's I, Geographic overlap) and Phylogenetic distance.
# To this analysis, we will use the tables (.csv) found in the folder "Tables_for_R".
######################################################################################################################################
# Load the tables (.csv) and convert these dataframes into distance matrix 
Croc_phylo <- read.csv(".../Tables_for_R/05_Phylo_distance.csv", row.names = 1)
class(Croc_phylo) # Needs to be "dataframe"
Croc_phylo_matrix <- as.dist(Croc_phylo) # Convert dataframe to distance matrix
any(is.na(Croc_phylo_matrix)) # Check that there are no NA data (should return "FALSE")
class(Croc_phylo_matrix) # Need to be "dist"

Croc_Schoener_D <- read.csv(".../Tables_for_R/06_Schoener_D.csv", row.names = 1)
Croc_Schoener_D_matrix <- as.dist(Croc_Schoener_D)
any(is.na(Croc_Schoener_D_matrix))

Croc_Hellinger_I <- read.csv(".../Tables_for_R/07_Hellinger_I.csv", row.names = 1)
Croc_Hellinger_I_matrix <- as.dist(Croc_Hellinger_I)
any(is.na(Croc_Hellinger_I_matrix)) 

Croc_Geo <- read.csv(".../Tables_for_R/08_Geographic.csv", row.names = 1)
class(Croc_Geo)
Croc_Geo_matrix <- as.dist(Croc_Geo)
any(is.na(Croc_Geo_matrix)) 

# Perform the Mantel test using the previous distance matrices (use package "vegan")
Mantel_phylo_D <- mantel(Croc_phylo_matrix, Croc_Schoener_D_matrix, method = "spearman", permutations = 1000, na.rm = FALSE)
Mantel_Phylo_I <- mantel(Croc_phylo_matrix, Croc_Hellinger_I_matrix, method = "spearman", permutations = 1000, na.rm = FALSE)
Mantel_Phylo_Geo <- mantel(Croc_phylo_matrix, Croc_Geo_matrix, method = "spearman", permutations = 1000, na.rm = FALSE)

Mantel_phylo_D
Mantel_Phylo_I
Mantel_Phylo_Geo

# Extract the permutation results and observed statistic for each test
observed_D <- Mantel_phylo_D$statistic
permuted_D <- Mantel_phylo_D$perm

observed_I <- Mantel_Phylo_I$statistic
permuted_I <- Mantel_Phylo_I$perm

observed_Geo <- Mantel_Phylo_Geo$statistic
permuted_Geo <- Mantel_Phylo_Geo$perm

# Two-tailed p-value calculation (both upper and lower)
# Schoener's D
p_value_D <- min(mean(permuted_D >= observed_D), mean(permuted_D <= observed_D))

# Hellinger I
p_value_I <- min(mean(permuted_I >= observed_I), mean(permuted_I <= observed_I))

# Geographic overlap
p_value_Geo <- min(mean(permuted_Geo >= observed_Geo), mean(permuted_Geo <= observed_Geo))

# Print the two-tailed p-values for each test
cat("Two-tailed p-value for Schoener's D:", p_value_D, "\n")
cat("Two-tailed p-value for Hellinger I:", p_value_I, "\n")
cat("Two-tailed p-value for Geographic overlap:", p_value_Geo, "\n")

# Figure 3.7 A, B, C.
mantel_stats_D <- Mantel_phylo_D$statistic  # Mantel statistic for Schoener's D
permutation_stats_D <- Mantel_phylo_D$perm  # Permutation r values for Schoener's D
min(permutation_stats_D)
max(permutation_stats_D)

mantel_stats_I <- Mantel_Phylo_I$statistic  # Mantel statistic for Hellinger I
permutation_stats_I <- Mantel_Phylo_I$perm  # Permutation r values for Hellinger I
min(permutation_stats_I)
max(permutation_stats_I)

mantel_stats_Geo <- Mantel_Phylo_Geo$statistic  # Mantel statistic for Geographic overlap
permutation_stats_Geo <- Mantel_Phylo_Geo$perm  # Permutation r values for Geographic overlap
min(permutation_stats_Geo)
max(permutation_stats_Geo)

# Figure 3.7 - A - Geographic overlap
Figure_3_7_A <- ggplot(data.frame(permutation_stats = permutation_stats_Geo), aes(x = permutation_stats)) +
  geom_density(aes(y = ..density.. / max(..density..)), fill = "gray", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mantel_stats_Geo), color = "black", linetype = "dashed", size = 1.5) +
  annotate("text", x = -0.55, y = 1, label = "A", hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = 0.35, y = 0.95, label = paste("Observed rho (ρ) =", round(mantel_stats_Geo, 3), "\n",
                                                     "p-value =", round(p_value_Geo, 3)),
           fontface = "bold.italic",  color = "darkgray", size = 11, hjust = 0.5, vjust = 0.5) +
  labs(x = "Permutation ρ (Phylogenetic distance - Geographic overlap)", y = "Density") +
  theme_classic(base_size = 17) +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20),  
        axis.text.y = element_text(size = 20),  
        panel.grid = element_blank(),  
        plot.title = element_blank(),  
        plot.subtitle = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1), limits = c(-0.55, 0.5)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))
Figure_3_7_A
ggsave(".../Figures/Figure_3_7_A.jpg", plot = Figure_3_7_A, width = 15, height = 10, dpi = 300)

# Figure 3.7 - B - Schoener's D
Figure_3_7_B <- ggplot(data.frame(permutation_stats = permutation_stats_D), aes(x = permutation_stats)) +
  geom_density(aes(y = ..density.. / max(..density..)), fill = "gray", color = "black", alpha = 0.7) + 
  geom_vline(aes(xintercept = mantel_stats_D), color = "black", linetype = "dashed", size = 1.5) +
  annotate("text", x = -0.55, y = 1, label = "B",  hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = 0.35, y = 0.95, label = paste("Observed rho (ρ) =", format(round(mantel_stats_D, 3), scientific = FALSE) , "\n",
                                                     "p-value =", round(p_value_D, 3)),
           fontface = "bold.italic", color = "darkgray", size = 11, hjust = 0.5, vjust = 0.5) +
  labs(x = "Permutation ρ (Phylogenetic distance -  Schoener's D)", y = "Density") +
  theme_classic(base_size = 17) +  
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20),  
        axis.text.y = element_text(size = 20),  
        panel.grid = element_blank(),  
        plot.title = element_blank(),  
        plot.subtitle = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1), limits = c(-0.55, 0.5)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))
Figure_3_7_B
ggsave(".../Figures/Figure_3_7_B.jpg", plot = Figure_3_7_B, width = 15, height = 10, dpi = 300)

# Figure 3.7 - C - Hellinger's I
Figure_3_7_C <- ggplot(data.frame(permutation_stats = permutation_stats_I), aes(x = permutation_stats)) +
  geom_density(aes(y = ..density.. / max(..density..)), fill = "gray", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mantel_stats_I), color = "black", linetype = "dashed", size = 1.5) +
  annotate("text", x = -0.55, y = 1, label = "C", hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = 0.35, y = 0.95, label = paste("Observed rho (ρ) =", round(mantel_stats_I, 3), "\n",
                                                     "p-value =", round(p_value_I, 3)),
           fontface = "bold.italic", color = "darkgray", size = 11, hjust = 0.5, vjust = 0.5) +
  labs(x = "Permutation ρ (Phylogenetic distance -  Hellinger's I)", y = "Density") +
  theme_classic(base_size = 17) +  
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20),  
        axis.text.y = element_text(size = 20),  
        panel.grid = element_blank(),  
        plot.title = element_blank(),  
        plot.subtitle = element_text(size = 18)) +
  scale_x_continuous(breaks = seq(-0.5, 0.5, by = 0.1), limits = c(-0.55, 0.5)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))  
Figure_3_7_C
ggsave(".../Figures/Figure_3_7_C.jpg", plot = Figure_3_7_C, width = 15, height = 10, dpi = 300)

######################################################################################################################################
# ANALYSIS OF PHYLOGENETIC INDEPENDENT CONTRASTS
# Correlation between Niche overlap (Scchoener's D, Hellinger's I, Geographic overlap) and Phylogenetic distance.
# To this analysis, we will use the tables (.csv) found in the folder "Tables_for_R".
######################################################################################################################################
# Load data
Croc_tree <- read.nexus("E:/THESIS_CH3/Tree/Oak/Croc_tree_nexus.tre")
Croc_phylo <- read.csv("E:/THESIS_CH3/Tables/05_Phylo_distance.csv", row.names = 1)
Croc_Schoener_D <- read.csv("E:/THESIS_CH3/Tables/06_Schoener_D.csv", row.names = 1)
Croc_Hellinger_I <- read.csv("E:/THESIS_CH3/Tables/07_Hellinger_I.csv", row.names = 1)
Croc_Geo <- read.csv("E:/THESIS_CH3/Tables/08_Geographic.csv", row.names = 1)

# Modify row names of the matrices to match the tree by replacing underscores with spaces
rownames(Croc_phylo) <- gsub("_", " ", rownames(Croc_phylo))
rownames(Croc_Schoener_D) <- gsub("_", " ", rownames(Croc_Schoener_D))
rownames(Croc_Hellinger_I) <- gsub("_", " ", rownames(Croc_Hellinger_I))
rownames(Croc_Geo) <- gsub("_", " ", rownames(Croc_Geo))

# Verify that the species are properly matched between the data and the phylogeny.
species_match <- match(rownames(Croc_Schoener_D), Croc_tree$tip.label)
all(!is.na(species_match))  # Should return TRUE if all species match

# Loop through each row and set diagonal values to NA
for (i in 1:nrow(Croc_Schoener_D)) {
  Croc_Schoener_D[i, i] <- NA
  Croc_Hellinger_I[i, i] <- NA
  Croc_phylo[i, i] <- NA
  Croc_Geo[i, i] <- NA
  }

# Calculate the average Schoener's D and average Hellinger's I for each species
Croc_Schoener_D_avg <- rowMeans(Croc_Schoener_D, na.rm = TRUE)
Croc_Hellinger_I_avg <- rowMeans(Croc_Hellinger_I, na.rm = TRUE)
Croc_phylo_avg <- rowMeans(Croc_phylo, na.rm = TRUE)
Croc_Geo_avg <- rowMeans(Croc_Geo, na.rm = TRUE)

length(Croc_Schoener_D_avg)  # Should match the number of species in Croc_tree$tip.label
length(Croc_Hellinger_I_avg)  
length(Croc_phylo_avg)
length(Croc_Geo_avg)

# Prepare a dataframe for PIC analysis
species_names <- Croc_tree$tip.label

analysis_data <- data.frame(Species = species_names,
                            Avg_Schoener_D = Croc_Schoener_D_avg[species_names],
                            Avg_Hellinger_I = Croc_Hellinger_I_avg[species_names],
                            Avg_Geo = Croc_Geo_avg[species_names], 
                            Avg_phylo_dist = Croc_phylo_avg[species_names])
analysis_data

# Step 6: Perform PIC Analysis.
is.ultrametric(Croc_tree)  # Should return TRUE

pic_Schoener_D <- pic(analysis_data$Avg_Schoener_D, Croc_tree)
pic_Hellinger_I <- pic(analysis_data$Avg_Hellinger_I, Croc_tree)
pic_phylo_dist <- pic(analysis_data$Avg_phylo_dist, Croc_tree)
pic_Geo <- pic(analysis_data$Avg_Geo, Croc_tree)

pic_data <- data.frame(PIC_Schoener_D = pic_Schoener_D,
                       PIC_Hellinger_I = pic_Hellinger_I,
                       PIC_Geo = pic_Geo,
                       PIC_phylo_dist = pic_phylo_dist)
pic_data

# Analyze relation between PIC Phylogenetic distance vs PIC overlaps (D, I, Geo)
model_PIC_phylo_schoener <- lm(PIC_Schoener_D ~ PIC_phylo_dist, data = pic_data)
model_PIC_phylo_hellinger <- lm(PIC_Hellinger_I ~ PIC_phylo_dist, data = pic_data)
model_PIC_phylo_geo <- lm(PIC_Geo ~ PIC_phylo_dist, data = pic_data)

# Extract coefficients and R^2
coeffs_PIC_schoener <- coef(model_PIC_phylo_schoener)
r_squared_PIC_schoener <- summary(model_PIC_phylo_schoener)$r.squared

coeffs_PIC_hellinger <- coef(model_PIC_phylo_hellinger)
r_squared_PIC_hellinger <- summary(model_PIC_phylo_hellinger)$r.squared

coeffs_PIC_geo <- coef(model_PIC_phylo_geo)
r_squared_PIC_geo <- summary(model_PIC_phylo_geo)$r.squared

# Spearman correlation tests
cor_PIC_phylo_schoener <- cor.test(pic_data$PIC_Schoener_D, pic_data$PIC_phylo_dist, method = "spearman")
cor_PIC_phylo_Hellinger <- cor.test(pic_data$PIC_Hellinger_I, pic_data$PIC_phylo_dist, method = "spearman")
cor_PIC_phylo_Geo <- cor.test(pic_data$PIC_Geo, pic_data$PIC_phylo_dist, method = "spearman")

# Output the correlation results
cor_PIC_phylo_schoener
cor_PIC_phylo_Hellinger
cor_PIC_phylo_Geo

# Figure 3.7 - D - Relationship between PIC Phylogenetic distance and PIC Geographic overlap
Figure_3_7_D <- ggplot(pic_data, aes(x = PIC_phylo_dist, y = PIC_Geo)) +
  geom_point(size = 3, color = "black") +  
  geom_smooth(method = "lm", se = FALSE, color = "darkgray", size = 1, linetype = "dashed") + 
  annotate("text", x = -5, y = 1.4, label = "D", hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = -2.5, y = 1.1, label = paste("rho (ρ) = ", round(cor_PIC_phylo_Geo$estimate, 3), "\n",
                                                    "p-value = ", round(cor_PIC_phylo_Geo$p.value, 3)),
           hjust = 0.5, vjust = 0.5, size = 11, fontface = "bold.italic", lineheight = 0.8, parse = FALSE, color = "darkgray") +
  labs(x = "PIC - Phylogenetic patristic distance", y = "PIC - Geographic overlap") +  
  theme_classic() + 
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(-5, 1, by = 1), limits = c(-5,1)) +
  scale_y_continuous(breaks = seq(-1, 1.4, by = 0.4), limits = c(-1, 1.4))
Figure_3_7_D
ggsave(".../Figures/Figure_3_7_D.jpg", plot = Figure_3_7_D, width = 15, height = 10, dpi = 300)

# Figure 3.7 - E - Relationship between PIC Phylogenetic distance and PIC Schoener's D
Figure_3_7_E <- ggplot(pic_data, aes(x = PIC_phylo_dist, y = PIC_Schoener_D)) +
  geom_point(size = 3, color = "black") +  
  geom_smooth(method = "lm", se = FALSE, color = "darkgray", size = 1, linetype = "dashed") +  
  annotate("text", x = -5, y = 0.03, label = "E", hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = -2.5, y = 0.025, label = paste("rho (ρ) = ", round(cor_PIC_phylo_schoener$estimate, 3), "\n",
                                                      "p-value = ", round(cor_PIC_phylo_schoener$p.value, 3)),
           hjust = 0.5, vjust = 0.5, size = 11, fontface = "bold.italic", lineheight = 0.8, parse = FALSE, color = "darkgray") +
  labs(x = "PIC - Phylogenetic patristic distance", y = "PIC - Schoener's D") +  
  theme_classic() +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(-5, 1, by = 1), limits = c(-5,1)) +
  scale_y_continuous(breaks = seq(-0.02, 0.03, by = 0.01), limits = c(-0.025,0.031))
Figure_3_7_E
ggsave(".../Figures/Figure_3_7_E.jpg", plot = Figure_3_7_E, width = 15, height = 10, dpi = 300)

# Figure 3.7 - F - Relationship between PIC Phylogenetic distance and PIC Hellinger's I
Figure_3_7_F <- ggplot(pic_data, aes(x = PIC_phylo_dist, y = PIC_Hellinger_I)) +
  geom_point(size = 3, color = "black") +  
  geom_smooth(method = "lm", se = FALSE, color = "darkgray", size = 1, linetype = "dashed") + 
  annotate("text", x = -5, y = 0.08, label = "F", hjust = 0.5, vjust = 0.5, size = 20, fontface = "bold") +
  annotate("text", x = -2.5, y = 0.07, label = paste("rho (ρ) = ", round(cor_PIC_phylo_Hellinger$estimate, 3), "\n",
                                                     "p-value = ", round(cor_PIC_phylo_Hellinger$p.value, 3)),
           hjust = 0.5, vjust = 0.5, size = 11, fontface = "bold.italic", lineheight = 0.8, parse = FALSE, color = "darkgray") +
  labs(x = "PIC - Phylogenetic patristic distance", y = "PIC - Hellinger's I") +  
  theme_classic() +
  theme(axis.text = element_text(color = "black"), 
        axis.title = element_text(face = "bold"),
        axis.title.x = element_text(size = 34, face = "bold", vjust = -1.5, margin = margin(b = 25)),
        axis.title.y = element_text(size = 34, face = "bold", vjust = 4, margin = margin(l = 25)),
        axis.text.x = element_text(size = 20),  # Change x-axis tick size
        axis.text.y = element_text(size = 20),  # Change y-axis tick size
        panel.grid = element_blank()) +
  scale_x_continuous(breaks = seq(-5, 1, by = 1), limits = c(-5,1)) +
  scale_y_continuous(breaks = seq(-0.06, 0.08, by = 0.02), limits = c(-0.06,0.08))
Figure_3_7_F
ggsave(".../Figures/Figure_3_7_F.jpg", plot = Figure_3_7_F, width = 15, height = 10, dpi = 300)

# Arrangement for the final figure included in the manuscript

Figure_3_7_top <- grid.arrange(Figure_3_7_A, Figure_3_7_B, Figure_3_7_C, ncol = 3)
Figure_3_7_bottom <- grid.arrange(Figure_3_7_D, Figure_3_7_E, Figure_3_7_F, ncol = 3)

# Final plot - Figure 3.7
Figure_3_7_final <- grid.arrange(Figure_3_7_top, nullGrob(), Figure_3_7_bottom, nrow = 3, heights = c(1, 0.1, 1))
ggsave(".../Figures/Figure_3_7_final.jpg", plot = Figure_3_7_final, width = 44, height = 23, dpi = 400)

######################################################################################################################################
# END
######################################################################################################################################