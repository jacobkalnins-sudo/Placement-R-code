## Filtering SNP data from DArT ----
## Example script July 2026 by J OHare


## Script info and package installations ----
## This script uses "dartR" to filter the data for quality,
## including the topics covered in BIOL2110 prac "Data Quality in Genomics"

## Note: the values will need some tweaking here and there,
## as each dataset is unique and you will need to investigate variations.
## Eg increasing the "call rate" values too high (being too strict) or too low (being too lenient)
## can filter out the true signal or introduce noise that is not biologically relevant

## Hence - there is no "one-size-fits-all"
## we often explore a few versions of the filtering before deciding
## We are trying to balance best quality data that we can get, with retaining as many inds and SNPs as possible

## Okay - into the data!!
## The script has bits of info, for more details, check the tutorials linked in here:
# https://github.com/green-striped-gecko/dartRverse

## 1 Loading Packages and set working directory: ----

# install first as necessary: 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("LEA")
# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("snpStats")


# The dartRverse website above also gives the best and thorough installation help
install.packages("dartRverse")
install.packages("dartR.data")
install.packages("dartR.base")
install.packages("HardyWeinberg")

## load basic packages needed ----

library(dartRverse)
library(HardyWeinberg)
library(SNPRelate)
library(dplyr)
library(directlabels)
library(ggplot2)
library(ggrepel)

## Note: 
## 2 set working directory to a folder where you have all the data saved
# Click options above: 
# Session > Set Working Directory > XX
#                                   To source file location (if this R code is saved in a folder near the dataset) 
#                                   OR 
#                                   Choose directory (to interactively choose the location of the dataset you'll use)


## 3 prep data for standard base filtering ----
# import dataset and code the pops using the metadta file (eg in "OrderAppendix.." folder)

## Note that first collum of numbers is for the first set of data and the seccond for the seccond set. When only 1 set it is the seccond set

unfiltered <- gl.read.dart(
  filename = "Report_DE23-8931_7_moreOrders_SNP_mapping_2.csv",
  ind.metafile = "RedGums_DE23-8931_metadata_merged.csv")


unfiltered
# <--- AS YOU GO THROUGH, UPDATE THE ANNOTATIONS FOR NUMBER OF INDIVIDUALS AND SNPS

# 560 individuals and 39403  SNPs in a 2 row format using Report_DE23-8113_SNP_1.csv 

#inspect

nInd(unfiltered)
#560, 1269
nLoc(unfiltered)
#39403, 23819
nPop(unfiltered)
#21, 57
popNames(unfiltered)
# "1_QLD_WNR"     "10_QLD_BRG"    "11_QLD_WON"    "12_QLD_TMF_P3" "13_QLD_BRV"    "14_QLD_MOO"    "15_QLD_WAR"    "16_NSW_SEL"    "17_NSW_NEW"    "18_NSW_WIN"    "19_NSW_CAL_2"  "2_QLD_HEL" "20_NSW_YUR"    "21_VIC_GIP"    "3_QLD_MIT"     "4_QLD_BER_P1"  "5_QLD_BUR"     "6_QLD_MID"     "7_QLD_MIR"     "8_QLD_BIL"     "9_QLD_WID" 
# 1_QLD_WNR"               "10_QLD_BRG"              "11_QLD_WON"              "12_QLD_TMF_P3"           "13_QLD_BRV"              "15_QLD_WAR"             "16_NSW_SEL"              "17_NSW_NEW"              "18_NSW_WIN"              "19_NSW_CAL_2"            "2_QLD_HEL"               "20_NSW_YUR"             "21_VIC_GIP"              "3_QLD_MIT"               "4_QLD_BER_P1"            "5_QLD_BUR"               "6_QLD_MID"               "7_QLD_MIR"              "8_QLD_BIL"               "9_QLD_WID"               "ABG"                     "Barrington_River"        "Batemans_Bay_Nelligen"   "Bulga"                  "Bungonia_NP"             "Camira_SF"               "Castle_Forbes_Hermitage" "Corlette_Point"          "Crescent_Head"           "Doonside_Crescent"      "Dungog"                  "Glenugie"                "Jackdaggery"             "Lansdowne_NR"            "Lansdowne_Reserve"       "Leets_Vale"             "Leeville"                "Lowther"                 "Merringo_Beach"          "Mount_Pleasant"          "Muswellbrook"            "Nundle_SF"              "Old_Station_Creek"       "Olive_Gap"               "pop1"                    "Purrungully_Woodland"    "Putty_Rd"                "Ramornie_Rd"            "Scheyville_NP"           "South_Pambula"           "Spring_Gully"            "Tooloom"                 "Tunglebung"              "Wollombi_Brook"         "Wombeyan_Caves_Rd"       "Worrigee_NR"             "Yarravel"               
indNames(unfiltered)
## 4, over 1000
# Ensure populations are treated as categorical variables
pop(unfiltered) <- as.factor(as.character(pop(unfiltered)))

# Examine number of individuals per population
popdf <- data.frame(pop = pop(unfiltered))

popdf %>%
  count(pop)
#pop   n
#1                1_QLD_WNR  28
#2               10_QLD_BRG  25
#3               11_QLD_WON  28
#4            12_QLD_TMF_P3  21
#5               13_QLD_BRV  25
#6               15_QLD_WAR  27
#7               16_NSW_SEL  25
#8               17_NSW_NEW  28
#9               18_NSW_WIN  28
#10            19_NSW_CAL_2  13
#11               2_QLD_HEL  22
#12              20_NSW_YUR  13
#13              21_VIC_GIP  19
#14               3_QLD_MIT  26
#15            4_QLD_BER_P1  20
#16               5_QLD_BUR  26
#17               6_QLD_MID  28
#18               7_QLD_MIR  29
#19               8_QLD_BIL  24
#20               9_QLD_WID  24
#21                     ABG   6
#22        Barrington_River   6
#23   Batemans_Bay_Nelligen   6
#24                   Bulga   1
#25             Bungonia_NP   5
#26               Camira_SF   6
#27 Castle_Forbes_Hermitage   5
#28          Corlette_Point   6
#29           Crescent_Head   6
#30       Doonside_Crescent   5
#31                  Dungog  12
#32                Glenugie   6
#33             Jackdaggery   2
#34            Lansdowne_NR   6
#35       Lansdowne_Reserve   7
#36              Leets_Vale   1
#37                Leeville   6
#38                 Lowther   6
#39          Merringo_Beach   6
#40          Mount_Pleasant   6
#41            Muswellbrook  11
#42               Nundle_SF   1
#43       Old_Station_Creek   1
#44               Olive_Gap   2
#45                    pop1 595
#46    Purrungully_Woodland   6
#47                Putty_Rd   2
#48             Ramornie_Rd   6
#49           Scheyville_NP   5
#50           South_Pambula   6
#51            Spring_Gully   1
#52                 Tooloom   6
#53              Tunglebung   5
#54          Wollombi_Brook   6
#55       Wombeyan_Caves_Rd   6
#56             Worrigee_NR   5
#57                Yarravel  16

## 5 remove any inds that we know we dont want 
# (wrong species, contamination, duplicates of the same individual etc)
## In the Euc example, it is extra species we dont want
## 


inds_to_keep <- indNames(unfiltered)[unfiltered$other$ind.metrics$Genetic_ID == "Eucalyptus tereticornis"]
# inds_to_keep <- indNames(unfiltered)[unfiltered$other$ind.metrics$id == "Eucalyptus tereticornis"]
unfiltered2 <- gl.keep.ind(unfiltered, ind.list=inds_to_keep)


# unfiltered2 <- unfiltered
# 813 inds, 674 inds

popNames(unfiltered2)
nPop(unfiltered2)
# 67 pops, 56 pops


## start of basic filtering (filter 1-7) ----

## FILTER 1: REMOVE MONOMORPHIC LOCI
filter1 <- gl.filter.monomorphs(unfiltered2)
nInd(filter1)
nLoc(filter1)
# 560, 674
# 39403, 20755 

## FILTER 2: REPRODUCIBILITY
filter2 <- gl.filter.reproducibility(filter1, threshold = 0.98)
filter2
nLoc(filter2)
# 32242 with 0.98, 17055 with 0.98

filter2 <- gl.filter.reproducibility(
  filter1,
  threshold = 0.99
)
filter2
nLoc(filter2)

# 17055 with 0.98 2nd set

# 27447 loci with 0.99 (recommended as better than 1.0 by DArT, due to missing data), 14404 with 0.99 2nd set

# 23715 with 1.0, 14404 with 1.0

# Hint: always look at the automatic graph that shows the distribution of the variable you are filtering
# In this case, if there are no datapoints at the 0.99 or 0.98 mark in the graph, 
# then filtering to 0.99 wont do anything different to the undesired 1.0 filter level
# I suspect for the Euv that 0.99 or 0.98 will work okay



## FILTER 3: READ DEPTH 
# remove loci that have low coverage (unreliable for heterozygote genotyping) 
# & remove loci with high coverage (could be paralogous loci)

# 10 as lower threshold should pick up most common alleles,
#lowest threshold we are comfortable with is 8.

#explore distribution to check the upper limit 
gl.report.rdepth(filter2) # scroll up in the console to see the summary in this output
# min 2.5, 2.5
# mean 22.82405, 22.75829 
# max 362.5, 266.5 

#the automatic plot should also give some insight into rdepth distribution for your data

mean <- mean(filter2@other$loc.metrics$rdepth)
#22.82405, 22.75829
sd <- sd(filter2@other$loc.metrics$rdepth)
#26.00427, 24.6931
upper <- mean+3*sd
#100.8369, 96.8376

filter3 <- gl.filter.rdepth(filter2, lower = 8, upper = upper)
filter3
# 560  individuals, 674
# 18282 SNPs at lower = 8, 9097
# 8490 SNPs at lower = 9 - 2nd set
# 15893  SNPs at lower = 10 <- did this, 7942

## FILTER 4: CALL RATE
filter4.1 <- gl.filter.callrate(filter3, method = "loc", threshold = 0.3, recalc = TRUE)
filter4.1
# 13170, 7,522


#Filtering Callrate by Individuals
filter4.2 <- gl.filter.callrate(filter4.1, method = "ind", threshold = 0.3, recalc = T, mono.rm = TRUE)
filter4.2
# 517, 674 

#locus
filter4.3 <- gl.filter.callrate(filter4.2, method = "loc", threshold = 0.6, recalc = TRUE)
filter4.3
# 0.7 6979, 3,782
# 0.6 8330  = chosen, 4,474


# ind
filter4.4 <-  gl.filter.callrate(filter4.3, method = "ind", threshold = 0.5, recalc = TRUE, mono.rm = TRUE)
filter4.4
# 511, 670
# 8318, 4470

#locus
filter4.5 <- gl.filter.callrate(filter4.4, method = "loc", threshold = 0.8, recalc = TRUE)
filter4.5
# 0.8 5638       = chosen for now, 3092
# 0.9 3904  , 2244

# ind
filter4.6 <-  gl.filter.callrate(filter4.5, method = "ind", threshold = 0.70, recalc = TRUE, mono.rm = TRUE)
filter4.6
# @70:  = chosen
# 500, 659

#locus
filter4.7 <- gl.filter.callrate(filter4.6, method = "loc", threshold = 0.9, recalc = TRUE)
filter4.7
# 0.9 4047, 2276

# ind
filter4.8 <-  gl.filter.callrate(filter4.7, method = "ind", threshold = 0.80, recalc = TRUE, mono.rm = TRUE)
filter4.8
# 496, 653


popdf <- as.data.frame(filter4.8$pop)
names(popdf) <- c("pop")
popdf %>% 
  count(pop)

#pop  n
#1      1_QLD_WNR 28
#2     10_QLD_BRG 25
#3     11_QLD_WON 27
#4  12_QLD_TMF_P3 19
#5     13_QLD_BRV 23
#6     14_QLD_MOO 26
#7     15_QLD_WAR 27
#8     16_NSW_SEL 23
#9     17_NSW_NEW 27
#10    18_NSW_WIN 27
#11  19_NSW_CAL_2 16
#12     2_QLD_HEL 21
#13    20_NSW_YUR 13
#14    21_VIC_GIP 25
#15     3_QLD_MIT 24
#16  4_QLD_BER_P1 16
#17     5_QLD_BUR 25
#18     6_QLD_MID 28
#19     7_QLD_MIR 29
#20     8_QLD_BIL 24
#21     9_QLD_WID 23

# 2nd set pop n
# 1                1_QLD_WNR 28
#2               10_QLD_BRG 25
#3               11_QLD_WON 27
#4            12_QLD_TMF_P3 19
#5               13_QLD_BRV 23
#6               15_QLD_WAR 27
#7               16_NSW_SEL 23
#8               17_NSW_NEW 27
#9               18_NSW_WIN 26
#10            19_NSW_CAL_2 12
#11               2_QLD_HEL 21
#12              20_NSW_YUR 13
#13              21_VIC_GIP 19
#14               3_QLD_MIT 24
#15            4_QLD_BER_P1 16
#16               5_QLD_BUR 24
#17               6_QLD_MID 28
#18               7_QLD_MIR 29
#19               8_QLD_BIL 24
#20               9_QLD_WID 23
#21                     ABG  6
#22        Barrington_River  6
#23   Batemans_Bay_Nelligen  6
#24                   Bulga  1
#25             Bungonia_NP  5
#26               Camira_SF  6
#27 Castle_Forbes_Hermitage  5
#28          Corlette_Point  6
#29           Crescent_Head  6
#30       Doonside_Crescent  5
#31                  Dungog 12
#32                Glenugie  6
#33             Jackdaggery  2
#34            Lansdowne_NR  6
#35       Lansdowne_Reserve  7
#36              Leets_Vale  1
#37                Leeville  6
#38                 Lowther  6
#39          Merringo_Beach  6
#40          Mount_Pleasant  6
#41            Muswellbrook 11
#42               Nundle_SF  1
#43       Old_Station_Creek  1
#44               Olive_Gap  2
#45    Purrungully_Woodland  6
#46                Putty_Rd  2
#47             Ramornie_Rd  6
#48           Scheyville_NP  5
#49           South_Pambula  6
#50            Spring_Gully  1
#51                 Tooloom  6
#52              Tunglebung  5
#53          Wollombi_Brook  6
#54       Wombeyan_Caves_Rd  6
#55             Worrigee_NR  5
#56                Yarravel 16
filter5 <- filter4.8

################ Decided to go with call rate 0.70 for individual, 0.80 for locus
# MAC 3

## FILTER 6: MINOR ALLELE COUNT
# Note: if You give a threshold >1, it is interpreted as MAC :)
filter6 <- gl.filter.maf(filter5, threshold = 3) #=MAF 0.0159574468085106    
filter6
# 3513, 1894    

## FILTER 7: SECONDARY SNPs
filter7 <- gl.filter.secondaries(filter6, method = "best")
filter7
# 2642, 1613
#       653 

# Save filtered data ----

# in dart format
output_name <- paste0(
  "PreHWE_",
  nPop(filter7), "_pops_",
  nInd(filter7), "_inds_",
  Loc(filter7), "_SNPs.csv"
)

gl.write.csv(
  filter7,
  output_name,
  getwd()
)

# our filename convention is to include the number of pops, inds, and snps 
# This helps with data tracking to know what step it was generated at or which version throughout the course of a project

## Saving Filter Stages to RDS:----
saveRDS(filter1, file = "filter1_monomorphic_removed_V2.RDS")
saveRDS(filter2, file = "filter2_reproducibility_V2.RDS")
saveRDS(filter3, file = "filter3_readdepth_V2.RDS")
saveRDS(filter5, file = "filter5_callrate_V2.RDS")
saveRDS(filter6, file = "filter6_MAC_V2.RDS")
saveRDS(filter7, file = "filter7_secondaries_V2.RDS")


## YOUR QUEST ----
## Complete a prelim check of the data by making a PCA or pcoa plot ----

# package dartR is a good option to use, and they will have tutorial documents online
# or here: https://github.com/green-striped-gecko/dartRverse

filter7 <- filter7_secondaries_V2


pca <- gl.pcoa(filter7)

gl.pcoa.plot(
  pca,
  filter7,
  xaxis = 1,
  yaxis = 2
)
# PCA axis 1, 8.3% X
# PCA axis 2, 2.1% Y

pca_data <- as.data.frame(pca$scores)
pca_data$IND <- filter7$ind.names
pca_data$PoP <- filter7$pop

saveRDS(pca_data, file = 'pca_data.RDS')
write.csv(pca_data,  'pca_data.csv')

##################### preping
filter7 <- filter7_secondaries_V2

# Run once to install any missing packages
packages <- c(
  "ggplot2", "dplyr", "ggrepel",
  "sf", "rnaturalearth", "rnaturalearthdata"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}


# Use the final filtered genetic dataset
genetic_data <- filter7

# Check that every individual has a collection site
if (anyNA(pop(genetic_data)) ||
    any(trimws(as.character(pop(genetic_data))) == "")) {
  stop("Some individuals are missing a collection-site assignment.")
}

# Use consistent site numbers and colours across the figures
site_key <- data.frame(
  Site = sort(unique(as.character(pop(genetic_data))))
)

site_key$Site_number <- seq_len(nrow(site_key))

site_colours <- setNames(
  grDevices::hcl.colors(nrow(site_key), palette = "Dynamic"),
  site_key$Site
)

write.csv(site_key, "Site_number_key.csv", row.names = FALSE)


############# PCA

# Use the same final filtered dataset for all figures
genetic_data <- filter7

individual_info <- data.frame(
  Individual = indNames(genetic_data),
  Site = as.character(pop(genetic_data))
)

if (anyNA(individual_info$Site) ||
    any(trimws(individual_info$Site) == "")) {
  stop("Some individuals are missing a collection-site assignment.")
}

# Extract coordinates
coordinates <- genetic_data@other$latlon

# If coordinates are stored in ind.metrics instead,
# replace the line above with the following.
# Substitute your actual metadata column names:
#
# coordinates <- data.frame(
#   lat = genetic_data@other$ind.metrics$Latitude,
#   lon = genetic_data@other$ind.metrics$Longitude
# )

if (is.null(coordinates) ||
    !all(c("lat", "lon") %in% names(coordinates))) {
  stop(
    "Coordinates were not found. Check ",
    "names(genetic_data@other$ind.metrics) ",
    "and supply the correct latitude and longitude columns."
  )
}

if (nrow(coordinates) != nInd(genetic_data)) {
  stop("Coordinate rows do not match the individuals in filter7.")
}

# Coordinate rows must follow the order of individuals
# in genetic_data.
individual_info$Latitude <- suppressWarnings(
  as.numeric(as.character(coordinates$lat))
)

individual_info$Longitude <- suppressWarnings(
  as.numeric(as.character(coordinates$lon))
)

# Check latitude and longitude separately:
# a site can still colour the PCA if latitude is known
# but longitude is missing.
individual_info <- individual_info %>%
  mutate(
    Valid_latitude =
      is.finite(Latitude) &
      between(Latitude, -90, 90),
    
    Valid_coordinates =
      Valid_latitude &
      is.finite(Longitude) &
      between(Longitude, -180, 180)
  )

# Assign a consistent number to each site
site_key <- data.frame(
  Site = sort(unique(individual_info$Site))
)

site_key$Site_number <- seq_len(nrow(site_key))

# Calculate site latitude for colouring and site location
# for mapping. Each site is represented by its mean location.
site_locations <- individual_info %>%
  group_by(Site) %>%
  summarise(
    n_retained = n(),
    
    n_with_coordinates = sum(Valid_coordinates),
    
    Latitude = if (any(Valid_latitude)) {
      mean(Latitude[Valid_latitude])
    } else {
      NA_real_
    },
    
    Map_latitude = if (any(Valid_coordinates)) {
      mean(Latitude[Valid_coordinates])
    } else {
      NA_real_
    },
    
    Longitude = if (any(Valid_coordinates)) {
      mean(Longitude[Valid_coordinates])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  left_join(site_key, by = "Site")

# Report sites with missing information
missing_latitude <- site_locations %>%
  filter(is.na(Latitude))

if (nrow(missing_latitude) > 0) {
  print(missing_latitude)
  warning("Sites without latitude will appear grey in the PCA.")
}

missing_coordinates <- site_locations %>%
  filter(n_with_coordinates == 0)

if (nrow(missing_coordinates) > 0) {
  print(missing_coordinates)
  warning("Sites without valid coordinate pairs cannot be mapped.")
}

valid_latitudes <- site_locations$Latitude[
  is.finite(site_locations$Latitude)
]

if (length(valid_latitudes) == 0) {
  stop("No valid site latitudes are available.")
}

latitude_limits <- range(valid_latitudes)

# Identical gradient and limits for every figure
latitude_scale <- function() {
  scale_colour_viridis_c(
    option = "C",
    limits = latitude_limits,
    na.value = "grey70",
    name = "Site latitude (°)"
  )
}


# ============================================================
# 3. Calculate PCA
# ============================================================

pca <- gl.pcoa(genetic_data, nfactors = 3)

pca_data <- data.frame(
  Individual = indNames(genetic_data),
  Site = as.character(pop(genetic_data)),
  PC1 = pca$scores[, 1],
  PC2 = pca$scores[, 2]
) %>%
  left_join(
    site_locations %>% select(Site, Site_number, Latitude),
    by = "Site"
  )

variance_explained <- 100 * pca$eig / sum(pca$eig)

x_label <- sprintf("PC1 (%.1f%%)", variance_explained[1])
y_label <- sprintf("PC2 (%.1f%%)", variance_explained[2])


# ============================================================
# 4. Individual-level PCA coloured by site latitude
# ============================================================

pca_individuals <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, colour = Latitude)
) +
  geom_point(size = 2, alpha = 0.8) +
  latitude_scale() +
  coord_equal() +
  labs(
    title = "Genetic PCA by collection-site latitude",
    x = x_label,
    y = y_label
  ) +
  theme_bw()

print(pca_individuals)


# ============================================================
# 5. Site-average PCA coloured by latitude
# ============================================================

site_pca_data <- pca_data %>%
  group_by(Site) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    n = n(),
    .groups = "drop"
  ) %>%
  left_join(
    site_locations %>% select(Site, Site_number, Latitude),
    by = "Site"
  )

pca_sites <- ggplot(
  site_pca_data,
  aes(x = PC1, y = PC2)
) +
  geom_point(
    aes(colour = Latitude, size = n),
    alpha = 0.9
  ) +
  geom_text_repel(
    aes(label = Site_number),
    colour = "black",
    seed = 123,
    size = 3,
    max.overlaps = Inf
  ) +
  latitude_scale() +
  scale_size_continuous(range = c(3, 8)) +
  coord_equal() +
  labs(
    title = "Average genetic PCA position by collection site",
    subtitle = "Numbers match the collection-site map",
    x = x_label,
    y = y_label,
    size = "Individuals"
  ) +
  theme_bw()

print(pca_sites)


# ============================================================
# 6. Collection-site map coloured by latitude
# ============================================================

# Map the sites represented in the final filtered dataset
map_data <- filter7

# Extract coordinates
coordinates <- map_data@other$latlon

# If coordinates are stored in ind.metrics instead,
# replace the line above with this, using your actual column names:
#
# coordinates <- data.frame(
#   lat = map_data@other$ind.metrics$Latitude,
#   lon = map_data@other$ind.metrics$Longitude
# )

# Check that coordinates are available
if (is.null(coordinates) ||
    !all(c("lat", "lon") %in% names(coordinates))) {
  stop(
    "Coordinates were not found in @other$latlon. ",
    "Check names(map_data@other$ind.metrics) for coordinate columns."
  )
}

if (nrow(coordinates) != nInd(map_data)) {
  stop("Coordinate rows do not match the individuals in filter7.")
}

# Combine coordinates with individual and site information.
# Stored coordinate rows must follow the individual order.
map_individuals <- data.frame(
  Individual = indNames(map_data),
  Site = as.character(pop(map_data)),
  Latitude = suppressWarnings(
    as.numeric(as.character(coordinates$lat))
  ),
  Longitude = suppressWarnings(
    as.numeric(as.character(coordinates$lon))
  )
)

if (anyNA(map_individuals$Site) ||
    any(trimws(map_individuals$Site) == "")) {
  stop("Some individuals are missing a collection-site assignment.")
}

# Identify valid coordinates
map_individuals <- map_individuals %>%
  mutate(
    Valid_coordinates =
      is.finite(Latitude) &
      is.finite(Longitude) &
      between(Latitude, -90, 90) &
      between(Longitude, -180, 180)
  )

# Summarise the retained individuals at each site
site_locations <- map_individuals %>%
  group_by(Site) %>%
  summarise(
    n_retained = n(),
    n_with_coordinates = sum(Valid_coordinates),
    Latitude = if (any(Valid_coordinates)) {
      mean(Latitude[Valid_coordinates])
    } else {
      NA_real_
    },
    Longitude = if (any(Valid_coordinates)) {
      mean(Longitude[Valid_coordinates])
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  left_join(site_key, by = "Site")

# Check that site numbering matches the PCA
if (anyNA(site_locations$Site_number)) {
  stop("Re-run the PCA setup using filter7 to update site_key.")
}

# Report sites that cannot be mapped
missing_sites <- site_locations %>%
  filter(n_with_coordinates == 0)

if (nrow(missing_sites) > 0) {
  print(missing_sites)
  warning("These sites need coordinates before they can be mapped.")
}

# Keep sites with usable coordinates for plotting
map_sites <- site_locations %>%
  filter(n_with_coordinates > 0)

if (nrow(map_sites) == 0) {
  stop("No valid coordinates are available.")
}

# Load the Australian basemap
australia <- rnaturalearth::ne_countries(
  country = "Australia",
  scale = "medium",
  returnclass = "sf"
)

# Set map boundaries around the sampling region
longitude_limits <- range(map_sites$Longitude) + c(-1, 1)
latitude_limits <- range(map_sites$Latitude) + c(-1, 1)

# Plot one point per site at its mean recorded location
site_map <- ggplot() +
  geom_sf(
    data = australia,
    fill = "grey95",
    colour = "grey55",
    linewidth = 0.3
  ) +
  geom_point(
    data = map_sites,
    aes(
      x = Longitude,
      y = Latitude,
      colour = Latitude
    ),
    size = 3
  ) +
  geom_text_repel(
    data = map_sites,
    aes(
      x = Longitude,
      y = Latitude,
      label = Site_number
    ),
    colour = "black",
    seed = 123,
    size = 3,
    max.overlaps = Inf
  ) +
  scale_colour_viridis_c(
    option = "C",
    name = "Site latitude (°)"
  ) +
  coord_sf(
    xlim = longitude_limits,
    ylim = latitude_limits,
    expand = FALSE
  ) +
  labs(
    title = "Eucalyptus tereticornis collection sites",
    subtitle = paste(
      "Sites retained after filtering;",
      "numbers match the Site number key"
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(site_map)

# Save the map and coordinate summary
ggsave(
  "Collection_site_map_filter7.png",
  plot = site_map,
  width = 10,
  height = 10,
  dpi = 300
)

write.csv(
  site_locations,
  "Collection_site_coordinates_filter7.csv",
  row.names = FALSE
)
# ============================================================
# 7. Save PCA figures and supporting tables
# ============================================================

ggsave(
  "PCA_individuals_by_latitude.png",
  plot = pca_individuals,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  "PCA_site_averages_by_latitude.png",
  plot = pca_sites,
  width = 10,
  height = 8,
  dpi = 300
)

write.csv(
  pca_data,
  "PCA_individual_scores_by_latitude.csv",
  row.names = FALSE
)

write.csv(
  site_pca_data,
  "PCA_site_averages_by_latitude.csv",
  row.names = FALSE
)

write.csv(
  site_locations,
  "Collection_site_coordinates_filter7.csv",
  row.names = FALSE
)

write.csv(
  site_key,
  "Site_number_key.csv",
  row.names = FALSE
)



##################################################

# ============================================================
# 8. Analysis settings and packages
# ============================================================
# Install once, if needed:
# install.packages(c("dartR.base", "dartR.spatial", "dartR.popgen",
#                    "dismo", "vegan", "geosphere", "ggplot2", "dplyr",
#                    "ggrepel", "sf", "rnaturalearth", "rnaturalearthdata"))
required <- c("dartR.base", "dartR.spatial", "dartR.popgen", "dismo",
              "vegan", "geosphere", "ggplot2", "dplyr", "ggrepel",
              "sf", "rnaturalearth", "rnaturalearthdata")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
library(dplyr)
library(ggplot2)
library(ggrepel)

out_dir <- "Spatial_results_noSite9"
dir.create(out_dir, showWarnings = FALSE)
out <- function(name) file.path(out_dir, name)
bin_width_km <- 100       # CHANGE THIS to adjust spatial distance classes
mantel_permutations <- 9999
spatial_reps <- 999       # Both permutations and bootstrap replicates
minimum_site_n <- 2      # Single-individual sites omitted from FST only
map_axis <- "PC1"         # Change to "PC2" to colour the map by PC2 scores
stopifnot(is.finite(bin_width_km), bin_width_km > 0,
          map_axis %in% c("PC1", "PC2"))



excluded_site <- "18_NSW_WIN"

# Remove every individual belonging to site 9
ids_to_remove <- adegenet::indNames(filter7)[
  as.character(adegenet::pop(filter7)) == excluded_site
]

if (length(ids_to_remove) > 0) {
  filter7 <- dartR.base::gl.drop.ind(
    filter7,
    ind.list = ids_to_remove
  )
}

# Update the object used for subsequent analyses
genetic_data <- filter7

# Confirm removal: should return FALSE
excluded_site %in% as.character(adegenet::pop(genetic_data))





genetic_data <- filter7
stopifnot(!excluded_site %in% as.character(pop(genetic_data)))
if (anyDuplicated(indNames(genetic_data))) stop("Individual IDs must be unique.")
individual_info <- data.frame(
  Individual = indNames(genetic_data), Site = as.character(pop(genetic_data)))
if (anyNA(individual_info$Site) || any(trimws(individual_info$Site) == "")) {
  stop("Some individuals have no site assignment.")
}
coordinates <- genetic_data@other$latlon
# If needed, replace the preceding line with your actual metadata columns:
# coordinates <- data.frame(lat = genetic_data@other$ind.metrics$Latitude,
#                           lon = genetic_data@other$ind.metrics$Longitude)
if (is.null(coordinates) || !all(c("lat", "lon") %in% names(coordinates)) ||
    nrow(coordinates) != nInd(genetic_data)) stop("Check coordinate columns and row alignment.")
# Match meaningful row names when available; otherwise assume genlight order,
# as in the original script. Never sort coordinate rows independently.
if (setequal(rownames(coordinates), individual_info$Individual)) {
  coordinates <- coordinates[individual_info$Individual, , drop = FALSE]
}
individual_info$Latitude <- suppressWarnings(as.numeric(as.character(coordinates$lat)))
individual_info$Longitude <- suppressWarnings(as.numeric(as.character(coordinates$lon)))
individual_info <- individual_info %>% mutate(
  Valid_latitude = is.finite(Latitude) & between(Latitude, -90, 90),
  Valid_coordinates = Valid_latitude & is.finite(Longitude) & between(Longitude, -180, 180))
write.csv(individual_info, out("Individual_coordinate_audit.csv"), row.names = FALSE)

# Retain original alphabetical map numbering, leaving a gap for the removed
# site. These labels are NOT the numeric prefix in the population name.
site_key <- data.frame(Site = sort(unique(c(individual_info$Site, excluded_site))))
site_key$Site_number <- seq_len(nrow(site_key))
site_key <- site_key[site_key$Site != excluded_site, ]
site_locations <- individual_info %>% group_by(Site) %>% summarise(
  n_retained = n(), n_with_coordinates = sum(Valid_coordinates),
  Latitude = if (any(Valid_latitude)) mean(Latitude[Valid_latitude]) else NA_real_,
  Map_latitude = if (any(Valid_coordinates)) mean(Latitude[Valid_coordinates]) else NA_real_,
  Longitude = if (any(Valid_coordinates)) mean(Longitude[Valid_coordinates]) else NA_real_,
  .groups = "drop") %>% left_join(site_key, by = "Site")
if (any(!individual_info$Valid_coordinates)) {
  warning("Missing/invalid coordinates: see audit. These individuals remain in PCA; spatial analyses exclude them as specified below.")
}
# ============================================================
# 9. Recompute PCA after excluding site 9
# ============================================================
pca <- dartR.base::gl.pcoa(genetic_data, nfactors = 3)
stopifnot(nrow(pca$scores) == nInd(genetic_data))
pca_data <- data.frame(Individual = indNames(genetic_data),
                       Site = as.character(pop(genetic_data)), PC1 = pca$scores[, 1], PC2 = pca$scores[, 2]) %>%
  left_join(site_locations %>% select(Site, Site_number, Latitude), by = "Site")
# Eigenvalues describe axes; individual/site scores describe their positions.
variance_explained <- 100 * pca$eig / sum(pca$eig)
x_label <- sprintf("PC1 (%.1f%%)", variance_explained[1])
y_label <- sprintf("PC2 (%.1f%%)", variance_explained[2])
site_pca_data <- pca_data %>% group_by(Site) %>% summarise(
  PC1 = mean(PC1), PC2 = mean(PC2), n = n(), .groups = "drop") %>%
  left_join(site_locations, by = "Site")
latitude_scale <- function() scale_colour_viridis_c(option = "C",
                                                    name = "Site latitude (°)", na.value = "grey70")
pca_individuals <- ggplot(pca_data, aes(PC1, PC2, colour = Latitude)) +
  geom_point(size = 2, alpha = 0.8) + latitude_scale() + coord_equal() +
  labs(title = "Genetic PCA by collection-site latitude", x = x_label, y = y_label) + theme_bw()
pca_sites <- ggplot(site_pca_data, aes(PC1, PC2)) +
  geom_point(aes(colour = Latitude, size = n)) +
  geom_text_repel(aes(label = Site_number), seed = 123, max.overlaps = Inf) +
  latitude_scale() + coord_equal() + theme_bw() +
  labs(title = "Site-average PCA", x = x_label, y = y_label, size = "Individuals")
print(pca_individuals)
print(pca_sites)
ggsave(out("PCA_individuals_by_latitude.png"), pca_individuals, width = 10, height = 8, dpi = 300)
ggsave(out("PCA_sites_by_latitude.png"), pca_sites, width = 10, height = 8, dpi = 300)

# ============================================================
# 10. Map coloured by mean site PC1 (or PC2) SCORE
# ============================================================
map_sites <- site_pca_data %>% filter(n_with_coordinates > 0)
# Extract coordinates
coordinates <- map_data@other$latlon

# If coordinates are stored in ind.metrics instead,
# replace the line above with this, using your actual column names:
#
# coordinates <- data.frame(
#   lat = map_data@other$ind.metrics$Latitude,
#   lon = map_data@other$ind.metrics$Longitude
# )

# Check that coordinates are available
if (is.null(coordinates) ||
    !all(c("lat", "lon") %in% names(coordinates))) {
  stop(
    "Coordinates were not found in @other$latlon. ",
    "Check names(map_data@other$ind.metrics) for coordinate columns."
  )
}

if (nrow(coordinates) != nInd(map_data)) {
  stop("Coordinate rows do not match the individuals in filter7.")
}

# Combine coordinates with individual and site information.
# Stored coordinate rows must follow the individual order.
map_individuals <- data.frame(
  Individual = indNames(map_data),
  Site = as.character(pop(map_data)),
  Latitude = suppressWarnings(
    as.numeric(as.character(coordinates$lat))
  ),
  Longitude = suppressWarnings(
    as.numeric(as.character(coordinates$lon))
  )
)

if (anyNA(map_individuals$Site) ||
    any(trimws(map_individuals$Site) == "")) {
  stop("Some individuals are missing a collection-site assignment.")
}

# Identify valid coordinates
map_individuals <- map_individuals %>%
  mutate(
    Valid_coordinates =
      is.finite(Latitude) &
      is.finite(Longitude) &
      between(Latitude, -90, 90) &
      between(Longitude, -180, 180)
  )

# Summarise the retained individuals at each site
site_locations <- map_individuals %>%
  group_by(Site) %>%
  summarise(
    n_retained = n(),
    n_with_coordinates = sum(Valid_coordinates),
    Latitude = if (any(Valid_coordinates)) {
      mean(Latitude[Valid_coordinates])
    } else {
      NA_real_
    },
    Longitude = if (any(Valid_coordinates)) {
      mean(Longitude[Valid_coordinates])
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  left_join(site_key, by = "Site")

# Check that site numbering matches the PCA
if (anyNA(site_locations$Site_number)) {
  stop("Re-run the PCA setup using filter7 to update site_key.")
}

# Report sites that cannot be mapped
missing_sites <- site_locations %>%
  filter(n_with_coordinates == 0)

if (nrow(missing_sites) > 0) {
  print(missing_sites)
  warning("These sites need coordinates before they can be mapped.")
}

# Keep sites with usable coordinates for plotting
map_sites <- site_locations %>%
  filter(n_with_coordinates > 0)

if (nrow(map_sites) == 0) {
  stop("No valid coordinates are available.")
}
# Recalculate PCA scores using the same dataset as the map
map_pca <- dartR.base::gl.pcoa(map_data, nfactors = 3)

map_scores <- data.frame(
  Site = as.character(adegenet::pop(map_data)),
  PC1 = map_pca$scores[, 1]
)

# Calculate one mean PC1 score per site
map_site_scores <- map_scores %>%
  group_by(Site) %>%
  summarise(
    Mean_PC1 = mean(PC1),
    .groups = "drop"
  )

# Preserve ALL map sites and attach their PC1 scores
map_sites <- map_sites %>%
  left_join(map_site_scores, by = "Site")

# Report how many sites are included
message(
  "Mapping ", nrow(map_sites), " sites out of ",
  nrow(site_locations), " retained sites."
)

# Missing scores remain visible as grey points
if (any(!is.finite(map_sites$Mean_PC1))) {
  warning("Some sites have missing PC1 scores and will appear grey.")
}
# Load the Australian basemap
australia <- rnaturalearth::ne_countries(
  country = "Australia",
  scale = "medium",
  returnclass = "sf"
)

# Set map boundaries around the sampling region
longitude_limits <- range(map_sites$Longitude) + c(-1, 1)
latitude_limits <- range(map_sites$Latitude) + c(-1, 1)

# Plot one point per site at its mean recorded location
site_map <- ggplot() +
  geom_sf(
    data = australia,
    fill = "grey95",
    colour = "grey55",
    linewidth = 0.3
  ) +
  geom_point(
    data = map_sites,
    aes(
      x = Longitude,
      y = Latitude,
      colour = Mean_PC1
    ),
    size = 3
  ) +
  geom_text_repel(
    data = map_sites,
    aes(
      x = Longitude,
      y = Latitude,
      label = Site_number
    ),
    colour = "black",
    seed = 123,
    size = 3,
    max.overlaps = Inf
  ) +
  scale_colour_viridis_c(option = "C", name = paste("Mean site", map_axis, "score")) +
  coord_sf(
    xlim = longitude_limits,
    ylim = latitude_limits,
    expand = FALSE
  ) +
  labs(
    title = "Collection sites coloured by mean PC1 scores",
    subtitle = paste(
      "Sites retained after filtering;",
      "numbers match the Site number key"
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

print(site_map)

ggsave(out("Collection_site_map_by_PC_score.png"), site_map, width = 10, height = 8, dpi = 300)
write.csv(site_key, out("Site_number_key.csv"), row.names = FALSE)
write.csv(pca_data, out("PCA_individual_scores.csv"), row.names = FALSE)
write.csv(site_pca_data, out("PCA_site_scores_and_coordinates.csv"), row.names = FALSE)

# Great-circle distances on the WGS84 ellipsoid, in METRES.
geo_matrix <- function(lon, lat, ids) {
  m <- geosphere::distm(cbind(lon, lat), fun = geosphere::distGeo)
  dimnames(m) <- list(ids, ids)
  diag(m) <- 0
  m
}

# ============================================================
# 11. SITE-LEVEL standard Mantel test: FST vs distance (km)
# ============================================================
# Retain sites with >=2 individuals AND a recorded site location.
# A site centroid uses only members with valid coordinate pairs; genotype
# estimates use all retained members of that site, including missing GPS.
site_audit <- site_locations %>% mutate(
  Included_Mantel = n_retained >= minimum_site_n & n_with_coordinates > 0)
write.csv(site_audit, out("Mantel_site_inclusion_audit.csv"), row.names = FALSE)
print(site_audit %>% filter(!Included_Mantel))
# Build the expected site list from eligible retained sites
ibd_sites <- as.character(
  site_audit$Site[site_audit$Included_Mantel]
)

# Only exclude the intended NSW site
excluded_site <- "18_NSW_WIN"
ibd_sites <- setdiff(ibd_sites, excluded_site)

# Rebuild the genetic dataset and FST matrix
ibd_data <- dartR.base::gl.keep.pop(
  genetic_data,
  pop.list = ibd_sites
)

fst_dist <- dartR.base::gl.fst.pop(
  ibd_data,
  nboots = 1
)

fst_matrix <- as.matrix(fst_dist)
site_order <- rownames(fst_matrix)

# Verify matching labels
stopifnot(
  !is.null(site_order),
  setequal(site_order, ibd_sites),
  !excluded_site %in% site_order
)
ibd_locations <- site_locations[match(site_order, site_locations$Site), ]
geographic_km <- geo_matrix(
  lon = ibd_locations$Longitude,
  lat = ibd_locations$Latitude,
  ids = site_order
) / 1000
stopifnot(identical(rownames(fst_matrix), rownames(geographic_km)))
if (any(!is.finite(fst_matrix[lower.tri(fst_matrix)]))) {
  stop("Undefined pairwise FST: inspect site sample sizes and shared genotyped loci.")
}
# Keep negative FST estimates as estimated; do not silently truncate them.
if (any(fst_matrix[lower.tri(fst_matrix)] < 0)) message("Negative FST estimates retained.")
if (sd(as.vector(as.dist(fst_matrix))) == 0 || sd(as.vector(as.dist(geographic_km))) == 0) {
  stop("Mantel requires variation in both distance matrices.")
}
set.seed(123)
mantel_result <- vegan::mantel(as.dist(fst_matrix), as.dist(geographic_km),
                               method = "pearson", permutations = mantel_permutations)
print(mantel_result)
# vegan tests a positive association (one-sided), appropriate to the stated IBD hypothesis.
write.csv(data.frame(r = unname(mantel_result$statistic), p = mantel_result$signif,
                     permutations = mantel_result$permutations, sites = length(site_order),
                     genetic_distance = "Pairwise FST", geographic_distance = "Untransformed km",
                     alternative = "Positive association"), out("Site_Mantel_summary.csv"), row.names = FALSE)
write.csv(fst_matrix, out("Site_pairwise_FST.csv"))
write.csv(geographic_km, out("Site_geographic_distance_km.csv"))
idx <- which(lower.tri(fst_matrix), arr.ind = TRUE)
ibd_pairs <- data.frame(Site1 = site_order[idx[, 1]], Site2 = site_order[idx[, 2]],
                        Distance_km = geographic_km[idx], FST = fst_matrix[idx])
write.csv(ibd_pairs, out("Site_IBD_pairs.csv"), row.names = FALSE)
ibd_plot <- ggplot(ibd_pairs, aes(Distance_km, FST)) + geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) + theme_bw() +
  labs(title = "Site-level isolation by distance",
       subtitle = sprintf("Mantel r = %.3f; permutation p = %.4g", mantel_result$statistic, mantel_result$signif),
       x = "Geographic distance (km)", y = expression(F[ST]))
# Line is descriptive only: pairwise points are not independent replicates.
print(ibd_plot)
ggsave(out("Site_Mantel_IBD.png"), ibd_plot, width = 9, height = 6, dpi = 300)
saveRDS(mantel_result, out("Site_Mantel_result.rds"))

# ============================================================
# 12. INDIVIDUAL spatial autocorrelation, 100 km classes
# ============================================================
# Use all individuals with valid GPS, including single-individual sites.
# Supply distance matrices to analyse ACROSS sites: passing the original
# genlight to this function would split the analysis by population.
valid_ids <- individual_info$Individual[individual_info$Valid_coordinates]
if (length(valid_ids) < 3) stop("Too few individuals with valid GPS.")
spatial_data <- dartR.base::gl.keep.ind(genetic_data, ind.list = valid_ids)
spatial_info <- individual_info[match(indNames(spatial_data), individual_info$Individual), ]
stopifnot(identical(indNames(spatial_data), spatial_info$Individual))
Dgeo_m <- geo_matrix(spatial_info$Longitude, spatial_info$Latitude, spatial_info$Individual)
Dgen <- as.matrix(dartR.base::gl.dist.ind(spatial_data, method = "Euclidean", verbose = 0))
if (is.null(rownames(Dgen)) || !setequal(rownames(Dgen), spatial_info$Individual)) {
  stop("Genetic distance labels do not match individual IDs.")
}
Dgen <- Dgen[spatial_info$Individual, spatial_info$Individual]
if (any(!is.finite(Dgen))) stop("Undefined genetic distances: check shared nonmissing loci.")
# Match dartR's genlight/Euclidean branch, which uses 1 - distance before
# autocorrelation. Do not feed allele-sharing similarity into this directly.
Dgen_for_sac <- 1 - Dgen
diag(Dgen_for_sac) <- 0
max_distance_m <- max(Dgeo_m)
if (max_distance_m <= 0) stop("Spatial autocorrelation needs distinct locations.")



bin_width_m <- bin_width_km * 1000



# Extend beyond the maximum so every pair is covered, even on a boundary.
breaks_m <- seq(0, (floor(max_distance_m / bin_width_m) + 1) * bin_width_m, by = bin_width_m)
# Negative infinitesimal first boundary includes distinct individuals at
# identical coordinates even in versions using an open lower bin boundary.
breaks_m[1] <- -1e-8
write.csv(data.frame(Lower_km = pmax(0, head(breaks_m, -1)) / 1000,
                     Upper_km = tail(breaks_m, -1) / 1000),
          out("Spatial_distance_classes_km.csv"), row.names = FALSE)
message("Distinct individuals sharing coordinates: ", sum(Dgeo_m[lower.tri(Dgeo_m)] == 0), " pairs.")
set.seed(123)
spatial_result <- dartR.spatial::gl.spatial.autoCorr(
  x = NULL, Dgen = as.dist(Dgen_for_sac), Dgeo = as.dist(Dgeo_m),
  bins = breaks_m, reps = spatial_reps, permutation = TRUE, bootstrap = TRUE,
  plot.out = FALSE, verbose = 2)
# dartR versions return either a list of per-population tables or a table.
if (is.data.frame(spatial_result)) {
  spatial_table <- spatial_result
} else if (is.list(spatial_result) && length(spatial_result) == 1 &&
           is.data.frame(spatial_result[[1]])) {
  spatial_table <- spatial_result[[1]]
} else stop("Unexpected autocorrelation result; inspect spatial_result and package version.")
stopifnot(nrow(spatial_table) == length(breaks_m) - 1)
# Explicit class bounds avoid relying on dartR's plot-axis unit formatting.
spatial_table$Lower_km <- pmax(0, head(breaks_m, -1)) / 1000
spatial_table$Upper_km <- tail(breaks_m, -1) / 1000
spatial_table$Midpoint_km <- (spatial_table$Lower_km + spatial_table$Upper_km) / 2
# Per-class p values are exploratory; also supply a BH-adjusted column.
spatial_table$p_BH <- p.adjust(spatial_table$p.one.tail, method = "BH")
print(spatial_table)
if (any(spatial_table$N == 0)) warning("Empty distance classes: consider a wider bin_width_km.")
write.csv(spatial_table, out("Individual_spatial_autocorrelation.csv"), row.names = FALSE)
saveRDS(spatial_result, out("Individual_spatial_autocorrelation.rds"))
sac_plot <- ggplot(spatial_table, aes(Midpoint_km, r)) +
  geom_ribbon(aes(ymin = L.r.null, ymax = U.r.null), fill = "grey70", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  geom_errorbar(aes(ymin = L.r, ymax = U.r), width = bin_width_km * 0.15, na.rm = TRUE) +
  theme_bw() + labs(title = "Individual-level genetic spatial autocorrelation",
                    subtitle = paste(bin_width_km, "km classes; ribbon: 95% null envelope; bars: 95% bootstrap intervals"))
print(sac_plot)

