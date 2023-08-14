suppressPackageStartupMessages({
  library(maptools)
  library(sp)
  library(robustbase)
  library(Rcpp)
  library(spatialreg)
  library(spData)
  library(Matrix)
  library(sf)
  library(GWmodel)
  library(dplyr)
  library(readr)
  library(forcats)
  library(stringr)
  library(lubridate)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(ranger)
  library(caret)
  library(lattice)
  library(h2o)
  library(pdp)
  library(vip)
  library(plyr)         ## Data management
  library(sp)           ## Spatial Data management
  library(spdep)        ## Spatial autocorrelation
  library(RColorBrewer) ## Visualization
  library(classInt)     ## Class intervals
  library(raster)       ## spatial data
  library(grid)         ## plot
  library(gridExtra)    ## Multiple plot
  library(ggplot2)      #  plotting
  library(tidyverse)    # data 
  library(SpatialML)    # Geographically weigted regression
  library(purrr)
  library(h2o)
  library(geojsonio)
  library(pdp)       # for partial dependence plots (and ICE curves)
  library(vip)       # for variable importance plots
  library(matrixStats)
})

library(readxl)
df <- read_csv(file_path,show_col_types = FALSE)

set.seed(123)
split <- sample(c("train", "test"), size = nrow(df), prob = c(0.8, 0.2), replace = TRUE)
df$split <- split

test.df<-df %>% 
  dplyr::select(FID_grid, x, y, Income, AreaProx,LUP_Other,Resi_lur,Edu_lur,Recre_lur,Public_lur,Health_lur,Other_lur,Comm_lur,
                Rd_Pri,Rd_Cyc,Rd_Ter,Rd_Pedes,T_Enable,T_Hinder,Rd_Sec,Weekend,Summer,Winter,Spring, Dist_center, trip_count) %>%
  filter(split == 'test')

valid.df<-df %>% 
  dplyr::select(FID_grid, x, y, Income, AreaProx,LUP_Other,Resi_lur,Edu_lur,Recre_lur,Public_lur,Health_lur,Other_lur,Comm_lur,
                Rd_Pri,Rd_Cyc,Rd_Ter,Rd_Pedes,T_Enable,T_Hinder,Rd_Sec,Weekend,Summer,Winter,Spring, Dist_center, trip_count) %>%
  filter(split == 'test')

train.df<-df %>% 
  dplyr::select(FID_grid, x, y, Income, AreaProx,LUP_Other,Resi_lur,Edu_lur,Recre_lur,Public_lur,Health_lur,Other_lur,Comm_lur,
                Rd_Pri,Rd_Cyc,Rd_Ter,Rd_Pedes,T_Enable,T_Hinder,Rd_Sec,Weekend,Summer,Winter,Spring, Dist_center, trip_count) %>%
  filter(split == 'train')

# Scale covariates
cols_sel <- c('Income','AreaProx','LUP_Other','Resi_lur','Edu_lur','Recre_lur','Public_lur','Health_lur','Other_lur','Comm_lur',
              'Rd_Pri','Rd_Cyc','Rd_Ter','Rd_Pedes','T_Enable','T_Hinder','Rd_Sec','Weekend','Summer','Winter','Spring',
              'Dist_center')

test.df[, cols_sel] = scale(test.df[, cols_sel])
valid.df[, cols_sel] = scale(valid.df[, cols_sel])
train.df[, cols_sel] = scale(train.df[, cols_sel])

h2o.init(nthreads = -1,max_mem_size ="48g",enable_assertions = FALSE)

cols_sel1 <- c('Income','AreaProx','LUP_Other','Resi_lur','Edu_lur','Recre_lur','Public_lur','Health_lur','Other_lur',
               'Comm_lur','Rd_Pri','Rd_Cyc','Rd_Ter','Rd_Pedes','T_Enable','T_Hinder','Rd_Sec','Weekend','Summer','Winter',
               'Spring','Dist_center','trip_count')

test.mf<-test.df[, cols_sel1] 
valid.mf<-valid.df[, cols_sel1]
train.mf<-train.df[, cols_sel1]

test.hex<-  as.h2o(test.mf)
valid.hex<-  as.h2o(valid.mf)
train.hex<-  as.h2o(train.mf)

# Local Random Forest
features <- as.data.frame(train.hex) %>%  dplyr::select(-trip_count)
response <- as.data.frame(train.hex) %>% pull(trip_count)

pred <- function(object, newdata)  {
  results <- as.vector(h2o.predict(object, as.h2o(newdata), local.w=1, global.w=0)) #weighted global as 0
  return(results)
}

# Optimal bandwidth
combined.df <- rbind(train.df, test.df)
Coords<-combined.df[ ,2:3]


# Optimality for modelling
results <- rf.mtry.optim(trip_count ~ Income+AreaProx+LUP_Other+Resi_lur+Edu_lur+Recre_lur+Public_lur+Health_lur+
                           Other_lur+Comm_lur+Rd_Pri+Rd_Cyc+Rd_Ter+Rd_Pedes+T_Enable+T_Hinder+Rd_Sec+Weekend+
                           Summer+Winter+Spring+Dist_center, combined.df)

# Set optimal values
mtry = 17
kernel = "adaptive"   # or fixed

bandwidth <- grf.bw(trip_count ~ Income+AreaProx+LUP_Other+Resi_lur+Edu_lur+Recre_lur+Public_lur+Health_lur+
                      Other_lur+Comm_lur+Rd_Pri+Rd_Cyc+Rd_Ter+Rd_Pedes+T_Enable+T_Hinder+Rd_Sec+Weekend+Summer+
                      Winter+Spring+Dist_center,
                    dataset = combined.df,
                    kernel="adaptive", coords=Coords, bw.min = 153, bw.max = 156, step = 1,
                    trees=500, mtry=mtry, importance="impurity", nthreads = 1, forests = TRUE, weighted = TRUE)

# Set optimal values
bw = bandwidth$Best.BW

# Generate GRF model
grf.model <- grf(trip_count ~ Income+AreaProx+LUP_Other+Resi_lur+Edu_lur+Recre_lur+Public_lur+Health_lur+
                   Other_lur+Comm_lur+Rd_Pri+Rd_Cyc+Rd_Ter+Rd_Pedes+T_Enable+T_Hinder+Rd_Sec+Weekend+
                   Summer+Winter+Spring+Dist_center,
                 bw=bw, importance = "impurity", nthreads = 1, dframe=combined.df,
                 ntree=120, mtry=2, kernel=kernel, forests = TRUE, coords=Coords, weighted=TRUE)

# Fetch the FID_grid and Local_R2 values from grf.model
df2 <- data.frame(FID_grid = combined.df$FID_grid, 
                  Coord = grf.model$Locations, 
                  Local_R2 = grf.model$LGofFit$LM_Rsq100, 
                  grf_Inc = grf.model$Local.Variable.Importance$Income, 
                  AreaProx = grf.model$Local.Variable.Importance$AreaProx, 
                  LUP_Other = grf.model$Local.Variable.Importance$LUP_Other,
                  Resi_lur = grf.model$Local.Variable.Importance$Resi_lur,
                  Edu_lur = grf.model$Local.Variable.Importance$Edu_lur,
                  Recre_lur = grf.model$Local.Variable.Importance$Recre_lur,
                  Public_lur = grf.model$Local.Variable.Importance$Public_lur,
                  Health_lur = grf.model$Local.Variable.Importance$Health_lur,
                  Other_lur = grf.model$Local.Variable.Importance$Other_lur,
                  Comm_lur = grf.model$Local.Variable.Importance$Comm_lur,
                  Rd_Pri = grf.model$Local.Variable.Importance$Rd_Pri,
                  Rd_Cyc = grf.model$Local.Variable.Importance$Rd_Cyc,
                  Rd_Ter = grf.model$Local.Variable.Importance$Rd_Ter,
                  Rd_Pedes = grf.model$Local.Variable.Importance$Rd_Pedes,
                  T_Enable = grf.model$Local.Variable.Importance$T_Enable,
                  T_Hinder = grf.model$Local.Variable.Importance$T_Hinder,
                  Rd_Sec = grf.model$Local.Variable.Importance$Rd_Sec,
                  Weekend = grf.model$Local.Variable.Importance$Weekend,
                  Summer = grf.model$Local.Variable.Importance$Summer,
                  Winter = grf.model$Local.Variable.Importance$Winter,
                  Spring = grf.model$Local.Variable.Importance$Spring,
                  Dist_center = grf.model$Local.Variable.Importance$Dist_center)

write.csv(df2,file=write_file)
