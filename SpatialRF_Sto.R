suppressPackageStartupMessages({
  library(spatialRF)
  library(kableExtra)
  library(tidyverse)
  library(randomForestExplainer)
  library(pdp)
  library(readr)
  library(dplyr)
})

df <- read_csv(file_path1, show_col_types = FALSE)
dfx <- data.frame(ColumnNames = c("FID_grid", "x", "y", "Comm_lur", "Income", "T_Enable", "Rd_Pedes", "Rd_Ter", "T_Hinder", "Resi_lur", "Health_lur", 
                                  "AreaProx", "Other_lur", "Edu_lur", "LUP_Other", "Recre_lur", "Rd_Cyc", "Public_lur", "Rd_Sec", "Rd_Pri",
                                  "Dist_center", "Summer", "Winter", "Spring", "Weekend", "trip_count"))
selected_columns <- dfx$ColumnNames
dfx <- as.data.frame(df[selected_columns])
dfx$FID_grid <- as.integer(dfx$FID_grid)

# set distance matrix
distance_matrix <- read.csv(file_path2)
distance_matrix <- as.matrix(distance_matrix)
colnames(distance_matrix) <- gsub("X", "", colnames(distance_matrix))
rownames(distance_matrix) <- distance_matrix[, "Orig"]
distance_matrix <- subset(distance_matrix, select = -c(Orig))
distance_matrix <- distance_matrix[1:(nrow(distance_matrix) - 1), 1:(ncol(distance_matrix) - 1)]
distance.matrix <- as.numeric(as.matrix(distance_matrix))
distance.matrix <- as.matrix(distance.matrix)
new_dim <- c(2990, 2990)
distance.matrix <- matrix(distance.matrix, nrow = new_dim[1], ncol = new_dim[2], byrow = TRUE)
new_row_names <- rownames(distance_matrix)
new_col_names <- colnames(distance_matrix)
rownames(distance.matrix) <- new_row_names
colnames(distance.matrix) <- new_col_names


dependent.variable.name <- "trip_count"
predictor.variable.names <- colnames(dfx)[4:25]
xy <- dfx[, c("x", "y")]
distance.thresholds <- c(50, 100, 200, 400, 700, 1000, 2000, 5000, 8000)
random.seed <- 123

# Get Morans-I figure
spatialRF::plot_training_df_moran(data = dfx,
                                 dependent.variable.name = dependent.variable.name,
                                 predictor.variable.names = predictor.variable.names,
                                 distance.matrix = distance.matrix,
                                 distance.thresholds = distance.thresholds,
                                 fill.color = viridis::viridis(100, option = "F", direction = -1),
                                 point.color = "gray40")

#Infl features
spatialRF::plot_training_df(
   data = dfx,
   dependent.variable.name = dependent.variable.name,
   predictor.variable.names = predictor.variable.names,
   ncol = 6, point.color = viridis::viridis(100, option = "F"), 
   line.color = "gray30")

# Finding promising variable interactions
interactions <- spatialRF::the_feature_engineer(
  data = dfx,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  xy = xy,
  importance.threshold = 0.65, #uses 40% best predictors
  cor.threshold = 0.70, #max corr between interactions and predictors
  seed = 123,
  repetitions = 100,
  verbose = TRUE
)

# Reducing multicollinearity in the predictors
preference.order <- c("Income","Resi_lur","Comm_lur","Dist_center")

predictor.variable.names <- spatialRF::auto_cor(
  x = dfx[, predictor.variable.names],
  cor.threshold = 0.65,
  preference.order = preference.order
) %>% 
  spatialRF::auto_vif(
    vif.threshold = 10,
    preference.order = preference.order
  )

dfx <- interactions$data
predictor.variable.names <- interactions$predictor.variable.names

# Fitting a non-spatial Random Forest model with rf()
model.non.spatial <- spatialRF::rf(
  data = dfx,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy,
  seed = 123,
  verbose = FALSE
)

#Residuals
spatialRF::plot_residuals_diagnostics(
  model.non.spatial,
  verbose = FALSE
)

# Global variable importance
spatialRF::plot_importance(
  model.non.spatial,
  verbose = FALSE
)

importance.dfx <- randomForestExplainer::measure_importance(
  model.non.spatial,
  measures = c("mean_min_depth", "no_of_nodes", "times_a_root", "p_value")
)

# Contribution of predictors to model transferability
model.non.spatial <- spatialRF::rf_importance(
  model = model.non.spatial
)

model.non.spatial$importance$per.variable %>% 
  ggplot2::ggplot() +
  ggplot2::aes(
    x = importance.oob,
    y = importance.cv
  ) + 
  ggplot2::geom_point(size = 3) + 
  ggplot2::theme_bw() +
  ggplot2::xlab("Importance (out-of-bag)") + 
  ggplot2::ylab("Contribution to transferability") + 
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "red4")

# Local variable importance
local.importance <- spatialRF::get_importance_local(model.non.spatial)

kableExtra::kbl(
  round(local.importance[1:10, 1:5], 0),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)

#adding coordinates
local.importance <- cbind(
  xy,
  local.importance
)

# Response curves and surfaces
spatialRF::plot_response_curves(
  model.non.spatial,
  quantiles = c(0.1, 0.5, 0.9),
  line.color = viridis::viridis(
    3, #same number of colors as quantiles
    option = "F", 
    end = 0.9
  ),
  ncol = 3,
  show.data = TRUE
)

spatialRF::plot_response_curves(
  model.non.spatial,
  quantiles = 0.5,
  ncol = 3
)

reponse.curves.dfx <- spatialRF::get_response_curves(model.non.spatial)
kableExtra::kbl(
  head(reponse.curves.dfx, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)


# Non-Spatial Model performance
spatialRF::print_performance(model.non.spatial)

# Spatial CV
model.non.spatial <- spatialRF::rf_evaluate(
  model = model.non.spatial,
  xy = xy,                  #data coordinates
  repetitions = 30,         #number of spatial folds
  training.fraction = 0.80, #training data fraction on each fold
  metrics = "r.squared",
  seed = 123,
  verbose = FALSE
)

spatialRF::plot_evaluation(model.non.spatial)
spatialRF::print_evaluation(model.non.spatial)

predicted <- stats::predict(
  object = model.non.spatial,
  data = dfx,
  type = "response"
)$predictions

model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  verbose = FALSE,
  seed = 123
)

# Plotting comparative variable importance
p1 <- spatialRF::plot_importance(
  model.non.spatial, 
  verbose = FALSE) + 
  ggplot2::ggtitle("Non-spatial model") 

p2 <- spatialRF::plot_importance(
  model.spatial,
  verbose = FALSE) + 
  ggplot2::ggtitle("Spatial model")
p1 | p2 

kableExtra::kbl(
  head(model.spatial$importance$per.variable, n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)



# Spatial predictors
p <- spatialRF::plot_optimization(model.spatial)


# Tuning Random Forest hyperparameters
model.spatial <- spatialRF::rf_spatial(
  model = model.non.spatial,
  method = "mem.moran.sequential", #default method
  ranger.arguments = list(
    mtry = 2,
    min.node.size = 20,
    num.trees = 500
  ),
  verbose = FALSE,
  seed = 123
)

# Repeating a model execution
model.spatial.repeat <- spatialRF::rf_repeat(
  model = model.spatial, 
  repetitions = 30,
  seed = 123,
  verbose = FALSE
)

# Spatial PDP
spatialRF::plot_response_curves(
  model.spatial.repeat, 
  quantiles = 0.5,
  ncol = 3
)

# Spatial performance
spatialRF::print_performance(model.spatial.repeat)

# Remove the zero distances from the distance matrix
zero_distances <- any(distance.matrix == 0)
if (zero_distances) {
  distance.matrix[distance.matrix == 0] <- 0.001
}

# Perform the spatial random forest analysis
model.full <- rf_spatial(
  data = dfx,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  xy = xy
)

# Comparing several models
comparison <- spatialRF::rf_compare(
  models = list(
    `Non-spatial` = model.non.spatial,
    `Spatial` = model.spatial
  ),
  xy = xy,
  repetitions = 30,
  training.fraction = 0.8,
  metrics = c("r.squared", "rmse"),
  seed = 123
)

x <- comparison$comparison.dfx %>% 
   dplyr::group_by(model, metric) %>% 
   dplyr::summarise(value = round(median(value), 3)) %>% 
   dplyr::arrange(metric) %>% 
   as.data.frame()
colnames(x) <- c("Model", "Metric", "Median")
kableExtra::kbl(
   x,
   format = "html"
) %>%
kableExtra::kable_paper("hover", full_width = F)


model.non.spatial <- spatialRF::rf_evaluate(
  model.non.spatial,
  xy = xy,
  metrics = "auc",
  verbose = FALSE
)
spatialRF::print_evaluation(model.non.spatial)

#single distance (0km by default)
mems <- spatialRF::mem(distance.matrix = distance.matrix)

#several distances
mems <- spatialRF::mem_multithreshold(
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds
)

kableExtra::kbl(
  head(mems[, 1:4], n = 10),
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)


mem.rank <- spatialRF::rank_spatial_predictors(
   distance.matrix = distance.matrix,
   spatial.predictors.dfx = mems,
   ranking.method = "moran"
)

mems <- mems[, mem.rank$ranking]
mems <- mem.rank$spatial.predictors.dfx


predictor.variable.names

#model definition
predictors <- c("Comm_lur", "Income", "T_Enable", "Rd_Pedes", "Rd_Ter", "T_Hinder", "Resi_lur", "Health_lur",
                "AreaProx", "Other_lur", "Edu_lur", "LUP_Other", "Recre_lur", "Rd_Cyc", "Public_lur", "Rd_Sec",
                "Rd_Pri", "Dist_center", "Summer", "Winter", "Spring", "Weekend", "Dist_center..pca..Other_lur",
                "Weekend..x..Summer", "Dist_center..x..AreaProx"  )

model.formula <- as.formula(
  paste(
    dependent.variable.name,
    " ~ ",
    paste(
      predictors,
      collapse = " + "
    )
  )
)

#scaling the data
model.data <- scale(dfx) %>% 
  as.data.frame()

#fitting the model
m <- lm(model.formula, data = model.data)

#Moran's I test of the residuals
moran.test <- spatialRF::moran(
  x = residuals(m),
  distance.matrix = distance.matrix,
  verbose = FALSE
)
moran.test$plot


comparison.dfx <- data.frame(
  Model = c("Non-spatial", "Spatial"),
  Predictors = c(length(predictors), length(predictors.i)),
  R_squared = round(c(summary(m)$r.squared, summary(m.i)$r.squared), 2),
  AIC = round(c(AIC(m), AIC(m.i)), 0),
  BIC = round(c(BIC(m), BIC(m.i)), 0),
  `Moran I` = round(c(moran.test$test$moran.i, moran.test.i$test$moran.i), 2)
)

kableExtra::kbl(
  comparison.dfx,
  format = "html"
) %>%
  kableExtra::kable_paper("hover", full_width = F)
