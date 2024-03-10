library(spdep)
library(dplyr)
data <- readRDS("../preprocessed_data/data123_spatial_DLPFC1.rds")

locations <- data$data_info$spatial_coordinate
count_data <- data$data %>% as.matrix()


#### All SVGs
all_svgs
count_data <- count_data[all_svgs, ]

moran_result <- tibble()

for(i in 1:length(all_svgs[11000:12000])){
  print(i)
  candicate_data <- locations %>%
    mutate(
      "gene" = count_data[i, ]
    )
  error <- try(
    result <- moranfast::moranfast(candicate_data$gene, candicate_data$x, candicate_data$y),
    silent = TRUE
  )
  # error <- try(
  #   result <- ape::Moran.I(x = count_data[i, ], weight = ozone.dists.inv),
  #   silent = TRUE
  # )
  if(is(error, "try-error")){
    next
  }
  moran_result <- moran_result %>% 
    rbind(tibble(
      "data" = "data123",
      "gene" = rownames(count_data)[i],
      "Moran" = result$observed
    ))
}


# ozone.dists <- as.matrix(dist(locations))
# ozone.dists.inv <- 1/ozone.dists
# diag(ozone.dists.inv) <- 0
# ape::Moran.I(candicate_data$gene, ozone.dists.inv)
# listw <- mat2listw(ozone.dists.inv, style = "W", zero.policy = TRUE)
# moran(x = candicate_data$gene, listw = listw, n = ncol(count_data), S0 = Szero(listw), zero.policy = TRUE)
# geary(x = candicate_data$gene, listw = listw, n = ncol(count_data), n1 = ncol(count_data) - 1, S0 = Szero(listw), zero.policy = TRUE)