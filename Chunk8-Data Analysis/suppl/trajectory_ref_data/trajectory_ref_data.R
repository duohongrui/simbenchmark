data_list <- list.files("../preprocessed_data/")

trajectory_plot <- list()

### iterate every reference data and visualize the trajectory
for(i in 1:152){
  ### read data
  data_name <- data_list[i]
  data <- readRDS(file.path("../preprocessed_data", data_name))
  ### data info
  data_info <- data$data_info
  ### ref data
  ref_data <- data$data
  if(!dynwrap::is_data_wrapper(ref_data)){
    next
  }else{
    message(paste0(data_info$id, "..."))
  }
  ### dimentionality reduction
  if(is.null(ref_data[["grouping"]])){
    ref_data <- dynwrap::add_grouping(ref_data, grouping = data_info$cluster_info)
  }
  dimred <- dyndimred::dimred_umap(ref_data$expression)
  cat("Umap...\n")
  p1 <- dynplot::plot_dimred(
    ref_data,
    dimred = dimred,
    expression_source = ref_data$expression, 
    grouping = ref_data$grouping,
    size_cells = 2,
    size_milestones = 1,
    arrow = grid::arrow(length = grid::unit(0.1, "inches"))
  )
  ### spatial data
  if(!is.null(data_info[["spatial_coordinate"]])){
    cat("Spatial coordinate...\n")
    p2 <- dynplot::plot_dimred(
      ref_data,
      dimred = data_info$spatial_coordinate,
      expression_source = ref_data$expression,
      grouping = ref_data$grouping,
      size_cells = 2,
      size_milestones = 1,
      arrow = grid::arrow(length = grid::unit(0.1, "inches"))
    )
    trajectory_plot[[data_info$id]] <- list(umap = p1, spatial = p2)
  }else{
    trajectory_plot[[data_info$id]] <- list(umap = p1)
  }
}
saveRDS(trajectory_plot, file = "./Chunk8-Data Analysis/suppl/trajectory_ref_data/trajectory_plot.rds")
