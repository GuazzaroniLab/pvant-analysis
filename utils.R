### Loading needed pkgs ----
if (!require("flowCore")) {
  BiocManager::install("flowCore")
  library(flowCore)
}

if (!require("ggcyto")) {
  BiocManager::install("ggcyto")
  library(ggcyto)
}

if (!require("flowStats")) {
  BiocManager::install("flowStats")
  library(flowStats)
}

if (!require("tidyverse")) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require("ggridges")) {
  install.packages("ggridges")
  library(ggridges)
}

if (!require("gt")) {
  install.packages("gt")
  library(gt)
}

### Define a manual color vector for coloring the plots consistently ----
basecolors <- RColorBrewer::brewer.pal(9, name = "Spectral") |>
  gt::adjust_luminance(steps = -0.1)

colorPalette <- basecolors[c(3:1, 7:9)]



### Defines a function for importing cytometry data ----
importCytometry <- function(input.path) {
  # listing all of the flow cytometry files in their data/ subdirectory
  # importing the Flow Cytometry data using flowCore
  flowFiles <- list.files(input.path,
                          pattern = ".fcs",
                          full.names = TRUE,
                          recursive = TRUE)

  # reading the data into a list
  flowData <- purrr::map(flowFiles, read.FCS)

  # the lines below name the elements of the flowData list using metadata information from the .fcs files themselves
  names(flowData) <- purrr::map_chr(.x = flowData, ~ {
    keyword(.x, keyword = c("EXPERIMENT NAME", "$SRC", "TUBE NAME")) |>
      unlist() |>
      glue::glue_collapse(sep = ".") |>
      stringr::str_replace(pattern = " ",
                           replacement = "-") })

  # finally, we join the individual flowFrames into a single flowSet object
  fs <- as(flowData, "flowSet")


  # lastly, we will store the metadata information in a tibble for later use (including colors, etc)
  metadata <- names(flowData) %>%
    tibble("SampleName" = .) %>%
    dplyr::filter(!str_detect(string = SampleName, pattern = "LT")) %>%
    mutate(vector = case_when(
      str_detect(SampleName, "SEVA") ~ "pSEVA231",
      str_detect(SampleName, "VANT") ~ "pVANT"),
      sampletype = case_when(
        str_detect(SampleName, "-SP") ~ "-*gfplva*",
        str_detect(SampleName, "(?<!SP)BRANCO") ~ "-∅",
        str_detect(SampleName, "GFP") ~ "-P*j100*-*gfplva*"),
      batch = lubridate::dmy(str_extract_all(SampleName, "^\\d{8}", simplify = TRUE)),
      replicate = glue::glue("R{as.numeric(as.factor(batch))}"),
      readableName = glue::glue("{vector}{sampletype}") %>%
        as_factor() %>%
        fct_relevel(c("pVANT-∅",
                      "pVANT-*gfplva*",
                      "pVANT-P*j100*-*gfplva*",
                      "pSEVA231-∅",
                      "pSEVA231-*gfplva*",
                      "pSEVA231-P*j100*-*gfplva*")
        ),
      light_color = case_when(
        readableName == "pSEVA231-∅" ~ colorPalette[1],
        readableName == "pSEVA231-*gfplva*" ~ colorPalette[2],
        readableName == "pSEVA231-P*j100*-*gfplva*" ~ colorPalette[3],
        readableName == "pVANT-∅" ~ colorPalette[4],
        readableName == "pVANT-*gfplva*" ~ colorPalette[5],
        readableName == "pVANT-P*j100*-*gfplva*" ~colorPalette[6],
      ),
      dark_color = gt::adjust_luminance(light_color, -1)
    )

  # returning the final result
  return <- dplyr::lst(flowSet = fs, names = names(flowData), metadata = metadata)

}


### Defines a function for retrieving cytometry data as a tibble ----
cyto2tibble <- function(fs, names, metadata) {
  fs_data_accessor <- glue::glue("fs@frames$`{names}`@exprs")
  fs_tibble <- purrr::map_dfr(fs_data_accessor,
                              ~ {
                                tibble::tibble(
                                  SampleName = stringr::str_extract(.x,
                                                                    pattern = "(?<=\\$`).+(?=`@)"),
                                  as_tibble(eval(parse(text = .x))))})  %>%
                                  left_join(metadata)
}



### Defines a function for calculating gate percentages ----

calcGated <- function(original_data, gated_data) {

  pre_counts <- count(original_data, SampleName, name = "pre_gating")
  post_counts <- count(gated_data, SampleName, name = "post_gating")

  metrics <- left_join(pre_counts, post_counts) %>% mutate(retained = post_gating/pre_gating)

  return(metrics)
}


### Defines custom plotting functions ----


vizParam <- function(replicate_choice, data) {

  replicate_choice <- replicate_choice

  fill_vector <- data$metadata %>%
    select(readableName, light_color) %>%
    distinct() %>%
    deframe()

  color_vector <- data$metadata %>%
    select(readableName, dark_color) %>%
    distinct() %>%
    deframe()

  return(list(replicate_choice = replicate_choice,
              fill_vector = fill_vector,
              color_vector = color_vector))
}


raincloudPlot <- function(data, visParam) {

  raincloud_plot <- data %>%
    dplyr::filter(replicate == visParam$replicate_choice ) %>%
    ggplot(aes(x = `FITC-A`, y = readableName, fill = readableName, color = readableName)) +
    coord_cartesian(xlim=c(-100, 1e5)) +
    scale_x_flowjo_biexp(expand = expansion(add = c(0, 100))) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.75)), limits = rev) +
    scale_fill_manual(values = visParam$fill_vector) +
    scale_color_manual(values = visParam$color_vector) +
    geom_density_ridges(scale = 0.5, position = position_nudge(y = 0.15), rel_min_height = 0.005) +
    geom_point(size = 0.55, alpha = 0.1, position = position_jitter(height = .125, width = 1)) +
    geom_boxplot(position = position_nudge(y = -0.2),
                 width = .1,
                 outlier.shape = NA) +
    labs(y = "", x = "Fluorescence intensity (A. U. )") +
    theme(


      panel.background = element_blank(),

      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = ggtext::element_markdown(size = rel(1.6), color = "black"), #

      axis.line.x = element_line(),
      axis.text.x = ggtext::element_markdown(size = rel(1.6), vjust = 0, color = "black"),
      axis.title.x = element_text(size = rel(1.6), margin = margin(20, 0, 0, 0), color = "black"),

      panel.grid.major.y = element_line(colour = "gray90", size = 0.5, linetype = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  return(raincloud_plot)
}


distributionPlot <- function(data, visParam, main, threshold, percentSize = 5.75) {

  percentages <- data %>%
    dplyr::filter(replicate == visParam$replicate_choice) %>%
    mutate(abv_thresh = if_else(`FITC-A` >= threshold, "yes", "no")) %>%
    group_by(readableName, abv_thresh) %>%
    tally() %>%
    mutate(total = sum(n),
           percentage = (n/total)*100,
           percentage_formatted = glue::glue("{round(percentage, 1)}%")) %>%
    dplyr::filter(abv_thresh == "yes") %>%
    select(readableName, percentage_formatted)



  plot <- data %>%
    dplyr::filter(replicate == visParam$replicate_choice ) %>%
    ggplot(aes(x = `FITC-A`, y = readableName, fill = readableName)) +
    coord_cartesian(xlim=c(-100, 1e5)) +
    geom_text(data = percentages,
              aes(y = readableName, x = threshold+15, label = percentage_formatted),
              position = position_nudge(y = 0.5),
              size = percentSize,
              hjust = 0) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "gray60") +
    scale_x_flowjo_biexp(expand = expansion(add = c(0, 100))) +
    scale_y_discrete(expand = expansion(add = c(0.2, 1)), limits = rev, position = "right") +
    scale_fill_manual(values = visParam$fill_vector) +
    geom_density_ridges(scale = 0.75, color = "black", size = 0.6, position = position_nudge(y = 0.2), rel_min_height = 0.01) +
    geom_boxplot(position = position_nudge(y = 0.1),
                 color = "black",
                 size = 0.4,
                 width = 0.1,
                 outlier.shape = NA) +
    labs(title = main, y = "", x = "Fluorescence intensity (A. U. )") +
    theme(

      plot.margin = margin(t = 0, r = 250, b = 0, l = 250),
      panel.background = element_blank(),
      plot.title = ggtext::element_markdown(size = rel(1.45), margin = margin(10, 0, 20, 0), hjust = 0.5),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y.right = ggtext::element_markdown(vjust = -1.5, size = rel(1.6), color = "black"), #

      axis.line.x = element_line(),
      axis.text.x = ggtext::element_markdown(size = rel(1.9), vjust = 0, color = "black"),
      axis.title.x = element_text(size = rel(1.6), margin = margin(20, 0, 0, 0), color = "black"),

      panel.grid.major.y = element_line(colour = "gray90", size = 0.5, linetype = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )

  return(plot)

}
