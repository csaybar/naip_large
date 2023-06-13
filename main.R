#remotes::install_github("r-spatial/rgee", branch = "ee_as_rast")
#remotes::install_github("r-spatial/rgeeExtra", branch = "ee_as_rast")

library(geojsonio)
library(rgeeExtra)
library(terra)
library(rgee)
library(sf)

source("utils.R")

ee_Initialize()

dataset <- read_sf("data/ROI.geojson")


for (index in 1:5000) {
    print(index)
    downloadNAIP(dataset[index,], output="/media/csaybar/0790BB3D255A0B7F/NAIPLARGE2/")
}

