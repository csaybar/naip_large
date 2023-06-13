library(lubridate)
library(dplyr)
library(terra)

downloadNAIP <- function(point, output) {
    # From sf to ee
    roi_id <- point$ROI
    ee_point <- sf_as_ee(point$geometry)
    
    # Get potential NAIP images
    hr_img <- ee$ImageCollection("USDA/NAIP/DOQQ") %>%
        ee$ImageCollection$filterBounds(ee_point) %>%
        ee$ImageCollection$filterDate("2010-01-01", "2024-12-31")
    
    # Get the dates
    hr_metadata <- ee_get_date_ic(hr_img)

    # Find the images to download
    to_download <- find_images(hr_metadata)
    hr_download <- hr_metadata[to_download,]

    if(sum(to_download) == 0) {
        return(0)
    } 
    
    # Create a tentative mask 
    ee_mask1 <- sapply(
        X = hr_download$id, 
        FUN = function(x) ee$Image(x)$unmask(0, sameFootprint = FALSE)$reduce("mean") == 0
    )
    ee_mask2 <- (ee_mask1[[1]] + ee_mask1[[2]])$eq(0)
    
    for (index in 1:2) {
        ee_img_id <- hr_download[index,]$id
        ee_img <- ee$Image(ee_img_id)$unmask(0, sameFootprint = FALSE)*ee_mask2

        # Obtain the CRS - Intersect the geometry with the UTM grid 
        ee_img_ref <- ee$ImageCollection("COPERNICUS/S2")$filterBounds(ee_point)$first()
        CRSS2 <- ee_img_ref$select("B2")$projection()$getInfo()$crs

        # Define the ROI
        roi <- st_transform(point, CRSS2) %>%
            st_buffer(1280, endCapStyle = "SQUARE") %>%
            sf_as_ee(proj = CRSS2)

        # Create the folder
        roi_folder <- sprintf("%s/%s/", output, roi_id)
        dir.create(roi_folder, recursive = TRUE, showWarnings = FALSE)

        # Start the download
        dsn1 <- sprintf("%s/%s__int8.tif", roi_folder, basename(ee_img_id))
        dsn2 <- gsub("__int8", "", dsn1)

        ee_as_rast(
            image = ee_img,
            region = roi$geometry(),
            scale = 1,
            via = "getDownloadURL",
            dsn = dsn1,
            crs = CRSS2
        )

        # clear tmp folder
        unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)

        # From float32 to int8
        writeRaster(
            terra::rast(dsn1),
            dsn2,
            overwrite = TRUE,
            datatype = "INT1U",
            gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES")
        )
        unlink(dsn1)
    }
}


find_images <- function(hr_metadata) {
    if(nrow(hr_metadata) < 2) {
        return(0)
    }

    # random pick
    idx <- sample(nrow(hr_metadata), 1)
    sample <- hr_metadata[idx,]
    years1 <- year(hr_metadata$time_start)
    years_diff <- abs(years1 - years1[idx]) * 6

    months1 <- month(hr_metadata$time_start)
    months_diff <- abs(months1 - months1[idx]) * 3

    # Calculate the distance between the points
    ref_points <- cbind(years_diff, months_diff)
    point <- c(years_diff[idx], months_diff[idx])
    dist <- sqrt(rowSums((ref_points - point)^2))
    high_idx <- which.max(dist)
    c(idx, high_idx)
}
