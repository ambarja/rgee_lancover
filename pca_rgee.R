#' authors: Fernando Prudencio y Antony Barja
#' source: https://www.l3harrisgeospatial.com/docs/principalcomponentanalysis.html

# PCA in rgee
ee_pca <-
  function(image, ee_feature, scale, nvar) {
    # collapse the bands of the image into a 1D array per pixel
    arrays <- image$toArray()
    # compute the covariance of the bands within the region
    covar <-
      arrays$reduceRegion(
        reducer = ee$Reducer$centeredCovariance(),
        geometry = ee_feature,
        scale = scale,
        maxPixels = 10**19
      )
    # get the 'array' covariance result and cast to an array
    # this represents the band-to-band covariance within the region
    covarArray <- ee$Array(covar$get("array"))
    # perform an eigen analysis and slice apart the values and vectors
    eigens <- covarArray$eigen()
    # this is a P-length vector of Eigenvalues
    eigenValues <- eigens$slice(1, 0, 1)
    # this is a PxP matrix with eigenvectors in rows
    eigenVectors <- eigens$slice(1, 1)
    # convert the array image to 2D arrays for matrix computations
    arrayImage <- arrays$toArray(1)
    # left multiply the image array by the matrix of eigenvectors
    pca <- ee$Image(eigenVectors)$matrixMultiply(arrayImage)
    # Turn the square roots of the Eigenvalues into a P-band image
    sdImage <-
      ee$Image(eigenValues$sqrt())$
        arrayProject(list(0))$
        arrayFlatten(list(sprintf("sd%1$s", 1:nvar)))
    # Turn the PCs into a P-band image, normalized by SD
    pca$arrayProject(list(0))$
      arrayFlatten(list(sprintf("pc%1$s", 1:nvar)))$
      divide(sdImage) %>%
      return()
  }

# standardize variables
ee_scale <-
  function(image, ee_feature, scale, namevar) {
    for (i in namevar) {
      # average value by band
      meanValues <-
        image$select(i)$
        reduceRegion(
          reducer = ee$Reducer$mean(),
          geometry = ee_feature,
          scale = scale,
          maxPixels = 10**19
        )
      meanImage <- ee$Image$constant(meanValues$values())
      # average value by band
      sdValues <-
        image$select(i)$
        reduceRegion(
          reducer = ee$Reducer$stdDev(),
          geometry = ee_feature,
          scale = scale,
          maxPixels = 10**19
        )
      sdImage <- ee$Image$constant(sdValues$values())
      # calculate standardized variables
      centred <- image$select(i)$subtract(meanImage)
      scaled <- centred$divide(sdImage)$rename(i)
      if (i == namevar[1]) {
        stack <- ee$Image(c(scaled))$toDouble()
      } else {
        stack <- ee$Image(c(stack, scaled))$toDouble()
      }
    }
    # return image
    return(stack)
  }

# Import eigenvectors
eVectors <- function(image, ee_feature, scale, nvar) {
  # collapse the bands of the image into a 1D array per pixel
  array <- image$toArray()
  # compute the covariance of the bands within the region
  covar <- array$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = ee_feature,
    scale = scale,
    maxPixels = 10**19
  )
  # get the 'array' covariance result and cast to an array
  # this represents the band-to-band covariance within the region
  covarArray <- ee$Array(covar$get("array"))
  # perform an eigen analysis and slice apart the values and vectors
  eigens <- covarArray$eigen()
  # this is a PxP matrix with eigenvectors in rows
  eigenVectors <- eigens$slice(1, 1)
  # get eigenVectors
  eigenVectors$getInfo() %>%
    map_dfr(
      ~ unlist(.x) %>%
        t() %>%
        as.data.frame(),
      .id = "eingvector"
    ) %>%
    rename_at(
      2:(nvar + 1), ~ sprintf("eVec%1$s", 1:nvar)
    ) %>%
    mutate_at(2:(nvar + 1), round, 3) %>%
    return()
}

# Import eigenvalues
eValues <- function(image, ee_feature, scale, nvar) {
  # collapse the bands of the image into a 1D array per pixel
  array <- image$toArray()
  # compute the covariance of the bands within the region
  covar <- array$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = ee_feature,
    scale = scale,
    maxPixels = 10**19
  )
  # get the 'array' covariance result and cast to an array
  # this represents the band-to-band covariance within the region
  covarArray <- ee$Array(covar$get("array"))
  # perform an eigen analysis and slice apart the values and vectors
  eigens <- covarArray$eigen()
  # this is a P-length vector of Eigenvalues
  eigenValues <- eigens$slice(1, 0, 1)
  # get eigenValues
  eigenValues$getInfo() %>%
    map_df(as.data.frame, .id = "eingvalue") %>%
    rename(values = `.x[[i]]`) %>%
    mutate(
      eingvalue = sprintf("eVal%1$s", 1:nvar),
      values = round(values, 2)
    ) %>%
    return()
}

# Importance of component
imporPCA <- function(image, ee_feature, scale, nvar) {
  # get eigenvalues
  eigenvalues <-
    eValues(
      image = image,
      ee_feature = ee_feature,
      scale = scale,
      nvar = nvar
    )
  eigenvalues %>%
    mutate(
      component = sprintf("pc%1$s", 1:nvar),
      variance = (values * 100 / sum(values)),
      cumulative = cumsum(variance)
    ) %>% return()
}

# Resample size
resample_size <- function(x) {
  reduce_size <- x$reproject(
    crs = nldc$projection(),
    scale = 100
  )
  return(reduce_size)
}
