# Original formula https://www.l3harrisgeospatial.com/docs/principalcomponentanalysis.html
# PCA in rgee
pca_rgee <- function(x){
  arrays <- x$toArray()
  covar  <- arrays$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = nc,
    scale = 100
  )
  
  covarArray = ee$Array(covar$get('array'))
  eigens = covarArray$eigen()
  eigenValues = eigens$slice(1, 0, 1)
  eigenVectors = eigens$slice(1, 1)
  arrayImage = arrays$toArray(1)
  pca = ee$Image(eigenVectors)$matrixMultiply(arrayImage)
  
  sd_img = ee$Image(eigenValues$sqrt())$
    arrayProject(list(0))$arrayFlatten(list(c("sd1","sd2","sd3","sd4","sd5","sd6")))

  pca_raster <- pca$
    arrayProject(list(0))$
    arrayFlatten(list(c("pc1","pc2","pc3","pc4","pc5","pc6")))$
    divide(sd_img)
  return(pca_raster) 
  
}

# Standar image
img_to_stand <- function(x){
  mean <- x$reduceRegion(
  reducer =  ee$Reducer$mean(),
  geometry = nc,
  scale =  100,
  maxPixels = 10**19
  )
  values_mean <- mean$values()
  img_mean <- ee$Image$constant(values_mean)
  centered <- x$subtract(img_mean)
  return(centered)
}  


# Importance of component 

impor_pca <- function(x){
  arrays <- x$toArray()
  covar  <- arrays$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = nc,
    scale = 100,
    maxPixels = 10**19
  )
  
  covarArray = ee$Array(covar$get('array'))
  eigens = covarArray$eigen()
  eigenValues = eigens$slice(1, 0, 1)
  eigenVectors = eigens$slice(1, 1)
  arrayImage = arrays$toArray(1)
  pca = ee$Image(eigenVectors)$matrixMultiply(arrayImage)
  
  sd_img = ee$Image(eigenValues$sqrt())$
    arrayProject(list(0))$arrayFlatten(list(c("sd1","sd2","sd3","sd4","sd5","sd6")))
  
  eing_values <- eigenValues$getInfo() %>%
    as.list() %>%
    map_df(.f = as.data.frame,.id = 'eingvalues') %>% 
    rename(value =`.x[[i]]`)
  
  var_comp <- eing_values %>% 
    mutate(Component = sprintf('PC%s',eingvalues),
           variance = (value*100/sum(value)),
           cumulative = cumsum(variance))
  
  return(var_comp)
  
}

eingvector_rgee <- function(x){
  arrays <- x$toArray()
  covar  <- arrays$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = nc,
    scale = 30,
    maxPixels = 10**19
  )
  
  covarArray = ee$Array(covar$get('array'))
  eigens = covarArray$eigen()
  eigenValues = eigens$slice(1, 0, 1)
  eigenVectors = eigens$slice(1, 1)
  
  table_eingvalues <- eigenVectors$getInfo() %>%
    map_dfr(~unlist(.x) %>% t() %>% as.data.frame(),.id = 'eingvector') %>% 
    rename_at(2:7,~c("nldc","imp" ,"tree","ndvi","tmax","tmin"))

  return(table_eingvalues)
}

eingvalues_rgee <- function(x){
  arrays <- x$toArray()
  covar  <- arrays$reduceRegion(
    reducer = ee$Reducer$centeredCovariance(),
    geometry = nc,
    scale = 30,
    maxPixels = 10**19
  )
  
  covarArray = ee$Array(covar$get('array'))
  eigens = covarArray$eigen()
  eigenValues = eigens$slice(1, 0, 1)
  eing_values <- eigenValues$getInfo() %>%
    map_df(as.data.frame,.id = 'Eig') %>% 
    rename(values =`.x[[i]]`) %>%
    mutate(Eig = paste0("Eig.",Eig))
  
  return(eing_values)
}

# Resample size 

# resample_size <- function(x){
#   reduce_size <- x$reproject(
#     crs = nldc$projection(),
#     scale = 100
#   )
#   return(reduce_size)
# }