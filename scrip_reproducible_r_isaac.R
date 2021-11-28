#video_3.----
library ( rgrass7 )
gisdbase <- 'grass-data-test' #Base de datos de GRASS GIS
wd <- getwd() #Directorio de trabajo
wd
## [1] "/home/jr/unidad-4-asignacion-1-procesos-fluviales"
loc <- initGRASS(gisBase = "/usr/lib/grass78/",
                 home = wd,
                 gisDbase = paste(wd, gisdbase, sep = '/'),
                 location = 'soco',
                 mapset = "PERMANENT",
                 override = TRUE)

#video_4 -----

gmeta ()
dem  <-  'datos-fuente/srtm_dem_cuenca_soco.tif'
execGRASS(
  cmd  =  'g.proj' ,
  flags  = c ('t','c'),
  georef=dem)
gmeta ()

execGRASS(
  cmd = 'r.in.gdal',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=dem,
    output='dem'
  )
)

execGRASS(
  cmd = 'g.region',
  parameters=list(
    raster = 'dem',
    align = 'dem'
  )
)

gmeta()

demext <- 'datos-fuente/srtm_dem_cuenca_soco.geojson'
execGRASS(
  cmd = 'v.in.ogr',
  flags=c('overwrite','quiet'),
  parameters=list(
    input=demext,
    output='dem_extent'
  )
)

# Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

execGRASS(
  cmd = 'g.extension',
  flags = 'c'
)


parseGRASS("r.in.gdal")

system('r.in.gdal --help')


#Video 5 Explorar datos espaciales básicos entre GRASS y R----


#Imprimir lista de mapas ráster y dentro del vector en la región / localización activa


execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

library(sp)
use_sp()
dem_sp <- readRAST('dem')
op <- par()
plot(dem_sp)

library(sf)
rutasoco <- 'datos-fuente/cuenca_soco.geojson'
soco <- st_read(rutasoco)

plot(dem_sp)
plot(soco, add=T, col='transparent', border='green', lwd=5);par(op[c('mfrow','mar')])

library(raster)
dem_r0 <- raster(dem_sp)
dem_r1 <- crop(dem_r0, soco)
dem_soco <- mask(dem_r1, soco)
plot(dem_soco, add=T, col='transparent', border='green', lwd=3)

summary(dem_soco)
hist(dem_soco)

#obtener variables de terreno básicas con el paquete raster dentro de R

pend_soco <- terrain(x = dem_soco, opt = 'slope', unit = 'degrees')
plot(pend_soco)
plot(soco, add=T, col='transparent', border='green', lwd=3)

summary(pend_soco)
hist(pend_soco)

#Obtener la misma variable de terreno con GRASS GIS

writeVECT(as_Spatial(soco), 'soco', v.in.ogr_flags='quiet')
execGRASS(
  "g.region",
  parameters=list(
    vector = "soco"
  )
)
execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'soco'
  )
)
execGRASS(
  cmd = 'r.slope.aspect',
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation='dem',
    slope='slope',
    aspect='aspect',
    pcurvature='pcurv',
    tcurvature='tcurv')
)
pend_soco_g <- readRAST('slope')
plot(pend_soco_g);par(op[c('mfrow','mar')])

summary(pend_soco_g)
summary(pend_soco)
gmeta()

execGRASS(
  "g.region",
  parameters=list(
    raster = "dem"
  )
)
execGRASS(
  "r.mask",
  flags = c('r','quiet')
)
gmeta()
#Limpiar archivo de bloqueo del conjunto de mapas de GRASS
unlink_.gislock()
#Video 6 ----

#lista de mapas cargados----
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

execGRASS(
  "r.watershed",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = "dem",
    accumulation = "accum-de-rwshed",
    stream = "stream-de-rwshed",
    drainage = "drainage-dir-de-rwshed",
    basin = 'basins',
    half_basin = 'half-basins',
    threshold = 80
  )
)
#Usar Spatial*
library(sp)
use_sp()
#Paquete manejo de los raster
library(raster)
#DEM
dem <- raster(readRAST('dem'))
#Basins
basins <- raster(readRAST('basins'))
#Stream network
stream <- raster(readRAST('stream-de-rwshed'))
stream3857 <- projectRaster(stream, crs = CRS("+init=epsg:3857"), method = 'ngb')
#Generar un vectorial de extensión de capa en EPSG:4326
e <- extent(stream)
e <- as(e, 'SpatialPolygons')
proj4string(e) <- CRS("+init=epsg:32619")
e <- spTransform(e, CRSobj = CRS("+init=epsg:4326"))

library(leaflet)
library(leafem)
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addRasterImage(dem, group='DEM', opacity = 0.5) %>%
  addRasterImage(
    ratify(basins),
    group='basins', opacity = 0.7,
    colors = sample(rep(RColorBrewer::brewer.pal(12, 'Set3'),1000))) %>% 
  addRasterImage(stream3857, project = F, group='str', opacity = 0.7, method = 'ngb', colors = 'blue') %>% 
  addLayersControl(
    overlayGroups = c('terrain','DEM','basins','str'),
    options = layersControlOptions(collapsed=FALSE)) %>% 
  addHomeButton(extent(e), 'Ver todo')

unlink_.gislock()

#video 7 -----
library ( rgrass7 )
gisdbase <- 'grass-data-test' #Base de datos de GRASS GIS
wd <- getwd() #Directorio de trabajo
wd
## [1] "/home/jr/unidad-4-asignacion-1-procesos-fluviales"
loc <- initGRASS(gisBase = "/usr/lib/grass78/",
                 home = wd,
                 gisDbase = paste(wd, gisdbase, sep = '/'),
                 location = 'soco',
                 mapset = "PERMANENT",
                 override = TRUE)
## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Obtener las coordenadas de la desembocadura de la cuenca de interés

library(mapview)
mapview(
  stream3857, method='ngb', col.regions = 'blue',
  legend = FALSE, label = FALSE, maxpixels =  910425
)


## Convertir las coordenadas lat/lon a EPSG:32619


my_trans <- function(coords = NULL) {
  require(sp)
  pt <- SpatialPoints(matrix(coords, ncol = 2), CRS("+init=epsg:4326"))
  foo <- spTransform(pt, CRSobj = CRS("+init=epsg:32619"))
  bar <- as.vector(coordinates(foo))
  return(bar)
}
soco_out <- my_trans(coords = c(-69.20322,18.45688))
soco_out
#este ultimo es el verdadero

## Extraer la cuenca de interés

execGRASS(
  "r.water.outlet",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'drainage-dir-de-rwshed',
    output = 'soco-basin',
    coordinates = soco_out
  )
)

## Convertir la cuenca a vectorial en GRASS

execGRASS(
  "r.to.vect",
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'soco-basin',
    output = 'soco_basin',
    type = 'area'
  )
)

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Traer a R la cuenca del soco
soco_bas <- readVECT('soco_basin')
soco_bas
plot(soco_bas)
soco_bas4326 <- spTransform(soco_bas, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain) %>%
  addRasterImage(stream, opacity = 0.7, method = 'ngb', colors = 'blue') %>% 
  addPolygons(data = soco_bas4326) %>% 
  leafem::addHomeButton(extent(soco_bas4326), 'Ver cuenca')
#video 8  Extraer una red drenaje con r.stream.extract. Visualizar con leaflet----

unlink_.gislock()

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)


## Usar la cuenca del arroyo Pantuflas como máscara


execGRASS(
  "r.mask",
  flags = c('verbose','overwrite','quiet'),
  parameters = list(
    vector = 'soco_basin'
  )
)

## Extraer la red de drenaje de la cuenca de interés


execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'dem',
    threshold = 80,
    stream_raster = 'soco-stream-de-rstr',
    stream_vector = 'soco_stream_de_rstr'
  )
)

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

soco_net <- readVECT('soco_stream_de_rstr', ignore.stderr = T)
soco_net
plot(soco_net)
soco_net4326 <- spTransform(soco_net, CRSobj = CRS("+init=epsg:4326"))
soco_net4326
soco_centroid <- coordinates(rgeos::gCentroid(soco_bas4326))
soco_centroid
soco_net_r <- raster(readRAST('soco-stream-de-rstr'))
soco_net_r
soco_net_r3857 <- projectRaster(soco_net_r, crs = CRS("+init=epsg:3857"), method = 'ngb')
soco_net_r3857
leaflet() %>% 
  setView(lng = soco_centroid[1], lat = soco_centroid[2], zoom = 11) %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addRasterImage(soco_net_r3857, opacity = 0.7, method = 'ngb', colors = 'grey20', group = 'str_raster') %>% 
  addPolylines(data = soco_net4326, weight = 3, opacity = 0.7, group = 'str_vect') %>% 
  leafem::addHomeButton(extent(soco_net4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','str_vect','str_raster'),
    options = layersControlOptions(collapsed=FALSE)) 

#Video 9 : Orden de red y razón de bifurcación explicados-----

#Solo explicacion

#Video 10: ----

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Crear mapa de dirección de flujo a partir de r.stream


execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'dem',
    threshold = 80,
    direction = 'drainage-dir-de-rstr'
  )
)

## Crear mapas de órdenes de red


execGRASS(
  "r.stream.order",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'soco-stream-de-rstr',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    accumulation = 'accum-de-rwshed',
    stream_vect = 'order_all',
    strahler = 'order-strahler',
    horton = 'order-horton',
    shreve = 'order-shreve',
    hack = 'order-hack-gravelius',
    topo = 'order-topology'
  )
)

## Mostrar lista nuevamente


execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Visualizar la red con leaflet

### Simbología única


order <- readVECT('order_all')



order4326 <- spTransform(order, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = order4326, weight = 3, opacity = 0.7, group = 'order',
    label = ~as.character(strahler),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = T,
                                style = list(
                                  "font-size" = "8px",
                                  "background" = "rgba(255, 255, 255, 0.5)",
                                  "background-clip" = "padding-box",
                                  "padding" = "1px"))) %>% 
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','order'),
    options = layersControlOptions(collapsed=FALSE))

### Simbología aplicando grosor según orden de red


leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = order4326, weight = order4326$strahler*1.5, opacity = 0.7, group = 'order',
    label = ~as.character(strahler),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = F)) %>% 
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','order'),
    options = layersControlOptions(collapsed=FALSE))

## Delimitar cuencas según orden de red de Strahler

### Obtener órdenes de red mínimo y máximo


#Estadísticas para obtener los valores mínimo y máximo del orden de red de Strahler
rinfo.ordstra <- execGRASS(
  'r.info',
  flags = 'r',
  parameters = list(
    map = 'order-strahler'
  )
)
#Órdenes de red mínimo y máximo
minmaxord <- as.numeric(
  stringr::str_extract_all(
    attributes(rinfo.ordstra)$resOut,
    "[0-9]+"
  )
)
minmaxord

### Delimitar cuencas, convertirlas de ráster a vectorial


sapply(
  min(minmaxord):max(minmaxord),
  function(x){
    execGRASS(
      "r.stream.basins",
      flags = c('overwrite','c','quiet'),
      parameters = list(
        direction = 'drainage-dir-de-rstr',
        stream_rast = 'order-strahler',
        cats = as.character(x),
        basins = paste0('r-stream-basins-',x)
      )
    )
    execGRASS(
      "r.to.vect",
      flags=c('overwrite','quiet'),
      parameters = list(
        input = paste0('r-stream-basins-',x),
        output = paste0('r_stream_basins_',x),
        type = 'area'
      )
    )
  }
)

### Representar las cuencas con leaflet

sapply(
  min(minmaxord):max(minmaxord),
  function(x){
    assign(
      paste0('orden', x),
      spTransform(readVECT(paste0('r_stream_basins_',x)), CRSobj = CRS("+init=epsg:4326")),
      envir = .GlobalEnv)
  }
)

paleta <- RColorBrewer::brewer.pal(12, 'Set3')
leaflet() %>% 
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolygons(data = orden4, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O4') %>% 
  addPolygons(data = orden3, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O3') %>%
  addPolygons(data = orden2, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O2') %>%
  addPolygons(data = orden1, stroke = T, weight = 2,
              color = ~paleta, fillOpacity = 0.4, group = 'O1') %>%
  addPolylines(
    data = order4326, weight = order4326$strahler*1.5,
    opacity = 0.7, group = 'str_order') %>%
  leafem::addHomeButton(extent(order4326), 'Ver todo') %>% 
  addLayersControl(
    overlayGroups = c('terrain','O1','O2','O3','O4','str_order'),
    options = layersControlOptions(collapsed=FALSE))

## Estadísticas de red resumidas por orden de red.


execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet','o'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    output = 'soco_stats.txt'
  )
)
file.show('soco_stats.txt')
d <- read.csv("soco_stats.txt", skip=1, header=TRUE)
plot(num_of_streams~order, data=d, log="y")
mod <- lm(log10(num_of_streams)~order, data=d)
abline(mod)
text(2, 20, 'logN=2.064-0.544u')
rb <- 1/10^mod$coefficients[[2]]
rb

## Estadísticas de red ampliadas


execGRASS(
  "r.stream.stats",
  flags = c('overwrite','quiet'),
  parameters = list(
    stream_rast = 'order-strahler',
    direction = 'drainage-dir-de-rstr',
    elevation = 'dem',
    output = 'soco_stats_expanded.txt'
  )
)
file.show('soco_stats_expanded.txt')

#Video 11-----


execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)


## Obtener coordenada
mapview(order, col.regions = 'blue', legend = FALSE)
 
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/lfp_network.R') 
#Cargada como función "LfpNetwork"
LfpNetwork(
  xycoords = my_trans(c(-69.19744, 18.45272)),
  suffix ='soco',
  stream_vect ='order_all',
  direction ='drainage-dir-de-rstr'
)

## Mostrar lista nuevamente


execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

##Representar con leaflet


lfp <- readVECT('LfpNetwork_lfp_all_final_soco')



lfp4326 <- spTransform(lfp, CRSobj = CRS("+init=epsg:4326"))
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(
    data = lfp4326, weight = 3, opacity = 0.7, group = 'order',
    label = ~as.character(cat),
    highlightOptions = highlightOptions(color = "white",
                                        weight = 5, bringToFront = F, opacity = 1),
    labelOptions = labelOptions(noHide = T,
                                style = list(
                                  "font-size" = "8px",
                                  "background" = "rgba(255, 255, 255, 0.5)",
                                  "background-clip" = "padding-box",
                                  "padding" = "1px"))) %>% 
  leafem::addHomeButton(extent(lfp4326), 'Ver todo')


## Exportar a KML

execGRASS(
  'v.out.ogr',
  flags = c('overwrite','quiet'),
  parameters = list(
    input = 'LfpNetwork_lfp_all_final_soco',
    output = 'lfp_kml.kml',
    format = 'KML',
    dsco = 'NameField=cat'
  )
)

## Obtención de perfiles longitudinales e índices de concavidad


devtools::source_url('https://raw.githubusercontent.com/geomorfologia-master/unidad-4-asignacion-1-procesos-fluviales/master/lfp_profiles_concavity.R')
soco_conv_prof <- LfpProfilesConcavity(
  xycoords = my_trans(c(-69.19744, 18.45272)),
  network = 'LfpNetwork_lfp_all_final_soco',
  prefix = 'soco',
  dem = 'dem',
  direction = 'drainage-dir-de-rstr',
  crs = '+init=epsg:32619',
  smns = 0.5,
  nrow = 3)

## Mostrar resultados


soco_conv_prof$profiles
soco_conv_prof$concavityindex
soco_conv_prof$dimensionlessprofiles

##  Tabla dx / dy, tanto en metros como adimensional. Útiles para construir perfiles por cuenta propia

soco_conv_prof$lengthzdata %>% tibble::as.tibble()
soco_conv_prof$lengthzdatadmnls %>% tibble::as.tibble()

##Video 12: GRASS GIS desde R: Parámetros de cuenca con r.basin -----

gmeta()
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Convertir a números enteros la extensión y la resolución del DEM


library(raster)
rutadem <- 'datos-fuente/srtm_dem_cuenca_soco.tif'
rawextent <- extent(raster(rutadem))
rawextent
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/integerextent.R')
devtools::source_url('https://raw.githubusercontent.com/geofis/rgrass/master/xyvector.R')
newextent <- intext(e = rawextent, r = 90, type = 'inner')
newextent
gdalUtils::gdalwarp(
  srcfile = 'datos-fuente/srtm_dem_cuenca_soco.tif',
  dstfile = 'srtm_dem_cuenca_socoint.tif',
  te = xyvector(newextent),
  tr = c(90,90),
  r = 'bilinear',
  overwrite = T
)

## Importar a sesión de GRASS


rutademint <- 'srtm_dem_cuenca_socoint.tif'
execGRASS(
  "g.proj",
  flags = c('t','c'),
  georef=rutademint)
gmeta()
execGRASS(
  "r.in.gdal",
  flags='overwrite',
  parameters=list(
    input=rutademint,
    output="demint"
  )
)

execGRASS(
  "g.region",
  parameters=list(
    raster = "demint",
    align = "demint"
  )
)


gmeta()
execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Generar red de drenaje para obtener coordenada posteriormente


execGRASS(
  "r.stream.extract",
  flags = c('overwrite','quiet'),
  parameters = list(
    elevation = 'demint',
    threshold = 80,
    stream_raster = 'stream-de-rstr-soco',
    stream_vector = 'stream_de_rstr_soco'
  )
)

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Obtener coordenada

library(sp)
use_sp()
library(mapview)
netw <- spTransform(
  readVECT('stream_de_rstr_soco'),
  CRSobj = CRS("+init=epsg:4326"))

mapview(netw, col.regions = 'blue', legend = FALSE)

source('https://raw.githubusercontent.com/geomorfologia-master/unidad-4-asignacion-1-procesos-fluviales/master/lfp_profiles_concavity.R')
outlet <- as.integer(my_trans(c(-69.20249, 18.45693)))

## Ejecutar`r.basin`


pref <- 'rbasin_soco'
execGRASS(
  "r.basin",
  flag ='overwrite',
  parameters=list(
    map='demint',
    prefix=pref,
    coordinates=outlet,
    threshold =80,
    dir='salidas-rbasin/soco'
  )
)

execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Cargar los vectoriales transformados a EPSG:4326 para visualizar en leaflet

rbnetw <- spTransform(
  readVECT('rbasin_soco_demint_network'),
  CRSobj = CRS("+init=epsg:4326"))
rbnetw
rbmain <- spTransform(
  readVECT('rbasin_soco_demint_mainchannel'),
  CRSobj = CRS("+init=epsg:4326"))
rbmain
rbbasin <- spTransform(
  readVECT('rbasin_soco_demint_basin'),
  CRSobj = CRS("+init=epsg:4326"))
rbbasin

library(leaflet)
leaflet() %>%
  addProviderTiles(providers$Stamen.Terrain, group = 'terrain') %>%
  addPolylines(data = rbnetw, weight = 3, opacity = 0.7) %>% 
  addPolylines(data = rbmain, weight = 3, opacity = 0.7, color = 'red') %>% 
  addPolygons(data = rbbasin) %>% 
  leafem::addHomeButton(extent(rbbasin), 'Ver cuenca')

## Explorar los parámetros de cuenca


library(readr)
rbsocopar1 <- read_csv("salidas-rbasin/soco/rbasin_soco_demint_parametersT.csv")
rbsocopar1 %>% tibble::as_tibble()
rbsocopar2 <- read_csv(
  "salidas-rbasin/soco/rbasin_soco_demint_parameters.csv",
  skip=2, col_names = c('Parameter', 'Value'))
rbsocopar2 %>% print(n=Inf)


##Video1 13:GRASS GIS desde R: Curva e integral hipsométrica

## Imprimir lista de mapas ráster y vectoriales dentro en la región/localización activa


execGRASS(
  'g.list',
  flags = 't',
  parameters = list(
    type = c('raster', 'vector')
  )
)

## Representar cuencas


library(sp)
use_sp()
library(mapview)
bas2 <- readVECT('r_stream_basins_2')
bas3 <- readVECT('r_stream_basins_3')


## Curva e integral hipsométrica

source('https://raw.githubusercontent.com/geomorfologia-master/unidad-4-asignacion-1-procesos-fluviales/master/integral_hypsometric_curve.R')
HypsoBasinsOrder2 <- HypsoIntCurve(
  basins = 'r_stream_basins_2',
  dem = 'dem',
  labelfield = 'cat',
  nrow = 5,
  labelsize = 3)

HypsoBasinsOrder2$HypsoInt
HypsoBasinsOrder2$HypsoCurve
mapview(bas2, zcol='cat', col.regions = 'blue', legend = FALSE) %>%
  addStaticLabels(label = bas2$cat)

HypsoBasinsOrder3 <- HypsoIntCurve(
  basins = 'r_stream_basins_3',
  dem = 'dem',
  labelfield = 'cat',
  nrow = 1,
  labelsize = 4
)


HypsoBasinsOrder3$HypsoInt
HypsoBasinsOrder3$HypsoCurve
mapview(bas3, zcol='cat', col.regions = 'blue', legend = FALSE) %>%
  addStaticLabels(label = bas3$cat)