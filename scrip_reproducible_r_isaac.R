#video_3.

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

#video_4

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


#Video 5 Explorar datos espaciales básicos entre GRASS y R


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
