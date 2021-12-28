library(tidyverse)
library(sf)
library(lwgeom)
library(sfnetworks)
library(tidygraph)
library(gtfstools)
library(mapview)
library(nngeo)
# downloaded gtfs from 7.7.2021
Sys.setlocale(locale = "hebrew")
# load gtfs
spo_gtfs <- read_gtfs("D:/Downloads/gtfs.zip",encoding = "UTF-8")
# create lines
lines <- spo_gtfs$routes %>% 
  # at first, only for Dan Beer sheva
  filter(agency_id == 32) %>% 
  select(route_id,route_desc) %>% 
  # join with trips
  left_join(spo_gtfs$trips, by = "route_id") %>% 
  # get only distinct values
  select(route_desc,shape_id) %>% 
  distinct() %>% 
  # join with shapes
  left_join(spo_gtfs$shapes, by = "shape_id") %>% 
  # arrange points in each line by sequence 
  arrange(route_desc,shape_id,shape_pt_sequence) %>% 
  # create linestring for each line
  group_by(route_desc,shape_id) %>% 
  nest() %>% 
  mutate(line = st_sfc(map(data, ~.x %>% 
                      select(shape_pt_lon ,shape_pt_lat) %>%
                      as.matrix() %>% 
                      st_linestring() ),crs=4326),
         # transform to ITM
         line = st_transform(line,2039)) %>% 
  select(-data) %>% 
  ungroup() %>% 
  st_sf()
# better intersections
# segs <- st_segments(lines$line) %>% st_as_sf()
# segs %>% 
#   distinct() %>% 
#   as_sfnetwork() %>% 
#   convert(to_spatial_smooth) %>% 
#   st_as_sf() 

# find all intersections between lines
intersections <- 
  # explode all lines
  st_segments(lines$line) %>% 
  st_as_sf() %>% 
  # get only distinct segments
  distinct() %>% 
  # parse as sfnetoworks object
  as_sfnetwork() %>% 
  # subdivide segments that intersect (assuming planar graph)
  convert(to_spatial_subdivision) %>% 
  # smooth useless nodes
  convert(to_spatial_smooth) %>% 
  st_as_sf() 
# mapview(intersections) + lines
# qwer <- as_sfnetwork(lines) %>% 
  # convert(to_spatial_subdivision) %>% 
  # convert(to_spatial_smooth)
# qwer %>% st_as_sf() %>% mapview()
# points <- st_difference(lines,lines) %>% 
#   st_cast("MULTILINESTRING") %>% 
#   st_cast("LINESTRING") %>% 
#   mutate(start = st_startpoint(line),
#          end = st_endpoint(line)) %>% 
#   st_drop_geometry() %>% 
#   select(start,end) %>% 
#   gather(type,geom)
# intersections <- distinct(points %>% select(geom) ) %>% st_sf()

# get all relevant stops, by line makat
stops_routes <- spo_gtfs$routes %>% 
  # at first, only for Dan Beer sheva
  filter(agency_id == 32) %>% 
  # join change until getting to stops
  left_join(spo_gtfs$trips, by = "route_id") %>% 
  left_join(spo_gtfs$stop_times, by = "trip_id") %>% 
  left_join(spo_gtfs$stops,by = "stop_id") %>% 
  # get only distinct stops per line
  select(route_id,route_desc,shape_id,stop_id,stop_sequence,stop_lat,stop_lon)  %>% 
  distinct() %>% 
  st_as_sf(coords = c("stop_lon","stop_lat"),crs = 4326) %>% 
  st_transform(2039)
# line_ends <- lines %>% 
#   mutate(start = st_startpoint(line),
#          end = st_endpoint(line)) %>% 
#   st_drop_geometry() %>% 
#   select(start,end) %>% 
#   gather(type,geom)%>% 
#   select(geom) %>% 
#   distinct() %>% 
#   st_sf(crs = 2039) 
# line_ends %>% mapview()

# demonstration for one line
# get all makats
makats <- lines$route_desc
# use first makat
map(makats[1],function(x){
  line <- lines %>% filter(route_desc == x)
  stops <- stops_routes %>% filter(route_desc == x)
  line %>% 
    as_sfnetwork() %>% 
    st_network_blend(stops,10) %>% 
    st_network_blend(stops_routes,10) %>% 
    st_network_blend(intersections,0.1)
}) %>% 
  `[[`(1) -> qqq
qqq %>% activate(nodes) %>% st_as_sf() -> q
qqq %>% activate(edges) %>% st_as_sf() -> w
mapview(q) + w + intersections + stops_routes
mapview(lines) + intersections
lines %>% filter(route_desc == makats[1]) %>% mapview() + stops_routes
