library(tidyverse)
library(sf)
library(lwgeom)
library(nngeo)
library(sfnetworks)
library(tidygraph)
library(gtfstools)
library(mapview)
library(patchwork)
library(job)
Sys.setlocale(locale = "hebrew")


spo_gtfs <- read_gtfs("C:/Users/yehuda/Downloads/gtfs.zip",encoding = "UTF-8")
bs_pol <- spo_gtfs$routes %>% 
  # at first, only for Dan Beer sheva
  filter(agency_id == 32) %>% 
  select(route_id,route_desc) %>% 
  # join with trips
  left_join(spo_gtfs$trips, by = "route_id") %>% 
  left_join(spo_gtfs$stop_times, by = "trip_id") %>% 
  distinct(stop_id) %>% 
  left_join(spo_gtfs$stops, by = "stop_id") %>% 
  st_as_sf(coords = c("stop_lon","stop_lat"),crs = 4326) %>% 
  summarise() %>% 
  st_convex_hull() 

bs_gtfs <- filter_by_sf(spo_gtfs,bs_pol)
bs_gtfs$shapes <- bs_gtfs$shapes %>% 
  group_by(shape_id) %>% 
  arrange(shape_id,shape_pt_sequence)%>% 
  mutate(shape_pt_sequence = row_number()) %>% 
  ungroup() %>% 
  data.table::as.data.table()
# TODO document and add names
# prepareing segments from shapes
canon_segments1 <- bs_gtfs$shapes %>% 
  group_by(shape_id) %>% 
  arrange(shape_id,shape_pt_sequence) %>% 
  mutate(shape_pt_lat_prev = lag(shape_pt_lat),
         shape_pt_lon_prev = lag(shape_pt_lon)) %>% 
  ungroup() 
# finding all distinct segments
canon_segments2 <- canon_segments1 %>% 
  select(-shape_id,-shape_pt_sequence) %>% 
  filter(!is.na(shape_pt_lat_prev)) %>% 
  distinct() %>% 
  mutate(fid = row_number())
# joining both
canon_segments3 <- canon_segments1 %>% 
  left_join(canon_segments2, by = c("shape_pt_lat","shape_pt_lon","shape_pt_lat_prev","shape_pt_lon_prev"))
# distinct segments as geometry
q2 <- canon_segments2 %>% 
  mutate(line = pmap(list(shape_pt_lon_prev,shape_pt_lat_prev,shape_pt_lon,shape_pt_lat),
                     function(x1,y1,x2,y2){
                       st_linestring(matrix(c(x1,y1,x2,y2),ncol=2,byrow = T))
                     }) %>% st_sfc(crs = 4326)
  ) %>% 
  st_sf() %>% 
  st_transform(2039) %>% 
  select(fid) %>% 
  mutate(buf = st_buffer(line,-10, singleSide = T))
q3 <- bs_gtfs$stops %>% 
  st_as_sf(coords = c("stop_lon","stop_lat"),crs = 4326) %>% 
  st_transform(2039) %>% 
  select(stop_id)
q31 <- q3 %>% 
  select() %>% 
  distinct() %>% 
  mutate(canon_stop = row_number()) %>% 
  st_join(q3)
q32 <- q31 %>% 
  select(canon_stop) %>% 
  mutate(canon_stop = as.character(canon_stop)) %>% 
  distinct() 
# q4 <- st_join(q32,q2,st_is_within_distance,10)
q4 <- st_join(q32,q2 %>% st_set_geometry("buf"))
nr = 100
q5 <- q4 %>% filter(is.na(fid))
q6 <- q5 %>% slice(0)
i <- 11
while(nr > 0){
  q6 <- bind_rows(q6,q5 %>% 
                    st_join(q2 %>% select(-buf) ,st_is_within_distance,i) %>% 
                    select(-fid.x) %>% 
                    rename(fid = fid.y) %>% 
                    filter(!is.na(fid)))
  q5 <- q5 %>% 
    st_join(q2%>% select(-buf),st_is_within_distance,i) %>%
    select(-fid.x) %>% 
    rename(fid = fid.y) %>% 
    filter(is.na(fid))
  nr <- q5  %>% nrow()
  
  i <-  i + 1
  print(i)
}
q7 <- bind_rows(q4,q6) %>% select(-line)


shapes <- canon_segments3$shape_id %>% unique()
qqq <- map_df(shapes,function(x){
  q <- canon_segments3 %>% 
    # filter certain shape_id with fid of canon segments
    filter(shape_id == x) %>% 
    # join to segments geometry
    left_join(q2, by = "fid") %>% 
    # filter those without fid - the first row
    filter(!is.na(fid)) %>% 
    # join segments and canon stops
    left_join(q7 ,by = "fid") %>%
    # turn to table
    st_sf()
  print(unique(q$shape_id))
  data.frame(shape_id = unique(q$shape_id),max= max(q$shape_pt_sequence),nr = nrow(q))})
# two problematic routes
qqq %>% filter(max != nr + 1)
shapes1 <- shapes[!shapes %in%  c("70893","72568")]






segmentize <- function(x){
  newtst1 <- canon_segments3 %>% 
    # filter certain shape_id with fid of canon segments
    filter(shape_id == x) %>% 
    # join to segments geometry
    left_join(q2, by = "fid") %>% 
    # filter those without fid - the first row
    filter(!is.na(fid)) %>% 
    # join segments and canon stops
    left_join(q7 ,by = "fid") %>%
    # turn to table
    st_sf() 
  # not taking into account recurring stop
  
  
  newtst2 <- newtst1 %>% 
    mutate(dist = map2_dbl(geometry,line,~st_distance(.x,.y))) %>% 
    group_by(canon_stop) %>%
    mutate(canon_stop = ifelse(dist == min(dist), canon_stop,NA)) %>% 
    ungroup() %>% 
    arrange(canon_stop) %>% 
    mutate(geometry = ifelse(is.na(canon_stop),NA, geometry) %>% st_sfc(crs = 2039)) %>% 
    arrange(shape_pt_sequence) %>% 
    select(shape_id,shape_pt_sequence,fid,canon_stop,geometry) %>% 
    mutate(perpendicular = map2(geometry,line,function(x,y){
      st_nearest_points(x,y) %>% `[[`(1)
    })%>% st_sfc(crs = 2039),
    start = st_startpoint(line),
    end = st_endpoint(line),
    # should take into account 
    splitter = map(perpendicular,function(x){
      st_endpoint(x) %>%`[[`(1)
    })%>% st_sfc(crs = 2039),
    split1 = map2(start,splitter,function(x,y){
      st_union(x,y) %>% st_cast("LINESTRING")
    }) %>% st_sfc(crs = 2039),
    split2 = map2(end,splitter,function(x,y){
      st_union(y,x) %>% st_cast("LINESTRING")
    }) %>% st_sfc(crs = 2039)
    )
  tst_cuts <- newtst2 %>% 
    filter(!is.na(canon_stop)) %>% 
    st_set_geometry("splitter") %>% 
    select(shape_id,canon_stop,shape_pt_sequence)
  tst_lines_first <- newtst2 %>% 
    filter(!is.na(canon_stop)) %>% 
    st_set_geometry("split1") %>% 
    select(shape_id,canon_stop,shape_pt_sequence) %>% 
    mutate(geom  = lead(split1),
           type = "second") %>% 
    st_set_geometry("geom") %>% 
    select(-split1)
  tst_lines_second <- newtst2 %>% 
    filter(!is.na(canon_stop)) %>% 
    st_set_geometry("split2") %>% 
    select(shape_id,canon_stop,shape_pt_sequence) %>% 
    mutate(geom  = split2,
           type = "first") %>% 
    st_set_geometry("geom") %>% 
    select(-split2)
  
  tst_lines_regular <- newtst2 %>% 
    filter(is.na(canon_stop)) %>% 
    select(shape_id,canon_stop,shape_pt_sequence) %>% 
    mutate(geom  = line,
           type = "regular") %>% 
    st_set_geometry("geom") %>% 
    select(-line)
  
  bind_rows(tst_lines_first,tst_lines_second,tst_lines_regular) %>% 
    filter(!st_is_empty(geom)) %>% 
    arrange(shape_pt_sequence) %>% 
    group_by(canon_stop) %>% 
    mutate(rn = ifelse(!is.na(canon_stop),row_number(),NA)) %>%
    ungroup() %>% 
    fill(canon_stop) %>% 
    fill(rn) %>% 
    mutate(type = factor(type,levels = c("first","regular","second")),
           canon_stop = factor(canon_stop,levels = unique(canon_stop)),
           oper_canon_stop = paste0(canon_stop,"_",ifelse(rn%%2 == 0,rn-1,rn)),
           oper_canon_stop = factor(oper_canon_stop,levels = unique(oper_canon_stop))) %>% 
    arrange(oper_canon_stop,type,shape_pt_sequence) %>% 
    mutate(mat = map2(geom,type,~if(.y != "second"){st_coordinates(.x)[1,1:2] %>% as.matrix() %>% t()}else{st_coordinates(.x)[,1:2] %>% as.matrix()})) %>% 
    st_drop_geometry() %>% 
    group_by(oper_canon_stop) %>% 
    summarise(tst = list(do.call(rbind,mat)) ) %>% 
    mutate(tst1 = st_sfc(map(tst,st_linestring),crs = 2039)) %>% 
    select(-tst) %>% 
    st_sf() %>%  
    # to be fixed when adding start points
    filter(oper_canon_stop != "NA_NA") %>% 
    # to be fixed - make sure no single point lines exist
    mutate(npts = map_int(tst1,npts)) %>% 
    filter(npts > 1) %>% 
    as_sfnetwork() %E>% st_as_sf() 
  # mapview(zcol = "oper_canon_stop") + tst_cuts + stops

}
i = 1
all_routes <- map(shapes1,function(x){
  i<<-i+1
  options(warn = -1)
  print(x)
  tim <- Sys.time()
  res <- segmentize(x)
  print(Sys.time() - tim)
  return(res)
})



newtst1 <- canon_segments3 %>% 
  # filter certain shape_id with fid of canon segments
  filter(shape_id == "112677") %>% 
  # join to segments geometry
  left_join(q2, by = "fid") %>% 
  # filter those without fid - the first row
  filter(!is.na(fid)) %>% 
  # join segments and canon stops
  left_join(q7 ,by = "fid") %>%
  # turn to table
  st_sf() 
# not taking into account recurring stop


newtst2 <- newtst1 %>% 
  mutate(dist = map2_dbl(geometry,line,~st_distance(.x,.y))) %>% 
  group_by(canon_stop) %>%
  mutate(canon_stop = ifelse(dist == min(dist), canon_stop,NA)) %>% 
  ungroup() %>% 
  arrange(canon_stop) %>% 
  mutate(geometry = ifelse(is.na(canon_stop),NA, geometry) %>% st_sfc(crs = 2039)) %>% 
  arrange(shape_pt_sequence) %>% 
  select(shape_id,shape_pt_sequence,fid,canon_stop,geometry) %>% 
  mutate(perpendicular = map2(geometry,line,function(x,y){
    st_nearest_points(x,y) %>% `[[`(1)
  })%>% st_sfc(crs = 2039),
  start = st_startpoint(line),
  end = st_endpoint(line),
  splitter = map(perpendicular,function(x){
  st_endpoint(x) %>%`[[`(1)
  })%>% st_sfc(crs = 2039),
  dist1 = map2_dbl(start,splitter,st_distance),
  dist2 = map2_dbl(end,splitter,st_distance)) %>% {
    a <- (.) %>% st_set_geometry("geometry")
    b <- (.) %>% st_set_geometry("splitter")
    mapview((.),zcol = "shape_pt_sequence") + mapview(a,zcol = "canon_stop") + mapview(b,zcol = "canon_stop")
  }
  
  # ,
  # should take into account 
  # splitter = map(perpendicular,function(x){
    # st_endpoint(x) %>%`[[`(1)
  # })%>% st_sfc(crs = 2039),
  split1 = map2(start,splitter,function(x,y){
    st_union(x,y) %>% st_cast("LINESTRING")
  }) %>% st_sfc(crs = 2039),
  split2 = map2(end,splitter,function(x,y){
    st_union(y,x) %>% st_cast("LINESTRING")
  }) %>% st_sfc(crs = 2039)
  )
tst_cuts <- newtst2 %>% 
  filter(!is.na(canon_stop)) %>% 
  st_set_geometry("splitter") %>% 
  select(shape_id,canon_stop,shape_pt_sequence)
tst_lines_first <- newtst2 %>% 
  filter(!is.na(canon_stop)) %>% 
  st_set_geometry("split1") %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = lead(split1),
         type = "second") %>% 
  st_set_geometry("geom") %>% 
  select(-split1)
tst_lines_second <- newtst2 %>% 
  filter(!is.na(canon_stop)) %>% 
  st_set_geometry("split2") %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = split2,
         type = "first") %>% 
  st_set_geometry("geom") %>% 
  select(-split2)

tst_lines_regular <- newtst2 %>% 
  filter(is.na(canon_stop)) %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = line,
         type = "regular") %>% 
  st_set_geometry("geom") %>% 
  select(-line)

bind_rows(tst_lines_first,tst_lines_second,tst_lines_regular) %>% 
  filter(!st_is_empty(geom)) %>% 
  arrange(shape_pt_sequence) %>% 
  group_by(canon_stop) %>% 
  mutate(rn = ifelse(!is.na(canon_stop),row_number(),NA)) %>%
  ungroup() %>% 
  fill(canon_stop) %>% 
  fill(rn) %>% 
  mutate(type = factor(type,levels = c("first","regular","second")),
         canon_stop = factor(canon_stop,levels = unique(canon_stop)),
         oper_canon_stop = paste0(canon_stop,"_",ifelse(rn%%2 == 0,rn-1,rn)),
         oper_canon_stop = factor(oper_canon_stop,levels = unique(oper_canon_stop))) %>% 
  arrange(oper_canon_stop,type,shape_pt_sequence) %>% 
  mutate(mat = map2(geom,type,~if(.y != "second"){st_coordinates(.x)[1,1:2] %>% as.matrix() %>% t()}else{st_coordinates(.x)[,1:2] %>% as.matrix()})) %>% 
  st_drop_geometry() %>% 
  group_by(oper_canon_stop) %>% 
  summarise(tst = list(do.call(rbind,mat))) %>% 
  mutate(tst1 = st_sfc(map(tst,st_linestring),crs = 2039)) %>% 
  select(-tst) %>% 
  st_sf() %>%  
  # to be fixed when adding start points
  filter(oper_canon_stop != "NA_NA") %>% 
  mutate(npts = map_int(tst1,npts)) %>% 
  filter(npts > 1) %>% 
  as_sfnetwork() %E>% st_as_sf() %>% 
  mapview(zcol = "oper_canon_stop") + tst_cuts + stops




canon_segments3 %>% 
  # filter certain shape_id with fid of canon segments
  filter(shape_id == "112677") %>% 
  left_join(q2, by = "fid") %>% 
  st_sf() %>% 
  filter(!is.na(fid)) %>% 
  arrange(shape_pt_sequence) %>% 
  summarise() %>% 
  st_line_merge() %>% 
  st_buffer(10,singleSide = T) %>% 
  mapview() + mapview(q2[q2$fid %in% canon_segments3[canon_segments3$shape_id == "112677",]$fid, ])

  st_sfc(crs=4326) %>% 
  st_transform(2039) %>% 
  %>% 
  {
    mapview((.)) + st_intersection(q7,(.))
  } 

  summarise(line = map2(shape_pt_lon,shape_pt_lat,~st_linestring(matrix(c(.x,.y),byrow = T,ncol = 2)))) 
  # join to segments geometry
  left_join(q2, by = "fid") %>% 
  