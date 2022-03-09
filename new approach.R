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
# TODO if this is fucked, go back to prod1
bs_lines <- get_trip_geometry(bs_gtfs) %>% 
  left_join(bs_gtfs$trips, by = "trip_id") %>% 
  select(shape_id) %>%
  group_by(shape_id) %>% 
  mutate(rn = row_number()) %>% 
  filter(rn == 1) %>%
  ungroup() %>% 
  select(-rn) %>% 
  st_sf() %>% 
  st_transform(2039) 
intersections <- 
  # explode all lines
  st_segments(bs_lines$geometry) %>% 
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

intersections_ready <- intersections %>% 
  mutate(own_id = paste0("inter_",.tidygraph_node_index)) %>% 
  rename(geometry =x) %>% 
  select(own_id)
start_points <- bs_lines %>% 
  mutate(geometry = st_startpoint(geometry),
         own_id = paste0("start_",shape_id)) %>% 
  select(-shape_id)
end_points <-  bs_lines %>% 
  mutate(geometry = st_endpoint(geometry),
         own_id = paste0("end_",shape_id)) %>% 
  select(-shape_id)
# TODO think about how these should be integrated
supple_geom <- bind_rows(start_points,end_points,intersections_ready) %>% 
  select() %>% 
  distinct() %>% 
  mutate(suppl_id = paste0("supple_",row_number()))
supple_to_own_id <- supple_geom %>% 
  st_join(bind_rows(start_points,end_points,intersections_ready))

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
# to be used in splitter context
# q41 <- st_join(supple_geom,q2,st_is_within_distance,0.1)
q41 <- st_join(supple_geom,q2 %>% st_set_geometry("buf")) %>% 
  rename(canon_stop = suppl_id)

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
# q7 <- bind_rows(q4,q6,q41) %>% select(-line)



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
qqq %>% filter(max != nr + 1)
bs_gtfs$trips %>% filter(shape_id == "72568")
# test a certain shape_id
tst1 <- canon_segments3 %>% 
  # filter certain shape_id with fid of canon segments
  filter(shape_id == "121612") %>% 
  # join to segments geometry
  left_join(q2, by = "fid") %>% 
  # filter those without fid - the first row
  filter(!is.na(fid)) %>% 
  # join segments and canon stops
  left_join(q7 ,by = "fid") %>%
  # turn to table
  st_sf() 
# manipulation in the stops level - each stop should appear only once in each 
tst2 <- tst1 %>% 
  # select relavant columns
  select(shape_id,shape_pt_sequence,fid,canon_stop,geometry) %>%
  # why do we get more than 1 canon stop?
  # because many segments are close enough. the solution is defining only those closest as connected with stop. 
  # create distance from each stop to each segment, in order to clarify one cut for each stop
  mutate(dist = map2_dbl(geometry,line,~st_distance(.x,.y))) %>% 
  group_by(canon_stop) %>%
  mutate(canon_stop = ifelse(dist == min(dist), canon_stop,NA)) %>% 
  # done to prevent bug of sf
  ungroup() %>% 
  arrange(canon_stop) %>% 
  mutate(geometry = ifelse(is.na(canon_stop),NA, geometry) %>% st_sfc(crs = 2039)) %>% 
  arrange(shape_pt_sequence)

# why do we get 2 or more shape sequence id?
# two stops or more are connected to one segemnt. 
# cases: 
# 1. 2 or more have canon stop, because there 2 or more stops which are closest to the segment. needs 2 ot more cuts. TODO
# 2. only one has canon stop. delete the other row
# 3. none has canon stop.

tst3 <- tst2 %>% 
  # solves 3
  select(-dist) %>% 
  distinct() %>% 
  # solves 2
  arrange(shape_pt_sequence) %>% 
  group_by(shape_pt_sequence) %>% 
  mutate(n = n(), 
         filled = sum(!is.na(canon_stop)),
         rn = row_number()) %>% 
  filter(!(filled ==1& n >1& rn > 1)) %>% 
  select(-n,-filled,-rn) %>% 
  ungroup()
tst4 <- tst3 %>% 
  mutate(shape_pt_sequence = sprintf("%05d", shape_pt_sequence)) %>% 
  group_by(shape_pt_sequence) %>% 
  mutate(rn = row_number(),
         new_shape_pt_sequence = paste0(shape_pt_sequence,"_",rn))
tst5 <- tst4 %>% 
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

tst51 <- tst5 %>% 
  filter(n() == 1)
# solves 1  
tst52 <- tst5 %>% 
  filter(n() > 1) %>% 
  mutate(lag = lag(splitter),
         lead = lead(splitter),
         split1 = pmap(list(rn,split1,lag,splitter),function(x1,x2,x3,x4){
           if(x1 == 1){
             x2 
           }else{
             st_union(x3,x4) %>% st_cast("LINESTRING")
           }
         }) %>% st_sfc(crs = 2039),
         split2 = pmap(list(rn,split2,splitter,lead),function(x1,x2,x3,x4){
           if(x1 == 2){
             x2 
           }else{
             st_union(x3,x4) %>% st_cast("LINESTRING")
           }
         }) %>% st_sfc(crs = 2039)) 
tst6 <- bind_rows(tst51,tst52) %>% arrange(new_shape_pt_sequence)
tst6 %>% mapview() + mapview(st_set_geometry(tst6,"geometry"))
  filter(!is.na(canon_stop)) %>% 
  st_set_geometry("split1") %>% 
  mapview()
newtst1 <- canon_segments3 %>% 
  # filter certain shape_id with fid of canon segments
  filter(shape_id == "121612") %>% 
  # join to segments geometry
  left_join(q2, by = "fid") %>% 
  # filter those without fid - the first row
  filter(!is.na(fid)) %>% 
  # join segments and canon stops
  left_join(q7 ,by = "fid") %>%
  # turn to table
  st_sf() 


newtst2 <- newtst1 %>% 
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
  }) %>% st_sfc(crs = 2039),
  split2lead = lead(split2)
  )
  
  


# TOOO make sure there is only one shape_pt_sequence
# TODO continue from here
mapview(tst2,zcol = "shape_pt_sequence")

  #  TODO this should be changed conceptually. we filter out stops that should not be filtered. 
  # I should implemnet a method that enables more than one splitter per segment
  # if multiple segemnts are joined to one canon stop - choose only closestm
  # else - create local sf network
  group_by(fid) %>% 
  filter(dist == min(dist) | is.na(dist)) %>% 
  ungroup() %>% 
  
tst21 <- tst1 %>% filter(canon_stop %>% is.na() %>% `!`())
tst22 <- tst1 %>% filter(canon_stop %>% is.na())


tst3 <- tst21 %>% 
  mutate(perpendicular = map2(geometry,line,function(x,y){
    st_nearest_points(x,y) %>% `[[`(1)
  })%>% st_sfc(crs = 2039),
  # should take into account 
  splitter = map(perpendicular,function(x){
    st_endpoint(x) %>%`[[`(1)
  })%>% st_sfc(crs = 2039),
  split1 = map2(line,splitter,function(x,y){
    st_union(st_startpoint(x),y) %>% st_cast("LINESTRING") %>%`[[`(1)
  }) %>% st_sfc(crs = 2039),
  split1lead = lead(split1),
  split2 = map2(line,splitter,function(x,y){
    st_union(y,st_endpoint(x)) %>% st_cast("LINESTRING")
  }) %>% st_sfc(crs = 2039),
  split2lead = lead(split2)
  )
tst_cuts <- tst3 %>% 
  st_set_geometry("splitter") %>% 
  select(shape_id,canon_stop,shape_pt_sequence)
tst_lines_second <- tst3 %>% 
  st_set_geometry("split1") %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = split1,
         type = "second") %>% 
  st_set_geometry("geom") %>% 
  select(-split1)
tst_lines_first <- tst3 %>% 
  st_set_geometry("split2lead") %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = split2lead,
         type = "first") %>% 
  st_set_geometry("geom") %>% 
  select(-split2lead)
tst_lines_regular <- tst22 %>% 
  select(shape_id,canon_stop,shape_pt_sequence) %>% 
  mutate(geom  = line,
         type = "regular") %>% 
  st_set_geometry("geom") %>% 
  select(-line)
bind_rows(tst_lines_first,tst_lines_second,tst_lines_regular) %>% 
  # filter(!st_is_empty(geom)) %>% 
  arrange(shape_pt_sequence) %>%
  fill(canon_stop) %>% 
  mutate(type = factor(type,levels = c("first","regular","second"))) %>% 
  arrange(canon_stop,type,shape_pt_sequence) %>% 
  group_by(canon_stop) %>% 
  summarise() %>% 
  mapview(zcol = "canon_stop")+tst_cuts



sfnet1 <- canon_segments3 %>% 
  filter(shape_id == "114341") %>% 
  left_join(q2, by = "fid") %>% 
  filter(!is.na(fid)) %>% 
  left_join(q7 ,by = "fid") %>%
  st_sf() %>% 
  select(shape_id,shape_pt_sequence,fid,canon_stop,geometry) %>%
  mutate(dist = map2_dbl(geometry,line,~st_distance(.x,.y))) %>% 
  group_by(canon_stop) %>%
  mutate(canon_stop = ifelse(dist == min(dist), canon_stop,NA)) %>% 
  ungroup() %>% 
  mutate(geometry = ifelse(!is.na(canon_stop), geometry,NA) %>% st_sfc(crs = 2039))  %>% 
  filter(canon_stop %>% is.na() %>% `!`()) %>% 
  mutate(perpendicular = map2(geometry,line,function(x,y){
    st_nearest_points(x,y) %>% `[[`(1)
  })%>% st_sfc(crs = 2039)) %>% 
  st_set_geometry("geometry") %>% 
  st_drop_geometry() %>% 
  gather(type,geom, line,perpendicular) %>% 
  mutate(geom = st_sfc(geom,crs = 2039)) %>% 
  st_sf() %>% 
  group_by(fid) %>% 
  nest() %>% 
  mutate(sfn = map(data,~as_sfnetwork(.x) %>% convert(to_spatial_subdivision) #%>% 
                     # filter(centrality_degree() > 1)
                   ),
         nds = map(sfn,st_as_sf)) #%>% 
  # ungroup() %>% 
  # unnest(nds) %>%
  # st_sf()
  #  TODO this should be changed conceptually. we filter out stops that should not be filtered. 
  # I should implemnet a method that enables more than one splitter per segment
  # if multiple segemnts are joined to one canon stop - choose only closestm
  # else - create local sf network
  group_by(fid) %>% 
  filter(dist == min(dist) | is.na(dist)) %>% 
  ungroup() 
  
  
  canon_segments3 %>% 
  filter(shape_id == "114341") %>% 
  left_join(q2, by = "fid") %>% 
  filter(!is.na(fid)) %>% 
  left_join(q7 ,by = "fid") %>%
  st_sf() %>% 
  select(shape_id,shape_pt_sequence,fid,canon_stop,geometry) %>% 
  mutate(dist = map2_dbl(geometry,line,~st_distance(.x,.y))) %>% 
  group_by(canon_stop) %>% 
  mutate(canon_stop = ifelse(dist == min(dist), canon_stop,NA)) 
  
  
  filter(dist == min(dist) | is.na(dist)) %>% filter(fid == 4872)
  ungroup() %>% 
  filter(canon_stop %>% is.na() %>% `!`()) %>% 
  
  mutate(perpendicular = map2(geometry,line,function(x,y){
    st_nearest_points(x,y) %>% `[[`(1)
  })%>% st_sfc(crs = 2039))

sfnet1 %>% 
  filter(canon_stop == 1372) 
sfnet1 %>% 
  st_set_geometry("geometry") %>% 
  st_drop_geometry() %>% 
  gather(type,geom, line,perpendicular) %>% 
  mutate(geom = st_sfc(geom,crs = 2039)) %>% 
  st_sf() %>% 
  group_by(fid) %>% 
  nest() %>% 
  mutate(sfn = map(data,~as_sfnetwork(.x) %>% convert(to_spatial_subdivision)),
         nds = map(sfn,st_as_sf)) %>% 
  # unnest(nds) %>%
   %>% 
  pull(sfn) %>%
  `[[`(1) %>%
  autoplot()
  st_sf() %>% 
  mapview()

# # ~st_nearest_points(.x,.y)%>% `[[`(1)) %>% st_sfc(crs = 2039),
# # splitter = map(perpendicular,~st_endpoint(.x) %>%`[[`(1))%>% st_sfc(crs = 2039),
# # split1 = map2(line,splitter,~st_union(st_startpoint(.x),.y) %>% st_cast("LINESTRING") %>%`[[`(1))%>%st_sfc(crs = 2039),
# # split1lead = lead(split1))#,
# # split2 = map2(line,splitter,~st_union(.y,st_endpoint(.x)) %>% st_cast("LINESTRING")) %>% st_sfc(crs = 2039)) %>% 
# tst3 <- tst2 %>% 
#   # select(-dist,-perpendicular,-splitter) %>%
#   fill(canon_stop) %>% 
#   group_by(canon_stop) %>% 
#   mutate(lrn = row_number(),
#          ismax = lrn == max(lrn)) #%>% 
# # filter(lrn == max(lrn)) %>% 
# ungroup() %>% 
#   mutate(megred_ending = pmap(list(line,split1lead,ismax),function(.x,.y,.z){
#     # if(!.z){
#     #   NA
#     # }
#     if(length(.y) == 0){
#       .x
#     }else{
#       st_union(.x,.y) %>% st_line_merge()
#     }
#   }
#   )%>% st_sfc(crs = 2039))
# 
# tsts1 <- tst %>% 
#   filter(lrn == max(lrn)) 
# tst3 %>% 
#   st_set_geometry("perpendicular") %>% 
#   mapview(color  = "red")  + 
#   mapview(tst3$line) + 
#   mapview(tst3$geometry,col.regions = "blue") + 
#   mapview(tst3$megred_ending,color = "black")
# mapview(tst$split1lead)
# mapview(tst) + q7[q7$fid %in% tst$fid,]
# todor::todor()
