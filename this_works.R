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
bs_stops <- bs_gtfs$stops %>% 
  st_as_sf(coords = c("stop_lon","stop_lat"),crs = 4326) %>% 
  st_transform(2039) %>% 
  select(stop_id)
intersections <- 
  # explode all lines
  st_segments(bs_lines$geometry) %>% 
  st_as_sf() %>% 
  # get only distinct segments
  distinct()  %>% 
  as_sfnetwork(directed = F) %>% 
  # subdivide segments that intersect (assuming planar graph)
  convert(to_spatial_subdivision) %>% 
  # smooth useless nodes
  convert(to_spatial_smooth) %>%
  mutate(inter_id = paste0("inter_",row_number()),
         degree = centrality_degree(mode = "all")) %>% 
  st_as_sf() 


joined <- bs_stops %>% 
  st_join(bs_lines %>% 
            st_buffer(-10,singleSide = T) ,left = FALSE) %>% 
  arrange(stop_id)
joined_inter <- intersections %>% 
  st_join(bs_lines %>% 
            st_buffer(-10,singleSide = T) ,left = FALSE)
res <- map(1:nrow(bs_lines),function(.x){
  print(bs_lines[.x,]$shape_id)
  one <- joined %>% filter(shape_id == bs_lines[.x,]$shape_id)
  two <- joined_inter%>% filter(shape_id == bs_lines[.x,]$shape_id)
  net <- bs_lines[.x,] %>% 
    as_sfnetwork(directed = T) %>% 
    convert(to_spatial_subdivision) %>%
    st_network_blend(two)
  rel_rows <- net %E>%
    st_as_sf() %>% 
    mutate(rn = row_number()) %>% 
    st_buffer(-10,singleSide = T) %>% 
    filter(apply(st_intersects(geometry,one,sparse = F),1,any)) 
  rel_idx <- net %E>%
    st_as_sf() %>% 
    mutate(rn = row_number()) %>% 
    st_buffer(-10,singleSide = T) %>% 
    st_intersects(one,sparse = F) %>% 
    apply( 1,which)
  
  nets <- net %E>%
    st_as_sf() %>% 
    mutate(rn = row_number(),
           idx = rel_idx) %>% 
    filter(rn %in% rel_rows$rn) %>% 
    mutate(net1 = map2(geometry,rn, ~.x %>% st_sfc(crs = 2039) %>% st_sf() %>% as_sfnetwork(directed = T) %E>% mutate(rn = .y)),
           net2 = map2(net1,idx, function(x,y){
             st_network_blend(x,one[y,]) %>% mutate(rn2 = row_number())
           })) %>% 
    pull(net2)
  stops <- bind_rows(map(nets,~.x %N>% st_as_sf())) %>% 
    filter(!is.na(stop_id)) %>% 
    distinct()
  final_net <- net %E>%
    mutate(rn = row_number()) %>%
    filter(!rn %in% rel_rows$rn) %E>%
    bind_graphs(bind_graphs(nets)) %>%
    as_sfnetwork(edges_as_lines = T,directed = T) %E>% 
    st_as_sf() %>%
    arrange(rn,rn2) %>% 
    as_sfnetwork(directed = T)
  route <- final_net %E>% 
    st_as_sf() %>% 
    mutate(cs = cumsum(st_length(geometry))) 
  cuts <- final_net %>% 
    st_join(stops) %>%
    st_join(two) %>% 
    mutate(rn = row_number()) %>% 
    as_tibble() %>% 
    as_tibble() %>%
    select(geometry,stop_id,inter_id,rn)
  # phase_1 <- route %>% 
  #   left_join(cuts,by = c("from" = "rn")) %>% 
  #   left_join(cuts,by = c("to" = "rn")) %>%
  #   fill(shape_id,.direction = "downup") %>% 
  #   mutate(origin = ifelse(!is.na(stop_id.x),stop_id.x,ifelse(!is.na(inter_id.x),inter_id.x,NA)),
  #          destination = ifelse(!is.na(stop_id.y),stop_id.y,ifelse(!is.na(inter_id.y),inter_id.y,NA)),
  #          origin = ifelse(row_number() == 1,paste0("start_",shape_id),origin),
  #          destination = ifelse(row_number() == nrow(.),paste0("end_",shape_id),destination)) %>% 
  #   select(origin,destination,from,to) %>% 
  #   fill(origin) %>% 
  #   fill(destination,.direction = "up") 
  # idx_list <- phase_1 %>% 
  #   group_by(origin,destination) %>% 
  #   group_indices() %>% unique()
  # phase_2 <- phase_1 %>%   
  #   group_by(origin,destination) %>% 
  #   mutate(gg = which(cur_group_id() == idx_list)) %>% 
  #   group_by(origin,destination,gg) %>% 
  #   summarise(do_union=T)
  # phase_21 <- phase_2 %>% 
  #   filter(st_geometry_type(geometry.x) == "MULTILINESTRING")
  # if(nrow(phase_21) == 0){
  #   phase_3 <- phase_2
  # }else{
  #   phase_3 <- bind_rows(phase_2 %>% 
  #                          filter(!st_geometry_type(geometry.x) == "MULTILINESTRING"),
  #                        phase_2 %>% 
  #                          filter(st_geometry_type(geometry.x) == "MULTILINESTRING") %>% 
  #                          st_line_merge())
  # }
  # phase_3 %>% 
  #   ungroup() %>% 
  #   arrange(gg) %>% 
  #   # st_buffer(-10,singleSide = T) %>%
  #   mutate(rn3 = row_number())
  # 
  phase_1 <- route %>%
    left_join(cuts,by = c("from" = "rn")) %>%
    left_join(cuts,by = c("to" = "rn")) %>%
    fill(shape_id,.direction = "downup") %>%
    mutate(origin = ifelse(!is.na(stop_id.x),stop_id.x,ifelse(!is.na(inter_id.x),inter_id.x,NA)),
           destination = ifelse(!is.na(stop_id.y),stop_id.y,ifelse(!is.na(inter_id.y),inter_id.y,NA)),
           origin = ifelse(row_number() == 1,paste0("start_",shape_id),origin),
           destination = ifelse(row_number() == nrow(.),paste0("end_",shape_id),destination)) %>%
    select(origin,destination,from,to) %>%
    fill(origin) %>%
    fill(destination,.direction = "up") %>% 
    mutate(lag1 = lag(origin),
           lag2 = lag(destination),
           counter =coalesce(ifelse(origin == lag1 & destination == lag2,0,1),1) %>% cumsum()) %>% 
    select(-lag1,-lag2)
  
  phase_2 <- phase_1 %>%
    group_by(origin,destination,counter) %>%
    summarise(do_union=T)
  phase_21 <- phase_2 %>%
    filter(st_geometry_type(geometry.x) == "MULTILINESTRING")
  if(nrow(phase_21) == 0){
    phase_3 <- phase_2
  }else{
    phase_3 <- bind_rows(phase_2 %>%
                           filter(!st_geometry_type(geometry.x) == "MULTILINESTRING"),
                         phase_2 %>%
                           filter(st_geometry_type(geometry.x) == "MULTILINESTRING") %>%
                           st_line_merge())
  }
  phase_3 %>%
    ungroup() %>%
    arrange(counter) %>%
    # st_buffer(-10,singleSide = T) %>%
    mutate(rn3 = row_number())
})
names(res) <- bs_lines$shape_id
with_shapes <- map(names(res),~res[[.x]] %>% mutate(shape_id = .x)) %>% bind_rows()
starts <- with_shapes %>% st_startpoint() %>% st_sf()
mapview(with_shapes) + mapview(starts)

with_shapes %>% 
  group_by(shape_id) %>% 
  mutate(lead1 = lead(origin)) %>% 
  filter(destination != lead1)



with_shapes %>% filter(shape_id == "124301",gg %in% 9:16) %>% ggplot()+ geom_sf() + facet_wrap(~gg)
stops_in_routes <- bs_gtfs$routes %>% 
  left_join(bs_gtfs$trips,by = "route_id") %>% 
  left_join(bs_gtfs$stop_times,by = "trip_id") %>% 
  select(route_id,route_desc,shape_id,stop_id,agency_id) %>% 
  distinct()  
# group_by(route_id,route_desc,shape_id) %>% 
# summarise(stop_id = list(stop_id))
# remove train
to_remove <- stops_in_routes %>% filter(agency_id== 2) %>% distinct(shape_id) %>% pull(shape_id)
with_shapes_2 <- with_shapes %>% filter(!shape_id %in% to_remove)

with_shapes_3 <- with_shapes_2 %>% 
  mutate(length = st_length(geometry.x)) %>% 
  as_tibble() %>% 
  group_by(shape_id) %>% 
  mutate(length= as.numeric(length),
         cs = replace_na(lag(cumsum(length)),0)) %>% 
  left_join(stops_in_routes,by = c("origin" = "stop_id","shape_id")) %>% 
  left_join(stops_in_routes,by = c("destination" = "stop_id","shape_id"))

makats <- bs_gtfs$routes %>% 
  filter(!agency_id== 2) %>% 
  left_join(bs_gtfs$trips,by = "route_id") %>% 
  select(route_id,route_desc,shape_id) %>% 
  distinct()
prod1 <- makats %>% 
  left_join(with_shapes_3,by = "shape_id") %>%
  select(-route_id,-route_id.x,-route_id.y,-agency_id.x,-agency_id.y,-rn3) %>% 
  mutate(start = str_detect(origin, "start"),
         end = str_detect(destination, "end"),
         is_stop_in_orig = coalesce(route_desc == route_desc.x,FALSE),
         is_stop_in_dest = coalesce(route_desc == route_desc.y,FALSE),
         is_orig_intersection = str_detect(origin,"inter"),
         is_dest_intersection = str_detect(destination,"inter"),
         orig_cut = ifelse(start | is_stop_in_orig, "first",
                           ifelse(is_orig_intersection,"third","second")),
         dest_cut = ifelse(start | is_stop_in_dest, "first",
                           ifelse(is_dest_intersection,"third","second")))
prod1 %>% 
  filter(shape_id == "124301") %>% 
  select(origin,destination) %>% 
  mutate(orig_cut_lead = lead(origin),
         clear = destination == orig_cut_lead) %>% 
  View()
# first order: first(name), last(name), stop of line
# second order: stop not of line
# third order: intersections

with_shapes %>% 
  filter(!shape_id %in% to_remove) %>% 
  filter(!st_geometry_type(geometry.x) == "LINESTRING") %>% 
  st_line_merge()
slice(10) %>% 
  mapview()
as_tibble() %>% 
  left_join(stops_in_routes,by = c("origin" = "stop_id","shape_id")) %>% 
  left_join(stops_in_routes,by = c("destination" = "stop_id","shape_id")) %>% 
  group_by(rn3,shape_id) 


one <- joined %>% filter(shape_id == bs_lines[1,]$shape_id)
two <- joined_inter%>% filter(shape_id == bs_lines[1,]$shape_id)
net <- bs_lines[1,] %>%
  as_sfnetwork(directed = T) %>%
  convert(to_spatial_subdivision) %>%
  st_network_blend(two)
net %E>% 
  mutate(rn = row_number()) %>% 
  filter(rn %in% c(1,2)) %N>% 
  filter(!node_is_isolated()) %E>%{
    (.) %E>%
      autoplot()  +
      aes(color =(.) %>% as_tibble() %>% pull(rn)) -
      ggspatial::annotation_map_tile(zoom = 17)[[1]]
  }
net %E>% 
  mutate(rn = row_number()) %>% 
  filter(rn %in% c(19,23)) %N>% 
  filter(!node_is_isolated()) %E>% 
  st_as_sf() %>% ggplot(aes(color = factor(rn))) + geom_sf()
rel_rows <- net %E>%
  st_as_sf() %>%
  mutate(rn = row_number()) %>%
  st_buffer(-10,singleSide = T) %>%
  filter(apply(st_intersects(geometry,one,sparse = F),1,any))
rel_idx <- net %E>%
  st_as_sf() %>%
  mutate(rn = row_number()) %>%
  st_buffer(-10,singleSide = T) %>%
  st_intersects(one,sparse = F) %>%
  apply( 1,which)

nets <- net %E>%
  st_as_sf() %>%
  mutate(rn = row_number(),
         idx = rel_idx) %>%
  filter(rn %in% rel_rows$rn) %>%
  mutate(net1 = map2(geometry,rn, ~.x %>% st_sfc(crs = 2039) %>% st_sf() %>% as_sfnetwork(directed = T) %E>% mutate(rn = .y)),
         net2 = map2(net1,idx, function(x,y){
           st_network_blend(x,one[y,]) %>% mutate(rn2 = row_number())
         })) %>%
  pull(net2)

`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}


stops <- bind_rows(map(nets,~.x %N>% st_as_sf())) %>%
  filter(!is.na(stop_id)) %>%
  distinct()
final_net <- net %E>%
  mutate(rn = row_number()) %>%
  filter(!rn %in% rel_rows$rn) %E>%
  bind_graphs(bind_graphs(nets)) %>%
  as_sfnetwork(edges_as_lines = T,directed = T) %E>%
  st_as_sf() %>%
  arrange(rn,rn2) %>%
  as_sfnetwork(directed = T)
route <- final_net %E>%
  st_as_sf() %>%
  mutate(cs = cumsum(st_length(geometry)))
cuts <- final_net %>%
  st_join(stops) %>%
  st_join(two) %>%
  mutate(rn = row_number()) %>%
  as_tibble() %>%
  as_tibble() %>%
  select(geometry,stop_id,inter_id,rn)
bind_graphs(nets[c(1,6)]) %>% 
  as_sfnetwork(edges_as_lines = T) %>% 
  autoplot() - ggspatial::annotation_map_tile(zoom = 17)[[1]] + geom_sf(data = cuts[20:22,] %>% st_sf(),mapping = aes(color = rn))
phase_1 <- route %>%
  left_join(cuts,by = c("from" = "rn")) %>%
  left_join(cuts,by = c("to" = "rn")) %>%
  fill(shape_id,.direction = "downup") %>%
  mutate(origin = ifelse(!is.na(stop_id.x),stop_id.x,ifelse(!is.na(inter_id.x),inter_id.x,NA)),
         destination = ifelse(!is.na(stop_id.y),stop_id.y,ifelse(!is.na(inter_id.y),inter_id.y,NA)),
         origin = ifelse(row_number() == 1,paste0("start_",shape_id),origin),
         destination = ifelse(row_number() == nrow(.),paste0("end_",shape_id),destination)) %>%
  select(origin,destination,from,to) %>%
  fill(origin) %>%
  fill(destination,.direction = "up") %>% 
  mutate(lag1 = lag(origin),
         lag2 = lag(destination),
         counter =coalesce(ifelse(origin == lag1 & destination == lag2,0,1),1) %>% cumsum()) %>% 
  select(-lag1,-lag2)

phase_2 <- phase_1 %>%
  group_by(origin,destination,counter) %>%
  summarise(do_union=T)
phase_21 <- phase_2 %>%
  filter(st_geometry_type(geometry.x) == "MULTILINESTRING")
if(nrow(phase_21) == 0){
  phase_3 <- phase_2
}else{
  phase_3 <- bind_rows(phase_2 %>%
                         filter(!st_geometry_type(geometry.x) == "MULTILINESTRING"),
                       phase_2 %>%
                         filter(st_geometry_type(geometry.x) == "MULTILINESTRING") %>%
                         st_line_merge())
}
phase_3 %>%
  ungroup() %>%
  arrange(counter) %>%
  # st_buffer(-10,singleSide = T) %>%
  mutate(rn3 = row_number()) %>%
  mapview(zcol = "rn3") +
  mapview(st_as_sf(cuts),col.regions = "red")



# idx_list <- phase_1 %>%
#   group_by(origin,destination) %>%
#   group_indices() %>% unique()
# 
# phase_1 %>%
#   group_by(origin,destination) %>%
#   mutate(gg = which(cur_group_id() == idx_list)) %>%
#   ungroup() %>%
#   mutate(mat = map(geometry.x,st_coordinates)) %>%
#   group_by(origin,destination,gg)  %>%
#   mutate(rn  = row_number(),
#          is_max = rn == max(rn),
#          mat = map2(mat,is_max,~if(.y){
#            .x
#          } else{
#            head(.x,-1)
#          })) %>%
#   st_drop_geometry() %>%
#   summarise(mat = list(mat)) %>%
#   ungroup() %>% View()
#   mutate(mat = map(mat,~rbind(.x) %>% bind_rows()))
#   arrange(gg)
# 
# 
# 
# phase_2 <- phase_1 %>%
#   group_by(origin,destination) %>%
#   mutate(gg = which(cur_group_id() == idx_list)) %>%
#   group_by(origin,destination,gg) %>%
#   summarise(do_union=T)
# phase_21 <- phase_2 %>%
#   filter(st_geometry_type(geometry.x) == "MULTILINESTRING")
# if(nrow(phase_21) == 0){
#   phase_3 <- phase_2
# }else{
#   phase_3 <- bind_rows(phase_2 %>%
#                          filter(!st_geometry_type(geometry.x) == "MULTILINESTRING"),
#                        phase_2 %>%
#                          filter(st_geometry_type(geometry.x) == "MULTILINESTRING") %>%
#                          st_line_merge())
# }
# phase_3 %>%
#   ungroup() %>%
#   arrange(gg) %>%
#   # st_buffer(-10,singleSide = T) %>%
#   mutate(rn3 = row_number()) %>%
#   mapview(zcol = "rn3") +
#   mapview(st_as_sf(cuts),col.regions = "red")
# mapview(two,col.regions = "yellow")
# 
# 
# 
# 
# 
# 

