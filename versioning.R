pacman::p_load(tidyverse,sf,lwgeom,sfnetworks,tidygraph,gtfstools,mapview,job,furrr,nngeo,assertr)
# library(tidyverse)
# library(sf)
# library(lwgeom)
# library(nngeo)
# library(sfnetworks)
# library(tidygraph)
# library(gtfstools)
# library(mapview)
# library(job)
early_gtfs <- read_gtfs("data/gtfs_010121.zip",encoding = "UTF-8") %>%  filter_by_agency_id("32")
late_gtfs <- read_gtfs("data/gtfs_010122.zip",encoding = "UTF-8")%>%  filter_by_route_id(early_gtfs$routes$route_id)
all_process <- function(gtfsq){
  linesq <-function(gtfs,crs){
    get_trip_geometry(gtfs,file = "shapes") %>% 
      # join trips to their geometry
      left_join(gtfs$trips, by = "trip_id") %>% 
      # get only one geometry per shape_id
      select(shape_id) %>%
      group_by(shape_id) %>% 
      mutate(rn = row_number()) %>% 
      filter(rn == 1) %>%
      ungroup() %>% 
      select(-rn) %>% 
      st_sf() %>% 
      st_transform(crs)
  } 
  stops <- function(gtfs,crs){
    gtfs$stops %>% 
      st_as_sf(coords = c("stop_lon","stop_lat"),crs = 4326) %>% 
      st_transform(crs) %>% 
      select(stop_id)
  }
  intersections <- function(lines){
    # explode all lines
    st_segments(lines$geometry) %>% 
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
  }
  # for each stop, join a relevant line
  joined_stops <- function(stops,lines){
    stops %>% 
      st_join(lines %>% 
                st_buffer(-10,singleSide = T) ,left = FALSE) %>% 
      arrange(stop_id) 
  }
  # for each intesection, join a relevant line
  joined_inter <- function(intersections,lines){
    intersections %>% 
      st_join(lines %>% 
                st_buffer(-10,singleSide = T) ,left = FALSE) 
  }
  filter_shape_id <- function(obj,rel_shape_id){
    obj %>% filter(shape_id == rel_shape_id)
  }
  
  
  segmentize_line <- function(line,stops,intersections){
    print(line$shape_id)
    # blend intersections to sfnetwork
    net <- line %>% 
      as_sfnetwork(directed = T) %>% 
      convert(to_spatial_subdivision) %>%
      st_network_blend(intersections)
    # filter segments where there are stops
    rel_rows <- net %E>%
      st_as_sf() %>% 
      mutate(rn = row_number()) %>% 
      st_buffer(-10,singleSide = T) %>% 
      filter(apply(st_intersects(geometry,stops,sparse = F),1,any)) 
    # get stops ids in order to blend into network
    rel_idx <- net %E>%
      st_as_sf() %>% 
      mutate(rn = row_number()) %>% 
      st_buffer(-10,singleSide = T) %>% 
      st_intersects(stops,sparse = F) %>% 
      apply( 1,which)
    nets <- net %E>%
      st_as_sf() %>% 
      mutate(rn = row_number(),
             idx = rel_idx) %>% 
      filter(rn %in% rel_rows$rn) %>% 
      mutate(
        # create an sfnetwork for each segment according to row number
        net1 = map2(geometry,rn, ~.x %>% 
                      st_sfc(crs = 2039) %>% 
                      st_sf() %>% 
                      as_sfnetwork(directed = T) %E>% 
                      mutate(rn = .y)),
        # blend stops into each sfnetwork
        net2 = map2(net1,idx, function(x,y){
          st_network_blend(x,stops[y,]) %>% 
            mutate(rn2 = row_number())
        })) %>% 
      pull(net2)
    # get all joined net stops
    stops <- bind_rows(map(nets,~.x %N>% st_as_sf())) %>% 
      filter(!is.na(stop_id)) %>% 
      distinct()
    # bind all graphs together, including those segments without stops
    final_net <- net %E>%
      mutate(rn = row_number()) %>%
      filter(!rn %in% rel_rows$rn) %E>%
      bind_graphs(bind_graphs(nets)) %>%
      # convert to sfnetwork
      as_sfnetwork(edges_as_lines = T,directed = T) %E>% 
      st_as_sf() %>%
      # arrange by order
      arrange(rn,rn2) %>% 
      as_sfnetwork(directed = T)
    route <- final_net %E>% 
      st_as_sf() %>% 
      mutate(cs = cumsum(st_length(geometry))) 
    # create cuts
    cuts <- final_net %>% 
      st_join(stops) %>%
      st_join(intersections) %>% 
      mutate(rn = row_number()) %>% 
      as_tibble() %>% 
      as_tibble() %>%
      select(geometry,stop_id,inter_id,rn)
    # process route and cuts
    phase_1 <- route %>%
      # join all cuts
      left_join(cuts,by = c("from" = "rn")) %>%
      left_join(cuts,by = c("to" = "rn")) %>%
      # fill all shape_ids to prevent errors
      fill(shape_id,.direction = "downup") %>%
      # give each origin and destination an id: stop id if stop, intersection id if intersection, start\end if edges
      mutate(origin = ifelse(!is.na(stop_id.x),stop_id.x,ifelse(!is.na(inter_id.x),inter_id.x,NA)),
             destination = ifelse(!is.na(stop_id.y),stop_id.y,ifelse(!is.na(inter_id.y),inter_id.y,NA)),
             origin = ifelse(row_number() == 1,paste0("start_",shape_id),origin),
             destination = ifelse(row_number() == nrow(.),paste0("end_",shape_id),destination)) %>%
      select(origin,destination,from,to) %>%
      fill(origin) %>%
      fill(destination,.direction = "up") %>%
      # preserve order
      mutate(lag1 = lag(origin),
             lag2 = lag(destination),
             counter =coalesce(ifelse(origin == lag1 & destination == lag2,0,1),1) %>% cumsum()) %>% 
      select(-lag1,-lag2)
    # join created lines, some have no named edge
    phase_2 <- phase_1 %>%
      group_by(origin,destination,counter) %>%
      summarise(do_union=T)
    # prevent bugs if no multilinestring
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
  }
  gtfsq_lines <- gtfsq %>% linesq(2039)
  gtfsq_stops <- gtfsq %>% stops(2039)
  gtfsq_intersections <- gtfsq_lines %>% intersections()
  gtfsq_joined_stops <- joined_stops(gtfsq_stops,gtfsq_lines)
  gtfsq_joined_intersections <- joined_inter(gtfsq_intersections,gtfsq_lines)
  print("segmentizing")
  res <- imap(gtfsq_lines$shape_id,function(chosen_shape_id,idx){
    out <- tryCatch(
      {
        line <- filter_shape_id(gtfsq_lines,chosen_shape_id)
        rel_stops <- filter_shape_id(gtfsq_joined_stops,chosen_shape_id)
        rel_intersections <- filter_shape_id(gtfsq_joined_intersections,chosen_shape_id)
        segmentize_line(line,rel_stops,rel_intersections)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)},
      finally = print(paste("line number",idx))
    )
    return(out)
  })
  print("remove na")
  names(res) <- gtfsq_lines$shape_id
  resna <- res[is.na(res)]
  resnotna <- res[!is.na(res)]
  with_shapes <- map(names(resnotna),~resnotna[[.x]] %>% mutate(shape_id = .x)) %>% bind_rows()
  starts <- with_shapes %>% st_startpoint() %>% st_sf()
  routes_to_remove <- names(resna)
  stops_in_routes <- gtfsq$routes %>% 
    left_join(gtfsq$trips,by = "route_id") %>% 
    left_join(gtfsq$stop_times,by = "trip_id") %>% 
    select(route_id,route_desc,shape_id,stop_id,agency_id) %>% 
    filter(!shape_id %in% routes_to_remove) %>% 
    distinct()  
  # remove israel rail - useless
  to_remove <- stops_in_routes %>% filter(agency_id== 2) %>% distinct(shape_id) %>% pull(shape_id)
  with_shapes <- with_shapes %>% filter(!shape_id %in% to_remove)
  
  with_shapes <- with_shapes %>% 
    # for each shape_id, add cumulative sum of length to know distance in line
    mutate(length = st_length(geometry.x)) %>% 
    as_tibble() %>% 
    group_by(shape_id) %>% 
    mutate(length= as.numeric(length),
           cs = replace_na(lag(cumsum(length)),0)) %>% 
    left_join(stops_in_routes,by = c("origin" = "stop_id","shape_id")) %>% 
    left_join(stops_in_routes,by = c("destination" = "stop_id","shape_id"))
  # getting relevent makatim without israel rail
  
  makats <- gtfsq$routes %>% 
    filter(!agency_id== 2,route_id %in% unique(stops_in_routes$route_id) ) %>% 
    left_join(gtfsq$trips,by = "route_id") %>% 
    select(route_id,route_desc,shape_id) %>% 
    distinct()
  
  prod1 <- makats %>% 
    left_join(with_shapes,by = "shape_id") %>%
    select(-route_id,-route_id.x,-route_id.y,-agency_id.x,-agency_id.y,-rn3) %>% 
    mutate(start = str_detect(origin, "start"),
           end = str_detect(destination, "end"),
           is_stop_in_orig = coalesce(route_desc == route_desc.x,FALSE),
           is_stop_in_dest = coalesce(route_desc == route_desc.y,FALSE),
           is_orig_intersection = str_detect(origin,"inter"),
           is_dest_intersection = str_detect(destination,"inter"),
           orig_cut = ifelse(start | is_stop_in_orig, "first",
                             ifelse(is_orig_intersection,"third","second")),
           dest_cut = ifelse(end | is_stop_in_dest, "first",
                             ifelse(is_dest_intersection,"third","second")))
  prod2 <- prod1 %>% 
    select(-c(route_desc.x:is_dest_intersection)) %>% 
    mutate(segment_type = ifelse(orig_cut == "third" | dest_cut == "third", "third",
                                 ifelse(orig_cut =="second"|dest_cut =="second","second","first")),
           counter_first = orig_cut == "first",
           counter_second = orig_cut == "second",
           counter_third = orig_cut == "third"
    ) %>% 
    group_by(shape_id) %>% 
    mutate(counter_seg_first = cumsum(counter_first)) %>% 
    ungroup() %>% 
    group_by(shape_id,counter_seg_first) %>% 
    mutate(counter_seg_second = cumsum(counter_second),
           has_second = sum(counter_second) > 0) %>% 
    ungroup() %>% 
    group_by(shape_id,counter_seg_first,counter_seg_second) %>% 
    mutate(counter_seg_third = cumsum(counter_third),
           has_third = sum(counter_third)>0) %>% 
    mutate(geometry.x = st_transform(geometry.x, 4326))
  
  final <- list()
  final$routes <- makats %>% 
    mutate(route_id=as.integer(route_id),
           shape_id = shape_id) %>% 
    select(route_id,shape_id,route_desc)
  st_x = function(x) round(st_coordinates(x)[,1],6)
  st_y = function(x) round(st_coordinates(x)[,2],6)
  final$cuts <- bind_rows(
    prod2 %>% ungroup() %>% 
      select(origin,geometry.x) %>% 
      mutate(geom = st_startpoint(geometry.x)) %>%
      select(-geometry.x) %>% 
      st_sf() %>% 
      rename(cut_id = origin),
    prod2 %>% 
      ungroup() %>% 
      select(destination,geometry.x) %>% 
      mutate(geom = st_endpoint(geometry.x)) %>%
      select(-geometry.x) %>% 
      st_sf() %>% 
      rename(cut_id = destination)
  ) %>% distinct() %>% 
    mutate(intersection = str_detect(cut_id,"inter"),
           start = str_detect(cut_id,"start"),
           end = str_detect(cut_id,"end"),
           type = ifelse(intersection,"inter",ifelse(start,"start",ifelse(end,"end","stop"))),
           lon = st_x(geom),
           lat = st_y(geom)) %>% 
    st_drop_geometry() %>% 
    select(cut_id,type,lon,lat)
  
  final$segments <- prod2 %>% 
    ungroup() %>% 
    select(geometry = geometry.x,length) %>% 
    mutate(geometry = map_chr(geometry,st_as_text,digits = 8)) %>% 
    distinct(geometry,.keep_all = T) %>% 
    mutate(seg_UID = row_number()) %>% 
    select(seg_UID,length,geometry)
  
  
  tmp <- prod2 %>%
    mutate(geometry.x = map_chr(geometry.x,st_as_text,digits = 8)) %>% 
    left_join(final$segments, by = c("geometry.x"="geometry"))
  rm(prod2)
  final$first_order_cuts <- tmp %>% 
    ungroup() %>% 
    select(shape_id,
           first_order_seg_id = counter_seg_first,
           from_cut = origin,
           to_cut = destination) %>% 
    group_by(shape_id,first_order_seg_id) %>% 
    mutate(from_cut = first(from_cut),
           to_cut = last(to_cut)
    ) %>% 
    distinct()
  final$second_order_cuts <- tmp %>% 
    ungroup() %>% 
    filter(has_second) %>% 
    select(shape_id,
           first_order_seg_id = counter_seg_first,
           second_order_seg_id = counter_seg_second,
           from_cut = origin,
           to_cut = destination) %>% 
    distinct() %>% 
    group_by(shape_id,first_order_seg_id,second_order_seg_id) %>% 
    mutate(from_cut = first(from_cut),
           to_cut = last(to_cut)
    ) %>% 
    distinct()  
  final$third_order_cuts <- tmp %>% 
    ungroup() %>% 
    filter(has_third) %>% 
    select(shape_id,
           first_order_seg_id = counter_seg_first,
           second_order_seg_id = counter_seg_second,
           third_order_seg_id = counter_seg_third,
           from_cut = origin,
           to_cut = destination,
           has_second) %>% 
    mutate(second_order_seg_id = ifelse(has_second,second_order_seg_id,NA_integer_)) %>% 
    select(-has_second) 
  final$route_order <- tmp %>% 
    ungroup() %>% 
    select(shape_id, 
           serial_cut_number = counter, 
           first_order_seg_id= counter_seg_first,
           second_order_seg_id= counter_seg_second,
           third_order_seg_id= counter_seg_third,
           seg_UID,
           distance_from_origin = cs,
           serial_cut_number =counter,
           has_second,
           has_third) %>% 
    mutate(second_order_seg_id = ifelse(has_second,second_order_seg_id,NA_integer_),
           third_order_seg_id = ifelse(has_third,third_order_seg_id,NA_integer_)) %>% 
    select(-has_second,-has_third)
  return(final)
}

job({early_gtfs_processed = all_process(early_gtfs)},import = "all")
job({late_gtfs_processed = all_process(late_gtfs)},import = "all")
# check for NA
early_gtfs_processed$segments %>% 
  assert(function(x)!is.na(x),everything())
early_gtfs_processed$cuts %>% 
  assert(function(x)!is.na(x),everything())
early_gtfs_processed$routes %>% 
  assert(function(x)!is.na(x),everything())
early_gtfs_processed$route_order %>% 
  assert(function(x)!is.na(x),shape_id,serial_cut_number,seg_UID,distance_from_origin)
early_gtfs_processed$first_order_cuts %>% 
  assert(function(x)!is.na(x),everything())
early_gtfs_processed$second_order_cuts %>% 
  assert(function(x)!is.na(x),everything())
early_gtfs_processed$third_order_cuts %>% 
  assert(function(x)!is.na(x),-second_order_seg_id)

late_gtfs_processed$segments %>% 
  assert(function(x)!is.na(x),everything())
late_gtfs_processed$cuts %>% 
  assert(function(x)!is.na(x),everything())
late_gtfs_processed$routes %>% 
  assert(function(x)!is.na(x),everything())
late_gtfs_processed$route_order %>% 
  assert(function(x)!is.na(x),shape_id,serial_cut_number,seg_UID,distance_from_origin)
late_gtfs_processed$first_order_cuts %>% 
  assert(function(x)!is.na(x),everything())
late_gtfs_processed$second_order_cuts %>% 
  assert(function(x)!is.na(x),everything())
late_gtfs_processed$third_order_cuts %>% 
  assert(function(x)!is.na(x),-second_order_seg_id)
# length sum is a bit different, in the milimetric level
early_gtfs_processed$route_order %>% 
  left_join(early_gtfs_processed$segments, by = "seg_UID") %>% 
  group_by(shape_id) %>% 
  mutate(cs = cumsum(lag(length) %>% coalesce(0))) %>% 
  filter(abs(distance_from_origin - cs)>0.01) %>% 
  st_as_sf(wkt = "geometry",crs=4326) %>% 
  mapview()

join_between_segments <- early_gtfs_processed$segments %>% 
  full_join(late_gtfs_processed$segments,by = "geometry") 
join_between_segments %>% 
  filter(is.na(seg_UID.x)|is.na(seg_UID.y)) 
  mutate(q=is.na(seg_UID.x)) %>% 
  st_as_sf(wkt = "geometry",crs=4326) %>% 
  mapview(zcol = "q")


early_gtfs_processed$second_order_cuts %>% 
  ungroup() %>% 
  add_count(from_cut,to_cut) %>% 
  filter(n >1)
early_gtfs_processed$route_order

