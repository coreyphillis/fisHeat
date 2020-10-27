
# Table 1. Parameters of Ratkowsky model derived at full ration by Perry et al 2015

b <-    0.388 # allometric growth exponent
b_se <- 0.025 # SE of b allometric exponent estimate

d <-    0.415 # shape parameter
d_se <- 0.037 # SE of shape parameter

g <-    0.315 # slope of shape parameter
g_se <- 0.044 # SE of shape parameter

TL <-    1.833 # lower temp threshold where growth goes to zero
TL_se <- 0.631 # SE of lower temp threshold where growth goes to zero

Topt <- 19.019 # Temp that maximizes growth
TU <-   24.918 # Upper Temp threshold
c <-     6.016 # specific growth rate of a 1-gram fish at Temp that maximizes growth

# Table 3 from Manhard et al 2018. DATA IS FOR COHO.
# how do these parameters relate to the above?
alphaT <- 2.483 # intercept of TM
BetaT_coho <-  0.306 # slope of TM from Manhard
BetaT <-  0.4624 # Fudged value from Brett's spreadsheet model
Betag <-  0.315 # chinook value (coho value = 0.142)

# SIT size classes
# small: < 42mm, medium: 42-72, large: 72 - 110, very large: >110
small_fl <- c(28,42) # assuming 28mm for fl at swim up. need reference
medium_fl <- c(42, 72)
large_fl <- c(72, 110)
vlarge_fl <- c(110, Inf)

# R_fxn Manhard
FORAGING_TIME <- 360 # in minutes
DRY_MASS <- 2.595e-5 # Carson's dry weight estimate for Daphnia pulex (see Zooplankton dry weigth data sheet)

# Type II functional response from Haskell et al 2017 using sub-yearling Chinook and Daphnia
max_consumption_rate = 29.858     # C_max = max consumption rate (Beta 0)
prey_density_at_half_max_consumption = 4.271    # P_Chalf = prey density at which consumption reaches half of its maximum (Beta 1)



# # Mass ~ FL parameters
# a <- 44.199
# f <-  0.3034

# Eqn 18: change in FL per change in Time
# dFL <- function(Temp, R){
#     Omega <- Omega_fxn(Temp)
#     a*(M0^b + (Omega*b*t/100))^f/(b-a*M0^f)
# }


# MWTS Detailed Approach and Model Equations

# eq 9 and 10
mass_fxn <- function(M0, omega, t){
  Mt = (M0 ^ b + ((omega * b * t)/100))^(1/b)
  return(Mt)
}

# Eqn 10
Omega_fxn <- function(Temp, R){
  g <- g_fxn(R)
  Topt <- Topt_fxn(R)
  TU <- TU_fxn(g, Topt, R)
  Omega = d*(Temp - TL) * (1 - exp(g * (Temp - TU)))
  return(Omega)
}

# Eqn 11
TU_fxn <- function(g, Topt, R){
  g <- g_fxn(R)
  TU <- Topt + ((log(1 + g*(Topt-TL)))/g)
  return(TU)
}

# Eqn 12
Topt_fxn <- function(R){
  Topt = exp(alphaT + BetaT * R)
  return(Topt)
}

# Eqn 13
g_fxn <- function(R){
  g = Betag * R
  return(g)
}


R_fxn <- function(foraging_time = FORAGING_TIME, # in minutes
                  consumption_rate, # from TypeII
                  dry_mass = DRY_MASS,
                  Temperature,
                  M0,
                  ...){
  # Eqn 3 Manhard
  R_max = (-8.82 + 2.5 * log(1.8 * Temperature + 32))*(0.17 * M0 ^-0.33)
  #C = (C_max * P)/(P_Chalf + P)
  R_dry = consumption_rate * foraging_time * dry_mass
  #R_dry = typeII(prey_density = P) * foraging_time * dry_mass
  #return(R_max)
  #return(R_dry)
  return(R_dry/R_max)
}


#' Type II functional feeding response
#'
#' Type II functional response from Haskell et al 2017 using sub-yearling Chinook and Daphnia
#'
#' @param C_max max consumption rate (Beta 0)
#' @param prey_density Prey density (individuals per liter)
#' @param P_Chalf prey density at which consumption reaches half of its maximum (Beta 1)
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
typeII <- function(C_max = max_consumption_rate,
                   prey_density,
                   P_Chalf = prey_density_at_half_max_consumption,
                   ...){

  consumption_rate = (C_max * prey_density)/(P_Chalf + prey_density)
  return(consumption_rate)

}



#' Calculate new mass
#'
#' Calculate the new mass of a fish from its old mass after a specficed number of days, prey density, and temperature
#'
#' @param df a dataframe with the columns "temperature" and "prey_density"
#' @param old_mass initial mass of fish in grams
#' @param day day on which mass should be calculated
#'
#' @return Vector
#' @export
#'
#' @examples
newmass_fxn <- function(df, old_mass, day) {

  # consumption rate given available prey
  C_prey = typeII(prey_density = df$prey_density)

  # Ration given consumption rate, fish mass, and water temperature
  R_temp = R_fxn(consumption_rate = C_prey, M0 = old_mass, Temperature = df$temperature)

  # Omega given ration and water temperature
  Omega_temp = Omega_fxn(Temp = df$temperature, R = R_temp)

  # new mass after t days given previous mass and Omega
  new_mass = mass_fxn(M0 = old_mass, omega = Omega_temp, t = day)
  return(new_mass)
}


# decay function ----------------------------------------------------------

# decay_fxn <- function(t, # time, but in this case distance
#                       N0 = 60, # initial quantity (prey density)
#                       lambda = 1 # decay constant
#                       ){
#     N_t = N0 * exp(-1*lambda *t)
#     return(N_t)
# }



# L:W regression functions ----------------------------------------------------------

#' Convert fork length to mass
#'
#' @param fl fork length to convert
#' @param a coefficient a
#' @param b coefficient b
#'
#' @return
#' @export
#'
#' @examples
L2W_fxn <- function(fl,
                    a = 0.005,
                    b = 3.44) {
  # Chinook biomass coefficients from Pascale EDI Data
  # Original study is Kimmerer et al., 2015
  # coefficients are for cm to gram

  fl_cm = fl/10 # converting FL to cm

  mass = a * (fl_cm) ^ b

  return(mass)
}

#' Convert mass to fork length
#'
#' @param mass mass to convert
#' @param a coefficient
#' @param b coefficient
#'
#' @return
#' @export
#'
#' @examples
W2L_fxn <- function(mass,
                    a = 0.005,
                    b = 3.44) {
  # Chinook biomass coefficients from Pascale EDI Data
  # Original study is Kimmerer et al., 2015
  # coefficients are for cm to gram

  fl_cm = (mass/a)^(1/b)

  return(fl_cm *10)

}





#' Create fish population
#'
#' Function to create fish population with specified abundance and
#' fork length distribution.
#'
#' @param nfish number of fish in the population
#' @param init_fl Named list ('min', 'max', 'mode') of initial fork lengths when
#'   distribution of fish for fl_distribution = 'triangle' (default)
#' @param fl_distribution 'triangle' uses rtriangle() distribution.
#'   'discrete' creates three discrete size classes ('small', 'medium', 'large')
#'   according to psizeclass
#' @param psizeclass proportion in each size class when fl_distribution = 'discrete'
#' @param ...
#'
#' @return list
#' @export
#'
#' @examples

create_fish_pop <- function(nfish = 5000,
                            init_fl = list(
                              'min' =  28,
                              'max' =  80,
                              'mode'  = 38),
                            fl_distribution = 'triangle',
                            psizeclass = c(0.3, 0.5, 0.2), # proportion in each size class if fl_distribution = 'discrete'
                            ...) {

  # SIT size classes
  # small: < 42mm, medium: 42-72, large: 72 - 110, very large: >110
  # small_fl <- c(28,42) # assuming 28mm for fl at swim up. need reference
  # medium_fl <- c(42, 72)
  # large_fl <- c(72, 110)
  # vlarge_fl <- c(110, Inf)

  small_rng <- L2W_fxn(small_fl)
  medium_rng <- L2W_fxn(medium_fl)
  large_rng <- L2W_fxn(large_fl)
  vlarge_rng <- L2W_fxn(vlarge_fl)

  if(fl_distribution == 'triangle'){
    # fill data frame with initial fork lengths from a triangle distribution and convert FL to mass in same step
    size_dist <- data.frame(
      init_size = L2W_fxn(triangle::rtriangle(n = nfish, a = init_fl$min, b = init_fl$max, c = init_fl$mode)),
      size_class = NA)

    size_dist$size_class <- ifelse(dplyr::between(size_dist$init_size, small_rng[1], small_rng[2]), 'small',
                                   ifelse(dplyr::between(size_dist$init_size, medium_rng[1], medium_rng[2]), 'medium',
                                          ifelse(dplyr::between(size_dist$init_size, large_rng[1], large_rng[2]), 'large',
                                                 ifelse(dplyr::between(size_dist$init_size, vlarge_rng[1], vlarge_rng[2]), 'very large', NA))))

    size_dist$size_class <- factor(size_dist$size_class, levels = c('small', 'medium', 'large', 'very large'))


    notverylarge <- sum(size_dist$size_class != 'very large', na.rm = TRUE)


    psmall <- sum(size_dist$size_class == 'small')/notverylarge
    pmedium <- sum(size_dist$size_class == 'medium')/notverylarge
    plarge <- sum(size_dist$size_class == 'large')/notverylarge

    psize <- c(psmall, pmedium, plarge)

    fish_init <- data.frame(size_class = c('small', 'medium', 'large'),
                            init_size = c(mean(size_dist$init_size[size_dist$size_class == 'small'], na.rm = TRUE),
                                          mean(size_dist$init_size[size_dist$size_class == 'medium'], na.rm = TRUE),
                                          mean(size_dist$init_size[size_dist$size_class == 'large'], na.rm = TRUE)),
                            psize = psize)

    fish <- dplyr::filter(size_dist, size_class != 'very large')

  }

  if(fl_distribution == 'discrete'){
    stopifnot(sum(psizeclass) == 1)

    # init_fl for recovering initial fl in discrete distribution mode
    init_fl <- list('small' = mean(small_fl),
                    'medium' = mean(medium_fl),
                    'large' = mean(large_fl),
                    'vlarge' = mean(vlarge_fl[1]))

    # initial masses of small, medium, and large size class fish for discrete sizes

    init_mass_s <- L2W_fxn(mean(small_fl))
    init_mass_m <- L2W_fxn(mean(medium_fl))
    init_mass_l <- L2W_fxn(mean(large_fl))
    init_mass_vl <- L2W_fxn(vlarge_fl[1])

    # data frame for discrete initial sizes
    size_dist <- data.frame('size_class' = factor(c('small', 'medium', 'large', 'very large'),
                                                  levels = c('small', 'medium', 'large', 'very large')),
                            'init_size' = c(init_mass_s, init_mass_m, init_mass_l, init_mass_vl))


    # discrete and triangle distributions w/o very large size class, w/ number of rows equal to nfish to assign
    fish <- dplyr::filter(size_dist, size_class != 'very large') %>%
      sample_n(size = nfish, replace = TRUE, weight = psizeclass)

    fish_init <- data.frame(size_class = c('small', 'medium', 'large'),
                            init_size = c(mean(fish$init_size[fish$size_class == 'small'], na.rm = TRUE),
                                          mean(fish$init_size[fish$size_class == 'medium'], na.rm = TRUE),
                                          mean(fish$init_size[fish$size_class == 'large'], na.rm = TRUE)),
                            psize = psizeclass)

  }

  fish_out <-list(
    'nfish' = nfish,
    'init_fl' = init_fl,
    'fl_distribution' = fl_distribution,
    'fish_init' = fish_init,
    'fish' = fish
  )
  return(fish_out)

}



# function for filling a reach with fish ----------------------------------

fill_reach <- function(reach_df,
                       fish_list,
                       territory = list(
                         'small' = 1,
                         'medium' = 3,
                         'large' = 9
                       ),
                       terr_max = 9,
                       spatial_distribution = 'random',
                       size_hierarchy = FALSE
) {
  nfish <- fish_list$nfish
  psize <- fish_list$fish_init$psize
  fish_init <- fish_list$fish_init

  fish <- fish_list$fish %>%
    mutate(territory = case_when(size_class == 'small' ~ territory$small,
                                 size_class == 'medium' ~ territory$medium,
                                 size_class == 'large' ~ territory$large))

  fish_id <- paste0('fish',1:terr_max)

  # add additional columns with name fish_id for fish up to territory max
  tmp_df <- cbind(reach_df, as.data.frame(
    array(NA,
          dim = c(nrow(reach_df), terr_max),
          dimnames = list(NULL, fish_id))))

  # an array to keep track of size class, initial size, and territory area
  fish_array_all <- array(data.matrix(tmp_df),
                          dim = c(nrow(tmp_df), ncol(tmp_df), 3),
                          dimnames = list(NULL, colnames(tmp_df), colnames(fish)))

  # subset initial day to fill reach initial conditions and
  # exclude ag input (reach_width < 0) from reach accessible to fish
  fish_array <- fish_array_all[which(fish_array_all[,'day',1] == 0 & fish_array_all[,'reach_width',1] > 0),,]



  if(size_hierarchy == TRUE) fish <- arrange((fish), desc(init_size))

  # initialize while loop
  n_small <- n_medium <- n_large <- n_vlarge <- 0
  fish_counter <- 0

  # the while loop assigns fish in APPROXIMATE proportion to psize. Large fish will be undersampled in scenarios where
  # size hierarchy rules are not in place.


  while(fish_counter < nfish){ # continue filling reach_cells while fish are available
    fish_counter = fish_counter + 1
    #print(fish_counter) #  print for debugging where loops are getting hung up

    # sample a fish to assign to a reach according to proportion size classes in the population

    fish_i <- fish[fish_counter, ]


    #print(fish_i$size_class)


    # skipping to next iteration of while loop if all inidividuals in a size class have been assigned
    # if(fish_i$size_class == 'small' && n_small >= nfish*psize[1]) next
    # if(fish_i$size_class == 'medium' && n_medium >= nfish*psize[2]) next
    # if(fish_i$size_class == 'large' && n_large >= nfish*psize[3]) next

    # draw a reach_cell that has available territory at random
    # array of only cells w/ territory <= fish_i territory
    available_cells <- fish_array[which(terr_max - rowSums(fish_array[, fish_id, 'territory'], na.rm = TRUE) >= fish_i$territory),'reach_cell','territory']
    if(length(available_cells) == 0) next

    if(spatial_distribution == 'random'){
      cell_i = available_cells[sample(length(available_cells), size = 1)]
    }

    if(spatial_distribution == 'IFD'){
      # find cells with max prey density
      maxprey <- max(fish_array[which(fish_array[,'reach_cell',1] %in% available_cells), 'prey_density', 1])
      # array of cells with max prey density
      available_cells_maxprey <- fish_array[which(fish_array[, 'reach_cell',1] %in% available_cells & fish_array[,'prey_density',1] == maxprey),'reach_cell','territory']
      # sample randomly one of the available max prey density cells
      # ?sample
      # when x is length 1, sample defaults to 1:x,
      # below is a workaround to avoid this behavior
      cell_i = available_cells_maxprey[sample(length(available_cells_maxprey), size = 1)]
    }

    # cell density prior to fish placement
    # cell_density = sum(fish_array[cell_i, fish_id, 'territory'], na.rm = TRUE) # indexing columns fish1:fish9
    # if(cell_density + fish_i$territory > terr_max) next # if territory is full, go to next iteration in while loop
    #

    for(fish_j in fish_id){
      # break out of for loop if no territory for fish_i
      if(sum(fish_array[which(fish_array[,'reach_cell','territory'] == cell_i), fish_id, "territory"], na.rm = TRUE) >= terr_max) break
      #print(fish_j)
      if(is.na(fish_array[which(fish_array[,'reach_cell' , 'size_class'] == cell_i),fish_j,1])){

        fish_array[which(fish_array[,'reach_cell' , 'size_class'] == cell_i),fish_j, 'size_class'] = fish_i$size_class
        fish_array[which(fish_array[,'reach_cell' , 'size_class'] == cell_i),fish_j, 'init_size'] = fish_i$init_size
        fish_array[which(fish_array[,'reach_cell' , 'size_class'] == cell_i),fish_j, 'territory'] = fish_i$territory
        # fish_i has been assigned to fish_j, break out of for loop
        break
      }

      next
    }

    # keeping track of number of fish in each size class that have been assigned
    n_small = sum(fish_array[, fish_id, 'size_class'] == 1, na.rm = TRUE)
    n_medium = sum(fish_array[, fish_id, 'size_class'] == 2, na.rm = TRUE)
    n_large = sum(fish_array[, fish_id, 'size_class'] == 3, na.rm = TRUE)
    n_vlarge = sum(fish_array[, fish_id, 'size_class'] == 4, na.rm = TRUE)

    # keeping track of total number of fish assigned
    fish_assigned = sum(!is.na(fish_array[, fish_id, 1]))

    fish_init$init_n <- c(n_small, n_medium, n_large) #, n_vlarge)
    cat(fish_assigned, ' of ', nfish,' available fish assigned;   ','    small fish:', fish_init$init_n[1], '  medium fish:', fish_init$init_n[2], '  large fish:', fish_init$init_n[3], '\r')
    flush.console()

    # control rules to break out of while loop if certain conditions are met:

    # break out of while loop if all cells are at max density
    if(min(rowSums(fish_array[, fish_id, 'territory'], na.rm = TRUE)) == terr_max) break

    # break if all of the small fish have been assigned and only small territory is available
    if(n_small == nfish*psize[1] &&
       min(rowSums(fish_array[, fish_id, 'territory'], na.rm = TRUE)) < territory$medium) break

    # break if all of the small and medium fish have been assigned and only small or medium territory is available
    if(n_small == nfish*psize[1] && n_medium == nfish*psize[2] &&
       min(rowSums(fish_array[, fish_id, 'territory'], na.rm = TRUE)) < territory$large) break

    # if the fish did not get assigned, return to the top of the while loop and try again
    #fish_counter <- if_else(fish_assigned == fish_counter, fish_counter + 1, fish_counter)

  }

  paste0(fish_assigned, ' fish were assigned from the ', fish_counter, ' selected and ', nfish, ' total available')
  fish_init$init_psize <- round(fish_init$init_n/fish_assigned, 2)


  fish_array_all[which(fish_array_all[,'day',1] == 0 & fish_array_all[,'reach_width',1] > 0),,] <- fish_array

  reach_out <- list('fish_available' = nfish,
                    'fish_assigned' = fish_assigned,
                    'fish_init' = fish_init,
                    'fish_array' = fish_array_all)

  return(reach_out)

}


# baseline reach function -------------------------------------------------

# create a reach w/o subsidy from a filled reach
basereach <- function(filled_reach, bkgrd_prey = river_zoop, bkgrd_temp = river_temp) {
  filled_reach$fish_array[,'prey_density',] = bkgrd_prey
  filled_reach$fish_array[,'temperature',] = bkgrd_temp
  return(filled_reach)
}

# function for growing fish in a reach ------------------------------------


# Really Sloooow on my computer. >5 minutes for 900 cell reach w/ max 9 fish
# How to speed up??? grouped reach_cell calculations should be able to run in parallel

grow_sum_fish <- function(filled_reach_list, ...) {

  full_df <- as.data.frame(filled_reach_list$fish_array[,, 'init_size'])
  occupied_cells <- full_df %>%
    dplyr::filter(day == 0, !is.na(fish1)) %>%
    pull(reach_cell)

  df <- dplyr::filter(full_df, reach_cell %in% occupied_cells)

  jdays <- unique(df$day)

  # get column names of fish_id
  fish_id <- dplyr::select(df, starts_with('fish')) %>% colnames()

  for(i in 1:max(df$reach_cell)){
    for(j in jdays[-1]){  # over j days, skipping day 0 (i.e. initial mass)
      df[which(df$reach_cell == i & df$day == j), fish_id] <- newmass_fxn(
        df[which(df$reach_cell == i & df$day == j),],
        df[which(df$reach_cell == i & df$day == jdays[match(j, jdays)-1]), fish_id],
        day = j - jdays[match(j, jdays)-1])
    }
  }

  return(df)

}



# Make reach from delt model output ---------------------------------------


make_reach <- function(reach_in,
                       ag_zoop,
                       river_zoop,
                       ag_temp,
                       river_temp,
                       vdays) # a vector of days
{

  reach_tmp <- data.frame(reach_length = reach_in$x.coordinate,
                          reach_width = reach_in$y.coordinate,
                          reach_depth = reach_in$water.depth..m.,
                          reach_cell = 1:nrow(reach_in),
                          prey_density = (reach_in$zoop * ag_zoop + (1-reach_in$zoop) * river_zoop),
                          temperature = (reach_in$zoop * ag_temp + (1- reach_in$zoop) * river_temp),
                          p = round(reach_in$zoop,2))

  reach_df <- do.call('rbind', replicate(length(vdays)+1, reach_tmp, simplify = FALSE))

  # day 0 for initial mass, then subsequent vdays where mass will be calculated
  reach_df$day <- rep(c(0, vdays), each = max(reach_df$reach_cell))

  return(reach_df)
}


# calculate number of cells that are subsidized ---------------------------

reach_subsidized <- function(filled_reach_list) {
  full_df <- as.data.frame(filled_reach_list$fish_array[,, 'territory']) %>% dplyr::filter(day == 0) %>% dplyr::filter(reach_width > 0)

  n_cells <- max(full_df$reach_cell)
  bkgrd <- min(full_df$prey_density)

  # prey density in cells w/ at least 1 fish
  one_fish <- full_df[which(!is.na(full_df$fish1)),]$prey_density

  # prey density in cells w/ at least 2 fish
  multi_fish <- full_df[which(!is.na(full_df$fish2)), ]$prey_density

  # prey density in cells w/ 1 large fish that occupies all of the territory
  large_fish <- full_df[which(na.exclude(full_df$fish1 == 9)), ]$prey_density



  df_out <- data.frame("x_bkgrd" =
                         c(1.1, 1.5, 2, 10, 100, 1000),
                       'p_reach' = NA,
                       'p_1_fish' = NA,
                       'p_multi_fish' = NA)

  for(i in 1:nrow(df_out)){
    df_out$p_reach[i] = round(sum(full_df$prey_density > (bkgrd * df_out$x_bkgrd[i]))/n_cells, 3) *100
    df_out$p_1_fish[i] = round(sum(one_fish > (bkgrd * df_out$x_bkgrd[i]))/length(one_fish), 3) *100
    df_out$p_multi_fish[i] = round((sum(multi_fish > (bkgrd * df_out$x_bkgrd[i])) + sum(large_fish > (bkgrd * df_out$x_bkgrd[i])))/(length(multi_fish) + length(large_fish)), 3) * 100

    sum(full_df$prey_density > (bkgrd * df_out$x_bkgrd[i]))/n_cells
  }
  return(df_out)
}
