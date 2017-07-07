#' Lakshadweep Archipelago coral fish community reassembly data (expanded)
#'
#' A list with three data frames containing presence-absence data for the
#' reassembly proccess of coral fish communities in three atolls (Agatti, Kadmat
#' and Kavaratti) of the Lakshadweep Archipelago in India. These data contains a number of 
#' replicates (transects) per sampling time. It is in this respect that expands
#' \code{alonso} community data (see \code{island} R package).
#'
#' Surveys were conducted from 2000 to 2013 in order to follow community
#' reassembly after a coral mass mortality event in the relatively unfished
#' Lakshadweep Archipelago. For most years, transects were taken in four 
#' locations per atoll. Although there might be some underlying heterogeneity, 
#' these transects are approximately taken as true replicates.   
#'
#' @format A list with three dataframes from the 3 different 
#'   atoll. The dataframe has in columns: \describe{ \item{Species}{Name
#'   of the species found} \item{Atoll}{Atoll surveyed} \item{Guild}{Feeding  
#'   strategy of the surveyed species} \item{Presence-absence data}{Several columns 
#'   corresponding to the year in which the surveys were done. Year repetition means
#'   repeated sampling of the same atoll at the same time. Presences are represented
#'   by 1 and true absencestrue or undetected presences by 0.} }
#'
#' @note Detectability per transect results to be of about 0.5, which means that the 
#' parameter 'Detectability' per atoll goes up to almost 0.94 if four transects per 
#' sampling time are taken (1-0.5^4). 
#'
#' @source Alonso, D., Pinyol-Gallemi, A., Alcoverro T. and Arthur, R.. (2015)
#'   Fish community reassembly after a coral mass mortality: higher trophic
#'   groups are subject to increased rates of extinction. \emph{Ecology
#'   Letters}, \bold{18}, 451--461.
#'
#' @name lakshadweep
NULL
