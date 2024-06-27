##' Process netCDF4 Gridded Data into Format of the fingerprint() Function
##' 
##' This function detects the signal factors on the observed data via total 
##' least square linear regression model.
##'
##' @param datafile path to the netCDF4 gridded datafile to be processed
##' @param variable the climate variable to be extracted
##' @param region the longitude and latitude boundary for selected region, 
##'   should match the format of IPCC AR6 regions, the lon and lat of the vertices
##' @param target.year vector of length 2, the starting and ending year of the 
##'   selected time period for D&A analysis
##' @param average number of years for average on each gridbox, default is 5-year average
##' @param reference vector of length 2, the starting and ending year of reference 
##'   time period for computing anomalies
##' @param regridding whether the grid box should be regridded. Specify the size of 
##'   the grid box, e.g., c(40, 30) for \eqn{40^\circ \times 30^\circ} grid box.
##'   If no regridding, leave empty
##' @return a dataset of the processed gridded climate variables for Y, 
##'   Xtilde or control runs
##' @author Yan Li
##' @keywords NetCDF fingerprinting
##' @import ncdf4
##' @importFrom utils head tail
##' @importFrom CFtime CFtime as_timestamp CFparse
##' @importFrom sp point.in.polygon
##' @export fpPrep

fpPrep <- function(datafile, variable, region = "GL", target.year,
                   average = 5, reference = c(1961, 1990),
                   regridding = NULL) {
  ## match args
  Call <- match.call()
  
  observation.file <- nc_open(datafile)
  # ncvar_get(observation, "latitude")
  # ncvar_get(observation, "longitude")
  
  ## get the array of observations
  observations <- ncvar_get(observation.file, variable)
  
  ## get the time
  time <- ncvar_get(observation.file,"time")
  tunits <- ncatt_get(observation.file,"time","units")
  ## convert time to CFtime class
  cftime <- CFtime(tunits$value, calendar = "proleptic_gregorian", time)
  ## get character-string times
  timestamps <- as_timestamp(cftime)
  ## parse the string into date components
  time_cf <- CFparse(cftime, timestamps)
  
  ## get the data year
  data.year <- c(head(time_cf$year, 1), tail(time_cf$year, 1))
  ## get the target year
  starting.year  <- target.year[1]; ending.year  <- target.year[2];
  
  if(ending.year > data.year[2]) {
    stop("The selected year exceeds the data timestamps")
  }
  
  ## get the target data order
  time_bnd <- seq(from = (starting.year - data.year[1]) * 12 + 1,
                  to = (ending.year - data.year[1] + 1) * 12)
  ## get the global data
  observation.GL <- read.obs(observations, time_bnd = time_bnd, year.step = average)
  # dim(observation.GL)
  # save(observation.GL, file=file.path(outDir, "observations.rdata"))
  ## get the longitude and latitude
  longitude <- ncvar_get(observation.file,"longitude")
  latitude <- ncvar_get(observation.file,"latitude")
  ## set the grid step
  if(is.null(regridding)) {
    step <- NULL
  } else {
    step <- c(regridding[1] / diff(longitude)[1],
              regridding[2] / diff(latitude)[1])
  }
  
  ## get the region
  if(region == "GL") {
    if(is.null(regridding)) {
      region <- rbind(c(-180, -90),
                      c(-180, 90),
                      c(180, 90),
                      c(180, -90))
    } else {
      region <- c(-180, 180, -90, 90)
    }
  }
  ## get the output
  output <- procDat(data.arr = observation.GL,
                    step = step,
                    region = region,
                    longitude = longitude, latitude = latitude)
  return(output)
}

#### function for processing the regional and global dataset
procDat <- function(data.arr, step = c(40 / 5, 30 / 5), 
                    region = c(-180, 180, -90, 90),
                    longitude, latitude) {
  ## expand the grid for longitude and latitude
  exp.grid <- expand.grid(longitude, latitude)
  ## output
  output <- NULL
  if(is.null(step)) {
    ## get the indicator for the data point
    Ind <- which(point.in.polygon(point.x = exp.grid[, 1], point.y = exp.grid[, 2],
                                      pol.x = region[, 1], pol.y = region[, 2]) == 1)
    ## take the average over the data
    for(i in 1:dim(data.arr)[3]) {
      tmp <- data.arr[, , i][Ind]
      if(length(Ind) == 112) {
        tmp <- colMeans(matrix(tmp, nrow = 2), na.rm = TRUE)
      }
      output <- cbind(output, tmp)
    }
  } else {
    Ind.row <-which(longitude > region[1] & longitude < region[2])
    Ind.col <-which(latitude > region[3] & latitude < region[4])
    for(i in 1:dim(data.arr)[3]) {
      tmp <- data.arr[Ind.row, Ind.col, i]
      tmp.vector <- NULL
      for(p in seq(1, nrow(tmp), by = step[1])) {
        for(q in seq(1, ncol(tmp), by = step[2])) {
          tmp.matrix <- tmp[p:(p+step[1] - 1), q:(q + step[2] - 1)]
          if(sum(is.na(tmp.matrix)) == length(tmp.matrix)) {
            tmp.vector <- c(tmp.vector, NA)
          } else {
            tmp.vector <- c(tmp.vector, mean(tmp.matrix, na.rm = TRUE))
          }
        }
      }
      output <- cbind(output, tmp.vector)
    }
  }
  return(output)
}

##############################################################################
############# inner functions 
##############################################################################
read.obs <- function(input, time_bnd, year.step, month = 12) {
  temp <- input[, , time_bnd]
  ## get annual means
  tmp <- vector("list", (length(time_bnd) / month))
  for(i in 1 : (length(time_bnd) / month)) {
    for(j in (i - 1) * month + 1:month) {
      tmp[[i]] <- cbind(tmp[[i]], as.vector(temp[, , j]))
    }
  }
  
  annual.list <- sapply(tmp,
                        function(inp) {
                          apply(inp, 1,
                                function(x) {
                                  if(sum(is.na(x)) > 3) {
                                    NA
                                  } else {
                                    mean(x, na.rm = TRUE)
                                  }
                                })
                        })
  
  bound <-   cbind(seq(1, length(time_bnd) / month, by = year.step),
                   seq(0, length(time_bnd) / month, by = year.step)[-1])
  
  ## get dimension
  slice.dim <- dim(temp)[3] / (year.step * month)
  row.dim <- dim(temp)[1]
  col.dim <- dim(temp)[2]
  output <- array(NA, dim = c(row.dim, col.dim, slice.dim))
  for(i in 1:nrow(bound)) {
    tmp.mat <- annual.list[, bound[i, 1] : bound[i, 2]]
    tmp.vec <- apply(tmp.mat, 1,
                     function(x) {
                       if(sum(is.na(x)) > 2) {
                         NA
                       } else {
                         mean(x, na.rm = TRUE)
                       }
                     })
    output[, , i] <- matrix(tmp.vec, nrow = row.dim, ncol = col.dim)
  }
  output
}

# dTrend <- function(x) {
#   timestep <- 1:length(x)
#   residuals(lm(x ~ timestep))
# }

rmModmean <- function(forcings) {
  models <-unique(colnames(forcings))
  output <- NULL
  for(mod.name in models) {
    m <- sum(colnames(forcings) == mod.name)
    tmp <- t(scale(t(forcings[, which(colnames(forcings) == mod.name)]),
                   scale = FALSE))[, - 1] * sqrt(m / (m - 1))
    output <- cbind(output, tmp)
  }
  output
}