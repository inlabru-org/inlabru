# Old methods needed to import mexdolphins data

as.spatial.dsdata <- function(dset, cnames, crs) {
  as.detection <- function(...) {
    UseMethod("as.detection")
  }
  detdata <- function(...) {
    UseMethod("detdata")
  }
  segdata <- function(...) {
    UseMethod("segdata")
  }

  as.detection.dsdata <- function(dsdata, ...) {
    return(as.detection.effort(dsdata$effort, ...))
  }


  detdata.dsdata <- function(data, detection = NULL, ...) {
    if (is.null(detection)) {
      return(data$effort[as.detection(data)[, "start"], ])
    } else {
      return(data$effort[detection[, "start"], ])
    }
  }

  segdata.dsdata <- function(data, ...) {
    return(data$effort[is.na(data$effort$det), ])
  }


  # Extract detection pointers from effort data
  #
  # @aliases as.detection.effort
  # @export
  # @param effort Effort data set
  # @return sighting sightings
  # @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

  as.detection.effort <- function(effort) {
    idx <- which(!(is.na(effort$det)))
    det <- data.frame(start = idx, end = idx)
    class(det) <- c("detection", "data.frame")
    return(det)
  }


  # Convert dsdata into Spatial* objects
  #
  # @aliases as.spatial.dsdata
  # @export
  # @param dsdata dsdata object
  # @param cnames Column names of the coordinates, e.g. c("lon","lat")
  # @param crs Coordinate reference system, e.g. CRS("+proj=longlat")
  # @return list with spatial points, spatial lines and a mesh
  # @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

  # as.spatial.dsdata = function(dset, cnames, crs)

  # Extract detections and convert to SpatialPointsDataFrame
  points <- SpatialPointsDataFrame(
    coords = as.data.frame(detdata(dset))[, cnames],
    data = as.data.frame(detdata(dset))[, setdiff(names(detdata(dset)), cnames)],
    proj4string = crs
  )

  # Extract effort and convert to SpatialLinesDataFrame
  sp <- as.data.frame(segdata(dset)[, paste0("start.", cnames)])
  ep <- as.data.frame(segdata(dset)[, paste0("end.", cnames)])
  colnames(sp) <- cnames
  colnames(ep) <- cnames
  lilist <- lapply(seq_len(nrow(sp)), function(k) {
    Lines(list(Line(rbind(sp[k, ], ep[k, ]))), ID = k)
  })

  spl <- SpatialLines(lilist, proj4string = crs)
  df <- as.data.frame(segdata(dset))[, setdiff(names(segdata(dset)), c(paste0("start.", cnames), paste0("end.", cnames)))]
  rownames(df) <- seq_len(nrow(df))
  samplers <- SpatialLinesDataFrame(spl, data = df)

  # Return detections ("points") and effort ("samplers")
  list(points = points, samplers = samplers, mesh = dset$mesh)
}
