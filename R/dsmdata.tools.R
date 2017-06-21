# iDistance project
# R functions to convert data formatted for dsm to a data format for iDistance

gap.in.segments.f <- function(seg=NULL,geometry="euc") {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  #
  # Identify where there is a gap in the effort along a transect (so the segments do not join together)
  # Create an indicator, gap=0 can join to previous segment, gap=1 cannot join
  # To get distance between points, if geometry='geo' uses columns latitude/longitude and great circle distances,
  #    if geometry='euc' uses x/y and trig 
  
  # Segments
  num.seg <- dim(seg)[1]
  # Initialise indicator variable
  gap <- rep(0,num.seg)
  
  for (i in 1:(num.seg-1)) {
    # Only need to check segments if the same transect
    if (seg$Transect.Label[i]==seg$Transect.Label[i+1]) {
      # Check distance between segments
      if (geometry=="euc") dist <- euc.distance.f(seg$x[i],seg$y[i],seg$x[i+1],seg$y[i+1])
      if (geometry=="geo") dist <- geo.distance.f(seg$longitude[i],seg$latitude[i],seg$longitude[i+1],seg$latitude[i+1])
      # Check the length of the effort (half of one segment + half of the other)
      effort.dist <- (seg$Effort[i]/2) + (seg$Effort[i+1]/2) 
      if (dist>effort.dist) gap[i+1] <- 1
    } # End if
  } # End of segments
  
  gap
  
}

# ---------------------------------------------------------- 

define.blocks.f <- function(seg=NULL,covar.col=NULL,geometry="euc") {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  # samoe for: covar.col=covariate.columns
  #
  # Define blocks - adjoining segments which can be combined because the covariates are the same
  # covar.cols = number of the column in the dataframe containing covariates
  #              thius can be a list of columns (i.e. c(2,6,7)) 
  
  num.covar <- length(covar.col)
  if (is.na(covar.col)) print("No covariates used to combined segments - using transects")
  
  # Number of segments
  num.seg <- dim(seg)[1]
  
  # Initialise things
  seg$Block.Label <- rep(NA,num.seg)
  j <- 1
  # First segment in data is block 1
  seg$Block.Label[1] <- j 
  
  # See if there is a gap in search effort along transects
  gap <- gap.in.segments.f(seg=seg,geometry=geometry)
  
  # Loop through remaining segments
  for (i in 2:num.seg) {
    # Update block if transect changes from previous segment
    if (seg$Transect.Label[i] != seg$Transect.Label[i-1]) j <- j + 1
    # Update block if there is a gap from previous segment  
    if (gap[i]==1) j <- j + 1
    # If same transect, then update block if covars change from previous segment
    if (seg$Transect.Label[i] == seg$Transect.Label[i-1]) {
      chg <- 0
      if (!is.na(covar.col)) {
        # Loop through all covariates
        for (k in 1:num.covar) {
          if ( seg[i,covar.col[k]] != seg[i-1,covar.col[k]] ) chg <- 1
        } # End of covars
      } # End of if covars
      j <- j + chg
    } # End of same transect
    seg$Block.Label[i] <- j 
  } # End of segments
  
  seg
}

# -------------------------------------------------------

get.blocks.f <- function(seg=NULL,geometry="euc") {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  #
  # Create a dataset for the blocks
  # Segments must contain column called Block.Label, Effort
  # geometry defines how the geometry is measured, either euclidean space (x, y) or geometic coords (lon, lat)
  
  name.blocks <- unique(seg$Block.Label)
  num.blocks <- length(name.blocks)
  
  blocks <- NULL
  for (i in 1:num.blocks) {
    temp <- seg[seg$Block.Label==name.blocks[i], ]
    num.seg <- dim(temp)[1]
    # Retain first line of each block
    first <- temp[1, ]
    # Add on end of segments
    if (geometry=="euc") {
      first$end.x <- temp$end.x[num.seg]
      first$end.y <- temp$end.y[num.seg]
    }
    if (geometry=="geo") {
      first$end.lon <- temp$end.lon[num.seg]
      first$end.lat <- temp$end.lat[num.seg]
    }
    # Total effort for block
    first$Effort <- sum(temp$Effort)
    if (i==1) blocks <- first
    if (i>1) blocks <- rbind(blocks,first)
    
  } # End of blocks
  
  # Tidy up - get rid of unnecessary columns
  exc.labels <- c("quadrant","angle","what.angle","x","y","latitude","longitude")
  col.names <- names(blocks)
  col.names <- !is.element(col.names,exc.labels)
  #print(col.names)
  blocks <- blocks[ ,col.names]
  
  blocks
  
} 

# -------------------------------------------------------

add.labels.to.obs.f <- function(dists=NULL,obs=NULL,seg=NULL) {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  # same for: dists=distances,obs=obsservations
  #
  # Add segment and block labels to observations
  # distance data and observation data MUST be in the same order
  
  # Check same number of observations in each dataframe
  num.dists <- dim(dists)[1]
  if (num.dists!=dim(obs)[1]) print("Perp distance data and observation data different number of records") 
  
  dists$Sample.Label <- obs$Sample.Label
  dists$Block.Label <- rep(NA,num.dists)
  
  # Now add blocks - don't merge because merge sorts the data
  for (i in 1:num.dists) {
    temp <- seg[seg$Sample.Label==dists$Sample.Label[i], ]
    # Should only have one record - check
    if (dim(temp)[1]>1) print(paste("More than one segment chosen for observation, object=",dists$object[i]))
    dists$Block.Label[i] <- temp$Block.Label[1]
  } # End of observations
  
  dists
  
}

# -------------------------------------------------------

combine.dsmdata.f <- function(blocks=NULL,dists=NULL) {
  #
  # NOTE: old problematic parameteriation: blocks=blockdata,dists=distdata
  #
  # Combine segments and perp distances into one dataframe for iDistance
  
  num.blocks <- dim(blocks)[1]
  
  # Get names of columns in distance data that also occur in blocks
  cols.dists <- names(dists)
  
  # Get rid of labels that may be repeated in both dataframes
  exc.labels <- c("Sample.Label","Block.Label","Transect.Label","Effort")
  cols.dist <- !is.element(cols.dists,exc.labels)
  one.dists.small <- dists[1,cols.dist]
  # Create missing distance data
  one.dists.small[1, ] <- NA
  
  # Loop through all blocks adding sightings if there are any in block
  k <- 0
  for (i in 1:num.blocks) {
    # Get any sightings for the block
    one.block <- blocks[i, ]
    temp <- dists[dists$Block.Label==blocks$Block.Label[i], ]
    num.temp <- dim(temp)[1]
    # Sightings exist for block so add info together
    if (num.temp>0) {
      # Get rid of columns that will be repeated
      dists.small <- temp[ ,cols.dist]
      for (j in 1:num.temp) {
        block.dist <- cbind(one.block,dists.small[j, ])
        if (j==1) comb <- rbind(cbind(one.block, one.dists.small), block.dist)
        if (j>1) comb <- rbind(comb,block.dist)
      } # End of sightings
    } 
    # No sightings for block and so combine with missing distance data
    if (num.temp==0) comb <- cbind(one.block,one.dists.small)
    # Add blocks together
    if (i==1) all <- comb
    if (i>1) all <- rbind(all,comb)
  } # End of blocks
  
  # Tidy up - get rid of unnecessary columns
  exc.labels <- c("quadrant","angle","what.angle","Sample.Label")
  col.names <- names(all)
  col.names <- !is.element(col.names,exc.labels)
  all <- all[ ,col.names]
  
  # Rename columns
  col1 <- match("Transect.Label",names(all))
  col2 <- match("Block.Label",names(all))
  col3 <- match("object",names(all))
  names(all)[c(col1,col2,col3)] <- c("trans","seg","det")
  
  all
  
}

# -------------------------------------------------------

get.direction.unit.f <- function(data=NULL,is.blocks=T,geometry="euc") {
  
  # NOTE: old param: data=data
  
  # Get the quadrant and angle (from 0 to 360o) using blocks (is.block=T) or transects (is.block=F)
  # NEED TO CHECK THAT THIS IS CORRECT FOR UNITS IN QUADRANT = 7
  
  if (is.blocks) {
    data$Unit <- data$Block.Label
    name.label <- "Block.Label"
  }
  else {
    data$Unit <- data$Transect.Label
    name.label <- "Transect.Label"
  }
  
  if (geometry=="euc") {
    data$new.x <- data$x
    data$new.y <- data$y
  } 
  if (geometry=="geo") {
    data$new.x <- data$longitude
    data$new.y <- data$latitude
  } 
  
  name.unit <- unique(data$Unit)
  num.unit <- length(name.unit)
  
  unit <- NULL
  for (i in 1:num.unit) {
    temp <- data[data$Unit==name.unit[i], ]
    num.temp <- dim(temp)[1]
    unit$Unit[i] <- as.character(name.unit[i])
    quad <- get.quadrant.f(temp$new.x[1],temp$new.y[1],temp$new.x[num.temp],temp$new.y[num.temp])
    unit$quadrant[i] <- quad
    diff.x <- temp$new.x[num.temp] - temp$new.x[1]
    diff.y <- temp$new.y[num.temp] - temp$new.y[1]
    #  what.angle <- atan2(diff.y,diff.x) * (180/pi)
    what.angle <- atan(diff.y/diff.x) * (180/pi)
    if (quad==1) unit$angle[i] <- 0
    if (quad==2) unit$angle[i] <- 90
    if (quad==3) unit$angle[i] <- 180
    if (quad==4) unit$angle[i] <- 270
    if (quad==5) unit$angle[i] <- 90 - abs(what.angle)
    if (quad==6) unit$angle[i] <- 90 + abs(what.angle)
    if (quad==7) unit$angle[i] <- 180 + abs(what.angle)
    if (quad==8) unit$angle[i] <- 270 + abs(what.angle)
    unit$what.angle[i] <- what.angle
  } # End of units
  
  unit <- data.frame(unit)
  
  # Rename unit
  names(unit)[1] <- name.label
  
  unit
  
}

# -------------------------------------------------------

get.direction.segment.f <- function(data=NULL,geometry="euc") {
  
  # NOTE: old param: data=data
  
  # Get the quadrant and angle (from 0 to 360o clockwise) for each segment
  # Last segment in transect is assumed to be in same direction as penultimate segment.
  # Quadrant=Direction of travel; 1=N, 2=E, 3=S, 4=W, 5=NE, 6=SE, 7=SW, 8=NW
  
  if (geometry=="euc") {
    data$new.x <- data$x
    data$new.y <- data$y
  } 
  if (geometry=="geo") {
    data$new.x <- data$longitude
    data$new.y <- data$latitude
  } 
  
  num.unit <- dim(data)[1]
  
  # Initialise new vars
  data$quadrant <- rep(NA,num.unit)
  data$angle <- rep(NA,num.unit)
  #data$what.angle <- rep(NA,num.unit)
  # Direction in reverse - this will be used to refine start and end points (not done at present) 
  data$quadrant.r <- rep(NA,num.unit)
  data$angle.r <- rep(NA,num.unit)
  
  for (i in 1:(num.unit-1)) {
    # Get pair of segments
    j <- i + 1
    temp <- data[i:j, ]
    # Get quadrant on a compass
    quad <- get.quadrant.f(temp$new.x[1],temp$new.y[1],temp$new.x[2],temp$new.y[2])
    data$quadrant[i] <- quad
    diff.x <- temp$new.x[2] - temp$new.x[1]
    diff.y <- temp$new.y[2] - temp$new.y[1]
    data$angle[i] <- what.angle.f(dy=diff.y,dx=diff.x,quad=quad)
    # Check if points are on the same transect
    if (temp$Transect.Label[1]!=temp$Transect.Label[2]) {
      # If not the same transect, then check if point is on the same transect as the previous point
      # If so, must be the last point on the transect so set values to be the same as last but one point on transect
      if (data$Transect.Label[i]==data$Transect.Label[i-1]) {
        data$quadrant[i] <- data$quadrant[i-1]
        data$angle[i] <- data$angle[i-1]
        #      data$what.angle[i] <- data$what.angle[i-1]
      }
      if (data$Transect.Label[i]!=data$Transect.Label[i-1]) {
        print(paste("Only one segment in transect",data$Transect.Label[i]))
      }
    }
  } # End of units
  # Make last segment of data same as last but one
  if (data$Transect.Label[num.unit]==data$Transect.Label[num.unit-1]) {
    data$quadrant[num.unit] <- data$quadrant[num.unit-1]
    data$angle[num.unit] <- data$angle[num.unit-1]
    #  data$what.angle[num.unit] <- data$what.angle[num.unit-1]
  }
  if (data$Transect.Label[num.unit]!=data$Transect.Label[num.unit-1]) print("Last segment on its own")
  
  # Now do the same again only for points going in reverse
  # THIS BIT NEEDS TO BE CHECKED!
  for (i in num.unit:2) {
    # Get pair of segments
    j <- i - 1
    temp <- rbind(data[i, ],data[j, ])
    # Get quadrant on a compass
    quad <- get.quadrant.f(temp$new.x[1],temp$new.y[1],temp$new.x[2],temp$new.y[2])
    data$quadrant.r[i] <- quad
    diff.x <- temp$new.x[2] - temp$new.x[1]
    diff.y <- temp$new.y[2] - temp$new.y[1]
    data$angle.r[i] <- what.angle.f(dy=diff.y,dx=diff.x,quad=quad)
    # Check if the same transect
    if (temp$Transect.Label[1]!=temp$Transect.Label[2]) {
      # Check if same as previous transect
      if (data$Transect.Label[i]==data$Transect.Label[i-1]) {
        data$quadrant.r[i] <- data$quadrant.r[i-1]
        data$angle.r[i] <- data$angle.r[i-1]
        #      data$what.angle[i] <- data$what.angle[i-1]
      }
      if (data$Transect.Label[i]!=data$Transect.Label[i+1]) {
        print(paste("Only one segment in transect",data$Transect.Label[i]))
      }
    }
  } # End of units
  
  # Make first segment of data same as second but one
  if (data$Transect.Label[1]==data$Transect.Label[2]) {
    data$quadrant.r[1] <- data$quadrant.r[2]
    data$angle.r[1] <- data$angle.r[2]
  }
  
  # Tidy up 
  exc.labels <- c("new.x","new.y")
  col.names <- names(data)
  col.names <- !is.element(col.names,exc.labels)
  data <- data[ ,col.names]
  
  data
}


# -------------------------------------------------------  

start.end.points.segments.f <- function(seg=NULL,use.tran=FALSE,tran=NULL,geometry="euc") {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  # same for: tran=transect.quadrant
  #
  # Calculate the start and end points of a block of segments
  # Start and end points based on half the length of the first and last segments in a block
  # Direction is given by direction of segment (use.tran=FALSE) or transect (use.tran=TRUE)
  #   If use.tran=TRUE must supply dataframe containing direction of transects
  #   If use.tran=FALSE, segments must contain angle and quadrant of travel
  
  num.seg <- dim(seg)[1]
  
  for (i in 1:num.seg) {
    # Get the angle and quadrant for segment from transects, otherwise info already in segments
    if (use.tran) {
      temp <- tran[tran$Transect.Label==seg$Transect.Label[i], ]
      seg$quadrant[i] <- temp$quadrant[1]
      seg$angle[i] <- temp$angle[1]
    }
    seg.half.len <- seg$Effort[i]/2
    if (geometry=="euc") {
      if (seg$quadrant[i]==1) {
        seg$start.x[i] <- seg$x[i] 
        seg$start.y[i] <- seg$y[i] - seg.half.len
        seg$end.x[i] <- seg$x[i] 
        seg$end.y[i] <- seg$y[i] + seg.half.len
      }
      if (seg$quadrant[i]==2) {
        seg$start.x[i] <- seg$x[i] - seg.half.len
        seg$start.y[i] <- seg$y[i] 
        seg$end.x[i] <- seg$x[i] + seg.half.len
        seg$end.y[i] <- seg$y[i] 
      }
      if (seg$quadrant[i]==3) {
        seg$start.x[i] <- seg$x[i] 
        seg$start.y[i] <- seg$y[i] + seg.half.len
        seg$end.x[i] <- seg$x[i] 
        seg$end.y[i] <- seg$y[i] - seg.half.len
      }
      if (seg$quadrant[i]==4) {
        seg$start.x[i] <- seg$x[i] + seg.half.len
        seg$start.y[i] <- seg$y[i] 
        seg$end.x[i] <- seg$x[i] - seg.half.len
        seg$end.y[i] <- seg$y[i] 
      }
      # Now get width (tri.side element 1) and height (2) of triangles
      if (seg$quadrant[i]==5) {
        angle <- 90 - seg$angle[i]
        tri.sides <- get.triangle.sides.f(seg.len=seg.half.len,angle=angle)
        seg$start.x[i] <- seg$x[i] - tri.sides[2]
        seg$start.y[i] <- seg$y[i] - tri.sides[1]
        seg$end.x[i] <- seg$x[i] + tri.sides[2]
        seg$end.y[i] <- seg$y[i] + tri.sides[1]
      }
      if (seg$quadrant[i]==6) {
        angle <- 180.0 - seg$angle[i]
        tri.sides <- get.triangle.sides.f(seg.len=seg.half.len,angle=angle)
        seg$start.x[i] <- seg$x[i] - tri.sides[1]
        seg$start.y[i] <- seg$y[i] + tri.sides[2]
        seg$end.x[i] <- seg$x[i] + tri.sides[1]
        seg$end.y[i] <- seg$y[i] - tri.sides[2]
      }
      if (seg$quadrant[i]==7) {
        angle <- seg$angle[i] - 180.0
        tri.sides <- get.triangle.sides.f(seg.len=seg.half.len,angle=angle)
        seg$start.x[i] <- seg$x[i] + tri.sides[1]
        seg$start.y[i] <- seg$y[i] + tri.sides[2]
        seg$end.x[i] <- seg$x[i] - tri.sides[1]
        seg$end.y[i] <- seg$y[i] - tri.sides[2]
      }
      if (seg$quadrant[i]==8) {
        angle <- 360 - seg$angle[i]
        tri.sides <- get.triangle.sides.f(seg.len=seg.half.len,angle=angle)
        seg$start.x[i] <- seg$x[i] + tri.sides[1]
        seg$start.y[i] <- seg$y[i] - tri.sides[2]
        seg$end.x[i] <- seg$x[i] - tri.sides[1]
        seg$end.y[i] <- seg$y[i] + tri.sides[2]
      }
      
    } # End of if euc
    if (geometry=="geo") {
      print("This option is not implemented yet - use geometry=euc!")
    }
  } # End of segments
  
  # Change name of x,y cols
  if (geometry=="euc") {
    x.col <- match("x",names(seg))
    y.col <- match("y",names(seg))
    names(seg)[c(x.col,y.col)] <- c("mid.x","mid.y")
  }
  
  seg
  
} # End of function

# -----------------------------------

get.quadrant.f <- function(start.x,start.y,end.x,end.y,tol=0.0000001) {
  # Get the quadrant of the points
  
  x.diff <- (start.x - end.x)
  y.diff <- (start.y - end.y)
  
  if (abs(x.diff)<tol & abs(y.diff)<tol) print("Segments on top of each other!")
  
  quad <- NA
  
  if (abs(x.diff)<tol) {
    # Points N
    if (end.y>start.y) quad <- 1
    # S
    if (end.y<start.y) quad <- 3
  }
  if (abs(y.diff)<tol) {
    # Points E 
    if (end.x>start.x) quad <- 2
    # W
    if (end.x<start.x) quad <- 4
  }
  # NE
  if (end.x>start.x & end.y>start.y) quad <- 5
  # SE
  if (end.x>start.x & end.y<start.y) quad <- 6
  # SW
  if (end.x<start.x & end.y<start.y) quad <- 7
  # NW 
  if (end.x<start.x & end.y>start.y) quad <- 8
  
  if (is.na(quad)) print("Quadrant not assigned")
  
  quad
}


# -----------------------------------------------

euc.distance.f <- function(x1,y1,x2,y2) {
  # Calculate the Euclidean distance between two points, (x1,y1) and (x2,y2)
  
  x.diff <- abs(x1 - x2)
  y.diff <- abs(y1 - y2)
  
  hypot <- sqrt(x.diff^2 + y.diff^2)
  
  hypot
}

# -------------------------------------------------

geo.distance.f <- function (lon1,lat1,lon2,lat2) {
  # This function calculates distance in km between two points following 
  # a great circle route. Points defined as (lon1,lat1) and (lon2,lat2)
  # Coded from Jeff Laake's Excel Geometry Functions
  
  rad <- pi/180
  nm2km <- 1.852
  
  if ( (lat1 == lat2) & (lon1 == lon2) ) posdist <- 0
  else {
    rlat1 <- lat1*rad
    rlat2 <- lat2*rad
    
    rlon <- (lon2 - lon1) * rad
    posdist <- 60 * (1/rad) * acos( sin(rlat1) * sin(rlat2) + cos(rlat1) * cos(rlat2) * cos(rlon) )
  }
  
  # Convert to km
  posdist <- posdist * nm2km
  posdist
}   

# ---------------------------------------------------------------

get.triangle.sides.f <- function(seg.len=NULL,angle=NULL) {
  
  # NOTE: old parameterization: seg.len=half.segment.length,angle=angle
  
  # Calculate the height and length of triangle given length of hypotenuse and angle
  
  deg2rad <- pi/180
  
  x <- seg.len * sin(angle*deg2rad)
  y <- seg.len * cos(angle*deg2rad)
  
  triangle <- c(x,y)
  
  triangle
  
}

# ------------------------------------------------------------------

generate.obs.location.f <- function(seg=NULL,dists=NULL,geometry="euc",do.plot=F) {
  #
  # NOTE the parameter seg used to be "seg=segments", which caused CRAN compatibility issues
  # same for dists=distance.data
  #
  # Generate coords for an observation if not contained in observations
  # Location is randomly generated along segment and random side of segment at perp distance of sighting
  # Side of line; 1=above and -1=below segment
  # Plots segment, position along segment (green dot) and at perp distance (red dot) from line
  
  deg2rad <- pi/180
  
  # Get the number of sightings
  num.sgt <- dim(dists)[1]
  
  # Create dataframe for coordinates
  new.sgt <- NULL
  new.sgt$x <- rep(NA,num.sgt)
  new.sgt$y <- rep(NA,num.sgt)
  new.sgt <- data.frame(new.sgt)
  
  if (do.plot)  stop("do.plot code has been removed")#par(ask=T)
  for (i in 1:num.sgt) {
    temp <- seg[seg$Sample.Label==dists$Sample.Label[i], ]
    # Randomly choose location along line
    new.coords <- get.point.along.segment.f(temp$start.x[1],temp$start.y[1],temp$end.x[1],temp$end.y[1],quad=temp$quadrant[1],seg.angle=temp$angle[1])
    new.x <- new.coords[1]
    new.y <- new.coords[2]
    # Choose side of line at random
    what.side <- sample(c(-1,1), 1)
    # Sort out coords of sighting
    pd <- dists$distance[i]
    sgt.coords <- get.coords.f(quad=temp$quadrant[1],alpha=temp$angle[1],new.x=new.x,new.y=new.y,pd=pd,side=what.side)
    # Save coordinates of sighting
    new.sgt$x[i] <- sgt.coords[1]
    new.sgt$y[i] <- sgt.coords[2]
    # Plot segment and location if requested
    if (do.plot) {
      stop("do.plot=TRUE. This code has been removed.")
      # plot.segment.f(temp$start.x,temp$start.y,temp$end.x,temp$end.y)
      # points(new.x,new.y,col=3,pch=16)
      # points(sgt.coords[1],sgt.coords[2],pch=19,col=2)
      # segments(new.x,new.y,sgt.coords[1],sgt.coords[2],lty=2)
      # print(paste("Sighting",i," perp. dist=",format(pd)))
      # title(temp$Sample.Label[1])
    }
    
  } # End of sgts
  
  new.sgt
  
}

# ---------------------------------

get.point.along.segment.f <- function(x1,y1,x2,y2,quad=NULL,seg.angle) {
  
  # NOTE: old parameterization: quad=quadrant
  
  # Get a random point along the segment (at which object was detected)
  # Segment end points defined by (x1,y1) and (x2,y2)
  # Direction of travel defined by quadrant and angle
  
  deg2rad <- pi/180
  if (quad==1) {
    new.x <- x1
    new.y <- runif(1,min=y1,max=y2)
  }
  if (quad==2) {
    new.x <- runif(1,min=x1,max=x2)
    new.y <- y1
  }
  if (quad==3) {
    new.x <- x1
    new.y <- runif(1,min=y2,max=y1)
  }
  if (quad==2) {
    new.x <- runif(1,min=x2,max=x1)
    new.y <- y1
  }
  if (quad==5) {
    theta <- 90 - seg.angle
    diffx <- x2 - x1
    diffy <- y2 - y1
    # Get length of hypotenuse
    hyp <- get.hypot.f(diffx,diffy)
    # Get random length along hypotenuse
    smallhyp <- runif(1,0,hyp)
    smallx <- cos(theta*deg2rad) * smallhyp
    smally <- sin(theta*deg2rad) * smallhyp
    new.x <- x1 + smallx
    new.y <- y1 + smally
  }
  if (quad==6) {
    theta <- 180 - seg.angle
    diffx <- x2 - x1
    diffy <- y1 - y2
    # Get length of hypotenuse
    hyp <- get.hypot.f(diffx,diffy)
    # Get random length along hypotenuse
    smallhyp <- runif(1,0,hyp)
    smallx <- sin(theta*deg2rad) * smallhyp
    smally <- cos(theta*deg2rad) * smallhyp
    new.x <- x1 + smallx
    new.y <- y1 - smally
  }
  if (quad==7) {
    theta <- seg.angle - 180
    diffx <- x2 - x1
    diffy <- y1 - y2
    # Get length of hypotenuse
    hyp <- get.hypot.f(diffx,diffy)
    # Get random length along hypotenuse
    smallhyp <- runif(1,0,hyp)
    smallx <- sin(theta*deg2rad) * smallhyp
    smally <- cos(theta*deg2rad) * smallhyp
    new.x <- x1 - smallx
    new.y <- y1 - smally
  }
  if (quad==8) {
    theta <- seg.angle - 270
    diffx <- x1 - x2
    diffy <- y1 - y2
    # Get length of hypotenuse
    hyp <- get.hypot.f(diffx,diffy)
    # Get random length along hypotenuse
    smallhyp <- runif(1,0,hyp)
    smallx <- cos(theta*deg2rad) * smallhyp
    smally <- sin(theta*deg2rad) * smallhyp
    new.x <- x1 - smallx
    new.y <- y1 + smally
  }
  
  new.coords <- c(new.x,new.y)
  new.coords
  
}

# --------------------------------------------------------------

get.coords.f <- function(quad=NULL,alpha=NULL,new.x=NULL,new.y=NULL,pd=NULL,side=NULL) {
  
  # NOTE old parameterization:
  # quad=quadrant,alpha=angle,new.x=x.reference,new.y=y.reference,pd=perp.dist,side=side.of.segment
  
  # Get the coordinates of the sighting based on 
  # 1. the reference point (new.x,new.y)
  # 2. the perp distance of the sighting
  # 3. which side of line sighting is 
  # 4. quadrant is quadrant on compass
  # 5. angle is angle from N going clockwise
  
  deg2rad <- pi/180
  
  # Quadrants 1 and 3 (i.e. heading directly N or S, respectively)
  if (quad==1 | quad==3) {
    sgt.x <- new.x + (side * pd)
    sgt.y <- new.y
  }
  # Quadrants 2 and 4 (i.e. heading directly E or W, respectively)
  if (quad==2 | quad==4) {
    sgt.x <- new.x 
    sgt.y <- new.y + (side * pd)
  }
  # Quadrant 5 (NE)
  if (quad==5) {
    theta <- 90 - alpha
    x1 <- sin(theta*deg2rad) * pd
    y1 <- sqrt(pd^2 - x1^2)
    if (side==1) {
      sgt.x <- new.x - x1
      sgt.y <- new.y + y1
    }
    if (side==(-1)) {
      sgt.x <- new.x + x1
      sgt.y <- new.y - y1
    }
  } # End of quad=5
  # Quadrant 6 (SE)
  if (quad==6) {
    theta <- alpha - 90
    x1 <- pd/tan(deg2rad*theta)
    y1 <- sin(theta*deg2rad) * x1
    x2 <- sqrt(pd^2 - y1^2)
    if (side==1) {
      sgt.x <- new.x + x2
      sgt.y <- new.y + y1
    }
    if (side==(-1)) {
      sgt.x <- new.x - x2
      sgt.y <- new.y - y1
    }
  } # End of quad=6
  # Quadrant 7 (SW)
  if (quad==7) {
    theta <- 270 - alpha
    x1 <- sin(theta*deg2rad) * pd
    y1 <- sqrt(pd^2 - x1^2)
    if (side==1) {
      sgt.x <- new.x - x1
      sgt.y <- new.y + y1
    }
    if (side==(-1)) {
      sgt.x <- new.x + x1
      sgt.y <- new.y - y1
    }
  } # End of quad=7
  # Quadrant 8 (NW)
  if (quad==8) {
    theta <- 360 - alpha
    x1 <- pd/tan(deg2rad*theta)
    x2 <- sin(theta*deg2rad) * x1
    y1 <- sqrt(pd^2 - x2^2)
    if (side==1) {
      sgt.x <- new.x + x2
      sgt.y <- new.y + y1
    }
    if (side==(-1)) {
      sgt.x <- new.x - x2
      sgt.y <- new.y - y1
    }
  } # End of quad=8
  
  sgt.coord <- c(sgt.x,sgt.y)
  sgt.coord
} 

# ------------------------------

get.hypot.f <- function(side1,side2) {
  # Get the length of the hypotenuse of triangle with sides of length side1 and side2
  hyp <- sqrt(side1^2 + side2^2)
  hyp
}

# -------------------------------

# plot.segment.f <- function(x1,y1,x2,y2) {
#   # Plots segment defined by (x1,y1) and (x2,y2)
#   
#   minx <- min(x1,x2)
#   maxx <- max(x1,x2)
#   miny <- min(y1,y2)
#   maxy <- max(y1,y2)
#   
#   plot(x1,y1,xlim=c(minx,maxx),ylim=c(miny,maxy),asp=1)
#   points(x2,y2,pch=2)
#   segments(x1,y1,x2,y2)
#   
# }

# ---------------------------------------------------------------------

what.angle.f <- function(dy=NULL,dx=NULL,quad=NULL) {
  
  # NOTE old parameterization: dy=diff.y,dx=diff.x,quad=quadrant
  
  # Get the angle between two points (from 0 to 360) using trig and the direction (quadrant)
  
  what.angle <- atan(dy/dx) * (180/pi)
  if (quad==1) angle <- 0
  if (quad==2) angle <- 90
  if (quad==3) angle <- 180
  if (quad==4) angle <- 270
  if (quad==5) angle <- 90 - abs(what.angle)
  if (quad==6) angle <- 90 + abs(what.angle)
  if (quad==7) angle <- 180 + (90 - abs(what.angle))
  if (quad==8) angle <- 270 + abs(what.angle)
  
  angle
}



