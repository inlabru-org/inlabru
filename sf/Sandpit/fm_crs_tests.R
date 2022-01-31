# Testing for fmesher functions using sf objects.
# For testing while editing code.  Eventually these will be the
# basis for proper package tests.  For now should work
# standalone by sourcing the relevant fmesher_*.R scripts directly

library(inlabru)
library(INLA)
library(sf)
source(here::here("R", "fmesher_crs.R"))

#### testing for fm_as_inla_mesh_segment.sf ###

## sfc_POINT ##

# compare inla.mesh.segment with matrix input
# to fm_as_inla_mesh_segment with sf point input

# matrix version
loc.bnd = matrix(c(0,0, 1,0, 1,1, 0,1), 4, 2, byrow=TRUE)
segm.bnd = inla.mesh.segment(loc.bnd, 
                             idx = seq_len(nrow(loc.bnd)),  
                             is.bnd = TRUE)

# note, this returns different indexing if idx is not set manually here. 
# Setting this to seq_len(nrow(loc.bnd)) makes the test work.

# From inla.mesh.segment documentation:
#   If ‘is.bnd==TRUE’, defaults to linking all
#   the points in ‘loc’, as ‘c(1:nrow(loc),1L)’, otherwise
#   ‘1:nrow(loc)’.
#
# However, fm_as_inla_mesh_segment.SpatialPoints uses seq_len(nrow(loc)) even if
# is.bnd == TRUE
#
# I am not sure if this is the desired behaviour.

# sf version
loc.sf = st_as_sf(as.data.frame(loc.bnd),
                  coords = c(1,2))

segm.bnd.sf = fm_as_inla_mesh_segment(loc.sf)

identical(segm.bnd.sf, segm.bnd)

