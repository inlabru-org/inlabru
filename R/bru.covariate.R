#
# @aliases covariate
# @name covariate
# @export
# @param mesh
# @param values
# @param mesh.coords
# @param time.coords
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

covariate = function(data, predictor = NULL, mesh = NULL){
  if ( inherits(data, "SpatialPointsDataFrame") ) {

    cv = covdata.import(as.data.frame(data), 
                        deparse(substitute(predictor)), 
                        data = list(mesh = mesh, mesh.coords = coordnames(data)))
    cv
  }
}

# Construct a function that evaluates a covariate given a location
#
# @aliases evaluator
# @name evaluator
# @export
# @param covdata A covariate data set
# @return A function that returns covariate values for a given location
# @examples \\dontrun{}
# @author Fabian E. Bachl <\email{f.e.bachl@@bath.ac.uk}>

evaluator = function(cdata, method = NULL, ...){
  
  if ( class(cdata)[1] == "covdata") { 
    if ( is.null(method) ) { method = get.value }
    return( function(...){ method(cdata, data.frame(...)) } )
  }  
  
  else if (class(cdata)[1] == "SpatialPolygonsDataFrame") {
    if ( is.null(method) ) { method = shapefile.to.covariate }
    return( method(cdata, ...) )
  }
  
}