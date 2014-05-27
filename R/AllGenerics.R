setGeneric(
           name = "RUVg",
           def = function(x, cIdx, k, drop=0, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8) {
             standardGeneric("RUVg")
           }
           )

setGeneric(
           name = "RUVs",
           def = function(x, cIdx, k, scIdx, round=TRUE, epsilon=1, tolerance=1e-8) {
             standardGeneric("RUVs")
           }
           )

setGeneric(
           name = "RUVr",
           def = function(x, cIdx, k, residuals, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8) {
             standardGeneric("RUVr")
           }
           )

