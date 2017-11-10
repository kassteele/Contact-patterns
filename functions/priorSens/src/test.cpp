#define _USE_MATH_DEFINES

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector coord_correc_func_cpp(double theta, List coordinate_correction){

  NumericVector result = NumericVector::create(1,1);

  double xpos = as<double>(coordinate_correction["xpos"]);
  double xneg = as<double>(coordinate_correction["xneg"]);
  double ypos = as<double>(coordinate_correction["ypos"]);
  double yneg = as<double>(coordinate_correction["yneg"]);

  if (theta <= -M_PI/2)
  {
    result[0] = xneg;
    result[1] = yneg;
  }
  else if (theta <= 0)
  {
    result[0] = xpos;
    result[1] = yneg;
  }
  else if (theta <= M_PI/2)
  {
    result[0] = xpos;
    result[1] = ypos;
  }
  else
  {
    result[0] = xneg;
    result[1] = ypos;
  }

   return result;

}
