#ifndef __MYFLOAT_H__
#define __MYFLOAT_H__

//Defines the floating-point type myfloat using Boost extended-precision libraries.
//Precision (number of decimal digits in mantissa) is set at compile time and can
// be changed below or with a compiler flag.

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
//#include <boost/geometry/geometries/point.hpp>
//#include <boost/geometry/core/cs.hpp>
//#include <boost/geometry/algorithms/transform.hpp>
//#include <boost/geometry.hpp>

using boost::math::itrunc;
using boost::math::tools::epsilon;
using boost::math::spherical_harmonic_r;
using namespace boost::math::constants;

//number of digits in myfloat type
// Setting precision to 0 gives unlimited precision
#ifndef PRECISION
//#warning "Using default precision value: 20 digits. Set at compile time with flag -DPRECISION=num_digits"
#define PRECISION 20
#endif
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float< PRECISION > > myfloat;

#endif