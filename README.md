## Purpose
Calculation of Martin Luscher's zeta function defined in [1].

Functions are fully templated and can be called with any floating-point type. Currently Boost cpp_dec_float is used for many decimal digits of precision (set at compile time), but this can be swapped out easily for any built in or extended precision type.

The form of the derivatives of the heat kernel is from [2] eqs. 109,110.

Numerical advice is taken from [3] Appendix B, and unless otherwise stated all references in the source code are to there.

## Dependencies
- Boost [multiprecision](http://www.boost.org/doc/libs/1_63_0/libs/multiprecision/doc/html/index.html)
- Boost [math/special_functions](http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/special.html)
- Boost [math/tools/precision](http://www.boost.org/doc/libs/1_63_0/boost/math/tools/precision.hpp)
- [GMP](https://gmplib.org/)
- [MPFR](http://www.mpfr.org/)


## References
[1]: [M. Lüscher, Nucl. Phys. B354, 531](http://inspirehep.net/record/300613) ([PDF](http://www-library.desy.de/preparch/desy/postpr/1990/desy90-131.pdf))

[2]: [K. Rummukainen, Steven A. Gottlieb, Nucl. Phys B450, 397](http://inspirehep.net/record/393935) ([arXiv](https://arxiv.org/abs/hep-lat/9503028))

[3]: Pelissier, Craig S., Ph.D., The George Washington University, 2013, 134; 3544183 (Dissertation) ([link](http://pqdtopen.proquest.com/doc/1220841203.html?FMT=ABS))
