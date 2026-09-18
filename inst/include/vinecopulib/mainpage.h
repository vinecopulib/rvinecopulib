/** @file mainpage.h
 *
 */

/** @mainpage vinecopulib: a C++ library for vine copula modeling

@authors Thomas Nagler and Thibault Vatter

@section what-vc What are vine copulas?

Vine copulas are a flexible class of dependence models consisting of bivariate
building blocks (see e.g.,
[Aas et al., 2009](https://mediatum.ub.tum.de/doc/1083600/1083600.pdf)).
You can find a comprehensive list of publications and other materials on
[vine-copula.org](http://www.statistics.ma.tum.de/en/research/vine-copula-models/).


@section what-vcl What is vinecopulib?

[vinecopulib](https://vinecopulib.github.io/vinecopulib/) is a header-only C++
library for vine copula models based on
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It provides
high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are
- a stand-alone C++ library with interfaces to both R and Python,
- a sleeker and more modern API,
- shorter runtimes and lower memory consumption, especially in high dimensions,
- nonparametric and multi-parameter families.


@section Contact
If you have any questions regarding the library, feel free to
[open an issue](https://github.com/vinecopulib/vinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.

*/
