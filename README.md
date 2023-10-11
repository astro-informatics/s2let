# S2LET: Fast wavelet transforms on the sphere

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://astro-informatics.github.io/s2let/
[conan-img]: https://img.shields.io/badge/ConanCenter-C%20Package-red.svg
[conan-url]: https://conan.io/center/s2let
[pypi-img]: https://badge.fury.io/py/pys2let.svg
[pypi-url]: https://badge.fury.io/py/pys2let
[codefactor-img]: https://www.codefactor.io/repository/github/astro-informatics/s2let/badge/main
[codefactor-url]: https://www.codefactor.io/repository/github/astro-informatics/s2let/overview/main

[![][docs-img]][docs-url]
[![][conan-img]][conan-url]
[![][pypi-img]][pypi-url]
![CMake Build](https://github.com/astro-informatics/s2let/workflows/CMake%20Build/badge.svg)
![Python Build](https://github.com/astro-informatics/s2let/workflows/Python%20Build/badge.svg)

## DESCRIPTION

S2LET provides functionality to perform fast and exact scale-discretised
wavelet transforms on the sphere.

## INSTALLATION

The python package, **pys2let**, is available on [pypi](https://pypi.org/project/pys2let/) and can be installed with:
 
 ```bash
 pip install pys2let
 ```

Alternatively, it can be installed from a local clone of the repository for development purposes by

 ```bash
 pip install -e .[dev]
 ```

The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html):

Both can be installed using pip:

```bash
pip install "conan<1" cmake
```

Then **S2LET** can be compiled with:

```bash
git clone https://github.com/astro-informatics/s2let.git
mkdir s2let/build && cd s2let/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -Dconan_deps=ON  -Dcfitsio=ON ..
make
make install
```

The above will also download all necessary dependencies.

Instructions for installing other languages can be found [here](https://astro-informatics.github.io/s2let/).


## DOCUMENTATION

Further documentation is available [here](https://astro-informatics.github.io/s2let/).

Usage for the python package is also given in the package docstring.


## REFERENCING

If you use **S2LET** for work that results in publication, please reference [http://github.com/astro-informatics/s2let](https://github.com/astro-informatics/s2let/) and cite the relevant academic paper(s):

- Y. Wiaux, J. D. McEwen,  P. Vandergheynst, O. Blanc,
  [Exact reconstruction with directional wavelets on the sphere](http://www.jasonmcewen.org/publication/wiaux-2007-sdw/), Mon. Not. Roy. Astron. Soc., 388(2):770-788, 2008. ([ArXiv](http://arxiv.org/abs/arXiv:0712.3519))
  | [DOI](http://dx.doi.org/10.1111/j.1365-2966.2008.13448.x))
- B. Leistedt, J. D. McEwen, P. Vandergheynst and Y. Wiaux, [S2LET: A code to perform fast wavelet analysis on the sphere](http://www.jasonmcewen.org/publication/leistedt-s-2-let-axisym/), Astronomy & Astrophysics, 558(A128):1-9, 2013 (http://arxiv.org/abs/1211.1680">ArXiv</a> | [DOI](http://dx.doi.org/10.1051/0004-6361/201220729)
- J. D. McEwen,  B. Leistedt, M. Büttner, H. V. Peiris, Y. Wiaux, [Directional spin wavelets on the sphere](http://www.jasonmcewen.org/publication/mcewen-s-2-let-spin/), IEEE Trans. Signal Proc., submitted, 2015 ([ArXiv](http://arxiv.org/abs/1509.06749)
-  J. D. McEwen, M. Price, [Ridgelet transform on the sphere](http://www.jasonmcewen.org/publication/mcewen-s-2-let-ridgelets/), 27th European Signal Processing Conference (EUSIPCO), 2019 ([ArXiv](http://arxiv.org/abs/1510.01595v1) | [DOI](http://dx.doi.org/10.23919/EUSIPCO.2019.8903034))
- J. Y. H. Chan, B. Leistedt, T. D. Kitching, J. D. McEwen, [Second-generation curvelets on the sphere](http://www.jasonmcewen.org/publication/chan-s-2-let-curvelets/), IEEE Trans. Signal Proc., 65(1):5-14, 2017 ([ArXiv](http://arxiv.org/abs/1511.05578) | [DOI](http://dx.doi.org/10.1109/TSP.2016.2600506))
- J. D. McEwen,  C. Durastanti, Y. Wiaux, [Localisation of directional scale-discretised wavelets on the sphere](http://www.jasonmcewen.org/publication/mcewen-s-2-let-localisation/), Applied Comput. Harm. Anal., 44(1), 59-88, 2018 ([ArXiv](http://arxiv.org/abs/1509.06749) | [DOI](http://dx.doi.org/10.1016/j.acha.2016.03.009))

You may also like to consider citing the following papers on which the fast algorithms of S2LET are based:
- J. D. McEwen, M. B&uuml;ttner, B. Leistedt, H. V. Peiris, Y. Wiaux, [A novel sampling theorem on the rotation group](http://ieeexplore.ieee.org/document/7298431/), IEEE Sig. Proc. Let., 22(12):2425-2429, 2015 ([ArXiv](http://arxiv.org/abs/1508.03101) | [DOI](http://dx.doi.org/10.1109/LSP.2015.2490676))
- J. D. McEwen and Y. Wiaux, <a href="http://www.jasonmcewen.org/publication/mcewen-so-3/">A
 novel sampling theorem on the sphere</a>, IEEE Trans. Signal Proc., 59, 5876-5887, 2011 ([ArXiv](http://arxiv.org/abs/1110.6298)
 | [DOI](http://dx.doi.org/10.1109/TSP.2011.2166394))




## LICENSE

S2LET is released under the GPL-3 license.  For further details see 
[LICENSE.txt](https://github.com/astro-informatics/s2let/blob/main/LICENSE).

## AUTHORS

**S2LET** was initially developed by Boris Leistedt, Martin Büttner, and [Jason McEwen](http://www.jasonmcewen.org/) but significant contributors have since been made by a number of [others](https://github.com/astro-informatics/s2let/graphs/contributors).
