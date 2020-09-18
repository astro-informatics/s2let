# S2LET: Fast wavelet transforms on the sphere

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://astro-informatics.github.io/s2let/
[bintray-img]: https://img.shields.io/bintray/v/mdavezac/AstroFizz/s2let:AstroFizz?label=C%20package
[bintray-url]: https://bintray.com/mdavezac/AstroFizz/s2let:AstroFizz/2.1.0:stable/link
[pypi-img]: https://badge.fury.io/py/pys2let.svg
[pypi-url]: https://badge.fury.io/py/pys2let

[![][docs-img]][docs-url]
[![][bintray-img]][bintray-url]
[![][pypi-img]][pypi-url]

## DESCRIPTION

S2LET provides functionality to perform fast and exact scale-discretised
wavelet transforms on the sphere.

## AUTHORS
- [B. Leistedt](www.ixkael.com/blog)
- M. Buttner
- [J. Y. H. Chan](www.mssl.ucl.ac.uk/~yhjc2/)
- [J. D. McEwen](www.jasonmcewen.org)

## REFERENCES

- J. D. McEwen, M. Büttner, B. Leistedt, H. V. Peiris and Y. Wiaux, 
  "Directional spin wavelets on the sphere", IEEE TSP, submitted, 2015 
  ([arXiv](http://arxiv.org/abs/1509.06749)).

- J. D. McEwen, M. Büttner, B. Leistedt, H. V. Peiris and Y. Wiaux, 
  "A novel sampling theorem on the rotation group", IEEE TSP, 2015
  ([arXiv](http://arxiv.org/abs/1508.03101)|[DOI](http://dx.doi.org/10.1109/LSP.2015.2490676)).

- B. Leistedt, J. D. McEwen, P. Vandergheynst and Y. Wiaux,
  "S2LET: A code to perform fast wavelet analysis on the sphere", 
  Astronomy & Astrophysics, 558, A128, 2013
  ([arXiv](http://arxiv.org/abs/1211.1680)|[DOI](http://dx.doi.org/10.1051/0004-6361/201220729)).

- J. D. McEwen, "Ridgelet transform on the sphere",
  IEEE TSP, submitted, 2015
  ([arXiv](http://arxiv.org/abs/1510.01595v1)).
     
- J. D. McEwen and Y. Wiaux, "A novel sampling theorem on the sphere",
- IEEE Trans. Sig. Proc., 59(12):5876-5887, 2011
  ([arXiv](http://arxiv.org/abs/1110.6298)|[DOI](http://dx.doi.org/10.1109/TSP.2011.2166394)).


- J. Y. H. Chan, B. Leidtedt, T.D. Kitching and J. D. McEwen, 
  "Second-generation curvelets on the sphere", 
  IEEE TSP, 65(1):5-14, 2017
  ([arXiv](http://arxiv.org/abs/1511.05578)|[DOI](http://dx.doi.org/10.1109/TSP.2016.2600506)).

- J. D. McEwen, G. Puy, J.-Ph. Thiran, P. Vandergheynst, D. Van De Ville, and Y. Wiaux,
  "Sparse image reconstruction on the sphere: implications of a new sampling theorem"
  IEEE Trans. Image Proc., 22(6):2275-2285, 2013
  ([arXiv](http://arxiv.org/abs/1205.1013)|[DOI](http://dx.doi.org/10.1109/TIP.2013.2249079)).

## INSTALLATION
The python package can be installed from the cloud with ``pip install pys2let``
or from a local repository for development with `pip install -e .[dev]`.

The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html):

Both can be installed using pip:

```bash
pip install conan cmake
```

Then s2let can be compiled with:

```bash
git clone https://github.com/astro-informatics/s2let.git
mkdir s2let/build && cd s2let/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -Dconan_deps=ON  -Dcfitsio=ON ..
make
make install
```

The above will also download all necessary dependencies.

Instructions for installing other languages can be found in docs/index.html.
