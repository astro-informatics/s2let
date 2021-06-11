from skbuild import setup

cmake_args = [
    "-DBUILD_TESTING:BOOL=OFF",
    "-Dconan_deps:BOOL=ON",
    "-Dcfitsio:BOOL=ON",
    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
]

build_requirements = [
    "setuptools",
    "wheel",
    "scikit-build",
    "cmake>=3.12",
    "ninja",
    "cython",
    "conan",
    "pip!=20.0.0,!=20.0.1",
]

setup(
    name="pys2let",
    version="2.2.3",
    author=["Boris Leistedt", "Martin Büttner", "Jennifer Chan", "Jason McEwen"],
    install_requires=["numpy"],
    extras_require={
        "dev": build_requirements + ["pytest", "black"],
        "plots": ["scipy", "healpy"],
    },
    description="Fast spin spherical transforms",
    url="http://astro-informatics.github.io/s2let/",
    package_dir={"pys2let": "src/main/pys2let"},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-2",
    packages=["pys2let"],
)
