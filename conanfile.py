from conans import CMake, ConanFile


class S2letConan(ConanFile):
    name = "s2let"
    version = "2.2.2"
    license = "GPL-2.0"
    url = "https://github.com/astro-informatics/s2let"
    homepage = "https://github.com/astro-informatics/s2let"
    description = "Fast wavelet transforms on the sphere"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"fPIC": [True, False], "with_cfitsio": [True, False]}
    default_options = {"fPIC": True, "with_cfitsio": True}
    generators = "cmake"
    exports_sources = [
        "src/main/c/*",
        "src/test/c/*",
        "include/s2let/*",
        "CMakeLists.txt",
        "cmake/*.cmake",
    ]

    def configure(self):
        if self.settings.compiler == "Visual Studio":
            del self.options.fPIC
        self.options["so3"].fPIC = self.options.fPIC
        if self.options.with_cfitsio:
            self.options["cfitsio"].fPIC = self.options.fPIC
            self.options["cfitsio"].shared = False
        del self.settings.compiler.libcxx

    def requirements(self):
        location = "astro-informatics/stable" if self.in_local_cache else "user/testing"
        self.requires(f"so3/1.3.1@{location}")
        if self.options.with_cfitsio:
            self.requires("cfitsio/3.490")

    @property
    def cmake(self):
        if not hasattr(self, "_cmake"):
            self._cmake = CMake(self)
            self._cmake.definitions["tests"] = True
            self._cmake.definitions["conan_deps"] = True
            self._cmake.definitions["python"] = False
            self._cmake.definitions["fPIC"] = self.options.fPIC
            self._cmake.configure(build_folder="build")
        return self._cmake

    def build(self):
        from pathlib import Path

        path = Path(self.source_folder)
        build = Path(self.source_folder) / "build"
        build.mkdir(exist_ok=True)
        (path / "conanbuildinfo.cmake").rename(path / "build" / "conanbuildinfo.cmake")
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["s2let"]
