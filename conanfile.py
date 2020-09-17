from conans import ConanFile, CMake


class S2letConan(ConanFile):
    name = "s2let"
    version = "2.1.0"
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
        "include/*",
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

    def requirements(self):
        self.requires("so3/1.2.1@AstroFizz/stable")
        if self.options.with_cfitsio:
            self.requires("cfitsio/3.480")

    @property
    def cmake(self):
        if not hasattr(self, "_cmake"):
            self._cmake = CMake(self)
            self._cmake.definitions["tests"] = True
            self._cmake.definitions["conan_deps"] = True
            self._cmake.definitions["python"] = False
            self._cmake.definitions["fPIC"] = self.options.fPIC
            self._cmake.configure(source_folder=".")
        return self._cmake

    def build(self):
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["s2let"]
