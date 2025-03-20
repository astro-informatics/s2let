include(FetchContent)

FetchContent_Declare(
  CMocka
  GIT_REPOSITORY https://gitlab.com/cmocka/cmocka.git
  GIT_TAG cmocka-1.1.5)

set(WITH_STATIC_LIB
    ON
    CACHE BOOL "CMocka: Build with a static library" FORCE)
set(WITH_CMOCKERY_SUPPORT
    OFF
    CACHE BOOL "CMocka: Install a cmockery header" FORCE)
set(PICKY_DEVELOPER
    OFF
    CACHE BOOL "CMocka: Build with picky developer flags" FORCE)
FetchContent_MakeAvailable("cmocka")
