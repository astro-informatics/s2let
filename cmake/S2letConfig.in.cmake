set(S2LET_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(astro-informatics-so3 REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/sshtTargets.cmake")
set(S2LET_LIBRARIES s2let::s2let)

check_required_components(S2LET)