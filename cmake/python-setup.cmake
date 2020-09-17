macro(find_python)
  set(OLD_CMAKE_FIND_FRAMEWORK "${CMAKE_FIND_FRAMEWORK}")
  set(CMAKE_FIND_FRAMEWORK LAST)
  find_package(Python3 REQUIRED COMPONENTS Development Interpreter)
  set(CMAKE_FIND_FRAMEWORK "${OLD_CMAKE_FIND_FRAMEWORK}")
  unset(OLD_CMAKE_FIND_FRAMEWORK)

  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "Path to python executable")
endmacro()

function(setup_skbuild)
  if(SKBUILD)
    return()
  endif()
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import skbuild; print(skbuild.__file__)"
    OUTPUT_VARIABLE SKBUILD_LOCATION
    RESULT_VARIABLE SKBUILD_FOUND
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT SKBUILD_FOUND EQUAL 0 OR NOT SKBUILD_LOCATION)
    message(FATAL_ERROR "Could not find scikit-build")
  else()
    set(SKBUILD_FOUND True)
    get_filename_component(SKBUILD_LOCATION "${SKBUILD_LOCATION}" DIRECTORY)
  endif()
  message(STATUS "Found skbuild at ${SKBUILD_LOCATION}")
  list(APPEND CMAKE_MODULE_PATH "${SKBUILD_LOCATION}/resources/cmake")
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} PARENT_SCOPE)
endfunction()
