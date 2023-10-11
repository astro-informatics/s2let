#ifndef S2LET_VERSION_H
#define S2LET_VERSION_H
inline const char *s2let_version_string() { return "@PROJECT_VERSION@"; }
inline const char *s2let_info() {
  return "package:\n"
         "  name: S2LET\n"
         "  description: Scale-discretised axisymmetric and direction wavelets\n"
         "  authors:\n"
         "      - Boris Leistedt\n"
         "      - Martin Büttner\n"
         "      - Jennifer Chan\n"
         "      - Jason McEwen\n"
         "  license: GPL-3\n"
         "  url: https://astro-informatics.github.io/s2let\n"
         "  version: @PROJECT_VERSION@\n";
};
// clang-format off
inline int s2let_version_major() { return @PROJECT_VERSION_MAJOR@; }
inline int s2let_version_minor() { return @PROJECT_VERSION_MINOR@; }
inline int s2let_version_patch() { return @PROJECT_VERSION_PATCH@; }
#define S2LET_VERSION_MAJOR @PROJECT_VERSION_MAJOR@
#define S2LET_VERSION_MINOR @PROJECT_VERSION_MINOR@
#define S2LET_VERSION_PATCH @PROJECT_VERSION_PATCH@
// clang-format on
#endif
