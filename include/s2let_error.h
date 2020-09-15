#ifndef S2LET_ERROR
#define S2LET_ERROR

#include <stdio.h>
#include <ssht/ssht.h>

#ifdef __cplusplus
extern "C" {
#endif

// Put this macro in a block so that it can be used with single-line
// if-statements.
#define S2LET_ERROR_GENERIC(comment)                  \
{                                                                       \
  printf("ERROR: %s.\n", comment);                  \
  printf("ERROR: %s <%s> %s %s %s %d.\n",               \
     "Occurred in function",                    \
       __PRETTY_FUNCTION__,                     \
       "of file", __FILE__,                     \
       "on line", __LINE__);                    \
  exit(1);                                                              \
}

#define S2LET_ERROR_MEM_ALLOC_CHECK(pointer)              \
  if(pointer == NULL) {                         \
    S2LET_ERROR_GENERIC("Memory allocation failed")           \
  }

#ifdef __cplusplus
}
#endif
#endif
