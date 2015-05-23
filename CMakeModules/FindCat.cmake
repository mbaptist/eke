# CAT_INCLUDE_DIR = cat.h
# CAT_LIBRARIES = cat.so
# CAT_FOUND = true if cat is found

IF(CAT_INCLUDE_DIRS)
  FIND_PATH(CAT_INCLUDE_DIR cat.h  ${CAT_INCLUDE_DIRS})
  FIND_LIBRARY(CAT_LIBRARY cat ${CAT_LIBRARY_DIRS})
ELSE(CAT_INCLUDE_DIRS)
  FIND_PATH(CAT_INCLUDE_DIR cat.h $ENV{HOME}/fcup/projects/cat/install/include/cat)
  FIND_LIBRARY(CAT_LIBRARIES cat $ENV{HOME}/fcup/projects/cat/install/lib)
ENDIF(CAT_INCLUDE_DIRS)

SET(CAT_FOUND FALSE)
IF(CAT_INCLUDE_DIR AND CAT_LIBRARIES)
  MESSAGE(STATUS "CAT_INCLUDE_DIR=${CAT_INCLUDE_DIR}")
  MESSAGE(STATUS "CAT_LIBRARIES=${CAT_LIBRARIES}")
  SET(CAT_FOUND TRUE)
ENDIF()

MARK_AS_ADVANCED(
   CAT_INCLUDE_DIR
   CAT_LIBRARIES
   CAT_FOUND
)

