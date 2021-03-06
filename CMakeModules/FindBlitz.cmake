# BLITZ_INCLUDE_DIR = blitz.h
# BLITZ_LIBRARIES = libblitz.so
# BLITZ_FOUND = true if BLITZ is found


IF(BLITZ_PREFIX)
  set(BLITZ_INCLUDE_DIRS "${BLITZ_PREFIX}/include/" ${BLITZ_INCLUDE_DIRS})
  set(BLITZ_LIBRARY_DIRS "${BLITZ_PREFIX}/lib64/" ${BLITZ_LIBRARY_DIRS})
ENDIF(BLITZ_PREFIX)

IF(BLITZ_PREFIX)
  FIND_PATH(BLITZ_INCLUDE_DIR NAMES blitz/blitz.h PATHS ${BLITZ_INCLUDE_DIRS} NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH)
  FIND_LIBRARY(BLITZ_LIBRARIES blitz ${BLITZ_LIBRARY_DIRS} NO_DEFAULT_PATH)
ELSE(BLITZ_PREFX)
  FIND_PATH(BLITZ_INCLUDE_DIR NAMES blitz/blitz.h PATHS ${BLITZ_INCLUDE_DIRS})
  FIND_LIBRARY(BLITZ_LIBRARIES blitz ${BLITZ_LIBRARY_DIRS})	
ENDIF(BLITZ_PREFIX)

SET(BLITZ_FOUND 0)
IF(BLITZ_INCLUDE_DIR AND BLITZ_LIBRARIES)
  SET(BLITZ_FOUND 1)
ENDIF(BLITZ_INCLUDE_DIR AND BLITZ_LIBRARIES)

IF(CMAKE_COMPILER_IS_GNUCXX)
  FIND_PATH(BLITZ_PLATFORM_INCLUDE_DIR blitz/gnu/bzconfig.h ${BLITZ_INCLUDE_DIR})
ENDIF(CMAKE_COMPILER_IS_GNUCXX)
IF(INTEL_COMPILER)
  FIND_PATH(BLITZ_PLATFORM_INCLUDE_DIR blitz/intel/bzconfig.h ${BLITZ_INCLUDE_DIR})
ENDIF(INTEL_COMPILER)

IF(NOT BLITZ_PLATFORM_INCLUDE_DIR)
  FIND_PATH(BLITZ_LIBRARY_DIR libblitz.so ${BLITZ_LIBRARY_DIRS})
  IF(BLITZ_LIBRARY_DIR)
    IF(CMAKE_COMPILER_IS_GNUCXX)
      FIND_PATH(BLITZ_PLATFORM_INCLUDE_DIR blitz/gnu/bzconfig.h "${BLITZ_LIBRARY_DIR}/blitz/include")
    ENDIF(CMAKE_COMPILER_IS_GNUCXX)
    IF(INTEL_COMPILER)
      FIND_PATH(BLITZ_PLATFORM_INCLUDE_DIR blitz/intel/bzconfig.h "${BLITZ_LIBRARY_DIR}/blitz/include")
    ENDIF(INTEL_COMPILER)
  ENDIF(BLITZ_LIBRARY_DIR)
ENDIF(NOT BLITZ_PLATFORM_INCLUDE_DIR)

SET(BLITZ_FOUND 0)
IF(BLITZ_INCLUDE_DIR AND BLITZ_LIBRARIES AND BLITZ_PLATFORM_INCLUDE_DIR)
  SET(BLITZ_FOUND 1)
ENDIF(BLITZ_INCLUDE_DIR AND BLITZ_LIBRARIES AND BLITZ_PLATFORM_INCLUDE_DIR)


#message("BLITZ_FOUND: " ${BLITZ_FOUND})
#message("BLITZ_INCLUDE_DIR: " ${BLITZ_INCLUDE_DIR})
#message("BLITZ_PLATFORM_INCLUDE_DIR: " ${BLITZ_PLATFORM_INCLUDE_DIR})
#message("BLITZ_LIBRARIES: " ${BLITZ_LIBRARIES})


MARK_AS_ADVANCED(
   BLITZ_INCLUDE_DIR
   BLITZ_LIBRARIES
   BLITZ_FOUND
)

IF(NOT ${BLITZ_INCLUDE_DIR} STREQUAL ${BLITZ_PLATFORM_INCLUDE_DIR})
  MARK_AS_ADVANCED(
    BLITZ_PLATFORM_INCLUDE_DIR
  )
ENDIF(NOT ${BLITZ_INCLUDE_DIR} STREQUAL ${BLITZ_PLATFORM_INCLUDE_DIR})

