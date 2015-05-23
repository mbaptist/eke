set(TARGET_1 pb_colloid)
set(TARGET_2 pb_double_plane)
set(TARGET_3 p_plane_capacitor)
set(TARGET_4 p_spherical_capacitor)

set(SRC_DIR "${CMAKE_SOURCE_DIR}/src")

set(USE_BLITZ true)
#set(USE_BLITZ false)

if (USE_BLITZ)
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-DUSE_BLITZ")
  find_package (Blitz)
  if (BLITZ_FOUND)
    include_directories(${BLITZ_INCLUDE_DIR})
    if(BLITZ_PLATFORM_INCLUDE_DIR)
      include_directories(${BLITZ_PLATFORM_INCLUDE_DIR})
    endif(BLITZ_PLATFORM_INCLUDE_DIR)
    set(LIBS ${LIBS} ${BLITZ_LIBRARIES})
  endif (BLITZ_FOUND)
else (USE_BLITZ)
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-DUSE_CAT")
  find_package (Cat)
  if (CAT_FOUND)
    include_directories(${CAT_INCLUDE_DIR})
    set(LIBS ${LIBS} ${CAT_LIBRARIES})
  endif (CAT_FOUND)
endif (USE_BLITZ)
  
include(FindPythonLibs)

find_package (PythonLibs)
if(PYTHONLIBS_FOUND)
  include_directories(${PYTHON_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${PYTHON_LIBRARIES})
endif(PYTHONLIBS_FOUND)

include_directories("${SRC_DIR}")

set(HEADERS

"${SRC_DIR}/coord.hpp"
"${SRC_DIR}/geometry.hpp"
"${SRC_DIR}/grid.hpp"
"${SRC_DIR}/io.hpp"
"${SRC_DIR}/maggs.hpp"
"${SRC_DIR}/random.hpp"
"${SRC_DIR}/types.hpp"

)

set(SOURCES

"${SRC_DIR}/coord.cpp"
"${SRC_DIR}/io.cpp"
"${SRC_DIR}/maggs.cpp"

)

add_executable(${TARGET_1} ${SOURCES} "${SRC_DIR}/pb_colloid.cpp")
target_link_libraries(${TARGET_1} ${LIBS})
add_executable(${TARGET_2} ${SOURCES} "${SRC_DIR}/pb_double_plane.cpp")
target_link_libraries(${TARGET_2} ${LIBS})
add_executable(${TARGET_3} ${SOURCES} "${SRC_DIR}/p_plane_capacitor.cpp")
target_link_libraries(${TARGET_3} ${LIBS})
add_executable(${TARGET_4} ${SOURCES} "${SRC_DIR}/p_spherical_capacitor.cpp")
target_link_libraries(${TARGET_4} ${LIBS})

install (TARGETS ${TARGET_1} ${TARGET_2} ${TARGET_3} ${TARGET_4} DESTINATION bin)

