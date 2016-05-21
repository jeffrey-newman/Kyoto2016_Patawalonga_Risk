PROJECT(PatawalongaRisk4Kyoto2016)
cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

# find_package(CUDA REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem program_options system graph)
find_package(GDAL)

set(AHHZ_RASTER_UTIL "Not Specified" CACHE PATH "Path to Alex Hagen-Zanker's Raster Utility Include direction (select above blink)")
set(AHHZ_RASTER_ITERATOR "Not Specified" CACHE PATH "Path to Alex Hagen-Zanker's Raster Iterator Include direction (select folder above blink)")
set(DAMAGE_SCANNER_SRC "Not Specified" CACHE PATH "Path to Damage Scanner sources")
set(CCCALM_SRC "Not Specified" CACHE PATH "Path to CCCALM sources")
set(INUNDATE_SRC "Not Specified" CACHE PATH "Path to inundate landscape sources")

include_directories(${GDAL_INCLUDE_DIR} )
include_directories(${Boost_INCLUDE_DIRS} )
include_directories(${AHHZ_RASTER_UTIL} )
include_directories(${AHHZ_RASTER_ITERATOR} )
include_directories(${DAMAGE_SCANNER_SRC} )
include_directories(${CCCALM_SRC} )
include_directories(${INUNDATE_SRC} )
include_directories("${INUNDATE_SRC}/eigen")


SET (Sources
		main.cpp
	)

SET (Headers

	)

file(GLOB_RECURSE ExternalDependancies
 "${AHHZ_RASTER_UTIL}/*.h"
 "${AHHZ_RASTER_UTIL}/*.hpp"
 "${AHHZ_RASTER_UTIL}/*.cpp" 
 "${AHHZ_RASTER_ITERATOR}/*.h" 
 "${AHHZ_RASTER_ITERATOR}/*.hpp" 
 "${AHHZ_RASTER_ITERATOR}/*.cpp"
 )
 
 SET (OtherModules 
 
 "${DAMAGE_SCANNER_SRC}/CalcDamage.cpp" 
 "${DAMAGE_SCANNER_SRC}/CalcDamage.hpp" 
 "${DAMAGE_SCANNER_SRC}/Types.hpp" 
 
 "${CCCALM_SRC}/CCCALM.cpp" 
 "${CCCALM_SRC}/Data.cpp" 
 "${CCCALM_SRC}/Initialise.cpp"
 "${CCCALM_SRC}/Simulate.cpp"
 "${CCCALM_SRC}/Finalise.cpp"
 "${CCCALM_SRC}/CCCALM.hpp"
 "${CCCALM_SRC}/Data.hpp"
 "${CCCALM_SRC}/Simulate.hpp"
 "${CCCALM_SRC}/Initialise.hpp"
 
 "${INUNDATE_SRC}/InundateLandscape.cpp"
 "${INUNDATE_SRC}/Neighbourhood.cpp"
 "${INUNDATE_SRC}/Neighbourhood_Impl.cpp"
 "${INUNDATE_SRC}/GraphSearches.h"
 "${INUNDATE_SRC}/IsChannel.h"
 "${INUNDATE_SRC}/PrintGraphsToFile.h"
 "${INUNDATE_SRC}/RasterCoordinates.h"
 "${INUNDATE_SRC}/ReadGraphsFromFile.h"
 "${INUNDATE_SRC}/ReadInControlsAndGuages.h"
 "${INUNDATE_SRC}/remove_processed.h"
 "${INUNDATE_SRC}/Types.h"
 "${INUNDATE_SRC}/Neighbourhood.h"
 "${INUNDATE_SRC}/GraphSearches.cpp"
 "${INUNDATE_SRC}/IsChannel.cpp"
 "${INUNDATE_SRC}/PrintGraphsToFile.cpp"
 "${INUNDATE_SRC}/RasterCoordinates.cpp"
 )





#LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)

#cuda_add_executable(cccalcm ${Sources} ${Headers})
	
#set_property(TARGET ccalcm PROPERTY CXX_STANDARD 11)
add_executable(pataRisk ${Sources} ${Headers} ${ExternalDependancies} ${OtherModules})
target_link_libraries(pataRisk ${GDAL_LIBRARY} ${Boost_LIBRARIES})
