#set(SOURCES
#        SphericalHarmonics.cpp
#        SphericalBessels.cpp
#        Waters.cpp
#        AtomisticModel.cpp
#        functions_WETSAXS.cpp
#        Fit.cpp
#        )

#set(HEADERS
#        SphericalHarmonics.h
#        SphericalBessels.h
#        Waters.h
#        AtomisticModel.h
#        functions_WETSAXS.h
#        Fit.h
#        )

list(APPEND SOURCES
        SphericalHarmonics.cpp
        SphericalBessels.cpp
        Waters.cpp
        AtomisticModel.cpp
        functions_WETSAXS.cpp
        Fit.cpp
        SBESJY.cpp
        )

list(APPEND HEADERS
        SphericalHarmonics.h
        SphericalBessels.h
        Waters.h
        AtomisticModel.h
        functions_WETSAXS.h
        Fit.h
        Result.h
        SBESJY.h
        )

add_library(${CMAKE_PROJECT_NAME}_lib STATIC ${SOURCES} ${HEADERS})

target_link_libraries(${CMAKE_PROJECT_NAME}_lib PUBLIC ${SASTools_LIB})

#target_include_directories(${CMAKE_PROJECT_NAME}_lib PUBLIC include)