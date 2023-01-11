set (SOURCES_LIB
    EquationOfState.cpp
    FileIO.cpp
    FreeStream.cpp
    FreestreamMilne.cpp
    FSUtilities.cpp
    HydroValidity.cpp
    InitialConditions.cpp
    LandauMatch.cpp
    Memoryf.cpp
    FSUtilities.cpp
    )

set (SOURCES
    RunWrapper.cpp
    )

add_library(freestream_lib SHARED ${SOURCES_LIB})
add_executable(RunWrapper ${SOURCES})
target_link_libraries (freestream_lib ${LIBS})
target_link_libraries (RunWrapper freestream_lib ${LIBS})

install(TARGETS RunWrapper DESTINATION ${CMAKE_HOME_DIRECTORY})
