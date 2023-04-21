cmake_minimum_required(VERSION 3.18.4)
find_package(PkgConfig REQUIRED)

# Check if a C compiler is explicitly stated
if (NOT DEFINED CMAKE_C_COMPILER)
    # Set the compilers based upon the PETSc package
    pkg_get_variable(PETSC_C_COMPILER PETSc ccompiler)
    set(CMAKE_C_COMPILER ${PETSC_C_COMPILER})

    pkg_get_variable(PETSC_C_FLAGS PETSc cflags_extra)
    string (REPLACE ";" " " PETSC_C_FLAGS_STR "${PETSC_C_FLAGS}")
    set(CMAKE_C_FLAGS_INIT ${PETSC_C_FLAGS_STR})
endif ()

# Check if a CXX compiler is explicitly stated
if (NOT DEFINED CMAKE_CXX_COMPILER)
    # Set the compilers based upon the PETSc package
    pkg_get_variable(PETSC_CXX_COMPILER PETSc cxxcompiler)
    set(CMAKE_CXX_COMPILER ${PETSC_CXX_COMPILER})

    pkg_get_variable(PETSC_CXX_FLAGS PETSc cxxflags_extra)
    string (REPLACE ";" " " PETSC_CXX_FLAGS_STR "${PETSC_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_INIT ${PETSC_CXX_FLAGS_STR})
endif ()

# Set the project details
project(petscTest VERSION 0.8.30)

# Load the Required 3rd Party Libaries
pkg_check_modules(PETSc REQUIRED IMPORTED_TARGET GLOBAL PETSc)

add_executable(example example.c)
target_link_libraries(example PUBLIC         PkgConfig::PETSc)
#
#add_executable(particlesInDm particlesInDm.c)
#target_link_libraries(particlesInDm PUBLIC         PkgConfig::PETSc)
#
#add_executable(dmLocate dmLocate.c)
#target_link_libraries(dmLocate PUBLIC         PkgConfig::PETSc)
#
#add_executable(petscSf petscSf.c)
#target_link_libraries(petscSf PUBLIC         PkgConfig::PETSc)
#
#add_executable(petscSfPlayground petscSfPlayground.cpp)
#target_link_libraries(petscSfPlayground PUBLIC         PkgConfig::PETSc)
#
#add_executable(gmshVolumeCheck gmshVolumeCheck.c)
#target_link_libraries(gmshVolumeCheck PUBLIC         PkgConfig::PETSc)
#
#add_executable(extchemfield extchemfield.c)
#target_link_libraries(extchemfield PUBLIC         PkgConfig::PETSc)
#
#add_executable(extchem extchem.c)
#target_link_libraries(extchem PUBLIC         PkgConfig::PETSc)
#


#
#add_executable(gMshDistribute gmshDistribute.c)
#target_link_libraries(gMshDistribute PUBLIC         PkgConfig::PETSc)
#
#add_executable(gMshExtrude gmshExtrude.c)
#target_link_libraries(gMshExtrude PUBLIC         PkgConfig::PETSc)
#
#add_executable(gmshIO gmshIO.c)
#target_link_libraries(gmshIO PUBLIC         PkgConfig::PETSc)
#
#add_executable(subDmTest subDmTest.c)
#target_link_libraries(subDmTest PUBLIC         PkgConfig::PETSc)
#
#add_executable(gmshAreaCheck gmshAreaCheck.c)
#target_link_libraries(gmshAreaCheck PUBLIC         PkgConfig::PETSc)
#
#add_executable(periodicMeshOutput periodicMeshOutput.c)
#target_link_libraries(periodicMeshOutput PUBLIC         PkgConfig::PETSc)
#
#
#
## XDMF
#add_executable(xdmfGen xdmfGen.c)
#target_link_libraries(xdmfGen PUBLIC         PkgConfig::PETSc)

# Include system specific hacks
#if ("${APPLE}" AND (${CMAKE_C_COMPILER_ID} STREQUAL "AppleClang"))
    # check for homebrew gfortran and get path for libstdc++.dylib
#    execute_process(COMMAND gfortran --print-file-name=libstdc++.dylib OUTPUT_VARIABLE LIBSTDCPP_PATH)
#
#    # convert to an absolute path and get the directory
#    get_filename_component(LIBSTDCPP_PATH ${LIBSTDCPP_PATH} ABSOLUTE)
#    get_filename_component(LIBSTDCPP_PATH ${LIBSTDCPP_PATH} DIRECTORY)
#
#    target_link_directories(particlesInDm PUBLIC ${LIBSTDCPP_PATH})
#    target_link_directories(example PUBLIC ${LIBSTDCPP_PATH})
#    target_link_directories(dmLocate PUBLIC ${LIBSTDCPP_PATH})
#    target_link_directories(petscSf PUBLIC ${LIBSTDCPP_PATH})
#    target_link_directories(gmshVolumeCheck PUBLIC ${LIBSTDCPP_PATH})
#target_link_directories(petscSfPlayground PUBLIC ${LIBSTDCPP_PATH})
#target_link_directories(extchemfield PUBLIC ${LIBSTDCPP_PATH})
#target_link_directories(extchem PUBLIC ${LIBSTDCPP_PATH})

#FIND_LIBRARY(TCHEM_LIBRARY
#        NAMES
#        tchem
#        HINTS
#        ${CMAKE_FIND_ROOT_PATH}
#        ${PETSc_LIBRARY_DIRS}
#        PATHS
#        ${CMAKE_FIND_ROOT_PATH}
#        ${PETSc_LIBRARY_DIRS}
#        )
#if(NOT TCHEM_LIBRARY)
#    message(FATAL_ERROR "Cannot find TChem library.  Please reconfigure PETSc with --download-tchem flag." )
#else()
#    add_library(TChem::TChem UNKNOWN IMPORTED)
#    set_target_properties(
#            TChem::TChem
#            PROPERTIES
#            IMPORTED_LOCATION ${TCHEM_LIBRARY})
#endif()
#target_link_libraries(extchemfield PUBLIC TChem::TChem)
#target_link_libraries(extchem PUBLIC TChem::TChem)

#endif ()
