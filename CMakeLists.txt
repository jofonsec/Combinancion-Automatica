######################################################################################
### 0) Definir las propiedades de la aplicación
######################################################################################
# check cmake version compatibility
# Comprobar la version de cmake 
cmake_minimum_required(VERSION 2.6)

# Here define your project name
# Aquí se define el nombre de tu proyecto
project(Automatic)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -pthread")
#SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -pthread")
#SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pthread")

# Here define the name and the version of your package
# Aqui se define el nombre y la version del problema
SET(PACKAGE_NAME "Experimentos" CACHE STRING "Package name" FORCE)
SET(PACKAGE_VERSION "1.0" CACHE STRING "Package version" FORCE)

# Aquí se define la ubicación de ParadisEO #IMPORTANTE#
set(PARADISEO_ROOT /home/jfonseca/GitHub/paradiseo)

# regular expression checking
#include_regular_expression("^.*$" "^$")

# set a language for the entire project.
enable_language(CXX)
enable_language(C)

######################################################################################

######################################################################################
### 2) include the sources
### 2) Incluir los codigos de la biblioteca y modulos
######################################################################################
######################################
### include required modules & utilities
#####################################################################################
include(CMakeBackwardCompatibilityCXX)

include(FindDoxygen)

include(CheckLibraryExists)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/module")
find_package(Paradiseo) #El orden es importante

include_directories(${PARADISEO_INCLUDE_DIR})

######################################################################################

#####################################################################################
### Manage the build type
### Administrar el tipo de construcción
#####################################################################################

# the user should choose the build type on windows environments,excepted under cygwin (default=none)
set(CMAKE_DEFAULT_BUILD_TYPE "Release" CACHE STRING "Variable that stores the default CMake build type" FORCE)

find_program(MEMORYCHECK_COMMAND
    NAMES purify valgrind
    PATHS
    "/usr/local/bin /usr/bin [HKEY_LOCAL_MACHINE\\SOFTWARE\\Rational Software\\Purify\\Setup;InstallFolder]"
    DOC "Path to the memory checking command, used for memory error detection.") 
       

if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DCMAKE_VERBOSE_MAKEFILE=ON)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

add_definitions(-DDEPRECATED_MESSAGES)

######################################################################################
### 3) Link the librairies for your executable
### 3) Encadenar las sub-carpetas para los ejecutables
######################################################################################

add_subdirectory(src)
add_subdirectory(application)
add_subdirectory(datos)

######################################################################################

