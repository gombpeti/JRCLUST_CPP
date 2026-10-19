# 3.24 for CMAKE_CUDA_ARCHITECTURES=native.
cmake_minimum_required(VERSION 3.24)

# Must be set BEFORE project(): enabling the CUDA language fills this in with a
# conservative default (sm_52), and a later assignment would be ignored.
# "native" compiles for exactly the GPU(s) in the build machine, e.g. sm_120 on
# an RTX 5070 Ti (which needs CUDA >= 12.8) or sm_89 on an RTX 4060 Ti.
# Override with -DCMAKE_CUDA_ARCHITECTURES=... to build for another machine
# (86 for a 3090, 75 for a 2080, "89;120" for a fat binary covering both).
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
endif()

project(jrclust
    VERSION 0.1.0
    DESCRIPTION "JRCLUST - native C++/CUDA spike sorter"
    LANGUAGES CXX CUDA)

# ------------------------------------------------------------------ options --
option(JRC_BUILD_GUI "Build the Qt curation GUI (Stage 2)" OFF)

# ---------------------------------------------------------------- standards --
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# One output dir keeps the .exe next to the Qt DLLs windeployqt will drop later.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(CUDAToolkit REQUIRED)

# ------------------------------------------------------------------- core ----
# NOTE: src/cuda/kernels/*.cu are deliberately NOT listed here. They are the
# unmodified JRCLUST kernel sources and are #included by DensityPeak.cu, so
# compiling them as separate translation units would duplicate the symbols.
add_library(jrc_core STATIC
    src/core/MemoryMappedFile.cpp
    src/core/Config.cpp
    src/core/Recording.cpp
    src/core/Probe.cpp
    src/detect/Filter.cpp
    src/detect/PeakDetector.cpp
    src/detect/SpikeVeto.cpp
    src/detect/DetectController.cpp
    src/features/FeatureExtractor.cpp
    src/sort/Clustering.cpp
    src/sort/PostProcess.cpp
    src/sort/UnitClassifier.cpp
    src/io/ResultWriter.cpp
    src/cuda/DensityPeak.cu
)

target_include_directories(jrc_core PUBLIC src)

target_link_libraries(jrc_core PUBLIC CUDA::cudart)

set_target_properties(jrc_core PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    POSITION_INDEPENDENT_CODE ON)

if(MSVC)
    # /utf-8: the .prm files and our comments carry non-ASCII (e.g. micrometre).
    target_compile_options(jrc_core PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:/utf-8 /W3 /permissive->)
    # std::getenv is fine for reading JRC_THREADS; MSVC's _s variants are not
    # portable and this code also builds on Linux.
    target_compile_definitions(jrc_core PUBLIC NOMINMAX WIN32_LEAN_AND_MEAN
                                               _CRT_SECURE_NO_WARNINGS)
else()
    target_compile_options(jrc_core PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra>)
endif()

# --------------------------------------------------------------- executable --
add_executable(jrclust src/main.cpp)
target_link_libraries(jrclust PRIVATE jrc_core)

if(MSVC)
    target_compile_options(jrclust PRIVATE /utf-8 /W3)
endif()

# ------------------------------------------------------------------- report --
message(STATUS "")
message(STATUS "  JRCLUST C++/CUDA  ${PROJECT_VERSION}")
message(STATUS "  build type ....... ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ compiler ..... ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  CUDA toolkit ..... ${CUDAToolkit_VERSION}")
message(STATUS "  CUDA arch ........ ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "  Qt GUI ........... ${JRC_BUILD_GUI} (Stage 2)")
message(STATUS "")

