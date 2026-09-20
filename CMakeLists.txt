cmake_minimum_required(VERSION 3.16)
project(raw_root_analysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree TreePlayer Hist Graf Gpad Imt)

add_executable(analyse_raw
    analyse_raw.cpp
    AngularCoincidenceHistograms.cpp
    RawAnalysis.cpp
    DetectorMatchingStatistics.cpp
    DetectorHistograms.cpp
    FoldStatistics.cpp
    GammaCoincidenceHistograms.cpp
    GermaniumCalibration.cpp
    GermaniumConditionHistograms.cpp
    GermaniumDetectorHistograms.cpp
    GateStatistics.cpp
    IndividualDetectorHistograms.cpp
    MultiplicityStatistics.cpp
    RunCalibrationManager.cpp
)
target_link_libraries(analyse_raw PRIVATE
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Imt
)

target_compile_options(analyse_raw PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
