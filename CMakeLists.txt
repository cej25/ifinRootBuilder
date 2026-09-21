cmake_minimum_required(VERSION 3.16)
project(raw_root_analysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree TreePlayer Hist Graf Gpad Imt)
find_package(Threads REQUIRED)

set(HISTOGRAM_SOURCES
    AngularCoincidenceHistograms.cpp
    DetectorMatchingStatistics.cpp
    DetectorHistograms.cpp
    FoldStatistics.cpp
    GammaCoincidenceHistograms.cpp
    GermaniumConditionHistograms.cpp
    GermaniumDetectorHistograms.cpp
    GateStatistics.cpp
    IndividualDetectorHistograms.cpp
    MultiplicityStatistics.cpp
    ProgressReporter.cpp
    RunningTimeMap.cpp
)

add_executable(analyse_raw
    analyse_raw.cpp
    RawAnalysis.cpp
    GermaniumCalibration.cpp
    RunCalibrationManager.cpp
    ${HISTOGRAM_SOURCES}
)

add_executable(build_analysis_tree
    build_analysis_tree.cpp
    AnalysisTreeBuilder.cpp
    GermaniumCalibration.cpp
    RunCalibrationManager.cpp
)

add_executable(analyse_tree
    analyse_tree.cpp
    AnalysisTreeAnalysis.cpp
    ${HISTOGRAM_SOURCES}
)

set(ROOT_ANALYSIS_LIBRARIES
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Imt
    Threads::Threads
)

foreach(target analyse_raw build_analysis_tree analyse_tree)
    target_link_libraries(${target} PRIVATE ${ROOT_ANALYSIS_LIBRARIES})
    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
    )
endforeach()
