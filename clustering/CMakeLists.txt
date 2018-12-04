# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(DAQWork)
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net XMLIO)
find_package(Boost REQUIRED)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIRS})

add_library(Clustering SHARED WireCluster.cc WireHit.cc OpticalCluster.cc OpticalHit.cc ClusterEngine.cc ClusterEnergyEstimator.cc Clustering.cc SmartTrigger.cc SmartTriggerToyThrower.cc BiTrigger.cc)
target_link_libraries(Clustering ${ROOT_LIBRARIES})

add_library(Display SHARED PhotonTimingDisplay.cc PhotonWirePositionDisplay.cc WirePositionTimingDisplay.cc)
target_link_libraries(Display ${ROOT_LIBRARIES})

add_library(SimpleOptHit SHARED HitSimplePlot.cc)
target_link_libraries(SimpleOptHit ${ROOT_LIBRARIES})

add_executable(RunDAQClustering RunDAQClustering.cxx)
target_link_libraries(RunDAQClustering Clustering)

add_executable(AnalyseOpticalCluster AnalyseOpticalCluster.cxx)
target_link_libraries(AnalyseOpticalCluster Clustering)

add_executable(GetEffBackRate GetEffBackRate.cxx)
target_link_libraries(GetEffBackRate Clustering)

add_executable(RunEventDiplay RunEventDisplay.cxx)
target_link_libraries(RunEventDiplay Display)

add_executable(AnalyseDAQClustering AnalyseDAQClustering.cxx)
target_link_libraries(AnalyseDAQClustering Clustering)

add_executable(RunSimpleHitInfo RunSimpleHitInfo.cxx)
target_link_libraries(RunSimpleHitInfo SimpleOptHit)

add_executable(PDSTiming PDSTiming.cxx)
target_link_libraries(PDSTiming ${ROOT_LIBRARIES})

add_executable(NeutronInfo NeutronInfo.cxx)
target_link_libraries(NeutronInfo ${ROOT_LIBRARIES})

add_executable(RunSmartTrigger RunSmartTrigger.cxx)
target_link_libraries(RunSmartTrigger Clustering)

add_executable(ThrowToySmartTrigger ThrowToySmartTrigger.cxx)
target_link_libraries(ThrowToySmartTrigger Clustering)

add_executable(RunBiTrigger RunBiTrigger.cxx)
target_link_libraries(RunBiTrigger Clustering)

add_executable(CreateReweightingScheme CreateReweightingScheme.cxx)
target_link_libraries(CreateReweightingScheme ${ROOT_LIBRARIES})

add_executable(AnalyseWireClustering_SolarNu AnalyseWireClustering_SolarNu.cxx)
target_link_libraries(AnalyseWireClustering_SolarNu ${ROOT_LIBRARIES})

add_executable(EfficiencyCalculator EfficiencyCalculator.cxx)
target_link_libraries(EfficiencyCalculator Clustering)

add_executable(DumpSNRate DumpRate.cxx)
target_link_libraries(DumpSNRate ${ROOT_LIBRARIES})

add_executable(RunNeutronEff RunNeutronEfficiency.cxx)
target_link_libraries(RunNeutronEff Clustering)

add_executable(SingleNeutronEfficiency SingleNeutronEfficiency.cxx)
target_link_libraries(SingleNeutronEfficiency Clustering)

# add_executable(OpHitTrigger OpHitTrigger.C)
# target_link_libraries(OpHitTrigger ${ROOT_LIBRARIES})
