cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(coal_fluctu CXX C)

add_executable(coal_fluctu coal_fluctu.cpp)
target_compile_features(coal_fluctu PUBLIC cxx_std_11)

find_package(libcloudph++ REQUIRED)
#find_package(synth_turb REQUIRED)
find_package(OpenMP REQUIRED)

target_link_libraries(coal_fluctu PRIVATE clphxx::cloudphxx_lgrngn)
#target_link_libraries(coal_fluctu PRIVATE synth_turb)
target_link_libraries(coal_fluctu PRIVATE OpenMP::OpenMP_CXX)
