cmake_minimum_required(VERSION 3.12)

project(TREW)

set(CMAKE_CXX_STANDARD 20)

find_package(argparse REQUIRED)
find_package(Boost REQUIRED)
find_package(unordered_dense REQUIRED)
find_package(ZLIB REQUIRED)
find_package(TBB REQUIRED)

add_executable(trew
        src/trew.cpp
        src/kmer.cpp
        src/kmer.h
)

add_executable(trew_test
        test/test.cpp
        src/kmer.cpp
        src/kmer.h
)

target_link_libraries(trew PRIVATE ZLIB::ZLIB TBB::tbb Boost::boost argparse::argparse unordered_dense::unordered_dense)
target_link_libraries(trew_test PRIVATE ZLIB::ZLIB TBB::tbb Boost::boost unordered_dense::unordered_dense)