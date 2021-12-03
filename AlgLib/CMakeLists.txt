set(CMAKE_CXX_STANDARD 20)

add_library(AlgLib)

target_sources(AlgLib
        PRIVATE
        AlgLib.cpp
        ElGamal.cpp
        PUBLIC
        AlgLib.h
        ElGamal.h
        )

target_include_directories(AlgLib
        PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
        )

add_subdirectory(BigIntLib)

target_link_libraries(AlgLib
        PUBLIC
        BigInt
        )