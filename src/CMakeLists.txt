include_directories(${PROJECT_SOURCE_DIR}/inc)

# Enable OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(FieldObj FieldStructure.cpp)
target_link_libraries(FieldObj       PUBLIC 
                                            OpenMP::OpenMP_CXX
                                            $<$<CXX_COMPILER_ID:GNU>:m>)
target_include_directories(FieldObj PUBLIC ${INCLUDE_HEADERS})


add_library(icsObj initial_distributions.cpp)
add_library(PanelObj Panel.cpp)
# target_include_directories(icsObj PUBLIC ${INCLUDE_HEADERS} ${Boost_INCLUDE})

add_library(AMRObj_cpu  AMRStructure.cpp 
                    AMRStructure_io.cpp
                    AMRStructure_amr.cpp
                    # AMRStructure_field.cpp
                    # AMRStructure_step.cpp
                    # AMRStructure_interp.cpp
                    )

target_link_libraries(AMRObj_cpu    PUBLIC  PanelObj
                                            FieldObj
                                            icsObj
                                            OpenMP::OpenMP_CXX
                                            $<$<CXX_COMPILER_ID:NVIDIA>:m>)
target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
                                                ${INCLUDE_HEADERS})


# Set additional flags for NVIDIA compilers if ENABLE_GPU_BUILD is on
if(ENABLE_GPU_BUILD AND CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -acc -Minfo=accel -gpu=managed") # -gpu=managed -ta=tesla:managed -gpu=nordc
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -acc")
else()
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "Build OpenACC GPU version" FORCE)
endif()


if(ENABLE_GPU_BUILD)
    add_library(FieldObj_gpu FieldStructure.cpp)
    target_compile_definitions(FieldObj_gpu PRIVATE OPENACC_ENABLED)
    target_link_libraries(FieldObj_gpu  PUBLIC  
                                                OpenMP::OpenMP_CXX
                                                OpenACC_CXX
                                                $<$<CXX_COMPILER_ID:NVIDIA>:m>) #GNU
    # target_compile_options(FieldObj_gpu PRIVATE -O0)
    target_include_directories(FieldObj_gpu PUBLIC ${INCLUDE_HEADERS})

#    set(TRGT AMRObj_gpu)
#    add_library(${TRGT}  AMRStructure.cpp 
#                        AMRStructure_io.cpp
#                        AMRStructure_amr.cpp
#                        AMRStructure_field.cpp
#                        AMRStructure_step.cpp
#                        AMRStructure_interp.cpp)
#    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
#    target_link_libraries(AMRObj_gpu    PUBLIC
#                                        OpenMP::OpenMP_CXX
#                                        MPI::MPI_CXX
#                                        OpenACC_CXX
#                                        $<$<CXX_COMPILER_ID:GNU>:m>)
#    target_include_directories(AMRObj_cpu PUBLIC    ${PROJECT_SOURCE_DIR}/inc_eigen/eigen-3.3.7
#                                                    ${INCLUDE_HEADERS})
endif()