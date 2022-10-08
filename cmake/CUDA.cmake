option (ENABLE_CUDA "Build CUDA-based QEngine type" OFF)

find_package (CUDA)
if (NOT CUDA_FOUND)
    set(ENABLE_CUDA OFF)
endif ()

message ("CUDA Support is: ${ENABLE_CUDA}")

if (ENABLE_CUDA)
    enable_language(CUDA)
    target_compile_definitions(qrack PUBLIC ENABLE_CUDA=1)
    target_include_directories (qrack PUBLIC ${PROJECT_BINARY_DIR} ${CUDA_INCLUDE_DIRS})
    target_compile_options (qrack PUBLIC ${CUDA_COMPILATION_OPTIONS})
    
    target_link_libraries (qrack_pinvoke ${CUDA_LIBRARIES})
    target_link_libraries (unittest ${CUDA_LIBRARIES})
    target_link_libraries (benchmarks ${CUDA_LIBRARIES})
    target_link_libraries (qrack_cl_precompile ${CUDA_LIBRARIES})
    target_link_libraries (grovers ${CUDA_LIBRARIES})
    target_link_libraries (grovers_lookup ${CUDA_LIBRARIES})
    target_link_libraries (ordered_list_search ${CUDA_LIBRARIES})
    target_link_libraries (quantum_perceptron ${CUDA_LIBRARIES})
    target_link_libraries (quantum_associative_memory ${CUDA_LIBRARIES})
    target_link_libraries (shors_factoring ${CUDA_LIBRARIES})
    target_link_libraries (pearson32 ${CUDA_LIBRARIES})
    target_link_libraries (teleport ${CUDA_LIBRARIES})
    target_link_libraries (qneuron_classification ${CUDA_LIBRARIES})
    
    # Add the CUDA objects to the library
    target_sources (qrack PRIVATE
        src/common/cudaengine.cu
        src/common/qengine.cu
        src/qengine/cuda.cu
        )
endif(ENABLE_CUDA)
