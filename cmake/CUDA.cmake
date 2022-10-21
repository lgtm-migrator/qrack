option (ENABLE_CUDA "Build CUDA-based QEngine type" OFF)

find_package (CUDA)
if (NOT CUDA_FOUND OR ENABLE_OPENCL)
    set(ENABLE_CUDA OFF)
endif ()

message ("CUDA Support is: ${ENABLE_CUDA}")

if (ENABLE_CUDA)
    enable_language(CUDA)
    target_compile_definitions(qrack PUBLIC ENABLE_CUDA=1)
    target_include_directories (qrack PUBLIC ${PROJECT_BINARY_DIR} ${CUDA_INCLUDE_DIRS})
    target_compile_options (qrack PUBLIC ${CUDA_COMPILATION_OPTIONS})
    
    target_link_libraries (qrack ${CUDA_LIBRARIES})
    set_target_properties(qrack PROPERTIES CUDA_ARCHITECTURES native)
    target_link_libraries (qrack_pinvoke ${CUDA_LIBRARIES})
    set_target_properties(qrack_pinvoke PROPERTIES CUDA_ARCHITECTURES native)
    if (NOT ENABLE_EMIT_LLVM)
        target_link_libraries (unittest ${CUDA_LIBRARIES})
        set_target_properties(unittest PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (benchmarks ${CUDA_LIBRARIES})
        set_target_properties(benchmarks PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (qrack_cl_precompile ${CUDA_LIBRARIES})
        set_target_properties(qrack_cl_precompile PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (quantum_associative_memory ${CUDA_LIBRARIES})
        set_target_properties(quantum_associative_memory PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (teleport ${CUDA_LIBRARIES})
        set_target_properties(teleport PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (qneuron_classification ${CUDA_LIBRARIES})
        set_target_properties(qneuron_classification PROPERTIES CUDA_ARCHITECTURES native)
        target_link_libraries (cosmology ${CUDA_LIBRARIES})
        set_target_properties(cosmology PROPERTIES CUDA_ARCHITECTURES native)
        if (ENABLE_ALU)
            target_link_libraries (grovers ${CUDA_LIBRARIES})
            set_target_properties(grovers PROPERTIES CUDA_ARCHITECTURES native)
            target_link_libraries (grovers_lookup ${CUDA_LIBRARIES})
            set_target_properties(grovers_lookup PROPERTIES CUDA_ARCHITECTURES native)
            target_link_libraries (ordered_list_search ${CUDA_LIBRARIES})
            set_target_properties(ordered_list_search PROPERTIES CUDA_ARCHITECTURES native)
            target_link_libraries (shors_factoring ${CUDA_LIBRARIES})
            set_target_properties(shors_factoring PROPERTIES CUDA_ARCHITECTURES native)
            target_link_libraries (pearson32 ${CUDA_LIBRARIES})
            set_target_properties(pearson32 PROPERTIES CUDA_ARCHITECTURES native)
            target_link_libraries (quantum_perceptron ${CUDA_LIBRARIES})
            set_target_properties(quantum_perceptron PROPERTIES CUDA_ARCHITECTURES native)
        endif (ENABLE_ALU)
    endif (NOT ENABLE_EMIT_LLVM)
    
    # Add the CUDA objects to the library
    target_sources (qrack PRIVATE
        src/common/cudaengine.cu
        src/common/qengine.cu
        src/qengine/cuda.cu
        src/qhybrid.cpp
        src/qunitmulti.cpp
        )
endif(ENABLE_CUDA)
