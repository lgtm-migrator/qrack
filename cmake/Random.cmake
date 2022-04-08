include(CheckIncludeFileCXX)
option (SEED_DEVRAND "Seeding PRNG from /dev/random instead of fallback to system clock" ON)
if (SEED_DEVRAND)
    CHECK_INCLUDE_FILE_CXX("sys/random.h" SEED_DEVRAND)
    if (APPLE OR MSVC)
        set(SEED_DEVRAND OFF)
    endif (APPLE OR MSVC)
    if (NOT SEED_DEVRAND)
        message ("/dev/random seed is NOT available on this system.")
    endif (NOT SEED_DEVRAND)
endif (SEED_DEVRAND)
message ("Seeding PRNG from /dev/random (instead of fallback to system clock): ${SEED_DEVRAND}")

option (ENABLE_RDRAND "Use RDRAND hardware random number generation, if available" ON)
if (ENABLE_RDRAND)
    set(QRACK_COMPILE_OPTS ${QRACK_COMPILE_OPTS} -mrdrnd)
    target_compile_definitions(qrack PUBLIC ENABLE_RDRAND=1)
endif (ENABLE_RDRAND)
message ("Try RDRAND is: ${ENABLE_RDRAND}")

option (ENABLE_RNDFILE "Get random numbers from ~/.qrack/rng directory" OFF)
if (ENABLE_RNDFILE)
    target_compile_definitions(qrack PUBLIC ENABLE_RNDFILE=1)
endif (ENABLE_RNDFILE)
message ("RNDFILE is: ${ENABLE_RNDFILE}")

option (ENABLE_DEVRAND "Get random numbers from /dev/urandom" OFF)
if (ENABLE_DEVRAND)
    CHECK_INCLUDE_FILE_CXX("sys/random.h" ENABLE_DEVRAND)
    if (NOT ENABLE_DEVRAND)
        message(FATAL_ERROR "Requested RNG from /dev/urandom, but header is not available on this system!")
    endif (NOT ENABLE_DEVRAND)

    target_compile_definitions(qrack PUBLIC ENABLE_DEVRAND=1)
endif (ENABLE_DEVRAND)
message ("All RNG direct from /dev/urandom is: ${ENABLE_DEVRAND}")

if ((NOT ENABLE_DEVRAND) AND ENABLE_RNDFILE)
    target_sources (qrack PRIVATE
        src/common/rdrandwrapper.cpp
        )
endif ((NOT ENABLE_DEVRAND) AND ENABLE_RNDFILE)
