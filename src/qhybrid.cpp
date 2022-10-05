//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2021. All rights reserved.
//
// QPager breaks a QEngine instance into pages of contiguous amplitudes.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "qfactory.hpp"

#include <thread>

namespace Qrack {

QHybrid::QHybrid(bitLenInt qBitCount, bitCapInt initState, qrack_rand_gen_ptr rgp, complex phaseFac, bool doNorm,
    bool randomGlobalPhase, bool useHostMem, int64_t deviceId, bool useHardwareRNG, bool useSparseStateVec,
    real1_f norm_thresh, std::vector<int64_t> devList, bitLenInt qubitThreshold, real1_f sep_thresh)
    : QEngine(qBitCount, rgp, doNorm, randomGlobalPhase, useHostMem, useHardwareRNG, norm_thresh)
    , useRDRAND(useHardwareRNG)
    , isSparse(useSparseStateVec)
    , separabilityThreshold(sep_thresh)
    , devID(deviceId)
    , phaseFactor(phaseFac)
    , deviceIDs(devList)
{
    if (qubitThreshold) {
        gpuThresholdQubits = qubitThreshold;
    } else {
        const bitLenInt gpuQubits =
            log2(OCLEngine::Instance().GetDeviceContextPtr(devID)->GetPreferredConcurrency()) + 1U;
        const bitLenInt cpuQubits = (GetStride() <= ONE_BCI) ? 0U : (log2(GetStride() - ONE_BCI) + 1U);
        gpuThresholdQubits = gpuQubits < cpuQubits ? gpuQubits : cpuQubits;
    }

    pagerThresholdQubits = log2(OCLEngine::Instance().GetDeviceContextPtr(devID)->GetMaxAlloc() / sizeof(complex));

#if ENABLE_ENV_VARS
    if (getenv("QRACK_SEGMENT_GLOBAL_QB")) {
        const bitLenInt segmentGlobalQb = (bitLenInt)std::stoi(std::string(getenv("QRACK_SEGMENT_GLOBAL_QB")));
        pagerThresholdQubits = ((pagerThresholdQubits - segmentGlobalQb) > gpuThresholdQubits)
            ? (pagerThresholdQubits - segmentGlobalQb)
            : (gpuThresholdQubits + 1U);
    }
#endif

    isGpu = (qubitCount >= gpuThresholdQubits);
    isPager = (qubitCount > pagerThresholdQubits);

    std::vector<QInterfaceEngine> engines;
    if (isPager) {
        engines.push_back(QINTERFACE_QPAGER);
    }
    engines.push_back(isGpu ? QINTERFACE_OPENCL : QINTERFACE_CPU);

    engine = std::dynamic_pointer_cast<QEngine>(CreateQuantumInterface(engines, qubitCount, initState, rand_generator,
        phaseFactor, doNormalize, randGlobalPhase, useHostRam, devID, useRDRAND, isSparse, (real1_f)amplitudeFloor,
        deviceIDs, pagerThresholdQubits, separabilityThreshold));
}

QEnginePtr QHybrid::MakeEngine(bool isOpenCL)
{
    QEnginePtr toRet =
        std::dynamic_pointer_cast<QEngine>(CreateQuantumInterface(isOpenCL ? QINTERFACE_OPENCL : QINTERFACE_CPU, 0U, 0U,
            rand_generator, phaseFactor, doNormalize, randGlobalPhase, useHostRam, devID, useRDRAND, isSparse,
            (real1_f)amplitudeFloor, deviceIDs, pagerThresholdQubits, separabilityThreshold));
    toRet->SetQubitCount(qubitCount);
    toRet->SetConcurrency(GetConcurrencyLevel());
    return toRet;
}
} // namespace Qrack
