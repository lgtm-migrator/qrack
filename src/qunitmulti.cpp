//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2021. All rights reserved.
//
// QUnitMulti is a multiprocessor variant of QUnit.
// QUnit maintains explicit separability of qubits as an optimization on a QEngine.
// See https://arxiv.org/abs/1710.05867
// (The makers of Qrack have no affiliation with the authors of that paper.)
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "qfactory.hpp"

namespace Qrack {

QUnitMulti::QUnitMulti(std::vector<QInterfaceEngine> eng, bitLenInt qBitCount, bitCapInt initState,
    qrack_rand_gen_ptr rgp, complex phaseFac, bool doNorm, bool randomGlobalPhase, bool useHostMem, int64_t deviceID,
    bool useHardwareRNG, bool useSparseStateVec, real1_f norm_thresh, std::vector<int64_t> devList,
    bitLenInt qubitThreshold, real1_f sep_thresh)
    : QUnit(eng, qBitCount, initState, rgp, phaseFac, doNorm, randomGlobalPhase, useHostMem, -1, useHardwareRNG,
          useSparseStateVec, norm_thresh, devList, qubitThreshold, sep_thresh)
    , isQEngineOCL(false)
{
#if ENABLE_ENV_VARS
    isRedistributing = (bool)getenv("QRACK_ENABLE_QUNITMULTI_REDISTRIBUTE");
#else
    isRedistributing = false;
#endif

    for (size_t i = 0U; i < engines.size(); i++) {
        if ((engines[i] == QINTERFACE_CPU) || (engines[i] == QINTERFACE_HYBRID)) {
            break;
        }
        if (engines[i] == QINTERFACE_OPENCL) {
            isQEngineOCL = true;
            break;
        }
    }
    if (engines.back() == QINTERFACE_QPAGER) {
        isQEngineOCL = true;
    }

    if (qubitThreshold) {
        thresholdQubits = qubitThreshold;
    } else {
        const bitLenInt gpuQubits =
            log2(OCLEngine::Instance().GetDeviceContextPtr(devID)->GetPreferredConcurrency()) + 1U;
        const bitLenInt cpuQubits = (GetStride() <= ONE_BCI) ? 0U : (log2(GetStride() - ONE_BCI) + 1U);
        thresholdQubits = gpuQubits < cpuQubits ? gpuQubits : cpuQubits;
    }

    std::vector<DeviceContextPtr> deviceContext = OCLEngine::Instance().GetDeviceContextPtrVector();
    defaultDeviceID = (deviceID < 0) ? OCLEngine::Instance().GetDefaultDeviceID() : (size_t)deviceID;

    const size_t devCount = devList.size() ? devList.size() : deviceContext.size();
    for (size_t i = 0; i < devCount; ++i) {
        DeviceInfo deviceInfo;
        deviceInfo.id =
            devList.size() ? ((devList[0U] < 0) ? OCLEngine::Instance().GetDefaultDeviceID() : (size_t)devList[i]) : i;
        deviceList.push_back(deviceInfo);
    }
    if (!devList.size()) {
        std::swap(deviceList[0U], deviceList[defaultDeviceID]);
    }

    for (size_t i = 0U; i < deviceList.size(); ++i) {
        deviceList[i].maxSize = deviceContext[deviceList[i].id]->GetMaxAlloc();
    }

    if (!devList.size()) {
        std::sort(deviceList.begin() + 1U, deviceList.end(), std::greater<DeviceInfo>());
    }
}

QInterfacePtr QUnitMulti::MakeEngine(bitLenInt length, bitCapInt perm)
{
    size_t deviceId = defaultDeviceID;
    uint64_t sz = OCLEngine::Instance().GetActiveAllocSize(deviceId);

    for (size_t i = 0U; i < deviceList.size(); ++i) {
        uint64_t tSz = OCLEngine::Instance().GetActiveAllocSize(deviceList[i].id);
        if (sz > tSz) {
            sz = tSz;
            deviceId = deviceList[i].id;
        }
    }

    // Suppress passing device list, since QUnitMulti occupies all devices in the list
    QInterfacePtr toRet = CreateQuantumInterface(engines, length, perm, rand_generator, phaseFactor, doNormalize,
        randGlobalPhase, useHostRam, deviceId, useRDRAND, isSparse, (real1_f)amplitudeFloor, std::vector<int64_t>{},
        thresholdQubits, separabilityThreshold);

    return toRet;
}

std::vector<QEngineInfo> QUnitMulti::GetQInfos()
{
    // Get shard sizes and devices
    std::vector<QInterfacePtr> qips;
    std::vector<QEngineInfo> qinfos;

    for (auto&& shard : shards) {
        if (shard.unit && (std::find(qips.begin(), qips.end(), shard.unit) == qips.end())) {
            qips.push_back(shard.unit);
            const size_t deviceIndex = std::distance(
                deviceList.begin(), std::find_if(deviceList.begin(), deviceList.end(), [&](DeviceInfo di) {
                    return di.id == (shard.unit->GetDevice() < 0) ? OCLEngine::Instance().GetDefaultDeviceID()
                                                                  : (size_t)shard.unit->GetDevice();
                }));
            qinfos.push_back(QEngineInfo(shard.unit, deviceIndex));
        }
    }

    // We distribute in descending size order:
    std::sort(qinfos.rbegin(), qinfos.rend());

    return qinfos;
}

void QUnitMulti::RedistributeQEngines()
{
    // Only redistribute if the env var flag is set and NOT a null string.
    // No need to redistribute, if there is only 1 device
    if (deviceList.size() <= 1U) {
        return;
    }

    // Get shard sizes and devices
    std::vector<QEngineInfo> qinfos = GetQInfos();

    std::vector<bitCapInt> devSizes(deviceList.size(), 0U);

    for (size_t i = 0U; i < qinfos.size(); ++i) {
        // We want to proactively set OpenCL devices for the event they cross threshold.
        if (!isRedistributing &&
            !((qinfos[i].unit->GetMaxQPower() <= 2U) ||
                (!isQEngineOCL && (qinfos[i].unit->GetQubitCount() < thresholdQubits)) ||
                qinfos[i].unit->isClifford())) {
            continue;
        }

        // If the original OpenCL device has equal load to the least, we prefer the original.
        int64_t deviceID = qinfos[i].unit->GetDevice();
        int64_t devIndex = qinfos[i].deviceIndex;
        bitCapInt sz = devSizes[devIndex];

        // If the original device has 0 determined load, don't switch the unit.
        if (sz) {
            // If the default OpenCL device has equal load to the least, we prefer the default.
            if (devSizes[0U] < sz) {
                deviceID = deviceList[0U].id;
                devIndex = 0U;
                sz = devSizes[0U];
            }

            // Find the device with the lowest load.
            for (size_t j = 0U; j < deviceList.size(); ++j) {
                if ((devSizes[j] < sz) && ((devSizes[j] + qinfos[i].unit->GetMaxQPower()) <= deviceList[j].maxSize)) {
                    deviceID = deviceList[j].id;
                    devIndex = j;
                    sz = devSizes[j];
                }
            }

            // Add this unit to the device with the lowest load.
            qinfos[i].unit->SetDevice(deviceID);
        }

        // Update the size of buffers handles by this device.
        devSizes[devIndex] += qinfos[i].unit->GetMaxQPower();
    }
}
} // namespace Qrack
