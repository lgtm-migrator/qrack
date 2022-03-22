//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2021. All rights reserved.
//
// QBinaryDecision tree is an alternative approach to quantum state representation, as
// opposed to state vector representation. This is a compressed form that can be
// operated directly on while compressed. Inspiration for the Qrack implementation was
// taken from JKQ DDSIM, maintained by the Institute for Integrated Circuits at the
// Johannes Kepler University Linz:
//
// https://github.com/iic-jku/ddsim
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#include "qbdt_node.hpp"
#include "qfactory.hpp"

namespace Qrack {

QBdt::QBdt(std::vector<QInterfaceEngine> eng, bitLenInt qBitCount, bitCapInt initState, qrack_rand_gen_ptr rgp,
    complex phaseFac, bool doNorm, bool randomGlobalPhase, bool useHostMem, int deviceId, bool useHardwareRNG,
    bool useSparseStateVec, real1_f norm_thresh, std::vector<int> ignored, bitLenInt qubitThreshold, real1_f sep_thresh)
    : QInterface(qBitCount, rgp, doNorm, useHardwareRNG, randomGlobalPhase, doNorm ? norm_thresh : ZERO_R1)
    , engines(eng)
    , devID(deviceId)
    , root(NULL)
{
#if ENABLE_PTHREAD
    SetConcurrency(std::thread::hardware_concurrency());
#endif
    SetPermutation(initState);
}

QInterfacePtr QBdt::MakeQInterface(bitLenInt qbCount)
{
    return CreateQuantumInterface(engines, qbCount, 0U, rand_generator, ONE_CMPLX, doNormalize, randGlobalPhase, false,
        devID, hardware_rand_generator != NULL, false, amplitudeFloor);
}

QStabilizerPtr QBdt::MakeQStabilizer(bitLenInt qbCount, bitCapInt perm)
{
    return std::dynamic_pointer_cast<QStabilizer>(
        CreateQuantumInterface({ QINTERFACE_STABILIZER }, qbCount, perm, rand_generator, ONE_CMPLX, doNormalize,
            randGlobalPhase, false, devID, hardware_rand_generator != NULL, false, amplitudeFloor));
}

void QBdt::FallbackMtrx(const complex* mtrx, bitLenInt target)
{
    if (!root->branches[0]) {
        root = root->PopSpecial();
    }

    Swap(0U, target);
    Mtrx(mtrx, 0U);
    Swap(0U, target);
}

void QBdt::FallbackMCMtrx(
    const complex* mtrx, const bitLenInt* controls, bitLenInt controlLen, bitLenInt target, bool isAnti)
{
    const bitLenInt gateQubitCount = controlLen + 1U;
    const bitCapInt gatePower = pow2(gateQubitCount);
    for (bitCapInt i = 0; i < gatePower; i++) {
        QBdtNodeInterfacePtr leaf = root;
        QBdtNodeInterfacePtr parent = NULL;
        for (bitLenInt j = 0; j < gateQubitCount; j++) {
            if (IS_NORM_0(leaf->scale)) {
                break;
            }
            if (!leaf->branches[0]) {
                leaf = leaf->PopSpecial();
                if (!j) {
                    root = leaf;
                } else {
                    parent->branches[SelectBit(i, (j - 1U))] = leaf;
                }
            }
            parent = leaf;
            leaf = leaf->branches[SelectBit(i, j)];
        }
    }

    std::unique_ptr<bitLenInt[]> lControls(new bitLenInt[controlLen]);
    for (bitLenInt i = 0U; i < controlLen; i++) {
        lControls[i] = i;
        Swap(i, controls[i]);
    }
    Swap(controlLen, target);

    ApplyControlledSingle(mtrx, lControls.get(), controlLen, controlLen, isAnti);

    Swap(controlLen, target);
    for (bitLenInt i = 0U; i < controlLen; i++) {
        Swap(i, controls[i]);
    }
}

void QBdt::SetPermutation(bitCapInt initState, complex phaseFac)
{
    if (!qubitCount) {
        return;
    }

    if (stateVec) {
        stateVec->SetPermutation(initState);
        return;
    }

    if (phaseFac == CMPLX_DEFAULT_ARG) {
        if (randGlobalPhase) {
            real1_f angle = Rand() * 2 * PI_R1;
            phaseFac = complex((real1)cos(angle), (real1)sin(angle));
        } else {
            phaseFac = ONE_CMPLX;
        }
    }

    root = std::make_shared<QBdtQStabilizerNode>(phaseFac, MakeQStabilizer(qubitCount, initState));
}

QInterfacePtr QBdt::Clone()
{
    QBdtPtr copyPtr = std::make_shared<QBdt>(qubitCount, 0, rand_generator, ONE_CMPLX, doNormalize, randGlobalPhase,
        false, -1, (hardware_rand_generator == NULL) ? false : true, false, (real1_f)amplitudeFloor);

    ResetStateVector();

    copyPtr->root = root ? root->ShallowClone() : NULL;

    return copyPtr;
}

template <typename Fn> void QBdt::GetTraversal(Fn getLambda)
{
    for (bitCapInt i = 0; i < maxQPower; i++) {
        QBdtNodeInterfacePtr leaf = root;
        complex scale = leaf->scale;
        bitLenInt j;
        for (j = 0; j < qubitCount; j++) {
            if (IS_NORM_0(scale) || !leaf->branches[0]) {
                break;
            }
            leaf = leaf->branches[SelectBit(i, j)];
            scale *= leaf->scale;
        }

        if (!IS_NORM_0(scale) && (j < qubitCount)) {
            scale *= NODE_TO_QINTERFACE(leaf)->GetAmplitude(i >> j);
        }

        getLambda((bitCapIntOcl)i, scale);
    }
}
template <typename Fn> void QBdt::SetTraversal(Fn setLambda)
{
    root = std::make_shared<QBdtNode>();

    for (bitCapInt i = 0; i < maxQPower; i++) {
        QBdtNodeInterfacePtr leaf = root;
        QBdtNodeInterfacePtr parent = NULL;
        for (bitLenInt j = 0; j < qubitCount; j++) {
            if (!leaf->branches[0]) {
                leaf = leaf->PopSpecial();
                if (!j) {
                    root = leaf;
                } else {
                    parent->branches[SelectBit(i, (j - 1U))] = leaf;
                }
            } else {
                leaf->Branch();
            }
            parent = leaf;
            leaf = leaf->branches[SelectBit(i, j)];
        }
        setLambda((bitCapIntOcl)i, leaf);
    }

    root->PopStateVector(qubitCount);
    root->Prune(qubitCount);
}
void QBdt::GetQuantumState(complex* state)
{
    GetTraversal([state](bitCapIntOcl i, complex scale) { state[i] = scale; });
}
void QBdt::GetQuantumState(QInterfacePtr eng)
{
    GetTraversal([eng](bitCapIntOcl i, complex scale) { eng->SetAmplitude(i, scale); });
}
void QBdt::SetQuantumState(const complex* state)
{
    if (stateVec) {
        return stateVec->SetQuantumState(state);
    }

    const bitLenInt qbCount = qubitCount;
    SetTraversal([qbCount, state](bitCapIntOcl i, QBdtNodeInterfacePtr leaf) { leaf->scale = state[i]; });
}
void QBdt::SetQuantumState(QInterfacePtr eng)
{
    if (stateVec) {
        stateVec = eng->Clone();
        return;
    }

    const bitLenInt qbCount = qubitCount;
    SetTraversal([qbCount, eng](bitCapIntOcl i, QBdtNodeInterfacePtr leaf) { leaf->scale = eng->GetAmplitude(i); });
}
void QBdt::GetProbs(real1* outputProbs)
{
    GetTraversal([outputProbs](bitCapIntOcl i, complex scale) { outputProbs[i] = norm(scale); });
}

real1_f QBdt::SumSqrDiff(QBdtPtr toCompare)
{
    if (this == toCompare.get()) {
        return ZERO_R1;
    }

    // If the qubit counts are unequal, these can't be approximately equal objects.
    if (qubitCount != toCompare->qubitCount) {
        // Max square difference:
        return ONE_R1;
    }

    ResetStateVector();
    toCompare->ResetStateVector();

    complex projection = ZERO_CMPLX;
    for (bitCapInt i = 0; i < maxQPower; i++) {
        QBdtNodeInterfacePtr leaf1 = root;
        QBdtNodeInterfacePtr leaf2 = toCompare->root;
        complex scale1 = leaf1->scale;
        complex scale2 = leaf2->scale;
        bitLenInt j;
        for (j = 0; j < qubitCount; j++) {
            if (IS_NORM_0(scale1)) {
                break;
            }
            leaf1 = leaf1->branches[SelectBit(i, j)];
            scale1 *= leaf1->scale;
        }
        if (j < qubitCount) {
            continue;
        }
        for (j = 0; j < qubitCount; j++) {
            if (IS_NORM_0(scale2)) {
                break;
            }
            leaf2 = leaf2->branches[SelectBit(i, j)];
            scale2 *= leaf2->scale;
        }
        if (j < qubitCount) {
            continue;
        }
        projection += conj(scale2) * scale1;
    }

    return ONE_R1 - clampProb(norm(projection));
}

complex QBdt::GetAmplitude(bitCapInt perm)
{
    if (stateVec) {
        return stateVec->GetAmplitude(perm);
    }

    QBdtNodeInterfacePtr leaf = root;
    complex scale = leaf->scale;
    bitLenInt j;
    for (j = 0; j < qubitCount; j++) {
        if (IS_NORM_0(scale) || !leaf->branches[0]) {
            break;
        }
        leaf = leaf->branches[SelectBit(perm, j)];
        scale *= leaf->scale;
    }

    if (!IS_NORM_0(scale) && (j < qubitCount)) {
        scale *= NODE_TO_QINTERFACE(leaf)->GetAmplitude(perm >> j);
    }

    return scale;
}

bitLenInt QBdt::Compose(QBdtPtr toCopy, bitLenInt start)
{
    ResetStateVector();
    toCopy->ResetStateVector();

    root->InsertAtDepth(toCopy->root, start, toCopy->qubitCount);
    SetQubitCount(qubitCount + toCopy->qubitCount);

    return start;
}

QInterfacePtr QBdt::Decompose(bitLenInt start, bitLenInt length)
{
    QBdtPtr dest = std::make_shared<QBdt>(qubitCount, length, rand_generator, ONE_CMPLX, doNormalize, randGlobalPhase,
        false, -1, (hardware_rand_generator == NULL) ? false : true, false, (real1_f)amplitudeFloor);

    Decompose(start, dest);

    return dest;
}

void QBdt::DecomposeDispose(bitLenInt start, bitLenInt length, QBdtPtr dest)
{
    ResetStateVector();

    if (dest) {
        dest->ResetStateVector();
        dest->root = root->RemoveSeparableAtDepth(start, length);
    } else {
        root->RemoveSeparableAtDepth(start, length);
    }
    SetQubitCount(qubitCount - length);

    root->Prune(qubitCount);
}

real1_f QBdt::Prob(bitLenInt qubit)
{
    if (stateVec) {
        return stateVec->Prob(qubit);
    }
    if (!root->branches[0]) {
        return NODE_TO_QINTERFACE(root)->Prob(qubit);
    }

    const bitCapInt qPower = pow2(qubit);

    std::map<QInterfacePtr, real1> qiProbs;

    real1 oneChance = ZERO_R1;
    for (bitCapInt i = 0; i < qPower; i++) {
        QBdtNodeInterfacePtr leaf = root;
        complex scale = leaf->scale;
        bitLenInt j;
        for (j = 0; j < qubit; j++) {
            if (IS_NORM_0(scale) || !leaf->branches[0]) {
                break;
            }
            leaf = leaf->branches[SelectBit(i, j)];
            scale *= leaf->scale;
        }

        if (IS_NORM_0(scale)) {
            continue;
        }

        if ((j == qubit) && leaf->branches[1]) {
            oneChance += norm(scale * leaf->branches[1]->scale);
            continue;
        }

        // Phase effects don't matter, for probability expectation.
        QInterfacePtr qi = NODE_TO_QINTERFACE(leaf);
        if (qiProbs.find(qi) == qiProbs.end()) {
            qiProbs[qi] = sqrt(NODE_TO_QINTERFACE(leaf)->Prob(qubit - j));
        }
        oneChance += norm(scale * qiProbs[qi]);
    }

    return clampProb(oneChance);
}

real1_f QBdt::ProbAll(bitCapInt perm)
{
    if (stateVec) {
        return stateVec->ProbAll(perm);
    }

    QBdtNodeInterfacePtr leaf = root;
    complex scale = leaf->scale;
    bitLenInt j;
    for (j = 0; j < qubitCount; j++) {
        if (IS_NORM_0(scale) || !leaf->branches[0]) {
            break;
        }
        leaf = leaf->branches[SelectBit(perm, j)];
        scale *= leaf->scale;
    }

    if (!IS_NORM_0(scale) && (j < qubitCount)) {
        scale *= NODE_TO_QINTERFACE(leaf)->GetAmplitude(perm >> j);
    }

    return clampProb(norm(scale));
}

bool QBdt::ForceM(bitLenInt qubit, bool result, bool doForce, bool doApply)
{
    if (stateVec) {
        return stateVec->ForceM(qubit, result, doForce, doApply);
    }
    if (!root->branches[0]) {
        return NODE_TO_QINTERFACE(root)->ForceM(qubit, result, doForce, doApply);
    }

    const real1_f oneChance = Prob(qubit);
    if (oneChance >= ONE_R1) {
        result = true;
    } else if (oneChance <= ZERO_R1) {
        result = false;
    } else if (!doForce) {
        result = (Rand() <= oneChance);
    }

    if (!doApply) {
        return result;
    }

    const bitCapInt qPower = pow2(qubit);
    root->scale = GetNonunitaryPhase();

    std::set<QInterfacePtr> qis;

    for (bitCapInt i = 0; i < qPower; i++) {
        QBdtNodeInterfacePtr leaf = root;
        bitLenInt j;
        for (j = 0; j < qubit; j++) {
            if (IS_NORM_0(leaf->scale) || !leaf->branches[0]) {
                break;
            }
            leaf->Branch();
            leaf = leaf->branches[SelectBit(i, j)];
        }

        if (IS_NORM_0(leaf->scale)) {
            continue;
        }

        if ((j < qubit) || !leaf->branches[1]) {
            QInterfacePtr qi = NODE_TO_QINTERFACE(leaf);
            if (qis.find(qi) == qis.end()) {
                NODE_TO_QINTERFACE(leaf)->ForceM(qubit - j, result, false, true);
                qis.insert(qi);
            }
            continue;
        }

        leaf->Branch();

        QBdtNodeInterfacePtr& b0 = leaf->branches[0];
        QBdtNodeInterfacePtr& b1 = leaf->branches[1];

        if (result) {
            if (IS_NORM_0(b1->scale)) {
                throw std::runtime_error("QBdt::ForceM() forced 0 probability!");
            }
            b0->SetZero();
            b1->scale /= abs(b1->scale);
        } else {
            if (IS_NORM_0(b0->scale)) {
                throw std::runtime_error("QBdt::ForceM() forced 0 probability!");
            }
            b0->scale /= abs(b0->scale);
            b1->SetZero();
        }
    }

    root->Prune(qubit);

    return result;
}

bitCapInt QBdt::MAll()
{
    if (stateVec) {
        return stateVec->MAll();
    }

    bitCapInt result = 0;
    QBdtNodeInterfacePtr leaf = root;
    bool bitResult = false;
    for (bitLenInt i = 0; i < qubitCount; i++) {
        if (!leaf->branches[0]) {
            result |= NODE_TO_QINTERFACE(leaf)->MAll() << i;
            break;
        }

        real1_f oneChance = clampProb(norm(leaf->branches[1]->scale));
        if (oneChance >= ONE_R1) {
            bitResult = true;
        } else if (oneChance <= ZERO_R1) {
            bitResult = false;
        } else {
            bitResult = (Rand() <= oneChance);
        }

        if (bitResult) {
            leaf->branches[0]->SetZero();
            leaf->branches[1]->scale = ONE_CMPLX;
            leaf = leaf->branches[1];
            result |= pow2(i);
        } else {
            leaf->branches[0]->scale = ONE_CMPLX;
            leaf->branches[1]->SetZero();
            leaf = leaf->branches[0];
        }
    }

    SetPermutation(result);

    return result;
}

#define IS_REAL_0(r) (abs(r) <= FP_NORM_EPSILON)
#define IS_CTRLED_CLIFFORD(top, bottom)                                                                                \
    ((IS_REAL_0(std::real(top)) || IS_REAL_0(std::imag(top))) && (IS_SAME(top, bottom) || IS_SAME(top, -bottom)))
#define IS_CLIFFORD(mtrx)                                                                                              \
    (IS_SAME(mtrx[0], mtrx[1]) || IS_SAME(mtrx[0], -mtrx[1]) || IS_SAME(mtrx[0], I_CMPLX * mtrx[1]) ||                 \
        IS_SAME(mtrx[0], -I_CMPLX * mtrx[1])) &&                                                                       \
        (IS_SAME(mtrx[0], mtrx[2]) || IS_SAME(mtrx[0], -mtrx[2]) || IS_SAME(mtrx[0], I_CMPLX * mtrx[2]) ||             \
            IS_SAME(mtrx[0], -I_CMPLX * mtrx[2])) &&                                                                   \
        (IS_SAME(mtrx[0], mtrx[3]) || IS_SAME(mtrx[0], -mtrx[3]) || IS_SAME(mtrx[0], I_CMPLX * mtrx[3]) ||             \
            IS_SAME(mtrx[0], -I_CMPLX * mtrx[3]))

void QBdt::ApplySingle(const complex* mtrx, bitLenInt target)
{
    if (IS_NORM_0(mtrx[1]) && IS_NORM_0(mtrx[2]) && (randGlobalPhase || IS_NORM_0(ONE_CMPLX - mtrx[0])) &&
        IS_NORM_0(mtrx[0] - mtrx[3])) {
        return;
    }

    if (stateVec) {
        stateVec->Mtrx(mtrx, target);
        return;
    }

    if (!root->branches[0]) {
        try {
            NODE_TO_QINTERFACE(root)->Mtrx(mtrx, target);
        } catch (const std::domain_error&) {
            FallbackMtrx(mtrx, target);
        }

        return;
    }

    if (!IS_CLIFFORD(mtrx)) {
        FallbackMtrx(mtrx, target);
        return;
    }

    const bitCapInt qPower = pow2(target);

    std::map<QInterfacePtr, bitLenInt> qis;
    std::set<QBdtNodeInterfacePtr> qns;

#if ENABLE_COMPLEX_X2
    const complex2 mtrxCol1(mtrx[0], mtrx[2]);
    const complex2 mtrxCol2(mtrx[1], mtrx[3]);
#endif

    par_for_qbdt(0, qPower, [&](const bitCapInt& i, const int& cpu) {
        QBdtNodeInterfacePtr leaf = root;
        bitLenInt j;
        QBdtNodeInterfacePtr parent = NULL;
        for (j = 0; j < target; j++) {
            if (IS_NORM_0(leaf->scale)) {
                // WARNING: Mutates loop control variable!
                return (bitCapInt)(pow2(target - j) - ONE_BCI);
            }
            if (!leaf->branches[0]) {
                leaf = leaf->PopSpecial();
                if (!j) {
                    root = leaf;
                } else {
                    parent->branches[SelectBit(i, j)] = leaf;
                }
            }

            leaf->Branch();
            parent = leaf;
            leaf = leaf->branches[SelectBit(i, target - (j + 1U))];
        }

        if (IS_NORM_0(leaf->scale)) {
            return (bitCapInt)0U;
        }

#if ENABLE_COMPLEX_X2
        leaf->Apply2x2(mtrxCol1, mtrxCol2, qubitCount - target);
#else
        leaf->Apply2x2(mtrx, qubitCount - target);
#endif
        qns.insert(leaf);

        return (bitCapInt)0U;
    });

    root->Prune(target);
}

void QBdt::ApplyControlledSingle(
    const complex* mtrx, const bitLenInt* controls, bitLenInt controlLen, bitLenInt target, bool isAnti)
{
    if (stateVec) {
        if (isAnti) {
            stateVec->MACMtrx(controls, controlLen, mtrx, target);
        } else {
            stateVec->MCMtrx(controls, controlLen, mtrx, target);
        }
        return;
    }

    if (!root->branches[0]) {
        try {
            QInterfacePtr qi = NODE_TO_QINTERFACE(root);
            if (isAnti) {
                qi->MACMtrx(controls, controlLen, mtrx, target);
            } else {
                qi->MCMtrx(controls, controlLen, mtrx, target);
            }
        } catch (const std::domain_error&) {
            FallbackMCMtrx(mtrx, controls, controlLen, target, isAnti);
        }

        return;
    }

    if ((controlLen > 1U) ||
        !((IS_NORM_0(mtrx[1]) && IS_NORM_0(mtrx[2]) && IS_CTRLED_CLIFFORD(mtrx[0], mtrx[1])) ||
            (IS_NORM_0(mtrx[0]) && IS_NORM_0(mtrx[3]) && IS_CTRLED_CLIFFORD(mtrx[1], mtrx[2])))) {
        FallbackMCMtrx(mtrx, controls, controlLen, target, isAnti);
        return;
    }

    std::vector<bitLenInt> controlVec(controlLen);
    std::copy(controls, controls + controlLen, controlVec.begin());
    std::sort(controlVec.begin(), controlVec.end());
    const bool isSwapped = target < controlVec.back();
    if (isSwapped) {
        Swap(target, controlVec.back());
        std::swap(target, controlVec.back());
    }

    bitCapInt controlMask = 0U;
    for (bitLenInt c = 0U; c < controlLen; c++) {
        controlMask |= pow2(target - (controlVec[c] + 1U));
    }

    const bitCapInt controlPerm = isAnti ? 0U : controlMask;
    const bitCapInt qPower = pow2(target);
#if ENABLE_COMPLEX_X2
    const complex2 mtrxCol1(mtrx[0], mtrx[2]);
    const complex2 mtrxCol2(mtrx[1], mtrx[3]);
#endif

    std::map<QInterfacePtr, bitLenInt> qis;
    std::set<QBdtNodeInterfacePtr> qns;

    par_for_qbdt(0, qPower, [&](const bitCapInt& i, const int& cpu) {
        if ((i & controlMask) != controlPerm) {
            return (bitCapInt)(controlMask - ONE_BCI);
        }

        QBdtNodeInterfacePtr leaf = root;
        QBdtNodeInterfacePtr parent = NULL;
        bitLenInt j;
        for (j = 0; j < target; j++) {
            if (IS_NORM_0(leaf->scale)) {
                // WARNING: Mutates loop control variable!
                return (bitCapInt)(pow2(target - j) - ONE_BCI);
            }
            if (!leaf->branches[0]) {
                leaf = leaf->PopSpecial();
                if (!j) {
                    root = leaf;
                } else {
                    parent->branches[SelectBit(i, j)] = leaf;
                }
            }

            leaf->Branch();
            parent = leaf;
            leaf = leaf->branches[SelectBit(i, target - (j + 1U))];
        }

        if (IS_NORM_0(leaf->scale)) {
            return (bitCapInt)0U;
        }

#if ENABLE_COMPLEX_X2
        leaf->Apply2x2(mtrxCol1, mtrxCol2, qubitCount - target);
#else
        leaf->Apply2x2(mtrx, qubitCount - target);
#endif
        qns.insert(leaf);

        return (bitCapInt)0U;
    });

    root->Prune(target);
    // Undo isSwapped.
    if (isSwapped) {
        Swap(target, controlVec.back());
    }
}

void QBdt::Mtrx(const complex* mtrx, bitLenInt target) { ApplySingle(mtrx, target); }

void QBdt::MCMtrx(const bitLenInt* controls, bitLenInt controlLen, const complex* mtrx, bitLenInt target)
{
    if (!controlLen) {
        Mtrx(mtrx, target);
    } else if (IS_NORM_0(mtrx[1]) && IS_NORM_0(mtrx[2])) {
        MCPhase(controls, controlLen, mtrx[0], mtrx[3], target);
    } else if (IS_NORM_0(mtrx[0]) && IS_NORM_0(mtrx[3])) {
        MCInvert(controls, controlLen, mtrx[1], mtrx[2], target);
    } else {
        ApplyControlledSingle(mtrx, controls, controlLen, target, false);
    }
}

void QBdt::MACMtrx(const bitLenInt* controls, bitLenInt controlLen, const complex* mtrx, bitLenInt target)
{

    if (!controlLen) {
        Mtrx(mtrx, target);
    } else if (IS_NORM_0(mtrx[1]) && IS_NORM_0(mtrx[2])) {
        MACPhase(controls, controlLen, mtrx[0], mtrx[3], target);
    } else if (IS_NORM_0(mtrx[0]) && IS_NORM_0(mtrx[3])) {
        MACInvert(controls, controlLen, mtrx[1], mtrx[2], target);
    } else {
        ApplyControlledSingle(mtrx, controls, controlLen, target, true);
    }
}

void QBdt::MCPhase(
    const bitLenInt* controls, bitLenInt controlLen, complex topLeft, complex bottomRight, bitLenInt target)
{
    if (!controlLen) {
        Phase(topLeft, bottomRight, target);
        return;
    }

    const complex mtrx[4] = { topLeft, ZERO_CMPLX, ZERO_CMPLX, bottomRight };
    if (!IS_NORM_0(ONE_CMPLX - topLeft)) {
        ApplyControlledSingle(mtrx, controls, controlLen, target, false);
        return;
    }

    if (IS_NORM_0(ONE_CMPLX - bottomRight)) {
        return;
    }

    std::unique_ptr<bitLenInt[]> lControls = std::unique_ptr<bitLenInt[]>(new bitLenInt[controlLen + 1U]);
    std::copy(controls, controls + controlLen, lControls.get());
    lControls[controlLen] = target;
    std::sort(lControls.get(), lControls.get() + controlLen + 1U);
    target = lControls[controlLen];

    ApplyControlledSingle(mtrx, lControls.get(), controlLen, target, false);
}

void QBdt::MCInvert(
    const bitLenInt* controls, bitLenInt controlLen, complex topRight, complex bottomLeft, bitLenInt target)
{
    if (!controlLen) {
        Invert(topRight, bottomLeft, target);
        return;
    }

    const complex mtrx[4] = { ZERO_CMPLX, topRight, bottomLeft, ZERO_CMPLX };
    if (!IS_NORM_0(ONE_CMPLX - topRight) || !IS_NORM_0(ONE_CMPLX - bottomLeft)) {
        ApplyControlledSingle(mtrx, controls, controlLen, target, false);
        return;
    }

    std::vector<bitLenInt> controlVec(controlLen);
    std::copy(controls, controls + controlLen, controlVec.begin());
    std::sort(controlVec.begin(), controlVec.end());

    if (controlVec.back() < target) {
        ApplyControlledSingle(mtrx, controls, controlLen, target, false);
        return;
    }

    H(target);
    MCPhase(controls, controlLen, ONE_CMPLX, -ONE_CMPLX, target);
    H(target);
}
} // namespace Qrack
