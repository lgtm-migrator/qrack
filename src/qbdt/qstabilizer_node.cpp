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
#include "qbdt_qstabilizer_node.hpp"

#define IS_SAME_AMP(a, b) (norm((a) - (b)) <= (REAL1_EPSILON * REAL1_EPSILON))

namespace Qrack {
bool QBdtQStabilizerNode::isEqual(QBdtNodeInterfacePtr r)
{
    if (!r) {
        return false;
    }

    if (this == r.get()) {
        return true;
    }

    if (!IS_SAME_AMP(scale, r->scale)) {
        return false;
    }

    if (IS_NORM_0(scale)) {
        return true;
    }

    if (r->branches[0]) {
        // "this" node is "special," but "r" is not.

        // TODO: Numerically compare this case, since QStabilizerHybrid coalesces "that" equivalence with "this,"
        // stabilizer, by replacing "that" with stabilizer.

        return false;
    }

    QStabilizerPtr rReg = std::dynamic_pointer_cast<QBdtQStabilizerNode>(r)->qReg;

    if (qReg.get() == rReg.get()) {
        return true;
    }

    if (qReg->ApproxCompare(rReg)) {
        qReg = rReg;
        return true;
    }

    return false;
}

bool QBdtQStabilizerNode::isEqualUnder(QBdtNodeInterfacePtr r)
{
    if (!r) {
        return false;
    }

    if (this == r.get()) {
        return true;
    }

    if (IS_NORM_0(scale)) {
        return IS_NORM_0(r->scale);
    }

    if (r->branches[0]) {
        // "this" node is "special," but "r" is not.

        // TODO: Numerically compare this case, since QStabilizerHybrid coalesces "that" equivalence with "this,"
        // stabilizer, by replacing "that" with stabilizer.

        return false;
    }

    QStabilizerPtr rReg = std::dynamic_pointer_cast<QBdtQStabilizerNode>(r)->qReg;

    if (qReg.get() == rReg.get()) {
        return true;
    }

    if (qReg->ApproxCompare(rReg)) {
        qReg = rReg;
        return true;
    }

    return false;
}

QBdtNodeInterfacePtr QBdtQStabilizerNode::PopSpecial()
{
    // NOTE: Stabilizer amplitudes are all-or-nothing equal superpositions in Qrack API and should not need rounding.

    const bitCapInt maxQPower = qReg->GetMaxQPower();
    complex amp0;
    bitCapInt perm;
    for (perm = 0U; perm < maxQPower; perm++) {
        amp0 = qReg->GetAmplitude(perm);
        if (amp0 != ZERO_CMPLX) {
            break;
        }
    }
    const bool isAmp1 = (bool)(perm & 1U);
    complex amp1;
    if (isAmp1) {
        amp1 = amp0;
        amp0 = ZERO_CMPLX;
    } else {
        // 0 index bit is definitely 0, so +1.
        amp1 = qReg->GetAmplitude(perm | 1U);
    }
    const real1 len = (real1)sqrt(norm(amp0) + norm(amp1));
    amp0 /= len;
    amp1 /= len;

    QBdtNodeInterfacePtr q[2];

    QStabilizerPtr qReg0 = std::dynamic_pointer_cast<QStabilizer>(qReg->Clone());
    qReg0->ForceM(0U, false, true);
    qReg0->Dispose(0U, 1U, 0U);
    q[0] = std::make_shared<QBdtQStabilizerNode>(amp0, qReg0);
    q[0]->Prune();

    QStabilizerPtr qReg1 = std::dynamic_pointer_cast<QStabilizer>(qReg->Clone());
    qReg0->ForceM(0U, true, true);
    qReg0->Dispose(0U, 1U, 1U);
    q[1] = std::make_shared<QBdtQStabilizerNode>(amp1, qReg1);
    q[1]->Prune();

    return std::make_shared<QBdtNode>(scale, q);
}

void QBdtQStabilizerNode::Normalize(bitLenInt depth)
{
    if (!depth) {
        return;
    }

    if (IS_NORM_0(scale)) {
        SetZero();
        return;
    }

    if (qReg) {
        qReg->UpdateRunningNorm();
        qReg->NormalizeState();
    }
}

void QBdtQStabilizerNode::Branch(bitLenInt depth)
{
    if (!depth) {
        return;
    }

    if (IS_NORM_0(scale)) {
        SetZero();
        return;
    }

    if (qReg) {
        qReg = std::dynamic_pointer_cast<QStabilizer>(qReg->Clone());
    }
}

void QBdtQStabilizerNode::Prune(bitLenInt depth)
{
    if (IS_NORM_0(scale)) {
        SetZero();
        return;
    }

    const real1_f phaseArg = qReg->FirstNonzeroPhase();
    qReg->NormalizeState(REAL1_DEFAULT_ARG, REAL1_DEFAULT_ARG, -phaseArg);
    scale *= std::polar(ONE_R1, (real1)phaseArg);
}

void QBdtQStabilizerNode::InsertAtDepth(QBdtNodeInterfacePtr b, bitLenInt depth, const bitLenInt& size)
{
    if (IS_NORM_0(scale)) {
        return;
    }

    if (depth) {
        throw std::runtime_error("QBdtQStabilizerNode::InsertAtDepth() not implemented for nonzero depth!");
    }

    QBdtQStabilizerNodePtr bEng = std::dynamic_pointer_cast<QBdtQStabilizerNode>(b);
    qReg->Compose(bEng->qReg, 0U);
}

QBdtNodeInterfacePtr QBdtQStabilizerNode::RemoveSeparableAtDepth(bitLenInt depth, const bitLenInt& size)
{
    if (!size || (IS_NORM_0(scale))) {
        return NULL;
    }

    QBdtQStabilizerNodePtr toRet = std::dynamic_pointer_cast<QBdtQStabilizerNode>(ShallowClone());
    toRet->scale /= abs(toRet->scale);

    if (!qReg) {
        return toRet;
    }

    toRet->qReg = std::dynamic_pointer_cast<QStabilizer>(qReg->Decompose(depth, size));

    return toRet;
}

#if ENABLE_COMPLEX_X2
void QBdtQStabilizerNode::PushSpecial(const complex2& mtrxCol1, const complex2& mtrxCol2, QBdtNodeInterfacePtr& b1)
{
    // const complex mtrx[4] = { mtrxCol1.c[0], mtrxCol2.c[0], mtrxCol1.c[1], mtrxCol2.c[1] };
#else
void QBdtQStabilizerNode::PushSpecial(const complex* mtrx, QBdtNodeInterfacePtr& b1)
{
#endif
    // TODO: QBdtQStabilizerNodes should "pop off" the lowest index stabilizer qubit into QBDT qubits until the gate can
    // be completed.
    throw std::out_of_range("QBdtQStabilizerNode::PushSpecial() not yet implemented!");
}
} // namespace Qrack
