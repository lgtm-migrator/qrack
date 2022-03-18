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

        QBdtNodeInterfacePtr l = ShallowClone();
        l->PopSpecial();

        return (l->branches[0] == r->branches[0]) && (l->branches[1] == r->branches[1]);
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

        QBdtNodeInterfacePtr l = ShallowClone();
        l->PopSpecial();

        return (l->branches[0] == r->branches[0]) && (l->branches[1] == r->branches[1]);
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
    // NOTE: This is basically the quantum teleportation algorithm.
    QStabilizerPtr q[2];

    // We clone the current stabilizer state.
    q[0] = std::dynamic_pointer_cast<QStabilizer>(qReg->Clone());

    // We add a qubit, set to |0>.
    q[0]->AddQbAt0();

    // We set the new qubit in superposition.
    q[0]->H(0);

    // We clone the NEW CLONE, with the superposed bit.
    q[1] = std::dynamic_pointer_cast<QStabilizer>(q[0]->Clone());

    // If we act an X gate just on the |1> permutation new branch, it's a CNOT in the QBdt tree.
    q[1]->X(0);

    // We have a Bell pair shared half-and-half across the simulation types.
    // (By common convention, stabilizer is "Alice" and QBdt is "Bob.")

    // Alice entangles her bit with the Bell pair.
    q[0]->CNOT(1, 0);
    q[1]->CNOT(1, 0);
    q[0]->H(1);
    q[1]->H(1);

    // At this point, Alice measures both her bits.

    // We measure Alice's Bell pair half...
    const real1_f probQb0 = (q[0]->Prob(0) + q[1]->Prob(0)) / 2U;
    bool m0;
    if (probQb0 <= FP_NORM_EPSILON) {
        m0 = false;
    } else if ((ONE_R1 - probQb0) <= FP_NORM_EPSILON) {
        m0 = true;
    } else {
        m0 = q[0]->Rand() <= probQb0;
    }
    q[0]->ForceM(0U, m0);
    q[1]->ForceM(0U, m0);

    // We measure Alice's original bit...
    const real1_f probQb1 = (q[0]->Prob(1) + q[1]->Prob(1)) / 2U;
    bool m1;
    if (probQb1 <= FP_NORM_EPSILON) {
        m1 = false;
    } else if ((ONE_R1 - probQb1) <= FP_NORM_EPSILON) {
        m1 = true;
    } else {
        m1 = q[1]->Rand() <= probQb1;
    }
    q[0]->ForceM(1U, m1);
    q[1]->ForceM(1U, m1);

    // Bob finishes teleportation based on Alice's measurement.
    complex amps[2] = { SQRT1_2_R1, SQRT1_2_R1 };
    if (m0) {
        std::swap(amps[0], amps[1]);
    }
    if (m1) {
        amps[1] *= -ONE_CMPLX;
    }

    // Clean up the measured bits.
    q[0]->Dispose(0U, 2U, 0U);
    q[1]->Dispose(0U, 2U, 0U);

    // Initialize and prune the sub-tree, and send it back up the caller.

    QBdtNodeInterfacePtr qn[2];
    qn[0] = std::make_shared<QBdtQStabilizerNode>(amps[0], q[0]);
    qn[1] = std::make_shared<QBdtQStabilizerNode>(amps[1], q[1]);

    qn[0]->Prune();
    qn[1]->Prune();

    QBdtNodePtr toRet = std::make_shared<QBdtNode>(scale, qn);
    toRet->Prune();

    return toRet;
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
} // namespace Qrack
