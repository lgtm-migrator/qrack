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
    QBdtNodeInterfacePtr qn[2];

    // We clone the current stabilizer state.
    q[0] = std::dynamic_pointer_cast<QStabilizer>(qReg->Clone());

    // First, let's divert Z-basis eigenstates, for simplicity.
    const bool isZ = q[0]->IsSeparableZ(0);
    if (isZ) {
        // The bit is in Z-basis eigenstate. We can measure without destroying any information.
        const bool m0 = q[0]->M(0);
        // We no longer need the original qubit.
        q[0]->Dispose(0U, 1U, m0);

        if (m0) {
            qn[0] = std::make_shared<QBdtQStabilizerNode>();
            qn[1] = std::make_shared<QBdtQStabilizerNode>(ONE_CMPLX, q[0]);
        } else {
            qn[0] = std::make_shared<QBdtQStabilizerNode>(ONE_CMPLX, q[0]);
            qn[1] = std::make_shared<QBdtQStabilizerNode>();
        }

        qn[0]->Prune();
        qn[1]->Prune();

        QBdtNodePtr toRet = std::make_shared<QBdtNode>(scale, qn);
        toRet->Prune();

        return toRet;
    }

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
    // Since we diverted all separable cases above, we know her original bit was in a maximally mixed state.
    // Hence, we can "force the measurement," since we're simulating.
    // If Alice measures both bits to be 0, Bob already has the teleported state.

    // We measure Alice's Bell pair half...
    const bool m0 = q[0]->M(0);
    q[1]->ForceM(0, m0);

    // We measure Alice's original bit...
    const bool m1 = q[0]->M(1);
    q[1]->ForceM(1, m1);

    // Clean up the measured bits.
    q[0]->Dispose(0U, 2U, 0U);
    q[1]->Dispose(0U, 2U, 0U);

    // Bob's Bell pair half starts in maximal superposition, without phase.
    complex amps[2] = { SQRT1_2_R1, SQRT1_2_R1 };
    if (m0) {
        // This acts an X gate on Bob's bit.
        // Bob's tree node amplitudes are unchanged, but the branches are swapped.
        q[0].swap(q[1]);
    }
    if (m1) {
        // This acts a Z gate on Bob's bit.
        amps[1] *= -ONE_CMPLX;
    }

    // Initialize and prune the sub-tree, and send it back up the caller.
    qn[0] = std::make_shared<QBdtQStabilizerNode>(SQRT1_2_R1, q[0]);
    qn[1] = std::make_shared<QBdtQStabilizerNode>(SQRT1_2_R1, q[1]);

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
