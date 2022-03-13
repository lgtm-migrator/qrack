//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2021. All rights reserved.
//
// Adapted from:
//
// CHP: CNOT-Hadamard-Phase
// Stabilizer Quantum Computer Simulator
// by Scott Aaronson
// Last modified June 30, 2004
//
// Thanks to Simon Anders and Andrew Cross for bugfixes
//
// https://www.scottaaronson.com/chp/
//
// Daniel Strano and the Qrack contributers appreciate Scott Aaronson's open sharing of the CHP code, and we hope that
// vm6502q/qrack is one satisfactory framework by which CHP could be adapted to enter the C++ STL. Our project
// philosophy aims to raise the floor of decentralized quantum computing technology access across all modern platforms,
// for all people, not commercialization.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#pragma once

#include "qbdt.hpp"

namespace Qrack {

class QStabilizerExtended;
typedef std::shared_ptr<QStabilizerExtended> QStabilizerExtendedPtr;

class QStabilizerExtended {
protected:
    QBdtPtr qbdt;

public:
    QStabilizerExtended(bitLenInt qubitCount, bitLenInt magicQubitCount = 1U, bitCapInt perm = 0U);

    bitLenInt GetQubitCount() { return qbdt->GetQubitCount(); }

    bitCapInt GetMaxQPower() { return qbdt->GetMaxQPower(); }

    void SetPermutation(const bitCapInt& perm) { qbdt->SetPermutation(perm); }

public:
    /// Apply an X (or NOT) gate to target
    virtual void X(bitLenInt t) { qbdt->X(t); }
    /// Apply a Pauli Y gate to target
    virtual void Y(bitLenInt t) { qbdt->Y(t); }
    /// Apply a phase gate (|0>->|0>, |1>->-|1>, or "Z") to qubit b
    virtual void Z(bitLenInt t) { qbdt->Z(t); }
    /// Apply a Hadamard gate to target
    virtual void H(bitLenInt t) { qbdt->H(t); }
    /// Apply a phase gate (|0>->|0>, |1>->i|1>, or "S") to qubit b
    virtual void S(bitLenInt t) { qbdt->S(t); }
    /// Apply an inverse phase gate (|0>->|0>, |1>->-i|1>, or "S adjoint") to qubit b
    virtual void IS(bitLenInt t) { qbdt->IS(t); }
    /// Apply a phase gate (square root of S, or "T") to qubit b
    virtual void T(bitLenInt t) { qbdt->T(t); }
    /// Apply an inverse phase gate (square root of S adjoint, or "T adjoint") to qubit b
    virtual void IT(bitLenInt t) { qbdt->IT(t); }
    /// Apply a CNOT gate with control and target
    virtual void CX(bitLenInt c, bitLenInt t) { qbdt->CNOT(c, t); }
    /// Apply a CY gate with control and target
    virtual void CY(bitLenInt c, bitLenInt t) { qbdt->CY(c, t); }
    /// Apply a CZ gate with control and target
    virtual void CZ(bitLenInt c, bitLenInt t) { qbdt->CZ(c, t); }
    /// Apply an AntiCNOT gate with control and target
    virtual void AntiCX(bitLenInt c, bitLenInt t) { qbdt->AntiCNOT(c, t); }
    /// Apply an AntiCY gate with control and target
    virtual void AntiCY(bitLenInt c, bitLenInt t) { qbdt->AntiCY(c, t); }
    /// Apply an AntiCZ gate with control and target
    virtual void AntiCZ(bitLenInt c, bitLenInt t) { qbdt->AntiCZ(c, t); }
    /// Swap two qubit indices
    virtual void Swap(bitLenInt q1, bitLenInt q2) { qbdt->Swap(q1, q2); }
    /// Swap two qubit indices. If they are different, apply a phase factor of i
    virtual void ISwap(bitLenInt q1, bitLenInt q2) { qbdt->ISwap(q1, q2); }
    /// Measure a qubit
    virtual bool M(bitLenInt t) { return qbdt->M(t); }
    /// Measure all qubits
    virtual bitCapInt MAll() { return qbdt->MAll(); }
};
} // namespace Qrack
