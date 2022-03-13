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

#include "qstabilizerextended.hpp"
#include "qfactory.hpp"

namespace Qrack {

QStabilizerExtended::QStabilizerExtended(bitLenInt qubitCount, bitLenInt magicQubitCount, bitCapInt perm)
{
    const bitLenInt stabilizerQubitCount = qubitCount - magicQubitCount;
    qbdt = std::dynamic_pointer_cast<QBdt>(
        CreateQuantumInterface({ QINTERFACE_BDT }, magicQubitCount, perm & (pow2(magicQubitCount) - ONE_BCI)));
    qbdt->Attach(std::dynamic_pointer_cast<QStabilizer>(
        CreateQuantumInterface({ QINTERFACE_STABILIZER }, stabilizerQubitCount, perm >> magicQubitCount)));
}
} // namespace Qrack
