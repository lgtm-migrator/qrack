//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano 2017. All rights reserved.
//
// This is a header-only, quick-and-dirty, multithreaded, universal quantum register
// simulation, allowing (nonphysical) register cloning and direct measurement of
// probability and phase, to leverage what advantages classical emulation of qubits
// can have.
//
// Licensed under the GNU General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/gpl-3.0.en.html
// for details.

#include <stdint.h>
#include <math.h>
#include <algorithm>
#include <ctime>
#include <random>
#include <stdexcept>
#include <memory>
#include <atomic>
#include <thread>
#include <future>

#ifdef __APPLE__
    #include <OpenCL/cl.hpp>
#else
    #include <CL/cl.hpp>
#endif

#include "complex16simd.hpp"

//#define Complex16 std::complex<double>
#define Complex16 Complex16Simd
#define bitLenInt uint8_t
#define bitCapInt uint64_t
#define bitsInByte 8

namespace Qrack {
	template <class BidirectionalIterator>
	void reverse (BidirectionalIterator first, BidirectionalIterator last, bitCapInt stride);

	template <class BidirectionalIterator>
	void rotate (BidirectionalIterator first, BidirectionalIterator middle, BidirectionalIterator last,  bitCapInt stride);

	/// "Qrack::OCLSingleton" manages the single OpenCL context
	/** "Qrack::OCLSingleton" manages the single OpenCL context. */
	class OCLSingleton{
		public:
			///Get a pointer to the Instance of the singleton. (The instance will be instantiated, if it does not exist yet.) 
			static OCLSingleton* Instance();
			///If this is the first time instantiating the OpenCL context, you may specify platform number and device number.
			static OCLSingleton* Instance(int plat, int dev);
			///Get a pointer to the OpenCL context
			cl::Context* GetContextPtr();
			///Get a pointer to the OpenCL queue
			cl::CommandQueue* GetQueuePtr();
			///Get a pointer to the Apply2x2 function kernel
			cl::Kernel* GetApply2x2Ptr();
			///Get a pointer to the ROL function kernel
			cl::Kernel* GetROLPtr();
			///Get a pointer to the ROR function kernel
			cl::Kernel* GetRORPtr();
			///Get a pointer to the ADD function kernel
			cl::Kernel* GetADDPtr();
			///Get a pointer to the SUB function kernel
			cl::Kernel* GetSUBPtr();
			///Get a pointer to the ADDC function kernel
			cl::Kernel* GetADDCPtr();
			///Get a pointer to the SUBC function kernel
			cl::Kernel* GetSUBCPtr();

		private:
			std::vector<cl::Platform> all_platforms;
			cl::Platform default_platform;
			std::vector<cl::Device> all_devices;
			cl::Device default_device;
			cl::Context context;
			cl::Program program;
			cl::CommandQueue queue;
			cl::Kernel apply2x2;
			cl::Kernel rol;
			cl::Kernel ror;
			cl::Kernel add;
			cl::Kernel sub;
			cl::Kernel addc;
			cl::Kernel subc;

			OCLSingleton(); // Private so that it can  not be called
			OCLSingleton(int plat, int dev); // Private so that it can  not be called
			OCLSingleton(OCLSingleton const&);             // copy constructor is private
			OCLSingleton& operator=(OCLSingleton const&);  // assignment operator is private
			static OCLSingleton* m_pInstance;

			void InitOCL(int plat, int dev);
	};
	/// The "Qrack::CoherentUnit" class represents one or more coherent quantum processor registers		
	/** The "Qrack::CoherentUnit" class represents one or more coherent quantum processor registers, including primitive bit logic gates and (abstract) opcodes-like methods. */
	class CoherentUnit {
		public:
			///Initialize a coherent unit with qBitCount number of bits, all to |0> state.
			CoherentUnit(bitLenInt qBitCount);
			///Initialize a coherent unit with qBitCount number pf bits, to initState unsigned integer permutation state
			CoherentUnit(bitLenInt qBitCount, bitCapInt initState);
			///PSEUDO-QUANTUM Initialize a cloned register with same exact quantum state as pqs
			CoherentUnit(const CoherentUnit& pqs);

			///Get the count of bits in this register
			int GetQubitCount();
			///PSEUDO-QUANTUM Output the exact quantum state of this register as a permutation basis array of complex numbers
			void CloneRawState(Complex16* output);
			///Generate a random double from 0 to 1
			double Rand();
			///Set |0>/|1> bit basis pure quantum permutation state, as an unsigned int
			void SetPermutation(bitCapInt perm);
			///Set arbitrary pure quantum state, in unsigned int permutation basis
			void SetQuantumState(Complex16* inputState);
			///Combine (a copy of) another CoherentUnit with this one, after the last bit index of this one.
			/** Combine (a copy of) another CoherentUnit with this one, after the last bit index of this one. (If the programmer doesn't want to "cheat," it is left up to them to delete the old coherent unit that was added. */
			void Cohere(CoherentUnit &toCopy);
			///Minimally decohere a set of contigious bits from the full coherent unit.
			/** Minimally decohere a set of contigious bits from the full coherent unit. The length of this coherent unit is reduced by the length of bits decohered, and the bits removed are output in the destination CoherentUnit pointer. The destination object must be initialized to the correct number of bits, in 0 permutation state. */
			void Decohere(bitLenInt start, bitLenInt length, CoherentUnit& destination);
			void Dispose(bitLenInt start, bitLenInt length);

			//Logic Gates:
			///"AND" compare two bits in CoherentUnit, and store result in outputBit
			void AND(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
			///"OR" compare two bits in CoherentUnit, and store result in outputBit
			void OR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
			///"XOR" compare two bits in CoherentUnit, and store result in outputBit
			void XOR(bitLenInt inputBit1, bitLenInt inputBit2, bitLenInt outputBit);
			/// Doubly-controlled not
			void CCNOT(bitLenInt control1, bitLenInt control2, bitLenInt target);
			/// "Anti-doubly-controlled not" - Apply "not" if control bits are both zero, do not apply if either control bit is one.
			void AntiCCNOT(bitLenInt control1, bitLenInt control2, bitLenInt target);
			///Controlled not
			void CNOT(bitLenInt control, bitLenInt target);
			///"Anti-controlled not" - Apply "not" if control bit is zero, do not apply if control bit is one.
			void AntiCNOT(bitLenInt control, bitLenInt target);
			///Hadamard gate
			void H(bitLenInt qubitIndex);
			///Measurement gate
			bool M(bitLenInt qubitIndex);
			///PSEUDO-QUANTUM Direct measure of bit probability to be in |1> state
			double Prob(bitLenInt qubitIndex);
			///PSEUDO-QUANTUM Direct measure of full register probability to be in permutation state
			double ProbAll(bitCapInt fullRegister);
			///PSEUDO-QUANTUM Direct measure of all bit probabilities in register to be in |1> state
			void ProbArray(double* probArray);
			///"Phase shift gate" - Rotates as e^(-i*\theta/2) around |1> state 
			void R1(double radians, bitLenInt qubitIndex);
			///Dyadic fraction "phase shift gate" - Rotates as e^(i*(M_PI * numerator) / denominator) around |1> state
			/** Dyadic fraction "phase shift gate" - Rotates as e^(i*(M_PI * numerator) / denominator) around |1> state. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS AND LACKS DIVISION BY A FACTOR OF TWO. */
			void R1Dyad(int numerator, int denominator, bitLenInt qubitIndex);
			///x axis rotation gate - Rotates as e^(-i*\theta/2) around Pauli x axis 
			void RX(double radians, bitLenInt qubitIndex);
			///Dyadic fraction x axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli x axis
			/** Dyadic fraction x axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli x axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS AND LACKS DIVISION BY A FACTOR OF TWO. */
			void RXDyad(int numerator, int denominator, bitLenInt qubitIndex);
			///y axis rotation gate - Rotates as e^(-i*\theta/2) around Pauli y axis 
			void RY(double radians, bitLenInt qubitIndex);
			///Dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis
			/** Dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS AND LACKS DIVISION BY A FACTOR OF TWO. */
			void RYDyad(int numerator, int denominator, bitLenInt qubitIndex);
			///z axis rotation gate - Rotates as e^(-i*\theta/2) around Pauli z axis 
			void RZ(double radians, bitLenInt qubitIndex);
			///Dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis
			/** Dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS AND LACKS DIVISION BY A FACTOR OF TWO. */
			void RZDyad(int numerator, int denominator, bitLenInt qubitIndex);
			///Set individual bit to pure |0> (false) or |1> (true) state
			void SetBit(bitLenInt qubitIndex1, bool value);
			///Swap values of two bits in register
			void Swap(bitLenInt qubitIndex1, bitLenInt qubitIndex2);
			///NOT gate, which is also Pauli x matrix
			void X(bitLenInt qubitIndex);
			///Apply Pauli Y matrix to bit
			void Y(bitLenInt qubitIndex);
			///Apply Pauli Z matrix to bit
			void Z(bitLenInt qubitIndex);
			///Controlled "phase shift gate"
			/** Controlled "phase shift gate" - if control bit is true, rotates target bit as e^(-i*\theta/2) around |1> state */
			void CR1(double radians, bitLenInt control, bitLenInt target);
			///Controlled dyadic fraction "phase shift gate"
			/** Controlled "phase shift gate" - if control bit is true, rotates target bit as e^(-i*\theta/2) around |1> state */
			void CR1Dyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
			///Controlled x axis rotation
			/** Controlled x axis rotation - if control bit is true, rotates as e^(-i*\theta/2) around Pauli x axis */
			void CRX(double radians, bitLenInt control, bitLenInt target);
			///Controlled dyadic fraction x axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli x axis
			/** Controlled dyadic fraction x axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli x axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS. */
			void CRXDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
			///Controlled y axis rotation
			/** Controlled y axis rotation - if control bit is true, rotates as e^(-i*\theta) around Pauli y axis */
			void CRY(double radians, bitLenInt control, bitLenInt target);
			///Controlled dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis
			/** Controlled dyadic fraction y axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli y axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS. */
			void CRYDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
			///Controlled z axis rotation
			/** Controlled z axis rotation - if control bit is true, rotates as e^(-i*\theta) around Pauli z axis */
			void CRZ(double radians, bitLenInt control, bitLenInt target);
			///Controlled dyadic fraction z axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli z axis
			/** Controlled dyadic fraction z axis rotation gate - Rotates as e^(i*(M_PI * numerator) / denominator) around Pauli z axis. NOTE THAT DYADIC OPERATION ANGLE SIGN IS REVERSED FROM RADIAN ROTATION OPERATORS. */
			void CRZDyad(int numerator, int denominator, bitLenInt control, bitLenInt target);
			///Apply controlled Pauli Y matrix to bit
			void CY(bitLenInt control, bitLenInt target);
			///Apply controlled Pauli Z matrix to bit
			void CZ(bitLenInt control, bitLenInt target);

			//Single register instructions:
			///"AND" compare two bit ranges in CoherentUnit, and store result in range starting at output
			void AND(bitLenInt inputStart1, bitLenInt inputStart2, bitLenInt outputStart, bitLenInt length);
			///"OR" compare two bit ranges in CoherentUnit, and store result in range starting at output
			void OR(bitLenInt inputStart1, bitLenInt inputStart2, bitLenInt outputStart, bitLenInt length);
			///"XOR" compare two bit ranges in CoherentUnit, and store result in range starting at output
			void XOR(bitLenInt inputStart1, bitLenInt inputStart2, bitLenInt outputStart, bitLenInt length);
			///Arithmetic shift left, with last 2 bits as sign and carry
			void ASL(bitLenInt shift, bitLenInt start, bitLenInt length);
			///Arithmetic shift right, with last 2 bits as sign and carry
			void ASR(bitLenInt shift, bitLenInt start, bitLenInt length);
			///Logical shift left, filling the extra bits with |0>
			void LSL(bitLenInt shift, bitLenInt start, bitLenInt length);
			///Logical shift right, filling the extra bits with |0>
			void LSR(bitLenInt shift, bitLenInt start, bitLenInt length);
			/// "Circular shift left" - shift bits left, and carry last bits.
			void ROL(bitLenInt shift, bitLenInt start, bitLenInt length);
			/// "Circular shift right" - shift bits right, and carry first bits.
			void ROR(bitLenInt shift, bitLenInt start, bitLenInt length);
			///Add integer (without sign)
			void INC(bitCapInt toAdd, bitLenInt start, bitLenInt length);
			///Subtract integer (without sign)
			void DEC(bitCapInt toSub, bitLenInt start, bitLenInt length);
			///Add two quantum integers
			/** Add integer of "length" bits in "inStart" to integer of "length" bits in "inOutStart," and store result in "inOutStart." */
			void ADD(const bitLenInt inOutStart, const bitLenInt inStart, const bitLenInt length);
			///Add two quantum integers with carry bit
			/** Add integer of "length" bits in "inStart" to integer of "length" bits in "inOutStart," and store result in "inOutStart." Get carry value from bit at "carryIndex" and place end result into this bit. */
			void ADDC(const bitLenInt inOutStart, const bitLenInt inStart, const bitLenInt length, const bitLenInt carryIndex);
			///Subtract two quantum integers
			/** Subtract integer of "length" bits in "toSub" from integer of "length" bits in "inOutStart," and store result in "inOutStart." */
			void SUB(const bitLenInt inOutStart, const bitLenInt toSub, const bitLenInt length);
			///Subtract two quantum integers with carry bit
			/** Subtract integer of "length" - 1 bits in "toSub" from integer of "length" - 1 bits in "inOutStart," and store result in "inOutStart." Get carry value from bit at "carryIndex" and place end result into this bit. */
			void SUBC(const bitLenInt inOutStart, const bitLenInt toSub, const bitLenInt length, const bitLenInt carryIndex);
			/// Quantum Fourier Transform - Apply the quantum Fourier transform to the register
			void QFT(bitLenInt start, bitLenInt length);

		private:
			double runningNorm;
			bitLenInt qubitCount;
			bitCapInt maxQPower;
			std::unique_ptr<Complex16[]> stateVec;

			std::default_random_engine rand_generator;
			std::uniform_real_distribution<double> rand_distribution;

			OCLSingleton* clObj;
			cl::CommandQueue queue;
			cl::Buffer stateBuffer;
			cl::Buffer cmplxBuffer;
			cl::Buffer ulongBuffer;
			cl::Buffer nrmBuffer;
			cl::Buffer maxBuffer;

			void Apply2x2(bitCapInt offset1, bitCapInt offset2, const Complex16* mtrx,
					const bitLenInt bitCount, const bitCapInt* qPowersSorted, bool doApplyNorm, bool doCalcNorm);
			void ApplySingleBit(bitLenInt qubitIndex, const Complex16* mtrx, bool doCalcNorm);
			void ApplyControlled2x2(bitLenInt control, bitLenInt target, const Complex16* mtrx, bool doCalcNorm);
			void ApplyAntiControlled2x2(bitLenInt control, bitLenInt target, const Complex16* mtrx, bool doCalcNorm);
			void Carry(bitLenInt integerStart, bitLenInt integerLength, bitLenInt carryBit);
			void InitOCL();
			void ReInitOCL();
			void NormalizeState();
			void Reverse(bitLenInt first, bitLenInt last);
			void UpdateRunningNorm();
	};
}
