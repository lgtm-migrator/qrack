#pragma once
#include "common/qrack_types.hpp"

namespace Qrack {

struct MpsShard;
typedef std::shared_ptr<MpsShard> MpsShardPtr;

struct MpsShard {
    complex gate[4U];

    MpsShard()
    {
        gate[0U] = ONE_CMPLX;
        gate[1U] = ZERO_CMPLX;
        gate[2U] = ZERO_CMPLX;
        gate[3U] = ONE_CMPLX;
    }

    MpsShard(const complex* g) { std::copy(g, g + 4, gate); }

    MpsShardPtr Clone() { return std::make_shared<MpsShard>(gate); }

    void Compose(const complex* g)
    {
        complex o[4U];
        std::copy(gate, gate + 4U, o);
        mul2x2((complex*)g, o, gate);
    }

    bool IsPhase() { return (norm(gate[1U]) <= FP_NORM_EPSILON) && (norm(gate[2U]) <= FP_NORM_EPSILON); }

    bool IsInvert() { return (norm(gate[0U]) <= FP_NORM_EPSILON) && (norm(gate[3U]) <= FP_NORM_EPSILON); }

    bool IsIdentity() { return IsPhase() && (norm(gate[0U] - gate[3U]) <= FP_NORM_EPSILON); }

    bool IsX(bool randGlobalPhase = true)
    {
        return IsInvert() && (norm(gate[1U] - gate[2U]) <= FP_NORM_EPSILON) &&
            (randGlobalPhase || (norm(ONE_CMPLX - gate[1U]) <= FP_NORM_EPSILON));
    }

    bool IsY(bool randGlobalPhase = true)
    {
        return IsInvert() && (norm(gate[1U] + gate[2U]) <= FP_NORM_EPSILON) &&
            (randGlobalPhase || (norm(ONE_CMPLX + gate[1U]) <= FP_NORM_EPSILON));
    }

    bool IsZ(bool randGlobalPhase = true)
    {
        return IsPhase() && (norm(gate[0U] + gate[3U]) <= FP_NORM_EPSILON) &&
            (randGlobalPhase || (norm(ONE_CMPLX - gate[0U]) <= FP_NORM_EPSILON));
    }
};

} // namespace Qrack
