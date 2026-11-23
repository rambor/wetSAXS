//
// Created by xos81802 on 08/10/2026.
//

#ifndef WATER_SHELL_H
#define WATER_SHELL_H
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <sastools/vector3.h>

#include "water_shell_data.h"

namespace water_shell {

// Rotation matrix stored as its three columns, so R*p = c0*p.x + c1*p.y + c2*p.z.
struct Rotation {
    vector3 c0, c1, c2;
    inline vector3 apply(const vector3& p) const { return c0 * p.x + c1 * p.y + c2 * p.z; }
};

// Independent, well-mixed generator for stream `stream` (thread id or task id).
inline std::mt19937 makeRng(std::uint64_t baseSeed, std::uint64_t stream) {
    std::seed_seq seq{std::uint32_t(baseSeed), std::uint32_t(baseSeed >> 32),
                      std::uint32_t(stream), std::uint32_t(stream >> 32)};
    return std::mt19937(seq);
}

// Uniformly random rotation (Shoemake's random unit quaternion).
inline Rotation randomRotation(std::mt19937 & rng) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double u1 = U(rng), u2 = U(rng), u3 = U(rng);
    const double twoPi = 2.0 * M_PI;
    const double a = std::sqrt(1.0 - u1), b = std::sqrt(u1);
    const double w = a * std::sin(twoPi * u2), x = a * std::cos(twoPi * u2);
    const double y = b * std::sin(twoPi * u3), z = b * std::cos(twoPi * u3);
    return {
        vector3(float(1 - 2 * (y * y + z * z)), float(2 * (x * y + z * w)),     float(2 * (x * z - y * w))),
        vector3(float(2 * (x * y - z * w)),     float(1 - 2 * (x * x + z * z)), float(2 * (y * z + x * w))),
        vector3(float(2 * (x * z + y * w)),     float(2 * (y * z - x * w)),     float(1 - 2 * (x * x + y * y))),
    };
}

inline int numConstellations() { return water_shell_data::kNumConstellations; }

// Number of points in constellation `which` (use to reserve buffers up front).
inline int constellationSize(int which) {
    return water_shell_data::kOffsets[which + 1] - water_shell_data::kOffsets[which];
}

// Deterministic placement: constellation `which`, rotated by R, translated to `center`.
// Fills `out` (cleared first; its capacity is reused, so no allocation after the first call).
inline void placeWaterShell(const vector3& center, const Rotation& R, int which,
                            std::vector<vector3>& out, float rMax = 14.0f) {
    using namespace water_shell_data;
    const float rMax2 = rMax * rMax;
    out.clear();
    for (int i = kOffsets[which]; i < kOffsets[which + 1]; ++i) {
        const vector3 p(simde_mm_load_ps(kPoints[i]));  // one aligned load, w = 0
        if (p.sqlength() > rMax2) break;                 // sorted by distance
        out.push_back(center + R.apply(p));
    }
}

// Random constellation (or `which` if >= 0) in a random orientation, written into `out`.
// Use this form inside parallel loops.
inline void randomWaterShell(const vector3& center, std::mt19937 & rng, std::vector<vector3>& out,
                             float rMax = 14.0f, int which = -1) {
    if (which < 0 || which >= water_shell_data::kNumConstellations)
        which = std::uniform_int_distribution<int>(0, water_shell_data::kNumConstellations - 1)(rng);
    placeWaterShell(center, randomRotation(rng), which, out, rMax);
}

// Convenience form that returns a new vector (allocates on every call).
inline std::vector<vector3> randomWaterShell(const vector3 & center, std::mt19937 & rng,
                                             float rMax = 14.0f, int which = -1) {
    std::vector<vector3> out;
    randomWaterShell(center, rng, out, rMax, which);
    return out;
}

}  // namespace water_shell
#endif //WATER_SHELL_H
