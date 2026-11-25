//
// Solvent-fluctuation correction to buffer-subtracted scattering.
//
// The mean solvent density is handled by Babinet's principle in the amplitudes (A_x, A_w0).
// What remains is the scattering of the density FLUCTUATIONS of the bulk water that the particle
// and its hydration sites displace. With chi(r) the excluded region (protein + hydration sites):
//
//   Delta I_fl(q) = <|F_w(q)|^2> * ( -Sigma(q) - N_x S(q) )      buffer subtracted as is
//   Delta I_fl(q) = <|F_w(q)|^2> * ( -Sigma(q) )                 buffer scaled by (1 - phi)
//
//   S(q)     = 1 + 4 pi rho  int r^2 h(r) sinc(qr) dr,            h = g_OO - 1
//   Sigma(q) = rho^2 int D(r) h(r) sinc(qr) dr
//   D(r)     = 4 pi r^2 V_x (1 - gamma(r))                         surface deficit of the region
//
// gamma(r) is the probability that a point a distance r from a random point of chi is also in chi.
// It is estimated by Monte Carlo on a grid of chi (solvent-excluded surface with a 1.4 A probe).
// Returned intensities are in absolute units (Debye): divide by 4 pi to add to sum_lm |A_lm|^2.
//
#ifndef WETSAXS_SOLVENTFLUCTUATION_H
#define WETSAXS_SOLVENTFLUCTUATION_H

#include <vector>
#include <cstdint>

class AtomisticModel;
class Waters;

class SolventFluctuation {

public:
    SolventFluctuation() = default;

    // centers in the same (centered) frame as the model; radii of the spheres making up the region
    void buildExcludedRegion(const std::vector<float> & x,
                             const std::vector<float> & y,
                             const std::vector<float> & z,
                             const std::vector<float> & radii,
                             float probe = 1.4f,
                             float spacing = 0.5f);

    // excluded region of a hydrated model: protein atoms (sphere of each atom's volume) + hydration sites
    void buildExcludedRegion(AtomisticModel & model, Waters & waters, float probe = 1.4f, float spacing = 0.5f);

    // Monte Carlo estimate of gamma(r) on [0, rmax]
    void computeSurfaceDeficit(unsigned long samples = 4000000UL, unsigned int seed = 1, float rmax = 16.0f, float binWidth = 0.25f);

    float structureFactor(float q) const;   // S(q) of bulk water
    float sigma(float q) const;             // surface integral
    float deltaI(float q, bool bufferScaled) const;  // absolute units

    float getExcludedVolume() const { return volume; }       // A^3
    float getExcludedWaters() const { return volume*rhoMolecules; } // N_x

    // water g_OO(r) used throughout (r < 2.5 A returns 0); S(0) = 0.064
    static float gOO(float r);

    static constexpr float rhoMolecules = 0.03342f; // water molecules per A^3
    static constexpr float sphereRadiusWater = 1.926f; // radius of a 29.9 A^3 sphere

private:
    // grid of the excluded region
    std::vector<uint8_t> grid;
    int nx = 0, ny = 0, nz = 0;
    float h = 0.5f;
    float ox = 0, oy = 0, oz = 0; // origin (center of voxel 0,0,0)
    float volume = 0.0f;

    // gamma(r) at bin centers and the fine integration grid
    std::vector<float> rBins, gammaBins;
    std::vector<float> rFine, weightFine; // weightFine = D(r) h(r) dr (trapezoid), q-independent
    std::vector<float> rS, weightS;       // r^2 h(r) dr for S(q)

    inline long index(int i, int j, int k) const { return ((long)i*ny + j)*nz + k; }
    inline bool inside(float px, float py, float pz) const;
    void prepareIntegrals();
};

#endif //WETSAXS_SOLVENTFLUCTUATION_H
