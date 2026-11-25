#include "SolventFluctuation.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <deque>
#include "functions_WETSAXS.h"
#include "AtomisticModel.h"
#include "Waters.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
    inline float sincf_(float x){ return (std::fabs(x) < 1e-6f) ? 1.0f - x*x/6.0f : std::sin(x)/x; }
}

/*
 * g_OO(r): wetSAXS g_OO_SQ_matched plus a small term at the first minimum (3.3 A) so that
 * S(0) = rho kT kappa_T = 0.064 and S(q) stays positive at low q (min ~0.04).
 */
float SolventFluctuation::gOO(float r){

    if (r < 2.5f) return 0.0f;

    auto G = [r](float a, float r0, float s){ const float d = r - r0; return a*std::exp(d*d*s); };

    float g = 1.0f
        + G( 1.6022288043f, 2.7992424956f,  -23.8463006847317f)
        + G(-0.5044352646f, 2.512087913132f, -82.5754762463057f)
        + G( 2.835902295f,  4.5106758451f,   -1.16992522513f)
        + G(-2.713383622f,  4.506014238f,    -0.999784342623f)
        + G( 0.07666115f,   6.689911741f,    -2.979445942709f)
        + G(-0.881242444f,  2.287563421f,   -17.1019985978f)
        + G(-0.01531494f,   7.62330643f,     -2.0926947568f)    // g_OO_SQ_matched
        + G( 0.01031f,      3.3f,            -1.388889f);        // compressibility: S(0) = 0.064
    return g;
}

/*
 * Excluded region = solvent-excluded surface (SES) of spheres with the given radii:
 * a voxel is excluded unless a probe sphere that touches no sphere, and is connected to the
 * outside, overlaps it. Internal cavities are therefore filled.
 */
void SolventFluctuation::buildExcludedRegion(const std::vector<float> & x,
                                             const std::vector<float> & y,
                                             const std::vector<float> & z,
                                             const std::vector<float> & radii,
                                             float probe, float spacing){

    const size_t n = x.size();
    if (n == 0 || y.size() != n || z.size() != n || radii.size() != n)
        throw std::invalid_argument("SolventFluctuation::buildExcludedRegion: inconsistent input");

    h = spacing;
    float rmax = *std::max_element(radii.begin(), radii.end());
    const float pad = rmax + 2.0f*probe + 2.0f*h;

    float lx = *std::min_element(x.begin(), x.end()) - pad, hx = *std::max_element(x.begin(), x.end()) + pad;
    float ly = *std::min_element(y.begin(), y.end()) - pad, hy = *std::max_element(y.begin(), y.end()) + pad;
    float lz = *std::min_element(z.begin(), z.end()) - pad, hz = *std::max_element(z.begin(), z.end()) + pad;

    ox = lx; oy = ly; oz = lz;
    nx = (int)std::ceil((hx - lx)/h) + 1;
    ny = (int)std::ceil((hy - ly)/h) + 1;
    nz = (int)std::ceil((hz - lz)/h) + 1;

    const long total = (long)nx*ny*nz;
    std::vector<uint8_t> inflated(total, 0);

    // 1. inflate every sphere by the probe radius
    for(size_t a = 0; a < n; a++){
        const float R = radii[a] + probe, R2 = R*R;
        const int ci = (int)std::lround((x[a]-ox)/h), cj = (int)std::lround((y[a]-oy)/h), ck = (int)std::lround((z[a]-oz)/h);
        const int span = (int)std::ceil(R/h) + 1;
        for(int i = std::max(0, ci-span); i <= std::min(nx-1, ci+span); i++){
            const float dx = ox + i*h - x[a];
            for(int j = std::max(0, cj-span); j <= std::min(ny-1, cj+span); j++){
                const float dy = oy + j*h - y[a];
                const float dxy = dx*dx + dy*dy;
                if (dxy > R2) continue;
                for(int k = std::max(0, ck-span); k <= std::min(nz-1, ck+span); k++){
                    const float dz = oz + k*h - z[a];
                    if (dxy + dz*dz <= R2) inflated[index(i,j,k)] = 1;
                }
            }
        }
    }

    // 2. probe centres reachable from the outside (flood fill over non-inflated voxels from the border)
    std::vector<uint8_t> outside(total, 0);
    std::deque<long> queue;
    auto seed = [&](int i, int j, int k){ long id = index(i,j,k); if (!inflated[id] && !outside[id]){ outside[id] = 1; queue.push_back(id);} };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ seed(i,j,0); seed(i,j,nz-1); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ seed(i,0,k); seed(i,ny-1,k); }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ seed(0,j,k); seed(nx-1,j,k); }

    while(!queue.empty()){
        const long id = queue.front(); queue.pop_front();
        const int k = (int)(id % nz), j = (int)((id / nz) % ny), i = (int)(id / ((long)ny*nz));
        if (i > 0) seed(i-1,j,k);
        if (i < nx-1) seed(i+1,j,k);
        if (j > 0) seed(i,j-1,k);
        if (j < ny-1) seed(i,j+1,k);
        if (k > 0) seed(i,j,k-1);
        if (k < nz-1) seed(i,j,k+1);
    }

    // 3. excluded = everything not swept by an outside probe; only probe centres next to the
    //    inflated surface can reach excluded voxels, so stamp only those
    grid.assign(total, 1);
    const float P2 = probe*probe;
    const int pspan = (int)std::ceil(probe/h) + 1;

    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){
        const long id = index(i,j,k);
        if (!outside[id]) continue;
        grid[id] = 0;
        bool boundary = (i>0 && inflated[index(i-1,j,k)]) || (i<nx-1 && inflated[index(i+1,j,k)]) ||
                        (j>0 && inflated[index(i,j-1,k)]) || (j<ny-1 && inflated[index(i,j+1,k)]) ||
                        (k>0 && inflated[index(i,j,k-1)]) || (k<nz-1 && inflated[index(i,j,k+1)]);
        if (!boundary) continue;
        for(int a = std::max(0,i-pspan); a <= std::min(nx-1,i+pspan); a++){
            const float dx = (a-i)*h;
            for(int b = std::max(0,j-pspan); b <= std::min(ny-1,j+pspan); b++){
                const float dy = (b-j)*h;
                for(int c = std::max(0,k-pspan); c <= std::min(nz-1,k+pspan); c++){
                    const float dz = (c-k)*h;
                    if (dx*dx+dy*dy+dz*dz <= P2) grid[index(a,b,c)] = 0;
                }
            }
        }
    }

    long count = 0;
    for(auto v : grid) count += v;
    volume = (float)count*h*h*h;
}

void SolventFluctuation::buildExcludedRegion(AtomisticModel & model, Waters & waters, float probe, float spacing){

    auto & pdb = model.getPDBModel();
    const unsigned int nAtoms = pdb.getTotalCoordinates();
    const float * px = pdb.getCenteredXVec();
    const float * py = pdb.getCenteredYVec();
    const float * pz = pdb.getCenteredZVec();
    const float * pv = pdb.getAtomicVolumeVec();

    std::vector<float> x, y, z, r;
    x.reserve(nAtoms + waters.getTotalWaters());
    y.reserve(x.capacity()); z.reserve(x.capacity()); r.reserve(x.capacity());

    for(unsigned int n = 0; n < nAtoms; n++){
        x.push_back(px[n]); y.push_back(py[n]); z.push_back(pz[n]);
        r.push_back(std::cbrt(3.0f*pv[n]/(4.0f*(float)M_PI)));   // sphere of the atom's volume
    }

    const auto & sites = waters.getHydration();
    for(unsigned int n = 0; n < waters.getTotalWaters(); n++){
        x.push_back(sites[n].x); y.push_back(sites[n].y); z.push_back(sites[n].z);
        r.push_back(sphereRadiusWater);
    }

    buildExcludedRegion(x, y, z, r, probe, spacing);
}

inline bool SolventFluctuation::inside(float px, float py, float pz) const {
    const int i = (int)std::lround((px-ox)/h), j = (int)std::lround((py-oy)/h), k = (int)std::lround((pz-oz)/h);
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return false;
    return grid[index(i,j,k)] != 0;
}

void SolventFluctuation::computeSurfaceDeficit(unsigned long samples, unsigned int seed, float rmax, float binWidth){

    if (grid.empty()) throw std::runtime_error("SolventFluctuation: build the excluded region first");

    std::vector<long> cells;
    for(long id = 0; id < (long)grid.size(); id++) if (grid[id]) cells.push_back(id);

    const int nb = (int)std::ceil(rmax/binWidth);
    std::vector<double> hits(nb, 0.0), tries(nb, 0.0);

    #pragma omp parallel
    {
        std::vector<double> lh(nb, 0.0), lt(nb, 0.0);
        #ifdef _OPENMP
        const unsigned int tid = (unsigned int)omp_get_thread_num();
        #else
        const unsigned int tid = 0;
        #endif
        std::mt19937_64 rng(seed + 7919u*tid);
        std::uniform_int_distribution<long> pickCell(0, (long)cells.size()-1);
        std::uniform_real_distribution<float> uni(0.0f, 1.0f);
        std::normal_distribution<float> gauss(0.0f, 1.0f);

        #pragma omp for schedule(static)
        for(long s = 0; s < (long)samples; s++){
            const long id = cells[pickCell(rng)];
            const int k = (int)(id % nz), j = (int)((id / nz) % ny), i = (int)(id / ((long)ny*nz));
            const float px = ox + (i + uni(rng) - 0.5f)*h;
            const float py = oy + (j + uni(rng) - 0.5f)*h;
            const float pz = oz + (k + uni(rng) - 0.5f)*h;

            float vx = gauss(rng), vy = gauss(rng), vz = gauss(rng);
            const float norm = std::sqrt(vx*vx + vy*vy + vz*vz);
            const float r = uni(rng)*rmax;
            const float f = r/norm;

            const int b = std::min(nb-1, (int)(r/binWidth));
            lt[b] += 1.0;
            if (inside(px + vx*f, py + vy*f, pz + vz*f)) lh[b] += 1.0;
        }

        #pragma omp critical
        for(int b = 0; b < nb; b++){ hits[b] += lh[b]; tries[b] += lt[b]; }
    }

    rBins.resize(nb); gammaBins.resize(nb);
    for(int b = 0; b < nb; b++){
        rBins[b] = (b + 0.5f)*binWidth;
        gammaBins[b] = (tries[b] > 0) ? (float)(hits[b]/tries[b]) : 1.0f;
    }
    prepareIntegrals();
}

void SolventFluctuation::prepareIntegrals(){

    const float dr = 0.005f, rEnd = rBins.back();
    const int n = (int)(rEnd/dr);

    rFine.resize(n); weightFine.resize(n);
    for(int i = 0; i < n; i++){
        const float r = (i + 0.5f)*dr;
        float gam;                                  // linear interpolation of gamma, gamma(0) = 1
        if (r <= rBins.front()){
            gam = 1.0f + (gammaBins.front() - 1.0f)*(r/rBins.front());
        } else {
            auto it = std::upper_bound(rBins.begin(), rBins.end(), r);
            if (it == rBins.end()) { gam = gammaBins.back(); }
            else {
                const size_t b = (size_t)(it - rBins.begin());
                const float t = (r - rBins[b-1])/(rBins[b] - rBins[b-1]);
                gam = gammaBins[b-1] + t*(gammaBins[b] - gammaBins[b-1]);
            }
        }
        const float D = 4.0f*(float)M_PI*r*r*volume*(1.0f - gam);
        rFine[i] = r;
        weightFine[i] = rhoMolecules*rhoMolecules*D*(gOO(r) - 1.0f)*dr;
    }

    const float rSend = 20.0f;
    const int ns = (int)(rSend/dr);
    rS.resize(ns); weightS.resize(ns);
    for(int i = 0; i < ns; i++){
        const float r = (i + 0.5f)*dr;
        rS[i] = r;
        weightS[i] = 4.0f*(float)M_PI*rhoMolecules*r*r*(gOO(r) - 1.0f)*dr;
    }
}

float SolventFluctuation::structureFactor(float q) const {
    double s = 1.0;
    for(size_t i = 0; i < rS.size(); i++) s += weightS[i]*sincf_(q*rS[i]);
    return (float)s;
}

float SolventFluctuation::sigma(float q) const {
    double s = 0.0;
    for(size_t i = 0; i < rFine.size(); i++) s += weightFine[i]*sincf_(q*rFine[i]);
    return (float)s;
}

float SolventFluctuation::deltaI(float q, bool bufferScaled) const {
    const float fw = functions_WETSAXS::asf(99, q);
    const float bulk = bufferScaled ? 0.0f : getExcludedWaters()*structureFactor(q);
    return fw*fw*(-sigma(q) - bulk);
}
