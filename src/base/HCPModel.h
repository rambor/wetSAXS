//
// Created by xos81802 on 20/07/2026.
//

#ifndef HCPMODEL_H
#define HCPMODEL_H
//
// Created by xos81802 on 18/07/2026.
//
#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>
#include <sastools/utils.h>
#include <sastools/vector3.h>

struct IVec3 {
    long i, j, k;

    bool operator==(const IVec3& other) const {
        return i == other.i && j == other.j && k == other.k;
    }
};

// Hash function so IVec3 can be used as an unordered_map / unordered_set key.
struct IVec3Hash {
    std::size_t operator()(const IVec3& v) const noexcept {
        // Simple combine of three integer hashes.
        std::size_t h1 = std::hash<long>{}(v.i);
        std::size_t h2 = std::hash<long>{}(v.j);
        std::size_t h3 = std::hash<long>{}(v.k);
        std::size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

class HCPModel {
    const float AVOGADRO = 6.02214076e23;          // 1/mol
    vector3 a1, a2, a3;
    float a_length;
    std::vector<vector3>coordinates;
    std::vector<IVec3> perPointSite;               // one per input point (duplicates allowed)
    std::vector<IVec3> fccModel;                   // distinct sites actually used

    std::unordered_map<IVec3, std::vector<int>, IVec3Hash> siteToPoints; // site -> point indices

public:
    HCPModel(){
        float radius = 2.8/2; // spacing between points in lattice, not radius
        float sqrt3 = sqrtf(3);
        float inv_rsqrt3= 1.0f/(radius * sqrt3);

        float a = 2.8;
        float c = a*sqrtf(8.0f/3.0f);

        a1 = vector3{ a, 0.0, 0.0 };
        a2 = vector3{ a/2.0, a*std::sqrtf(3.0f)*0.5f, 0.0 };
        a3 = vector3{ 0.0, 0.0, c };

        basisFrac[0] = Vec3{0.0, 0.0, 0.0};
        basisFrac[1] = Vec3{1.0/3.0, 1.0/3.0, 1.0/2.0};

        M.m[0][0]=a1.x; M.m[0][1]=a2.x; M.m[0][2]=a3.x;
        M.m[1][0]=a1.y; M.m[1][1]=a2.y; M.m[1][2]=a3.y;
        M.m[2][0]=a1.z; M.m[2][1]=a2.z; M.m[2][2]=a3.z;

        /*
        float radius = 2.8/2; // spacing between points in lattice, not radius
        float sqrt3 = sqrtf(3);
        float inv_rsqrt3= 1.0f/(radius * sqrt3);

        // create HCP grid model
        float inv_radius = 1.0/radius;
        const vector3 a1 = vector3(0.5*radius, -0.5*sqrt3*radius, 0);
        const vector3 a2 = vector3(0.5*radius,  0.5*sqrt3*radius, 0);
        const vector3 a3 = vector3(0, 0, (float)(2.0f*sqrt(6)/3.0f*radius)); // seems fine, spreads out

        const vector3 basis2 = a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f;
        */
    };


    float number_density_per_Angstrom3(float density_g_cm3=1.0,
                            float molar_mass=18.015){

        float mol_per_cm3 = density_g_cm3 / molar_mass;
        return mol_per_cm3 * AVOGADRO/1e24;
        //# 1 cm^3 = 1e24 A^3
    }

    float simple_cubic_spacing(float density_g_cm3=1.0,
                          float molar_mass=18.015) {
        /*
         * Lattice spacing (nm) for a simple cubic lattice that reproduces
         * the target number density: a = (1 / n)^(1/3)
        */
        float n = number_density_per_Angstrom3(density_g_cm3, molar_mass);
        a_length = 1.0f/cbrtf(n);
    }


    // Convert an integer FCC lattice index (i, j, k) to real-space coordinates.
    vector3 toReal(const IVec3& idx) {
        double half = a_length / 2.0;
        return { static_cast<float>(idx.i * half), static_cast<float>(idx.j * half), static_cast<float>(idx.k * half) };
    }

    // Returns the FCC lattice point (conventional cubic cell edge = a,
    // e.g. one corner of the cell at the origin) closest to point
    // Return the single closest FCC lattice site to p, as an exact integer
    // triple (i, j, k) in units of a/2, where (i + j + k) is even.
    Vec3 nearestHCPPoint(const vector3 & p, const HCPLattice& lat) {
        double bestDist = std::numeric_limits<double>::infinity();
        vector3 best{0,0,0};

        for (int s = 0; s < 2; ++s) {
            vector3 shifted = p - lat.basisCartesian(s);
            vector3 frac = lat.Minv.apply(shifted);  // real-valued (n1, n2, n3)

            long lo[3] = {
                (long)std::floor(frac.x),
                (long)std::floor(frac.y),
                (long)std::floor(frac.z)
            };

            // All 8 floor/ceil combinations around the real-valued fractional
            // coordinate -- guarantees we don't miss the true nearest point
            // due to the non-orthogonal a1/a2 basal vectors.
            for (int bx = 0; bx < 2; ++bx) {
                for (int by = 0; by < 2; ++by) {
                    for (int bz = 0; bz < 2; ++bz) {
                        vector3 n{
                            (double)(lo[0] + bx),
                            (double)(lo[1] + by),
                            (double)(lo[2] + bz)
                        };
                        vector3 candidate = lat.M.apply(n) + lat.basisCartesian(s);
                        double d = dist(p, candidate);
                        if (d < bestDist) {
                            bestDist = d;
                            best = candidate;
                        }
                    }
                }
            }
        }

        return best;
    }

    void addToModel(vector3 & point) {
        // determine nearest FCC lattice point and put into perPointSite
        auto possibility = nearestFCCIndex(point);

        auto find_it = siteToPoints.find(possibility);

        if (find_it == siteToPoints.end()) { // if not found, add it
            perPointSite.push_back(possibility);
            siteToPoints[perPointSite.back()].push_back(perPointSite.size() - 1); // map keeps an un-duplicated list
            fccModel.push_back(siteToPoints.find(perPointSite.back())->first);
        }
    }

    void convertLatticeModelToXYZCoordinates() {
        coordinates.clear();

        std::random_device rd;
        std::mt19937 gen(rd());

        for (auto & mm : fccModel) {
            coordinates.emplace_back(toReal(mm) + addRandomNoiseToModel(gen) );
        }
        // water O-to-O distance 2.8 A, random displacements should be 0 to 0.15 in any direction
        // if two waters move maximally in same direction, the distance will be 2.5 instead of 2.8 (average)
        // liquid water correlation function of water shows closest distance to be 2.5

        // to calculate excluded volume scattering
        // spherical coordinates
        // spherical harmonics and bessels
        // calculate nearest neighbor weights
        // combine for amplitudes.
    }


    /*
     * create random displacement vector of 0 to 0.4 Angstroms to mimic water
     */
    vector3 addRandomNoiseToModel(std::mt19937& rng) {

            std::normal_distribution<double> normal(0.0, 1.0); // can be positive or negative
            double x = normal(rng), y = normal(rng), z = normal(rng);
            double norm = std::sqrt(x * x + y * y + z * z);

            if (norm < 1e-12) { x = 1.0; y = 0.0; z = 0.0; norm = 1.0; }
            x /= norm; y /= norm; z /= norm;

            std::uniform_real_distribution<float> magDist(0, 0.15);
            float mag = magDist(rng);

            return { static_cast<float>(x) * mag, static_cast<float>(y) * mag, static_cast<float>(z) * mag};
    }

    std::vector<vector3> getCoordinates(){ return coordinates;}

    void writeFCCModelToPDB() {
        // write out bead model
        std::string nameOf = "FCC_water_model.pdb";
        //const char * outputFileName = nameOf.c_str();
        FILE * pFile;
        pFile = fopen(nameOf.c_str(), "w");

        unsigned int residue_index=1;
        std::vector<std::string> chains(25);
        chains[0] = "A";
        chains[1] = "B";
        chains[2] = "C";
        chains[3] = "D";
        chains[4] = "E";
        chains[5] = "F";
        chains[6] = "G";
        chains[7] = "H";
        chains[8] = "I";
        chains[9] = "K";
        chains[10] = "L";
        chains[11] = "M";
        chains[12] = "N";
        chains[13] = "O";
        chains[14] = "P";
        chains[15] = "Q";
        chains[16] = "R";
        chains[17] = "S";
        chains[18] = "T";
        chains[19] = "U";
        chains[20] = "V";
        chains[21] = "W";
        chains[22] = "X";
        chains[23] = "Y";
        chains[24] = "J";

        std::string chainID = chains[0];

        //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
        //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        unsigned int residueCnt = 1;
        unsigned locale = 0;

        for (auto & mm : fccModel){

            vector3 b = toReal(mm);

            if (residueCnt%9999 == 0){ // reinitialize counter
                residueCnt = 1;
                residue_index = 1;
                locale += 1;
                if (locale > 0){
                    chainID = chains[locale];
                }
            }

            printAtomLine(pFile, residueCnt, chainID, residue_index, b.x, b.y, b.z );

            residue_index += 1;
            residueCnt++;
        }

        fclose(pFile);
    }


    void writeWaterExVolModelToPDB() {

        if (coordinates.size() < 1) {
            throw std::runtime_error("writing file without FCC coordinates populated from fccModel");
        }

        // write out bead model
        std::string nameOf = "exvol_water_model.pdb";
        //const char * outputFileName = nameOf.c_str();
        FILE * pFile;
        pFile = fopen(nameOf.c_str(), "w");

        unsigned int residue_index=1;
        std::vector<std::string> chains(25);
        chains[0] = "A";
        chains[1] = "B";
        chains[2] = "C";
        chains[3] = "D";
        chains[4] = "E";
        chains[5] = "F";
        chains[6] = "G";
        chains[7] = "H";
        chains[8] = "I";
        chains[9] = "K";
        chains[10] = "L";
        chains[11] = "M";
        chains[12] = "N";
        chains[13] = "O";
        chains[14] = "P";
        chains[15] = "Q";
        chains[16] = "R";
        chains[17] = "S";
        chains[18] = "T";
        chains[19] = "U";
        chains[20] = "V";
        chains[21] = "W";
        chains[22] = "X";
        chains[23] = "Y";
        chains[24] = "J";

        std::string chainID = chains[0];

        //fprintf(pFile, "REMARK BEAD RADIUS %.3f\nREMARK CONTACTS PER BEAD %i\nREMARK VOLUME UPPER: %i LOWER: %i\n", this->bead_radius );
        //  fprintf(pFile, "REMARK TEMP RANGE %.3f => %.4E\nREMARK TOTAL TEMPERATURE STEPS %i\nREMARK TOTAL BEADS %i\n", "ATOM", i+1, "CA", "ALA", "A", residue_index.c_str(), currentBead->getX(), currentBead->getY(), currentBead->getZ() );
        unsigned int residueCnt = 1;
        unsigned locale = 0;

        for (auto & b : coordinates){

            if (residueCnt%9999 == 0){ // reinitialize counter
                residueCnt = 1;
                residue_index = 1;
                locale += 1;
                if (locale > 0){
                    chainID = chains[locale];
                }
            }

            printAtomLine(pFile, residueCnt, chainID, residue_index, b.x, b.y, b.z );

            residue_index += 1;
            residueCnt++;
        }

        fclose(pFile);
    }
};

#endif //HCPMODEL_H
