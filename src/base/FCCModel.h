//
// Created by xos81802 on 18/07/2026.
//

#ifndef FCCMODEL_H
#define FCCMODEL_H
#include <cmath>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sastools/utils.h>
#include <sastools/vector3.h>

#include "AtomisticModel.h"
#include "SphericalBessels.h"
#include "SphericalHarmonics.h"



class FCCModel {
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
            seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2); // Bitwise XOR combination
            seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    const float AVOGADRO = 6.02214076e23;          // 1/mol
    float a_length = 2.8;
    float half = 1.4;

    std::vector<vector3>coordinates;
    std::unordered_set<IVec3, IVec3Hash> indexed_lattice;

public:
    FCCModel(){
        this->simple_cubic_spacing();
    };


    float number_density_per_Angstrom3(float density_g_cm3=1.0,
                            float molar_mass=18.015){

        float mol_per_cm3 = density_g_cm3 / molar_mass;
        return mol_per_cm3 * AVOGADRO/1e24;
        //# 1 cm^3 = 1e24 A^3
    }

    void simple_cubic_spacing(float density_g_cm3=1.0,
                          float molar_mass=18.015) {
        /*
         * Lattice spacing (nm) for a simple cubic lattice that reproduces
         * the target number density: a = (1 / n)^(1/3)
        */
        float n = number_density_per_Angstrom3(density_g_cm3, molar_mass);
        a_length = 1.0f/cbrtf(n);
        half = a_length / 2.0;
    }


    // Convert an integer FCC lattice index (i, j, k) to real-space coordinates.
    vector3 toReal(const IVec3& idx) {
        return { static_cast<float>(idx.i * half), static_cast<float>(idx.j * half), static_cast<float>(idx.k * half) };
    }

    // Returns the FCC lattice point (conventional cubic cell edge = a,
    // e.g. one corner of the cell at the origin) closest to point
    // Return the single closest FCC lattice site to p, as an exact integer
    // triple (i, j, k) in units of a/2, where (i + j + k) is even.
    IVec3 nearestFCCIndex(const vector3 & p) {

        double u = p.x / half;
        double v = p.y / half;
        double w = p.z / half;

        long i = std::lround(u);
        long j = std::lround(v);
        long k = std::lround(w);

        const double du = u - static_cast<double>(i);
        const double dv = v - static_cast<double>(j);
        const double dw = w - static_cast<double>(k);

        long sum = i + j + k;
        long parity = ((sum % 2) + 2) % 2;

        if (parity != 0) {
            double adu = std::fabs(du);
            double adv = std::fabs(dv);
            double adw = std::fabs(dw);

            if (adu >= adv && adu >= adw) {
                i += (du >= 0.0) ? 1 : -1;
            } else if (adv >= adu && adv >= adw) {
                j += (dv >= 0.0) ? 1 : -1;
            } else {
                k += (dw >= 0.0) ? 1 : -1;
            }
        }

        return IVec3{ i, j, k };
    }


    void addToModel(vector3 & point) {
        // determine nearest FCC lattice point and put into perPointSite
        auto possibility = nearestFCCIndex(point);

        if (indexed_lattice.find(possibility) == indexed_lattice.end()) { // if not found, add it
            indexed_lattice.insert({possibility});
        }
    }

    void addPDBModel(AtomisticModel & model) {

        const float * pXv = model.getPDBModel().getCenteredXVec();
        const float * pYv = model.getPDBModel().getCenteredYVec();
        const float * pZv = model.getPDBModel().getCenteredZVec();

        auto total = model.getPDBModel().getTotalCoordinates();
        for (int i = 0; i < total; i++) {
            auto vec = vector3(pXv[i], pYv[i], pZv[i]);
            addToModel(vec);
        }
    }

    void addWaterModel(Waters & model) {
        auto & waters = model.getHydration();

        for (auto & w : waters) {
            auto vec = vector3(w.x, w.y, w.z);
            addToModel(vec);
        }
    }

    void convertLatticeModelToXYZCoordinates() {
        coordinates.clear();

        std::random_device rd;
        std::mt19937 gen(rd());

        for (auto & mm : indexed_lattice) {
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

            std::uniform_real_distribution<float> magDist(0, 0.15); // adds jitter up to 0.15 Angstroms
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

        std::random_device rd;
        std::mt19937 gen(rd());

        for (auto & mm : indexed_lattice){

            vector3 b = toReal(mm) + addRandomNoiseToModel(gen);

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



#endif //FCCMODEL_H
