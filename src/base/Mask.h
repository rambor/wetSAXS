//
// Created by xos81802 on 26/07/2026.
//

#ifndef MASK_H
#define MASK_H

#include <cmath>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "Waters.h"

class Mask {
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

    float a_length = 3.48; // Angstroms
    double half;

    std::unordered_set<IVec3, IVec3Hash> indexed_mask_lattice;
    unsigned int total_lattice_points;
    std::unordered_map<IVec3, vector3, IVec3Hash> lattice_r_theta_phi;

    // need to map particle pt -> lattice points with weights
    std::vector<float> rvalues;
    std::vector<float> thetas; // could be aligned memory
    std::vector<float> phis;

    std::vector<IVec3> lattice_order;
    std::vector<float> weights;
    std::vector<float> distances;
    std::vector<unsigned int> neighbors_offset;
    std::vector<unsigned int> lattice_indices_ordered; // potential limit to total number of atoms/waters that can be stored (~21m atoms)

    std::vector<float> cross_Ap_Am;
    std::vector<float> cross_Ax_Am;

    SphericalHarmonics she;
    SphericalBessels sbj;

public:
    Mask(){
        half = a_length / 2.0;
    };


    // Returns the FCC lattice point (conventional cubic cell edge = a,
    // e.g. one corner of the cell at the origin) closest to point
    // Return the single closest FCC lattice site to p, as an exact integer
    // triple (i, j, k) in units of a/2, where (i + j + k) is even.
    IVec3 nearestFCCIndex(const vector3 & p) {
        const double half = a_length / 2.0;

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



    /*
     * these are waters using Waters algorithm against the input PDB structure
     * The waters represent the hydration points which will be used to make the mask
     */
    void readInWaterModel(Waters & waterModel) {

        auto & waters = waterModel.getHydration();

        for (auto & w : waters) {
            auto possibility = nearestFCCIndex(vector3(w.x, w.y, w.z));

            auto find_it = indexed_mask_lattice.find(possibility); // check if water has been already added, if not

            if (find_it == indexed_mask_lattice.end()) { // if not found, add it
                // add the water and the associated block surrounding it
                indexed_mask_lattice.insert({possibility});
                // this water is tangent point for a sphere at distance to the specified water
                //
            }

            // use the water point to add layer of water at 2.81 and 2*2.81 away at a tangent plane
            vector3 point = vector3(w.x, w.y, w.z); // hydrated water position
            vector3 normalized = point.normalize();
            vector3 extended = normalized*2.81 + point;
            pointsTangentToPlaneAtDistanceR(extended, 2.81);
            vector3 extended2 = normalized*(2*2.81) + point; // need points surrounding this parallel to the plane defined at water
            pointsTangentToPlaneAtDistanceR(extended2, 2.81);

            // if water is internal cavity, need method to build the mask
            // above method is too outward looking

        }
    }

    // given a point, assuming it is on the surface of a sphere, calculate the lattice points that are tangent
    // to the sphere at fixed radial distance.
    void pointsTangentToPlaneAtDistanceR(const vector3 & p, float radial_distance=2.81) {

        float dis = p.length();
        if (dis < 1e-14) throw std::runtime_error("tangent point at origin");

        vector3 u_vec = p.normalize();

        vector3 a_vec = std::abs(u_vec.z) < 0.9 ? vector3(0, 0, 1) : vector3(1, 0, 0);
        vector3 e1 = u_vec.cross(a_vec).normalize();
        vector3 e2 = u_vec.cross(e1);

        for (int i =0; i<8; i++) {
            float theta = M_PI*i*0.25f; // first layer in the plane surrounding the water
            vector3 v = radial_distance*(std::cosf(theta)*e1 + std::sinf(theta)*e2) + p; // only keep if not in list

            auto possibility = nearestFCCIndex(v);
            auto find_it = indexed_mask_lattice.find(possibility);
            if (find_it == indexed_mask_lattice.end()) {
                indexed_mask_lattice.insert({possibility});
            }

            v *= 2;
            v -= p;

            possibility = nearestFCCIndex(v);
            find_it = indexed_mask_lattice.find(possibility);
            if (find_it == indexed_mask_lattice.end()) {
                indexed_mask_lattice.insert({possibility});
            }

            // second layer around the plane
            theta = M_PI*(2*i+1)/8.0f;
            v = radial_distance*(std::cosf(theta)*e1 + std::sinf(theta)*e2) + p; // only keep if not in list
            possibility = nearestFCCIndex(v);
            find_it = indexed_mask_lattice.find(possibility);
            if (find_it == indexed_mask_lattice.end()) {
                indexed_mask_lattice.insert({possibility});
            }
        }
    }


    // remove points that overlap with protein model
    // should really be points that are within 1.4 Angstroms + van der Waal's radii
    // then for each point in Particle -> identify/associate the neighboring waters
    // calculate distances, calculate weight
    void createParticleMask(AtomisticModel & model) {

        auto totalPDBAtoms = model.getPDBModel().getTotalCoordinates();

        const float * pX = model.getPDBModel().getCenteredXVec();
        const float * pY = model.getPDBModel().getCenteredYVec();
        const float * pZ = model.getPDBModel().getCenteredZVec();

        float * gaussian_radii = model.getPDBModel().getAtomicGaussianRadii();

        for (unsigned int i=0; i<totalPDBAtoms; i++) { // go over each atom in the structure - remove any hydration waters that are too close

            vector3 p = vector3(pX[i], pY[i], pZ[i]);
            auto atom_position = nearestFCCIndex(p);

            auto find_it = indexed_mask_lattice.find(atom_position);
            if (find_it != indexed_mask_lattice.end()) { // if atom_position overlaps with existing lattice point, remove the lattice point
                // remove the indexed_mask_lattice
                indexed_mask_lattice.erase(find_it);
            }

            float radial_limit = 1.41f+gaussian_radii[i]; // remove points that are within the van der Waals radii
            int search_range = std::ceil(radial_limit/half) + 2;

            for (int h=-search_range; h<search_range; h++) {
                int test_h = h + atom_position.i;
                for (int j=-search_range; j<search_range; j++) {
                    int test_j = j + atom_position.j;
                    for (int k=-search_range; k<search_range; k++) {
                        int test_k = k + atom_position.k;

                        if ((test_h + test_j + test_k) % 2 == 0) {
                            float diff = (p - toReal({test_h, test_j, test_k})).length();
                            if ( diff < radial_limit) {
                                find_it = indexed_mask_lattice.find({test_h, test_j, test_k});
                                if (find_it != indexed_mask_lattice.end()) { // an isolated point will have no neighbors
                                    indexed_mask_lattice.erase(find_it);
                                }
                            }
                        }
                    }
                }
            }
        }

        total_lattice_points = indexed_mask_lattice.size(); // these are my remaining hydration waters to assess
        unsigned int counter = 0;
        for (auto & ml : indexed_mask_lattice) {
            auto vec = toReal({ml.i, ml.j, ml.k});
            float rval = vec.length();
            lattice_r_theta_phi.emplace(IVec3{ml.i, ml.j, ml.k}, vector3(rval,acosf(vec.z/rval),(vec.x == 0 || vec.x == 0.0f) ? 0.0f : atan2f(vec.y, vec.x)));
            lattice_order.push_back({ml.i, ml.j, ml.k});
            counter++;
        }

        weights.clear();
        lattice_indices_ordered.clear();

        neighbors_offset.reserve(totalPDBAtoms);
        neighbors_offset[0] = 0;

        // determine water neighbors for each protein atom
        for (unsigned int i=0; i<totalPDBAtoms; i++) {
            vector3 p = vector3(pX[i], pY[i], pZ[i]);
            auto atom_position = nearestFCCIndex(p);
            int search_range = std::ceil(10/half) + 1;

            // neighbor count
            int counter = 0;

            for (int h=-search_range; h<search_range; h++) {
                int test_h = h + atom_position.i;
                for (int j=-search_range; j<search_range; j++) {
                    int test_j = j + atom_position.j;
                    for (int k=-search_range; k<search_range; k++) {
                        int test_k = k + atom_position.k;

                        if ((test_h + test_j + test_k) % 2 == 0) {
                                auto find_it = indexed_mask_lattice.find({test_h, test_j, test_k});
                                if (find_it != indexed_mask_lattice.end()) { // an isolated point will have no neighbors
                                    float diff = (p - toReal({test_h, test_j, test_k})).length();
                                    if (diff < 10 && diff > 2.3) { // relevant neighbors
                                        // need to map protein atom to lattice indices
                                        auto it = std::find(lattice_order.begin(), lattice_order.end(), IVec3{test_h, test_j, test_k});
                                        lattice_indices_ordered.push_back(it - lattice_order.begin()); // store the indexed lattice point
                                        weights.push_back(calculateWaterWeight(diff));
                                        distances.push_back(diff); // for sin (qr)/qr // just need to track the atomic form factor
                                        counter++;
                                    }
                                }
                        }
                    }
                }
            } // end of lattice search

            neighbors_offset[i+1] = counter + neighbors_offset[i]; // allows me to select which form factor
        }

        // mask is created but need the
    }

    const float * getMaskCross_Ap_Am() { return cross_Ap_Am.data(); }
    const float * getMaskCross_Ax_Am() { return cross_Ax_Am.data(); }

    void calculateMaskCrossTerm(unsigned int totalqvalues, std::vector < float > & qvalues, AtomisticModel & model) {

        cross_Ap_Am.resize(totalqvalues);
        cross_Ax_Am.resize(totalqvalues);

        auto totalPDBAtoms = model.getPDBModel().getTotalCoordinates(); // total atoms in the PDB model (EXCLUDES hydration)
        auto pAtomicNumbers = model.getPDBModel().getAdjustedASFNumbers(); // does not include hydration

        // Grab the atom type list once, outside the parallel region
        std::vector<int> atomTypes;
        for (auto & atom : model.getAtomList()) {
          atomTypes.push_back(atom);
        }

        const size_t distSize = distances.size();


#pragma omp parallel
        {
            // Thread-private copies — each thread gets its own map and scratch buffer
            std::map<int, float> asf;
            for (auto & atom : atomTypes) {
                asf.insert(std::pair<int, float>(atom, 0.0f));
            }
            std::vector<float> qr(distSize);


#pragma omp for schedule(dynamic)
            for (unsigned int qvalue_index = 0; qvalue_index < totalqvalues; qvalue_index++) {

                float qvalue = qvalues[qvalue_index];

                for (auto & ff : asf) {
                    ff.second = functions_WETSAXS::asf(ff.first, qvalue);
                }


                for (size_t i = 0; i < distSize; i++) {
                    qr[i] = qvalue * distances[i];
                }

                float temp_cross_Ap_Am = 0.0f;
                float temp_cross_Ax_Am = 0.0f;
                //
                for (unsigned int i = 0; i < totalPDBAtoms; i++) {
                    float sinc = 0.0f;
                    for (unsigned int v = neighbors_offset[i]; v < neighbors_offset[i + 1]; v++) {
                        float qrVal = qr[v];
                        sinc += sinf(qrVal) / qrVal * weights[v];
                    }

                    temp_cross_Ap_Am += asf[pAtomicNumbers[i]] * sinc;
                    temp_cross_Ax_Am += sinc;
                }

                float waterASF  = asf[99];
                float waterASF2 = waterASF * waterASF; // 2*f_w^2
                // No race here: each thread writes a distinct index
                cross_Ap_Am[qvalue_index] = waterASF * temp_cross_Ap_Am;
                cross_Ax_Am[qvalue_index] = waterASF2 * temp_cross_Ax_Am;
            }
        } // end parallel

        // for (unsigned int i=0; i<totalqvalues; i++) {
        //     std::cout << qvalues[i] << " " << cross_Ap_Am[i] << "  " << cross_Ax_Am[i] << std::endl;
        // }
        // does this become a collection of diffuse 10-Angstroms spheres?
    }


    void calculateMaskCrossTerm(float q_value, AtomisticModel & model, float * cross_Ap_Am, float * cross_Ax_Am) {

        auto totalAtoms = model.getTotalAtoms();
        auto totalPDBAtoms = model.getPDBModel().getTotalCoordinates();
        auto pAtomicNumbers = model.getPDBModel().getAdjustedASFNumbers(); // does not include the hydration

        std::map<int, float> asf; // for each q-value, this is calculated over the number of unique atoms

        for(auto & atom : model.getAtomList()){ // water is element 99
            asf.insert(std::pair<int, float>(atom, functions_WETSAXS::asf(atom, q_value)));
        }

        float waterASF = asf[99];
        float waterASF2 = waterASF*waterASF;

        std::vector<float> qr(distances.size());
        for (unsigned int i=0; i<distances.size(); i++) {
            qr[i] = q_value * distances[i];
        }

        float sinc, * pQr;
        float form_factor;
        *cross_Ap_Am = 0;
        *cross_Ax_Am = 0;

        for (unsigned int i=0; i<totalAtoms; i++) {
            // form factor for atom i
            form_factor = ( i < totalPDBAtoms  ? asf[pAtomicNumbers[i]] : waterASF);

            sinc = 0;
            for (unsigned int v =neighbors_offset[i]; v < neighbors_offset[i+1]; v++) { // over all atoms including hydration waters
                pQr = &qr[v];
                sinc += sinf(*pQr)/(*pQr) * weights[v];
            }
            *cross_Ap_Am += 2*waterASF*form_factor*sinc;
            *cross_Ax_Am += 2*waterASF2*sinc;
        }
        // two sums need to be calculated, A_p x A_m and A_x x A_m
    }

    // Convert an integer FCC lattice index (i, j, k) to real-space coordinates.
    vector3 toReal(const IVec3& idx) {
        return { static_cast<float>(idx.i * half), static_cast<float>(idx.j * half), static_cast<float>(idx.k * half) };
    }

    void writeMaskToPDB() {
        // write out bead model
        std::string nameOf = "mask_model.pdb";
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

        for (auto & mm : indexed_mask_lattice){

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

    float calculateWaterWeight(const float dis) {

        float a_1 = 2.4 ;
        float r_1 = 2.8 ;
        float s_1 = 0.5/(0.28*0.28);
        float a_2 = 0.5 ;
        float r_2 = 4.5 ;
        float s_2 = 0.5/(0.55*0.55);
        float a_3 = 0.2 ;
        float r_3 = 7.1 ;
        float s_3 = 0.5/(0.9*0.9) ;
        float normalization = 1.0/17.5;
        // unsigned int num;
        float minLimit = 2.49*2.49;
        float maxLimit = 10.1*10.1;
        // P(r) distribution of water => r^2*( g(r)-1 )
        // maximum value at P(2.8) is 1, normalized to peak value
        return normalization*dis*dis*((1 + a_1*expf(-std::powf(dis-r_1,2)*s_1) + a_2*expf(-std::powf(dis-r_2,2)*s_2) + a_3*expf(-std::powf(dis-r_3,2)*s_3))*(1-expf(-8*std::powf(dis-2.3,2))) - 1);
    }
};



#endif //MASK_H
