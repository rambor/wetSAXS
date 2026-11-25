// IonCloud: charged groups at the buffer pH and their counterion atmosphere as form-factor pseudo-atoms.
//
// Charges (Henderson-Hasselbalch, average pKa values of residues in folded proteins, Grimsley et al. 2009,
// Protein Sci. 18, 247; Arg from Fitch et al. 2015):
//   Asp 3.5 (OD1/OD2)   Glu 4.2 (OE1/OE2)   His 6.6 (ND1/NE2)   Cys 6.8 (SG, free thiols only)
//   Tyr 10.3 (OH)       Lys 10.5 (NZ)       Arg 12.0 (CZ)       N-terminus 7.7 (N)   C-terminus 3.3 (O/OXT)
//   nucleic-acid phosphate -1 (P); 5'-terminal phosphate (OP3) a further -1 with pKa 6.5
//   metal ions in the model (MG, CA, ZN, MN, NA, K, ...) carry their formal charge.
//
// Counterion atmosphere: each charged group k is neutralised by a pseudo-atom placed along the local outward
// normal n_k (from the centroid of atoms within 10 A to the group) at the distance a + lambda_D, the centroid
// of a linearised Gouy-Chapman layer next to a charged surface:
//     r_ion,k = r_k + (a + lambda_D) n_k
// It carries the extra electrons of the ions that neutralise the group (BufferModel::cloudElectronsPerCharge):
//     E_k = |q_k| [ f dN_counter - (1 - f) dN_co ]
// with f = 0.71 for nucleic-acid phosphates (Pabit et al. 2010) and 0.5 for protein groups (linear screening),
// and the form factor F(q) = E_k exp(-q^2 sigma^2 / 2).
// No parameter of the cloud is fitted: a, sigma and f are fixed physical settings.
//
// The amplitudes are added to the particle amplitudes A_p (AtomisticModel), since the cloud's electrons
// are already counted against the buffer they displace.

#ifndef WETSAXS_IONCLOUD_H
#define WETSAXS_IONCLOUD_H

#include <string>
#include <vector>
#include <sastools/vector3.h>
#include <sastools/PDBModel.h>

class BufferModel;

struct ChargeSite {
    vector3 position;      // centred model coordinates
    float charge;          // fractional, in e
    bool nucleic;
    std::string label;     // e.g. "ASP A 12"
};

struct IonCloudSettings {
    float closestApproach = 3.0f;    // a, Angstrom
    float width = 2.5f;              // sigma of the Gaussian form factor, Angstrom
    float nucleicFraction = 0.71f;   // counterion share of the neutralisation, phosphates
    float proteinFraction = 0.50f;   // counterion share, protein groups
    float minCharge = 0.02f;         // groups with |q| below this are ignored
};

class IonCloud {

    std::vector<ChargeSite> sites;
    std::vector<vector3> positions;  // pseudo-atoms, centred coordinates
    std::vector<float> electrons;    // E_k
    float sigma = 2.5f;
    float debyeLength = 0.0f;
    std::vector<std::string> report;

public:
    IonCloud() = default;
    IonCloud(PDBModel & model, const BufferModel & buffer, const IonCloudSettings & settings = IonCloudSettings());

    static std::vector<ChargeSite> chargeSites(PDBModel & model, double pH, float minCharge = 0.02f);

    /*
     * Adds sum_k E_k F(q) j_l(q r_k) Y_lm(r_k) to the partial amplitudes, in the layout of
     * AtomisticModel: index (lmax+1)^2 * q + l*l + l + m.
     */
    void addPartialAmplitudes(unsigned int lmax, const std::vector<float> & qvalues,
                              float * almReal, float * almImag) const;

    const std::vector<ChargeSite> & getSites() const { return sites; }
    const std::vector<vector3> & getPositions() const { return positions; }
    const std::vector<float> & getElectrons() const { return electrons; }
    float netCharge() const;
    float totalElectrons() const;
    const std::vector<std::string> & getReport() const { return report; }
    void writePseudoAtoms(const std::string & filename, const vector3 * centeringVector) const;
};

#endif //WETSAXS_IONCLOUD_H
