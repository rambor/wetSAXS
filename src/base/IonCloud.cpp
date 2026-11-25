#include "IonCloud.h"
#include "BufferModel.h"
#include "SphericalBessels.h"

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>
#include <algorithm>

namespace {

    std::string trim(const std::string & s){
        const auto b = s.find_first_not_of(" \t");
        if (b == std::string::npos) return "";
        const auto e = s.find_last_not_of(" \t");
        return s.substr(b, e - b + 1);
    }

    // fraction protonated of a group with this pKa
    inline double protonated(double pKa, double pH){ return 1.0/(1.0 + std::pow(10.0, pH - pKa)); }

    struct Titratable { const char * atoms[2]; double pKa; int chargeProtonated; };   // charge = q_p * protonated + (q_p - 1) * (1 - protonated)

    const std::map<std::string, Titratable> & sideChains(){
        static const std::map<std::string, Titratable> t = {
                {"ASP", {{"OD1", "OD2"}, 3.5, 0}},
                {"GLU", {{"OE1", "OE2"}, 4.2, 0}},
                {"HIS", {{"ND1", "NE2"}, 6.6, 1}},
                {"CYS", {{"SG",  ""   }, 6.8, 0}},
                {"TYR", {{"OH",  ""   }, 10.3, 0}},
                {"LYS", {{"NZ",  ""   }, 10.5, 1}},
                {"ARG", {{"CZ",  ""   }, 12.0, 1}},
        };
        return t;
    }

    const std::map<std::string, int> & metalIons(){
        static const std::map<std::string, int> t = {
                {"MG", 2}, {"CA", 2}, {"ZN", 2}, {"MN", 2}, {"FE", 2}, {"CO", 2}, {"NI", 2}, {"CU", 2},
                {"SR", 2}, {"BA", 2}, {"CD", 2}, {"NA", 1}, {"K", 1}, {"RB", 1}, {"CS", 1}, {"LI", 1},
                {"TL", 1}, {"CL", -1}, {"BR", -1}, {"IOD", -1},
        };
        return t;
    }

    bool isNucleotide(const std::string & r){
        return r == "rA" || r == "rG" || r == "rU" || r == "rC" || r == "DA" || r == "DG" || r == "DT" || r == "DC" ||
               r == "A" || r == "G" || r == "U" || r == "C" || r == "DU";
    }

    bool isAminoAcid(const std::string & r){
        static const char * aa[] = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET",
                                    "PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","HID","HIE","HIP","CYX","SEC"};
        for (const char * a : aa) if (r == a) return true;
        return false;
    }
}

std::vector<ChargeSite> IonCloud::chargeSites(PDBModel & model, double pH, float minCharge){

    const unsigned int n = model.getTotalCoordinates();
    const float * X = model.getCenteredXVec();
    const float * Y = model.getCenteredYVec();
    const float * Z = model.getCenteredZVec();
    const std::string * resName = model.getResiduesVector();
    const std::string * atomName = model.getPointerToAtomTypes();
    auto resIt = model.getResIDIterator();
    auto chainIt = model.getChainIDIterator();

    // residues in file order: key (chain, resid, name) -> atom name -> index
    using Key = std::tuple<std::string, int, std::string>;
    std::vector<Key> order;
    std::map<Key, std::map<std::string, unsigned int>> residues;
    std::map<std::string, std::pair<int,int>> proteinChainEnds;   // chain -> first, last residue number

    for (unsigned int i = 0; i < n; i++){
        const std::string res = trim(resName[i]);
        const Key key{*(chainIt + i), *(resIt + i), res};
        auto & atoms = residues[key];
        if (atoms.empty()) order.push_back(key);
        atoms.emplace(trim(atomName[i]), i);
        if (isAminoAcid(res)){
            auto it = proteinChainEnds.find(std::get<0>(key));
            if (it == proteinChainEnds.end()) proteinChainEnds[std::get<0>(key)] = {std::get<1>(key), std::get<1>(key)};
            else { it->second.first = std::min(it->second.first, std::get<1>(key)); it->second.second = std::max(it->second.second, std::get<1>(key)); }
        }
    }

    auto pos = [&](unsigned int i){ return vector3(X[i], Y[i], Z[i]); };
    auto centroid = [&](const std::map<std::string, unsigned int> & atoms, const char * a, const char * b, vector3 & out){
        auto ia = atoms.find(a);
        if (ia == atoms.end()) return false;
        out = pos(ia->second);
        if (b[0] != '\0'){
            auto ib = atoms.find(b);
            if (ib != atoms.end()) out = (out + pos(ib->second))*0.5f;
        }
        return true;
    };

    // free cysteines: SG with no other SG within 2.5 A
    std::vector<vector3> sulfurs;
    for (const auto & key : order){
        if (std::get<2>(key) != "CYS") continue;
        auto it = residues[key].find("SG");
        if (it != residues[key].end()) sulfurs.push_back(pos(it->second));
    }

    std::vector<ChargeSite> out;
    auto addSite = [&](const vector3 & p, double q, bool nucleic, const Key & key, const std::string & what){
        if (std::fabs(q) < minCharge) return;
        char label[64];
        std::snprintf(label, sizeof(label), "%s %s %d %s", std::get<2>(key).c_str(), std::get<0>(key).c_str(), std::get<1>(key), what.c_str());
        out.push_back({p, (float)q, nucleic, label});
    };

    for (const auto & key : order){
        const std::string & res = std::get<2>(key);
        const auto & atoms = residues[key];
        vector3 p;

        if (isNucleotide(res)){
            auto ip = atoms.find("P");
            if (ip != atoms.end()){
                double q = -1.0;
                if (atoms.count("OP3") || atoms.count("O3P")) q -= 1.0 - protonated(6.5, pH);   // 5'-terminal phosphate
                addSite(pos(ip->second), q, true, key, "phosphate");
            }
            continue;
        }

        if (isAminoAcid(res)){
            const std::string base = (res == "HID" || res == "HIE" || res == "HIP") ? "HIS" : (res == "CYX" ? "" : res);
            auto t = sideChains().find(base);
            if (t != sideChains().end() && centroid(atoms, t->second.atoms[0], t->second.atoms[1], p)){
                bool skip = false;
                if (base == "CYS"){
                    int partners = 0;
                    for (const auto & s : sulfurs){ const vector3 d = s - p; if (d.length() > 0.1f && d.length() < 2.5f) partners++; }
                    skip = partners > 0;
                }
                if (!skip){
                    const double f = protonated(t->second.pKa, pH);
                    addSite(p, t->second.chargeProtonated - (1.0 - f), false, key, "side chain");
                }
            }
            const auto ends = proteinChainEnds[std::get<0>(key)];
            if (std::get<1>(key) == ends.first && centroid(atoms, "N", "", p))
                addSite(p, protonated(7.7, pH), false, key, "N-terminus");
            if (std::get<1>(key) == ends.second && centroid(atoms, "OXT", "O", p))
                addSite(p, -(1.0 - protonated(3.3, pH)), false, key, "C-terminus");
            continue;
        }

        auto m = metalIons().find(res);
        if (m != metalIons().end() && atoms.size() == 1)
            addSite(pos(atoms.begin()->second), m->second, false, key, "ion");
    }
    return out;
}


IonCloud::IonCloud(PDBModel & model, const BufferModel & buffer, const IonCloudSettings & settings) : sigma(settings.width) {

    sites = chargeSites(model, buffer.getPH(), settings.minCharge);
    debyeLength = (float)buffer.debyeLength();
    const float shift = settings.closestApproach + std::min(debyeLength, 100.0f);

    const unsigned int n = model.getTotalCoordinates();
    const float * X = model.getCenteredXVec();
    const float * Y = model.getCenteredYVec();
    const float * Z = model.getCenteredZVec();

    // electrons per unit charge, by group sign and kind
    const double eNucNeg = buffer.cloudElectronsPerCharge(-1, settings.nucleicFraction);
    const double eNucPos = buffer.cloudElectronsPerCharge(+1, settings.nucleicFraction);
    const double eProNeg = buffer.cloudElectronsPerCharge(-1, settings.proteinFraction);
    const double eProPos = buffer.cloudElectronsPerCharge(+1, settings.proteinFraction);

    vector3 center(0, 0, 0);   // coordinates are centred; used only for buried groups
    int buried = 0;
    const float r2 = 100.0f;   // 10 A neighbourhood for the normal

    for (const auto & s : sites){
        vector3 sum(0, 0, 0); int count = 0;
        for (unsigned int i = 0; i < n; i++){
            const float dx = X[i] - s.position.x, dy = Y[i] - s.position.y, dz = Z[i] - s.position.z;
            if (dx*dx + dy*dy + dz*dz < r2){ sum += vector3(X[i], Y[i], Z[i]); count++; }
        }
        vector3 normal = s.position - ((count > 0) ? sum/(float)count : center);
        if (normal.length() < 0.5f){ normal = s.position - center; buried++; }
        if (normal.length() < 1e-3f) normal = vector3(0, 0, 1);
        normal = normal/normal.length();

        positions.push_back(s.position + normal*shift);
        const double perCharge = (s.charge < 0) ? (s.nucleic ? eNucNeg : eProNeg) : (s.nucleic ? eNucPos : eProPos);
        electrons.push_back((float)(std::fabs(s.charge)*perCharge));
    }

    int nNuc = 0; for (const auto & s : sites) if (s.nucleic) nNuc++;
    char b[160];
    report.push_back("Ion cloud :: " + std::to_string(sites.size()) + " charged groups (" + std::to_string(nNuc) + " phosphates)");
    std::snprintf(b, sizeof(b), "Net charge of model (e) :: %+.2f", netCharge()); report.push_back(b);
    std::snprintf(b, sizeof(b), "Pseudo-atoms at a + lambda_D = %.1f + %.1f A along the local normal, width %.1f A", settings.closestApproach, debyeLength, sigma); report.push_back(b);
    std::snprintf(b, sizeof(b), "Cloud electrons :: %+.1f (e per charge: phosphate %.2f, protein -/+ %.2f/%.2f)", totalElectrons(), eNucNeg, eProNeg, eProPos); report.push_back(b);
    if (buried > 0) report.push_back("Groups without a defined normal (buried) :: " + std::to_string(buried));
}

float IonCloud::netCharge() const { double q = 0; for (const auto & s : sites) q += s.charge; return (float)q; }
float IonCloud::totalElectrons() const { double e = 0; for (float v : electrons) e += v; return (float)e; }


void IonCloud::addPartialAmplitudes(unsigned int lmax, const std::vector<float> & qvalues,
                                    float * almReal, float * almImag) const {
    const unsigned int N = (unsigned int)positions.size();
    if (N == 0) return;
    const unsigned int L1 = lmax + 1, perQ = L1*L1;

    std::vector<float> r(N), th(N), ph(N);
    for (unsigned int k = 0; k < N; k++){
        const vector3 & p = positions[k];
        r[k] = p.length();
        th[k] = (r[k] > 0) ? std::acos(std::max(-1.0f, std::min(1.0f, p.z/r[k]))) : 0.0f;
        ph[k] = std::atan2(p.y, p.x);
    }

    // Y_lm for m >= 0 (same convention as SphericalHarmonics: boost, complex)
    std::vector<float> yR((size_t)perQ*N), yI((size_t)perQ*N);
    for (int l = 0; l <= (int)lmax; l++)
        for (int m = 0; m <= l; m++)
            for (unsigned int k = 0; k < N; k++){
                const auto y = boost::math::spherical_harmonic(l, m, (double)th[k], (double)ph[k]);
                const size_t idx = (size_t)(l*l + l + m)*N + k;
                yR[idx] = (float)y.real(); yI[idx] = (float)y.imag();
            }

    std::vector<float> qv(qvalues);
    SphericalBessels sbj(lmax, (unsigned int)qv.size(), N, qv, &r);
    const float * pSBJ = sbj.getpSphericalBessels();

    for (size_t q = 0; q < qv.size(); q++){
        const float ff = std::exp(-0.5f*qv[q]*qv[q]*sigma*sigma);
        for (int l = 0; l <= (int)lmax; l++){
            const size_t bIdx = (size_t)L1*N*q + (size_t)N*l;
            const size_t base = (size_t)perQ*q + l*l + l;
            for (int m = 0; m <= l; m++){
                double sR = 0, sI = 0;
                const size_t yIdx = (size_t)(l*l + l + m)*N;
                for (unsigned int k = 0; k < N; k++){
                    const double w = electrons[k]*ff*pSBJ[bIdx + k];
                    sR += w*yR[yIdx + k];
                    sI += w*yI[yIdx + k];
                }
                almReal[base + m] += (float)sR;
                almImag[base + m] += (float)sI;
                if (m > 0){   // A_{l,-m} = (-1)^m conj(A_{l,m})
                    if (m & 1){ almReal[base - m] += (float)(-sR); almImag[base - m] += (float)sI; }
                    else      { almReal[base - m] += (float)sR;    almImag[base - m] += (float)(-sI); }
                }
            }
        }
    }
}

void IonCloud::writePseudoAtoms(const std::string & filename, const vector3 * c) const {
    FILE * f = std::fopen(filename.c_str(), "w");
    if (!f) return;
    std::fprintf(f, "REMARK  ion-cloud pseudo-atoms; occupancy column = extra electrons, B column = group charge\n");
    for (size_t k = 0; k < positions.size(); k++){
        const vector3 p = positions[k] + *c;
        std::fprintf(f, "HETATM%5d  ION ION X%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", (int)(k + 1) % 100000, (int)(k + 1) % 10000,
                     p.x, p.y, p.z, std::max(-99.99f, std::min(999.99f, electrons[k])), sites[k].charge);
    }
    std::fprintf(f, "END\n");
    std::fclose(f);
}
