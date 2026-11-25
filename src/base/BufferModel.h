// BufferModel: electron density of the buffer (c_x) from a plain-text buffer file.
//
// File format, one entry per line, '#' starts a comment:
//
//   temperature 25 C            # or K; default 25 C if absent
//   NaCl        137   mM        # <component> <value> <unit>
//   KCl         2.7   mM
//   sucrose     1     %w/v
//   density     1.0081          # optional, measured buffer density in g/cm^3 at the stated temperature
//   myAdditive  20 mM  Z=64 MW=128.1 V=95.0   # custom component: electrons, g/mol, partial molar volume cm^3/mol
//
// Units: M, mM, uM, g/L, mg/mL, %w/v, %v/v (%v/v only for liquids with a tabulated pure density).
//
// c_x is the buffer electron density relative to 0.334 e/A^3, the reference used by the
// excluded-volume form factors.  Without a measured density:
//     rho_e = rho_e,water(T) + sum_i c_i N_A (Z_i - rho_e,water(T) * Vbar_i)
// With a measured density, the electron count per litre is taken directly from the composition
// and the remaining mass (water), which avoids partial-volume errors.
//
// Partial molar volumes are infinite-dilution values near 25 C and are approximate; they are
// used without temperature or concentration corrections.
//
// Ions (for the counterion cloud):
//   pH 7.5                      # default 7.0 if absent
// Every component is split into free ions at the stated pH (salts fully dissociated; weak-acid buffers
// such as HEPES, Tris, MES, MOPS, phosphate by Henderson-Hasselbalch, assuming titration with NaOH / HCl).
// From these follow the ionic strength I = 1/2 sum c_i z_i^2, the Debye length
//     lambda_D = sqrt(eps0 eps_r(T) k T / (2 N_A e^2 I))      (3.04 A / sqrt(I) at 25 C)
// and the electron contrast of each ion against the water it displaces,
//     dN_i = Z_i - rho_w(T) V_i,abs
// with ABSOLUTE single-ion partial molar volumes (conventional values with V(H+) = -5.4 cm^3/mol), e.g.
// Na+ 13.7, K+ 16.0, Mg2+ 27.8, Cl- 5.1 extra electrons.

#ifndef WETSAXS_BUFFERMODEL_H
#define WETSAXS_BUFFERMODEL_H

#include <string>
#include <vector>
#include <map>

struct BufferComponent {
    float electrons;      // Z per formula unit
    float molarMass;      // g/mol
    float partialVolume;  // cm^3/mol (infinite dilution, ~25 C)
    float liquidDensity;  // g/cm^3 of the pure liquid, 0 if not a liquid (needed for %v/v)
};

struct IonSpecies {
    float electrons;      // Z of the ion
    int   charge;         // z
    float volumeConv;     // conventional partial molar volume, cm^3/mol (V(H+) = 0)
};

struct IonEntry {
    std::string name;
    double molar;
    IonSpecies props;
};

struct BufferEntry {
    std::string name;
    double molar;         // mol/L after unit conversion
    BufferComponent props;
};

class BufferModel {

    double temperatureC = 25.0;
    double pH = 7.0;
    bool pHGiven = false;
    std::map<std::string, double> ionMolar;          // free ions, mol/L
    static const std::map<std::string, IonSpecies> & ionTable();
    void speciate(const std::string & key, double molar);
    double measuredDensity = -1.0;   // g/cm^3, < 0 if not given
    std::vector<BufferEntry> entries;
    std::vector<std::string> report;

    static const std::map<std::string, BufferComponent> & table();
    static std::string lower(const std::string & s);
    double toMolar(double value, const std::string & unit, const BufferComponent & props, const std::string & name) const;

public:
    static constexpr double referenceDensity = 0.334;  // e/A^3, reference for c_x

    explicit BufferModel(const std::string & filename);

    static double waterMassDensity(double tempC);       // g/cm^3, Tanaka et al. (2001)
    static double waterElectronDensity(double tempC);   // e/A^3

    double getTemperature() const { return temperatureC; }
    double electronDensity() const;                     // e/A^3 of the buffer
    double getCx() const { return electronDensity()/referenceDensity; }
    const std::vector<std::string> & getReport() const { return report; }

    double getPH() const { return pH; }
    bool isPHGiven() const { return pHGiven; }
    std::vector<IonEntry> getIons() const;
    double ionicStrength() const;                       // mol/L
    double debyeLength() const;                         // Angstrom (infinite if no ions)
    static double waterPermittivity(double tempC);      // eps_r, Malmberg & Maryott (1956)
    double ionContrast(const IonSpecies & ion) const;   // extra electrons vs displaced water
    /*
     * Extra electrons of the counterion atmosphere per unit of neutralised charge of a group with the given
     * sign (-1: anionic group such as a phosphate).  A fraction counterFraction of the neutralisation is an
     * excess of counterions, the rest a depletion of co-ions; within each class the share of species i is
     * c_i z_i^2 / sum c z^2 (linear screening).
     *   E = f sum_counter w_i dN_i/|z_i|  -  (1 - f) sum_co w_j dN_j/|z_j|
     */
    double cloudElectronsPerCharge(int sign, double counterFraction) const;
};

#endif //WETSAXS_BUFFERMODEL_H
