#include "BufferModel.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

namespace {
    constexpr double AVOGADRO = 6.02214076e23;
    constexpr double WATER_MW = 18.01528;   // g/mol
    constexpr double WATER_Z  = 10.0;

    std::string fmt(const char * f, double v){ char b[64]; std::snprintf(b, sizeof(b), f, v); return b; }
}

/*
 * Electrons, molar mass (g/mol), partial molar volume (cm^3/mol, infinite dilution, ~25 C),
 * pure-liquid density (g/cm^3, 0 if not needed).
 * Values marked ~ are approximate; add or override entries with Z= MW= V= in the buffer file.
 * Keys are lower case.
 */
const std::map<std::string, BufferComponent> & BufferModel::table(){
    static const std::map<std::string, BufferComponent> t = {
            // salts
            {"nacl",       {28.f,  58.44f,  16.6f, 0.f}},
            {"kcl",        {36.f,  74.55f,  26.8f, 0.f}},
            {"licl",       {20.f,  42.39f,  17.0f, 0.f}},
            {"nabr",       {46.f, 102.89f,  23.5f, 0.f}},
            {"kbr",        {54.f, 119.00f,  33.7f, 0.f}},
            {"mgcl2",      {46.f,  95.21f,  14.5f, 0.f}},
            {"cacl2",      {54.f, 110.98f,  17.8f, 0.f}},
            {"mgso4",      {58.f, 120.37f,  -7.2f, 0.f}},
            {"rbcl",       {54.f, 120.92f,  31.9f, 0.f}},
            {"cscl",       {72.f, 168.36f,  39.2f, 0.f}},
            {"mncl2",      {61.f, 125.84f,  18.1f, 0.f}},   // ~
            {"nan3",       {32.f,  65.01f,  27.5f, 0.f}},   // ~ sodium azide
            {"sodiumazide",{32.f,  65.01f,  27.5f, 0.f}},   // ~
            {"(nh4)2so4",  {70.f, 132.14f,  49.8f, 0.f}},
            {"nh42so4",    {70.f, 132.14f,  49.8f, 0.f}},
            // phosphates
            {"nah2po4",    {60.f, 119.98f,  27.9f, 0.f}},   // ~
            {"na2hpo4",    {70.f, 141.96f,   3.0f, 0.f}},   // ~
            {"kh2po4",     {68.f, 136.09f,  38.1f, 0.f}},   // ~
            {"k2hpo4",     {86.f, 174.18f,  23.3f, 0.f}},   // ~
            // buffers
            {"tris",       {66.f, 121.14f,  82.6f, 0.f}},
            {"tris-hcl",   {84.f, 157.60f, 101.0f, 0.f}},   // ~
            {"trishcl",    {84.f, 157.60f, 101.0f, 0.f}},   // ~
            {"hepes",      {128.f,238.30f, 168.0f, 0.f}},   // ~
            {"mes",        {104.f,195.24f, 135.0f, 0.f}},   // ~
            {"pipes",      {160.f,302.37f, 190.0f, 0.f}},   // ~
            {"mops",       {112.f,209.26f, 147.0f, 0.f}},   // ~
            {"imidazole",  {36.f,  68.08f,  56.0f, 0.f}},   // ~
            // reducing agents, chelators
            {"dtt",        {82.f, 154.25f, 118.0f, 0.f}},   // ~
            {"tcep",       {150.f,286.65f, 170.0f, 0.f}},   // ~ (TCEP-HCl)
            {"edta",       {174.f,336.21f, 140.0f, 0.f}},   // ~ (Na2H2EDTA)
            // co-solvents, sugars, denaturants
            {"glycerol",   {50.f,  92.09f,  70.9f, 1.261f}},
            {"ethyleneglycol", {34.f, 62.07f, 54.6f, 1.113f}},
            {"dmso",       {42.f,  78.13f,  71.0f, 1.100f}},
            {"sucrose",    {182.f,342.30f, 211.3f, 0.f}},
            {"trehalose",  {182.f,342.30f, 207.0f, 0.f}},   // ~
            {"glucose",    {96.f, 180.16f, 112.0f, 0.f}},
            {"urea",       {32.f,  60.06f,  44.2f, 0.f}},
            {"gdnhcl",     {50.f,  95.53f,  59.0f, 0.f}},   // ~
            {"peg",        {24.f,  44.05f,  36.7f, 0.f}},   // per ethylene-oxide monomer, use %w/v or g/L
    };
    return t;
}

/*
 * Free ions: electrons, charge, conventional partial molar volume (cm^3/mol, 25 C, V(H+) = 0;
 * Millero 1971, Marcus 2011).  Values marked ~ are estimates; ions of minor buffer components.
 */
const std::map<std::string, IonSpecies> & BufferModel::ionTable(){
    static const std::map<std::string, IonSpecies> t = {
            {"Li+",     { 2.f,  1,  -0.88f}},
            {"Na+",     {10.f,  1,  -1.21f}},
            {"K+",      {18.f,  1,   9.02f}},
            {"Rb+",     {36.f,  1,  14.07f}},
            {"Cs+",     {54.f,  1,  21.42f}},
            {"NH4+",    {10.f,  1,  17.86f}},
            {"Mg2+",    {10.f,  2, -21.17f}},
            {"Ca2+",    {18.f,  2, -17.85f}},
            {"Mn2+",    {23.f,  2, -17.7f}},
            {"Cl-",     {18.f, -1,  17.83f}},
            {"Br-",     {36.f, -1,  24.71f}},
            {"N3-",     {22.f, -1,  27.0f}},    // ~
            {"SO4 2-",  {50.f, -2,  13.98f}},
            {"H2PO4-",  {50.f, -1,  29.1f}},
            {"HPO4 2-", {50.f, -2,   5.4f}},    // ~
            {"HEPES-",  {128.f,-1, 165.f}},     // ~
            {"MES-",    {104.f,-1, 130.f}},     // ~
            {"MOPS-",   {112.f,-1, 145.f}},     // ~
            {"PIPES-",  {160.f,-1, 185.f}},     // ~
            {"PIPES 2-",{160.f,-2, 180.f}},     // ~
            {"TrisH+",  {66.f,  1,  83.2f}},
            {"ImH+",    {36.f,  1,  56.f}},     // ~ imidazolium
            {"Gdm+",    {32.f,  1,  41.2f}},
            {"H2EDTA 2-",{154.f,-2, 140.f}},    // ~
    };
    return t;
}

/*
 * Free ions of one component at the buffer pH.  Weak-acid buffers are assumed titrated to the stated pH
 * with NaOH (acids) or HCl (bases).  Non-ionic components (sugars, glycerol, urea, ...) give no ions.
 */
void BufferModel::speciate(const std::string & key, double c){
    auto add = [&](const char * ion, double n){ ionMolar[ion] += n*c; };
    const double T = temperatureC;
    auto acidFraction = [&](double pKa){ return 1.0/(1.0 + std::pow(10.0, pKa - pH)); };   // deprotonated
    auto baseFraction = [&](double pKa){ return 1.0/(1.0 + std::pow(10.0, pH - pKa)); };   // protonated

    if (key == "nacl")       { add("Na+",1); add("Cl-",1); }
    else if (key == "kcl")   { add("K+",1);  add("Cl-",1); }
    else if (key == "licl")  { add("Li+",1); add("Cl-",1); }
    else if (key == "rbcl")  { add("Rb+",1); add("Cl-",1); }
    else if (key == "cscl")  { add("Cs+",1); add("Cl-",1); }
    else if (key == "nabr")  { add("Na+",1); add("Br-",1); }
    else if (key == "kbr")   { add("K+",1);  add("Br-",1); }
    else if (key == "mgcl2") { add("Mg2+",1); add("Cl-",2); }
    else if (key == "cacl2") { add("Ca2+",1); add("Cl-",2); }
    else if (key == "mncl2") { add("Mn2+",1); add("Cl-",2); }
    else if (key == "mgso4") { add("Mg2+",1); add("SO4 2-",1); }
    else if (key == "(nh4)2so4" || key == "nh42so4") { add("NH4+",2); add("SO4 2-",1); }
    else if (key == "gdnhcl") { add("Gdm+",1); add("Cl-",1); }
    else if (key == "nan3" || key == "sodiumazide") { add("Na+",1); add("N3-",1); }
    else if (key == "nah2po4" || key == "na2hpo4" || key == "kh2po4" || key == "k2hpo4"){
        const double a = acidFraction(7.21);
        add("H2PO4-", 1.0 - a); add("HPO4 2-", a);
        add(key[0] == 'n' ? "Na+" : "K+", (key.find("2hpo4") != std::string::npos) ? 2.0 : 1.0);
    }
    else if (key == "hepes") { const double a = acidFraction(7.48 - 0.014*(T - 25.0)); add("HEPES-", a); add("Na+", a); }
    else if (key == "mes")   { const double a = acidFraction(6.10 - 0.011*(T - 25.0)); add("MES-", a);   add("Na+", a); }
    else if (key == "mops")  { const double a = acidFraction(7.20 - 0.013*(T - 25.0)); add("MOPS-", a);  add("Na+", a); }
    else if (key == "pipes") { const double a = acidFraction(6.76 - 0.009*(T - 25.0));
                               add("PIPES-", 1.0 - a); add("PIPES 2-", a); add("Na+", 1.0 + a); }
    else if (key == "tris" || key == "tris-hcl" || key == "trishcl") {
                               const double b = baseFraction(8.07 - 0.028*(T - 25.0)); add("TrisH+", b); add("Cl-", b); }
    else if (key == "imidazole") { const double b = baseFraction(6.95); add("ImH+", b); add("Cl-", b); }
    else if (key == "edta")  { add("H2EDTA 2-", 1); add("Na+", 2); }
    // tcep, dtt, sugars, polyols, urea, peg, dmso: treated as non-ionic (TCEP-HCl is a minor additive)
}

std::string BufferModel::lower(const std::string & s){
    std::string o(s);
    std::transform(o.begin(), o.end(), o.begin(), [](unsigned char c){ return std::tolower(c); });
    return o;
}

double BufferModel::waterMassDensity(double T){
    // Tanaka, Girard, Davis, Peuto & Bignell (2001) Metrologia 38, 301; T in Celsius, result g/cm^3
    const double a1 = -3.983035, a2 = 301.797, a3 = 522528.9, a4 = 69.34881, a5 = 999.974950;
    const double rho = a5*(1.0 - (T + a1)*(T + a1)*(T + a2)/(a3*(T + a4)));
    return rho*1.0e-3;
}

double BufferModel::waterElectronDensity(double T){
    return waterMassDensity(T)*AVOGADRO*WATER_Z/WATER_MW*1.0e-24; // e/A^3
}

double BufferModel::toMolar(double value, const std::string & unitIn, const BufferComponent & p, const std::string & name) const {
    const std::string unit = lower(unitIn);
    if (unit == "m")                     return value;
    if (unit == "mm")                    return value*1.0e-3;
    if (unit == "um" || unit == "µm")    return value*1.0e-6;
    if (unit == "g/l" || unit == "mg/ml") return value/p.molarMass;
    if (unit == "%w/v" || unit == "%")   return value*10.0/p.molarMass;        // g per 100 mL
    if (unit == "%v/v") {
        if (p.liquidDensity <= 0) throw std::runtime_error("buffer: %v/v needs a liquid density for " + name + ", use %w/v or g/L");
        return value*10.0*p.liquidDensity/p.molarMass;
    }
    throw std::runtime_error("buffer: unknown unit '" + unitIn + "' for " + name);
}

BufferModel::BufferModel(const std::string & filename){

    std::ifstream in(filename);
    if (!in) throw std::runtime_error("buffer: cannot open " + filename);

    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)){
        ++lineNo;
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ss(line);
        std::string name;
        if (!(ss >> name)) continue;
        const std::string key = lower(name);

        if (key == "temperature" || key == "temp" || key == "t"){
            double v; std::string u = "C";
            if (!(ss >> v)) throw std::runtime_error("buffer line " + std::to_string(lineNo) + ": temperature value missing");
            ss >> u;
            temperatureC = (lower(u) == "k") ? v - 273.15 : v;
            continue;
        }

        if (key == "ph"){
            if (!(ss >> pH) || pH < 0.0 || pH > 14.0)
                throw std::runtime_error("buffer line " + std::to_string(lineNo) + ": pH value missing or out of range");
            pHGiven = true;
            continue;
        }

        if (key == "density"){
            if (!(ss >> measuredDensity) || measuredDensity <= 0.5 || measuredDensity > 2.0)
                throw std::runtime_error("buffer line " + std::to_string(lineNo) + ": density must be in g/cm^3");
            continue;
        }

        if (key == "water" || key == "h2o") continue;

        double value; std::string unit;
        if (!(ss >> value >> unit))
            throw std::runtime_error("buffer line " + std::to_string(lineNo) + ": expected <component> <value> <unit>");

        BufferComponent props{0.f, 0.f, 0.f, 0.f};
        bool found = false;
        const auto it = table().find(key);
        if (it != table().end()){ props = it->second; found = true; }

        // optional custom properties Z= MW= V= rho=
        std::string tok; bool custom = false;
        while (ss >> tok){
            const auto eq = tok.find('=');
            if (eq == std::string::npos) continue;
            const std::string k = lower(tok.substr(0, eq));
            const float v = std::stof(tok.substr(eq + 1));
            if (k == "z") props.electrons = v;
            else if (k == "mw") props.molarMass = v;
            else if (k == "v") props.partialVolume = v;
            else if (k == "rho") props.liquidDensity = v;
            custom = true;
        }

        if (!found && !(custom && props.electrons > 0 && props.molarMass > 0))
            throw std::runtime_error("buffer line " + std::to_string(lineNo) + ": unknown component '" + name +
                                     "', give Z=<electrons> MW=<g/mol> V=<cm^3/mol>");

        entries.push_back({name, toMolar(value, unit, props, name), props});
    }

    for (const auto & e : entries) speciate(lower(e.name), e.molar);

    // report
    report.push_back("Buffer file :: " + filename);
    report.push_back("Temperature (C) :: " + fmt("%.1f", temperatureC));
    report.push_back("Water electron density (e/A^3) :: " + fmt("%.5f", waterElectronDensity(temperatureC)));
    const double rhoW = waterElectronDensity(temperatureC);
    for (const auto & e : entries){
        const double perMolecule = e.props.electrons - rhoW*e.props.partialVolume*1.0e24/AVOGADRO;
        const double dRho = e.molar*AVOGADRO*1.0e-27*perMolecule;
        report.push_back("  " + e.name + " :: " + fmt("%.4g M", e.molar) + "  d(c_x) " + fmt("%+.5f", dRho/referenceDensity));
    }
    if (measuredDensity > 0)
        report.push_back("Measured density (g/cm^3) :: " + fmt("%.5f", measuredDensity) + " (used for c_x)");
    report.push_back("Buffer electron density (e/A^3) :: " + fmt("%.5f", electronDensity()));
    report.push_back("c_x (relative to 0.334 e/A^3) :: " + fmt("%.4f", getCx()));
    report.push_back("pH :: " + fmt("%.2f", pH) + (pHGiven ? "" : " (default, not given in buffer file)"));
    for (const auto & ion : getIons())
        report.push_back("  ion " + ion.name + " :: " + fmt("%.4g M", ion.molar) + "  extra electrons " + fmt("%+.1f", ionContrast(ion.props)));
    report.push_back("Ionic strength (M) :: " + fmt("%.4f", ionicStrength()));
    if (ionicStrength() > 0){
        report.push_back("Debye length (A) :: " + fmt("%.2f", debyeLength()));
        report.push_back("Counterion cloud, e per charge (anionic group, f = 0.71) :: " + fmt("%.2f", cloudElectronsPerCharge(-1, 0.71)));
        report.push_back("Counterion cloud, e per charge (protein, f = 0.50: -, +) :: " + fmt("%.2f", cloudElectronsPerCharge(-1, 0.5))
                         + ", " + fmt("%.2f", cloudElectronsPerCharge(+1, 0.5)));
    }
}

std::vector<IonEntry> BufferModel::getIons() const {
    std::vector<IonEntry> out;
    for (const auto & kv : ionMolar){
        if (kv.second <= 0.0) continue;
        out.push_back({kv.first, kv.second, ionTable().at(kv.first)});
    }
    return out;
}

double BufferModel::ionicStrength() const {
    double I = 0.0;
    for (const auto & ion : getIons()) I += 0.5*ion.molar*ion.props.charge*ion.props.charge;
    return I;
}

double BufferModel::waterPermittivity(double T){
    // Malmberg & Maryott (1956) J. Res. NBS 56, 1; T in Celsius
    return 87.740 - 0.40008*T + 9.398e-4*T*T - 1.410e-6*T*T*T;
}

double BufferModel::debyeLength() const {
    const double I = ionicStrength();
    if (I <= 0.0) return std::numeric_limits<double>::infinity();
    constexpr double eps0 = 8.8541878128e-12, kB = 1.380649e-23, e = 1.602176634e-19;
    const double T = temperatureC + 273.15;
    const double lam = std::sqrt(eps0*waterPermittivity(temperatureC)*kB*T/(2.0*AVOGADRO*e*e*I*1000.0)); // m
    return lam*1.0e10;
}

double BufferModel::ionContrast(const IonSpecies & ion) const {
    constexpr double VH = -5.4;                                     // absolute V(H+), cm^3/mol
    const double vAbs = ion.volumeConv + ion.charge*VH;             // cm^3/mol
    return ion.electrons - waterElectronDensity(temperatureC)*vAbs*1.0e24/AVOGADRO;
}

double BufferModel::cloudElectronsPerCharge(int sign, double f) const {
    double sc = 0.0, so = 0.0, wc = 0.0, wo = 0.0;   // weighted sums and normalisation, counter / co-ions
    for (const auto & ion : getIons()){
        const int z = ion.props.charge;
        const double w = ion.molar*z*z;
        const double perCharge = ionContrast(ion.props)/std::abs(z);
        if (z*sign < 0) { sc += w*perCharge; wc += w; }
        else            { so += w*perCharge; wo += w; }
    }
    const double counter = (wc > 0) ? sc/wc : 0.0;
    const double co      = (wo > 0) ? so/wo : 0.0;
    return f*counter - (1.0 - f)*co;
}

double BufferModel::electronDensity() const {

    const double rhoW = waterElectronDensity(temperatureC);

    if (measuredDensity > 0){
        // electrons per litre: solutes + the remaining mass as water
        double soluteMass = 0.0, soluteElectrons = 0.0;          // per litre
        for (const auto & e : entries){
            soluteMass      += e.molar*e.props.molarMass;
            soluteElectrons += e.molar*e.props.electrons;
        }
        const double waterMass = measuredDensity*1000.0 - soluteMass; // g per litre
        const double electronsPerLitre = soluteElectrons + waterMass/WATER_MW*WATER_Z; // mol of electrons per litre
        return electronsPerLitre*AVOGADRO*1.0e-27;
    }

    double rho = rhoW;
    for (const auto & e : entries){
        const double perMolecule = e.props.electrons - rhoW*e.props.partialVolume*1.0e24/AVOGADRO;
        rho += e.molar*AVOGADRO*1.0e-27*perMolecule;
    }
    return rho;
}
