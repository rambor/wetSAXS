// Copyright (c) 2022.
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//
// Created by Robert Rambo on 04/05/2022.
//

#include <random>
#include "Waters.h"
#include "sastools/Residues.h"

Waters::Waters() {
    /*
     * closest observed waters in Angstroms (ignoring hydrogens)
     */
    minima["ALA"] = 2.85256f;
    minima["ARG"] = 2.59187f;
    minima["ASN"] = 2.75131f;
    minima["ASP"] = 2.42051f;
    minima["CYS"] = 2.69287f;
    minima["DA"] = 2.74323f; // DNA
    minima["DC"] = 3.26352f; // DNA
    minima["DG"] = 2.6976f;  // DNA
    minima["DT"] = 3.17513f; // DNA
    minima["GLN"] = 2.90458f;
    minima["GLU"] = 2.35418f;
    minima["GLY"] = 2.69029f;
    minima["HIS"] = 2.58119f;
    minima["ILE"] = 2.81091f;
    minima["LEU"] = 2.9619f;
    minima["LYS"] = 2.17605f;
    minima["MET"] = 3.29484f;
    minima["PHE"] = 2.95653f;
    minima["PRO"] = 3.38124f;
    minima["SER"] = 2.5156f;
    minima["THR"] = 2.43802f;
    minima["TRP"] = 2.93771f;
    minima["TYR"] = 2.45725f;
    minima["VAL"] = 2.78789f;
    minima["rA"] = 2.14787f; // RNA
    minima["rC"] = 2.73046f; // RNA
    minima["rG"] = 2.83705f; // RNA
    minima["rU"] = 2.24522f; // RNA

    gen = std::mt19937(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    randomIndex = std::uniform_int_distribution<int> (0,360); // guaranteed unbiased
    randomBeta = std::uniform_int_distribution<int> (0,180); // guaranteed unbiased
    createBackbonesAndSideChains();
}

/*
 * no terminal OXT is included in the library which means that OXT will never be hydrated unless we make a special generic model
 *
 */
void Waters::createBackbonesAndSideChains(){
    // need number of atoms per backbone type
    // create protein side chains
    // specify only the atoms used in the alignment
    // total models used :: 257
    centeringAtoms["PROTEIN"] = "CA";  // use this for backbone water placements
    sideChains["PROTEIN"].insert({"N",Coords(1.1187,0.9331,-0.1188,"N", 1.0)});
    //sideChains["PROTEIN"].insert({"CA",Coords(0.0000,0.0000,0.0000,"CA", 1.0)});
    sideChains["PROTEIN"].insert({"C",Coords(0.1553,-1.0269,-1.1186,"C", 1.0)});
    sideChains["PROTEIN"].insert({"O",Coords(0.7187,-0.7198,-2.1729,"O", 1.0)});

    centeringAtoms["ALA"] = "CA";
    sideChains["ALA"].insert({"N",Coords(0.0934,-1.3677,0.4904,"N", 1.0)});
//    sideChains["ALA"].insert({"CA",Coords(0.0000,0.0000,0.0000,"CA", 1.0)});
    sideChains["ALA"].insert({"C",Coords(0.4932,0.1092,-1.4387,"C", 1.0)});
//    sideChains["ALA"].insert({"O",Coords(-0.1131,0.8131,-2.2526,"O", 1.0)});
    sideChains["ALA"].insert({"CB",Coords(0.7925,0.9368,0.9100,"CB", 1.0)});

    centeringAtoms["ARG"] = "NE";
    // sideChains["ARG"].insert({"N",Coords(3.2083,-4.6138,-2.7459,"N", 1.0)});
    // sideChains["ARG"].insert({"CA",Coords(3.2915,-3.4222,-1.9157,"CA", 1.0)});
    // sideChains["ARG"].insert({"C",Coords(3.4279,-3.7978,-0.4377,"C", 1.0)});
    // sideChains["ARG"].insert({"O",Coords(3.5730,-2.9236,0.3994,"O", 1.0)});
    // sideChains["ARG"].insert({"CB",Coords(2.0374,-2.5445,-2.0272,"CB", 1.0)});
    // sideChains["ARG"].insert({"CG",Coords(1.8765,-1.3183,-1.1369,"CG", 1.0)});
    sideChains["ARG"].insert({"CD",Coords(0.4452,-1.2514,-0.5805,"CD", 1.0)});
//    sideChains["ARG"].insert({"NE",Coords(0.0000,0.0000,0.0000,"NE", 1.0)});
    sideChains["ARG"].insert({"CZ",Coords(-1.1522,0.5704,-0.3704,"CZ", 1.0)});
    // sideChains["ARG"].insert({"NH1",Coords(-1.9218,0.0082,-1.2996,"NH1", 1.0)});
    // sideChains["ARG"].insert({"NH2",Coords(-1.5062,1.7393,0.2166,"NH2", 1.0)});

    centeringAtoms["ASN"] = "CG";
    // sideChains["ASN"].insert({"N",Coords(-2.0986,2.0972,1.0121,"N", 1.0)});
    // sideChains["ASN"].insert({"CA",Coords(-0.7692,1.8155,1.5366,"CA", 1.0)});
    // sideChains["ASN"].insert({"C",Coords(-0.1819,2.9533,2.3462,"C", 1.0)});
    // sideChains["ASN"].insert({"O",Coords(1.0052,2.8983,2.7000,"O", 1.0)});
    sideChains["ASN"].insert({"CB",Coords(0.1740,1.4294,0.4037,"CB", 1.0)});
//    sideChains["ASN"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["ASN"].insert({"OD1",Coords(-0.1975,-0.8648,0.8475,"OD1", 1.0)});
    sideChains["ASN"].insert({"ND2",Coords(0.0308,-0.2644,-1.2917,"ND2", 1.0)});

    centeringAtoms["ASP"] = "CG";
    // sideChains["ASP"].insert({"N",Coords(-1.2444,-2.7299,1.5292,"N", 1.0)});
    // sideChains["ASP"].insert({"CA",Coords(0.1049,-2.2248,1.2979,"CA", 1.0)});
    // sideChains["ASP"].insert({"C",Coords(1.0277,-2.7905,2.3707,"C", 1.0)});
    // sideChains["ASP"].insert({"O",Coords(0.5786,-3.1001,3.4954,"O", 1.0)});
    sideChains["ASP"].insert({"CB",Coords(0.2153,-0.6890,1.3211,"CB", 1.0)});
    // sideChains["ASP"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["ASP"].insert({"OD1",Coords(0.5809,-0.4630,-1.0191,"OD1", 1.0)});
    sideChains["ASP"].insert({"OD2",Coords(-0.7002,1.0317,0.0065,"OD2", 1.0)});


    centeringAtoms["CYS"] = "CB";
    sideChains["CYS"].insert({"N",Coords(1.5321,-1.7611,-0.7573,"N", 1.0)});
    sideChains["CYS"].insert({"CA",Coords(0.2722,-1.4998,-0.0708,"CA", 1.0)});
    // sideChains["CYS"].insert({"C",Coords(-0.9243,-2.1641,-0.7358,"C", 1.0)});
    // sideChains["CYS"].insert({"O",Coords(-2.0315,-2.1004,-0.1902,"O", 1.0)});
    //sideChains["CYS"].insert({"CB",Coords(0.0000,0.0000,0.0000,"CB", 1.0)});
    sideChains["CYS"].insert({"SG",Coords(1.0822,0.8001,1.1822,"SG", 1.0)});


    centeringAtoms["GLN"] = "CD";
    // sideChains["GLN"].insert({"N",Coords(-2.6467,-2.2398,1.7537,"N", 1.0)});
    // sideChains["GLN"].insert({"CA",Coords(-1.2391,-2.4590,1.4061,"CA", 1.0)});
    // sideChains["GLN"].insert({"C",Coords(-0.7808,-3.8851,1.7245,"C", 1.0)});
    // sideChains["GLN"].insert({"O",Coords(0.1971,-4.1298,2.4386,"O", 1.0)});
    sideChains["GLN"].insert({"CB",Coords(-0.3600,-1.3967,2.0682,"CB", 1.0)});
    sideChains["GLN"].insert({"CG",Coords(-0.5534,-0.0286,1.4149,"CG", 1.0)});
//    sideChains["GLN"].insert({"CD",Coords(0.0000,0.0000,0.0000,"CD", 1.0)});
    sideChains["GLN"].insert({"OE1",Coords(1.1969,0.1844,-0.1943,"OE1", 1.0)});
    sideChains["GLN"].insert({"NE2",Coords(-0.8625,-0.1958,-0.9918,"NE2", 1.0)});


    centeringAtoms["GLU"] = "CG";
    // sideChains["GLU"].insert({"N",Coords(2.0942,2.7098,-1.7139,"N", 1.0)});
    // sideChains["GLU"].insert({"CA",Coords(0.9054,2.0917,-1.1372,"CA", 1.0)});
    // sideChains["GLU"].insert({"C",Coords(0.4140,2.9160,0.0413,"C", 1.0)});
    // sideChains["GLU"].insert({"O",Coords(-0.7699,3.2412,0.1406,"O", 1.0)});
    sideChains["GLU"].insert({"CB",Coords(1.1965,0.6662,-0.6674,"CB", 1.0)});
//    sideChains["GLU"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["GLU"].insert({"CD",Coords(0.3832,-1.1904,0.8559,"CD", 1.0)});
    sideChains["GLU"].insert({"OE1",Coords(1.5954,-1.4469,1.0136,"OE1", 1.0)});
    sideChains["GLU"].insert({"OE2",Coords(-0.5298,-1.8671,1.3755,"OE2", 1.0)});


    // need GLY
    centeringAtoms["GLY"] = "CA";
    sideChains["GLY"].insert({"N",Coords(1.1611,-0.4046,-0.7683,"N", 1.0)});
//    sideChains["GLY"].insert({"CA",Coords(0.0000,0.0000,0.0000,"CA", 1.0)});
    sideChains["GLY"].insert({"C",Coords(-1.2843,-0.0080,-0.8130,"C", 1.0)});
    sideChains["GLY"].insert({"O",Coords(-2.3491,-0.3576,-0.3003,"O", 1.0)});


    centeringAtoms["HIS"] = "CG";
    // sideChains["HIS"].insert({"N",Coords(-1.5829,2.5757,0.4693,"N", 1.0)});
    // sideChains["HIS"].insert({"CA",Coords(-0.1697,2.5467,0.1347,"CA", 1.0)});
    // sideChains["HIS"].insert({"C",Coords(0.2410,3.8496,-0.5515,"C", 1.0)});
    // sideChains["HIS"].insert({"O",Coords(1.0592,3.8588,-1.4697,"O", 1.0)});
    sideChains["HIS"].insert({"CB",Coords(0.1702,1.3167,-0.7131,"CB", 1.0)});
    //sideChains["HIS"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["HIS"].insert({"ND1",Coords(0.2003,-1.2070,-0.6361,"ND1", 1.0)});
    sideChains["HIS"].insert({"CD2",Coords(-0.3740,-0.3075,1.2675,"CD2", 1.0)});
    sideChains["HIS"].insert({"CE1",Coords(-0.0282,-2.1962,0.2080,"CE1", 1.0)});
    sideChains["HIS"].insert({"NE2",Coords(-0.3845,-1.6804,1.3699,"NE2", 1.0)});

    centeringAtoms["ILE"] = "CB";
    // sideChains["ILE"].insert({"N",Coords(-1.6383,0.1718,-1.8553,"N", 1.0)});
    // sideChains["ILE"].insert({"CA",Coords(-0.6073,-0.6638,-1.2511,"CA", 1.0)});
    // sideChains["ILE"].insert({"C",Coords(-1.1846,-2.0365,-0.9070,"C", 1.0)});
    // sideChains["ILE"].insert({"O",Coords(-0.6128,-3.0759,-1.2403,"O", 1.0)});
    //sideChains["ILE"].insert({"CB",Coords(0.0000,0.0000,0.0000,"CB", 1.0)});
    sideChains["ILE"].insert({"CG1",Coords(0.7479,1.2920,-0.3700,"CG1", 1.0)});
    sideChains["ILE"].insert({"CG2",Coords(0.9165,-0.9780,0.7215,"CG2", 1.0)});
    sideChains["ILE"].insert({"CD1",Coords(1.1515,2.1283,0.8242,"CD1", 1.0)});

    // redo LEU and GLY

    centeringAtoms["LEU"] = "CG"; // needs to be fixed
    // sideChains["LEU"].insert({"N",Coords(0.8701,-2.8051,0.7512,"N", 1.0)});
    // sideChains["LEU"].insert({"CA",Coords(-0.5093,-2.5295,0.3813,"CA", 1.0)});
    // sideChains["LEU"].insert({"C",Coords(-1.1358,-3.7171,-0.3471,"C", 1.0)});
    // sideChains["LEU"].insert({"O",Coords(-0.4388,-4.5856,-0.8956,"O", 1.0)});
    sideChains["LEU"].insert({"CB",Coords(-0.5799,-1.3089,-0.5401,"CB", 1.0)});
    //sideChains["LEU"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["LEU"].insert({"CD1",Coords(-0.1551,1.0985,-1.0440,"CD1", 1.0)});
    sideChains["LEU"].insert({"CD2",Coords(-0.6712,0.3944,1.3053,"CD2", 1.0)});


    centeringAtoms["LYS"] = "CD";
    // sideChains["LYS"].insert({"N",Coords(0.5853,-2.4583,-2.5890,"N", 1.0)});
    // sideChains["LYS"].insert({"CA",Coords(0.0654,-2.8255,-1.2860,"CA", 1.0)});
    // sideChains["LYS"].insert({"C",Coords(-0.0970,-4.3208,-1.1432,"C", 1.0)});
    // sideChains["LYS"].insert({"O",Coords(-0.9530,-4.7966,-0.3811,"O", 1.0)});
    sideChains["LYS"].insert({"CB",Coords(-1.1839,-2.0317,-0.9070,"CB", 1.0)});
    sideChains["LYS"].insert({"CG",Coords(-1.0242,-0.5329,-0.9685,"CG", 1.0)});
    //sideChains["LYS"].insert({"CD",Coords(0.0000,0.0000,0.0000,"CD", 1.0)});
    sideChains["LYS"].insert({"CE",Coords(-0.0052,1.4982,0.1670,"CE", 1.0)});
    sideChains["LYS"].insert({"NZ",Coords(0.2681,2.2209,-1.0662,"NZ", 1.0)});


    centeringAtoms["MET"] = "SD";
    // sideChains["MET"].insert({"N",Coords(3.4249,2.9208,1.1906,"N", 1.0)});
    // sideChains["MET"].insert({"CA",Coords(2.5418,2.3814,2.2257,"CA", 1.0)});
    // sideChains["MET"].insert({"C",Coords(3.2950,2.0850,3.5198,"C", 1.0)});
    // sideChains["MET"].insert({"O",Coords(2.7732,2.3184,4.6224,"O", 1.0)});
    sideChains["MET"].insert({"CB",Coords(1.8683,1.1312,1.6668,"CB", 1.0)});
    sideChains["MET"].insert({"CG",Coords(0.8341,1.4690,0.5953,"CG", 1.0)});
//    sideChains["MET"].insert({"SD",Coords(0.0000,0.0000,0.0000,"SD", 1.0)});
    sideChains["MET"].insert({"CE",Coords(-1.6309,0.6898,-0.4097,"CE", 1.0)});

    centeringAtoms["PHE"] = "CG";
    // sideChains["PHE"].insert({"N",Coords(-2.0971,0.1681,3.1625,"N", 1.0)});
    // sideChains["PHE"].insert({"CA",Coords(-1.1948,0.6495,2.1272,"CA", 1.0)});
    // sideChains["PHE"].insert({"C",Coords(-1.8703,1.7271,1.2819,"C", 1.0)});
    // sideChains["PHE"].insert({"O",Coords(-1.2926,2.7869,1.0485,"O", 1.0)});
    sideChains["PHE"].insert({"CB",Coords(-0.7024,-0.4881,1.2289,"CB", 1.0)});
    //sideChains["PHE"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["PHE"].insert({"CD1",Coords(-0.6891,0.1344,-1.1926,"CD1", 1.0)});
    sideChains["PHE"].insert({"CD2",Coords(1.3390,0.3661,0.0420,"CD2", 1.0)});
    sideChains["PHE"].insert({"CE1",Coords(-0.0572,0.6086,-2.3268,"CE1", 1.0)});
    sideChains["PHE"].insert({"CE2",Coords(1.9768,0.8363,-1.0932,"CE2", 1.0)});
    sideChains["PHE"].insert({"CZ",Coords(1.2759,0.9601,-2.2788,"CZ", 1.0)});

    centeringAtoms["PRO"] = "CG";
    // sideChains["PRO"].insert({"N",Coords(-1.4526,0.7688,1.6635,"N", 1.0)});
    sideChains["PRO"].insert({"CA",Coords(-2.1501,0.9466,0.3866,"CA", 1.0)});
    // sideChains["PRO"].insert({"C",Coords(-3.6355,0.5827,0.4349,"C", 1.0)});
    // sideChains["PRO"].insert({"O",Coords(-4.4327,1.1538,-0.3158,"O", 1.0)});
    sideChains["PRO"].insert({"CB",Coords(-1.3889,0.0068,-0.5542,"CB", 1.0)});
    sideChains["PRO"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["PRO"].insert({"CD",Coords(-0.2013,0.0027,1.4968,"CD", 1.0)});


    centeringAtoms["SER"] = "CB";
    // sideChains["SER"].insert({"N",Coords(1.5855,-0.5545,1.7881,"N", 1.0)});
    sideChains["SER"].insert({"CA",Coords(0.2598,-0.0332,1.5035,"CA", 1.0)});
    // sideChains["SER"].insert({"C",Coords(-0.7981,-0.8978,2.1694,"C", 1.0)});
    // sideChains["SER"].insert({"O",Coords(-0.5383,-2.0310,2.5868,"O", 1.0)});
    sideChains["SER"].insert({"CB",Coords(0.0000,0.0000,0.0000,"CB", 1.0)});
    sideChains["SER"].insert({"OG",Coords(0.0341,-1.3122,-0.5297,"OG", 1.0)});


    centeringAtoms["THR"] = "CB";
    // sideChains["THR"].insert({"N",Coords(-1.1494,-1.8000,-1.1772,"N", 1.0)});
    sideChains["THR"].insert({"CA",Coords(-0.4695,-1.4639,0.0585,"CA", 1.0)});
    // sideChains["THR"].insert({"C",Coords(-1.4715,-1.6589,1.1886,"C", 1.0)});
    // sideChains["THR"].insert({"O",Coords(-2.5986,-1.1507,1.0951,"O", 1.0)});
    sideChains["THR"].insert({"CB",Coords(0.0000,0.0000,0.0000,"CB", 1.0)});
    sideChains["THR"].insert({"OG1",Coords(0.7438,0.2066,-1.2093,"OG1", 1.0)});
    sideChains["THR"].insert({"CG2",Coords(0.8668,0.3498,1.1868,"CG2", 1.0)});


    centeringAtoms["TRP"] = "CD2";
    // sideChains["TRP"].insert({"N",Coords(3.2349,-2.3998,-2.4844,"N", 1.0)});
    // sideChains["TRP"].insert({"CA",Coords(2.4877,-2.0593,-1.2736,"CA", 1.0)});
    // sideChains["TRP"].insert({"C",Coords(1.1415,-2.7784,-1.2450,"C", 1.0)});
    // sideChains["TRP"].insert({"O",Coords(0.7417,-3.3528,-0.2213,"O", 1.0)});
    sideChains["TRP"].insert({"CB",Coords(2.2914,-0.5370,-1.1779,"CB", 1.0)});
    sideChains["TRP"].insert({"CG",Coords(1.4216,-0.1501,-0.0084,"CG", 1.0)});
    sideChains["TRP"].insert({"CD1",Coords(1.8352,0.1124,1.2656,"CD1", 1.0)});
//    sideChains["TRP"].insert({"CD2",Coords(0.0000,0.0000,0.0000,"CD2", 1.0)});
    sideChains["TRP"].insert({"NE1",Coords(0.7605,0.3918,2.0699,"NE1", 1.0)});
    sideChains["TRP"].insert({"CE2",Coords(-0.3781,0.3420,1.3146,"CE2", 1.0)});
    sideChains["TRP"].insert({"CE3",Coords(-0.9929,-0.1115,-0.9752,"CE3", 1.0)});
    sideChains["TRP"].insert({"CZ2",Coords(-1.6962,0.5756,1.6744,"CZ2", 1.0)});
    sideChains["TRP"].insert({"CZ3",Coords(-2.2991,0.1131,-0.6109,"CZ3", 1.0)});
    sideChains["TRP"].insert({"CH2",Coords(-2.6404,0.4550,0.7008,"CH2", 1.0)});


    centeringAtoms["TYR"] = "CG";
    // sideChains["TYR"].insert({"N",Coords(0.9805,-0.2460,2.8742,"N", 1.0)});
    // sideChains["TYR"].insert({"CA",Coords(1.1436,-1.2996,1.8757,"CA", 1.0)});
    // sideChains["TYR"].insert({"C",Coords(2.2857,-2.2019,2.1686,"C", 1.0)});
    // sideChains["TYR"].insert({"O",Coords(2.2109,-3.4114,1.8720,"O", 1.0)});
    sideChains["TYR"].insert({"CB",Coords(1.2498,-0.7285,0.4482,"CB", 1.0)});
//    sideChains["TYR"].insert({"CG",Coords(0.0000,0.0000,0.0000,"CG", 1.0)});
    sideChains["TYR"].insert({"CD1",Coords(-1.1086,-0.7560,-0.3483,"CD1", 1.0)});
    sideChains["TYR"].insert({"CD2",Coords(-0.0991,1.3498,-0.0528,"CD2", 1.0)});
    sideChains["TYR"].insert({"CE1",Coords(-2.2987,-0.2072,-0.7579,"CE1", 1.0)});
    sideChains["TYR"].insert({"CE2",Coords(-1.2567,1.9707,-0.4466,"CE2", 1.0)});
    sideChains["TYR"].insert({"CZ",Coords(-2.3402,1.1806,-0.8007,"CZ", 1.0)});
    sideChains["TYR"].insert({"OH",Coords(-3.4845,1.8182,-1.1980,"OH", 1.0)});

    centeringAtoms["VAL"] = "CB";
//    sideChains["VAL"].insert({"N",Coords(-1.3599,-0.3408,-2.0525,"N", 1.0)});
    sideChains["VAL"].insert({"CA",Coords(-0.8605,-0.9698,-0.8309,"CA", 1.0)});
//    sideChains["VAL"].insert({"C",Coords(-2.0322,-1.5009,-0.0069,"C", 1.0)});
//    sideChains["VAL"].insert({"O",Coords(-3.1992,-1.2106,-0.3005,"O", 1.0)});
//    sideChains["VAL"].insert({"CB",Coords(0.0000,0.0000,0.0000,"CB", 1.0)});
    sideChains["VAL"].insert({"CG1",Coords(1.1661,0.5228,-0.8471,"CG1", 1.0)});
    sideChains["VAL"].insert({"CG2",Coords(-0.8535,1.1470,0.5450,"CG2", 1.0)});


    /*
     * add nucleics - center on C1'
     *
     */
// total models used :: 52
     centeringAtoms["DA"] = "C1'";
    // sideChains["DA"].insert({"P",Coords(-3.7645,-1.0688,2.6759,"P", 1.0)});
    // sideChains["DA"].insert({"OP1",Coords(-4.5541,-2.3053,2.6174,"OP1", 1.0)});
    // sideChains["DA"].insert({"OP2",Coords(-3.8643,-0.1510,3.8605,"OP2", 1.0)});
    // sideChains["DA"].insert({"O5'",Coords(-2.2319,-1.4166,2.4782,"O5'", 1.0)});
    // sideChains["DA"].insert({"C5'",Coords(-1.8873,-2.4875,1.5795,"C5'", 1.0)});
    // sideChains["DA"].insert({"C4'",Coords(-0.6650,-2.0915,0.7937,"C4'", 1.0)});
    // sideChains["DA"].insert({"O4'",Coords(-0.9752,-1.0036,-0.1444,"O4'", 1.0)});
    // sideChains["DA"].insert({"C3'",Coords(0.4982,-1.5823,1.6379,"C3'", 1.0)});
    // sideChains["DA"].insert({"O3'",Coords(1.7263,-2.1416,1.1642,"O3'", 1.0)});
    // sideChains["DA"].insert({"C2'",Coords(0.3795,-0.0732,1.4974,"C2'", 1.0)});
    // sideChains["DA"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DA"].insert({"N9",Coords(-0.6130,1.2251,-0.4790,"N9", 1.0)});
    sideChains["DA"].insert({"C8",Coords(-1.4636,2.0587,0.1692,"C8", 1.0)});
    sideChains["DA"].insert({"N7",Coords(-1.8588,3.0910,-0.5103,"N7", 1.0)});
    sideChains["DA"].insert({"C5",Coords(-1.1999,2.9178,-1.7335,"C5", 1.0)});
    sideChains["DA"].insert({"C6",Coords(-1.2099,3.6953,-2.8941,"C6", 1.0)});
    sideChains["DA"].insert({"N6",Coords(-1.8913,4.8431,-3.0995,"N6", 1.0)});
    sideChains["DA"].insert({"N1",Coords(-0.4449,3.2480,-3.9149,"N1", 1.0)});
    sideChains["DA"].insert({"C2",Coords(0.2570,2.0875,-3.7503,"C2", 1.0)});
    sideChains["DA"].insert({"N3",Coords(0.3421,1.2592,-2.6870,"N3", 1.0)});
    sideChains["DA"].insert({"C4",Coords(-0.4379,1.7823,-1.7271,"C4", 1.0)});

    // total models used :: 112
    centeringAtoms["DC"] = "C1'";
    // sideChains["DC"].insert({"P",Coords(-0.4474,4.6050,-2.2079,"P", 1.0)});
    // sideChains["DC"].insert({"OP1",Coords(-1.6757,4.9981,-1.5240,"OP1", 1.0)});
    // sideChains["DC"].insert({"OP2",Coords(0.6569,5.5484,-2.3068,"OP2", 1.0)});
    // sideChains["DC"].insert({"O5'",Coords(0.1175,3.2648,-1.4968,"O5'", 1.0)});
    // sideChains["DC"].insert({"C5'",Coords(0.7337,3.3538,-0.1826,"C5'", 1.0)});
    // sideChains["DC"].insert({"C4'",Coords(1.2006,1.9961,0.2612,"C4'", 1.0)});
    // sideChains["DC"].insert({"O4'",Coords(0.0168,1.2651,0.6267,"O4'", 1.0)});
    // sideChains["DC"].insert({"C3'",Coords(1.9026,1.1248,-0.7996,"C3'", 1.0)});
    // sideChains["DC"].insert({"O3'",Coords(2.8917,0.3515,-0.1480,"O3'", 1.0)});
    // sideChains["DC"].insert({"C2'",Coords(0.7826,0.2045,-1.3015,"C2'", 1.0)});
    // sideChains["DC"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DC"].insert({"N1",Coords(-1.4216,-0.3703,-0.2564,"N1", 1.0)});
    sideChains["DC"].insert({"C2",Coords(-1.7531,-1.7115,-0.3550,"C2", 1.0)});
    sideChains["DC"].insert({"O2",Coords(-0.8577,-2.5387,-0.1646,"O2", 1.0)});
    sideChains["DC"].insert({"N3",Coords(-3.0360,-2.0414,-0.6024,"N3", 1.0)});
    sideChains["DC"].insert({"C4",Coords(-3.9842,-1.1211,-0.8061,"C4", 1.0)});
    sideChains["DC"].insert({"N4",Coords(-5.2294,-1.5482,-1.0929,"N4", 1.0)});
    sideChains["DC"].insert({"C5",Coords(-3.6604,0.2620,-0.7713,"C5", 1.0)});
    sideChains["DC"].insert({"C6",Coords(-2.3863,0.5664,-0.4736,"C6", 1.0)});

    // total models used :: 143
    centeringAtoms["DG"] = "C1'";
    // sideChains["DG"].insert({"P",Coords(-2.5379,0.5634,4.9725,"P", 1.0)});
    // sideChains["DG"].insert({"OP1",Coords(-1.5156,0.7990,6.0221,"OP1", 1.0)});
    // sideChains["DG"].insert({"OP2",Coords(-3.8839,1.1745,5.1092,"OP2", 1.0)});
    // sideChains["DG"].insert({"O5'",Coords(-1.9197,1.0194,3.5780,"O5'", 1.0)});
    // sideChains["DG"].insert({"C5'",Coords(-2.5930,0.7347,2.3523,"C5'", 1.0)});
    // sideChains["DG"].insert({"C4'",Coords(-1.7128,1.0956,1.1791,"C4'", 1.0)});
    // sideChains["DG"].insert({"O4'",Coords(-0.4583,0.3936,1.2818,"O4'", 1.0)});
    // sideChains["DG"].insert({"C3'",Coords(-2.2949,0.7558,-0.1847,"C3'", 1.0)});
    // sideChains["DG"].insert({"O3'",Coords(-2.9081,1.9608,-0.6620,"O3'", 1.0)});
    // sideChains["DG"].insert({"C2'",Coords(-1.0785,0.3764,-1.0147,"C2'", 1.0)});
    // sideChains["DG"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DG"].insert({"N9",Coords(0.3248,-1.4211,0.0702,"N9", 1.0)});
    sideChains["DG"].insert({"C8",Coords(1.5672,-1.9935,-0.0626,"C8", 1.0)});
    sideChains["DG"].insert({"N7",Coords(1.5506,-3.2915,0.0729,"N7", 1.0)});
    sideChains["DG"].insert({"C5",Coords(0.2159,-3.5957,0.3048,"C5", 1.0)});
    sideChains["DG"].insert({"C6",Coords(-0.4169,-4.8449,0.5345,"C6", 1.0)});
    sideChains["DG"].insert({"O6",Coords(0.0944,-5.9709,0.5815,"O6", 1.0)});
    sideChains["DG"].insert({"N1",Coords(-1.7870,-4.6969,0.7218,"N1", 1.0)});
    sideChains["DG"].insert({"C2",Coords(-2.4631,-3.5032,0.6947,"C2", 1.0)});
    sideChains["DG"].insert({"N2",Coords(-3.7877,-3.5683,0.8977,"N2", 1.0)});
    sideChains["DG"].insert({"N3",Coords(-1.8862,-2.3320,0.4853,"N3", 1.0)});
    sideChains["DG"].insert({"C4",Coords(-0.5546,-2.4519,0.3000,"C4", 1.0)});

    // total models used :: 51
    centeringAtoms["DT"] = "C1'";
    // sideChains["DT"].insert({"P",Coords(-0.0774,-4.6557,-2.3600,"P", 1.0)});
    // sideChains["DT"].insert({"OP1",Coords(-1.5541,-4.4989,-2.3340,"OP1", 1.0)});
    // sideChains["DT"].insert({"OP2",Coords(0.5769,-5.2460,-3.5550,"OP2", 1.0)});
    // sideChains["DT"].insert({"O5'",Coords(0.5674,-3.2264,-2.0775,"O5'", 1.0)});
    // sideChains["DT"].insert({"C5'",Coords(0.2712,-2.1139,-2.9221,"C5'", 1.0)});
    // sideChains["DT"].insert({"C4'",Coords(0.4275,-0.8186,-2.1592,"C4'", 1.0)});
    // sideChains["DT"].insert({"O4'",Coords(-0.5195,-0.7797,-1.0709,"O4'", 1.0)});
    // sideChains["DT"].insert({"C3'",Coords(1.7907,-0.6063,-1.4875,"C3'", 1.0)});
    // sideChains["DT"].insert({"O3'",Coords(2.6805,0.0424,-2.4022,"O3'", 1.0)});
    // sideChains["DT"].insert({"C2'",Coords(1.4564,0.3220,-0.3331,"C2'", 1.0)});
    // sideChains["DT"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DT"].insert({"N1",Coords(-0.1812,-0.7673,1.2451,"N1", 1.0)});
    sideChains["DT"].insert({"C2",Coords(-0.3937,-0.0703,2.4131,"C2", 1.0)});
    sideChains["DT"].insert({"O2",Coords(-0.4365,1.1481,2.4660,"O2", 1.0)});
    sideChains["DT"].insert({"N3",Coords(-0.5578,-0.8585,3.5256,"N3", 1.0)});
    sideChains["DT"].insert({"C4",Coords(-0.5292,-2.2386,3.5830,"C4", 1.0)});
    sideChains["DT"].insert({"O4",Coords(-0.6942,-2.8101,4.6575,"O4", 1.0)});
    sideChains["DT"].insert({"C5",Coords(-0.2987,-2.9034,2.3209,"C5", 1.0)});
    sideChains["DT"].insert({"C7",Coords(-0.2498,-4.3977,2.2849,"C7", 1.0)});
    sideChains["DT"].insert({"C6",Coords(-0.1384,-2.1439,1.2303,"C6", 1.0)});

    // total models used :: 148
    centeringAtoms["rA"] = "C1'";
//    sideChains["rA"].insert({"P",Coords(-0.2426,-0.4658,-5.2875,"P", 1.0)});
//    sideChains["rA"].insert({"OP1",Coords(-0.2850,-1.7679,-5.9889,"OP1", 1.0)});
//    sideChains["rA"].insert({"OP2",Coords(-1.2687,0.5539,-5.5850,"OP2", 1.0)});
//    sideChains["rA"].insert({"O5'",Coords(-0.2675,-0.7078,-3.7084,"O5'", 1.0)});
//    sideChains["rA"].insert({"C5'",Coords(0.4972,-1.7366,-3.0964,"C5'", 1.0)});
//    sideChains["rA"].insert({"C4'",Coords(0.2525,-1.6959,-1.6010,"C4'", 1.0)});
//    sideChains["rA"].insert({"O4'",Coords(0.8199,-0.4860,-1.0430,"O4'", 1.0)});
//    sideChains["rA"].insert({"C3'",Coords(-1.2121,-1.6295,-1.2087,"C3'", 1.0)});
//    sideChains["rA"].insert({"O3'",Coords(-1.7388,-2.9429,-1.2071,"O3'", 1.0)});
//    sideChains["rA"].insert({"C2'",Coords(-1.1230,-1.0185,0.1860,"C2'", 1.0)});
//    sideChains["rA"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["rA"].insert({"N9",Coords(-0.5477,1.2937,-0.3901,"N9", 1.0)});
    sideChains["rA"].insert({"C8",Coords(-0.8524,1.6942,-1.6626,"C8", 1.0)});
    sideChains["rA"].insert({"N7",Coords(-1.3653,2.8971,-1.7261,"N7", 1.0)});
    sideChains["rA"].insert({"C5",Coords(-1.4064,3.3028,-0.4045,"C5", 1.0)});
    sideChains["rA"].insert({"C6",Coords(-1.8468,4.4956,0.1873,"C6", 1.0)});
    sideChains["rA"].insert({"N6",Coords(-2.3410,5.5079,-0.5327,"N6", 1.0)});
    sideChains["rA"].insert({"N1",Coords(-1.7619,4.5996,1.5285,"N1", 1.0)});
    sideChains["rA"].insert({"C2",Coords(-1.2588,3.5695,2.2200,"C2", 1.0)});
    sideChains["rA"].insert({"N3",Coords(-0.8096,2.3954,1.7746,"N3", 1.0)});
    sideChains["rA"].insert({"C4",Coords(-0.9111,2.3291,0.4378,"C4", 1.0)});

    // total models used :: 161
    centeringAtoms["rC"] = "C1'";
//    sideChains["CYT"].insert({"P",Coords(-1.1349,-5.2231,0.1749,"P", 1.0)});
//    sideChains["CYT"].insert({"OP1",Coords(-0.2958,-6.2931,-0.4376,"OP1", 1.0)});
//    sideChains["CYT"].insert({"OP2",Coords(-1.5456,-5.3702,1.5762,"OP2", 1.0)});
//    sideChains["CYT"].insert({"O5'",Coords(-0.4386,-3.8143,-0.0160,"O5'", 1.0)});
//    sideChains["CYT"].insert({"C5'",Coords(0.0044,-3.3653,-1.2797,"C5'", 1.0)});
//    sideChains["CYT"].insert({"C4'",Coords(0.5261,-1.9538,-1.1909,"C4'", 1.0)});
//    sideChains["CYT"].insert({"O4'",Coords(-0.5404,-1.0606,-0.7726,"O4'", 1.0)});
//    sideChains["CYT"].insert({"C3'",Coords(1.6081,-1.7002,-0.1512,"C3'", 1.0)});
//    sideChains["CYT"].insert({"O3'",Coords(2.8956,-2.1506,-0.5522,"O3'", 1.0)});
//    sideChains["CYT"].insert({"C2'",Coords(1.5091,-0.1934,0.0374,"C2'", 1.0)});
//    sideChains["CYT"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["rC"].insert({"N1",Coords(-0.5794,-0.0537,1.3552,"N1", 1.0)});
    sideChains["rC"].insert({"C2",Coords(-0.4380,1.1006,2.1377,"C2", 1.0)});
    sideChains["rC"].insert({"O2",Coords(0.1336,2.0803,1.6314,"O2", 1.0)});
    sideChains["rC"].insert({"N3",Coords(-0.9473,1.0821,3.3972,"N3", 1.0)});
    sideChains["rC"].insert({"C4",Coords(-1.5532,-0.0221,3.8647,"C4", 1.0)});
    sideChains["rC"].insert({"N4",Coords(-2.0557,0.0193,5.1057,"N4", 1.0)});
    sideChains["rC"].insert({"C5",Coords(-1.7043,-1.2074,3.0753,"C5", 1.0)});
    sideChains["rC"].insert({"C6",Coords(-1.1942,-1.1866,1.8370,"C6", 1.0)});

    // total models used :: 252
    centeringAtoms["rG"] = "C1'";
//    sideChains["rG"].insert({"P",Coords(-1.6626,4.1606,-2.4487,"P", 1.0)});
//    sideChains["rG"].insert({"OP1",Coords(-3.0618,4.6169,-2.5771,"OP1", 1.0)});
//    sideChains["rG"].insert({"OP2",Coords(-0.6094,5.1543,-2.1564,"OP2", 1.0)});
//    sideChains["rG"].insert({"O5'",Coords(-1.5853,2.9735,-1.4027,"O5'", 1.0)});
//    sideChains["rG"].insert({"C5'",Coords(-2.5535,1.9397,-1.3839,"C5'", 1.0)});
//    sideChains["rG"].insert({"C4'",Coords(-2.1403,0.8663,-0.4218,"C4'", 1.0)});
//    sideChains["rG"].insert({"O4'",Coords(-1.0035,0.1568,-0.9717,"O4'", 1.0)});
//    sideChains["rG"].insert({"C3'",Coords(-1.6815,1.3777,0.9434,"C3'", 1.0)});
//    sideChains["rG"].insert({"O3'",Coords(-2.7537,1.4344,1.8689,"O3'", 1.0)});
//    sideChains["rG"].insert({"C2'",Coords(-0.5902,0.3991,1.3501,"C2'", 1.0)});
//    sideChains["rG"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["rG"].insert({"N9",Coords(1.1562,0.8347,-0.3917,"N9", 1.0)});
    sideChains["rG"].insert({"C8",Coords(1.2781,1.5484,-1.5683,"C8", 1.0)});
    sideChains["rG"].insert({"N7",Coords(2.4122,2.1736,-1.6966,"N7", 1.0)});
    sideChains["rG"].insert({"C5",Coords(3.1189,1.8356,-0.5333,"C5", 1.0)});
    sideChains["rG"].insert({"C6",Coords(4.4313,2.2242,-0.1164,"C6", 1.0)});
    sideChains["rG"].insert({"O6",Coords(5.2502,2.9430,-0.7322,"O6", 1.0)});
    sideChains["rG"].insert({"N1",Coords(4.7470,1.6313,1.1021,"N1", 1.0)});
    sideChains["rG"].insert({"C2",Coords(3.9006,0.8269,1.8321,"C2", 1.0)});
    sideChains["rG"].insert({"N2",Coords(4.3813,0.3582,2.9956,"N2", 1.0)});
    sideChains["rG"].insert({"N3",Coords(2.6746,0.4657,1.4568,"N3", 1.0)});
    sideChains["rG"].insert({"C4",Coords(2.3615,0.9979,0.2638,"C4", 1.0)});


    // total models used :: 116
    centeringAtoms["rU"] = "C1'";
//    sideChains["rU"].insert({"P",Coords(-4.9480,-0.9853,1.6931,"P", 1.0)});
//    sideChains["rU"].insert({"OP1",Coords(-5.3497,-2.3727,2.0213,"OP1", 1.0)});
//    sideChains["rU"].insert({"OP2",Coords(-5.0768,0.0621,2.7379,"OP2", 1.0)});
//    sideChains["rU"].insert({"O5'",Coords(-3.4634,-0.9704,1.1126,"O5'", 1.0)});
//    sideChains["rU"].insert({"C5'",Coords(-3.0815,-1.7942,0.0242,"C5'", 1.0)});
//    sideChains["rU"].insert({"C4'",Coords(-1.5997,-1.6964,-0.2354,"C4'", 1.0)});
//    sideChains["rU"].insert({"O4'",Coords(-1.2498,-0.3271,-0.5626,"O4'", 1.0)});
//    sideChains["rU"].insert({"C3'",Coords(-0.7034,-2.0424,0.9460,"C3'", 1.0)});
//    sideChains["rU"].insert({"O3'",Coords(-0.4743,-3.4285,1.0643,"O3'", 1.0)});
//    sideChains["rU"].insert({"C2'",Coords(0.5624,-1.2480,0.6653,"C2'", 1.0)});
//    sideChains["rU"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["rU"].insert({"N1",Coords(-0.1779,1.1400,0.9375,"N1", 1.0)});
    sideChains["rU"].insert({"C2",Coords(0.9527,1.8163,1.3541,"C2", 1.0)});
    sideChains["rU"].insert({"O2",Coords(2.0645,1.4433,1.0247,"O2", 1.0)});
    sideChains["rU"].insert({"N3",Coords(0.7407,2.9089,2.1765,"N3", 1.0)});
    sideChains["rU"].insert({"C4",Coords(-0.4840,3.3943,2.6058,"C4", 1.0)});
    sideChains["rU"].insert({"O4",Coords(-0.5393,4.3980,3.3216,"O4", 1.0)});
    sideChains["rU"].insert({"C5",Coords(-1.6018,2.6485,2.1151,"C5", 1.0)});
    sideChains["rU"].insert({"C6",Coords(-1.4154,1.5879,1.3199,"C6", 1.0)});

    // total models used :: 8
    centeringAtoms["rI"] = "C1'";
//    sideChains["rI"].insert({"P",Coords(-4.0551,1.3171,3.1999,"P", 1.0)});
//    sideChains["rI"].insert({"OP1",Coords(-3.9486,2.1747,4.3919,"OP1", 1.0)});
//    sideChains["rI"].insert({"OP2",Coords(-5.2299,1.4569,2.3068,"OP2", 1.0)});
//    sideChains["rI"].insert({"O5'",Coords(-2.7229,1.1124,2.3356,"O5'", 1.0)});
//    sideChains["rI"].insert({"C5'",Coords(-1.4998,1.5638,2.8689,"C5'", 1.0)});
//    sideChains["rI"].insert({"C4'",Coords(-0.3549,1.2799,1.9077,"C4'", 1.0)});
//    sideChains["rI"].insert({"O4'",Coords(-0.4740,-0.0083,1.3344,"O4'", 1.0)});
//    sideChains["rI"].insert({"C3'",Coords(-0.3305,2.2464,0.7367,"C3'", 1.0)});
//    sideChains["rI"].insert({"O3'",Coords(0.4047,3.4106,1.0309,"O3'", 1.0)});
//    sideChains["rI"].insert({"C2'",Coords(0.3508,1.4387,-0.3483,"C2'", 1.0)});
//    sideChains["rI"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["rI"].insert({"N9",Coords(-1.0519,-0.4535,-0.9274,"N9", 1.0)});
    sideChains["rI"].insert({"C8",Coords(-2.3850,-0.6097,-0.6518,"C8", 1.0)});
    sideChains["rI"].insert({"N7",Coords(-3.0087,-1.0380,-1.7730,"N7", 1.0)});
    sideChains["rI"].insert({"C5",Coords(-2.0913,-1.1540,-2.7578,"C5", 1.0)});
    sideChains["rI"].insert({"C6",Coords(-2.1748,-1.5459,-4.0847,"C6", 1.0)});
    sideChains["rI"].insert({"O6",Coords(-3.2560,-1.8795,-4.5686,"O6", 1.0)});
    sideChains["rI"].insert({"N1",Coords(-1.0320,-1.5629,-4.8593,"N1", 1.0)});
    sideChains["rI"].insert({"C2",Coords(0.1603,-1.1985,-4.3244,"C2", 1.0)});
    sideChains["rI"].insert({"N3",Coords(0.2461,-0.8119,-3.0273,"N3", 1.0)});
    sideChains["rI"].insert({"C4",Coords(-0.8582,-0.7868,-2.2432,"C4", 1.0)});

    // total models used :: 257
    centeringAtoms["RNA"] = "C4'";  // use this for backbone water placements
    sideChains["RNA"].insert({"P",Coords(-1.1345,-3.2655,-1.6251,"P", 1.0)});
//    sideChains["RNA"].insert({"OP1",Coords(-2.6256,-3.3246,-1.5605,"OP1", 1.0)});
//    sideChains["RNA"].insert({"OP2",Coords(-0.4449,-3.9022,-2.7792,"OP2", 1.0)});
    sideChains["RNA"].insert({"O5'",Coords(-0.6601,-1.7730,-1.4855,"O5'", 1.0)});
    sideChains["RNA"].insert({"C5'",Coords(-1.0386,-1.0170,-0.3430,"C5'", 1.0)});
    sideChains["RNA"].insert({"C4'",Coords(0.0000,0.0000,0.0000,"C4'", 1.0)});
    sideChains["RNA"].insert({"O4'",Coords(1.2588,-0.6633,0.2718,"O4'", 1.0)});
    sideChains["RNA"].insert({"C3'",Coords(0.3470,0.9650,-1.0993,"C3'", 1.0)});
    sideChains["RNA"].insert({"O3'",Coords(-0.6173,1.9918,-1.2240,"O3'", 1.0)});
    sideChains["RNA"].insert({"C2'",Coords(1.7234,1.4527,-0.6744,"C2'", 1.0)});
    sideChains["RNA"].insert({"C1'",Coords(2.3307,0.1507,-0.1754,"C1'", 1.0)});

    // any new residue requires set of centered coordinates - maybe averaged followed by a set of hydrating waters
    // total models used :: 126
    centeringAtoms["DNA"] = "C4'";
    sideChains["DNA"].insert({"P",Coords(3.2223,-0.4999,-2.1642,"P", 1.0)});
    // sideChains["NBCKBONE"].insert({"OP1",Coords(4.5557,-0.0301,-1.7101,"OP1", 1.0)});
    // sideChains["NBCKBONE"].insert({"OP2",Coords(3.0679,-1.8860,-2.6754,"OP2", 1.0)});
    sideChains["DNA"].insert({"O5'",Coords(2.1903,-0.2926,-0.9690,"O5'", 1.0)});
    sideChains["DNA"].insert({"C5'",Coords(0.8012,-0.5690,-1.1474,"C5'", 1.0)});
//    sideChains["NBCKBONE"].insert({"C4'",Coords(0.0000,0.0000,0.0000,"C4'", 1.0)});
    sideChains["DNA"].insert({"O4'",Coords(0.3568,1.3808,0.2139,"O4'", 1.0)});
    sideChains["DNA"].insert({"C3'",Coords(-1.4963,0.0244,-0.2605,"C3'", 1.0)});
    sideChains["DNA"].insert({"O3'",Coords(-2.1030,-1.2527,-0.0437,"O3'", 1.0)});
    sideChains["DNA"].insert({"C2'",Coords(-1.9832,1.1621,0.6170,"C2'", 1.0)});
    sideChains["DNA"].insert({"C1'",Coords(-0.7781,2.1041,0.6834,"C1'", 1.0)});


    std::vector<Coords> * pwater;
    // Create Vector for new Residue ALA
    waters.insert(std::pair<std::string, std::vector<Coords> >("ALA", std::vector<Coords>(6)));
    pwater = &waters["ALA"];
    (*pwater)[0].x = -1.299;
    (*pwater)[0].y = 3.036;
    (*pwater)[0].z = -0.996;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 2.406;
    (*pwater)[1].y = -1.987;
    (*pwater)[1].z = 2.319;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 3.658;
    (*pwater)[2].y = 2.160;
    (*pwater)[2].z = 2.457;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 3.997;
    (*pwater)[3].y = 1.773;
    (*pwater)[3].z = -0.600;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 1.511;
    (*pwater)[4].y = 3.980;
    (*pwater)[4].z = -0.907;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -2.309;
    (*pwater)[5].y = 2.363;
    (*pwater)[5].z = 1.490;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("ARG", std::vector<Coords>(5)));
    pwater = &waters["ARG"];
    (*pwater)[0].x = -4.215;
    (*pwater)[0].y = 1.614;
    (*pwater)[0].z = -1.893;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 1.672;
    (*pwater)[1].y = 1.344;
    (*pwater)[1].z = 2.044;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -1.769;
    (*pwater)[2].y = -2.329;
    (*pwater)[2].z = -2.744;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -3.882;
    (*pwater)[3].y = 2.991;
    (*pwater)[3].z = -0.349;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 0.130;
    (*pwater)[4].y = 3.322;
    (*pwater)[4].z = 2.075;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("ASN", std::vector<Coords>(5)));
    pwater = &waters["ASN"];
    (*pwater)[0].x = 0.351;
    (*pwater)[0].y = 1.952;
    (*pwater)[0].z = -3.291;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -0.334;
    (*pwater)[1].y = -2.669;
    (*pwater)[1].z = -2.559;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -1.654;
    (*pwater)[2].y = -2.994;
    (*pwater)[2].z = 0.572;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 0.738;
    (*pwater)[3].y = -3.629;
    (*pwater)[3].z = 0.620;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 0.002;
    (*pwater)[4].y = -0.874;
    (*pwater)[4].z = 3.565;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;


    waters.insert(std::pair<std::string, std::vector<Coords> >("CYS", std::vector<Coords>(9)));
    pwater = &waters["CYS"];
    (*pwater)[0].x = 2.904;
    (*pwater)[0].y = 0.654;
    (*pwater)[0].z = -1.647;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 2.324;
    (*pwater)[1].y = -1.853;
    (*pwater)[1].z = 3.007;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.974;
    (*pwater)[2].y = 3.168;
    (*pwater)[2].z = 2.989;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 0.674;
    (*pwater)[3].y = 1.355;
    (*pwater)[3].z = 4.655;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -1.097;
    (*pwater)[4].y = -0.322;
    (*pwater)[4].z = 3.649;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -1.829;
    (*pwater)[5].y = 2.984;
    (*pwater)[5].z = 0.655;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = 4.516;
    (*pwater)[6].y = 0.666;
    (*pwater)[6].z = 0.832;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = -0.219;
    (*pwater)[7].y = 3.317;
    (*pwater)[7].z = -0.993;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = 0.968;
    (*pwater)[8].y = 4.339;
    (*pwater)[8].z = 1.840;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("GLN", std::vector<Coords>(6)));
    pwater = &waters["GLN"];
    (*pwater)[0].x = -3.520;
    (*pwater)[0].y = 0.248;
    (*pwater)[0].z = -0.669;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 2.862;
    (*pwater)[1].y = 1.327;
    (*pwater)[1].z = -1.708;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 3.123;
    (*pwater)[2].y = -1.974;
    (*pwater)[2].z = 0.500;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 0.358;
    (*pwater)[3].y = 1.310;
    (*pwater)[3].z = -3.374;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 3.294;
    (*pwater)[4].y = 0.681;
    (*pwater)[4].z = 1.720;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -1.541;
    (*pwater)[5].y = 2.976;
    (*pwater)[5].z = 1.658;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("GLU", std::vector<Coords>(7)));
    pwater = &waters["GLU"];
    (*pwater)[0].x = -2.974;
    (*pwater)[0].y = -1.017;
    (*pwater)[0].z = 1.225;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 2.302;
    (*pwater)[1].y = -3.554;
    (*pwater)[1].z = 2.483;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 4.108;
    (*pwater)[2].y = -0.680;
    (*pwater)[2].z = 0.695;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -0.174;
    (*pwater)[3].y = -3.655;
    (*pwater)[3].z = 3.066;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 2.864;
    (*pwater)[4].y = 0.634;
    (*pwater)[4].z = 2.382;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = 2.919;
    (*pwater)[5].y = -2.008;
    (*pwater)[5].z = -1.562;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -2.312;
    (*pwater)[6].y = -2.774;
    (*pwater)[6].z = -0.723;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("GLY", std::vector<Coords>(7)));
    pwater = &waters["GLY"];
    (*pwater)[0].x = -0.934;
    (*pwater)[0].y = 3.308;
    (*pwater)[0].z = -0.981;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -2.630;
    (*pwater)[1].y = 0.158;
    (*pwater)[1].z = 2.597;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -0.357;
    (*pwater)[2].y = 2.664;
    (*pwater)[2].z = 2.262;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -5.227;
    (*pwater)[3].y = -0.291;
    (*pwater)[3].z = -0.657;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 1.769;
    (*pwater)[4].y = -1.353;
    (*pwater)[4].z = 2.649;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = 2.582;
    (*pwater)[5].y = 1.664;
    (*pwater)[5].z = 1.683;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -0.911;
    (*pwater)[6].y = -3.318;
    (*pwater)[6].z = 0.632;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("HIS", std::vector<Coords>(3)));
    pwater = &waters["HIS"];
    (*pwater)[0].x = 0.991;
    (*pwater)[0].y = -1.571;
    (*pwater)[0].z = -3.399;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -0.978;
    (*pwater)[1].y = -3.200;
    (*pwater)[1].z = 3.463;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -0.016;
    (*pwater)[2].y = -5.374;
    (*pwater)[2].z = 0.141;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("ILE", std::vector<Coords>(9)));
    pwater = &waters["ILE"];
    (*pwater)[0].x = 1.648;
    (*pwater)[0].y = 1.236;
    (*pwater)[0].z = -3.892;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -1.781;
    (*pwater)[1].y = -0.986;
    (*pwater)[1].z = 2.987;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.156;
    (*pwater)[2].y = -2.978;
    (*pwater)[2].z = -2.027;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -0.140;
    (*pwater)[3].y = 0.652;
    (*pwater)[3].z = 4.001;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 3.312;
    (*pwater)[4].y = -0.063;
    (*pwater)[4].z = 3.325;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -0.109;
    (*pwater)[5].y = 4.219;
    (*pwater)[5].z = -2.305;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -0.653;
    (*pwater)[6].y = -3.421;
    (*pwater)[6].z = 2.679;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = 3.299;
    (*pwater)[7].y = -0.860;
    (*pwater)[7].z = -2.231;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = 2.707;
    (*pwater)[8].y = -3.985;
    (*pwater)[8].z = 0.462;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("LEU", std::vector<Coords>(10)));
    pwater = &waters["LEU"];
    (*pwater)[0].x = -3.646;
    (*pwater)[0].y = -1.638;
    (*pwater)[0].z = 1.507;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 2.970;
    (*pwater)[1].y = -1.555;
    (*pwater)[1].z = -1.326;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -3.020;
    (*pwater)[2].y = 1.716;
    (*pwater)[2].z = 3.652;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -1.231;
    (*pwater)[3].y = -1.672;
    (*pwater)[3].z = 4.243;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -1.651;
    (*pwater)[4].y = 0.201;
    (*pwater)[4].z = -4.112;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -0.636;
    (*pwater)[5].y = 4.308;
    (*pwater)[5].z = -2.353;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -0.882;
    (*pwater)[6].y = 2.311;
    (*pwater)[6].z = -4.316;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = -0.918;
    (*pwater)[7].y = 0.950;
    (*pwater)[7].z = 4.608;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = 1.687;
    (*pwater)[8].y = 0.389;
    (*pwater)[8].z = 3.824;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;
    (*pwater)[9].x = 2.799;
    (*pwater)[9].y = 1.453;
    (*pwater)[9].z = 1.855;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("LYS", std::vector<Coords>(7)));
    pwater = &waters["LYS"];
    (*pwater)[0].x = 1.977;
    (*pwater)[0].y = 1.361;
    (*pwater)[0].z = 2.642;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -0.221;
    (*pwater)[1].y = 4.798;
    (*pwater)[1].z = 0.023;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.731;
    (*pwater)[2].y = 1.562;
    (*pwater)[2].z = -1.724;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -3.636;
    (*pwater)[3].y = 1.688;
    (*pwater)[3].z = 0.197;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -2.302;
    (*pwater)[4].y = 1.662;
    (*pwater)[4].z = 2.602;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -2.036;
    (*pwater)[5].y = 2.009;
    (*pwater)[5].z = -2.966;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = 0.688;
    (*pwater)[6].y = 0.292;
    (*pwater)[6].z = -3.656;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("MET", std::vector<Coords>(9)));
    pwater = &waters["MET"];
    (*pwater)[0].x = -2.428;
    (*pwater)[0].y = 1.752;
    (*pwater)[0].z = 2.321;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -1.221;
    (*pwater)[1].y = 4.306;
    (*pwater)[1].z = -0.630;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 0.826;
    (*pwater)[2].y = 4.194;
    (*pwater)[2].z = -2.322;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -1.214;
    (*pwater)[3].y = 4.326;
    (*pwater)[3].z = 1.230;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 2.971;
    (*pwater)[4].y = 2.749;
    (*pwater)[4].z = -2.308;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -2.650;
    (*pwater)[5].y = 2.670;
    (*pwater)[5].z = -3.041;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -3.396;
    (*pwater)[6].y = -0.465;
    (*pwater)[6].z = -2.672;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = -1.691;
    (*pwater)[7].y = 0.203;
    (*pwater)[7].z = -4.309;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = -0.415;
    (*pwater)[8].y = -1.120;
    (*pwater)[8].z = 3.294;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("PHE", std::vector<Coords>(4)));
    pwater = &waters["PHE"];
    (*pwater)[0].x = 4.604;
    (*pwater)[0].y = 2.031;
    (*pwater)[0].z = -3.309;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 0.855;
    (*pwater)[1].y = 1.799;
    (*pwater)[1].z = -5.651;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -3.008;
    (*pwater)[2].y = 0.030;
    (*pwater)[2].z = -4.091;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 4.596;
    (*pwater)[3].y = 0.799;
    (*pwater)[3].z = 1.311;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("PRO", std::vector<Coords>(9)));
    pwater = &waters["PRO"];
    (*pwater)[0].x = 3.300;
    (*pwater)[0].y = 0.129;
    (*pwater)[0].z = 0.569;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 0.764;
    (*pwater)[1].y = -3.314;
    (*pwater)[1].z = 1.001;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.576;
    (*pwater)[2].y = 0.655;
    (*pwater)[2].z = -2.310;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -0.651;
    (*pwater)[3].y = 3.420;
    (*pwater)[3].z = 0.470;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -1.155;
    (*pwater)[4].y = 1.346;
    (*pwater)[4].z = 4.627;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -3.651;
    (*pwater)[5].y = 0.550;
    (*pwater)[5].z = 2.048;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = -0.631;
    (*pwater)[6].y = -1.886;
    (*pwater)[6].z = 4.324;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = 2.006;
    (*pwater)[7].y = 2.785;
    (*pwater)[7].z = 0.922;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = -2.211;
    (*pwater)[8].y = -2.760;
    (*pwater)[8].z = 1.913;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("SER", std::vector<Coords>(3)));
    pwater = &waters["SER"];
    (*pwater)[0].x = -2.334;
    (*pwater)[0].y = -2.628;
    (*pwater)[0].z = -0.485;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 0.044;
    (*pwater)[1].y = -1.332;
    (*pwater)[1].z = -3.298;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.935;
    (*pwater)[2].y = -1.666;
    (*pwater)[2].z = -0.778;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("THR", std::vector<Coords>(3)));
    pwater = &waters["THR"];
    (*pwater)[0].x = 2.206;
    (*pwater)[0].y = 2.338;
    (*pwater)[0].z = -1.321;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -0.764;
    (*pwater)[1].y = 0.721;
    (*pwater)[1].z = -3.329;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 3.319;
    (*pwater)[2].y = -0.681;
    (*pwater)[2].z = -1.231;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("TRP", std::vector<Coords>(3)));
    pwater = &waters["TRP"];
    (*pwater)[0].x = 0.713;
    (*pwater)[0].y = 1.078;
    (*pwater)[0].z = 4.951;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 3.772;
    (*pwater)[1].y = 0.341;
    (*pwater)[1].z = 4.037;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 5.376;
    (*pwater)[2].y = -0.663;
    (*pwater)[2].z = 1.093;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("TYR", std::vector<Coords>(2)));
    pwater = &waters["TYR"];
    (*pwater)[0].x = -5.626;
    (*pwater)[0].y = 0.320;
    (*pwater)[0].z = -1.529;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -3.670;
    (*pwater)[1].y = 4.438;
    (*pwater)[1].z = -1.201;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;

    waters.insert(std::pair<std::string, std::vector<Coords> >("VAL", std::vector<Coords>(12)));
    pwater = &waters["VAL"];
    (*pwater)[0].x = 0.155;
    (*pwater)[0].y = -1.548;
    (*pwater)[0].z = -3.641;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 0.570;
    (*pwater)[1].y = 3.539;
    (*pwater)[1].z = 1.288;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -3.980;
    (*pwater)[2].y = 1.502;
    (*pwater)[2].z = 2.391;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -1.750;
    (*pwater)[3].y = 2.916;
    (*pwater)[3].z = -1.320;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -3.982;
    (*pwater)[4].y = -0.241;
    (*pwater)[4].z = -0.669;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = -0.039;
    (*pwater)[5].y = -4.864;
    (*pwater)[5].z = -0.782;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;
    (*pwater)[6].x = 0.664;
    (*pwater)[6].y = 2.017;
    (*pwater)[6].z = 3.741;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.000;
    (*pwater)[7].x = -1.161;
    (*pwater)[7].y = -3.591;
    (*pwater)[7].z = -2.365;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.000;
    (*pwater)[8].x = -3.481;
    (*pwater)[8].y = -2.195;
    (*pwater)[8].z = 0.087;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.000;
    (*pwater)[9].x = 3.295;
    (*pwater)[9].y = -3.313;
    (*pwater)[9].z = 0.894;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.000;
    (*pwater)[10].x = 3.661;
    (*pwater)[10].y = -1.352;
    (*pwater)[10].z = 1.707;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 1.000;
    (*pwater)[11].x = -1.567;
    (*pwater)[11].y = -1.723;
    (*pwater)[11].z = 3.969;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 1.000;

    // backbone
    waters.insert(std::pair<std::string, std::vector<Coords> >("PROTEIN", std::vector<Coords>(5)));
    pwater = &waters["PROTEIN"];
    (*pwater)[0].x = -2.277;
    (*pwater)[0].y = -2.609;
    (*pwater)[0].z = -1.341;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -1.422;
    (*pwater)[1].y = -2.658;
    (*pwater)[1].z = 1.802;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -1.333;
    (*pwater)[2].y = -0.593;
    (*pwater)[2].z = -4.259;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 1.445;
    (*pwater)[3].y = 3.189;
    (*pwater)[3].z = -1.527;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = 0.762;
    (*pwater)[4].y = 1.604;
    (*pwater)[4].z = -3.662;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;


    // ADE in RNA
    // Create Vector for new Residue A
    waters["rA"] = std::vector<Coords>(14);
    pwater = &waters["rA"];
    (*pwater)[0].x = -3.5162;
    (*pwater)[0].y = 4.3393;
    (*pwater)[0].z = -3.5739;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -3.5499;
    (*pwater)[1].y = 3.4190;
    (*pwater)[1].z = 6.3245;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -1.8203;
    (*pwater)[2].y = 12.1766;
    (*pwater)[2].z = 1.8823;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = -5.0946;
    (*pwater)[3].y = -1.9847;
    (*pwater)[3].z = -2.3524;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = -3.1875;
    (*pwater)[4].y = 8.0440;
    (*pwater)[4].z = 0.9046;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -0.5154;
    (*pwater)[5].y = -4.6238;
    (*pwater)[5].z = 1.9894;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 0.7010;
    (*pwater)[6].y = 0.6447;
    (*pwater)[6].z = 3.8475;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = -6.8347;
    (*pwater)[7].y = 3.5568;
    (*pwater)[7].z = -5.5333;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = -0.2917;
    (*pwater)[8].y = 11.3987;
    (*pwater)[8].z = 4.6815;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = -1.4407;
    (*pwater)[9].y = 9.4941;
    (*pwater)[9].z = -0.7456;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;
    (*pwater)[10].x = -2.4892;
    (*pwater)[10].y = 1.3303;
    (*pwater)[10].z = -7.4622;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 1.0000;
    (*pwater)[11].x = 3.9484;
    (*pwater)[11].y = -0.4974;
    (*pwater)[11].z = -1.1912;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 1.0000;
    (*pwater)[12].x = 3.4358;
    (*pwater)[12].y = 0.2977;
    (*pwater)[12].z = 1.2967;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 1.0000;
    (*pwater)[13].x = 2.9584;
    (*pwater)[13].y = -2.8113;
    (*pwater)[13].z = -5.4033;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 1.0000;


// Create Vector for new Residue G
    waters["rG"] = std::vector<Coords>(14);
    pwater = &waters["rG"];
    (*pwater)[0].x = 0.4793;
    (*pwater)[0].y = 5.4623;
    (*pwater)[0].z = -1.5606;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -5.0206;
    (*pwater)[1].y = 0.2953;
    (*pwater)[1].z = 2.2531;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 0.7592;
    (*pwater)[2].y = 4.2446;
    (*pwater)[2].z = -5.7303;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 1.5897;
    (*pwater)[3].y = -1.3761;
    (*pwater)[3].z = 4.1003;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 4.7960;
    (*pwater)[4].y = 1.4949;
    (*pwater)[4].z = 7.0876;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -4.3525;
    (*pwater)[5].y = -1.2387;
    (*pwater)[5].z = -5.1036;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = -4.9237;
    (*pwater)[6].y = 4.4384;
    (*pwater)[6].z = -1.9495;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = 0.4046;
    (*pwater)[7].y = 7.6985;
    (*pwater)[7].z = 1.3661;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = 2.9967;
    (*pwater)[8].y = 4.1110;
    (*pwater)[8].z = -3.6331;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = 5.5047;
    (*pwater)[9].y = 5.0157;
    (*pwater)[9].z = -2.3983;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;
    (*pwater)[10].x = 12.6400;
    (*pwater)[10].y = 4.1757;
    (*pwater)[10].z = 0.4429;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 1.0000;
    (*pwater)[11].x = 7.0910;
    (*pwater)[11].y = 4.0643;
    (*pwater)[11].z = 7.6615;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 1.0000;
    (*pwater)[12].x = 4.8785;
    (*pwater)[12].y = 7.9540;
    (*pwater)[12].z = -1.8751;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 1.0000;
    (*pwater)[13].x = 6.8735;
    (*pwater)[13].y = 9.2044;
    (*pwater)[13].z = -2.5838;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 1.0000;

    // Create Vector for new Residue U
    waters["rU"] = std::vector<Coords>(8);
    pwater = &waters["rU"];
    (*pwater)[0].x = -3.5019;
    (*pwater)[0].y = 0.7294;
    (*pwater)[0].z = 4.2101;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 0.5694;
    (*pwater)[1].y = -4.8909;
    (*pwater)[1].z = -0.5521;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -4.6876;
    (*pwater)[2].y = 5.1135;
    (*pwater)[2].z = 8.8948;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 4.0947;
    (*pwater)[3].y = -0.9770;
    (*pwater)[3].z = 0.1620;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 7.2768;
    (*pwater)[4].y = 0.5164;
    (*pwater)[4].z = 4.3032;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -8.4280;
    (*pwater)[5].y = 1.6611;
    (*pwater)[5].z = 2.3350;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = -3.2258;
    (*pwater)[6].y = -5.1515;
    (*pwater)[6].z = -0.9635;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = -1.1344;
    (*pwater)[7].y = 1.3280;
    (*pwater)[7].z = -3.8784;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;

    // Create Vector for new Residue C
    waters["rC"] = std::vector<Coords>(9);
    pwater = &waters["rC"];
    (*pwater)[0].x = -3.3560;
    (*pwater)[0].y = -2.4191;
    (*pwater)[0].z = 6.0745;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 3.5198;
    (*pwater)[1].y = -0.3978;
    (*pwater)[1].z = -3.3341;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -1.9399;
    (*pwater)[2].y = -5.5953;
    (*pwater)[2].z = -3.6108;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 4.0252;
    (*pwater)[3].y = -3.4458;
    (*pwater)[3].z = 5.4060;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = -0.0218;
    (*pwater)[4].y = -3.9137;
    (*pwater)[4].z = 4.1516;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 1.6532;
    (*pwater)[5].y = 3.6588;
    (*pwater)[5].z = -0.5053;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 3.6128;
    (*pwater)[6].y = 6.1194;
    (*pwater)[6].z = 3.9159;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = 2.6020;
    (*pwater)[7].y = -6.3601;
    (*pwater)[7].z = 0.2461;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = -3.3227;
    (*pwater)[8].y = 4.1352;
    (*pwater)[8].z = 9.7322;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;

    waters["rI"] = std::vector<Coords>(2);
    pwater = &waters["rI"];
    (*pwater)[0].x = 3.0091;
    (*pwater)[0].y = -0.2373;
    (*pwater)[0].z = -2.2087;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 3.0372;
    (*pwater)[1].y = 0.0299;
    (*pwater)[1].z = -5.3884;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;


    // RNA backbone
    // Create Vector for new Residue RNA backbone
    waters["RNA"] = std::vector<Coords>(10);
    pwater = &waters["RNA"];
    (*pwater)[0].x = 0.9132;
    (*pwater)[0].y = -2.2776;
    (*pwater)[0].z = -4.7812;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -3.7335;
    (*pwater)[1].y = -0.3780;
    (*pwater)[1].z = -3.4614;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 2.6709;
    (*pwater)[2].y = -0.9810;
    (*pwater)[2].z = 3.4014;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = -0.8079;
    (*pwater)[3].y = 3.6259;
    (*pwater)[3].z = 1.3553;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = -2.2406;
    (*pwater)[4].y = 5.2196;
    (*pwater)[4].z = -0.5451;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -2.2852;
    (*pwater)[5].y = -6.0268;
    (*pwater)[5].z = -4.0371;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = -3.3756;
    (*pwater)[6].y = -3.4238;
    (*pwater)[6].z = 1.6775;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = 1.3288;
    (*pwater)[7].y = 1.6198;
    (*pwater)[7].z = 3.3486;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = 4.3270;
    (*pwater)[8].y = 3.4083;
    (*pwater)[8].z = 1.1492;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = -4.3766;
    (*pwater)[9].y = 0.9460;
    (*pwater)[9].z = -0.1951;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;

    waters["DNA"] = std::vector<Coords>(5);
    pwater = &waters["DNA"];

    // (*pwater)[0].x = -0.928;
    // (*pwater)[0].y = 0.476;
    // (*pwater)[0].z = 3.976;
    // (*pwater)[0].type = "O";
    // (*pwater)[0].occ = 1.000;
    (*pwater)[0].x = -4.981;
    (*pwater)[0].y = 1.515;
    (*pwater)[0].z = -2.208;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -0.785;
    (*pwater)[1].y = 2.654;
    (*pwater)[1].z = -4.064;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 0.335;
    (*pwater)[2].y = -2.958;
    (*pwater)[2].z = -3.220;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 1.991;
    (*pwater)[3].y = 0.659;
    (*pwater)[3].z = 3.077;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -2.007;
    (*pwater)[4].y = -3.291;
    (*pwater)[4].z = 2.424;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;

    // (*pwater)[0].x = -4.9789;
    // (*pwater)[0].y = 1.5232;
    // (*pwater)[0].z = -2.3501;
    // (*pwater)[0].type = "O";
    // (*pwater)[0].occ = 1.0000;
    // (*pwater)[1].x = -1.0576;
    // (*pwater)[1].y = 1.0054;
    // (*pwater)[1].z = 3.8727;
    // (*pwater)[1].type = "O";
    // (*pwater)[1].occ = 1.0000;
    // (*pwater)[1].x = 0.0059;
    // (*pwater)[1].y = -3.0690;
    // (*pwater)[1].z = -3.3615;
    // (*pwater)[1].type = "O";
    // (*pwater)[1].occ = 1.0000;
    // (*pwater)[2].x = -0.8245;
    // (*pwater)[2].y = 2.7030;
    // (*pwater)[2].z = -4.1705;
    // (*pwater)[2].type = "O";
    // (*pwater)[2].occ = 1.0000;
    // (*pwater)[2].x = 1.9927;    // seems to be a phosphate hydration
    // (*pwater)[2].y = 0.7604;    // seems to be a phosphate hydration
    // (*pwater)[2].z = 3.0047;    // seems to be a phosphate hydration
    // (*pwater)[2].type = "O";    // seems to be a phosphate hydration
    // (*pwater)[2].occ = 1.0000;  // seems to be a phosphate hydration


    // Create Vector for new Residue DC
    waters["DC"] = std::vector<Coords>(4);
    pwater = &waters["DC"];
    (*pwater)[0].x = 1.680;
    (*pwater)[0].y = -3.298;
    (*pwater)[0].z = 0.528;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -7.066;
    (*pwater)[1].y = 0.612;
    (*pwater)[1].z = -1.301;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -3.754;
    (*pwater)[2].y = 3.737;
    (*pwater)[2].z = -0.020;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -5.793;
    (*pwater)[3].y = 3.105;
    (*pwater)[3].z = -1.303;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;

    // Create Vector for new Residue DG
    waters["DG"] = std::vector<Coords>(6);
    pwater = &waters["DG"];
    pwater = &waters["DG"];
    (*pwater)[0].x = -3.636;
    (*pwater)[0].y = -0.081;
    (*pwater)[0].z = -0.280;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -5.521;
    (*pwater)[1].y = -1.188;
    (*pwater)[1].z = 1.621;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 2.561;
    (*pwater)[2].y = -7.053;
    (*pwater)[2].z = 0.674;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = 3.991;
    (*pwater)[3].y = -4.637;
    (*pwater)[3].z = 0.302;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;
    (*pwater)[4].x = -7.020;
    (*pwater)[4].y = -4.227;
    (*pwater)[4].z = 1.984;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.000;
    (*pwater)[5].x = 4.403;
    (*pwater)[5].y = -0.088;
    (*pwater)[5].z = -0.295;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.000;

    // Create Vector for new Residue DA
    waters["DA"] = std::vector<Coords>(4);
    pwater = &waters["DA"];
    (*pwater)[0].x = -3.322;
    (*pwater)[0].y = 5.140;
    (*pwater)[0].z = 0.333;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = 1.631;
    (*pwater)[1].y = -1.073;
    (*pwater)[1].z = -3.270;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = -3.633;
    (*pwater)[2].y = 6.341;
    (*pwater)[2].z = -1.456;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
    (*pwater)[3].x = -1.965;
    (*pwater)[3].y = 2.416;
    (*pwater)[3].z = 3.304;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.000;


    // Create Vector for new Residue DT
    waters["DT"] = std::vector<Coords>(3);         // 3
    pwater = &waters["DT"];
    (*pwater)[0].x = -1.010;
    (*pwater)[0].y = -5.369;
    (*pwater)[0].z = 5.629;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.000;
    (*pwater)[1].x = -1.031;
    (*pwater)[1].y = 3.512;
    (*pwater)[1].z = 2.038;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.000;
    (*pwater)[2].x = 3.355;
    (*pwater)[2].y = -3.950;
    (*pwater)[2].z = 3.038;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.000;
}



/*
 * residues must be arranged as continuous blocks
 * waters are added in Cartesian coordinates
 * Water coordinates need to be converted to spherical using createSphericalCoordinateOfHydration();
 * Use this function to create a water model based on input Atomistic Model
 */
void Waters::hydrateAtomisticModel(AtomisticModel & aModel) {

    SASTOOLS_UTILS_H::logger("", "HYDRATING MODEL");

    tempWaters = new Coords[totalTempWaters];
    hydration.clear();

    auto & model = aModel.getPDBModel();

    unsigned int numAtoms = model.getTotalCoordinates();

    auto pResidues = model.getResiduesVector();

    auto pResIDs = model.getResid();

   // auto atomicNumbers = model.getAtomicNumberVec();

    auto xvalue = model.getCenteredXVec();
    auto yvalue = model.getCenteredYVec();
    auto zvalue = model.getCenteredZVec();

    auto atomType = model.getPointerToAtomTypes();

    unsigned int startHere=0;
    int atomsInResidue = 1;
    int currentID = pResIDs[0];
    std::string currentRes = pResidues[0];

    for (unsigned int i=1; i<numAtoms; i++){ // find block of atoms that belong to same RESID

        if (currentID == pResIDs[i] && currentRes == pResidues[i] && i != (numAtoms-1)){ // keep adding until resid and resname changes

            atomsInResidue++;

        } else {
            // hydrate
            hydrateResidueDirect(currentRes, currentID, atomsInResidue, startHere, atomType, xvalue, yvalue, zvalue);

            //reset
            currentID = pResIDs[i];
            currentRes = pResidues[i];
            atomsInResidue = 1;
            startHere = i;
        }
    }

    std::string residue_index;
    totalwaters = (unsigned int)hydration.size();

    std::clock_t startTime = std::clock();

    // check for waters too close to each other - average their positions
    // average nearest-neighbor Oxygen-Oxygen  distance in liquid water is approximately 2.80 to 2.85 Angstroms
    // most probable is 2.7 to 2.75 Angstrom

    float dwater = 2.65*2.65;

    for(int w=0; w<totalwaters; w++){ // sort the waters.  waters too close will be moved to end of the vector
        auto pW = &hydration[w];

        vector3 vecW = vector3(pW->x, pW->y, pW->z);

        float x_v = pW->x;
        float y_v = pW->y;
        float z_v = pW->z;
        float occ = pW->occ;
        unsigned int counter = 1;

        unsigned int k = w + 1;

        while(k < totalwaters){
            auto pK = &hydration[k];
            if ((vector3(pK->x, pK->y, pK->z) - vecW).sqlength() < dwater){

                x_v += pK->x;
                y_v += pK->y;
                z_v += pK->z;
                occ += pK->occ;
                counter += 1;

                std::iter_swap(hydration.begin()+k, hydration.begin()+(totalwaters-1));
                hydration.pop_back();

                totalwaters--;
            } else {
                k++;
            }
        }

        float inv = 1.0/(float)counter;
        pW->x = inv*x_v;
        pW->y = inv*y_v;
        pW->z = inv*z_v;
        pW->occ = inv*occ;

        //recheck if updated by averaging location
        if (counter > 1) {
            w -= 1;
        }
    }

    /*
     * prune waters too close to protein or nucleic acid
     * does not consider hydrogens
     */
    for (unsigned int i=0; i<numAtoms; i++){

        currentID = pResIDs[i];
        currentRes = pResidues[i];
        boost::trim(currentRes);

        //auto ptr = minima.find(currentRes);
        //float rlimit = 2.2;//(ptr != minima.end()) ? ptr->second : 2.7; // 2.2
        // could fail if residue not found - PDBModel should create the residue during instantiation

        float rlimit = 2.2562;//1.35 + SASTOOLS_RESIDUES_H::Residues::getResidue(currentRes)->getAtom(atomType[i])->getRadii();
        rlimit *= rlimit; // compare the square distance to reduce operations in the loop

        vector3 anchor(xvalue[i], yvalue[i], zvalue[i]);

        for(int w=0; w < totalwaters; w++){ // sort the waters.  waters too close will be moved to end of the vector
            auto pW = &hydration[w]; // these waters belong to a resid
            //float val = (anchor - vector3(pW->x, pW->y, pW->z)).length();

            if ((anchor - vector3(pW->x, pW->y, pW->z)).sqlength() < rlimit ){ // too close for hydrogen bond
                // reject
                std::iter_swap(hydration.begin()+w, hydration.begin()+(totalwaters-1));
                hydration.pop_back();
                w = (w==0) ? 0 : w-1 ;
                totalwaters--;
            }

            if (w >= totalwaters){
                break;
            }
        }
    }

    // for each water, what is the closest protein residue
    // reject waters that are too far from closest atom
    int w=0;
    while(w<totalwaters){
        auto pW = &hydration[w];
        float closest, min_dis = FLT_MAX;

        for (unsigned int i=0; i<numAtoms; i++){
            vector3 anchor(xvalue[i], yvalue[i], zvalue[i]);

            closest = (anchor - vector3(pW->x, pW->y, pW->z)).sqlength();
            if (closest < min_dis){
                min_dis = closest;
            }
        }

        if (min_dis > 16){ // reject water if too far
            std::cout << " Too far " << std::endl;
            std::iter_swap(hydration.begin()+w, hydration.begin()+(totalwaters-1));
            totalwaters--;
        } else {
            w++;
        }
    }

//    std::shuffle(hydration.begin(), hydration.begin()+totalwaters, gen);
//    totalwaters = totalwaters*0.791;

    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    std::cout << " Averaging time " << runtime << std::endl;
    std::cout << " Total waters " << hydration.size() << std::endl;

    totalwatersInExcludedVolume = hydration.size() - totalwaters;

    delete[] tempWaters;

    SASTOOLS_UTILS_H::logger("TOTAL WATERS IN HYDRATION", std::to_string(totalwaters));
    //SASTOOLS_UTILS_H::logger("TOTAL WATERS IN HYDRATION", std::to_string(totalwaters));

    createSphericalCoordinateOfHydration();
}


/*
 * residues must be arranged as continuous blocks
 * waters are added in Cartesian coordinates
 * Water coordinates need to be converted to spherical using createSphericalCoordinateOfHydration();
 * Use this function to create a water model based on input Atomistic Model
 */
void Waters::extractWatersFromAtomisticModel(AtomisticModel &model) {

    SASTOOLS_UTILS_H::logger("", "HYDRATING MODEL");

    boost::regex ifOAtom("^[ ]?O['A-Z0-9]?+");

    auto waterLines = model.getPDBModel().getWaterLines();
    auto centering_vector = model.getPDBModel().getCenteringVector();

    hydration.clear();

    boost::regex wat("HOH");
    float x, y , z;
    std::string atomType;
    int count = 1;
    for(auto & line : waterLines){
        // if water is oxygen,
        if (boost::regex_match(line.substr(17,3), wat)){

            x = std::strtof(line.substr(30,8).c_str(), nullptr);
            y = std::strtof(line.substr(38,8).c_str(), nullptr);
            z = std::strtof(line.substr(46,8).c_str(), nullptr);

            //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
            atomType = line.substr(12,4);// needs to be this way until atomic numbers are assigned
            //resiname = line.substr(17,3);

            if (boost::regex_search(atomType, ifOAtom)){
                hydration.emplace_back(Coords(x-centering_vector->x,y-centering_vector->y, z-centering_vector->z, atomType, 1.0, "HOH", count));
                count++;
            }
        }
    }

    totalwaters = (unsigned int)hydration.size();
    totalwatersInExcludedVolume = hydration.size() - totalwaters;

    SASTOOLS_UTILS_H::logger("TOTAL WATERS IN HYDRATION", std::to_string(totalwaters));
}

float Waters::getForwardScatteringOfHydrationModel(){
    double sumOcc = 0.0;
    for(unsigned int n=0; n < totalwaters; n++){
        sumOcc += hydration[n].occ;
    }
    return (float)sumOcc*functions_WETSAXS::asf_at_q_zero(99);
}


void Waters::hydrateResidueDirect(std::string residue, int resid, int atomsInResidue, unsigned int startAt,
                                  const std::string * atomType, const float *xpos, const float *ypos, const float *zpos){
    // create target and reference vectors
    // reference is already centered
    std::string tempAtom;

    if ( sideChains.find(residue) == sideChains.end() ) { // not found
//    if ( true ) { // not found
            // If residue is not found, hydrate using tetrahedron model
            // create tetrahedron around the atom, randomly rotate
            // do this for each atom in residue
            // remove clashes

        boost::algorithm::trim(residue);
        std::string note = "RESID " + std::to_string(resid) + " => _" + residue+"_";
        SASTOOLS_UTILS_H::logger("PSEUDO HYDRATING", note);

        boost::regex ifCarbon("^C[0-9]+?"); // match any character

        // collect atoms in the unknown residue, center, determine size of bounding box, increase size and form a rectangular grid
//        float sum_x=0, sum_y=0, sum_z=0;
        float at_x, at_y, at_z;
        // float min_x=FLT_MAX, max_x = -FLT_MAX;
        // float min_y=FLT_MAX, max_y = -FLT_MAX;
        // float min_z=FLT_MAX, max_z = -FLT_MAX;

        // random rotation
//        std::random_device rd;
//        std::mt19937 gen(rd());
//        auto convert = (float)(M_PI/180.0f);
//        std::uniform_int_distribution<int> randomIndex(0,360); // guaranteed unbiased
//        std::uniform_int_distribution<int> randomBeta(0,180); // guaranteed unbiased

        std::vector <vector3> beads;

        auto * pRes = SASTOOLS_RESIDUES_H::Residues::getResidue(residue);

        // remove waters too close and within range
        // needs to be van der Waals radii plus 1.4 Angstrom (radius of water?)
        for (int i=0; i<atomsInResidue; i++) {

            tempAtom = atomType[startAt + i];

            // if atom type does not match C # where # can be any integer
            if (!boost::regex_match(tempAtom, ifCarbon)){

                float radius = pRes->getAtom(tempAtom)->getRadii();

                // diameter of water is 2.8 Angstroms -
                // water molecule should 1.4
                // so water should be 1.4 + VDW from atom

                float a_edge = (1.4f + radius)/sqrtf(6.0)*4.0f;
                float half_a_edge = 0.5f*a_edge;
                float x_tetra = a_edge*sqrtf(3.0)/3.0f;
                float z_tetra = a_edge*sqrtf(6.0)/3.0f;
                float half_x_tetra = 0.5f*x_tetra;

                std::vector<vector3> tetrahedron(4);
                std::vector<vector3> transformed(4);
                tetrahedron[0] = vector3(x_tetra,0,0);
                tetrahedron[1] = vector3(-half_x_tetra,half_a_edge,0);
                tetrahedron[2] = vector3(half_x_tetra,-half_a_edge,0);
                tetrahedron[3] = vector3(0,0,z_tetra);

                at_x = xpos[startAt + i];
                at_y = ypos[startAt + i];
                at_z = zpos[startAt + i];
                vector3 cvec = vector3(at_x, at_y, at_z);

                // rotate and translate
                float alpha_x = randomIndex(gen)*convert;
                float gamma_z = randomIndex(gen)*convert;
                float beta_y = randomBeta(gen)*convert;

                rotate(alpha_x, beta_y, gamma_z,  tetrahedron, transformed, cvec);

                for(auto & tetra : transformed){
                    beads.emplace_back(vector3(tetra));
                }
            }
        }

        int totalBeads = beads.size();
        float rlimit = 2.81*2.81;

        // move waters that are too close to molecule
        for (int i=0; i<atomsInResidue; i++) {

            tempAtom = atomType[startAt + i];
            at_x = xpos[startAt + i];
            at_y = ypos[startAt + i];
            at_z = zpos[startAt + i];
            vector3 cvec = vector3(at_x, at_y, at_z);

            float radius = pRes->getAtom(tempAtom)->getGaussianRadii() + 1.4;
            rlimit = radius*radius*0.9;

            for(int w=0; w<totalBeads; w++){ // sort the waters.  waters too close will be moved to end of the vector
                auto pW = &beads[w];

                if ((cvec - vector3(pW->x, pW->y, pW->z)).sqlength() < rlimit ){ // too close for hydrogen bond
                    // reject
                    std::iter_swap(beads.begin()+w, beads.begin()+(totalBeads-1));
                    totalBeads--;
                    w--; // check newly swapped bead on next round
                }

                if (w >= totalBeads){
                    break;
                }
            }
        }

        // now merge any waters that are too close to each other
        float limit = 2.65*2.65;
        for(int w=0; w<totalBeads; w++){ // sort the waters.  waters too close will be moved to end of the vector
            auto pCurrent = &beads[w];
            vector3 tvec = vector3(*pCurrent);

            int count = 1;
            for(int k=w+1; k<totalBeads; k++){
                auto pNext = &beads[k];
                if ( ( *pNext - *pCurrent).sqlength() < limit){
                    //*pCurrent = ( *pNext + *pCurrent)*0.5;
                    tvec += *pNext;
                    std::iter_swap(beads.begin()+k, beads.begin()+(totalBeads-1));
                    count += 1;
                    totalBeads--;
                    k--;
                }
            }
            *pCurrent = tvec/(float)count;

            if (w >= totalBeads){
                break;
            }
        }

        for (unsigned int i=0; i < totalBeads; i++){
            auto & bead = beads[i];
            hydration.emplace_back( Coords(bead.x, bead.y, bead.z, "O", 1.0) );
        }

//        std::string nameOf = "model.pdb";
//        FILE * pFile = fopen(nameOf.c_str(), "w");
//        fprintf(pFile, "REMARK 265                   EDGE RADIUS : %.3f\n", 2.8);
//
//        fprintf(pFile, "TER\n");
//        // print beads
//        for (unsigned int i=0; i < totalBeads; i++){
//            auto & bead = beads[i];
//            int residue_index = (i + 1);
//            SASTOOLS_UTILS_H::printAtomLine(pFile, i+1, "B", residue_index, bead.x, bead.y, bead.z );
//        }
//        fclose(pFile);
//        fprintf(pFile, "TER\n");
//        fprintf(pFile, "END\n");
//        exit(0);

    } else { // residue in library hydrate with known hydration

        std::vector<float> tar_x = std::vector<float>();
        std::vector<float> tar_y = std::vector<float>();
        std::vector<float> tar_z = std::vector<float>();
        std::vector<float> ref_x = std::vector<float>();
        std::vector<float> ref_y = std::vector<float>();
        std::vector<float> ref_z = std::vector<float>();

        // get all atoms in the sideChain, includes backbone
        auto findIt = sideChains.find(residue);
        std::string centering_atom = centeringAtoms.find(residue)->second;

        // need to add backbone for RNA - center on C4' or protein centered on CA
        // && (residue == "DG" || residue == "DC")

        if (findIt != sideChains.end()){
            auto psidechain = &(findIt->second); // std::unordered_map<std::string, Coords> >sideChains
            const float * at_x, * at_y, * at_z;
            unsigned int centering_index=startAt;

            for (int i=0; i< atomsInResidue; i++) {

                tempAtom = atomType[startAt + i];
                boost::algorithm::trim (tempAtom);

                if (centering_atom == tempAtom){
                    centering_index = startAt + i;
                }

                auto pAtom = psidechain->find(tempAtom);
                if (pAtom != psidechain->end()){ //
                    // DA, DT, DG, DC = P, OP1, OP2 and sugar atoms
                    // Protein residues N, CA, CB, C, O
                    Coords *ptempCoord = &(*pAtom).second;
                    at_x = &xpos[startAt + i];
                    at_y = &ypos[startAt + i];
                    at_z = &zpos[startAt + i];

                    tar_x.push_back(*at_x); // atoms in residue in structure
                    tar_y.push_back(*at_y);
                    tar_z.push_back(*at_z);

                    ref_x.push_back(ptempCoord->x); // atoms in residue in library
                    ref_y.push_back(ptempCoord->y);
                    ref_z.push_back(ptempCoord->z);
                }
            }

            // create alignment vectors - center residue on the centering_atom
            const float * aveX = &xpos[centering_index];
            const float * aveY = &ypos[centering_index];
            const float * aveZ = &zpos[centering_index];

            // center the residue in the structure
            auto totalAtomsInResidue = tar_x.size() ;//(int)sccount;

            for(int i=0; i < totalAtomsInResidue;i++){
                tar_x[i] = tar_x[i] - *aveX;
                tar_y[i] = tar_y[i] - *aveY;
                tar_z[i] = tar_z[i] - *aveZ;
            }

            // rotate library residue into structure residue
            setRotationMatrix(totalAtomsInResidue, &ref_x[0], &ref_y[0], &ref_z[0], &tar_x[0], &tar_y[0], &tar_z[0]);

            //cout << "____HYDRATE___" << endl;
            int totalToAdd = placeWaters(residue, resid, *aveX, *aveY, *aveZ);
            // copy into vector
            for(int m=0; m<totalToAdd; m++){
                hydration.emplace_back( Coords(tempWaters[m]) );
            }
        }

        // if nucleic, need to hydrate backbone separately from residue
        std::string bck = assignBackbone(residue);

        if (bck.size() > 1) {
            ref_x.clear();
            ref_y.clear();
            ref_z.clear();
            tar_x.clear();
            tar_y.clear();
            tar_z.clear();

             findIt = sideChains.find(bck);
            centering_atom = centeringAtoms.find(bck)->second;

            auto pBackbone = &(findIt->second); // std::unordered_map<std::string, Coords> >sideChains
            // backbone will not contain the centering atom

            unsigned int centering_index=startAt, index;
            float * pTar_x, * pTar_y, * pTar_z;
            for (int i=0; i< atomsInResidue; i++) {

                index = startAt + i;
                tempAtom = atomType[index];
                boost::algorithm::trim (tempAtom);

                if (centering_atom == tempAtom){ // find index of the centering atom for backbone
                    centering_index = index;
                }

                auto pBackboneAtom = pBackbone->find(tempAtom);
                if (pBackboneAtom != pBackbone->end()){ //
                    // DA, DT, DG, DC = P, OP1, OP2 and sugar atoms
                    // Protein residues N, CA, CB, C, O
                    Coords *ptempCoord = &(*pBackboneAtom).second;

                    tar_x.push_back(xpos[index]); // atoms in residue in structure
                    tar_y.push_back(ypos[index]);
                    tar_z.push_back(zpos[index]);

                    ref_x.push_back(ptempCoord->x); // atoms in residue in library
                    ref_y.push_back(ptempCoord->y);
                    ref_z.push_back(ptempCoord->z);
                }
            }

            // create alignment vectors
            const float * aveX = &xpos[centering_index];
            const float * aveY = &ypos[centering_index];
            const float * aveZ = &zpos[centering_index];

            // center the residue in the structure
            auto totalAtomsInResidue = tar_x.size() ;//(int)sccount;

            for(int i=0; i < totalAtomsInResidue;i++){
                pTar_x = &tar_x[i];
                pTar_y = &tar_y[i];
                pTar_z = &tar_z[i];
                *pTar_x -= *aveX;
                *pTar_y -= *aveY;
                *pTar_z -= *aveZ;
                // tar_x[i] = tar_x[i] - *aveX;
                // tar_y[i] = tar_y[i] - *aveY;
                // tar_z[i] = tar_z[i] - *aveZ;
            }

            // rotate library residue into structure residue
            setRotationMatrix(totalAtomsInResidue, &ref_x[0], &ref_y[0], &ref_z[0], &tar_x[0], &tar_y[0], &tar_z[0]);
            int totalToAdd = placeWaters(bck, resid, *aveX, *aveY, *aveZ);
            // copy into vector
            for(int m=0; m<totalToAdd; m++){
                hydration.emplace_back( Coords(tempWaters[m]) );
            }
        }

    }
}


// reset array elements to zero
void Waters::resetTempWaters(){
    for (int i=0; i<totalTempWaters; i++){
        tempWaters[i].x = 0;
        tempWaters[i].y = 0;
        tempWaters[i].z = 0;
        tempWaters[i].type = "";
        tempWaters[i].occ = 0;
    }
}

/**
 * Determine Rotation Matrix for rotation of the reference residue with respect to the target residue (from input PDB)
 * coordinates for ref and tar must be centered
 *
 * @param totalAtomsInResidue
 * @param ref_x
 * @param ref_y
 * @param ref_z
 * @param tar_x
 * @param tar_y
 * @param tar_z
 */
void Waters::setRotationMatrix(int totalAtomsInResidue, float * ref_x, float * ref_y, float * ref_z, float * tar_x, float * tar_y, float * tar_z) {
    /*
     * for the 3x3 matrix
     * rotation of P into Q
     * calculation P^T*Q
     *
     * rotate reference into target
     *
     */
    float a11 =0;
    float a12 =0;
    float a13 =0;
    float a21 =0, a22=0, a23=0;
    float a31 =0, a32=0, a33=0;

    for(int i=0; i<totalAtomsInResidue; i++){
        a11 += ref_x[i]*tar_x[i];
        a12 += ref_x[i]*tar_y[i];
        a13 += ref_x[i]*tar_z[i];

        a21 += ref_y[i]*tar_x[i];
        a22 += ref_y[i]*tar_y[i];
        a23 += ref_y[i]*tar_z[i];

        a31 += ref_z[i]*tar_x[i];
        a32 += ref_z[i]*tar_y[i];
        a33 += ref_z[i]*tar_z[i];
    }

    float u11, u12, u13, u21, u22, u23, u31, u32, u33;
    float s11, s12, s13, s21, s22, s23, s31, s32, s33;
    float v11, v12, v13, v21, v22, v23, v31, v32, v33;

    // M = USV*
    svd3::svd(a11, a12, a13, a21, a22, a23, a31, a32, a33,
        u11, u12, u13, u21, u22, u23, u31, u32, u33,
        s11, s12, s13, s21, s22, s23, s31, s32, s33,
        v11, v12, v13, v21, v22, v23, v31, v32, v33);

    float d11, d12, d13, d21, d22, d23, d31, d32, d33;

    // calculate det(VU*) where * represents transpose
    svd3::multAB(v11, v12, v13, v21, v22, v23, v31, v32, v33,
           u11, u21, u31, u12, u22, u32, u13, u23, u33,
           d11, d12, d13, d21, d22, d23, d31, d32, d33);

    // calculate sign of determinant
    float det = d11*d22*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31 - d11*d23*d32 - d12*d21*d33;
    float d = (det < 0) ? -1.0f : 1.0f;

    svd3::multAB(v11, v12, v13, v21, v22, v23, v31, v32, v33,
           1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, d,
           d11, d12, d13, d21, d22, d23, d31, d32, d33);

    svd3::multAB(d11, d12, d13, d21, d22, d23, d31, d32, d33,
           u11, u21, u31, u12, u22, u32, u13, u23, u33,
           rotation[0], rotation[1], rotation[2],
           rotation[3], rotation[4], rotation[5],
           rotation[6], rotation[7], rotation[8]);

    // rotate reference and print out
//    float new_x, new_y, new_z;
//    int index = 1;
//    for(int i=0;i<totalAtomsInResidue; i++){
//        float x_val = ref_x[i];
//        float y_val = ref_y[i];
//        float z_val = ref_z[i];
//
//        new_x = rotation[0]*x_val + rotation[1]*y_val + rotation[2]*z_val;
//        new_y = rotation[3]*x_val + rotation[4]*y_val + rotation[5]*z_val;
//        new_z = rotation[6]*x_val + rotation[7]*y_val + rotation[8]*z_val;
//
//        printf("%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", index," CA ", "TRP", "D", "4", new_x, new_y, new_z );
//        index++;
//    }

}

// grab all the waters associated with the residue and perform rotation and translation
// this will populate tempWaters
int Waters::placeWaters(std::string residue, int resid, float delX, float delY, float delZ){

    std::vector<Coords> * pwater = &waters[residue];
    auto total = (unsigned int)pwater->size();

    float x_val, y_val, z_val;// new_x, new_y, new_z;

    resetTempWaters();

    for(unsigned int i=0; i<total; i++){ // for each water
        Coords * pC = &(*pwater)[i]; // original water
        x_val = pC->x;
        y_val = pC->y;
        z_val = pC->z;

        Coords * pT = &tempWaters[i];

        pT->x = rotation[0]*x_val + rotation[1]*y_val + rotation[2]*z_val + delX;
        pT->y = rotation[3]*x_val + rotation[4]*y_val + rotation[5]*z_val + delY;
        pT->z = rotation[6]*x_val + rotation[7]*y_val + rotation[8]*z_val + delZ;
        pT->type = "O";
        pT->occ = pC->occ;
        pT->resid = resid;
        pT->resname = residue;
    }

    return total;
    //writeWaterCoords(residue, total);
}


/**
 * write the waters used to determine the hydrated particle
 * chain W are waters that are included in the SAXS model
 * chain Y are the excluded waters that are too close
 * @param name
 */
void Waters::writeWatersToFile(std::string name) {
    std::string nameOf = name+"_w.pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    int count=1;
    for(unsigned int i=0; i < totalwaters; i++){
        Coords * pCoord = &hydration[i];
        std::string resid = std::to_string(count);
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  %3.2f100.00\n", "ATOM", count, " O  ", "HOH", "W",
                resid.c_str(),
                pCoord->x,
                pCoord->y,
                pCoord->z,
                pCoord->occ);
        count++;
    }

    for(unsigned int i=totalwaters; i < hydration.size(); i++){
        Coords * pCoord = &hydration[i];
        std::string resid = std::to_string(count);
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", count, " O  ", "HOH", "Y", resid.c_str(), pCoord->x, pCoord->y, pCoord->z);
        count++;
    }

    fclose(pFile);
}


/**
 * write the waters used to determine the hydrated particle
 * chain W are waters that are included in the SAXS model
 * chain Y are the excluded waters that are too close
 * @param name
 */
void Waters::translateAndWriteWatersToFile(std::string name, const vector3 * pCenVec) {
    std::string nameOf = name+"_w.pdb";
    FILE * pFile = fopen(nameOf.c_str(), "w");

    int count=1;
    for(unsigned int i=0; i < totalwaters; i++){
        Coords * pCoord = &hydration[i];
        std::string resid = std::to_string(count);

        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  %3.2f100.00\n", "ATOM", count, " O  ", "HOH", "W",
                resid.c_str(),
                pCoord->x + pCenVec->x,
                pCoord->y + pCenVec->y,
                pCoord->z + pCenVec->z,
                pCoord->occ);

        count++;
    }

    // 9999 is largest number for resid
    int index = 0;
    std::string chains[3] ={ "X", "Y", "Z"};
    std::string chainID = chains[index];

    for(unsigned int i=totalwaters; i < hydration.size(); i++){
        Coords * pCoord = &hydration[i];
        std::string resid = std::to_string(count);
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", count, " O  ", "HOH", chainID.c_str(),
                resid.c_str(),
                pCoord->x + pCenVec->x,
                pCoord->y + pCenVec->y,
                pCoord->z + pCenVec->z);
        count++;

        if (count % 9999 == 0 ){ // reset count and use new chain ID
            index += 1;
            chainID = chains[index];
            count = 1;
        }

    }

    fclose(pFile);
}


/**
 * must only be used after the atomistic model is hydrated
 */
void Waters::createSphericalCoordinateOfHydration() {

    if (totalwaters == 0)
        throw std::invalid_argument("no waters, run hydration first");

    rvalues.resize(totalwaters);
    phis.resize(totalwaters);
    thetas.resize(totalwaters);
    float xval, yval, zval, rval;

//    std::clock_t startTime = std::clock();
    for(unsigned int n=0; n < totalwaters; n++){
        Coords * pCoord = &hydration[n];

        xval = pCoord->x;
        yval = pCoord->y;
        zval = pCoord->z;

        rval = sqrtf(xval*xval + yval*yval +zval*zval);

        rvalues[n] = rval;
        thetas[n] = (rval > 0.0f) ? acosf(zval/rval) : 0.0f; // guard: acos(0/0) is NaN
        phis[n] = (xval == 0 || xval == 0.0f) ? 0.0f : atan2f(yval, xval); // need to check if this is a good way to test for 0
    }

//    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//    std::cout << "THREADS elapsed time= s :: "  << runtime << std::endl;

}


/*
 * Partial amplitudes of the hydration sites, two sets sharing one geometric sum
 *
 *   S_lm(q) = sum_n occ_n * j_l(q r_n) * Y_lm(theta_n, phi_n)
 *
 *   A_w   = f_w(q)  * S_lm   bound water: real molecular form factor of water (ASF 99)
 *   A_w0  = g_w(q)  * S_lm   bulk water displaced by the same sites (buffer side):
 *                            smooth Gaussian blob of 29.9 A^3, g_w(q) = 0.334*29.9*exp(-q^2 29.9^(2/3)/4pi)
 *
 * In Fit, A_xw0 = A_x (AtomisticModel excluded volume) + A_w0, and
 *   I(q) = < | A_p + c_o A_w - c_x A_xw0 |^2 >
 * At q = 0 both f_w and g_w are ~10 electrons, so the shell contrast is (c_o - c_x).
 * Both sets are stored with a POSITIVE sign; Babinet's minus sign is applied in Fit.
 *
 * Index of A_lm(q_i): (lmax+1)^2 * q_index + l*l + l + m
 */
void Waters::calculatePartialAmplitudes(unsigned int lmax,
                                        unsigned int totalqvalues,
                                        std::vector < float > & qvalues,
                                        bool recalculate){

    // if qvalues changed, need to recalculate spherical bessels
    if (recalculate){
        sbj.recalculate(totalqvalues, qvalues, &rvalues);
    } else {
        this->createSphericalHarmonics(lmax);
        this->createSphericalBessels(totalqvalues, qvalues);
    }

    const unsigned int partialsPerQ = (lmax+1)*(lmax+1);
    almRealSum.resize(partialsPerQ*totalqvalues);
    almImagSum.resize(partialsPerQ*totalqvalues);
    almRealDisplacedSum.resize(partialsPerQ*totalqvalues);
    almImagDisplacedSum.resize(partialsPerQ*totalqvalues);

    if (totalwaters == 0) {
        std::fill(almRealSum.begin(), almRealSum.end(), 0.0f);
        std::fill(almImagSum.begin(), almImagSum.end(), 0.0f);
        std::fill(almRealDisplacedSum.begin(), almRealDisplacedSum.end(), 0.0f);
        std::fill(almImagDisplacedSum.begin(), almImagDisplacedSum.end(), 0.0f);
        return;
    }

    float * const ptrR  = almRealSum.data();
    float * const ptrI  = almImagSum.data();
    float * const ptrR0 = almRealDisplacedSum.data();
    float * const ptrI0 = almImagDisplacedSum.data();

    const float * pSBJ = sbj.getpSphericalBessels();
    const float * pRealYlm = she.getpDataYLMReal();
    const float * pImagYlm = she.getpDataYLMImag();

    std::vector<float> occ(totalwaters);
    for (unsigned int n = 0; n < totalwaters; n++){
        occ[n] = hydration[n].occ;
    }

    const float rho0 = 0.334f, vWater = 29.9f;
    const float waterWidth = std::cbrt(vWater*vWater)/(4.0f*(float)M_PI);

    std::clock_t startTime = std::clock();
    unsigned int q_index = 0;

    for(auto & qvalue : qvalues){

        const float fw = functions_WETSAXS::asf(99, qvalue);                        // bound water
        const float gw = rho0*vWater*std::exp(-qvalue*qvalue*waterWidth);         // displaced bulk water

        const unsigned int ll_index = partialsPerQ*q_index;

        for(int l=0; l<=(int)lmax; l++){

            const unsigned long bessel_index = (unsigned long)(lmax+1)*totalwaters*q_index + (unsigned long)totalwaters*l;
            const unsigned long base_index = (unsigned long)(l*l)*totalwaters;
            const unsigned int partials_index = ll_index + l*l + l;

            // m = 0 : imaginary part of Y_l0 is zero
            {
                double sumR = 0.0;
                unsigned long ylm_index = base_index + (unsigned long)l*totalwaters;

                for (unsigned int n = 0; n < totalwaters; n++, ylm_index++) {
                    sumR += occ[n]*pSBJ[bessel_index + n]*pRealYlm[ylm_index];
                }

                ptrR[partials_index]  = fw*(float)sumR;  ptrI[partials_index]  = 0.0f;
                ptrR0[partials_index] = gw*(float)sumR;  ptrI0[partials_index] = 0.0f;
            }

            // m > 0, and m < 0 by symmetry of a real density: A_{l,-m} = (-1)^m conj(A_{l,m})
            for(int m=1; m<=l; m++){

                double sumR = 0.0, sumI = 0.0;
                unsigned long ylm_index = base_index + (unsigned long)(l+m)*totalwaters;

                for(unsigned int n=0; n<totalwaters; n++, ylm_index++){
                    const float b = occ[n]*pSBJ[bessel_index + n];
                    sumR += b*pRealYlm[ylm_index];
                    sumI += b*pImagYlm[ylm_index];
                }

                const float sR = (float)sumR, sI = (float)sumI;
                unsigned int idx = partials_index + m;
                ptrR[idx]  = fw*sR;  ptrI[idx]  = fw*sI;
                ptrR0[idx] = gw*sR;  ptrI0[idx] = gw*sI;

                idx = partials_index - m;
                if (m & 1){ // odd m: flip sign of the real part
                    ptrR[idx]  = -fw*sR;  ptrI[idx]  = fw*sI;
                    ptrR0[idx] = -gw*sR;  ptrI0[idx] = gw*sI;
                } else {    // even m: flip sign of the imaginary part
                    ptrR[idx]  = fw*sR;   ptrI[idx]  = -fw*sI;
                    ptrR0[idx] = gw*sR;   ptrI0[idx] = -gw*sI;
                }
            }
        }

        q_index++;
    }

    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
    logger("Total TIME (Waters)", formatNumber((float)runtime,8));
}


/*
 * must only be used after hydration has been set and spherical coordinates
 *
 */
void Waters::populateAmplitudesExcludedVolume(unsigned lmax, unsigned int totalqvalues, std::vector<float> & qvalues){

    if (rvalues.empty())
        throw std::invalid_argument("spherical coordinates not initialized");

    this->lmax = lmax;
    she = SphericalHarmonics(lmax, totalwaters);
    she.populateSHETable(thetas.data(), phis.data());
    sbj = SphericalBessels(lmax, totalqvalues, totalwaters, qvalues, &rvalues);

}

// thetas and phis need to be prepopulated
void Waters::createSphericalHarmonics(unsigned int lmax) {

    this->lmax = lmax;
    she = SphericalHarmonics(lmax, totalwaters);
    she.populateSHETable(thetas.data(), phis.data());

}

/*
 * must be run after SphericalHarmonics
 */
void Waters::createSphericalBessels(unsigned int totalqvalues, std::vector < float > & qvalues) {

    sbj = SphericalBessels(lmax, totalqvalues, totalwaters, qvalues, &rvalues);

}

std::string Waters::assignBackbone(std::string residue) {

    if (!(residue == "rA" || residue == "rG" || residue == "rU" || residue == "rC" || residue == "rI" ||
        residue == "DA" || residue == "DG" || residue == "DT" || residue == "DC")) {
        return "PROTEIN";
    }

    if (residue == "rA" || residue == "rG" || residue == "rU" || residue == "rC" || residue == "rI") {
        return "RNA";
    }

    if (residue == "DA" || residue == "DG" || residue == "DT" || residue == "DC") {
        return "DNA";
    }

    return "";
}

void Waters::rotate(float alpha, float beta, float gamma, std::vector<vector3> & tetrahedron,
                    std::vector<vector3> & transformedCoordinates, vector3 & translateToVector) {

    float cosAlpha = cos(alpha);
    float sinAlpha = sin(alpha);
    float cosBeta = cos(beta);
    float sinBeta = sin(beta);
    float cosGamma = cos(gamma);
    float sinGamma = sin(gamma);

    float x1 = cosGamma*cosBeta;
    float x2 = sinGamma*cosBeta;
    float x3 = -sinBeta;
    float y1 = -sinGamma*cosAlpha + cosGamma*sinBeta*sinAlpha;
    float y2 = cosGamma*cosAlpha + sinGamma*sinBeta*sinAlpha;
    float y3 = cosBeta*sinAlpha;
    float z1 = sinGamma*sinAlpha + cosGamma*sinBeta*cosAlpha;
    float z2 = -cosGamma*sinAlpha + sinGamma*sinBeta*cosAlpha;
    float z3 = cosBeta*cosAlpha;

    vector3 * tempVec1, *trans;

    for(unsigned int i=0; i<4; i++){
        tempVec1 = &tetrahedron[i];
        trans = &transformedCoordinates[i];

        (*trans).x = x1*tempVec1->x + y1*tempVec1->y + z1*tempVec1->z + translateToVector.x;
        (*trans).y = x2*tempVec1->x + y2*tempVec1->y + z2*tempVec1->z + translateToVector.y;
        (*trans).z = x3*tempVec1->x + y3*tempVec1->y + z3*tempVec1->z + translateToVector.z;
    }
}


