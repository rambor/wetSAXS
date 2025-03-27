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
    sideChains["ALA"].insert({"N",Coords(-1.1640f,0.7258,-0.7960f,"N",1.0)});
    sideChains["ALA"].insert({"CA",Coords(-0.0600f,-0.2272f,-0.6820f,"CA",1.0)});
    sideChains["ALA"].insert({"C",Coords(-0.0560f,-0.8512f,0.6900,"C",1.0)});
    sideChains["ALA"].insert({"O",Coords(-0.0200f,-0.1132f,1.6960,"O",1.0)});
    sideChains["ALA"].insert({"CB",Coords(1.3000,0.4658,-0.9080,"CB",1.0)});
    sideChains["ARG"].insert({"N",Coords(0.1116,0.6292,3.9827,"N",1.0)});
    sideChains["ARG"].insert({"CA",Coords(-0.3734,0.5782,2.5887,"CA",1.0)});
    sideChains["ARG"].insert({"C",Coords(-1.5014,-0.4348,2.4027,"C",1.0)});
    sideChains["ARG"].insert({"O",Coords(-2.5044,-0.1758,1.6977,"O",1.0)});
    sideChains["ARG"].insert({"CB",Coords(0.7716,0.2722,1.6077,"CB",1.0)});
    sideChains["ARG"].insert({"CG",Coords(0.3186,-0.0578,0.1957,"CG",1.0)});
    sideChains["ARG"].insert({"CD",Coords(1.5136,-0.3258,-0.7193,"CD",1.0)});
    sideChains["ARG"].insert({"NE",Coords(1.0156,-0.8408,-1.9943,"NE",1.0)});
    sideChains["ARG"].insert({"CZ",Coords(0.4216,-0.1118,-2.9413,"CZ",1.0)});
    sideChains["ARG"].insert({"NH1",Coords(0.2796,1.1982,-2.7983,"NH1",1.0)});
    sideChains["ARG"].insert({"NH2",Coords(-0.0534,-0.7308,-4.0223,"NH2",1.0)});
    sideChains["ASN"].insert({"N",Coords(1.2974,1.3450,-0.7467,"N",1.0)});
    sideChains["ASN"].insert({"CA",Coords(-0.0166,0.8460,-0.3387,"CA",1.0)});
    sideChains["ASN"].insert({"C",Coords(-0.9646,1.9630,0.0763,"C",1.0)});
    sideChains["ASN"].insert({"O",Coords(-2.1546,1.9320,-0.2758,"O",1.0)});
    sideChains["ASN"].insert({"CB",Coords(0.1304,-0.1990,0.7752,"CB",1.0)});
    sideChains["ASN"].insert({"CG",Coords(0.4834,-1.5700,0.2482,"CG",1.0)});
    sideChains["ASN"].insert({"OD1",Coords(0.0844,-1.9730,-0.8327,"OD1",1.0)});
    sideChains["ASN"].insert({"ND2",Coords(1.1404,-2.3440,1.0943,"ND2",1.0)});
    sideChains["ASP"].insert({"N",Coords(-1.2417,-0.4275,-1.2545,"N",1.0)});
    sideChains["ASP"].insert({"CA",Coords(0.1383,-0.4995,-0.7225,"CA",1.0)});
    sideChains["ASP"].insert({"C",Coords(0.7683,-1.8815,-0.8365,"C",1.0)});
    sideChains["ASP"].insert({"O",Coords(0.1183,-2.8975,-0.5645,"O",1.0)});
    sideChains["ASP"].insert({"CB",Coords(0.1813,-0.0125,0.7315,"CB",1.0)});
    sideChains["ASP"].insert({"CG",Coords(0.0463,1.5075,0.8485,"CG",1.0)});
    sideChains["ASP"].insert({"OD1",Coords(0.0303,2.2045,-0.1895,"OD1",1.0)});
    sideChains["ASP"].insert({"OD2",Coords(-0.0407,2.0065,1.9875,"OD2",1.0)});
    sideChains["CYS"].insert({"N",Coords(0.9042,-1.4382,0.9132,"N",1.0)});
    sideChains["CYS"].insert({"CA",Coords(0.4992,-0.0842,0.5522,"CA",1.0)});
    sideChains["CYS"].insert({"C",Coords(-0.9988,0.1158,0.7882,"C",1.0)});
    sideChains["CYS"].insert({"O",Coords(-1.7768,-0.5952,0.1482,"O",1.0)});
    sideChains["CYS"].insert({"CB",Coords(0.8222,0.1518,-0.9138,"CB",1.0)});
    sideChains["CYS"].insert({"SG",Coords(0.5502,1.8498,-1.4878,"SG",1.0)});
    sideChains["GLN"].insert({"N",Coords(-0.1012,-0.1694,2.0707,"N",1.0)});
    sideChains["GLN"].insert({"CA",Coords(-0.8012,-0.0514,0.7937,"CA",1.0)});
    sideChains["GLN"].insert({"C",Coords(-2.3272,-0.1464,0.9047,"C",1.0)});
    sideChains["GLN"].insert({"O",Coords(-3.0062,-0.1364,-0.1073,"O",1.0)});
    sideChains["GLN"].insert({"CB",Coords(-0.2542,-1.0954,-0.2033,"CB",1.0)});
    sideChains["GLN"].insert({"CG",Coords(1.2598,-0.9394,-0.4413,"CG",1.0)});
    sideChains["GLN"].insert({"CD",Coords(1.6238,0.4596,-0.9183,"CD",1.0)});
    sideChains["GLN"].insert({"OE1",Coords(1.1618,0.9066,-1.9733,"OE1",1.0)});
    sideChains["GLN"].insert({"NE2",Coords(2.4448,1.1726,-0.1253,"NE2",1.0)});
    sideChains["GLU"].insert({"N",Coords(2.0491,-0.8326,-1.1278,"N",1.0)});
    sideChains["GLU"].insert({"CA",Coords(0.5951,-1.0666,-1.1118,"CA",1.0)});
    sideChains["GLU"].insert({"C",Coords(0.0231,-1.1156,-2.5358,"C",1.0)});
    sideChains["GLU"].insert({"O",Coords(-1.1909,-1.0566,-2.7218,"O",1.0)});
    sideChains["GLU"].insert({"CB",Coords(-0.0919,0.0164,-0.3008,"CB",1.0)});
    sideChains["GLU"].insert({"CG",Coords(0.2891,-0.0566,1.1832,"CG",1.0)});
    sideChains["GLU"].insert({"CD",Coords(-0.3719,1.0654,2.0062,"CD",1.0)});
    sideChains["GLU"].insert({"OE1",Coords(0.1361,2.2234,1.9852,"OE1",1.0)});
    sideChains["GLU"].insert({"OE2",Coords(-1.4379,0.8224,2.6232,"OE2",1.0)});
    sideChains["GLY"].insert({"N",Coords(-1.6083,0.1225,-0.4023,"N",1.0)});
    sideChains["GLY"].insert({"CA",Coords(-0.3422,-0.5765,-0.5853,"CA",1.0)});
    sideChains["GLY"].insert({"C",Coords(0.8478,0.2265,-0.1123,"C",1.0)});
    sideChains["GLY"].insert({"O",Coords(1.1028,0.2275,1.0998,"O",1.0)});
    sideChains["HIS"].insert({"N",Coords(0.8684,1.2451,-1.3149,"N",1.0)});
    sideChains["HIS"].insert({"CA",Coords(1.2274,1.1471,0.1021,"CA",1.0)});
    sideChains["HIS"].insert({"C",Coords(2.7174,1.2321,0.2951,"C",1.0)});
    sideChains["HIS"].insert({"O",Coords(3.4644,0.6131,-0.4249,"O",1.0)});
    sideChains["HIS"].insert({"CB",Coords(0.7264,-0.1649,0.7441,"CB",1.0)});
    sideChains["HIS"].insert({"CG",Coords(-0.6746,-0.5669,0.3571,"CG",1.0)});
    sideChains["HIS"].insert({"ND1",Coords(-1.7776,0.2581,0.4891,"ND1",1.0)});
    sideChains["HIS"].insert({"CD2",Coords(-1.1626,-1.7529,-0.0999,"CD2",1.0)});
    sideChains["HIS"].insert({"CE1",Coords(-2.8676,-0.3889,0.0931,"CE1",1.0)});
    sideChains["HIS"].insert({"NE2",Coords(-2.5216,-1.6219,-0.2409,"NE2",1.0)});
    sideChains["ILE"].insert({"N",Coords(-0.1651,0.3650,-1.7521,"N",1.0)});
    sideChains["ILE"].insert({"CA",Coords(0.5189,0.2690,-0.4671,"CA",1.0)});
    sideChains["ILE"].insert({"C",Coords(1.6859,-0.7180,-0.6011,"C",1.0)});
    sideChains["ILE"].insert({"O",Coords(1.5049,-1.8220,-1.0861,"O",1.0)});
    sideChains["ILE"].insert({"CB",Coords(-0.4791,-0.1940,0.6349,"CB",1.0)});
    sideChains["ILE"].insert({"CG1",Coords(-1.3001,0.9870,1.1789,"CG1",1.0)});
    sideChains["ILE"].insert({"CG2",Coords(0.2439,-0.7170,1.8639,"CG2",1.0)});
    sideChains["ILE"].insert({"CD1",Coords(-2.0091,1.8300,0.2289,"CD1",1.0)});
    sideChains["LEU"].insert({"N",Coords(-1.8109,-0.5595,-0.3068,"N",1.0)});
    sideChains["LEU"].insert({"CA",Coords(-0.3899,-0.8955,-0.1888,"CA",1.0)});
    sideChains["LEU"].insert({"C",Coords(-0.2319,-2.0705,0.7382,"C",1.0)});
    sideChains["LEU"].insert({"O",Coords(-0.7119,-2.0485,1.8872,"O",1.0)});
    sideChains["LEU"].insert({"CB",Coords(0.4091,0.3015,0.3332,"CB",1.0)});
    sideChains["LEU"].insert({"CG",Coords(0.4201,1.5225,-0.5718,"CG",1.0)});
    sideChains["LEU"].insert({"CD1",Coords(1.3361,2.5965,0.0313,"CD1",1.0)});
    sideChains["LEU"].insert({"CD2",Coords(0.9791,1.1535,-1.9228,"CD2",1.0)});
    sideChains["LYS"].insert({"N",Coords(-1.8770,1.4420,1.2700,"N",1.0)});
    sideChains["LYS"].insert({"CA",Coords(-1.1030,0.2040,1.2830,"CA",1.0)});
    sideChains["LYS"].insert({"C",Coords(-0.2470,0.2160,2.5440,"C",1.0)});
    sideChains["LYS"].insert({"O",Coords(0.4870,1.1720,2.8220,"O",1.0)});
    sideChains["LYS"].insert({"CB",Coords(-0.2290,0.0820,0.0390,"CB",1.0)});
    sideChains["LYS"].insert({"CG",Coords(-1.0000,-0.3180,-1.2490,"CG",1.0)});
    sideChains["LYS"].insert({"CD",Coords(0.0090,-0.6300,-2.3730,"CD",1.0)});
    sideChains["LYS"].insert({"CE",Coords(1.3350,-1.2010,-1.7780,"CE",1.0)});
    sideChains["LYS"].insert({"NZ",Coords(2.6250,-0.9670,-2.5580,"NZ",1.0)});
    sideChains["MET"].insert({"N",Coords(-1.5350,1.3150,-0.3264,"N",1.0)});
    sideChains["MET"].insert({"CA",Coords(-1.2900,0.1670,0.5846,"CA",1.0)});
    sideChains["MET"].insert({"C",Coords(-1.4490,0.6100,2.0396,"C",1.0)});
    sideChains["MET"].insert({"O",Coords(-1.1510,1.7600,2.3906,"O",1.0)});
    sideChains["MET"].insert({"CB",Coords(0.1280,-0.3840,0.3846,"CB",1.0)});
    sideChains["MET"].insert({"CG",Coords(0.5180,-0.7680,-1.0464,"CG",1.0)});
    sideChains["MET"].insert({"SD",Coords(2.2530,-1.3070,-1.1254,"SD",1.0)});
    sideChains["MET"].insert({"CE",Coords(2.5260,-1.3930,-2.9014,"CE",1.0)});
    sideChains["PHE"].insert({"N",Coords(-0.0683,1.9033,-1.4177,"N",1.0)});
    sideChains["PHE"].insert({"CA",Coords(-0.7693,1.6893,-0.1377,"CA",1.0)});
    sideChains["PHE"].insert({"C",Coords(-1.1843,3.0113,0.4973,"C",1.0)});
    sideChains["PHE"].insert({"O",Coords(-0.4523,4.0003,0.3973,"O",1.0)});
    sideChains["PHE"].insert({"CB",Coords(0.0857,0.8453,0.8333,"CB",1.0)});
    sideChains["PHE"].insert({"CG",Coords(0.2457,-0.5797,0.3863,"CG",1.0)});
    sideChains["PHE"].insert({"CD1",Coords(-0.8303,-1.4677,0.4413,"CD1",1.0)});
    sideChains["PHE"].insert({"CD2",Coords(1.4687,-1.0317,-0.0977,"CD2",1.0)});
    sideChains["PHE"].insert({"CE1",Coords(-0.6763,-2.7917,0.0323,"CE1",1.0)});
    sideChains["PHE"].insert({"CE2",Coords(1.6217,-2.3467,-0.5127,"CE2",1.0)});
    sideChains["PHE"].insert({"CZ",Coords(0.5587,-3.2317,-0.4217,"CZ",1.0)});
    sideChains["PRO"].insert({"N",Coords(-1.0674,0.7484,-0.0024,"N",1.0)});
    sideChains["PRO"].insert({"CA",Coords(0.0706,0.4534,0.8926,"CA",1.0)});
    sideChains["PRO"].insert({"C",Coords(1.3416,0.9804,0.2916,"C",1.0)});
    sideChains["PRO"].insert({"O",Coords(1.5356,0.8334,-0.9164,"O",1.0)});
    sideChains["PRO"].insert({"CB",Coords(0.0756,-1.0796,0.9256,"CB",1.0)});
    sideChains["PRO"].insert({"CG",Coords(-0.4734,-1.5056,-0.3944,"CG",1.0)});
    sideChains["PRO"].insert({"CD",Coords(-1.4824,-0.4306,-0.7964,"CD",1.0)});
    sideChains["SER"].insert({"N",Coords(0.4822,0.6052,-1.4895,"N",1.0)});
    sideChains["SER"].insert({"CA",Coords(0.6842,0.4172,-0.0505,"CA",1.0)});
    sideChains["SER"].insert({"C",Coords(0.5632,-1.0448,0.4105,"C",1.0)});
    sideChains["SER"].insert({"O",Coords(0.1412,-1.9218,-0.3535,"O",1.0)});
    sideChains["SER"].insert({"CB",Coords(-0.3128,1.2982,0.7165,"CB",1.0)});
    sideChains["SER"].insert({"OG",Coords(-1.5578,0.6462,0.7665,"OG",1.0)});
    sideChains["THR"].insert({"N",Coords(0.0024,0.5187,-1.4923,"N",1.0)});
    sideChains["THR"].insert({"CA",Coords(0.0614,0.7137,-0.0463,"CA",1.0)});
    sideChains["THR"].insert({"C",Coords(1.5014,0.4577,0.4297,"C",1.0)});
    sideChains["THR"].insert({"O",Coords(2.3504,0.0357,-0.3643,"O",1.0)});
    sideChains["THR"].insert({"CB",Coords(-0.9406,-0.1963,0.7107,"CB",1.0)});
    sideChains["THR"].insert({"OG1",Coords(-0.5796,-1.5753,0.5337,"OG1",1.0)});
    sideChains["THR"].insert({"CG2",Coords(-2.3956,0.0457,0.2287,"CG2",1.0)});
    sideChains["TRP"].insert({"N",Coords(2.1504,2.3076,0.8213,"N",1.0)});
    sideChains["TRP"].insert({"CA",Coords(1.0854,1.9146,1.7013,"CA",1.0)});
    sideChains["TRP"].insert({"C",Coords(1.4564,2.1126,3.1693,"C",1.0)});
    sideChains["TRP"].insert({"O",Coords(0.6304,2.6086,3.9683,"O",1.0)});
    sideChains["TRP"].insert({"CB",Coords(0.6914,0.4656,1.4183,"CB",1.0)});
    sideChains["TRP"].insert({"CG",Coords(-0.1216,0.2876,0.1633,"CG",1.0)});
    sideChains["TRP"].insert({"CD1",Coords(-0.9616,1.1956,-0.4307,"CD1",1.0)});
    sideChains["TRP"].insert({"CD2",Coords(-0.2266,-0.9214,-0.6047,"CD2",1.0)});
    sideChains["TRP"].insert({"NE1",Coords(-1.5656,0.6346,-1.5477,"NE1",1.0)});
    sideChains["TRP"].insert({"CE2",Coords(-1.1436,-0.6754,-1.6447,"CE2",1.0)});
    sideChains["TRP"].insert({"CE3",Coords(0.3664,-2.1984,-0.4947,"CE3",1.0)});
    sideChains["TRP"].insert({"CZ2",Coords(-1.4976,-1.6584,-2.5897,"CZ2",1.0)});
    sideChains["TRP"].insert({"CZ3",Coords(0.0144,-3.1794,-1.4507,"CZ3",1.0)});
    sideChains["TRP"].insert({"CH2",Coords(-0.8776,-2.8934,-2.4787,"CH2",1.0)});
    sideChains["TYR"].insert({"N",Coords(-0.1382,-0.1286,2.4852,"N",1.0)});
    sideChains["TYR"].insert({"CA",Coords(0.8638,-0.9836,1.8112,"CA",1.0)});
    sideChains["TYR"].insert({"C",Coords(2.0308,-1.2416,2.7522,"C",1.0)});
    sideChains["TYR"].insert({"O",Coords(2.4558,-0.3266,3.4722,"O",1.0)});
    sideChains["TYR"].insert({"CB",Coords(1.4218,-0.2836,0.5772,"CB",1.0)});
    sideChains["TYR"].insert({"CG",Coords(0.3608,0.0394,-0.4238,"CG",1.0)});
    sideChains["TYR"].insert({"CD1",Coords(-0.3762,1.2174,-0.3068,"CD1",1.0)});
    sideChains["TYR"].insert({"CD2",Coords(0.0498,-0.8396,-1.4578,"CD2",1.0)});
    sideChains["TYR"].insert({"CE1",Coords(-1.3732,1.5204,-1.1828,"CE1",1.0)});
    sideChains["TYR"].insert({"CE2",Coords(-0.9772,-0.5376,-2.3668,"CE2",1.0)});
    sideChains["TYR"].insert({"CZ",Coords(-1.6722,0.6424,-2.2108,"CZ",1.0)});
    sideChains["TYR"].insert({"OH",Coords(-2.6462,0.9214,-3.1488,"OH",1.0)});
    sideChains["VAL"].insert({"N",Coords(-0.0281,1.3981,-1.2411,"N",1.0)});
    sideChains["VAL"].insert({"CA",Coords(0.1279,0.7201,0.0459,"CA",1.0)});
    sideChains["VAL"].insert({"C",Coords(1.5869,0.3021,0.1039,"C",1.0)});
    sideChains["VAL"].insert({"O",Coords(2.0039,-0.6279,-0.6001,"O",1.0)});
    sideChains["VAL"].insert({"CB",Coords(-0.8001,-0.5079,0.1679,"CB",1.0)});
    sideChains["VAL"].insert({"CG1",Coords(-0.6111,-1.2009,1.5399,"CG1",1.0)});
    sideChains["VAL"].insert({"CG2",Coords(-2.2791,-0.0839,-0.0161,"CG2",1.0)});

    /*
     * add nucleics - center on C1'
     *
     */
// total models used :: 52
     centeringAtoms["DA"] = "C1'";
//    sideChains["DA"].insert({"P",Coords(-3.7834,-1.1321,3.1716,"P", 1.0)});
//    sideChains["DA"].insert({"OP1",Coords(-3.9316,-2.3331,4.0326,"OP1", 1.0)});
//    sideChains["DA"].insert({"OP2",Coords(-3.9939,0.2170,3.7521,"OP2", 1.0)});
//    sideChains["DA"].insert({"O5'",Coords(-2.3354,-1.1802,2.5093,"O5'", 1.0)});
//    sideChains["DA"].insert({"C5'",Coords(-1.9156,-2.3251,1.7715,"C5'", 1.0)});
//    sideChains["DA"].insert({"C4'",Coords(-0.6777,-2.0131,0.9647,"C4'", 1.0)});
//    sideChains["DA"].insert({"O4'",Coords(-0.9781,-1.0233,-0.0494,"O4'", 1.0)});
//    sideChains["DA"].insert({"C3'",Coords(0.5373,-1.4670,1.7376,"C3'", 1.0)});
//    sideChains["DA"].insert({"O3'",Coords(1.6305,-2.0623,1.0310,"O3'", 1.0)});
//    sideChains["DA"].insert({"C2'",Coords(0.4758,0.0194,1.4374,"C2'", 1.0)});
//    sideChains["DA"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DA"].insert({"N9",Coords(-0.6163,1.2378,-0.4736,"N9", 1.0)});
    sideChains["DA"].insert({"C8",Coords(-1.4800,2.0752,0.1883,"C8", 1.0)});
    sideChains["DA"].insert({"N7",Coords(-1.8648,3.1067,-0.5244,"N7", 1.0)});
    sideChains["DA"].insert({"C5",Coords(-1.2085,2.9383,-1.7363,"C5", 1.0)});
    sideChains["DA"].insert({"C6",Coords(-1.1975,3.6935,-2.9224,"C6", 1.0)});
    sideChains["DA"].insert({"N6",Coords(-1.8961,4.8179,-3.0914,"N6", 1.0)});
    sideChains["DA"].insert({"N1",Coords(-0.4356,3.2469,-3.9432,"N1", 1.0)});
    sideChains["DA"].insert({"C2",Coords(0.2629,2.1172,-3.7755,"C2", 1.0)});
    sideChains["DA"].insert({"N3",Coords(0.3336,1.3211,-2.7130,"N3", 1.0)});
    sideChains["DA"].insert({"C4",Coords(-0.4344,1.7928,-1.7167,"C4", 1.0)});

    // total models used :: 98
    centeringAtoms["DC"] = "C1'";
//    sideChains["DC"].insert({"P",Coords(2.6615,-0.9099,2.7762,"P", 1.0)});
//    sideChains["DC"].insert({"OP1",Coords(2.7370,0.2994,3.5904,"OP1", 1.0)});
//    sideChains["DC"].insert({"OP2",Coords(2.0248,-2.1033,3.3144,"OP2", 1.0)});
//    sideChains["DC"].insert({"O5'",Coords(1.9100,-0.5379,1.3914,"O5'", 1.0)});
//    sideChains["DC"].insert({"C5'",Coords(0.4748,-0.3031,1.3932,"C5'", 1.0)});
//    sideChains["DC"].insert({"C4'",Coords(0.0000,0.0000,0.0000,"C4'", 1.0)});
//    sideChains["DC"].insert({"O4'",Coords(0.4376,1.3393,-0.2899,"O4'", 1.0)});
//    sideChains["DC"].insert({"C3'",Coords(0.5804,-0.8694,-1.1334,"C3'", 1.0)});
//    sideChains["DC"].insert({"O3'",Coords(-0.4206,-1.0303,-2.1199,"O3'", 1.0)});
//    sideChains["DC"].insert({"C2'",Coords(1.7182,-0.0076,-1.6954,"C2'", 1.0)});
//    sideChains["DC"].insert({"C1'",Coords(1.1034,1.3867,-1.5341,"C1'", 1.0)});
    sideChains["DC"].insert({"N1",Coords(2.1457,2.4507,-1.4607,"N1", 1.0)});
    sideChains["DC"].insert({"C2",Coords(2.5596,3.0597,-2.6338,"C2", 1.0)});
    sideChains["DC"].insert({"O2",Coords(1.9967,2.7290,-3.6808,"O2", 1.0)});
    sideChains["DC"].insert({"N3",Coords(3.5130,4.0095,-2.5649,"N3", 1.0)});
    sideChains["DC"].insert({"C4",Coords(4.1059,4.3450,-1.4145,"C4", 1.0)});
    sideChains["DC"].insert({"N4",Coords(5.0814,5.2734,-1.4526,"N4", 1.0)});
    sideChains["DC"].insert({"C5",Coords(3.7441,3.6935,-0.2047,"C5", 1.0)});
    sideChains["DC"].insert({"C6",Coords(2.7570,2.7863,-0.2907,"C6", 1.0)});

    // total models used :: 123
    centeringAtoms["DG"] = "C1'";
//    sideChains["DG"].insert({"P",Coords(-5.3892,1.4224,-0.6462,"P", 1.0)});
//    sideChains["DG"].insert({"OP1",Coords(-5.9068,0.8283,-1.9037,"OP1", 1.0)});
//    sideChains["DG"].insert({"OP2",Coords(-6.3033,1.6057,0.5091,"OP2", 1.0)});
//    sideChains["DG"].insert({"O5'",Coords(-4.1469,0.5485,-0.1688,"O5'", 1.0)});
//    sideChains["DG"].insert({"C5'",Coords(-3.3172,0.9848,0.9073,"C5'", 1.0)});
//    sideChains["DG"].insert({"C4'",Coords(-2.1103,0.0853,1.0314,"C4'", 1.0)});
//    sideChains["DG"].insert({"O4'",Coords(-1.3986,0.0552,-0.2214,"O4'", 1.0)});
//    sideChains["DG"].insert({"C3'",Coords(-1.1073,0.5050,2.0955,"C3'", 1.0)});
//    sideChains["DG"].insert({"O3'",Coords(-1.4250,-0.2662,3.2617,"O3'", 1.0)});
//    sideChains["DG"].insert({"C2'",Coords(0.2364,0.1060,1.5058,"C2'", 1.0)});
//    sideChains["DG"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DG"].insert({"N9",Coords(0.6112,1.0526,-0.8052,"N9", 1.0)});
    sideChains["DG"].insert({"C8",Coords(1.4635,0.8898,-1.8711,"C8", 1.0)});
    sideChains["DG"].insert({"N7",Coords(1.8294,2.0198,-2.4119,"N7", 1.0)});
    sideChains["DG"].insert({"C5",Coords(1.1835,2.9882,-1.6550,"C5", 1.0)});
    sideChains["DG"].insert({"C6",Coords(1.1912,4.4025,-1.7696,"C6", 1.0)});
    sideChains["DG"].insert({"O6",Coords(1.7907,5.1060,-2.5926,"O6", 1.0)});
    sideChains["DG"].insert({"N1",Coords(0.3981,4.9985,-0.7949,"N1", 1.0)});
    sideChains["DG"].insert({"C2",Coords(-0.3132,4.3239,0.1652,"C2", 1.0)});
    sideChains["DG"].insert({"N2",Coords(-1.0204,5.0800,1.0187,"N2", 1.0)});
    sideChains["DG"].insert({"N3",Coords(-0.3324,3.0070,0.2815,"N3", 1.0)});
    sideChains["DG"].insert({"C4",Coords(0.4322,2.4070,-0.6551,"C4", 1.0)});

    // total models used :: 48
    centeringAtoms["DT"] = "C1'";
//    sideChains["DT"].insert({"P",Coords(0.2077,-4.8901,1.8154,"P", 1.0)});
//    sideChains["DT"].insert({"OP1",Coords(0.5055,-4.2532,3.1236,"OP1", 1.0)});
//    sideChains["DT"].insert({"OP2",Coords(-0.6843,-6.0759,1.7626,"OP2", 1.0)});
//    sideChains["DT"].insert({"O5'",Coords(-0.3787,-3.7657,0.8507,"O5'", 1.0)});
//    sideChains["DT"].insert({"C5'",Coords(-1.5505,-3.0343,1.2124,"C5'", 1.0)});
//    sideChains["DT"].insert({"C4'",Coords(-1.5498,-1.6778,0.5460,"C4'", 1.0)});
//    sideChains["DT"].insert({"O4'",Coords(-0.4193,-0.9102,1.0101,"O4'", 1.0)});
//    sideChains["DT"].insert({"C3'",Coords(-1.4020,-1.6906,-0.9813,"C3'", 1.0)});
//    sideChains["DT"].insert({"O3'",Coords(-2.6956,-1.7804,-1.5878,"O3'", 1.0)});
//    sideChains["DT"].insert({"C2'",Coords(-0.7797,-0.3361,-1.2708,"C2'", 1.0)});
//    sideChains["DT"].insert({"C1'",Coords(0.0000,0.0000,0.0000,"C1'", 1.0)});
    sideChains["DT"].insert({"N1",Coords(1.4624,-0.1211,-0.1367,"N1", 1.0)});
    sideChains["DT"].insert({"C2",Coords(2.1765,0.9990,-0.4979,"C2", 1.0)});
    sideChains["DT"].insert({"O2",Coords(1.6572,2.0825,-0.7117,"O2", 1.0)});
    sideChains["DT"].insert({"N3",Coords(3.5320,0.8031,-0.5983,"N3", 1.0)});
    sideChains["DT"].insert({"C4",Coords(4.2244,-0.3724,-0.3805,"C4", 1.0)});
    sideChains["DT"].insert({"O4",Coords(5.4457,-0.3960,-0.5079,"O4", 1.0)});
    sideChains["DT"].insert({"C5",Coords(3.4109,-1.5070,-0.0078,"C5", 1.0)});
    sideChains["DT"].insert({"C7",Coords(4.0731,-2.8235,0.2471,"C7", 1.0)});
    sideChains["DT"].insert({"C6",Coords(2.0885,-1.3261,0.0936,"C6", 1.0)});

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
    sideChains["CYT"].insert({"N1",Coords(-0.5794,-0.0537,1.3552,"N1", 1.0)});
    sideChains["CYT"].insert({"C2",Coords(-0.4380,1.1006,2.1377,"C2", 1.0)});
    sideChains["CYT"].insert({"O2",Coords(0.1336,2.0803,1.6314,"O2", 1.0)});
    sideChains["CYT"].insert({"N3",Coords(-0.9473,1.0821,3.3972,"N3", 1.0)});
    sideChains["CYT"].insert({"C4",Coords(-1.5532,-0.0221,3.8647,"C4", 1.0)});
    sideChains["CYT"].insert({"N4",Coords(-2.0557,0.0193,5.1057,"N4", 1.0)});
    sideChains["CYT"].insert({"C5",Coords(-1.7043,-1.2074,3.0753,"C5", 1.0)});
    sideChains["CYT"].insert({"C6",Coords(-1.1942,-1.1866,1.8370,"C6", 1.0)});

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
    sideChains["DNA"].insert({"P",Coords(3.6962,-1.2041,-0.2176,"P", 1.0)});
//    sideChains["DNA"].insert({"OP1",Coords(3.8853,-2.1457,0.9082,"OP1", 1.0)});
//    sideChains["DNA"].insert({"OP2",Coords(4.3353,-1.4720,-1.5237,"OP2", 1.0)});
    sideChains["DNA"].insert({"O5'",Coords(2.1332,-1.0217,-0.4514,"O5'", 1.0)});
    sideChains["DNA"].insert({"C5'",Coords(1.2963,-0.5015,0.5897,"C5'", 1.0)});
    sideChains["DNA"].insert({"C4'",Coords(0.0000,0.0000,0.0000,"C4'", 1.0)});
    sideChains["DNA"].insert({"O4'",Coords(0.2309,1.2644,-0.6669,"O4'", 1.0)});
    sideChains["DNA"].insert({"C3'",Coords(-0.6141,-0.9136,-1.0622,"C3'", 1.0)});
    sideChains["DNA"].insert({"O3'",Coords(-2.0251,-0.7011,-1.0625,"O3'", 1.0)});
    sideChains["DNA"].insert({"C2'",Coords(-0.1427,-0.2829,-2.3562,"C2'", 1.0)});
    sideChains["DNA"].insert({"C1'",Coords(-0.3023,1.1710,-1.9746,"C1'", 1.0)});
//    sideChains["DG"].insert({"N9",Coords(0.3796,2.1420,-2.8215,"N9", 1.0)});
//    sideChains["DG"].insert({"C8",Coords(1.4450,1.9262,-3.6606,"C8", 1.0)});
//    sideChains["DG"].insert({"N7",Coords(1.8109,2.9989,-4.3081,"N7", 1.0)});
//    sideChains["DG"].insert({"C5",Coords(0.9217,3.9760,-3.8836,"C5", 1.0)});
//    sideChains["DG"].insert({"C6",Coords(0.8230,5.3454,-4.2393,"C6", 1.0)});
//    sideChains["DG"].insert({"O6",Coords(1.5182,5.9879,-5.0317,"O6", 1.0)});
//    sideChains["DG"].insert({"N1",Coords(-0.2221,5.9725,-3.5682,"N1", 1.0)});
//    sideChains["DG"].insert({"C2",Coords(-1.0633,5.3627,-2.6718,"C2", 1.0)});
//    sideChains["DG"].insert({"N2",Coords(-2.0145,6.1374,-2.1320,"N2", 1.0)});
//    sideChains["DG"].insert({"N3",Coords(-0.9770,4.0900,-2.3270,"N3", 1.0)});
//    sideChains["DG"].insert({"C4",Coords(0.0296,3.4616,-2.9684,"C4", 1.0)});

    std::vector<Coords> * pwater;
    // Create Vector for new Residue ALA
    waters.insert(std::pair<std::string, std::vector<Coords> >("ALA", std::vector<Coords>(12)));
    pwater = &waters["ALA"];
    (*pwater)[0].x = 1.8194;
    (*pwater)[0].y = -3.0286f;
    (*pwater)[0].z = 0.1559;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = 4.7257;
    (*pwater)[1].y = 0.1993;
    (*pwater)[1].z = -1.9795f;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -0.0286f;
    (*pwater)[2].y = -0.2394f;
    (*pwater)[2].z = 4.6083;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -0.2119f;
    (*pwater)[3].y = -1.7953f;
    (*pwater)[3].z = -3.8500f;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 3.1218;
    (*pwater)[4].y = 0.8221;
    (*pwater)[4].z = 2.2970;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 1.0030;
    (*pwater)[5].y = -3.1859f;
    (*pwater)[5].z = -1.9756f;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.6670;
    (*pwater)[6].x = 4.4477;
    (*pwater)[6].y = 2.5877;
    (*pwater)[6].z = -0.6527f;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -2.0723f;
    (*pwater)[7].y = -0.4076f;
    (*pwater)[7].z = 3.6552;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 4.2466;
    (*pwater)[8].y = -1.4169f;
    (*pwater)[8].z = 0.0168;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 2.0576;
    (*pwater)[9].y = 3.0443;
    (*pwater)[9].z = -3.5803f;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = 3.3228;
    (*pwater)[10].y = -1.6780f;
    (*pwater)[10].z = -3.2131f;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 2.3350;
    (*pwater)[11].y = -1.6908f;
    (*pwater)[11].z = 2.0430;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;

    // Create Vector for new Residue GLY
    waters.insert(std::pair<std::string, std::vector<Coords> >("GLY", std::vector<Coords>(10)));
    pwater = &waters["GLY"];
    (*pwater)[0].x = 1.6164;
    (*pwater)[0].y = -2.7153;
    (*pwater)[0].z = -1.1367;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 4.0055;
    (*pwater)[1].y = -0.3816;
    (*pwater)[1].z = 1.8228;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -2.1785;
    (*pwater)[2].y = -3.5583;
    (*pwater)[2].z = -0.4525;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 3.9719;
    (*pwater)[3].y = 0.7516;
    (*pwater)[3].z = -0.9554;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 1.8416;
    (*pwater)[4].y = -1.0375;
    (*pwater)[4].z = 3.5277;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -0.4796;
    (*pwater)[5].y = -1.7318;
    (*pwater)[5].z = 3.0685;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 0.6059;
    (*pwater)[6].y = 1.8712;
    (*pwater)[6].z = 3.5282;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 1.4779;
    (*pwater)[7].y = 0.8419;
    (*pwater)[7].z = 3.6920;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = -3.1990;
    (*pwater)[8].y = -2.6064;
    (*pwater)[8].z = 1.8467;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = -0.2441;
    (*pwater)[9].y = 0.2452;
    (*pwater)[9].z = -4.5276;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;


// Create Vector for new Residue ARG
    waters.insert(std::pair<std::string, std::vector<Coords> >("ARG", std::vector<Coords>(21) ));
    pwater = &waters["ARG"];

    (*pwater)[0].x = -0.9202;
    (*pwater)[0].y = -3.2030;
    (*pwater)[0].z = -5.1172;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = 3.7116;
    (*pwater)[1].y = 1.8026;
    (*pwater)[1].z = -2.1487;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = -1.4915;
    (*pwater)[2].y = 0.7648;
    (*pwater)[2].z = -6.0131;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -0.0641;
    (*pwater)[3].y = -3.5482;
    (*pwater)[3].z = -1.6271;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 4.1228;
    (*pwater)[4].y = 0.5081;
    (*pwater)[4].z = 1.7213;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 5.0371;
    (*pwater)[5].y = 0.3402;
    (*pwater)[5].z = -0.0500;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = -2.0521;
    (*pwater)[6].y = 2.5556;
    (*pwater)[6].z = -2.4945;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 3.6150;
    (*pwater)[7].y = -2.6461;
    (*pwater)[7].z = -0.9889;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 0.8890;
    (*pwater)[8].y = -0.7937;
    (*pwater)[8].z = -6.9492;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = -0.6157;
    (*pwater)[9].y = 4.1316;
    (*pwater)[9].z = 1.0697;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = -2.8522;
    (*pwater)[10].y = 0.2233;
    (*pwater)[10].z = -1.1125;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 3.0449;
    (*pwater)[11].y = -3.1499;
    (*pwater)[11].z = 0.7095;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = 0.5620;
    (*pwater)[12].y = 3.4120;
    (*pwater)[12].z = -1.4803;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 0.3330;
    (*pwater)[13].x = 2.7574;
    (*pwater)[13].y = 2.0884;
    (*pwater)[13].z = 3.2493;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.3330;
    (*pwater)[14].x = 0.1818;
    (*pwater)[14].y = -2.3448;
    (*pwater)[14].z = -7.3175;
    (*pwater)[14].type = "O";
    (*pwater)[14].occ = 0.3330;
    (*pwater)[15].x = 3.8552;
    (*pwater)[15].y = 1.1529;
    (*pwater)[15].z = -5.4908;
    (*pwater)[15].type = "O";
    (*pwater)[15].occ = 0.3330;
    (*pwater)[16].x = -2.6921;
    (*pwater)[16].y = -2.4749;
    (*pwater)[16].z = -0.1004;
    (*pwater)[16].type = "O";
    (*pwater)[16].occ = 0.3330;
    (*pwater)[17].x = 0.3297;
    (*pwater)[17].y = -3.4427;
    (*pwater)[17].z = 0.6405;
    (*pwater)[17].type = "O";
    (*pwater)[17].occ = 0.3330;
    (*pwater)[18].x = 2.3617;
    (*pwater)[18].y = 2.9007;
    (*pwater)[18].z = -4.6133;
    (*pwater)[18].type = "O";
    (*pwater)[18].occ = 0.3330;
    (*pwater)[19].x = 1.5189;
    (*pwater)[19].y = -3.5564;
    (*pwater)[19].z = -2.4670;
    (*pwater)[19].type = "O";
    (*pwater)[19].occ = 0.3330;
    (*pwater)[20].x = 0.3913;
    (*pwater)[20].y = 2.4426;
    (*pwater)[20].z = -5.9981;
    (*pwater)[20].type = "O";
    (*pwater)[20].occ = 0.3330;

    // ASN
// Create Vector for new Residue ASN
    waters.insert(std::pair<std::string, std::vector<Coords> >("ASN", std::vector<Coords>(9)));
    pwater = &waters["ASN"];

    (*pwater)[0].x = -2.8171;
    (*pwater)[0].y = -0.7220;
    (*pwater)[0].z = 0.0598;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 0.6216;
    (*pwater)[1].y = -4.3980;
    (*pwater)[1].z = 3.5030;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = 1.6212;
    (*pwater)[2].y = -4.3882;
    (*pwater)[2].z = -0.7999;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 1.5040;
    (*pwater)[3].y = -1.2586;
    (*pwater)[3].z = 3.8511;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.6670;
    (*pwater)[4].x = -1.3369;
    (*pwater)[4].y = -3.2494;
    (*pwater)[4].z = -3.1901;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -4.6528;
    (*pwater)[5].y = 2.3431;
    (*pwater)[5].z = 0.8011;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 0.1112;
    (*pwater)[6].y = 1.1631;
    (*pwater)[6].z = 4.1515;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 3.6856;
    (*pwater)[7].y = -0.8639;
    (*pwater)[7].z = -1.8972;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = -1.5136;
    (*pwater)[8].y = -3.8739;
    (*pwater)[8].z = 2.8804;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;

    // ASP
    // Create Vector for new Residue ASP
    waters["ASP"] = std::vector<Coords>(12);
    pwater = &waters["ASP"];

    (*pwater)[0].x = 1.6853;
    (*pwater)[0].y = 2.9411;
    (*pwater)[0].z = 3.6584;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = -2.5141;
    (*pwater)[1].y = 2.5170;
    (*pwater)[1].z = 2.7862;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -1.2643;
    (*pwater)[2].y = 0.3306;
    (*pwater)[2].z = 3.9048;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 3.4443;
    (*pwater)[3].y = -3.0228;
    (*pwater)[3].z = 1.3967;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 3.5437;
    (*pwater)[4].y = -0.2716;
    (*pwater)[4].z = -1.6293;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 3.8170;
    (*pwater)[5].y = -0.8801;
    (*pwater)[5].z = -0.1211;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 0.5566;
    (*pwater)[6].y = 4.6603;
    (*pwater)[6].z = 0.1878;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 1.9622;
    (*pwater)[7].y = -1.1092;
    (*pwater)[7].z = 3.0756;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 1.6527;
    (*pwater)[8].y = -2.9711;
    (*pwater)[8].z = 2.7722;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 2.6802;
    (*pwater)[9].y = 2.1310;
    (*pwater)[9].z = -1.2750;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = 0.3306;
    (*pwater)[10].y = 4.2698;
    (*pwater)[10].z = 2.7611;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 2.4564;
    (*pwater)[11].y = 3.5514;
    (*pwater)[11].z = 0.7012;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;

    // CYS
    // Create Vector for new Residue CYS
    waters["CYS"] = std::vector<Coords>(5);
    pwater = &waters["CYS"];
    (*pwater)[0].x = -4.1878f;
    (*pwater)[0].y = -1.6206f;
    (*pwater)[0].z = 1.0453;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = -4.8562f;
    (*pwater)[1].y = 2.4769;
    (*pwater)[1].z = -0.4143f;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -0.6391f;
    (*pwater)[2].y = 3.9497;
    (*pwater)[2].z = -0.2930f;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -3.2043f;
    (*pwater)[3].y = 0.2392;
    (*pwater)[3].z = -2.2271f;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 0.0960;
    (*pwater)[4].y = 3.1794;
    (*pwater)[4].z = 2.2657;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;


    // GLN
    // Create Vector for new Residue GLN
    waters["GLN"] = std::vector<Coords>(7);
    pwater = &waters["GLN"];
    (*pwater)[0].x = 3.6571;
    (*pwater)[0].y = 1.2530;
    (*pwater)[0].z = -3.7994;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -1.3633;
    (*pwater)[1].y = 2.6682;
    (*pwater)[1].z = -1.7547;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 4.8520;
    (*pwater)[2].y = -1.4280;
    (*pwater)[2].z = -2.1606;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 0.0716;
    (*pwater)[3].y = -1.6981;
    (*pwater)[3].z = -3.5993;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 3.7159;
    (*pwater)[4].y = 3.5766;
    (*pwater)[4].z = -1.3324;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -3.4612;
    (*pwater)[5].y = -2.0754;
    (*pwater)[5].z = -2.2215;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 5.4311;
    (*pwater)[6].y = 2.7412;
    (*pwater)[6].z = -0.9313;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;

    // GLU
// Create Vector for new Residue GLU
    waters["GLU"] = std::vector<Coords>(14);
    pwater = &waters["GLU"];
    (*pwater)[0].x = -2.9360;
    (*pwater)[0].y = -1.1168;
    (*pwater)[0].z = 0.6716;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -3.0990;
    (*pwater)[1].y = 1.2488;
    (*pwater)[1].z = -2.2489;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = -2.0902;
    (*pwater)[2].y = 2.7773;
    (*pwater)[2].z = 3.7612;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 0.1893;
    (*pwater)[3].y = 3.3134;
    (*pwater)[3].z = 4.3816;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -0.2553;
    (*pwater)[4].y = -3.9246;
    (*pwater)[4].z = 1.1594;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -2.8243;
    (*pwater)[5].y = -2.6349;
    (*pwater)[5].z = -4.2005;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 1.2205;
    (*pwater)[6].y = 3.6917;
    (*pwater)[6].z = 0.1794;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -0.7029;
    (*pwater)[7].y = 3.9004;
    (*pwater)[7].z = -0.9380;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 0.3579;
    (*pwater)[8].y = -1.1831;
    (*pwater)[8].z = 4.5715;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 2.7686;
    (*pwater)[9].y = 2.9763;
    (*pwater)[9].z = 1.3065;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = -1.6243;
    (*pwater)[10].y = -0.0216;
    (*pwater)[10].z = 4.8823;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 0.8202;
    (*pwater)[11].y = 0.9833;
    (*pwater)[11].z = 4.9212;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = -3.0688;
    (*pwater)[12].y = 3.9501;
    (*pwater)[12].z = 1.2100;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 0.3330;
    (*pwater)[13].x = -3.4749;
    (*pwater)[13].y = -0.2795;
    (*pwater)[13].z = -3.9040;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.3330;

    // HIS
    // Create Vector for new Residue HIS
    waters["HIS"] = std::vector<Coords>(7);
    pwater = &waters["HIS"];
    (*pwater)[0].x = -4.3518;
    (*pwater)[0].y = -3.1090;
    (*pwater)[0].z = -1.2904;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = 3.5241;
    (*pwater)[1].y = -0.4878;
    (*pwater)[1].z = 2.9943;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = 0.9139;
    (*pwater)[2].y = 3.2521;
    (*pwater)[2].z = 2.9615;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -1.4703;
    (*pwater)[3].y = 2.9121;
    (*pwater)[3].z = 1.2264;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -3.0409;
    (*pwater)[4].y = -4.1669;
    (*pwater)[4].z = 1.0670;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -5.9028;
    (*pwater)[5].y = -0.7166;
    (*pwater)[5].z = -0.4736;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 3.2458;
    (*pwater)[6].y = -2.3028;
    (*pwater)[6].z = -0.3388;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;

    // ILE
    // Create Vector for new Residue ILE
    waters["ILE"] = std::vector<Coords>(13);
    pwater = &waters["ILE"];
    (*pwater)[0].x = 3.5850;
    (*pwater)[0].y = -4.5614;
    (*pwater)[0].z = -1.4590;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = 3.7145;
    (*pwater)[1].y = -0.4266;
    (*pwater)[1].z = 1.8319;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = 1.4073;
    (*pwater)[2].y = -2.9962;
    (*pwater)[2].z = 4.0535;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -4.9984;
    (*pwater)[3].y = -0.5650;
    (*pwater)[3].z = 3.5030;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 1.4572;
    (*pwater)[4].y = -3.9516;
    (*pwater)[4].z = 1.7508;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -4.7012;
    (*pwater)[5].y = 1.5256;
    (*pwater)[5].z = 4.5302;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 2.9520;
    (*pwater)[6].y = 1.8150;
    (*pwater)[6].z = 2.3363;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -0.7989;
    (*pwater)[7].y = -2.2725;
    (*pwater)[7].z = 4.5129;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 0.0093;
    (*pwater)[8].y = -3.9403;
    (*pwater)[8].z = -0.0011;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = -0.4385;
    (*pwater)[9].y = 5.2123;
    (*pwater)[9].z = 1.0803;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = 2.0102;
    (*pwater)[10].y = -3.2656;
    (*pwater)[10].z = -3.5013;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = -3.4773;
    (*pwater)[11].y = -2.9465;
    (*pwater)[11].z = 2.5979;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = -0.0907;
    (*pwater)[12].y = -0.8461;
    (*pwater)[12].z = 5.4360;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 0.3330;

    // LEU
    // Create Vector for new Residue LEU
    waters["LEU"] = std::vector<Coords>(10);
    pwater = &waters["LEU"];
    (*pwater)[0].x = 2.6135;
    (*pwater)[0].y = -0.8435;
    (*pwater)[0].z = 3.0639;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = 2.9171;
    (*pwater)[1].y = -1.9274;
    (*pwater)[1].z = -1.3897;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 1.4662;
    (*pwater)[2].y = 1.4511;
    (*pwater)[2].z = -5.4427;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 4.4627;
    (*pwater)[3].y = 3.4386;
    (*pwater)[3].z = -1.3742;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -1.2462;
    (*pwater)[4].y = -0.0756;
    (*pwater)[4].z = 4.0308;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 3.9783;
    (*pwater)[5].y = -0.3993;
    (*pwater)[5].z = 1.3453;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 3.8264;
    (*pwater)[6].y = 1.3502;
    (*pwater)[6].z = 2.6594;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 2.1936;
    (*pwater)[7].y = -4.9280;
    (*pwater)[7].z = 1.4838;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 0.9350;
    (*pwater)[8].y = 5.6672;
    (*pwater)[8].z = 2.2716;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 3.4443;
    (*pwater)[9].y = 4.7909;
    (*pwater)[9].z = 1.7566;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;

    // LYS
// Create Vector for new Residue LYS
    waters["LYS"] = std::vector<Coords>(11);
    pwater = &waters["LYS"];
    (*pwater)[0].x = 2.9215;
    (*pwater)[0].y = -3.4968;
    (*pwater)[0].z = -3.0286;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = 2.6353;
    (*pwater)[1].y = -1.2531;
    (*pwater)[1].z = 1.7792;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -0.2691;
    (*pwater)[2].y = -4.1620;
    (*pwater)[2].z = -3.6566;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 1.8398;
    (*pwater)[3].y = 0.0435;
    (*pwater)[3].z = -6.0966;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -0.6216;
    (*pwater)[4].y = 1.1652;
    (*pwater)[4].z = -5.1751;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 4.4267;
    (*pwater)[5].y = 0.1670;
    (*pwater)[5].z = -2.1073;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 0.8504;
    (*pwater)[6].y = -3.9537;
    (*pwater)[6].z = -0.4357;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -0.3562;
    (*pwater)[7].y = -2.1610;
    (*pwater)[7].z = -4.9159;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 3.0896;
    (*pwater)[8].y = 1.2003;
    (*pwater)[8].z = -3.8057;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = -0.8425;
    (*pwater)[9].y = -1.9435;
    (*pwater)[9].z = -6.1917;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = 2.3883;
    (*pwater)[10].y = -1.6387;
    (*pwater)[10].z = -6.2405;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;

    // MET
    // Create Vector for new Residue MET
    waters["MET"] = std::vector<Coords>(2);
    pwater = &waters["MET"];
    (*pwater)[0].x = 1.2020;
    (*pwater)[0].y = 3.1766;
    (*pwater)[0].z = -3.0646;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = -1.7255;
    (*pwater)[1].y = -2.7419;
    (*pwater)[1].z = 2.0694;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;

    // PHE
// Create Vector for new Residue PHE
    waters["PHE"] = std::vector<Coords>(9);
    pwater = &waters["PHE"];
    (*pwater)[0].x = 5.1444;
    (*pwater)[0].y = -0.8998;
    (*pwater)[0].z = -0.3772;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = -2.5157;
    (*pwater)[1].y = 0.8654;
    (*pwater)[1].z = 2.9472;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -3.6910;
    (*pwater)[2].y = -3.8688;
    (*pwater)[2].z = 1.3881;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -4.0298;
    (*pwater)[3].y = -1.5403;
    (*pwater)[3].z = 2.3317;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = 0.7490;
    (*pwater)[4].y = -6.5949;
    (*pwater)[4].z = -1.4971;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 1.8581;
    (*pwater)[5].y = 0.3850;
    (*pwater)[5].z = 3.9751;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 0.6362;
    (*pwater)[6].y = -3.7793;
    (*pwater)[6].z = 3.2490;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 3.3106;
    (*pwater)[7].y = -4.1907;
    (*pwater)[7].z = 2.2425;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = -0.5794;
    (*pwater)[8].y = 6.9281;
    (*pwater)[8].z = 0.0063;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;

    // PRO
    // Create Vector for new Residue PRO
    waters["PRO"] = std::vector<Coords>(9);
    pwater = &waters["PRO"];
    (*pwater)[0].x = 2.3440;
    (*pwater)[0].y = -0.3463;
    (*pwater)[0].z = 3.3234;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 0.4841;
    (*pwater)[1].y = -4.4487;
    (*pwater)[1].z = 1.6156;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = 0.8795;
    (*pwater)[2].y = -2.4198;
    (*pwater)[2].z = 4.0414;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 0.2620;
    (*pwater)[3].y = -2.2055;
    (*pwater)[3].z = -4.0237;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -0.5291;
    (*pwater)[4].y = -3.3470;
    (*pwater)[4].z = 3.5450;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = -3.4629;
    (*pwater)[5].y = -3.1900;
    (*pwater)[5].z = -1.8817;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 4.1331;
    (*pwater)[6].y = -1.9238;
    (*pwater)[6].z = 0.5205;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 4.5634;
    (*pwater)[7].y = 2.5348;
    (*pwater)[7].z = 1.6607;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 1.5337;
    (*pwater)[8].y = -4.2773;
    (*pwater)[8].z = -1.1376;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;

    // SER
    // Create Vector for new Residue SER
    waters["SER"] = std::vector<Coords>(7);
    pwater = &waters["SER"];
    (*pwater)[0].x = 1.8881;
    (*pwater)[0].y = 0.5784;
    (*pwater)[0].z = 3.3323;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -2.8149;
    (*pwater)[1].y = 2.5513;
    (*pwater)[1].z = 2.5862;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = -2.7848;
    (*pwater)[2].y = -0.9318;
    (*pwater)[2].z = 2.2938;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.6670;
    (*pwater)[3].x = 0.8809;
    (*pwater)[3].y = -5.6659;
    (*pwater)[3].z = -0.3631;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -3.5865;
    (*pwater)[4].y = 2.8896;
    (*pwater)[4].z = -0.3212;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 1.3636;
    (*pwater)[5].y = 3.7513;
    (*pwater)[5].z = 2.6829;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 2.8889;
    (*pwater)[6].y = -3.5453;
    (*pwater)[6].z = 1.6646;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;

    // THR
    // Create Vector for new Residue THR
    waters["THR"] = std::vector<Coords>(10);
    pwater = &waters["THR"];
    (*pwater)[0].x = -0.1170;
    (*pwater)[0].y = 2.7503;
    (*pwater)[0].z = 2.6289;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -2.4063;
    (*pwater)[1].y = -2.9645;
    (*pwater)[1].z = 2.1163;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = 0.0747;
    (*pwater)[2].y = -0.7845;
    (*pwater)[2].z = 4.0407;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.6670;
    (*pwater)[3].x = 1.5311;
    (*pwater)[3].y = -3.1712;
    (*pwater)[3].z = 1.2931;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -3.8421;
    (*pwater)[4].y = 0.5495;
    (*pwater)[4].z = 3.2340;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;
    (*pwater)[5].x = 3.5293;
    (*pwater)[5].y = 0.4756;
    (*pwater)[5].z = 3.4535;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = -5.9016;
    (*pwater)[6].y = -0.6513;
    (*pwater)[6].z = 0.6366;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -0.3062;
    (*pwater)[7].y = -3.0772;
    (*pwater)[7].z = 2.4631;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 2.6645;
    (*pwater)[8].y = -1.7422;
    (*pwater)[8].z = 3.5524;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 2.6795;
    (*pwater)[9].y = 1.3298;
    (*pwater)[9].z = 2.3779;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;

    //TRP
    // Create Vector for new Residue TRP
    waters["TRP"] = std::vector<Coords>(4);
    pwater = &waters["TRP"];
    (*pwater)[0].x = -3.6537;
    (*pwater)[0].y = 1.5699;
    (*pwater)[0].z = -3.3903;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.3330;
    (*pwater)[1].x = -1.2314;
    (*pwater)[1].y = -3.6476;
    (*pwater)[1].z = -5.4470;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -0.0637;
    (*pwater)[2].y = 5.3013;
    (*pwater)[2].z = 1.0488;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -3.6294;
    (*pwater)[3].y = -5.0915;
    (*pwater)[3].z = -2.2316;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;

    // TYR
    // Create Vector for new Residue TYR
    waters["TYR"] = std::vector<Coords>(5);
    pwater = &waters["TYR"];
    (*pwater)[0].x = -3.7604;
    (*pwater)[0].y = 3.2507;
    (*pwater)[0].z = -2.9513;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = 2.0821;
    (*pwater)[1].y = -3.8333;
    (*pwater)[1].z = 0.3806;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = -3.1787;
    (*pwater)[2].y = -0.7734;
    (*pwater)[2].z = -4.8465;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.6670;
    (*pwater)[3].x = 4.8459;
    (*pwater)[3].y = -1.6081;
    (*pwater)[3].z = -0.2639;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;
    (*pwater)[4].x = -1.5064;
    (*pwater)[4].y = -2.9448;
    (*pwater)[4].z = -4.8543;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 0.3330;

    // VAL
    // Create Vector for new Residue VAL
    waters["VAL"] = std::vector<Coords>(4);
    pwater = &waters["VAL"];
    (*pwater)[0].x = 0.8419;
    (*pwater)[0].y = 1.6832;
    (*pwater)[0].z = 3.4416;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = -0.1380;
    (*pwater)[1].y = -4.7028;
    (*pwater)[1].z = 1.9950;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.3330;
    (*pwater)[2].x = 2.4260;
    (*pwater)[2].y = -3.3721;
    (*pwater)[2].z = -0.3480;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = -4.0055;
    (*pwater)[3].y = -1.2465;
    (*pwater)[3].z = 3.1386;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.3330;

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


    // Create Vector for new Residue DA
    waters["DA"] = std::vector<Coords>(3);
    pwater = &waters["DA"];
    (*pwater)[0].x = -3.5028;
    (*pwater)[0].y = 5.6436;
    (*pwater)[0].z = -0.4700;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 1.6761;
    (*pwater)[1].y = -1.0325;
    (*pwater)[1].z = -3.2875;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -2.0555;
    (*pwater)[2].y = 2.3420;
    (*pwater)[2].z = 3.4435;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;

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

    waters["DNA"] = std::vector<Coords>(2);
    pwater = &waters["DNA"];
    (*pwater)[0].x = 2.9540;
    (*pwater)[0].y = -2.0312;
    (*pwater)[0].z = -3.0254;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -2.9630;
    (*pwater)[1].y = 4.4956;
    (*pwater)[1].z = -0.4762;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;


    // Create Vector for new Residue DC
    waters["DC"] = std::vector<Coords>(4);
    pwater = &waters["DC"];
    (*pwater)[0].x = 5.5571;
    (*pwater)[0].y = 2.4535;
    (*pwater)[0].z = -0.5641;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -2.1049;
    (*pwater)[1].y = 3.7118;
    (*pwater)[1].z = -5.3136;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -1.7182;
    (*pwater)[2].y = 2.9993;
    (*pwater)[2].z = 1.6921;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 1.6478;
    (*pwater)[3].y = 3.1839;
    (*pwater)[3].z = 4.2544;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;

    // Create Vector for new Residue DG
    waters["DG"] = std::vector<Coords>(6);
    pwater = &waters["DG"];
    (*pwater)[0].x = -1.2622;
    (*pwater)[0].y = 2.0057;
    (*pwater)[0].z = 2.8183;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -3.1615;
    (*pwater)[1].y = 4.1614;
    (*pwater)[1].z = 2.9041;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 1.3153;
    (*pwater)[2].y = 4.4565;
    (*pwater)[2].z = -5.3356;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 3.2426;
    (*pwater)[3].y = 1.9755;
    (*pwater)[3].z = -4.6978;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 3.4048;
    (*pwater)[4].y = 4.9277;
    (*pwater)[4].z = -4.6732;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 2.2562;
    (*pwater)[5].y = -2.1508;
    (*pwater)[5].z = -3.1460;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;


    // Create Vector for new Residue DT
    waters["DT"] = std::vector<Coords>(3);         // 3
    pwater = &waters["DT"];
    (*pwater)[0].x = 0.2881;
    (*pwater)[0].y = 4.4010;
    (*pwater)[0].z = -0.4560;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 7.5308;
    (*pwater)[1].y = -2.0235;
    (*pwater)[1].z = -0.1088;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 0.6567;
    (*pwater)[2].y = -4.3684;
    (*pwater)[2].z = -0.1900;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
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

    auto atomicNumbers = model.getAtomicNumberVec();

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
    float dwater = 2.8*2.7;

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
    }

    /*
     * prune waters too close to protein
     * does not consider hydrogens
     */
    for (unsigned int i=0; i<numAtoms; i++){

        currentID = pResIDs[i];
        currentRes = pResidues[i];
        boost::trim(currentRes);

        // could fail if residue not found - PDBModel should create the residue during instantiation
        float rlimit = 2.2;//1.4 + SASTOOLS_RESIDUES_H::Residues::getResidue(currentRes)->getAtom(atomType[i])->getRadii();

        vector3 anchor(xvalue[i], yvalue[i], zvalue[i]);
        rlimit *= rlimit; // compare the square distance to reduce operations in the loop

        for(int w=0; w < totalwaters; w++){ // sort the waters.  waters too close will be moved to end of the vector
            auto pW = &hydration[w]; // these waters belong to a resid
            float val = (anchor - vector3(pW->x, pW->y, pW->z)).length();

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

        if (min_dis > 36){ // reject water if too far
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
    std::cout << " Total waters " << totalwaters << std::endl;

    totalwatersInExcludedVolume = hydration.size() - totalwaters;

    delete[] tempWaters;

    SASTOOLS_UTILS_H::logger("TOTAL WATERS IN HYDRATION", std::to_string(totalwaters));
    SASTOOLS_UTILS_H::logger("TOTAL WATERS IN HYDRATION", std::to_string(totalwaters));

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
    std::string atomType, resiname;
    int count = 1;
    for(auto & line : waterLines){
        // if water is oxygen,
        if (boost::regex_match(line.substr(17,3), wat)){

            x = std::strtof(line.substr(30,8).c_str(), nullptr);
            y = std::strtof(line.substr(38,8).c_str(), nullptr);
            z = std::strtof(line.substr(46,8).c_str(), nullptr);

            //Atom type taken from rows 13-16 needs to be converted to the actual atom, e.g., CB1 is C
            atomType = line.substr(12,4);// needs to be this way until atomic numbers are assigned
            resiname = line.substr(17,3);

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
    return totalwaters*functions_WETSAXS::asf_at_q_zero(99);
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
        float sum_x=0, sum_y=0, sum_z=0;
        float at_x, at_y, at_z;
        float min_x=FLT_MAX, max_x = -FLT_MAX;
        float min_y=FLT_MAX, max_y = -FLT_MAX;
        float min_z=FLT_MAX, max_z = -FLT_MAX;

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

        // need to add backbone for RNA - center on C4'
        //if (false){
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
                } // else {
//                    std::string msg = " _"+ tempAtom + "_ " + std::to_string(resid);
//                    SASTOOLS_UTILS_H::logger("ALIGNMENT ATOM TYPE NOT FOUND", msg);
//                    SASTOOLS_UTILS_H::logger("RENAME ATOM TYPE TO STANDARD", "i.e, C5M => C7");
//                }
            }

            // create alignment vectors
            float aveX = xpos[centering_index];
            float aveY = ypos[centering_index];
            float aveZ = zpos[centering_index];

            // center the residue in the structure
            auto totalAtomsInResidue = tar_x.size() ;//(int)sccount;

            for(int i=0; i < totalAtomsInResidue;i++){
                tar_x[i] = tar_x[i] - aveX;
                tar_y[i] = tar_y[i] - aveY;
                tar_z[i] = tar_z[i] - aveZ;
            }

            // rotate library residue into structure residue
            setRotationMatrix(totalAtomsInResidue, &ref_x[0], &ref_y[0], &ref_z[0], &tar_x[0], &tar_y[0], &tar_z[0]);
            //cout << "____HYDRATE___" << endl;
            int totalToAdd = placeWaters(residue, resid, aveX, aveY, aveZ);
            // copy into vector
            for(int m=0; m<totalToAdd; m++){
                hydration.emplace_back( Coords(tempWaters[m]) );
            }
        }

        // if nucleic, need to hydrate
        //if (false){
        if (residue == "rA" || residue == "rG" || residue == "rU" || residue == "rC" || residue == "rI" ||
                residue == "DA" || residue == "DG" || residue == "DT" || residue == "DC"){

            std::string bck = (residue == "DA" || residue == "DG" || residue == "DT" || residue == "DC") ? "DNA" : "RNA";

            findIt = sideChains.find(bck);
            centering_atom = centeringAtoms.find(bck)->second;

            auto pBackbone = &(findIt->second); // std::unordered_map<std::string, Coords> >sideChains
            unsigned int centering_index=startAt, index;

            for (int i=0; i< atomsInResidue; i++) {

                index = startAt + i;
                tempAtom = atomType[index];
                boost::algorithm::trim (tempAtom);

                if (centering_atom == tempAtom){
                    centering_index = index;
                }

                auto pAtom = pBackbone->find(tempAtom);
                if (pAtom != pBackbone->end()){ //
                    // DA, DT, DG, DC = P, OP1, OP2 and sugar atoms
                    // Protein residues N, CA, CB, C, O
                    Coords *ptempCoord = &(*pAtom).second;
                    index = startAt + i;
                    tar_x.push_back(xpos[index]); // atoms in residue in structure
                    tar_y.push_back(ypos[index]);
                    tar_z.push_back(zpos[index]);

                    ref_x.push_back(ptempCoord->x); // atoms in residue in library
                    ref_y.push_back(ptempCoord->y);
                    ref_z.push_back(ptempCoord->z);
                }
            }

            // create alignment vectors
            float aveX = xpos[centering_index];
            float aveY = ypos[centering_index];
            float aveZ = zpos[centering_index];

            // center the residue in the structure
            auto totalAtomsInResidue = tar_x.size() ;//(int)sccount;

            for(int i=0; i < totalAtomsInResidue;i++){
                tar_x[i] = tar_x[i] - aveX;
                tar_y[i] = tar_y[i] - aveY;
                tar_z[i] = tar_z[i] - aveZ;
            }

            // rotate library residue into structure residue
            setRotationMatrix(totalAtomsInResidue, &ref_x[0], &ref_y[0], &ref_z[0], &tar_x[0], &tar_y[0], &tar_z[0]);
            int totalToAdd = placeWaters(bck, resid, aveX, aveY, aveZ);
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
        thetas[n] = acosf(zval/rval);
        phis[n] = (xval == 0 || xval == 0.0f) ? 0.0f : atan2f(yval, xval); // need to check if this is a good way to test for 0
    }

//    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//    std::cout << "THREADS elapsed time= s :: "  << runtime << std::endl;

}


void Waters::calculatePartialAmplitudes(unsigned int lmax,
                                        unsigned int totalqvalues,
                                        std::vector < float > & qvalues,
                                        bool recalculate){

    if (recalculate){
        sbj.recalculate(totalqvalues, qvalues, &rvalues);
    } else {
        this->createSphericalHarmonics(lmax);
        this->createSphericalBessels(totalqvalues, qvalues);
    }


    // total partial amplitudes per q-value is (lmax+1)^2
    almRealSum.resize((lmax+1)*(lmax+1)*totalqvalues);
    almImagSum.resize((lmax+1)*(lmax+1)*totalqvalues);
    float * const ptrR = (totalwaters != 0) ? almRealSum.data() : nullptr;
    float * const ptrI = (totalwaters != 0) ? almImagSum.data() : nullptr;

    auto pSBJ = sbj.getpSphericalBessels();

    auto pRealYlm = she.getpDataYLMReal();
    auto pImagYlm = she.getpDataYLMImag();

    std::clock_t startTime = std::clock();

    unsigned int bessel_index, q_index=0, base_index;
    unsigned int long ylm_index;
    float bessel_product, asf_at_q, partialSumR, partialSumI;

    float occupancies [totalwaters];
    for (unsigned int i=0; i< totalwaters; i++){
        occupancies[i] = 1.0;//hydration[i].occ;
    }

    float * pR, *pI;

    for(auto & qvalue : qvalues){

        asf_at_q = functions_WETSAXS::asf(99, qvalue);

        int ll_index = (lmax+1)*(lmax+1)*q_index;// + l*l + l + m;

        for(int l=0; l<=lmax; l++){

            bessel_index = (lmax+1)*totalwaters*q_index + totalwaters*l;
            base_index = (l*l)*totalwaters;

            /*
             * for m = 0
             */
            // for given (l,m) calculate Spherical Harmonic
            partialSumR = 0.0f;
            partialSumI = 0.0f;

            ylm_index = base_index + l * totalwaters; // for m = 0

            for (unsigned int n = 0; n < totalwaters; n++) { // over all atoms

                bessel_product = pSBJ[bessel_index + n] * occupancies[n];

                // real term
                partialSumR += bessel_product * pRealYlm[ylm_index];
                //imaginary term
                partialSumI += bessel_product * pImagYlm[ylm_index];

                ++ylm_index;
            }

            // update for m = 0
            int partials_index = ll_index + l*l + l;

            ptrR[partials_index] = asf_at_q * partialSumR;
            ptrI[partials_index] = asf_at_q * partialSumI;

            for(int m=1; m<=l; m++){
                // for given (l,m) calculate Spherical Harmonic
                // she.calculateSHE(l,m,thetas.data(), phis.data());

                // over all atoms
                partialSumR = 0.0f;
                partialSumI = 0.0f;

                ylm_index = base_index + (l+m) * totalwaters;

                // for each (l,m) calculate partial sum at given q value
                for(unsigned int n=0;n<totalwaters;n++){

                    bessel_product = pSBJ[bessel_index + n] * occupancies[n];

                    // real term
                    partialSumR += bessel_product * pRealYlm[ylm_index];
                    //imaginary term
                    partialSumI += bessel_product * pImagYlm[ylm_index];

                    ++ylm_index;
                }

                int at_m = partials_index + m;

                pR = &ptrR[at_m];
                pI = &ptrI[at_m];

                *pR = asf_at_q * partialSumR;
                *pI = asf_at_q * partialSumI;

                at_m = partials_index - m;
                // update negative m values
                if (m & 1){
                    ptrR[at_m] = -*pR;
                    ptrI[at_m] = *pI;
                } else {
                    ptrR[at_m] = *pR;
                    ptrI[at_m] = -*pI;
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


