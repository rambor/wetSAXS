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
     * add nucleics
     */
    sideChains[" DA"].insert({"P",Coords(-0.4995,-3.9442,-2.3020,"P",1.0)});
    sideChains[" DA"].insert({"OP1",Coords(0.0725,-4.6312,-3.4950,"OP1",1.0)});
    sideChains[" DA"].insert({"OP2",Coords(-0.7975,-4.7732,-1.1100,"OP2",1.0)});
    sideChains[" DA"].insert({"O5'",Coords(0.4835,-2.7762,-1.8460,"O5'",1.0)});
    sideChains[" DA"].insert({"C5'",Coords(0.9215,-1.8032,-2.7890,"C5'",1.0)});
    sideChains[" DA"].insert({"C4'",Coords(1.6335,-0.6932,-2.0630,"C4'",1.0)});
    sideChains[" DA"].insert({"O4'",Coords(0.6665,0.0788,-1.3140,"O4'",1.0)});
    sideChains[" DA"].insert({"C3'",Coords(2.6635,-1.1682,-1.0330,"C3'",1.0)});
    sideChains[" DA"].insert({"O3'",Coords(3.7695,-0.2642,-1.0910,"O3'",1.0)});
    sideChains[" DA"].insert({"C2'",Coords(1.9405,-1.0372,0.2830,"C2'",1.0)});
    sideChains[" DA"].insert({"C1'",Coords(1.1765,0.2568,0.0160,"C1'",1.0)});
    sideChains[" DA"].insert({"N9",Coords(0.0775,0.5448,0.9120,"N9",1.0)});
    sideChains[" DA"].insert({"C8",Coords(-0.7025,-0.2962,1.6580,"C8",1.0)});
    sideChains[" DA"].insert({"N7",Coords(-1.6725,0.3278,2.2820,"N7",1.0)});
    sideChains[" DA"].insert({"C5",Coords(-1.5155,1.6708,1.9410,"C5",1.0)});
    sideChains[" DA"].insert({"C6",Coords(-2.1875,2.8438,2.3250,"C6",1.0)});
    sideChains[" DA"].insert({"N6",Coords(-3.2455,2.8758,3.1270,"N6",1.0)});
    sideChains[" DA"].insert({"N1",Coords(-1.7255,3.9998,1.8170,"N1",1.0)});
    sideChains[" DA"].insert({"C2",Coords(-0.6735,4.0048,1.0060,"C2",1.0)});
    sideChains[" DA"].insert({"N3",Coords(0.0415,2.9608,0.5540,"N3",1.0)});
    sideChains[" DA"].insert({"C4",Coords(-0.4265,1.8238,1.1210,"C4",1.0)});
    sideChains[" DC"].insert({"P",Coords(0.3991,-3.3471,1.7177,"P",1.0)});
    sideChains[" DC"].insert({"OP1",Coords(0.9611,-4.6721,1.3597,"OP1",1.0)});
    sideChains[" DC"].insert({"OP2",Coords(0.9711,-2.6571,2.8837,"OP2",1.0)});
    sideChains[" DC"].insert({"O5'",Coords(0.3981,-2.4031,0.4277,"O5'",1.0)});
    sideChains[" DC"].insert({"C5'",Coords(1.4211,-2.5301,-0.5893,"C5'",1.0)});
    sideChains[" DC"].insert({"C4'",Coords(1.9761,-1.1541,-0.9153,"C4'",1.0)});
    sideChains[" DC"].insert({"O4'",Coords(0.9131,-0.3201,-1.4373,"O4'",1.0)});
    sideChains[" DC"].insert({"C3'",Coords(2.5411,-0.4101,0.2437,"C3'",1.0)});
    sideChains[" DC"].insert({"O3'",Coords(3.7001,0.3649,-0.0823,"O3'",1.0)});
    sideChains[" DC"].insert({"C2'",Coords(1.4331,0.5289,0.6497,"C2'",1.0)});
    sideChains[" DC"].insert({"C1'",Coords(0.8141,0.8769,-0.6943,"C1'",1.0)});
    sideChains[" DC"].insert({"N1",Coords(-0.5749,1.2569,-0.6023,"N1",1.0)});
    sideChains[" DC"].insert({"C2",Coords(-0.8669,2.6019,-0.3673,"C2",1.0)});
    sideChains[" DC"].insert({"O2",Coords(0.0661,3.4039,-0.2613,"O2",1.0)});
    sideChains[" DC"].insert({"N3",Coords(-2.1789,2.9529,-0.2733,"N3",1.0)});
    sideChains[" DC"].insert({"C4",Coords(-3.1499,2.0489,-0.4053,"C4",1.0)});
    sideChains[" DC"].insert({"N4",Coords(-4.4019,2.4499,-0.3023,"N4",1.0)});
    sideChains[" DC"].insert({"C5",Coords(-2.8589,0.6779,-0.6273,"C5",1.0)});
    sideChains[" DC"].insert({"C6",Coords(-1.5629,0.3299,-0.7243,"C6",1.0)});
    sideChains[" DG"].insert({"P",Coords(0.3656,4.0151,-2.9567,"P",1.0)});
    sideChains[" DG"].insert({"OP1",Coords(-0.5834,5.0801,-3.4077,"OP1",1.0)});
    sideChains[" DG"].insert({"OP2",Coords(1.7496,4.3971,-2.5577,"OP2",1.0)});
    sideChains[" DG"].insert({"O5'",Coords(-0.2254,3.2081,-1.7297,"O5'",1.0)});
    sideChains[" DG"].insert({"C5'",Coords(-1.5884,2.7551,-1.7977,"C5'",1.0)});
    sideChains[" DG"].insert({"C4'",Coords(-1.8954,2.0131,-0.5307,"C4'",1.0)});
    sideChains[" DG"].insert({"O4'",Coords(-1.1984,0.7371,-0.5657,"O4'",1.0)});
    sideChains[" DG"].insert({"C3'",Coords(-1.4004,2.7131,0.7383,"C3'",1.0)});
    sideChains[" DG"].insert({"O3'",Coords(-2.4124,2.8251,1.7233,"O3'",1.0)});
    sideChains[" DG"].insert({"C2'",Coords(-0.2764,1.8331,1.2243,"C2'",1.0)});
    sideChains[" DG"].insert({"C1'",Coords(-0.6854,0.4631,0.7343,"C1'",1.0)});
    sideChains[" DG"].insert({"N9",Coords(0.3716,-0.5089,0.5693,"N9",1.0)});
    sideChains[" DG"].insert({"C8",Coords(1.6166,-0.3199,0.0343,"C8",1.0)});
    sideChains[" DG"].insert({"N7",Coords(2.3216,-1.4019,-0.0317,"N7",1.0)});
    sideChains[" DG"].insert({"C5",Coords(1.4806,-2.4009,0.4473,"C5",1.0)});
    sideChains[" DG"].insert({"C6",Coords(1.6906,-3.7999,0.6503,"C6",1.0)});
    sideChains[" DG"].insert({"O6",Coords(2.6906,-4.4669,0.3493,"O6",1.0)});
    sideChains[" DG"].insert({"N1",Coords(0.5936,-4.4299,1.2293,"N1",1.0)});
    sideChains[" DG"].insert({"C2",Coords(-0.5724,-3.7859,1.5713,"C2",1.0)});
    sideChains[" DG"].insert({"N2",Coords(-1.5194,-4.5709,2.1003,"N2",1.0)});
    sideChains[" DG"].insert({"N3",Coords(-0.8024,-2.4959,1.3583,"N3",1.0)});
    sideChains[" DG"].insert({"C4",Coords(0.2796,-1.8599,0.8483,"C4",1.0)});
    sideChains[" DT"].insert({"P",Coords(-1.8285,3.1266,0.1304,"P",1.0)});
    sideChains[" DT"].insert({"OP1",Coords(-2.9125,3.7846,0.8844,"OP1",1.0)});
    sideChains[" DT"].insert({"OP2",Coords(-0.4465,3.6167,0.3204,"OP2",1.0)});
    sideChains[" DT"].insert({"O5'",Coords(-1.8055,1.5766,0.4734,"O5'",1.0)});
    sideChains[" DT"].insert({"C5'",Coords(-2.9225,0.7286,0.2504,"C5'",1.0)});
    sideChains[" DT"].insert({"C4'",Coords(-2.5455,-0.6944,0.6204,"C4'",1.0)});
    sideChains[" DT"].insert({"O4'",Coords(-1.4535,-1.1464,-0.2066,"O4'",1.0)});
    sideChains[" DT"].insert({"C3'",Coords(-2.0705,-0.8824,2.0574,"C3'",1.0)});
    sideChains[" DT"].insert({"O3'",Coords(-2.6215,-2.1244,2.4834,"O3'",1.0)});
    sideChains[" DT"].insert({"C2'",Coords(-0.5705,-0.8974,1.9354,"C2'",1.0)});
    sideChains[" DT"].insert({"C1'",Coords(-0.3405,-1.5454,0.5914,"C1'",1.0)});
    sideChains[" DT"].insert({"N1",Coords(0.8845,-1.1274,-0.0946,"N1",1.0)});
    sideChains[" DT"].insert({"C2",Coords(1.7195,-2.0884,-0.6076,"C2",1.0)});
    sideChains[" DT"].insert({"O2",Coords(1.5835,-3.2934,-0.4006,"O2",1.0)});
    sideChains[" DT"].insert({"N3",Coords(2.7325,-1.6284,-1.4026,"N3",1.0)});
    sideChains[" DT"].insert({"C4",Coords(3.0215,-0.3154,-1.7026,"C4",1.0)});
    sideChains[" DT"].insert({"O4",Coords(4.0015,-0.0584,-2.3976,"O4",1.0)});
    sideChains[" DT"].insert({"C5",Coords(2.1115,0.6526,-1.1426,"C5",1.0)});
    sideChains[" DT"].insert({"C7",Coords(2.3805,2.1086,-1.3686,"C7",1.0)});
    sideChains[" DT"].insert({"C6",Coords(1.0825,0.2066,-0.4236,"C6",1.0)});
    sideChains[" rA"].insert({"P",Coords(2.0505,1.9103,-3.6766,"P",1.0)});
    sideChains[" rA"].insert({"OP1",Coords(1.8505,1.4763,-5.1006,"OP1",1.0)});
    sideChains[" rA"].insert({"OP2",Coords(2.8375,3.1193,-3.3926,"OP2",1.0)});
    sideChains[" rA"].insert({"O5'",Coords(0.5875,1.9923,-3.0616,"O5'",1.0)});
    sideChains[" rA"].insert({"C5'",Coords(0.2605,2.8263,-1.9476,"C5'",1.0)});
    sideChains[" rA"].insert({"C4'",Coords(-0.8545,2.1893,-1.1196,"C4'",1.0)});
    sideChains[" rA"].insert({"O4'",Coords(-0.4985,0.8523,-0.7396,"O4'",1.0)});
    sideChains[" rA"].insert({"C3'",Coords(-1.2065,2.9593,0.1624,"C3'",1.0)});
    sideChains[" rA"].insert({"O3'",Coords(-2.6415,3.0063,0.2684,"O3'",1.0)});
    sideChains[" rA"].insert({"C2'",Coords(-0.8455,2.0163,1.3064,"C2'",1.0)});
    sideChains[" rA"].insert({"O2'",Coords(-1.7575,2.0023,2.3794,"O2'",1.0)});
    sideChains[" rA"].insert({"C1'",Coords(-0.8705,0.6643,0.6004,"C1'",1.0)});
    sideChains[" rA"].insert({"N9",Coords(-0.0925,-0.4037,1.1904,"N9",1.0)});
    sideChains[" rA"].insert({"C8",Coords(0.9965,-0.3517,2.0134,"C8",1.0)});
    sideChains[" rA"].insert({"N7",Coords(1.4285,-1.5287,2.3914,"N7",1.0)});
    sideChains[" rA"].insert({"C5",Coords(0.5595,-2.4087,1.7684,"C5",1.0)});
    sideChains[" rA"].insert({"C6",Coords(0.4775,-3.8007,1.7754,"C6",1.0)});
    sideChains[" rA"].insert({"N6",Coords(1.3175,-4.5707,2.4504,"N6",1.0)});
    sideChains[" rA"].insert({"N1",Coords(-0.5055,-4.3757,1.0554,"N1",1.0)});
    sideChains[" rA"].insert({"C2",Coords(-1.3465,-3.5847,0.3704,"C2",1.0)});
    sideChains[" rA"].insert({"N3",Coords(-1.3725,-2.2587,0.2854,"N3",1.0)});
    sideChains[" rA"].insert({"C4",Coords(-0.3755,-1.7307,1.0214,"C4",1.0)});
    sideChains[" rC"].insert({"P",Coords(2.2922,2.6740,1.5097,"P",1.0)});
    sideChains[" rC"].insert({"OP1",Coords(2.0632,3.9129,2.3257,"OP1",1.0)});
    sideChains[" rC"].insert({"OP2",Coords(3.0612,1.5729,2.0677,"OP2",1.0)});
    sideChains[" rC"].insert({"O5'",Coords(0.8802,2.0759,1.0667,"O5'",1.0)});
    sideChains[" rC"].insert({"C5'",Coords(-0.0908,2.8959,0.4517,"C5'",1.0)});
    sideChains[" rC"].insert({"C4'",Coords(-1.1358,2.0319,-0.1853,"C4'",1.0)});
    sideChains[" rC"].insert({"O4'",Coords(-0.4788,1.2099,-1.1753,"O4'",1.0)});
    sideChains[" rC"].insert({"C3'",Coords(-1.7838,1.0140,0.7417,"C3'",1.0)});
    sideChains[" rC"].insert({"O3'",Coords(-2.8258,1.6180,1.5107,"O3'",1.0)});
    sideChains[" rC"].insert({"C2'",Coords(-2.3258,0.0190,-0.2703,"C2'",1.0)});
    sideChains[" rC"].insert({"O2'",Coords(-3.4908,0.5329,-0.9193,"O2'",1.0)});
    sideChains[" rC"].insert({"C1'",Coords(-1.1368,-0.0551,-1.2443,"C1'",1.0)});
    sideChains[" rC"].insert({"N1",Coords(-0.1638,-1.1071,-0.9003,"N1",1.0)});
    sideChains[" rC"].insert({"C2",Coords(-0.4498,-2.4031,-1.3273,"C2",1.0)});
    sideChains[" rC"].insert({"O2",Coords(-1.5098,-2.5811,-1.9583,"O2",1.0)});
    sideChains[" rC"].insert({"N3",Coords(0.4042,-3.4061,-1.0543,"N3",1.0)});
    sideChains[" rC"].insert({"C4",Coords(1.5272,-3.1541,-0.3783,"C4",1.0)});
    sideChains[" rC"].insert({"N4",Coords(2.3342,-4.1751,-0.1593,"N4",1.0)});
    sideChains[" rC"].insert({"C5",Coords(1.8512,-1.8291,0.0917,"C5",1.0)});
    sideChains[" rC"].insert({"C6",Coords(0.9782,-0.8471,-0.1933,"C6",1.0)});
    sideChains[" rG"].insert({"P",Coords(-1.8214,4.4453,0.1367,"P",1.0)});
    sideChains[" rG"].insert({"OP1",Coords(-3.0364,5.0523,0.7317,"OP1",1.0)});
    sideChains[" rG"].insert({"OP2",Coords(-0.4834,4.8353,0.6097,"OP2",1.0)});
    sideChains[" rG"].insert({"O5'",Coords(-1.8764,2.8613,0.2327,"O5'",1.0)});
    sideChains[" rG"].insert({"C5'",Coords(-3.0684,2.1593,-0.1283,"C5'",1.0)});
    sideChains[" rG"].insert({"C4'",Coords(-2.8154,0.6813,-0.1063,"C4'",1.0)});
    sideChains[" rG"].insert({"O4'",Coords(-1.7564,0.3243,-1.0533,"O4'",1.0)});
    sideChains[" rG"].insert({"C3'",Coords(-2.2214,0.2013,1.1967,"C3'",1.0)});
    sideChains[" rG"].insert({"O3'",Coords(-3.2694,0.1423,2.1577,"O3'",1.0)});
    sideChains[" rG"].insert({"C2'",Coords(-1.6644,-1.1607,0.7797,"C2'",1.0)});
    sideChains[" rG"].insert({"O2'",Coords(-2.7104,-2.1427,0.6647,"O2'",1.0)});
    sideChains[" rG"].insert({"C1'",Coords(-1.0614,-0.8077,-0.5863,"C1'",1.0)});
    sideChains[" rG"].insert({"N9",Coords(0.3596,-0.4777,-0.4893,"N9",1.0)});
    sideChains[" rG"].insert({"C8",Coords(0.9406,0.7613,-0.4313,"C8",1.0)});
    sideChains[" rG"].insert({"N7",Coords(2.2436,0.6993,-0.4023,"N7",1.0)});
    sideChains[" rG"].insert({"C5",Coords(2.5356,-0.6437,-0.4183,"C5",1.0)});
    sideChains[" rG"].insert({"C6",Coords(3.7526,-1.3067,-0.3733,"C6",1.0)});
    sideChains[" rG"].insert({"O6",Coords(4.8886,-0.8297,-0.3313,"O6",1.0)});
    sideChains[" rG"].insert({"N1",Coords(3.5996,-2.6657,-0.3763,"N1",1.0)});
    sideChains[" rG"].insert({"C2",Coords(2.3986,-3.3427,-0.4313,"C2",1.0)});
    sideChains[" rG"].insert({"N2",Coords(2.4426,-4.6577,-0.4563,"N2",1.0)});
    sideChains[" rG"].insert({"N3",Coords(1.2416,-2.7327,-0.4663,"N3",1.0)});
    sideChains[" rG"].insert({"C4",Coords(1.3816,-1.3947,-0.4603,"C4",1.0)});
    sideChains[" rU"].insert({"P",Coords(-3.2431,2.0657,-0.4487,"P",1.0)});
    sideChains[" rU"].insert({"OP1",Coords(-4.5471,2.2078,0.2153,"OP1",1.0)});
    sideChains[" rU"].insert({"OP2",Coords(-2.2751,3.1598,-0.3757,"OP2",1.0)});
    sideChains[" rU"].insert({"O5'",Coords(-2.5061,0.7487,0.0603,"O5'",1.0)});
    sideChains[" rU"].insert({"C5'",Coords(-3.0681,-0.5743,-0.1087,"C5'",1.0)});
    sideChains[" rU"].insert({"C4'",Coords(-2.0031,-1.5802,0.2773,"C4'",1.0)});
    sideChains[" rU"].insert({"O4'",Coords(-0.8541,-1.3462,-0.5647,"O4'",1.0)});
    sideChains[" rU"].insert({"C3'",Coords(-1.5691,-1.3483,1.7303,"C3'",1.0)});
    sideChains[" rU"].insert({"O3'",Coords(-1.3651,-2.6412,2.3523,"O3'",1.0)});
    sideChains[" rU"].insert({"C2'",Coords(-0.2071,-0.6863,1.5923,"C2'",1.0)});
    sideChains[" rU"].insert({"O2'",Coords(0.7069,-0.8663,2.6543,"O2'",1.0)});
    sideChains[" rU"].insert({"C1'",Coords(0.2849,-1.2363,0.2443,"C1'",1.0)});
    sideChains[" rU"].insert({"N1",Coords(1.2789,-0.4022,-0.4287,"N1",1.0)});
    sideChains[" rU"].insert({"C2",Coords(2.5039,-0.9902,-0.7467,"C2",1.0)});
    sideChains[" rU"].insert({"O2",Coords(2.7309,-2.1853,-0.6687,"O2",1.0)});
    sideChains[" rU"].insert({"N3",Coords(3.4539,-0.1223,-1.1467,"N3",1.0)});
    sideChains[" rU"].insert({"C4",Coords(3.3179,1.2308,-1.3167,"C4",1.0)});
    sideChains[" rU"].insert({"O4",Coords(4.3109,1.8878,-1.6157,"O4",1.0)});
    sideChains[" rU"].insert({"C5",Coords(2.0039,1.7498,-1.0637,"C5",1.0)});
    sideChains[" rU"].insert({"C6",Coords(1.0449,0.9288,-0.6417,"C6",1.0)});

    // any new residue requires set of centered coordinates - maybe averaged followed by a set of hydrating waters

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
    waters[" rA"] = std::vector<Coords>(16);
    pwater = &waters[" rA"];
    (*pwater)[0].x = -2.3151;
    (*pwater)[0].y = -1.9034;
    (*pwater)[0].z = 6.6707;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = -2.9978;
    (*pwater)[1].y = 5.1559;
    (*pwater)[1].z = 1.1633;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 3.8593;
    (*pwater)[2].y = 3.4099;
    (*pwater)[2].z = 3.7829;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 2.4351;
    (*pwater)[3].y = 6.5518;
    (*pwater)[3].z = 4.3852;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.6670;
    (*pwater)[4].x = -0.9107;
    (*pwater)[4].y = 7.1066;
    (*pwater)[4].z = 1.0449;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -2.7475;
    (*pwater)[5].y = 0.7178;
    (*pwater)[5].z = 3.7877;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = -4.6706;
    (*pwater)[6].y = 1.7970;
    (*pwater)[6].z = -1.1000;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.6670;
    (*pwater)[7].x = -0.4712;
    (*pwater)[7].y = 7.7070;
    (*pwater)[7].z = 3.9309;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 4.6818;
    (*pwater)[8].y = -1.5165;
    (*pwater)[8].z = 1.9604;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.6670;
    (*pwater)[9].x = 3.1355;
    (*pwater)[9].y = 1.3777;
    (*pwater)[9].z = 1.3600;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = -3.7317;
    (*pwater)[10].y = 2.2702;
    (*pwater)[10].z = 3.5079;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 5.2641;
    (*pwater)[11].y = 3.2647;
    (*pwater)[11].z = 1.5906;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = 1.6179;
    (*pwater)[12].y = 7.7585;
    (*pwater)[12].z = 3.1372;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 0.3330;
    (*pwater)[13].x = -0.9632;
    (*pwater)[13].y = -0.4734;
    (*pwater)[13].z = 8.7067;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.3330;
    (*pwater)[14].x = 6.0020;
    (*pwater)[14].y = -5.3822;
    (*pwater)[14].z = 2.2091;
    (*pwater)[14].type = "O";
    (*pwater)[14].occ = 0.3330;
    (*pwater)[15].x = 4.4519;
    (*pwater)[15].y = -2.2909;
    (*pwater)[15].z = -0.3282;
    (*pwater)[15].type = "O";
    (*pwater)[15].occ = 0.3330;

// Create Vector for new Residue G
    waters[" rG"] = std::vector<Coords>(21);
    pwater = &waters[" rG"];
    (*pwater)[0].x = -5.0474;
    (*pwater)[0].y = 0.6048;
    (*pwater)[0].z = 5.5357;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 0.1009;
    (*pwater)[1].y = 0.2850;
    (*pwater)[1].z = 6.3765;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -1.3157;
    (*pwater)[2].y = -4.3911;
    (*pwater)[2].z = -0.3592;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = -4.3172;
    (*pwater)[3].y = -4.2792;
    (*pwater)[3].z = 1.7921;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 6.8415;
    (*pwater)[4].y = 2.8652;
    (*pwater)[4].z = 2.9403;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 9.8814;
    (*pwater)[5].y = -1.0003;
    (*pwater)[5].z = -0.9952;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = -0.0358;
    (*pwater)[6].y = -5.3617;
    (*pwater)[6].z = -3.8965;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.6670;
    (*pwater)[7].x = 0.7130;
    (*pwater)[7].y = 2.0723;
    (*pwater)[7].z = 3.7056;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = 0.2772;
    (*pwater)[8].y = -6.5138;
    (*pwater)[8].z = -0.2537;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = 2.8151;
    (*pwater)[9].y = -7.9146;
    (*pwater)[9].z = 1.9377;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;
    (*pwater)[10].x = -6.8175;
    (*pwater)[10].y = 0.2770;
    (*pwater)[10].z = 1.1481;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = 3.3244;
    (*pwater)[11].y = 1.9211;
    (*pwater)[11].z = 3.2390;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 1.0000;
    (*pwater)[12].x = -5.4723;
    (*pwater)[12].y = -1.8342;
    (*pwater)[12].z = 1.4139;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 1.0000;
    (*pwater)[13].x = -1.2037;
    (*pwater)[13].y = -7.1516;
    (*pwater)[13].z = 2.7445;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.6670;
    (*pwater)[14].x = 8.5251;
    (*pwater)[14].y = -3.5069;
    (*pwater)[14].z = -3.9051;
    (*pwater)[14].type = "O";
    (*pwater)[14].occ = 0.3330;
    (*pwater)[15].x = 10.4743;
    (*pwater)[15].y = 0.6524;
    (*pwater)[15].z = 1.1108;
    (*pwater)[15].type = "O";
    (*pwater)[15].occ = 0.3330;
    (*pwater)[16].x = -0.1180;
    (*pwater)[16].y = 4.3302;
    (*pwater)[16].z = 5.7508;
    (*pwater)[16].type = "O";
    (*pwater)[16].occ = 0.3330;
    (*pwater)[17].x = 7.3654;
    (*pwater)[17].y = 1.0194;
    (*pwater)[17].z = 4.3941;
    (*pwater)[17].type = "O";
    (*pwater)[17].occ = 0.3330;
    (*pwater)[18].x = -6.2170;
    (*pwater)[18].y = -2.8016;
    (*pwater)[18].z = -1.1848;
    (*pwater)[18].type = "O";
    (*pwater)[18].occ = 0.3330;
    (*pwater)[19].x = 6.6659;
    (*pwater)[19].y = -1.5060;
    (*pwater)[19].z = 2.9142;
    (*pwater)[19].type = "O";
    (*pwater)[19].occ = 0.3330;
    (*pwater)[20].x = -6.8774;
    (*pwater)[20].y = -3.1975;
    (*pwater)[20].z = 4.3375;
    (*pwater)[20].type = "O";
    (*pwater)[20].occ = 0.3330;

// Create Vector for new Residue U
    waters[" rU"] = std::vector<Coords>(12);
    pwater = &waters[" rU"];
    (*pwater)[0].x = 2.6331;
    (*pwater)[0].y = 3.3829;
    (*pwater)[0].z = 2.0980;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 1.4436;
    (*pwater)[1].y = 2.2385;
    (*pwater)[1].z = 5.6663;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -2.9314;
    (*pwater)[2].y = 2.6933;
    (*pwater)[2].z = 4.6556;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.6670;
    (*pwater)[3].x = -2.4196;
    (*pwater)[3].y = -4.0484;
    (*pwater)[3].z = 3.7483;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.6670;
    (*pwater)[4].x = 2.2664;
    (*pwater)[4].y = -4.7059;
    (*pwater)[4].z = 1.4133;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 1.9806;
    (*pwater)[5].y = 5.1754;
    (*pwater)[5].z = 0.4960;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = -1.5934;
    (*pwater)[6].y = -1.3653;
    (*pwater)[6].z = 7.8245;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = 4.1326;
    (*pwater)[7].y = 5.9174;
    (*pwater)[7].z = 0.9400;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.6670;
    (*pwater)[8].x = -0.4232;
    (*pwater)[8].y = -5.4076;
    (*pwater)[8].z = 5.0181;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 1.6060;
    (*pwater)[9].y = -3.6626;
    (*pwater)[9].z = 3.3271;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.6670;
    (*pwater)[10].x = -4.7927;
    (*pwater)[10].y = -0.2808;
    (*pwater)[10].z = 4.6545;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = -5.3728;
    (*pwater)[11].y = 1.2300;
    (*pwater)[11].z = 4.5763;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;

// Create Vector for new Residue C
    waters[" rC"] = std::vector<Coords>(17);
    pwater = &waters[" rC"];
    (*pwater)[0].x = 1.3207;
    (*pwater)[0].y = -2.4963;
    (*pwater)[0].z = 3.4978;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -0.6073;
    (*pwater)[1].y = 2.3336;
    (*pwater)[1].z = 5.5621;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -5.4325;
    (*pwater)[2].y = 3.5253;
    (*pwater)[2].z = 1.9739;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = -1.8222;
    (*pwater)[3].y = -1.1333;
    (*pwater)[3].z = 5.6276;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = -4.9458;
    (*pwater)[4].y = 2.8218;
    (*pwater)[4].z = -0.4615;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = -5.7648;
    (*pwater)[5].y = 1.4829;
    (*pwater)[5].z = 5.7102;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 0.9915;
    (*pwater)[6].y = -1.1308;
    (*pwater)[6].z = 6.5248;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -4.6099;
    (*pwater)[7].y = -1.6508;
    (*pwater)[7].z = -2.1170;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = -5.3239;
    (*pwater)[8].y = -4.6550;
    (*pwater)[8].z = -2.9277;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = 4.2943;
    (*pwater)[9].y = -5.0980;
    (*pwater)[9].z = 3.2428;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.6670;
    (*pwater)[10].x = 4.9298;
    (*pwater)[10].y = -9.0100;
    (*pwater)[10].z = 1.5053;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = -2.5827;
    (*pwater)[11].y = 5.3573;
    (*pwater)[11].z = 2.5159;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = 1.0438;
    (*pwater)[12].y = 1.0153;
    (*pwater)[12].z = 6.7282;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 0.6670;
    (*pwater)[13].x = 3.3454;
    (*pwater)[13].y = -6.4521;
    (*pwater)[13].z = 3.2185;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.3330;
    (*pwater)[14].x = 3.5398;
    (*pwater)[14].y = -7.6049;
    (*pwater)[14].z = 0.6300;
    (*pwater)[14].type = "O";
    (*pwater)[14].occ = 0.3330;
    (*pwater)[15].x = -5.8325;
    (*pwater)[15].y = -0.7430;
    (*pwater)[15].z = -1.5070;
    (*pwater)[15].type = "O";
    (*pwater)[15].occ = 0.3330;
    (*pwater)[16].x = 4.9845;
    (*pwater)[16].y = -3.6332;
    (*pwater)[16].z = 0.9446;
    (*pwater)[16].type = "O";
    (*pwater)[16].occ = 0.3330;


// Create Vector for new Residue DA
    waters[" DA"] = std::vector<Coords>(11);
    pwater = &waters[" DA"];
    (*pwater)[0].x = 3.9244;
    (*pwater)[0].y = 2.8812;
    (*pwater)[0].z = -2.4520;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 0.6670;
    (*pwater)[1].x = 4.9959;
    (*pwater)[1].y = -1.2603;
    (*pwater)[1].z = 3.1532;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 0.3329;
    (*pwater)[2].y = -1.8341;
    (*pwater)[2].z = 4.3724;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.3330;
    (*pwater)[3].x = 5.3992;
    (*pwater)[3].y = 1.0315;
    (*pwater)[3].z = -3.7158;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 4.8103;
    (*pwater)[4].y = -3.4228;
    (*pwater)[4].z = -3.3444;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 2.3935;
    (*pwater)[5].y = 4.1370;
    (*pwater)[5].z = -0.2270;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 3.3967;
    (*pwater)[6].y = 6.9620;
    (*pwater)[6].z = 2.0570;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 7.4121;
    (*pwater)[7].y = 1.6179;
    (*pwater)[7].z = -4.1550;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 6.5471;
    (*pwater)[8].y = -3.8436;
    (*pwater)[8].z = 0.7699;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 0.3330;
    (*pwater)[9].x = 1.4105;
    (*pwater)[9].y = 8.1095;
    (*pwater)[9].z = 1.0448;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;
    (*pwater)[10].x = 7.3812;
    (*pwater)[10].y = -0.9858;
    (*pwater)[10].z = -4.0873;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;


// Create Vector for new Residue DC
    waters[" DC"] = std::vector<Coords>(12);
    pwater = &waters[" DC"];
    (*pwater)[0].x = 4.1711;
    (*pwater)[0].y = 2.3387;
    (*pwater)[0].z = -3.4908;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 4.9321;
    (*pwater)[1].y = 3.5497;
    (*pwater)[1].z = -0.5246;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 0.6670;
    (*pwater)[2].x = 2.3683;
    (*pwater)[2].y = 6.8861;
    (*pwater)[2].z = -1.6557;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = 2.1057;
    (*pwater)[3].y = 5.8861;
    (*pwater)[3].z = -5.1574;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 7.9822;
    (*pwater)[4].y = 1.4160;
    (*pwater)[4].z = -0.0039;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 4.4257;
    (*pwater)[5].y = 6.1779;
    (*pwater)[5].z = 0.4837;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 0.3330;
    (*pwater)[6].x = 8.1642;
    (*pwater)[6].y = -1.5826;
    (*pwater)[6].z = 1.8757;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = 1.2021;
    (*pwater)[7].y = -0.0178;
    (*pwater)[7].z = 4.7893;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = 4.5325;
    (*pwater)[8].y = -6.7553;
    (*pwater)[8].z = -2.0668;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = 3.4651;
    (*pwater)[9].y = 1.6782;
    (*pwater)[9].z = 4.0884;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 1.0000;
    (*pwater)[10].x = -3.3545;
    (*pwater)[10].y = 2.6379;
    (*pwater)[10].z = 4.1373;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.6670;
    (*pwater)[11].x = 5.1173;
    (*pwater)[11].y = -2.4093;
    (*pwater)[11].z = 3.4462;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;


// Create Vector for new Residue DG
    waters[" DG"] = std::vector<Coords>(14);
    pwater = &waters[" DG"];
    (*pwater)[0].x = -4.0134;
    (*pwater)[0].y = -4.1199;
    (*pwater)[0].z = 3.0242;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = -3.0232;
    (*pwater)[1].y = 6.5053;
    (*pwater)[1].z = 1.2605;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = 4.6464;
    (*pwater)[2].y = 1.2016;
    (*pwater)[2].z = 4.0452;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 1.0000;
    (*pwater)[3].x = -1.8470;
    (*pwater)[3].y = 3.9853;
    (*pwater)[3].z = 5.6001;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 1.0000;
    (*pwater)[4].x = 7.3519;
    (*pwater)[4].y = -6.2410;
    (*pwater)[4].z = 0.9215;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 5.1823;
    (*pwater)[5].y = -1.1410;
    (*pwater)[5].z = 2.3331;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 5.9589;
    (*pwater)[6].y = -8.4685;
    (*pwater)[6].z = 0.8028;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 0.3330;
    (*pwater)[7].x = -1.5398;
    (*pwater)[7].y = 6.5863;
    (*pwater)[7].z = 5.2076;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 1.0000;
    (*pwater)[8].x = -5.4612;
    (*pwater)[8].y = 0.7465;
    (*pwater)[8].z = 0.9575;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = -5.1384;
    (*pwater)[9].y = 3.8719;
    (*pwater)[9].z = -0.4148;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = -0.2801;
    (*pwater)[10].y = 5.5866;
    (*pwater)[10].z = 3.5442;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 0.3330;
    (*pwater)[11].x = -3.9657;
    (*pwater)[11].y = -0.7548;
    (*pwater)[11].z = 4.6576;
    (*pwater)[11].type = "O";
    (*pwater)[11].occ = 0.3330;
    (*pwater)[12].x = 4.1258;
    (*pwater)[12].y = 4.1804;
    (*pwater)[12].z = 3.8482;
    (*pwater)[12].type = "O";
    (*pwater)[12].occ = 1.0000;
    (*pwater)[13].x = -5.8264;
    (*pwater)[13].y = 3.3762;
    (*pwater)[13].z = 1.6902;
    (*pwater)[13].type = "O";
    (*pwater)[13].occ = 0.3330;


// Create Vector for new Residue DT
    waters[" DT"] = std::vector<Coords>(11);
    pwater = &waters[" DT"];
    (*pwater)[0].x = 2.0406;
    (*pwater)[0].y = -7.6477;
    (*pwater)[0].z = 1.5050;
    (*pwater)[0].type = "O";
    (*pwater)[0].occ = 1.0000;
    (*pwater)[1].x = 4.0482;
    (*pwater)[1].y = -7.8872;
    (*pwater)[1].z = 0.0667;
    (*pwater)[1].type = "O";
    (*pwater)[1].occ = 1.0000;
    (*pwater)[2].x = -5.5711;
    (*pwater)[2].y = -3.2725;
    (*pwater)[2].z = 2.7347;
    (*pwater)[2].type = "O";
    (*pwater)[2].occ = 0.6670;
    (*pwater)[3].x = 5.7182;
    (*pwater)[3].y = 1.7706;
    (*pwater)[3].z = 1.9058;
    (*pwater)[3].type = "O";
    (*pwater)[3].occ = 0.6670;
    (*pwater)[4].x = 7.4716;
    (*pwater)[4].y = 2.6815;
    (*pwater)[4].z = 0.2087;
    (*pwater)[4].type = "O";
    (*pwater)[4].occ = 1.0000;
    (*pwater)[5].x = 0.6752;
    (*pwater)[5].y = -0.6026;
    (*pwater)[5].z = 5.5922;
    (*pwater)[5].type = "O";
    (*pwater)[5].occ = 1.0000;
    (*pwater)[6].x = 3.4773;
    (*pwater)[6].y = 1.1435;
    (*pwater)[6].z = 3.2623;
    (*pwater)[6].type = "O";
    (*pwater)[6].occ = 1.0000;
    (*pwater)[7].x = -3.5632;
    (*pwater)[7].y = -5.9902;
    (*pwater)[7].z = 5.5048;
    (*pwater)[7].type = "O";
    (*pwater)[7].occ = 0.3330;
    (*pwater)[8].x = -3.5639;
    (*pwater)[8].y = -5.2633;
    (*pwater)[8].z = 1.1695;
    (*pwater)[8].type = "O";
    (*pwater)[8].occ = 1.0000;
    (*pwater)[9].x = 0.0230;
    (*pwater)[9].y = -2.7772;
    (*pwater)[9].z = 7.3155;
    (*pwater)[9].type = "O";
    (*pwater)[9].occ = 0.3330;
    (*pwater)[10].x = -5.5912;
    (*pwater)[10].y = -3.7620;
    (*pwater)[10].z = 5.7415;
    (*pwater)[10].type = "O";
    (*pwater)[10].occ = 1.0000;
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

//        if (i == (numAtoms-1)){
//            hydrateResidueDirect(currentRes, currentID, atomsInResidue, startHere, atomType, xvalue, yvalue, zvalue);
//        }
    }

    std::string residue_index;
    totalwaters = (unsigned int)hydration.size();

    //std::clock_t startTime = std::clock();
    // check for waters too close to each other - average their positions
    float dwater = 2.8*2.8;

    for(int w=0; w<totalwaters; w++){ // sort the waters.  waters too close will be moved to end of the vector
        auto pW = &hydration[w];
        vector3 vecW = vector3(pW->x, pW->y, pW->z);

        unsigned int k = w + 1;
        while(k < totalwaters){
            auto pK = &hydration[k];
            if ((vector3(pK->x, pK->y, pK->z) - vecW).sqlength() < dwater){
                pW->x = 0.5f*(pK->x + vecW.x);
                pW->y = 0.5f*(pK->y + vecW.y);
                pW->z = 0.5f*(pK->z + vecW.z);

                vecW = vector3(pW->x, pW->y, pW->z);
                std::iter_swap(hydration.begin()+k, hydration.begin()+(totalwaters-1));
                //std::cout << "Averaging water position " << std::endl;
                totalwaters--;
            } else {
                k++;
            }
        }
    }
//    double runtime = (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//    std::cout << " Averaging time " << runtime << std::endl;
    /*
     * prune waters too close to protein
     * does not consider hydrogens
     *
     */
    for (unsigned int i=0; i<numAtoms; i++){

        currentID = pResIDs[i];
        currentRes = pResidues[i];
        int atmnumber = atomicNumbers[i];
        float rlimit;

        auto it = minima.find(currentRes);
        if (it == minima.end()){
            rlimit = 2.2; //arbitrary
        } else {
            rlimit = it->second;
        }

        float xval = xvalue[i];
        float yval = yvalue[i];
        float zval = zvalue[i];

        vector3 anchor(xval, yval, zval);
        rlimit *= rlimit; // compare the square distance to reduce operations in the loop

        for(int w=0; w<totalwaters; w++){ // sort the waters.  waters too close will be moved to end of the vector
            auto pW = &hydration[w];

            if ((pW->resid != currentID) && (anchor - vector3(pW->x, pW->y, pW->z)).sqlength() < rlimit ){ // too close for hydrogen bond
                // reject
                std::iter_swap(hydration.begin()+w, hydration.begin()+(totalwaters-1));
                totalwaters--;
            }

            if (w >= totalwaters){
                break;
            }
        }
    }


    totalwatersInExcludedVolume = hydration.size() - totalwaters;

    delete[] tempWaters;

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
            atomType = line.substr(12,4) ;// needs to be this way until atomic numbers are assigned
            resiname = line.substr(17,3);
            if (model.getPDBModel().ifOxygen(atomType)){
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

    //float rotation[9];
    std::vector<float> tar_x = std::vector<float>();
    std::vector<float> tar_y = std::vector<float>();
    std::vector<float> tar_z = std::vector<float>();
    std::vector<float> ref_x = std::vector<float>();
    std::vector<float> ref_y = std::vector<float>();
    std::vector<float> ref_z = std::vector<float>();

    // apply rotation and add delx, dely, delz to waters

    if ( sideChains.find(residue) == sideChains.end() ) {
        // not found
        std::string note = "RESID " + std::to_string(resid) + " => " + residue;
        SASTOOLS_UTILS_H::logger("PSEUDO HYDRATING", note);

        boost::regex ifCarbon("^C[0-9]+?"); // match any character

        // collect atoms in the unknown residue, center, determine size of bounding box, increase size and form a rectangular grid
        float sum_x=0, sum_y=0, sum_z=0;
        float at_x, at_y, at_z;
        float min_x=FLT_MAX, max_x = -FLT_MAX;
        float min_y=FLT_MAX, max_y = -FLT_MAX;
        float min_z=FLT_MAX, max_z = -FLT_MAX;


        // random rotation
        std::random_device rd;
        std::mt19937 gen(rd());
        auto convert = (float)(M_PI/180.0f);
        std::uniform_int_distribution<int> randomIndex(0,360); // guaranteed unbiased
        std::uniform_int_distribution<int> randomBeta(0,180); // guaranteed unbiased

        std::vector <vector3> beads;

        // remove waters too close and within range
        // needs to be van der Waals radii plus 1.4 Angstrom (radius of water?)

        for (int i=0; i<atomsInResidue; i++) {

            tempAtom = atomType[startAt + i];

            // if atom type does not match C# where # can be any integer
            if (!boost::regex_match(tempAtom, ifCarbon)){

                float radius = cbrtf(0.75/M_PI*SASTOOLS_UTILS_H::residueToVolume(tempAtom, residue));

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

        for (int i=0; i<atomsInResidue; i++) {
            tempAtom = atomType[startAt + i];
            at_x = xpos[startAt + i];
            at_y = ypos[startAt + i];
            at_z = zpos[startAt + i];
            vector3 cvec = vector3(at_x, at_y, at_z);

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

            for(int k=w+1; k<totalBeads; k++){
                auto pNext = &beads[k];
                if ( ( *pNext - *pCurrent).sqlength() < limit){
                    *pCurrent = ( *pNext + *pCurrent)*0.5;
                    std::iter_swap(beads.begin()+k, beads.begin()+(totalBeads-1));
                    totalBeads--;
                    k--;
                }
            }

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

    } else {
        // get all atoms in the sideChain, includes backbone

        auto findIt = sideChains.find(residue);
        if (findIt != sideChains.end()){
            auto psidechain = &(findIt->second); // std::unordered_map<std::string, Coords> >sideChains
            float sum_x=0, sum_y=0, sum_z=0;
            float at_x, at_y, at_z;
            float sccount=0.0;
            for (int i=0; i<atomsInResidue; i++) {

                tempAtom = atomType[startAt + i];
                boost::algorithm::trim (tempAtom);

                auto pAtom = psidechain->find(tempAtom);
                if (pAtom != psidechain->end()){
                    Coords *ptempCoord = &(*pAtom).second;
                    at_x = xpos[startAt + i];
                    at_y = ypos[startAt + i];
                    at_z = zpos[startAt + i];

                    sum_x += at_x;
                    sum_y += at_y;
                    sum_z += at_z;

                    tar_x.push_back(at_x); // atoms in residue in structure
                    tar_y.push_back(at_y);
                    tar_z.push_back(at_z);

                    ref_x.push_back(ptempCoord->x); // atoms in residue in library
                    ref_y.push_back(ptempCoord->y);
                    ref_z.push_back(ptempCoord->z);
                    sccount += 1.0;
                } else {
                    std::string msg = " _"+ tempAtom + "_ " + std::to_string(resid);
                    SASTOOLS_UTILS_H::logger("ALIGNMENT ATOM TYPE NOT FOUND", msg);
//                    std::cout << "ATOM TYPES DEFINED FOR " << (*findIt).first << std::endl;
//                    for(auto fit : (*findIt).second){
//                        std::cout << "ATOM => " << " " << fit.first << std::endl;
//                    }
                }
            }

            // create alignment vectors
            float invSV = 1.0f/sccount;
            float aveX = invSV*sum_x;
            float aveY = invSV*sum_y;
            float aveZ = invSV*sum_z;

            // center the residue in the structure
            auto totalAtomsInResidue = (int)sccount;

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

        } else { // create water box and only take waters within 3 Angstroms of structure
            // if residue not found, must come up with an alternative, right now, it is ignored


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

//        new_x = rotation[0]*x_val + rotation[3]*y_val + rotation[6]*z_val + delX;
//        new_y = rotation[1]*x_val + rotation[4]*y_val + rotation[7]*z_val + delY;
//        new_z = rotation[2]*x_val + rotation[5]*y_val + rotation[8]*z_val + delZ;

//        new_x = rotation[0]*x_val + rotation[1]*y_val + rotation[2]*z_val + delX;
//        new_y = rotation[3]*x_val + rotation[4]*y_val + rotation[5]*z_val + delY;
//        new_z = rotation[6]*x_val + rotation[7]*y_val + rotation[8]*z_val + delZ;

        Coords * pT = &tempWaters[i];

        pT->x = rotation[0]*x_val + rotation[1]*y_val + rotation[2]*z_val + delX;
        pT->y = rotation[3]*x_val + rotation[4]*y_val + rotation[5]*z_val + delY;
        pT->z = rotation[6]*x_val + rotation[7]*y_val + rotation[8]*z_val + delZ;
        pT->type = "O";
        pT->occ = pC->occ;
        pT->resid = resid;
        pT->resname = residue;

//        tempWaters[i].x = new_x;
//        tempWaters[i].y = new_y;
//        tempWaters[i].z = new_z;
//        tempWaters[i].type = "O";
//        tempWaters[i].occ = (*pwater)[i].occ;
//        tempWaters[i].resid = resid;
//        tempWaters[i].resname = residue;
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
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", count, " O  ", "HOH", "W", resid.c_str(), pCoord->x, pCoord->y, pCoord->z);
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

        fprintf(pFile, "%-6s%5i %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", count, " O  ", "HOH", "W",
                resid.c_str(),
                pCoord->x + pCenVec->x,
                pCoord->y + pCenVec->y,
                pCoord->z + pCenVec->z);

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
        occupancies[i] = hydration[i].occ;
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


