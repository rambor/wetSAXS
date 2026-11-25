#include <iostream>
#include <boost/program_options.hpp>
#include <sastools/IofQData.h>
#include <filesystem>
#include "AtomisticModel.h"
#include "Waters.h"
#include "Fit.h"
#include "Ensemble.h"
#include "Aligner.h"
#include "BufferModel.h"
#include "IonCloud.h"
#include <memory>

namespace po = boost::program_options;

namespace {
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}

void checkCXRange(float cxlow, float cxhi){
    if ((cxlow < 0 || cxlow > cxhi)) {
        throw "ERROR: lower CX (cxlow) can not be NEGATIVE or greater than cxhi :: " + std::to_string(cxlow) + " < " + std::to_string(cxhi);
    }
}

void checkCX(float cx){
    if ((cx < 0)){
        throw "ERROR: cx can not be NEGATIVE";
    }
}

void checkCORange(float colow, float cohi){
    if ((colow < 0 || colow > cohi)) {
        throw "ERROR: lower C_o (colow) can not be NEGATIVE or greater than cohi :: " + std::to_string(colow) + " < " + std::to_string(cohi);
    }
}

// Comparator function
bool comp(std::string a, std::string b) {

    // first parameter should be greater than
    // second one (for decreasing order)
    return a > b;
}

int main(int argc, const char *const *argv) {

    std::string const FILES_KEY("pdb-files");
    IofQData iofqdata;

    std::cout << "********************                                      *********************\n" << std::endl;
    printf("                               wetSAXS version 1.0 2025 \n\n");
    printf("                           Authors: RP Rambo and FA Reyes\n\n");
    printf("                         Diamond Light Source Ltd, Didcot UK \n");
    printf("                      Lawrence Berkeley National Laboratory, USA \n\n");
    std::cout << "********************                                      *********************" << std::endl;
    std::cout << "    For RNA, PDB file should specify rA, rU, rG and rC to avoid confusion " << std::endl;
    std::cout << "    For DNA, PDB file should specify DA, DT, DG and DC to avoid confusion " << std::endl;
    std::cout << "    ADE and GUA are ambiguous and can specify either DNA or RNA           " << std::endl;
    std::cout << "********************                                      *********************" << std::endl;

    float cxlow = 0.90, cxhi = 1.10, dmax=0.0f, qmin = 0.001, qmax = 0.41;
    float cx = 1.0, colow = 0.5, cohi = 2.0, co = 1.09;

    bool convertToAngstrom = false, forceRNA = false, keepWaters = false, noChiFree = false, singleFit = false,
    writeWaters = false, invacuo=false, mixture=false, ensemble=false, background=false, bufferScaled=false, dmaxFilter = false,
    align_and_map = false;

    int divisions = 101;

    int ns = 237, rounds = 5;

    std::string datFile, atomFile, listOfFiles, bufferFile;
    bool noIonCloud = false;
    IonCloudSettings ionSettings;
    std::unique_ptr<BufferModel> bufferModel;

    // counterion atmosphere for one model: needs the buffer (ionic strength, pH, ion contrasts)
    auto makeIonCloud = [&](AtomisticModel & model) -> std::unique_ptr<IonCloud> {
        if (!bufferModel || noIonCloud || bufferModel->ionicStrength() <= 0.0) return nullptr;
        auto cloud = std::make_unique<IonCloud>(model.getPDBModel(), *bufferModel, ionSettings);
        for (const auto & line : cloud->getReport()) SASTOOLS_UTILS_H::logger("", line);
        model.attachIonCloud(cloud.get());
        if (writeWaters)
            cloud->writePseudoAtoms(model.getPDBModel().getFileStemName() + "_ion_cloud.pdb", model.getPDBModel().getCenteringVector());
        return cloud;
    };

    std::vector<std::string> filenames;

    std::string descText = "\n Usage: wetSAXS myPDBFile.pdb";
    boost::program_options::options_description generic(descText);

    // if no dat file present, simulate using values specified by cw and cx
    generic.add_options()
            ("help,h", "help message")
            ("dat,i", po::value<std::string>(&datFile))
            ("toNM", po::bool_switch(&convertToAngstrom), "Reduce bead radius defining the lattice")
            ("buffer,b", po::value<std::string>(&bufferFile), "Buffer file (components, concentrations, temperature; default 25 C): sets c_x and fixes it")
            ("noIonCloud", po::bool_switch(&noIonCloud), "Do not model the counterion atmosphere (it is modelled whenever --buffer is given)")
            ("ionApproach", po::value<float>(&ionSettings.closestApproach), "Ion cloud: distance of closest approach a, A; pseudo-atoms at a + Debye length (default: 3.0)")
            ("ionWidth", po::value<float>(&ionSettings.width), "Ion cloud: Gaussian width of each pseudo-atom, A (default: 2.5)")
            ("counterFraction", po::value<float>(&ionSettings.nucleicFraction), "Ion cloud: counterion share of the neutralisation of phosphates (default: 0.71)")
            ("contrast,x", po::value<float>(&cx), "Solvent density relative to 0.334 e/A^3, scales the excluded volume (default: 1.0)")
            ("occupancy,o", po::value<float>(&co), "Hydration-site density relative to bulk; shell contrast is co - cx (default: 1.09)")
            (",g", po::bool_switch(&background), "Include constant background in fit (default: false)")
            ("bufferScaled", po::bool_switch(&bufferScaled), "Buffer was scaled by (1 - volume fraction) before subtraction (default: false)")
            ("cxlo", po::value<float>(&cxlow),"lower bound for cx (default: 0.90)")
            ("cxhi", po::value<float>(&cxhi),"upper bound for cx (default: 1.10)")
            ("colo", po::value<float>(&colow),"lower bound for co (default: 0.5)")
            ("cohi", po::value<float>(&cohi),"upper bound for co (default: 2.0)")

            ("divisions,v", po::value<int>(&divisions),"number of divisions between lower and upper bounds for searching (default: 101)")

            ("qmin", po::value<float>(&qmin),"qmin")
            ("qmax", po::value<float>(&qmax),"qmax")
            ("ns,n", po::value<int>(&ns),"total number of points in simulation")
            ("dmax,d", po::value<float>(&dmax),"dmax used to determine the Shannon limit for fitting")
            ("forceRNA", po::bool_switch(&forceRNA), "force A, G, C, U residues to be RNA")
            ("keepWaters", po::bool_switch(&keepWaters), "keep modeled waters in PDB")
            ("invacuo", po::bool_switch(&invacuo), "calculate in vacuo scattering profile")
            ("noChiFree", po::bool_switch(&noChiFree), "Fit all the data without subsampling")
            (",w", po::bool_switch(&writeWaters), "write waters to file")
            ("rounds,r", po::value<int>(&rounds), "Number of rounds for chi-free fitting - odd number")
            (",m", po::bool_switch(&mixture), "mixture model")
            (",e", po::bool_switch(&ensemble), "Ensemble fitting")
            ("listOfFiles,p", po::value<std::string>(&listOfFiles), "Align and create density map from pdb list ")
            (",u", po::bool_switch(&dmaxFilter), "Exclude models that exceed dmax in IofQ file during ensemble fitting")
            ("atomFile,t", po::value<std::string>(&atomFile), "form-factor file")
            ;

//    po::positional_options_description positionalOptions;
//    positionalOptions.add("pdb-file", 1);
    po::options_description hidden("Hidden options");
    hidden.add_options()
            (FILES_KEY.c_str(), po::value(&filenames)->required(), "list of PDB files");

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);

    // And this one to display help
    po::options_description visible_options;
    visible_options.add(generic);

    po::positional_options_description pos;
    pos.add(FILES_KEY.c_str(), -1);

    po::variables_map vm;

    try {
        // Only parse the options, so we can catch the explicit `--files`
        auto parsed = po::command_line_parser(argc, argv)
                .options(cmdline_options)
                .positional(pos)
                .run();

       // po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
        // Make sure there were no non-positional `files` options
        for (auto const& opt : parsed.options) {
            if ((opt.position_key == -1) && (opt.string_key == FILES_KEY)) {
                throw po::unknown_option(FILES_KEY);
            }
        }

        po::store(parsed, vm);
        // checking options list
        if (vm.count("help") || !vm.count("pdb-files")) {
            std::cout << generic << "\n";
            return SUCCESS;
        }

        po::notify(vm);

        // buffer composition fixes c_x: only c_o (hydration excess c_o - c_x) is searched
        if (vm.count("buffer")){
            try {
                bufferModel = std::make_unique<BufferModel>(bufferFile);
                const BufferModel & buffer = *bufferModel;
                for (const auto & line : buffer.getReport()) std::cout << "  " << line << std::endl;
                if (vm.count("contrast")){
                    std::cerr << "ERROR: give either -x (contrast) or --buffer, not both" << std::endl;
                    return ERROR_IN_COMMAND_LINE;
                }
                if (vm.count("cxlo") || vm.count("cxhi")){
                    std::cout << "  WARNING: --cxlo/--cxhi ignored, c_x fixed from buffer file" << std::endl;
                }
                cx = (float)buffer.getCx();
                cxlow = cx;
                cxhi = cx;
            } catch (const std::exception & e){
                std::cerr << "ERROR: " << e.what() << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }
        }


         if (vm["toNM"].as<bool>()){
             convertToAngstrom = true;
         }

         // performs median-based SVD alignment using input file and suggested residue list
         if (vm.count("listOfFiles")){ // perform alignment
             //getcwd = boost::filesystem::path full_path(boost::filesystem::current_path());
             Aligner alignit(listOfFiles, boost::filesystem::current_path().string());
             alignit.align();
             std::cout << " Align file " << std::endl;
         }

         if (vm.count("occupancy") || vm.count("contrast")){

             if (co < 0){
                 std::cerr << "ERROR: Occupancy can not be NEGATIVE :: " << co << std::endl;
                 return ERROR_IN_COMMAND_LINE;
             }

             if (cx < 0){
                 std::cerr << "ERROR: Contrast Term (CX) can not be NEGATIVE :: " << cx << std::endl;
                 return ERROR_IN_COMMAND_LINE;
             }

             std::cout << "Occupancy " << co << std::endl;
             std::cout << "Contrast " << cx << std::endl;
         }
//
//
//
//         if (datFile.empty()){ // simulate using specified cx anc bf values
//
//             std::cout << "Dat file empty simulating in vacuo" << std::endl;
//             std::vector<float> qvalues;
//             float delta = (qmax-qmin)/(float)ns;
//             for(int i=0; i<ns; i++){
//                 qvalues.push_back((float)(i+1)*delta);
//             }
//
//             for(const std::string& filename : filenames){
//
//                 SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);
//                 AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
//
//                 // if dmax is not specified, must retrieve from PDB
//                 // dmax should really be based on the empirical value not model taken from IofQfile
//                 if (!vm.count("dmax")){
//                     dmax = model.getDmax();
//                 }
//
//                 SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));
//
//                 unsigned lmax = 2*(unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
//                 SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
//                 SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");
//                 model.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);
//
//                 if (invacuo) {
//
//                 } else { // simulate with waters
//                     Waters waters = Waters();
//
//                     if (keepWaters){
//
//                         waters.extractWatersFromAtomisticModel(model);
//
//                     } else {
//
//                         waters.hydrateAtomisticModel(model); // add waters to PDB
//                         waters.createSphericalCoordinateOfHydration();
//
//                         waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);
//                         std::string name = model.getPDBModel().getFileStemName();
//
//                         // translate water model back to original input coordinates
//                         waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
//                     }
//
//                     // create water model
//                     Fit simulated = Fit(lmax);
//                     simulated.simulate(cx, co, qvalues, model, waters);
//
//                 }
//             }
//
//             return 0;
//
//         }
//
//         if (!ensemble) { // fit data
//
//             iofqdata = IofQData(datFile, convertToAngstrom);
//             iofqdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format
//
//             if (iofqdata.getDmax() < 1 && !vm.count("dmax")) {
//                 std::cerr << "ERROR: DMAX of dataset must be specified, not found in IofQ file" << std::endl;
//                 return ERROR_IN_COMMAND_LINE;
//             }
//
//             if (iofqdata.getShannonBins() < 1) { // dmax is specified
//                 iofqdata.setDmax(dmax); // dmax is specific to the experimental data, ideally, no model should exceed this....
//             } else {
//                 dmax = iofqdata.getDmax();
//             }
//
//             // need to use qmax if specifz0ied otherwise stick to IofQ qmax
//             if (vm.count("qmax") > 0){
//                 qmax = vm["qmax"].as<float>() < iofqdata.getQmax() ? vm["qmax"].as<float>() : iofqdata.getQmax();
//             }
//
//             // truncate q-values of iofqdata to specified qmax
//             if (qmax < iofqdata.getQmax()){
//                 iofqdata.truncateToQmax(qmax);
//             }
//
//             unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
//             SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
//
//             // check that the cx and co are sensible
//             try {
//                 checkCORange(colow, cohi);
//                 checkCXRange(cxlow, cxhi);
//                 checkCX(cx);
//             } catch (std::string & fl) {
//                 std::cerr << fl << std::endl;
//                 return ERROR_IN_COMMAND_LINE;
//             }
//
//             std::vector<std::string> logfileOutput;
//
//             // chi free fitting is default
//             if (!noChiFree && !vm.count("contrast") && !vm.count("occupancy")){
//
//                 Fit fitting = Fit(lmax, divisions, cxlow, cxhi, colow, cohi);
//
//                 for(const std::string& filename : filenames){
//
//                     try{
//                         AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
//
//                         if (model.getPDBModel().getTotalResidues() < 1){
//                             throw (filename);
//                         }
//
//                         // if atomFile present load the atom descriptors
//                         if (vm.count("atomFile")){
//                             if (!model.getPDBModel().validateATOMTYPESFileFormat(atomFile)){
//                                 std::cerr << "ERROR: check input --atomFile : "  << atomFile << std::endl;
//                                 exit(0);
//                             }
//
//                             model.getPDBModel().updateAtomDescriptions(atomFile);
//                             exit(0);
//                         }
//
//                         SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
//                         SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");
//
//                         Waters waters = Waters();
//                         waters.hydrateAtomisticModel(model); // add waters to PDB
//                         waters.createSphericalCoordinateOfHydration();
//
//                         fitting.chiFreeSearch(rounds, iofqdata, model, waters);
//
//                         if (writeWaters){
//                             std::string name = model.getPDBModel().getFileStemName();
//                             model.getPDBModel().getCenteringVector();
//                             waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
//                         }
//
//                         logfileOutput.emplace_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));
//
//                     } catch (std::string & fl){
//                         SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
//                     }
//                 }
//
//             } else {
//
//                 iofqdata.setAllDataToWorkingSet();
//
//
//
//                 // no chi free, use same working set for all PDB files in fit
// //                iofqdata.makeWorkingSet(3);
//                 std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();
//                 // can fit with specified cx and b-factor or do a search
//                 // if fitting multiple PDB files to dataset with no chiFree, use same working set
//
//                 for(const std::string& filename : filenames){
//
//                     SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);
//
//                     // dmax should really be based on the empirical value not model taken from IofQfile
//
//                     try{
//
//                         AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
//
//                         if (model.getPDBModel().getTotalResidues() < 1){
//                             throw (filename);
//                         }
//
//                         SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
//                         SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");
//
//                         model.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);
//
//                         // create water model
//                         Waters waters = Waters();
//                         waters.hydrateAtomisticModel(model); // add waters to PDB
//                         waters.createSphericalCoordinateOfHydration();
//                         waters.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);
//
//                         if (writeWaters){
//                             std::string name = model.getPDBModel().getFileStemName();
//                             model.getPDBModel().getCenteringVector();
//                             waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
//                         }
//
//                         if (!vm.count("contrast")){ // no Chi-free fitting single working set only
//
//                             Fit fitting = Fit(lmax, divisions, cxlow, cxhi, colow, cohi);
//                             SASTOOLS_UTILS_H::logger("", "Grid Search of B-factor and Cx");
//                             fitting.search(iofqdata, model, waters);
//                             // append a summary file with the best fitted parameters so a user could
//                             logfileOutput.push_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));
//
//                         } else {
//
//                             Fit fitting = Fit(lmax);
//
//                             if (vm.count("contrast")) { // fit using fixed contrast and occupancy
//
//                                 // working set was set above to 3 times Shannon Number, but for fixed fitting should do all?
//                                 iofqdata.setAllDataToWorkingSet();
//                                 std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();
//                                 // recalculate Partials for the initial models
//                                 waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);
//                                 model.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues, true);
//
//                                 fitting.fitFixedCoAndCx(iofqdata, model, waters, co, cx);
//                                 // (qvalues, iofqdata.getFilename(), model);
//                                 fitting.writeBestModelToFile(qvalues, iofqdata.getFilename(), model);
//
//                             } else if (vm.count("occupancy")){ // search contrast with fixed bfactor
//
//                                 fitting.searchCX(iofqdata, model, waters, co, cxlow, cxhi, divisions);
//
//                             } else if (vm.count("contrast")){ // search bfactor with fixed contrast
//
//                                 //fitting.searchBFactor(iofqdata, model, waters, cx, upperB, divisions);
//                             }
//
//                             logfileOutput.emplace_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));
//                         }
//
//                     } catch (std::string & fl){
//                         SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
//                     }
//                 }
//             }
//
//             // write logfile
//             FILE * pFile = fopen("fit_summary.txt", "w");
//             fprintf(pFile, "  SCORE    CHI2   DW    CX  BFactor NAME\n");
//             for(auto & line : logfileOutput){
//                 fprintf(pFile, line.c_str());
//             }
//             fclose(pFile);
//
//
//         } else if (ensemble){
//
//             Ensemble ens = Ensemble(filenames, std::filesystem::current_path().string());
//
//             iofqdata = IofQData(datFile, convertToAngstrom);
//             iofqdata.extractData();
//
//             if (iofqdata.getDmax() < 1 && !vm.count("dmax")) {
//                 std::cerr << "ERROR: DMAX of dataset must be specified, not found in IofQ file" << std::endl;
//                 return ERROR_IN_COMMAND_LINE;
//             } else if (iofqdata.getShannonBins() < 1) { // dmax is specified
//                 iofqdata.setDmax(dmax); // dmax is specific to the experimental data, ideally, no model should exceed this....
//             } else {
//                 dmax = iofqdata.getDmax();
//             }
//
//             // need to use qmax if specified otherwise stick to IofQ qmax
//             if (vm.count("qmax") > 0){
//                 qmax = vm["qmax"].as<float>() < iofqdata.getQmax() ? vm["qmax"].as<float>() : iofqdata.getQmax();
//             }
//
//             // truncate q-values of iofqdata to specified qmax
//             if (qmax < iofqdata.getQmax()){
//                 iofqdata.truncateToQmax(qmax);
//             }
//
//             ens.ce_search(iofqdata, cx, co, dmaxFilter);
//             //ens.combinatorial_search(iofqdata, cx, bfactor);
//
//         }

    // } catch (boost::program_options::required_option& e){
    //
    //     //if (!vm.count("help")) {
    //         std::cerr << "ERROR: " << e.what() << std::endl;
    //     std::cerr << generic << std::endl;
    //         return ERROR_IN_COMMAND_LINE;
    //     //}

    } catch (boost::program_options::error& e){

        std::cerr << "ERROR:: " << e.what() << std::endl;
        return ERROR_IN_COMMAND_LINE;

    } catch (const std::invalid_argument& e){

        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr<<"Type "<<typeid(e).name()<<std::endl;

    }


    if (datFile.empty()){ // simulate using specified cx anc bf values

            std::cout << "Dat file empty simulating in vacuo" << std::endl;
            std::vector<float> qvalues;
            float delta = (qmax-qmin)/(float)ns;
            for(int i=0; i<ns; i++){
                qvalues.push_back((float)(i+1)*delta);
            }

            for(const std::string& filename : filenames){

                SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);
                AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
                std::unique_ptr<IonCloud> ionCloud = makeIonCloud(model);

                // if dmax is not specified, must retrieve from PDB
                // dmax should really be based on the empirical value not model taken from IofQfile
                if (!vm.count("dmax")){
                    dmax = model.getDmax();
                }

                SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));

                unsigned lmax = 2*(unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
                SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));
                SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");
                model.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);

                if (invacuo) {

                } else { // simulate with waters
                    Waters waters = Waters();

                    if (keepWaters){

                        waters.extractWatersFromAtomisticModel(model);

                    } else {

                        waters.hydrateAtomisticModel(model); // add waters to PDB
                        waters.createSphericalCoordinateOfHydration();

                        waters.calculatePartialAmplitudes(lmax, qvalues.size(), qvalues);
                        std::string name = model.getPDBModel().getFileStemName();

                        // translate water model back to original input coordinates
                        waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
                    }

                    // create water model
                    Fit simulated = Fit(lmax);
                    simulated.simulate(cx, co, qvalues, model, waters);

                }
            }

            return 0;

        }

        if (!ensemble) { // fit data

            iofqdata = IofQData(datFile, convertToAngstrom);
            iofqdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format

            if (iofqdata.getDmax() < 1 && !vm.count("dmax")) {
                std::cerr << "ERROR: DMAX of dataset must be specified, not found in IofQ file" << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }

            if (iofqdata.getShannonBins() < 1) { // dmax is specified
                iofqdata.setDmax(dmax); // dmax is specific to the experimental data, ideally, no model should exceed this....
            } else {
                dmax = iofqdata.getDmax();
            }

            // need to use qmax if specifz0ied otherwise stick to IofQ qmax
            if (vm.count("qmax") > 0){
                qmax = vm["qmax"].as<float>() < iofqdata.getQmax() ? vm["qmax"].as<float>() : iofqdata.getQmax();
            }

            // truncate q-values of iofqdata to specified qmax
            if (qmax < iofqdata.getQmax()){
                iofqdata.truncateToQmax(qmax);
            }

            unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
            SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));

            // check that the cx and co are sensible
            try {
                checkCORange(colow, cohi);
                checkCXRange(cxlow, cxhi);
                checkCX(cx);
            } catch (std::string & fl) {
                std::cerr << fl << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }

            std::vector<std::string> logfileOutput;

            // chi free fitting is default
            if (!noChiFree && !vm.count("contrast") && !vm.count("occupancy")){

                Fit fitting = Fit(lmax, divisions, cxlow, cxhi, colow, cohi);

                for(const std::string& filename : filenames){

                    try{
                        AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
                        std::unique_ptr<IonCloud> ionCloud = makeIonCloud(model);

                        if (model.getPDBModel().getTotalResidues() < 1){
                            throw (filename);
                        }

                        // if atomFile present load the atom descriptors
                        if (vm.count("atomFile")){
                            if (!model.getPDBModel().validateATOMTYPESFileFormat(atomFile)){
                                std::cerr << "ERROR: check input --atomFile : "  << atomFile << std::endl;
                                exit(0);
                            }

                            model.getPDBModel().updateAtomDescriptions(atomFile);
                            exit(0);
                        }

                        SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
                        SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");

                        Waters waters = Waters();
                        waters.hydrateAtomisticModel(model); // add waters to PDB
                        waters.createSphericalCoordinateOfHydration();

                        // l needed is ~ q*r for the outermost atom (r measured from the centre), not the Shannon number
                        const unsigned int modelLmax = (unsigned int)std::ceil(qmax*model.getPDBModel().getSMax()) + 7;
                        SASTOOLS_UTILS_H::logger("LMAX (model)", std::to_string(modelLmax));
                        fitting.setLmax(modelLmax);

                        // solvent fluctuations of the bulk water displaced by protein + hydration sites
                        SolventFluctuation solvent;
                        solvent.buildExcludedRegion(model, waters);
                        solvent.computeSurfaceDeficit();
                        SASTOOLS_UTILS_H::logger("Excluded volume (SES, A^3)", std::to_string((int)solvent.getExcludedVolume()));
                        fitting.setSolventFluctuation(&solvent, bufferScaled);

                        fitting.chiFreeSearch(rounds, iofqdata, model, waters);
                        fitting.setSolventFluctuation(nullptr, bufferScaled); // solvent goes out of scope

                        if (writeWaters){
                            std::string name = model.getPDBModel().getFileStemName();
                            model.getPDBModel().getCenteringVector();
                            waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
                        }

                        logfileOutput.emplace_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));

                    } catch (std::string & fl){
                        SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
                    }
                }

            } else {

                iofqdata.setAllDataToWorkingSet();

                // no chi free, use same working set for all PDB files in fit
//                iofqdata.makeWorkingSet(3);
                std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();
                // can fit with specified cx and b-factor or do a search

                for(const std::string& filename : filenames){

                    SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);

                    // dmax should really be based on the empirical value not model taken from IofQfile

                    try{

                        AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);
                        std::unique_ptr<IonCloud> ionCloud = makeIonCloud(model);

                        if (model.getPDBModel().getTotalResidues() < 1){
                            throw (filename);
                        }

                        SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
                        SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");

                        // l needed is ~ q*r for the outermost atom (r measured from the centre), not the Shannon number
                        const unsigned int modelLmax = (unsigned int)std::ceil(qmax*model.getPDBModel().getSMax()) + 7;
                        SASTOOLS_UTILS_H::logger("LMAX (model)", std::to_string(modelLmax));

                        model.calculatePartialAmplitudes(modelLmax, iofqdata.getTotalInWorkingSet(), qvalues);

                        // create water model
                        Waters waters = Waters();
                        waters.hydrateAtomisticModel(model); // add waters to PDB
                        waters.createSphericalCoordinateOfHydration();
                        waters.calculatePartialAmplitudes(modelLmax, iofqdata.getTotalInWorkingSet(), qvalues);

                        SolventFluctuation solvent;
                        solvent.buildExcludedRegion(model, waters);
                        solvent.computeSurfaceDeficit();
                        SASTOOLS_UTILS_H::logger("Excluded volume (SES, A^3)", std::to_string((int)solvent.getExcludedVolume()));

                        if (writeWaters){
                            std::string name = model.getPDBModel().getFileStemName();
                            model.getPDBModel().getCenteringVector();
                            waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
                        }

                        if (!vm.count("contrast")){ // no Chi-free fitting single working set only

                            Fit fitting = Fit(modelLmax, divisions, cxlow, cxhi, colow, cohi);
                            fitting.setSolventFluctuation(&solvent, bufferScaled);
                            SASTOOLS_UTILS_H::logger("", "Grid Search of B-factor and Cx");
                            fitting.search(iofqdata, model, waters);
                            // append a summary file with the best fitted parameters so a user could
                            logfileOutput.push_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));

                        } else {

                            Fit fitting = Fit(modelLmax);
                            fitting.setSolventFluctuation(&solvent, bufferScaled);

                            if (vm.count("contrast")) { // fit using fixed contrast and occupancy

                                // working set was set above to 3 times Shannon Number, but for fixed fitting should do all?
                                iofqdata.setAllDataToWorkingSet();
                                std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();
                                // recalculate Partials for the initial models
                                waters.calculatePartialAmplitudes(modelLmax, qvalues.size(), qvalues, true);
                                model.calculatePartialAmplitudes(modelLmax, qvalues.size(), qvalues, true);

                                fitting.fitFixedCoAndCx(iofqdata, model, waters, co, cx);
                                // (qvalues, iofqdata.getFilename(), model);
                                fitting.writeBestModelToFile(qvalues, iofqdata.getFilename(), model);

                            } else if (vm.count("occupancy")){ // search contrast with fixed bfactor

                                fitting.searchCX(iofqdata, model, waters, co, cxlow, cxhi, divisions);

                            } else if (vm.count("contrast")){ // search bfactor with fixed contrast

                                //fitting.searchBFactor(iofqdata, model, waters, cx, upperB, divisions);
                            }

                            logfileOutput.emplace_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));
                        }

                    } catch (std::string & fl){
                        SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
                    }
                }
            }

            // write logfile
            FILE * pFile = fopen("fit_summary.txt", "w");
            fprintf(pFile, "  SCORE    CHI2   DW    CX  BFactor NAME\n");
            for(auto & line : logfileOutput){
                fprintf(pFile, line.c_str());
            }
            fclose(pFile);


        } else if (ensemble){

            Ensemble ens = Ensemble(filenames, std::filesystem::current_path().string());

            iofqdata = IofQData(datFile, convertToAngstrom);
            iofqdata.extractData();

            if (iofqdata.getDmax() < 1 && !vm.count("dmax")) {
                std::cerr << "ERROR: DMAX of dataset must be specified, not found in IofQ file" << std::endl;
                return ERROR_IN_COMMAND_LINE;
            } else if (iofqdata.getShannonBins() < 1) { // dmax is specified
                iofqdata.setDmax(dmax); // dmax is specific to the experimental data, ideally, no model should exceed this....
            } else {
                dmax = iofqdata.getDmax();
            }

            // need to use qmax if specified otherwise stick to IofQ qmax
            if (vm.count("qmax") > 0){
                qmax = vm["qmax"].as<float>() < iofqdata.getQmax() ? vm["qmax"].as<float>() : iofqdata.getQmax();
            }

            // truncate q-values of iofqdata to specified qmax
            if (qmax < iofqdata.getQmax()){
                iofqdata.truncateToQmax(qmax);
            }

            ens.ce_search(iofqdata, cx, co, dmaxFilter);
            //ens.combinatorial_search(iofqdata, cx, bfactor);

        }

    return 0;



}


