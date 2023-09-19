#include <iostream>
#include <boost/program_options.hpp>
#include <sastools/IofQData.h>
#include "AtomisticModel.h"
#include "Waters.h"
#include "Fit.h"

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

void checkUpperB(float upperB){
    if ((upperB < 0 || upperB > 500)){
        throw "ERROR: Upper B-Factor can not be NEGATIVE or greater than 500:: " + std::to_string(upperB);
    }
}

int main(int argc, const char *const *argv) {

    std::string const FILES_KEY("pdb-files");
    IofQData iofqdata;

    std::cout << "********************                                      *********************\n" << std::endl;
    printf("                               wetSAXS version 1.0 2022 \n\n");
    printf("                           Authors: RP Rambo and FA Reyes\n\n");
    printf("                         Diamond Light Source Ltd, Didcot UK \n");
    printf("                      Lawrence Berkeley National Laboratory, USA \n\n");
    std::cout << "********************                                      *********************" << std::endl;
    std::cout << "    For RNA, PDB file should specify rA, rU, rG and rC to avoid confusion " << std::endl;
    std::cout << "    For DNA, PDB file should specify DA, DT, DG and DC to avoid confusion " << std::endl;
    std::cout << "    ADE and GUA are ambiguous and can specify either DNA or RNA           " << std::endl;
    std::cout << "********************                                      *********************" << std::endl;

    float cxlow = 0.8, cxhi = 1.5, dmax=0.0f, qmin = 0.001, qmax = 0.41;

    float bfactor = 0.98, cx = 1.108, upperB=87;

    bool convertToAngstrom = false, forceRNA = false, keepWaters = false, noChiFree = false, singleFit = false,
    writeWaters = false, invacuo=false;

    int divisions = 101;

    int ns = 237, rounds = 5;

    std::string datFile;
    std::vector<std::string> file_names;

    boost::program_options::options_description generic("\n Usage: wetSAXS myPDBFile.pdb");

    // if no dat file present, simulate using values specified by cw and cx
    generic.add_options()
            ("help,h", "help message")
            ("dat,i", po::value<std::string>(&datFile))
            ("toNM", po::bool_switch(&convertToAngstrom), "Reduce bead radius defining the lattice")
            ("bfactor,b", po::value<float>(&bfactor), "Use to specify fixed B-factor for simulation or fitting (default: 0.98)")
            ("contrast,x", po::value<float>(&cx), "Excluded volume contrast adjustment (default: 1.108)")
            ("cxlo", po::value<float>(&cxlow),"lower bound for excluded volume adjustment")
            ("cxhi", po::value<float>(&cxhi),"upper bound for excluded volume adjustment")
            ("upperB", po::value<float>(&upperB),"upper bound for B-factor during search")
            ("divisions,v", po::value<int>(&divisions),"number of divisions between lower and upper bounds for searching (default: 101)")

            ("qmin", po::value<float>(&qmin),"qmin")
            ("qmax", po::value<float>(&qmax),"qmax")
            ("ns,n", po::value<int>(&ns),"total number of points in simulation")
            ("dmax,d", po::value<float>(&dmax),"dmax used to determine the Shannon limit for fitting")
            ("forceRNA", po::bool_switch(&forceRNA), "force A, G, C, U residues to be RNA")
            ("keepWaters", po::bool_switch(&keepWaters), "keep modeled waters in PDB")
            ("invacuo", po::bool_switch(&invacuo), "calculate in vacuo scattering profile")
            ("noChiFree", po::bool_switch(&noChiFree), "Fit all the data without subsampling")
            (",w", po::bool_switch(&writeWaters), "Fit all the data without subsampling")
            ("rounds,r", po::value<int>(&rounds), "Number of rounds for chi-free fitting - odd number")
            ;

//    po::positional_options_description positionalOptions;
//    positionalOptions.add("pdb-file", 1);

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
            (FILES_KEY.c_str(), po::value(&file_names)->required(), "list of PDB files");

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
        po::notify(vm);

        if (vm["toNM"].as<bool>()){
            convertToAngstrom = true;
        }

        if (vm.count("bfactor") || vm.count("contrast")){

            if (bfactor < 0){
                std::cerr << "ERROR: B-Factor can not be NEGATIVE :: " << bfactor << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }

            if (cx < 0){
                std::cerr << "ERROR: Contrast Term (CX) can not be NEGATIVE :: " << cx << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }

            std::cout << "B-Factor " << bfactor << std::endl;
            std::cout << "Contrast " << cx << std::endl;


        }

        if (datFile.empty()){ // simulate using specified cx anc bf values

            std::cout << "Dat file empty simulating in vacuo" << std::endl;
            std::vector<float> qvalues;
            float delta = (qmax-qmin)/(float)ns;
            for(int i=0; i<ns; i++){
                qvalues.push_back((float)(i+1)*delta);
            }


            for(const std::string& filename : file_names){

                SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);
                AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);

                // if dmax is not specified, must retrieve from PDB
                // dmax should really be based on the empirical value not model taken from IofQfile
                if (!vm.count("dmax")){
                    dmax = model.getDmax();
                }

                SASTOOLS_UTILS_H::logger("DMAX SET", std::to_string(dmax));

                unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
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
                    simulated.simulate(cx, bfactor, qvalues, model, waters);

                }
            }

            return 0;

        } else { // fit

            iofqdata = IofQData(datFile, convertToAngstrom);
            iofqdata.extractData(); // assume dmax is not set as data from other facilities will not be compatible with Scatter format

            if (iofqdata.getDmax() < 1 && !vm.count("dmax")) {
                std::cerr << "ERROR: DMAX of dataset must be specified, not found in IofQ file" << std::endl;
                return ERROR_IN_COMMAND_LINE;
            } else if (iofqdata.getShannonBins() < 1) { // dmax is specified
                iofqdata.setDmax(dmax); // dmax is specific to the experimental data, ideally, no model should exceed this....
            } else {
                dmax = iofqdata.getDmax();
            }

            unsigned lmax = (unsigned int)((qmax*dmax)/SIMDE_MATH_PI) + 1;
            SASTOOLS_UTILS_H::logger("LMAX", std::to_string(lmax));


            try {
                checkUpperB(upperB);
                checkCXRange(cxlow, cxhi);
                checkCX(cx);
            } catch (std::string & fl) {
                std::cout << fl << std::endl;
                return ERROR_IN_COMMAND_LINE;
            }



            std::vector<std::string> logfileOutput;
            // chi free fitting is default
            if (!noChiFree && !vm.count("upperB") && !vm.count("contrast")){

                Fit fitting = Fit(lmax, divisions,cxlow,cxhi, upperB);

                for(const std::string& filename : file_names){

                    try{
                        AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);

                        if (model.getPDBModel().getTotalResidues() < 1){
                            throw (filename);
                        }

                        SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
                        SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");

                        Waters waters = Waters();
                        waters.hydrateAtomisticModel(model); // add waters to PDB
                        waters.createSphericalCoordinateOfHydration();

                        if (writeWaters){
                            std::string name = model.getPDBModel().getFileStemName();
                            model.getPDBModel().getCenteringVector();
                            waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
                        }

                        fitting.chiFreeSearch(rounds, iofqdata, model, waters);

                        logfileOutput.emplace_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));

                    } catch (std::string & fl){
                        SASTOOLS_UTILS_H::logger("Improper PDB file", "skipping..." + fl);
                    }
                }

            } else {
                // no chi free, use same working set for all PDB files in fit
                iofqdata.makeWorkingSet(3);
                std::vector<float> qvalues = iofqdata.getWorkingSetQvalues();
                // can fit with specified cx and b-factor or do a search
                // if fitting multiple PDB files to dataset with no chiFree, use same working set

                for(const std::string& filename : file_names){

                    SASTOOLS_UTILS_H::logger("FITTING MODEL FILE",filename);

                    // dmax should really be based on the empirical value not model taken from IofQfile

                    try{

                        AtomisticModel model = AtomisticModel(filename, forceRNA, keepWaters);

                        if (model.getPDBModel().getTotalResidues() < 1){
                            throw (filename);
                        }

                        SASTOOLS_UTILS_H::logger("Model DMAX", std::to_string(model.getDmax()));
                        SASTOOLS_UTILS_H::logger("", "Calculating Amplitudes MODEL");

                        model.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);

                        // create water model
                        Waters waters = Waters();
                        waters.hydrateAtomisticModel(model); // add waters to PDB
                        waters.createSphericalCoordinateOfHydration();
                        waters.calculatePartialAmplitudes(lmax, iofqdata.getTotalInWorkingSet(), qvalues);

                        if (writeWaters){
                            std::string name = model.getPDBModel().getFileStemName();
                            model.getPDBModel().getCenteringVector();
                            waters.translateAndWriteWatersToFile(name + "_hydration_model", model.getPDBModel().getCenteringVector());
                        }

                        if (!vm.count("upperB") && !vm.count("contrast")){ // no Chi-free fitting single working set only

                            Fit fitting = Fit(lmax, divisions,cxlow,cxhi, upperB);
                            SASTOOLS_UTILS_H::logger("", "Grid Search of B-factor and Cx");
                            fitting.search(iofqdata, model, waters);
                            // append a summary file with the best fitted parameters so a user could
                            logfileOutput.push_back(fitting.getScoreText(model.getPDBModel().getFileStemName()));

                        } else {

                            Fit fitting = Fit(lmax);

                            if (vm.count("bfactor") && vm.count("contrast")) { // fit using fixed bfactor and contrast

                                fitting.fitFixedCxandBFactor(iofqdata, model, waters, bfactor, cx);

                            } else if (vm.count("bfactor")){ // search contrast with fixed bfactor

                                fitting.searchCX(iofqdata, model, waters, bfactor, cxlow, cxhi, divisions);

                            } else if (vm.count("contrast")){ // search bfactor with fixed contrast

                                fitting.searchBFactor(iofqdata, model, waters, cx, upperB, divisions);
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

        }

    } catch (boost::program_options::required_option& e){

        if (!vm.count("help")) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return ERROR_IN_COMMAND_LINE;
        }

    } catch (boost::program_options::error& e){
        std::cerr << "ERROR: " << e.what() << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (const std::invalid_argument& e){
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr<<"Type "<<typeid(e).name()<<std::endl;
    }

    if (vm.count("help") || !vm.count("pdb-files")) {
        std::cout << generic << "\n";
        return SUCCESS;
    }


    return 0;



}


