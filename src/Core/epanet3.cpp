/* EPANET 3.1
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

///////////////////////////////////////////////
// Implementation of EPANET 3.1's API library  //
///////////////'''''///////////////////////////

// TO DO:
// - finish implementing all of the functions declared in EPANET3.H
// - provide a brief comment on what each function does

#include "epanet3.h"
#include "Core/Project.h"
#include "Core/datamanager.h"
#include "Core/constants.h"
#include "Core/error.h"
#include "Core/network.h"
#include "Core/hydengine.h"
#include "Core/hydbalance.h"
#include "Elements/valve.h"
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Elements/pumpcurve.h"
#include "Elements/control.h"
#include "Elements/junction.h"
#include "Elements/tank.h"
#include "Elements/link.h"
#include "Utilities/utilities.h"
#include "rwcggasolver.h"
#include "matrixsolver.h"
#include "linkparser.h"
#include "options.h"
#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>   //for debugging
#include <algorithm>
#include <time.h>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>


using namespace Epanet;
using namespace std;

#define project(p) ((Project *)p)

extern "C" {

//-----------------------------------------------------------------------------

int EN_getVersion(int* version)
{
    *version = VERSION;
    return 0;
}

namespace {

    struct ResultPreset {
        const char* key;     // .inp adinda aranan kelime (kucuk harf)
        const char* nodes;   // virgulle ayrilmis dugum kimlikleri
    };

    // Sira onemli: dugum varligina gore secimde daha ozel aglar once denenir.
    const ResultPreset RESULT_PRESETS[] = {
        { "hk-burst",  "612,3234" },
        { "nault2018", "J71,J70,J874" },
        { "nault2016", "J36,J35,J16,J24" },
        { "tnet1",     "N2,N3,N4,N5,N6,N7,N8" },
        { "onizuka",   "J1,J2,J3,J4,J5,J6,T1" },
        { "rpv",       "2,3" },
    };

    std::vector<std::string> splitCsv(const std::string& s) {
        std::vector<std::string> out;
        std::stringstream ss(s);
        std::string item;
        while (std::getline(ss, item, ',')) {
            std::size_t b = item.find_first_not_of(" \t\r\n");
            std::size_t e = item.find_last_not_of(" \t\r\n");
            if (b != std::string::npos) out.push_back(item.substr(b, e - b + 1));
        }
        return out;
    }

    std::string toLower(std::string s) {
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }

    // EN_getNodeIndex char* bekleyebilecegi icin degistirilebilir bir kopya verilir.
    bool findNode(Epanet::Project& p, const std::string& id, int& index) {
        std::string name = id;
        index = -1;
        return !name.empty() && EN_getNodeIndex(&name[0], &index, &p) == 0;
    }

    std::vector<std::string> chooseResultNodes(Epanet::Project& p, const std::string& inpStem,
        std::string& source) {
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
        const char* env = std::getenv("CGGA_RESULT_NODES");
        if (env && *env) {
            source = "CGGA_RESULT_NODES";
            return splitCsv(env);
        }
        const std::string stem = toLower(inpStem);
        for (const ResultPreset& r : RESULT_PRESETS) {
            if (stem.find(r.key) != std::string::npos) {
                source = std::string(r.key) + ", dosya adindan";
                return splitCsv(r.nodes);
            }
        }
        for (const ResultPreset& r : RESULT_PRESETS) {
            std::vector<std::string> names = splitCsv(r.nodes);
            bool all = true;
            int idx;
            for (const std::string& n : names) {
                if (!findNode(p, n, idx)) { all = false; break; }
            }
            if (all) {
                source = std::string(r.key) + ", dugumlerden";
                return names;
            }
        }
        source = "taninmadi";
        return std::vector<std::string>();
    }

} // namespace

//-----------------------------------------------------------------------------

int EN_runEpanet(const char* inpFile, const char* rptFile, const char* outFile) {
    std::cout << "\n... EPANET Version 3.0 with CGGA Adaptive Scheme\n";

    Epanet::Project p;
    int err = 0;

    // ------------------------------------------------------------------
    // Sonuc dosyasi: .inp ile ayni klasorde, "Results-<inp adi>.txt"
    // ------------------------------------------------------------------
    const std::string inpPath(inpFile);
    const std::size_t slash = inpPath.find_last_of("\\/");
    const std::string folder = (slash == std::string::npos) ? "" : inpPath.substr(0, slash + 1);
    std::string stem = (slash == std::string::npos) ? inpPath : inpPath.substr(slash + 1);
    const std::size_t dot = stem.find_last_of('.');
    if (dot != std::string::npos) stem = stem.substr(0, dot);
    const std::string resultsPath = folder + "Results-" + stem + ".txt";

    std::ofstream hadilan(resultsPath);
    if (!hadilan.is_open()) {
        std::cerr << "Error: Could not open " << resultsPath << " for writing\n";
    }

    // Zaman bicimi: varsayilan s:dd:ss.sss; CGGA_TIME_FORMAT=s ise saniye.
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
    const char* timeEnv = std::getenv("CGGA_TIME_FORMAT");
    const bool timeInSeconds = timeEnv && toLower(timeEnv) == "s";

    std::vector<std::string> outNames;   // sonuc dosyasina yazilan dugumler
    std::vector<int> outIndex;

    clock_t start_t = clock();

    double T_end = p.getTotalSimulationTime(); // Total simulation time

    for (;;)
    {
        // ... open the command line files and load network data
        if ((err = p.openReport(rptFile))) break;
        std::cout << "\n    Reading input file ...";
        if ((err = p.load(inpFile))) break;
        if ((err = p.openOutput(outFile))) break;
        p.writeSummary();

        // ... yazilacak dugumleri sec, indekslerini bir kez bul, basligi yaz
        std::string source;
        for (const std::string& n : chooseResultNodes(p, stem, source)) {
            int idx;
            if (findNode(p, n, idx)) {
                outNames.push_back(n);
                outIndex.push_back(idx);
            }
            else {
                std::cerr << "[CGGA] result node not found: " << n << "\n";
            }
        }
        std::cerr << "[CGGA] result nodes = ";
        for (std::size_t k = 0; k < outNames.size(); ++k) std::cerr << (k ? "," : "") << outNames[k];
        std::cerr << "  (" << source << ")\n";
        std::cerr << "[CGGA] results file = " << resultsPath << "\n";

        hadilan << "Time";
        for (const std::string& n : outNames) hadilan << "\t\tHead_" << n << "_(m)";
        hadilan << "\n";

        // ... step through each time period
        double t = 0;
        double tstep = 0;

        // Adaptive time steps
        double Delta_WH = 0.05, Delta_RWC = 1, Delta_QS = 5.0;
        tstep = Delta_QS;
        if (p.getHydraulicSolverType() == 2)
            tstep = Delta_WH;
        else if (p.getHydraulicSolverType() == 1)
            tstep = Delta_RWC;
        else
            tstep = Delta_QS;
        Epanet::Project::SolverMode currentMode = Epanet::Project::QUASI_STEADY;

        // ... initialize the solver
        std::cout << "\n    Initializing solver ...";
        if ((err = p.initSolver(false))) break;
        std::cout << "\n    ";

        do {
            std::cout << "\r    Solving network at "
                << Utilities::getTime(static_cast<double>(t + tstep)) << " hrs ...        ";
            if (p.getHydraulicSolverType() == 2) { // This is CGGA
                err = p.runSolver(&t);
                if (err) break;
                if (!err) err = p.advanceSolver(&tstep);
                if (err) break;

                if (timeInSeconds) hadilan << t;
                else               hadilan << Utilities::getTime(t);
                for (std::size_t k = 0; k < outIndex.size(); ++k) {
                    double h = 0.0;
                    EN_getNodeValue(outIndex[k], EN_HEAD, &h, &p);
                    hadilan << "\t\t" << h;
                }
                hadilan << "\n";
                hadilan.flush();
            }
            else if (p.getHydraulicSolverType() == 1) { // This is RWC GGA
                err = p.runSolver(&t);
                if (err) break;
                if (!err) err = p.advanceSolver(&tstep);
                if (err) break;
                std::cout << "Time: " << t << " s\n";
            }
            else { // This is quasi-steady
                err = p.runSolver(&t);
                if (err) break;
                if (!err) err = p.advanceSolver(&tstep);
                if (err) break;
            }
        } while (tstep > 0 && !err);

        break;
    }

    // ... simulation was successful
    if (!err)
    {
        clock_t end_t = clock();
        double cpu_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
        std::stringstream ss;
        ss << "\n  Simulation completed in ";
        p.writeMsg(ss.str());
        ss.str("");
        if (cpu_t < 0.001) ss << "< 0.001 sec.";
        else ss << std::setprecision(3) << cpu_t << " sec.";
        p.writeMsg(ss.str());

        std::cout << "\n    Writing report ...                           ";
        err = p.writeReport();
        std::cout << "\n    Simulation completed.                         \n";
        std::cout << "\n... EPANET completed in " << ss.str() << "\n";
    }

    if (err)
    {
        p.writeMsgLog();
        std::cout << "\n\n    There were errors. See report file for details.\n";
        return err;
    }
    return 0;
}


//-----------------------------------------------------------------------------

EN_Project EN_createProject()
{
    Epanet::Project* p = new Epanet::Project();
    return (EN_Project *)p;
}

//-----------------------------------------------------------------------------

int EN_deleteProject(EN_Project p)
{
    delete (Epanet::Project *)p;
    return 0;
}

//-----------------------------------------------------------------------------

int EN_loadProject(const char* fname, EN_Project p)
{
    return project(p)->load(fname);
}

//-----------------------------------------------------------------------------

int EN_saveProject(const char* fname, EN_Project p)
{
    return project(p)->save(fname);
}

//-----------------------------------------------------------------------------

int EN_clearProject(EN_Project p)
{
    project(p)->clear();
    return 0;
}

//-----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////
//  NOT SURE IF THIS METHOD WORKS CORRECTLY -- NEEDS TESTING  //
////////////////////////////////////////////////////////////////
int EN_cloneProject(EN_Project pClone, EN_Project pSource)
{
    if ( pSource == nullptr || pClone == nullptr ) return 102;
    int err = 0;
    std::string tmpFile;
    if ( Utilities::getTmpFileName(tmpFile) )
    {
        try
        {
            EN_saveProject(tmpFile.c_str(), pSource);
            EN_loadProject(tmpFile.c_str(), pClone);
        }
        catch (ENerror const& e)
        {
            static_cast<Epanet::Project*>(pSource)->writeMsg(e.msg);
            err = e.code;
  	    }
        catch (...)
        {
            err = 208; //Unspecified error
        }
        if ( err > 0 )
        {
            EN_clearProject(pClone);
        }
        remove(tmpFile.c_str());
        return err;
    }
    return 208;
}

//-----------------------------------------------------------------------------

/*int EN_runProject(EN_Project p)    // <<=============  TO BE COMPLETED
{
    return 0;
}

//-----------------------------------------------------------------------------

int EN_initSolver(int initFlows, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->initSolver(initFlows);
}

//-----------------------------------------------------------------------------

int EN_runSolver(double* t, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->runSolver(t);
}

//-----------------------------------------------------------------------------

int EN_advanceSolver(int *dt, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->advanceSolver(dt);
}

//-----------------------------------------------------------------------------

int EN_openOutputFile(const char* fname, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->openOutput(fname);
}

//-----------------------------------------------------------------------------

int EN_saveOutput(EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->saveOutput();
}

//-----------------------------------------------------------------------------

int EN_openReportFile(const char* fname, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->openReport(fname);
}

//-----------------------------------------------------------------------------

int EN_writeReport(EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->writeReport();
}

//-----------------------------------------------------------------------------

int EN_writeSummary(EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeSummary();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeResults(int t, EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeResults(t);
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeMsgLog(EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeMsgLog();
    return 0;
} // */

int EN_runProject(EN_Project p)    // <<=============  TO BE COMPLETED
{
    return 0;
}

//-----------------------------------------------------------------------------

int EN_initSolver(int initFlows, EN_Project p)
{
    return project(p)->initSolver(initFlows);
}

//-----------------------------------------------------------------------------

int EN_runSolver(double* t, EN_Project p)
{
    return project(p)->runSolver(t);
}

//-----------------------------------------------------------------------------

int EN_advanceSolver(double* dt, EN_Project p)
{
    return project(p)->advanceSolver(dt);
}

//-----------------------------------------------------------------------------

int EN_openOutputFile(const char* fname, EN_Project p)
{
    return project(p)->openOutput(fname);
}

//-----------------------------------------------------------------------------

int EN_saveOutput(EN_Project p)
{
    return project(p)->saveOutput();
}

//-----------------------------------------------------------------------------

int EN_openReportFile(const char* fname, EN_Project p)
{
    return project(p)->openReport(fname);
}

//-----------------------------------------------------------------------------

int EN_writeReport(EN_Project p)
{
    return project(p)->writeReport();
}

//-----------------------------------------------------------------------------

int EN_writeSummary(EN_Project p)
{
    project(p)->writeSummary();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeResults(int t, EN_Project p)
{
    project(p)->writeResults(t);
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeMsgLog(EN_Project p)
{
    project(p)->writeMsgLog();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_getCount(int element, int* result, EN_Project p)
{
    return DataManager::getCount(element, result, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeIndex(char* name, int* index, EN_Project p)
{
    return DataManager::getNodeIndex(name, index, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeId(int index, char* id, EN_Project p)
{
    return DataManager::getNodeId(index, id, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeType(int index, int* type, EN_Project p)
{
    return DataManager::getNodeType(index, type, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeValue(int index, int param, double* value, EN_Project p)
{
    return DataManager::getNodeValue(index, param, value, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkIndex(char* name, int* index, EN_Project p)
{
    return DataManager::getLinkIndex(name, index, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkId(int index, char* id, EN_Project p)
{
    return DataManager::getLinkId(index, id, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkType(int index, int* type, EN_Project p)
{
    return DataManager::getLinkType(index, type, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkNodes(int index, int* fromNode, int* toNode, EN_Project p)
{
    return DataManager::getLinkNodes(index, fromNode, toNode,
        static_cast<Epanet::Project*>(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkValue(int index, int param, double* value, EN_Project p)
{
   return DataManager::getLinkValue(index, param, value, static_cast<Epanet::Project*>(p)->getNetwork());
}

int EN_setLinkValue(int index, int param, double value, EN_Project p)
{
	return DataManager::setLinkValue(index, param, value, static_cast<Epanet::Project*>(p)->setNetwork());
}

int       EN_getTimeParam(int, long*, EN_Project)
{
    return 0;
}
}  // end of namespace
