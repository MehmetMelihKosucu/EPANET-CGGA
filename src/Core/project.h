/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file project.h
//! \brief Describes EPANET's Project class.


#ifndef PROJECT_H_
#define PROJECT_H_

#include "Core/network.h"
#include "Core/qualengine.h"
#include "Output/outputfile.h"
#include "hydengine.h"
#include "Elements/pipe.h"

#include <string>
#include <fstream>
#include <memory>

// Forward declarations
class Network;
class HydEngine;
class QualEngine;
class OutputFile;

namespace Epanet
{
    //!
    //! \class Project
    //! \brief Encapsulates a pipe network and its simulation engines.
    //!
    //! A project contains a description of the pipe network being analyzed
    //! and the engines and methods used to carry out the analysis. All
    //! methods applied to a project and its components can be done in a
    //! thread-safe manner.
	
    class Project
    {
      public:

        Project();
        ~Project();

        // Define the SolverMode enumeration
        enum SolverMode {
            QUASI_STEADY,
            RWC,
            WATER_HAMMER
        };

        // define thresholds:
        static constexpr double PHI_A_D = 0.5;   // dynamic compressibility threshold
        static constexpr double PHI_R_D = 0.25;  // relative dynamic threshold
        static constexpr double PHI_A_I = 0.1;   // inertial threshold
        static constexpr double PHI_R_I = 0.05;  // relative inertial threshold

        // define time steps:
        static constexpr double Delta_WH = 0.05; // Water Hammer step
        static constexpr double Delta_RWC = 1.0;  // RWC step
        static constexpr double Delta_QS = 10.0; // Quasi-steady step

		MatrixSolver* matrixSolver;  // Pointer to the matrix solver

        int   load(const char* fname);
        int   save(const char* fname);
        void  clear();

        int   initSolver(bool initFlows);
        int   runSolver(double* t);
        int   advanceSolver(double* dt);

        int   openOutput(const char* fname);
        int   saveOutput();

        int   openReport(const char* fname);
        void  writeSummary();
        void  writeResults(int t);
        int   writeReport();

        void  writeMsg(const std::string& msg);
        void  writeMsgLog(std::ostream& out);
        void  writeMsgLog();
        Network* getNetwork() { return &network; }
		Network* setNetwork() { return &network; }
		double pressureManagement(double Xm, double error, double delta_Xm, double Xm_Last, int t);
        int getTotalSimulationTime() const;
       
        SolverMode getCurrentMode() const { return currentMode; }
        void setCurrentMode(SolverMode mode) { currentMode = mode; }
		int Project::getHydraulicSolverType();

        // ... other members ...

        // For example, assume our Project has a Network member named 'network'
        // and containers to store the base (unsegmented) network.
        Network network;

      private:
        //Network        network;        //!< pipe network to be analyzed
        HydEngine      hydEngine;      //!< hydraulic simulation engine.
        QualEngine     qualEngine;     //!< water quality simulation engine.
        OutputFile     outputFile;     //!< binary output file for saved results.
        std::string    inpFileName;    //!< name of project's input file.
        std::string    outFileName;    //!< name of project's binary output file.
        std::string    tmpFileName;    //!< name of project's temporary binary output file.
        std::string    rptFileName;    //!< name of project's report file.
        std::ofstream  rptFile;        //!< reporting file stream.
        double tstep;  // Global timestep for the project
        
        SolverMode currentMode = SolverMode::QUASI_STEADY;

        // Project status conditions
        bool           networkEmpty;
        bool           hydEngineOpened;
        bool           qualEngineOpened;
        bool           outputFileOpened;
        bool           solverInitialized;
        bool           runQuality;
		
        void           finalizeSolver();
        void           closeReport();
    };
}
#endif
