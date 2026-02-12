/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

//! \file hydengine.h
//! \brief Describes the HydEngine class.

#ifndef HYDENGINE_H_
#define HYDENGINE_H_

#include <string>

#include "network.h"
#include "dualflowhistory.h"


//class Project;
class Network;
class HydSolver;
class MatrixSolver;

//! \class HydEngine
//! \brief Simulates extended period hydraulics.
//!
//! The HydEngine class carries out an extended period hydraulic simulation on
//! a pipe network, calling on its HydSolver object to solve the conservation of
//! mass and energy equations at each time step.

class HydEngine
{
  public:

    // Constructor/Destructor

    HydEngine();

    ~HydEngine();

    // Public Methods

    void   open(Network* nw);
    void   init(bool initFlows);
    int    solve(double* t);
    void   advance(double* tstep);
    void   close();

    int    getElapsedTime() { return currentTime; }
    double getPeakKwatts()  { return peakKwatts;  }
	double    currentTime;        //!< current simulation time (sec)
    void HydEngine::updateFlowHistory(double currentTime);
    double HydEngine::calculateMaxWaveTravel();

  private:
      
    //Project* project;  // Pointer to the Project instance
    Network* network;            //!< network being analyzed
    // Engine state
    enum EngineState {CLOSED, OPENED, INITIALIZED};
    EngineState engineState;
    double lastHistoryUpdateTime = -1.0;  // Track when we last updated history
    // Engine components

    HydSolver*     hydSolver;          //!< steady state or rwc unsteady hydraulic solver
    MatrixSolver*  matrixSolver;       //!< sparse matrix solver

    // Engine properties
	double 	       tstep; 			//!< current time step (sec)
    bool           saveToFile;         //!< true if results saved to file
    bool           halted;             //!< true if simulation has been halted
    double            startTime;          //!< starting time of day (sec)
    double            rptTime;            //!< current reporting time (sec)
    double            hydStep;            //!< hydraulic time step (sec)
    
    double            timeOfDay;          //!< current time of day (sec)
    double         peakKwatts;         //!< peak energy usage (kwatts)
    std::string    timeStepReason;     //!< reason for taking next time step

    void           initMatrixSolver();

    int            getTimeStep();
	void           pastJunction();
	void           pastLink();
    int            timeToPatternChange(double tstep);
    int            timeToActivateControl(double tstep);
    int            timeToCloseTank(double tstep);

    void           updateCurrentConditions();
    void           updateTanks();
    void           updatePatterns();
    void           updateEnergyUsage();

    bool           isPressureDeficient();
    int            resolvePressureDeficiency(int& trials);
    void           reportDiagnostics(int statusCode, int trials);
    
};

#endif
