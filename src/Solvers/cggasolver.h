#ifndef CGGASOLVER_H_
#define CGGASOLVER_H_

// Forward declarations
class Network;
class Node;
class Link;
class Junction;
class Pipe;
class PumpCurve;
class MatrixSolver;
class HydSolver;

#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <map>
#include <unordered_set>
#include "Solvers/hydsolver.h"
//#include "Core/hydengine.h"
#include "Core/hydbalance.h"
#include "Core/network.h"
#include "link.h"
#include "node.h"
#include "pipe.h"
#include <unordered_map>
#include <string>
#include "Elements/junction.h"
#include "Elements/valve.h"
#include "Elements/pump.h"
#include "Core/options.h"
#include <set>
#include <map>
#include <unordered_set>
#include "Elements/tank.h"
#include "dualflowhistory.h"
#include <fstream>

using Vector = std::vector<double>;

//! \class CGGASolver
//! \brief A hydraulic solver based on Comprehensive Global Gradient Algorithm.

class CGGASolver : public HydSolver
{
public:

    //! Constructor
    CGGASolver(Network* nw, MatrixSolver* ms);
		
    // Make sure we have a proper destructor
    ~CGGASolver();

    int solve(double& tstep, int& trials, double currentTime);

    // Thresholds (adjustable)
    double phiA_D = 0.008;    // Dynamic threshold [m]
    double phiR_D = 0.004;   // Dynamic threshold [-]
    double phiA_I = 0.001;  // 0.001; // 0.1;    // Inertial threshold [m]
    double phiR_I = 0.0005; // 0.0005; // 0.05;   // Inertial threshold [-]

    std::ofstream indicatorLog;

    // Maximum indicators from last computation (for external monitoring)
    double maxPhiA = 0.0;
    double maxPhiR = 0.0;

    // define time steps:
    static constexpr double Delta_WH = 0.05; // Water Hammer step
    static constexpr double Delta_RWC = 1.0;  // RWC step
    static constexpr double Delta_QS = 10.0; // Quasi-steady step

    // Define the SolverMode enumeration
    enum SolverMode {
        QUASI_STEADY = 0,
        RWC = 1,
        WATER_HAMMER = 2
    };

    SolverMode pendingMode;

    bool pendingModeActive;

    enum SimulationPhase {
        FRONT_END = 0,           // QS before event
		MIDWAY = 1,              // RWC/WH, event detected
        BACK_END = 2,            // QS after event settles
        NONE_PHASE = -1
    };

    SimulationPhase pendingPhase = NONE_PHASE;

    double lastModeChangeTime = 0.0;
    double lastEventTime = -1000.0;
    int modeStepsInCurrentPhase = 0;
    bool eventOngoing = false;

    bool whLockedOut = false;
	bool rwcLockedOut = false;

    double computeTimeStepForMode(SolverMode mode, double currentTime, double endTime);
    
    // Wave attenuation point structure
    struct WaveAttenuationPoint {
        double flow;      // Flow at this point
        double head;      // Head at this point
        double position;  // Relative position along pipe (0-1)
        double amplitude;  // Wave amplitude (for compatibility)

        // Default constructor
        WaveAttenuationPoint() : position(0.0), flow(0.0), head(0.0), amplitude(0.0) {}
        
        // Parameterized constructor
        WaveAttenuationPoint(double pos, double f, double h, double amp = 0.0) 
            : position(pos), flow(f), head(h), amplitude(amp) {}
    };

    struct IntegerSecondState {
        double time = -1.0;
        std::vector<double> nodeHeads;
        std::vector<double> pipeStartFlows;
        std::vector<double> pipeEndFlows;
        std::vector<double> pipeFlows;
        std::vector<double> pipeFrictionLoss;  // Store friction loss too
        bool isValid = false;
    };

    IntegerSecondState prevSecondState;
    IntegerSecondState currentSecondState;
    bool evaluateModeThisStep = false;

    // Logging helper
    void logIndicators(double currentTime, double dt, bool modeEval);

    std::map<Pipe*, std::tuple<double, double, double>> transitionMapCache;

    // Current state
    SolverMode currentMode;
    SolverMode previousMode = QUASI_STEADY;
    SimulationPhase simulationPhase = FRONT_END;

    double previousTimestep;

    // Forced step counters (prevents cycling)
    int forcedWHSteps;      // Current WH counter
    int forcedWHSteps0;     // Initial WH counter value (e.g., 10)
    int forcedRWCSteps;     // Current RWC counter
    int forcedRWCSteps0;    // Initial RWC counter value (e.g., 5)

    // Event tracking
    double eventStartTime;

    // Event detection
    bool detectHydraulicEvent(double currentTime);
    bool isEventOngoing(double currentTime);

    // Mode/Phase helpers
    const char* CGGASolver::getModeString(SolverMode mode);
    const char* CGGASolver::getPhaseString(SimulationPhase phase);

    // Indicator computation
    void computeFlowIndicators(double currentTime, double tstep);

    // Solver execution
    int executeWaterHammerSolve(double tstep, int& trials, double currentTime);
    int CGGASolver::executeRWCSolve(double tstep_advance, int& trials, double currentTime);
    int executeQuasiSteadySolve(double tstep, int& trials, double currentTime);

    std::vector<double> dQ_start;
    std::vector<double> dQ_end;

    // Accessor methods for pending mode state (used by HydEngine::advance)
    bool hasPendingModeChange() const { return pendingModeChange; }
    SolverMode getPendingMode() const { return pendingMode; }

    // Make currentMode accessible (if not already public)
    SolverMode getCurrentMode() const { return currentMode; }

    double minErrorNorm;      //!< Minimum observed error norm
    double dl;

private:

    double lamda;
    static const int ARF = 1;  // Adaptive Relaxation Factor
    static const int BRF = 2;  // Backward Relaxation Factor

    // General solver variables
    int    nodeCount;         //!< Number of network nodes
    int    linkCount;         //!< Number of network links
    int    hLossEvalCount;    //!< Number of head loss evaluations
    int    stepSizing;        //!< Newton step sizing method
    int    trialsLimit;       //!< Limit on the number of trials
    bool   reportTrials;      //!< Flag to report summary of each trial
    double tstep;             //!< Time step (sec)
    double errorNorm;         //!< Solution error norm
    double oldErrorNorm;      //!< Previous error norm
    
	double theta;             //!< Relaxation factor for step sizing
	double kappa;             //!< Temporal discretization parameter
	double fixedHead;         //!< Fixed head value for a node

    int CGGASolver::sign(double val);

    bool pendingModeChange = false;

    double minSegmentLength;

    HydBalance hydBalance;  // Instance of the HydBalance class
	
    std::vector<double> phiA; //!< Absolute flow indicators for each link
    std::vector<double> phiR; //!< Relative flow indicators for each link
    std::vector<int>    SC;   //!< Compressibility indicator for each link
    std::vector<int>    SI;   //!< Inertial effects indicator for each link

    // Hydraulic variables
    std::vector<double> dH;   //!< Head change at each node (ft)
    std::vector<double> dQ;   //!< Flow change in each link (cfs)
    std::vector<double> xQ;   //!< Node flow imbalances (cfs)

    // Iteration state storage (for line search)
    std::vector<double> iterHeads;
    std::vector<double> iterFlows;
    std::vector<double> iterStartFlows;
    std::vector<double> iterEndFlows;

    // Convergence parameters
    double headErrLimit;      //!< Allowable head error (ft)
    double flowErrLimit;      //!< Allowable flow error (cfs)
    double flowChangeLimit;   //!< Allowable flow change (cfs)
    double flowRatioLimit;    //!< Allowable total flow change / total flow
    
    // Water hammer specific error tracking
    double maxPressureWaveErr;
    double maxCharEqErr;

    void restoreState();
    double calculateTrialError(double lambda, double currentTime, double tstep);

    // Functions for solver initialization and convergence
    bool hasConverged();
    bool linksChangedStatus();

    // Line search state
    void storeIterationState();
    void restoreIterationState();

    // Utility
    std::string getModeString(SolverMode mode) const;

    // Functions to assemble and solve equations
	void setFixedGradeNodes();  //!< Adjust fixed grade status of specific nodes
   
    void setSteadyMatrixCoeffs();       //!< Set coefficients for hybrid equations
    void setSteadyLinkCoeffs(); 	 //!< Set link coefficients
    void setSteadyNodeCoeffs();	   //!< Set node coefficients
    void setSteadyValveCoeffs();	  //!< Set valve coefficients

    void setUnsteadyMatrixCoeffs(double currentTime);       //!< Set coefficients for hybrid equations
    void setUnsteadyLinkCoeffs(double currentTime); 	 //!< Set link coefficients
    void setUnsteadyNodeCoeffs(double currentTime);	   //!< Set node coefficients
    void setUnsteadyValveCoeffs(double currentTime);	  //!< Set valve coefficients
    
    void handleCompressibleFlowCoeffs(Link* link, Node* node1, Node* node2, double currentTime, int n1, int n2);
    void handleQuasiSteadyFlowCoeffs(Link* link, Node* node1, Node* node2, int n1, int n2);

    int findSteadyHeadChanges();       //!< Solve for nodal head changes in Steady State Condition
    void findSteadyFlowChanges();       //!< Solve for link flow changes in Steady State Condition

	// Functions that update the hydraulic solution for RWC Condition
    int    findRWCHeadChanges(double currentTime, double tstep, double rwcTheta);
    void   findRWCFlowChanges(double tstep, double rwcTheta);
    double findRWCStepSize(int trials, double currentTime, double tstep, double rwcTheta);
    void   updateRWCSolution(double lamda);
    double findRWCErrorNorm(double lamda, double currentTime, double tstep);
    void CGGASolver::setRWCMatrixCoeffs(double tstep);
    void CGGASolver::setRWCLinkCoeffs(double tstep);
	void CGGASolver::setRWCNodeCoeffs(double tstep);
	void CGGASolver::setRWCValveCoeffs();

    int findUnsteadyHeadChanges(double currentTime);
    void findUnsteadyFlowChanges (int currentMode, double currentTime);
    double findUnsteadyErrorNorm(double lamda, double currentTime, double tstep);

    // Functions for solution updates and error handling
    double findStepSize(int trials, double currentTime);
    //double findWHStepSize(int trials, double currentTime, double dl, double minErrorNorm);
    double CGGASolver::findWHStepSize(int trials, double currentTime, double dl, double minErrorNorm);
    int lambdaNumber;
    std::vector<double> Lambda;
	double findErrorNorm(double lamda, double currentTime, double tstep, double currentMode);
    void updateSolution(double lamda, bool useShadowMode = false );
    void reportTrial(int trials, double lamda);
   
    void setConvergenceLimits();

    std::vector<double> currentHeads;  // Current head values at nodes
    std::vector<double> currentFlows;  // Current flow values in links

    void CGGASolver::handleQuasiSteadyFlowChange(Link* link, double H_A, double H_B, int linkIndex);
    
    void CGGASolver::updateWHSolution(double lamda);

    // Shadow calculation vectors
    std::vector<double> shadowHeads;
    std::vector<double> shadowFlows;
    
    void CGGASolver::configurePipeForTimestep(Pipe* pipe, double tstep, Network* network);
    double* savedNodeHeads;     // Temporary storage for node heads
    double* savedLinkFlows;     // Temporary storage for link flows
    double* savedDH;            // Temporary storage for head changes
    double* savedDQ;            // Temporary storage for flow changes
    
    // State storage for adaptive relaxation factor
    std::vector<double> saved_nodeHeads;      // Store node heads
    std::vector<double> saved_linkFlows;      // Store link flows
    std::vector<double> saved_startFlows;     // Store pipe start flows
    std::vector<double> saved_endFlows;       // Store pipe end flows
    
    // Store the changes as well (important for reapplying with different lambda)
    std::vector<double> saved_dH;             // Store head changes
    std::vector<double> saved_dQ;             // Store flow changes
    std::vector<double> saved_dQ_start;       // Store start flow changes
    std::vector<double> saved_dQ_end;         // Store end flow changes
    double CGGASolver::findUnsteadyStepSize(int trials, double currentTime, double tstep, double prevErrorNorm, double dl);
};
#endif // CGGASOLVER_H