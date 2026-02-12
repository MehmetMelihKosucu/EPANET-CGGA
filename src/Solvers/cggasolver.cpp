// Ensure that all necessary headers are included
#include "cggasolver.h"
#include "matrixsolver.h"
#include "Core/network.h"
#include "Core/constants.h"
#include "Elements/control.h"
#include "Core/hydbalance.h"
#include "Core/options.h" 
#include "Elements/junction.h"
#include "Elements/tank.h"
#include "Elements/valve.h"
#include "Elements/node.h"
#include "Elements/link.h"
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Elements/pumpcurve.h"
#include "Elements/curve.h"
#include "Elements/reservoir.h"
#include "Elements/emitter.h"
#include "Models/headlossmodel.h"
#include "Models/leakagemodel.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>  // Required for std::accumulate
#include "sparspak.h"
#include "Models/demandmodel.h"
#include "sparspaksolver.h"
#include "project.h"
#include <unordered_map>
#include <string>
#include <set>
#include <unordered_set>
#include <queue>
#include "dualflowhistory.h"

#include <iomanip>
using namespace std;
static const string s_Trial = "    Trial ";
static const string s_IllConditioned = "  Hydraulic matrix ill-conditioned at node ";
static const string s_StepSize = "    Step Size   = ";
static const string s_TotalError = "    Error Norm  = ";
static const string s_HlossEvals = "    Head Loss Evaluations = ";
static const string s_HeadError = "    Head Error  = ";
static const string s_ForLink = " for Link ";
static const string s_FlowError = "    Flow Error  = ";
static const string s_AtNode = " at Node ";
static const string s_FlowChange = "    Flow Change = ";
static const string s_TotFlowChange = "    Total Flow Change Ratio = ";
static const string s_NodeLabel = "  Node ";
static const string s_FGChange = "    Fixed Grade Status changed to ";

//-----------------------------------------------------------------------------

// error norm threshold for status checking
static const double ErrorThreshold = 1e-6;
static const double Huge = numeric_limits<double>::max();

// step sizing enumeration
enum StepSizing { FULL, RELAXATION, LINESEARCH, BRF, ARF };

// Updated constructor with network state management initialization
CGGASolver::CGGASolver(Network* nw, MatrixSolver* ms) : HydSolver(nw, ms),
    hydBalance(),
    fixedHead(50.0),
    previousMode(QUASI_STEADY),    // ← Use enum value, not int
    currentMode(QUASI_STEADY),     // ← Add this if not present
    simulationPhase(FRONT_END)     // ← Add this if not present
{  

    if (!nw || !ms) {
        throw std::invalid_argument("Network or MatrixSolver cannot be null.");
    }

    // Store network and matrix solver pointers
    network = nw;
    matrixSolver = ms;

    nodeCount = nw->nodes.size();
    linkCount = nw->links.size();

    // Resize vectors to match network dimensions
    dH.resize(nodeCount, 0.0);    dQ.resize(linkCount, 0.0);    xQ.resize(nodeCount, 0.0);    
    dQ_start.resize(linkCount, 0.0);    dQ_end.resize(linkCount, 0.0);
    phiA.resize(linkCount, 0.0);
    phiR.resize(linkCount, 0.0);

    // Initialize convergence parameters
    errorNorm = std::numeric_limits<double>::max();
    headErrLimit = 0.1;  // Default head error limit
    flowErrLimit = 0.5;  // Default flow error limit
    flowChangeLimit = 0.1; // Default flow change limit

    savedNodeHeads = new double[nodeCount];
    savedLinkFlows = new double[linkCount];
    savedDH = new double[nodeCount];
    savedDQ = new double[linkCount];

    pendingMode = QUASI_STEADY;
    pendingModeActive = false;

    if (network->option(Options::STEP_SIZING) == "RELAXATION" )
        stepSizing = RELAXATION;
    else if (network->option(Options::STEP_SIZING) == "LINESEARCH" )
        stepSizing = LINESEARCH;
    else if (network->option(Options::STEP_SIZING) == "BRF")
        stepSizing = BRF;
    else if (network->option(Options::STEP_SIZING) == "ARF")
        stepSizing = ARF;
    else stepSizing = FULL;

    // Log initialization (optional)
    std::cout << "CGGASolver initialized with "
        << nodeCount << " active nodes and "
        << linkCount << " active links." << std::endl;
}

// Updated destructor to clean up all resources
CGGASolver::~CGGASolver() {
    // Clean up vectors
    dH.clear();
    dQ.clear();
    xQ.clear();
    SC.clear();
    SI.clear();
    phiA.clear();
    phiR.clear();

    delete[] savedNodeHeads;
    delete[] savedLinkFlows;
    delete[] savedDH;
    delete[] savedDQ;
    
}

//=============================================================================
// CGGA SOLVER 
//=============================================================================

int CGGASolver::solve(double& tstep_, int& trials, double currentTime) {

    trials = 1;
    int result = HydSolver::SUCCESSFUL;

    // Integration parameters
    theta = network->option(Options::TIME_WEIGHT);
    theta = std::max(0.5, std::min(theta, 1.0));
    kappa = network->option(Options::TEMP_DISC_PARA);
    kappa = std::max(0.0, std::min(kappa, 1.0));

    double simulationEndTime = network->option(Options::TOTAL_DURATION);

    // Check simulation end
    if (currentTime >= simulationEndTime - 1e-10) {
        network->msgLog << "\nSimulation reached end time";
        return HydSolver::SUCCESSFUL;
    }

    // ========================================================================
    // CHECK: Are we at an integer second?
    // Mode evaluation and indicator computation only happen at integer seconds
    // ========================================================================
    const double INTEGER_TOLERANCE = 1e-6;
    double fractionalPart = currentTime - std::floor(currentTime);
    bool isIntegerSecond = (fractionalPart < INTEGER_TOLERANCE) ||
        (fractionalPart > 1.0 - INTEGER_TOLERANCE);

    //=========================================================================
    // STEP 0: APPLY PENDING MODE CHANGE FROM PREVIOUS TIMESTEP
    //=========================================================================

    if (pendingModeChange) {
        SolverMode oldMode = currentMode;
        currentMode = pendingMode;
        pendingModeChange = false;

        network->msgLog << "\n==================================================";
        network->msgLog << "\n  APPLYING MODE CHANGE at t=" << currentTime;
        network->msgLog << "\n  " << oldMode << " -> " << currentMode;
        network->msgLog << "\n==================================================";

        if (currentMode == WATER_HAMMER) {
            forcedWHSteps = forcedWHSteps0;
        }
        else if (currentMode == RWC) {
            forcedRWCSteps = forcedRWCSteps0;
        }
    }

    //=========================================================================
    // STEP 1: INITIAL SOLVE AT t=0
    //=========================================================================

    const double initialTimeTolerance = 1e-9;
    if (currentTime < initialTimeTolerance) {
        network->msgLog << "\n==================================================";
        network->msgLog << "\n  INITIAL SOLVE at t=0 (Quasi-Steady)";
        network->msgLog << "\n==================================================";

        currentMode = QUASI_STEADY;
        simulationPhase = FRONT_END;
        whLockedOut = false;
		rwcLockedOut = false;
        pendingModeChange = false;

        tstep = computeTimeStepForMode(QUASI_STEADY, currentTime, simulationEndTime);
        result = executeQuasiSteadySolve(tstep, trials, currentTime);
        if (result != HydSolver::SUCCESSFUL) return result;

        // Initialize indicators to zero at t=0
        maxPhiA = 0.0;
        maxPhiR = 0.0;
        previousMode = QUASI_STEADY;
        previousTimestep = tstep;
        tstep_ = tstep;

        return result;
    }

    //=========================================================================
    // NON-INTEGER SECOND: Just execute current mode solver, no mode changes
    //=========================================================================

    if (!isIntegerSecond) {
        // Continue with current mode, no indicator computation or mode evaluation

        if (currentMode == WATER_HAMMER) {
            tstep = computeTimeStepForMode(WATER_HAMMER, currentTime, simulationEndTime);
            result = executeWaterHammerSolve(tstep, trials, currentTime);
            tstep_ = tstep;
        }
        else if (currentMode == RWC) {
            tstep = computeTimeStepForMode(RWC, currentTime, simulationEndTime);
            result = executeRWCSolve(tstep, trials, currentTime);
            tstep_ = tstep;
        }
        else {  // QUASI_STEADY
            tstep = computeTimeStepForMode(QUASI_STEADY, currentTime, simulationEndTime);
            result = executeQuasiSteadySolve(tstep, trials, currentTime);
            tstep_ = tstep;
        }

        previousTimestep = tstep;
        return result;
    }

    //=========================================================================
    // INTEGER SECOND: Full mode evaluation logic
    //=========================================================================

    //=========================================================================
    // STEP 2: EVENT DETECTION (only at integer seconds)
    //=========================================================================

    bool newEventDetected = detectHydraulicEvent(currentTime);
    bool eventOngoing = isEventOngoing(currentTime);

    // New event clears WH lockout (ratchet-down reset)
    if (newEventDetected && whLockedOut) {
        network->msgLog << "\n  NEW hydraulic event detected - clearing WH lockout";
        whLockedOut = false;
    }

    if (newEventDetected && rwcLockedOut) {
        network->msgLog << "\n  NEW hydraulic event detected - clearing RWC lockout";
        rwcLockedOut = false;
    }

    //=========================================================================
    // STEP 3: EVENT DETECTED - Enter MIDWAY, evaluate indicators for next mode
    //=========================================================================

    if ((simulationPhase == FRONT_END || simulationPhase == BACK_END) && newEventDetected) {

        network->msgLog << "\n==================================================";
        network->msgLog << "\n  EVENT DETECTED at t=" << currentTime;
        network->msgLog << "\n  Entering MIDWAY phase";
        network->msgLog << "\n==================================================";

        // Enter MIDWAY phase
        simulationPhase = MIDWAY;
        eventStartTime = currentTime;
        whLockedOut = false;

        // Complete THIS timestep with current mode (QS)
        tstep = computeTimeStepForMode(QUASI_STEADY, currentTime, simulationEndTime);
        result = executeQuasiSteadySolve(tstep, trials, currentTime);
        if (result != HydSolver::SUCCESSFUL) return result;

        previousTimestep = tstep;
        previousMode = QUASI_STEADY;

        // Compute indicators and let them decide the next mode
        computeFlowIndicators(currentTime, tstep);

        bool aboveDynamic = (maxPhiA > phiA_D) && (maxPhiR > phiR_D);
        bool aboveInertial = (maxPhiA >= phiA_I) && (maxPhiR >= phiR_I);

        if (aboveDynamic) {
            pendingModeChange = true;
            pendingMode = WATER_HAMMER;
            tstep_ = computeTimeStepForMode(WATER_HAMMER, currentTime + tstep, simulationEndTime);

            network->msgLog << "\n  Indicators above dynamic -> scheduling WH";
            network->msgLog << "\n  phiA=" << maxPhiA << ", phiR=" << maxPhiR;
        }
        else if (aboveInertial) {
            pendingModeChange = true;
            pendingMode = RWC;
            tstep_ = computeTimeStepForMode(RWC, currentTime + tstep, simulationEndTime);

            network->msgLog << "\n  Indicators above inertial -> scheduling RWC";
            network->msgLog << "\n  phiA=" << maxPhiA << ", phiR=" << maxPhiR;
        }
        else {
            // Indicators still low, stay in MIDWAY but continue QS for now
            tstep_ = tstep;

            network->msgLog << "\n  Indicators still low -> continuing QS in MIDWAY";
            network->msgLog << "\n  phiA=" << maxPhiA << ", phiR=" << maxPhiR;
        }

        return result;
    }

    //=========================================================================
    // STEP 4: MIDWAY MICROSIMULATIONS (mode evaluation at integer seconds)
    //=========================================================================

    if (simulationPhase == MIDWAY) {

        //---------------------------------------------------------------------
        // CASE A: Currently in WATER HAMMER mode
        //---------------------------------------------------------------------
        if (currentMode == WATER_HAMMER) {

            tstep = computeTimeStepForMode(WATER_HAMMER, currentTime, simulationEndTime);
            result = executeWaterHammerSolve(tstep, trials, currentTime);
            if (result != HydSolver::SUCCESSFUL) return result;

            // Compute indicators at integer second
            computeFlowIndicators(currentTime, tstep);
            forcedWHSteps--;

            // Check for WH -> RWC transition
            bool belowDynamic = (maxPhiA < phiA_D) || (maxPhiR < phiR_D);

            if (belowDynamic && forcedWHSteps <= 0) {

                network->msgLog << "\n==================================================";
                network->msgLog << "\n  [t=" << currentTime << "] Scheduling WH->RWC for next timestep";
                network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
                network->msgLog << "\n  WH now LOCKED OUT (ratchet-down active)";
                network->msgLog << "\n==================================================";

                // Schedule RWC for next timestep
                pendingModeChange = true;
                pendingMode = RWC;
                whLockedOut = true;

                // Return RWC timestep for next iteration
                tstep_ = computeTimeStepForMode(RWC, currentTime + tstep, simulationEndTime);
            }
            else {
                tstep_ = tstep;

                if (belowDynamic && forcedWHSteps > 0) {
                    network->msgLog << "\n  [WH] Indicators low but " << forcedWHSteps
                        << " forced steps remain";
                }
            }

            previousTimestep = tstep;
            previousMode = WATER_HAMMER;
            return result;
        }

        //---------------------------------------------------------------------
        // CASE B: Currently in RWC mode
        //---------------------------------------------------------------------
        else if (currentMode == RWC) {

            tstep = computeTimeStepForMode(RWC, currentTime, simulationEndTime);
            result = executeRWCSolve(tstep, trials, currentTime);
            if (result != HydSolver::SUCCESSFUL) return result;

            // Compute indicators at integer second
            computeFlowIndicators(currentTime, tstep);
            forcedRWCSteps--;

            bool aboveDynamic = (maxPhiA > phiA_D) && (maxPhiR > phiR_D);
            bool belowInertial = (maxPhiA < phiA_I) || (maxPhiR < phiR_I);

            //--- RWC -> WH transition (only if not locked out) ---
            if (aboveDynamic && !whLockedOut) {

                network->msgLog << "\n==================================================";
                network->msgLog << "\n  [t=" << currentTime << "] Scheduling RWC->WH for next timestep";
                network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
                network->msgLog << "\n==================================================";

                // Schedule WH for next timestep
                pendingModeChange = true;
                pendingMode = WATER_HAMMER;

                // Return WH timestep for next iteration
                tstep_ = computeTimeStepForMode(WATER_HAMMER, currentTime + tstep, simulationEndTime);

                previousTimestep = tstep;
                previousMode = RWC;
                return result;
            }

            //--- RWC -> QS transition (exit MIDWAY) ---
            if (!eventOngoing && belowInertial && forcedRWCSteps <= 0) {

                network->msgLog << "\n==================================================";
                network->msgLog << "\n  [t=" << currentTime << "] Scheduling RWC->QS (exiting MIDWAY)";
                network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
                network->msgLog << "\n  Entering BACK-END phase";
                network->msgLog << "\n  RWC now LOCKED OUT (ratchet-down active)";
                network->msgLog << "\n==================================================";

                // Schedule QS for next timestep
                pendingModeChange = true;
                pendingMode = QUASI_STEADY;
                simulationPhase = BACK_END;
                whLockedOut = false;
                rwcLockedOut = true;

                // Return QS timestep for next iteration
                tstep_ = computeTimeStepForMode(QUASI_STEADY, currentTime + tstep, simulationEndTime);
            }
            else {
                tstep_ = tstep;

                if (aboveDynamic && whLockedOut) {
                    network->msgLog << "\n  [RWC] Indicators high but WH LOCKED OUT";
                }
                if (eventOngoing && belowInertial) {
                    network->msgLog << "\n  [RWC] BOTH indicators below inertial but event still ongoing";
                }
            }

            previousTimestep = tstep;
            previousMode = RWC;
            return result;
        }

        //---------------------------------------------------------------------
        // CASE C: Currently in QUASI-STEADY mode during MIDWAY
        //---------------------------------------------------------------------
        else {
            tstep = computeTimeStepForMode(QUASI_STEADY, currentTime, simulationEndTime);
            result = executeQuasiSteadySolve(tstep, trials, currentTime);
            if (result != HydSolver::SUCCESSFUL) return result;

            // Compute indicators at integer second
            computeFlowIndicators(currentTime, tstep);

            bool aboveDynamic = (maxPhiA > phiA_D) && (maxPhiR > phiR_D);
            bool belowInertial = (maxPhiA < phiA_I) && (maxPhiR < phiR_I);
            bool aboveInertial = !belowInertial;

            if (aboveDynamic && !whLockedOut) {
                network->msgLog << "\n==================================================";
                network->msgLog << "\n  [t=" << currentTime << "] Scheduling QS->WH in MIDWAY";
                network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
                network->msgLog << "\n==================================================";

                // Schedule WH for next timestep
                pendingModeChange = true;
                pendingMode = WATER_HAMMER;

                // Return WH timestep for next iteration
                tstep_ = computeTimeStepForMode(WATER_HAMMER, currentTime + tstep, simulationEndTime);
            }
            else if (aboveInertial && !rwcLockedOut) {
                network->msgLog << "\n==================================================";
                network->msgLog << "\n  [t=" << currentTime << "] Scheduling QS->RWC in MIDWAY";
                network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
                network->msgLog << "\n==================================================";

                // Schedule RWC for next timestep
                pendingModeChange = true;
                pendingMode = RWC;

                // Return RWC timestep for next iteration
                tstep_ = computeTimeStepForMode(RWC, currentTime + tstep, simulationEndTime);
            }
            else if (!eventOngoing && belowInertial) {
                network->msgLog << "\n  [t=" << currentTime << "] QS->BACK-END (event complete, indicators low)";

                simulationPhase = BACK_END;
                whLockedOut = false;
                tstep_ = tstep;
            }
            else {
                tstep_ = tstep;
            }

            previousTimestep = tstep;
            previousMode = QUASI_STEADY;
            return result;
        }
    }

    //=========================================================================
    // STEP 5: FRONT-END / BACK-END (QUASI-STEADY)
    //=========================================================================

    tstep = computeTimeStepForMode(QUASI_STEADY, currentTime, simulationEndTime);
    result = executeQuasiSteadySolve(tstep, trials, currentTime);
    if (result != HydSolver::SUCCESSFUL) return result;

    // Compute indicators at integer second
    computeFlowIndicators(currentTime, tstep);

    bool aboveInertial = (maxPhiA >= phiA_I) && (maxPhiR >= phiR_I);

    if (aboveInertial && !rwcLockedOut) {

        network->msgLog << "\n==================================================";
        network->msgLog << "\n  [t=" << currentTime << "] Scheduling QS->RWC";
        network->msgLog << "\n  Indicators: phiA=" << maxPhiA << ", phiR=" << maxPhiR;
        network->msgLog << "\n  Entering MIDWAY phase";
        network->msgLog << "\n==================================================";

        // Schedule RWC for next timestep
        pendingModeChange = true;
        pendingMode = RWC;
        simulationPhase = MIDWAY;

        // Return RWC timestep for next iteration
        tstep_ = computeTimeStepForMode(RWC, currentTime + tstep, simulationEndTime);
    }
    else {
        tstep_ = tstep;
    }

    previousTimestep = tstep;
    previousMode = QUASI_STEADY;

    return result;
}


const char* CGGASolver::getModeString(SolverMode mode) {
    switch (mode) {
    case QUASI_STEADY: return "QUASI_STEADY";
    case RWC: return "RWC";
    case WATER_HAMMER: return "WATER_HAMMER";
    default: return "NONE";
    }
}

const char* CGGASolver::getPhaseString(SimulationPhase phase) {
    switch (phase) {
    case FRONT_END: return "FRONT_END";
    case MIDWAY: return "MIDWAY";
    case BACK_END: return "BACK_END";
    default: return "NONE";
    }
}

//=============================================================================
// COMPUTE TIMESTEP FOR MODE
//=============================================================================
double CGGASolver::computeTimeStepForMode(SolverMode mode, double currentTime, double endTime)
{
    double dt;

    switch (mode) {
    case WATER_HAMMER:
        dt = Delta_WH;  // e.g., 0.05s
        break;

    case RWC:
        dt = Delta_RWC;  // e.g., 1.0s
        break;

    case QUASI_STEADY:
    default:
        dt = Delta_QS;  // e.g., 1.0s or larger
        break;
    }

    // Don't exceed remaining simulation time
    if (currentTime + dt > endTime) {
        dt = endTime - currentTime;
    }

    return std::max(dt, 1e-6);  // Ensure positive timestep
}

bool CGGASolver::isEventOngoing(double currentTime) {
    /*
     * Checks if any boundary condition is still actively changing.
     * This is used to determine when the MIDWAY phase can end.
     *
     * Returns true if any valve/pump is still in transition.
     */

    const double settingTolerance = 1e-6;

    // Check if any valve is still moving
    for (Link* link : network->links) {
        if (link->type() == Link::VALVE) {
            Valve* valve = static_cast<Valve*>(link);

            if (valve->settingPattern != nullptr) {
                double currentSetting = valve->setting;
                double previousSetting = valve->pastSetting;

                // Valve is still moving if setting changed this step
                if (std::abs(currentSetting - previousSetting) > settingTolerance) {
                    return true;
                }
            }
        }
    }

    // Check if any pump is changing
    for (Link* link : network->links) {
        if (link->type() == Link::PUMP) {
            Pump* pump = static_cast<Pump*>(link);

            if (pump->speedPattern != nullptr) {
                double currentSpeed = pump->speed;
                double previousSpeed = pump->pastSpeed;

                if (std::abs(currentSpeed - previousSpeed) > settingTolerance) {
                    return true;
                }
            }
        }
    }

    return false;
}

void CGGASolver::computeFlowIndicators(double currentTime, double dt) {
    // ========================================================================
    // UNSTEADINESS INDICATOR COMPUTATION - Integer Second Evaluation Only
    //
    // Key insight: Indicators WILL oscillate during active transients because
    // the physical quantities (dH/dt, DQ/Dt) oscillate with pressure waves.
    // This is correct behavior. What matters for mode selection is whether
    // indicators are consistently ABOVE or BELOW the thresholds.
    //
    // Important: We compute friction loss F(Q) directly from pipe properties
    // and flow, NOT from pipe->hLoss which may have different meaning in
    // different solver modes.
    // ========================================================================

    const double g = 32.2;              // ft/s²
    const double FT_TO_M = 0.3048;      // Conversion factor
    const double MIN_DENOM = 1e-10;
    const double INTEGER_TOLERANCE = 1e-6;

    // Initialize vectors if needed
    if (phiA.size() != linkCount) {
        phiA.resize(linkCount, 0.0);
        phiR.resize(linkCount, 0.0);
    }

    // Reset mode evaluation flag
    evaluateModeThisStep = false;

    // ------------------------------------------------------------------------
    // CHECK: Are we at an exact integer second?
    // ------------------------------------------------------------------------
    double fractionalPart = currentTime - std::floor(currentTime);
    bool isIntegerSecond = (fractionalPart < INTEGER_TOLERANCE) ||
        (fractionalPart > 1.0 - INTEGER_TOLERANCE);

    if (!isIntegerSecond) {
        return;
    }

    // ------------------------------------------------------------------------
    // We are at an integer second - capture state and compute indicators
    // ------------------------------------------------------------------------
    int currentSecond = static_cast<int>(std::round(currentTime));

    // Shift previous state and capture current state
    prevSecondState = currentSecondState;

    currentSecondState.time = static_cast<double>(currentSecond);
    currentSecondState.isValid = true;

    int nodeCount = network->count(Element::NODE);
    if (currentSecondState.nodeHeads.size() != nodeCount) {
        currentSecondState.nodeHeads.resize(nodeCount, 0.0);
    }
    if (currentSecondState.pipeStartFlows.size() != linkCount) {
        currentSecondState.pipeStartFlows.resize(linkCount, 0.0);
        currentSecondState.pipeEndFlows.resize(linkCount, 0.0);
        currentSecondState.pipeFlows.resize(linkCount, 0.0);
    }

    // Store current node heads
    for (Node* node : network->nodes) {
        if (node && node->index >= 0 && node->index < nodeCount) {
            currentSecondState.nodeHeads[node->index] = node->head;
        }
    }

    // Store current pipe flows
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) continue;

        int idx = pipe->index;
        if (idx < 0 || idx >= linkCount) continue;

        currentSecondState.pipeStartFlows[idx] = pipe->startFlow;
        currentSecondState.pipeEndFlows[idx] = pipe->endFlow;
        currentSecondState.pipeFlows[idx] = pipe->flow;
    }

    // ------------------------------------------------------------------------
    // Check if we have valid history
    // ------------------------------------------------------------------------
    if (!prevSecondState.isValid && currentSecondState.time - prevSecondState.time >= 1.0) {
        // Initialize prevSecondState from past values stored on nodes/pipes
        prevSecondState.time = currentSecondState.time - 1.0;  // Assume 1 second ago
        prevSecondState.isValid = true;

        // Resize if needed
        if (prevSecondState.nodeHeads.size() != nodeCount) {
            prevSecondState.nodeHeads.resize(nodeCount, 0.0);
        }
        if (prevSecondState.pipeStartFlows.size() != linkCount) {
            prevSecondState.pipeStartFlows.resize(linkCount, 0.0);
            prevSecondState.pipeEndFlows.resize(linkCount, 0.0);
            prevSecondState.pipeFlows.resize(linkCount, 0.0);
        }

        // Fill prevSecondState with pastHead values from nodes
        for (Node* node : network->nodes) {
            if (node && node->index >= 0 && node->index < nodeCount) {
                prevSecondState.nodeHeads[node->index] = node->pastHead;
            }
        }

        // Fill prevSecondState with pastFlow values from pipes
        for (Link* link : network->links) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue;

            int idx = pipe->index;
            if (idx < 0 || idx >= linkCount) continue;

            prevSecondState.pipeStartFlows[idx] = pipe->pastStartFlow;
            prevSecondState.pipeEndFlows[idx] = pipe->pastEndFlow;
            prevSecondState.pipeFlows[idx] = pipe->pastFlow;
        }

        network->msgLog << "\n  [Indicators] Initialized prevSecondState from past values at t="
            << currentTime;
    }
    else if (!prevSecondState.isValid) {
        // dt < 1.0 and no valid history - can't compute meaningful indicators
        maxPhiA = 0.0;
        maxPhiR = 0.0;
        std::fill(phiA.begin(), phiA.end(), 0.0);
        std::fill(phiR.begin(), phiR.end(), 0.0);
        logIndicators(currentTime, dt, false);
        return;
    }

    // ------------------------------------------------------------------------
    // Compute indicators
    // ------------------------------------------------------------------------
    double timeDiff;
    timeDiff = currentSecondState.time - prevSecondState.time;

    if (timeDiff < dt)
    {
        timeDiff = dt;
    }
  
    maxPhiA = 0.0;
    maxPhiR = 0.0;

    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);

        if (!pipe) {
            if (link && link->index >= 0 && link->index < linkCount) {
                phiA[link->index] = 0.0;
                phiR[link->index] = 0.0;
            }
            continue;
        }

        int idx = pipe->index;
        if (idx < 0 || idx >= linkCount) continue;

        double L = pipe->length;
        double D = pipe->diameter;
        double A = PI * D * D / 4.0;
        double a = pipe->getWaveSpeed();

        if (L < MIN_DENOM || A < MIN_DENOM || a < MIN_DENOM) {
            phiA[idx] = 0.0;
            phiR[idx] = 0.0;
            continue;
        }

        // ====================================================================
        // COMPRESSIBILITY TERM
        // ====================================================================
        int n1 = pipe->fromNode->index;
        int n2 = pipe->toNode->index;

        double H1_current = currentSecondState.nodeHeads[n1];
        double H2_current = currentSecondState.nodeHeads[n2];
        double H1_prev = prevSecondState.nodeHeads[n1];
        double H2_prev = prevSecondState.nodeHeads[n2];

        double dH1_dt = (H1_current - H1_prev) / timeDiff;
        double dH2_dt = (H2_current - H2_prev) / timeDiff;

        double compressibilityTerm = (L / (2.0 * a)) * std::abs(dH1_dt + dH2_dt);

        // ====================================================================
        // INERTIA TERM
        // ====================================================================
        double DQ_Dt;

        double QA_current = currentSecondState.pipeStartFlows[idx];
        double QB_current = currentSecondState.pipeEndFlows[idx];
        double Q_current = currentSecondState.pipeFlows[idx];

        double QA_prev = prevSecondState.pipeStartFlows[idx];
        double QB_prev = prevSecondState.pipeEndFlows[idx];
        double Q_prev = prevSecondState.pipeFlows[idx];

        if (currentMode == WATER_HAMMER) {
            double dQA_dt = (QA_current - QA_prev) / timeDiff;
            double dQB_dt = (QB_current - QB_prev) / timeDiff;
            DQ_Dt = 0.5 * (std::abs(dQA_dt) + std::abs(dQB_dt));
        }
        else {
            DQ_Dt = std::abs((Q_current - Q_prev) / timeDiff);
        }

        double inertiaTerm = (L / (g * A)) * DQ_Dt;

        // ====================================================================
        // FRICTION LOSS F(Q) - COMPUTED DIRECTLY FROM PIPE PROPERTIES
        // 
        // This is critical: we must NOT use pipe->hLoss because it may have
        // different meanings in different solver modes. Instead, we compute
        // the steady-state friction loss directly using the resistance
        // coefficient and current flow.
        // ====================================================================
        double Q_avg;
        if (currentMode == WATER_HAMMER) {
            // Use average of start and end flows
            Q_avg = 0.5 * (std::abs(QA_current) + std::abs(QB_current));
        }
        else {
            Q_avg = std::abs(Q_current);
        }

        double frictionLoss_feet = 0.0;

        // Get the head loss model type from network options
        std::string headlossModel = network->option(Options::HEADLOSS_MODEL);

        if (headlossModel == "H-W") {
            // Hazen-Williams: hf = K * Q^1.852
            // pipe->resistance stores the H-W resistance coefficient
            double K = pipe->resistance;
            if (K > MIN_DENOM && Q_avg > MIN_DENOM) {
                frictionLoss_feet = K * std::pow(Q_avg, 1.852);
            }
        }
        else {
            // Darcy-Weisbach: hf = K * Q^2  or  hf = f * L / (2*g*D*A^2) * Q^2
            // pipe->resistance stores the D-W resistance coefficient
            double K = pipe->resistance;
            if (K > MIN_DENOM && Q_avg > MIN_DENOM) {
                frictionLoss_feet = K * Q_avg * Q_avg;
            }
            else {
                // Fallback: compute from friction factor
                double f = pipe->roughness;  // This might be Darcy friction factor
                if (f > 0 && D > 0 && A > 0) {
                    double K_calc = f * L / (2.0 * g * D * A * A);
                    frictionLoss_feet = K_calc * Q_avg * Q_avg;
                }
            }
        }

        // Ensure minimum friction loss to prevent φ_R = 1.0 exactly
        // Even a small flow produces some friction loss
        double minFrictionLoss_feet = 0.01;  // 0.01 ft minimum (~3mm)
        frictionLoss_feet = std::max(frictionLoss_feet, minFrictionLoss_feet);

        // ====================================================================
        // FINAL INDICATORS (convert to meters)
        // ====================================================================
        double phi_A_feet = compressibilityTerm + inertiaTerm;
        double phi_A = phi_A_feet * FT_TO_M;

        double frictionLoss = frictionLoss_feet * FT_TO_M;

        double phi_R = 0.0;
        double denominator = phi_A + frictionLoss;
        if (denominator > MIN_DENOM) {
            phi_R = phi_A / denominator;
        }

        // Store computed indicators
        phiA[idx] = phi_A;
        phiR[idx] = phi_R;
        pipe->phiA = phi_A;
        pipe->phiR = phi_R;

        if (phi_A > maxPhiA) maxPhiA = phi_A;
        if (phi_R > maxPhiR) maxPhiR = phi_R;
    }

    evaluateModeThisStep = true;
    logIndicators(currentTime, dt, true);
}


void CGGASolver::logIndicators(double currentTime, double dt, bool modeEval) {
    // ========================================================================
    // Logging helper - only logs at integer seconds
    // ========================================================================
    if (!indicatorLog.is_open()) {
        indicatorLog.open("C:\\EPANET-CGGA\\Networks\\Unsteadiness_Indicators.txt");
        indicatorLog << "Time_s" << "\t"
            << "Mode" << "\t"
            << "Max_phiA_m" << "\t"
            << "Max_phiR" << "\t"
            << "EvalMode" << std::endl;
    }

    indicatorLog << static_cast<int>(std::round(currentTime)) << "\t"
        << static_cast<int>(currentMode) << "\t"
        << maxPhiA << "\t"
        << maxPhiR << "\t"
        << (modeEval ? 1 : 0) << std::endl;
}

bool CGGASolver::detectHydraulicEvent(double currentTime) {
   
    const double settingTolerance = 1e-6;
    const double demandChangeTolerance = 0.01;  // 1% change threshold

    // Check valves
    for (Link* link : network->links) {
        if (link->type() == Link::VALVE) {
            Valve* valve = static_cast<Valve*>(link);

            // Check if setting is changing
            if (valve->settingPattern != nullptr) {
                double currentSetting = valve->setting;
                double previousSetting = valve->pastSetting;

                if (std::abs(currentSetting - previousSetting) > settingTolerance) {
                    network->msgLog << "\n  Event: Valve " << valve->index
                        << " setting change detected"
                        << " (" << previousSetting << " → " << currentSetting << ")";
                    return true;
                }
            }
        }
    }

    // Check pumps
    for (Link* link : network->links) {
        if (link->type() == Link::PUMP) {
            Pump* pump = static_cast<Pump*>(link);

            // Check status change
            if (pump->status != pump->previousStatus) {
                network->msgLog << "\n  Event: Pump " << pump->index
                    << " status change detected";
                return true;
            }

            // Check speed change
            if (pump->speedPattern != nullptr) {
                double currentSpeed = pump->speed;
                double previousSpeed = pump->pastSpeed;

                if (std::abs(currentSpeed - previousSpeed) > settingTolerance) {
                    network->msgLog << "\n  Event: Pump " << pump->index
                        << " speed change detected";
                    return true;
                }
            }
        }
    }

    // Check for significant demand changes
    /*
    for (Node* node : network->nodes) {
        if (node->type() == Node::JUNCTION) {
            double currentDemand = node->actualDemand;
            double previousDemand = node->pastDemand;
            double avgDemand = (currentDemand + previousDemand) / 2.0;

            if (avgDemand > 0 &&
                std::abs(currentDemand - previousDemand) / avgDemand > demandChangeTolerance) {
                return true;
            }
        }
    }
    */

    return false;
}

int CGGASolver::executeQuasiSteadySolve(double dt, int& trials, double currentTime) {

    network->msgLog << "\n--- Quasi-Steady Solve (t=" << currentTime << " s) ---";

    bool converged = false;
    bool statusChanged = true;
    double lambda = 1.0;

    setConvergenceLimits();
    trials = 1;

    while (trials <= trialsLimit) {
        oldErrorNorm = errorNorm;
        setFixedGradeNodes();

        if (statusChanged) {
            oldErrorNorm = findErrorNorm(0.0, currentTime, dt, currentMode);
            lambda = 1.0;
        }
        statusChanged = false;

        int errorCode = findSteadyHeadChanges();
        if (errorCode >= 0) {
            Node* node = network->node(errorCode);
            network->msgLog << "\n  QS ill-conditioned at node: " << node->name;
            return HydSolver::FAILED_ILL_CONDITIONED;
        }

        findSteadyFlowChanges();

        errorNorm = findStepSize(trials, currentTime);
        updateSolution(1.0, false);

        if (reportTrials) reportTrial(trials, lambda);

        converged = hasConverged();
        if (converged) {
            statusChanged = linksChangedStatus();
        }

        if (converged && !statusChanged) break;
        trials++;
    }

    network->msgLog << "\n  QS: " << trials << " iterations, error=" << errorNorm;

    return converged ? HydSolver::SUCCESSFUL : HydSolver::FAILED_NO_CONVERGENCE;
}


int CGGASolver::executeRWCSolve(double tstep_advance, int& trials, double currentTime) {
    
    network->msgLog << "\n--- RWC Solve (t=" << currentTime << " s) ---";
    network->msgLog << "\n    dt_advance=" << tstep_advance << "s";
    network->msgLog << "\n    dt_calc=" << previousTimestep << "s (from "
        << getModeString(previousMode) << ")";

    // Time advancement
    tstep = tstep_advance;

    // ADAPTIVE: Use previous mode's timestep for ALL RWC calculations
    double dt_calc = previousTimestep;

    bool converged = false;
    bool statusChanged = true;
    double lambda = 1.0;
    const double rwcTheta = 1.0;

    setConvergenceLimits();
    trials = 1;

    while (trials <= trialsLimit) {
        oldErrorNorm = errorNorm;
        setFixedGradeNodes();

        if (statusChanged) {
            oldErrorNorm = findRWCErrorNorm(currentTime, tstep_advance, rwcTheta);
            lambda = 1.0;
        }
        statusChanged = false;

        int errorCode = findRWCHeadChanges(currentTime, tstep_advance, rwcTheta);
        if (errorCode >= 0) {
            Node* node = network->node(errorCode);
            network->msgLog << "\n  RWC ill-conditioned at node: " << node->name;
            return HydSolver::FAILED_ILL_CONDITIONED;
        }

        findRWCFlowChanges(tstep_advance, rwcTheta);
        errorNorm = findRWCStepSize(trials, currentTime, tstep_advance, rwcTheta);
        updateRWCSolution(1.0);

        if (reportTrials) reportTrial(trials, lambda);

        converged = hasConverged();
        if (converged) {
            statusChanged = linksChangedStatus();
        }
        if (converged && !statusChanged) break;

        trials++;
    }

    network->msgLog << "\n  RWC: " << trials << " iterations, error=" << errorNorm;
    network->msgLog << "\n  Calculated with dt_calc=" << dt_calc << "s, advancing by " << tstep_advance << "s";

    return converged ? HydSolver::SUCCESSFUL : HydSolver::FAILED_NO_CONVERGENCE;
}

//  Find changes in nodal heads by solving a linearized system of equations.
int CGGASolver::findRWCHeadChanges(double currentTime, double tstep, double rwcTheta)
{
    // ... setup the coeff. matrix of the RWCGGA linearized system

    setRWCMatrixCoeffs(tstep);

    // ... temporarily use the head change array dH[] to store new heads

    double* h = &dH[0];

    // ... solve the linearized RWCGGA system for new nodal heads
    //     (matrixSolver returns a negative integer if it runs successfully;
    //      otherwise it returns the index of the row that caused it to fail.)

    int errorCode = matrixSolver->solve(nodeCount, h);
    if (errorCode >= 0) return errorCode;

    // ... save new heads as head changes

    for (int i = 0; i < nodeCount; i++)
    {
        dH[i] = h[i] - network->node(i)->head;
    }

    // ... return a negative number indicating that
    //     the matrix solver ran successfully

    return -1;
}

void CGGASolver::setRWCMatrixCoeffs(double tstep)
{
    memset(&xQ[0], 0, nodeCount * sizeof(double));
    matrixSolver->reset();
    setRWCLinkCoeffs(tstep);
    setRWCNodeCoeffs(tstep);
    setRWCValveCoeffs();
}

//  Compute matrix coefficients for link head loss gradients.
void CGGASolver::setRWCLinkCoeffs(double tstep)
{
    for (int j = 0; j < linkCount; j++)
    {
        // ... skip links with zero head gradient
        //     (e.g. active pressure regulating valves)

        Link* link = network->link(j);
        if (link->hGrad == 0.0) continue;

        // ... identify end nodes of link

        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        int n1 = node1->index;
        int n2 = node2->index;

        // ... update node flow balances

        xQ[n1] -= link->flow;
        xQ[n2] += link->flow;

        // ... a is contribution to coefficient matrix
        //     b is contribution to right hand side

        if (tstep == 0)
        {
            double a = 1.0 / link->hGrad;
            double b = a * link->hLoss;

            // ... update off-diagonal coeff. of matrix if both start and
            //     end nodes are not fixed grade

            if (!node1->fixedGrade && !node2->fixedGrade)
            {
                matrixSolver->addToOffDiag(j, -a);
            }

            // ... if start node has fixed grade, then apply a to r.h.s.
            //     of that node's row;

            if (node1->fixedGrade)
            {
                matrixSolver->addToRhs(n2, a * node1->head);
            }

            // ... otherwise add a to row's diagonal coeff. and
            //     add b to its r.h.s.

            else
            {
                matrixSolver->addToDiag(n1, a);
                matrixSolver->addToRhs(n1, b);
            }

            // ... do the same for the end node, except subtract b from r.h.s

            if (node2->fixedGrade)
            {
                matrixSolver->addToRhs(n1, a * node2->head);
            }
            else
            {
                matrixSolver->addToDiag(n2, a);
                matrixSolver->addToRhs(n2, -b);
            }
        }

        else
        {
            // a and b are upgraded according to RWC-GGA

            double a = 1.0 / ((link->hGrad) + (link->inertialTerm / (kappa * tstep)));
            double b = a * ((link->hLoss) + ((1 - kappa) / kappa) * (link->pastHloss - (node1->pastHead - node2->pastHead)) - (link->inertialTerm / (kappa * tstep)) * (link->pastFlow) - link->hGrad * link->flow) + link->flow;

            // ... update off-diagonal coeff. of matrix if both start and
            //     end nodes are not fixed grade

            if (!node1->fixedGrade && !node2->fixedGrade)
            {
                matrixSolver->addToOffDiag(j, -a);
            }

            // ... if start node has fixed grade, then apply a to r.h.s. of that node's row;

            if (node1->fixedGrade)
            {
                matrixSolver->addToRhs(n2, a * node1->head);
            }

            // ... otherwise add a to row's diagonal coeff. and add b to its r.h.s.

            else
            {
                matrixSolver->addToDiag(n1, a);
                matrixSolver->addToRhs(n1, b);
            }

            // ... do the same for the end node, except subtract b from r.h.s

            if (node2->fixedGrade)
            {
                matrixSolver->addToRhs(n1, a * node2->head);
            }
            else
            {
                matrixSolver->addToDiag(n2, a);
                matrixSolver->addToRhs(n2, -b);
            }
        }
    }
}

//  Compute matrix coefficients for dynamic tanks and external node outflows.
void  CGGASolver::setRWCNodeCoeffs(double tstep)
{
    for (int i = 0; i < nodeCount; i++)
    {
        // ... if node's head not fixed

        Node* node = network->node(i);
        if (!node->fixedGrade)
        {
            // ... for dynamic tanks, add area terms to row i
            //     of the head solution matrix & r.h.s. vector

            if (node->type() == Node::TANK && theta != 0.0)
            {
                Tank* tank = static_cast<Tank*>(node);

                if (tank->head == tank->pastHead)
                {
                    double a = tank->area / (theta * tstep);
                    matrixSolver->addToDiag(i, a);

                    a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
                    matrixSolver->addToRhs(i, a); //  */
                }

                else
                {
                    double a = (tank->area + ((tank->area - tank->pastArea) / (tank->head - tank->pastHead)) * (tank->head)) / (theta * tstep);
                    matrixSolver->addToDiag(i, a);

                    double b = (tank->pastArea) * tank->pastHead / (theta * tstep) + (1.0 - theta) * tank->pastOutflow / theta;
                    double c = b + ((tank->area - tank->pastArea) / (tank->head - tank->pastHead)) * (tank->head) / (theta * tstep);
                    matrixSolver->addToRhs(i, c); //
                }
            }

            // ... for junctions, add effect of external outflows

            else if (node->type() == Node::JUNCTION)
            {
                // ... update junction's net inflow
                xQ[i] -= node->outflow;
                matrixSolver->addToDiag(i, node->qGrad);
                matrixSolver->addToRhs(i, node->qGrad * node->head);
            }

            // ... add node's net inflow to r.h.s. row
            matrixSolver->addToRhs(i, (double)xQ[i]);
        }

        // ... if node has fixed head, force solution to produce it

        else
        {
            matrixSolver->setDiag(i, 1.0);
            matrixSolver->setRhs(i, node->head);
        }
    }
}

//-----------------------------------------------------------------------------

//  Compute matrix coefficients for pressure regulating valves.

void  CGGASolver::setRWCValveCoeffs()
{
    for (Link* link : network->links)
    {
        // ... skip links that are not active pressure regulating valves

        if (link->hGrad > 0.0) continue;

        // ... determine end node indexes of link

        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... add net inflow of downstream node of a PRV to the
        //     r.h.s. row of its upstream node

        if (link->isPRV())
        {
            matrixSolver->addToRhs(n1, (double)xQ[n2]);
        }

        // ... add net inflow of upstream node of a PSV to the
        //     r.h.s. row of its downstream node

        if (link->isPSV())
        {
            matrixSolver->addToRhs(n2, (double)xQ[n1]);
        }
    }
}

//-----------------------------------------------------------------------------

//  Find the changes in link flows resulting from a set of nodal head changes.

void CGGASolver::findRWCFlowChanges(double tstep, double rwcTheta)
{
    for (int i = 0; i < linkCount; i++)
    {
        // ... get link object and its end node indexes

        dQ[i] = 0.0;
        Link* link = network->link(i);
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... flow change for pressure regulating valves

        if (link->hGrad == 0.0)
        {
            if (link->isPRV()) dQ[i] = -xQ[n2] - link->flow;
            if (link->isPSV()) dQ[i] = xQ[n1] - link->flow;
            continue;
        }

        if (tstep == 0) // || network->link->type() == Link::VALVE || network->link->type() == Link::PUMP)
        {

            // ... apply GGA flow change formula:

            double dh = (link->fromNode->head + dH[n1]) -
                (link->toNode->head + dH[n2]);
            double dq = (link->hLoss - dh) / link->hGrad;

            // ... special case to prevent negative flow in constant HP pumps

            if (link->isHpPump() &&
                link->status == Link::LINK_OPEN &&
                dq > link->flow) dq = link->flow / 2.0;
            // ... save flow change

            dQ[i] = -dq;
        }
        else
        {
            // ... apply RWCGGA flow change formula:

            double dh = (link->fromNode->head + dH[n1]) - (link->toNode->head + dH[n2]);
            double dhpast = (link->fromNode->pastHead) - (link->toNode->pastHead);
            double pastTerms = ((1 - kappa) / kappa) * (link->pastHloss - dhpast);
            //double flows = (link->inertialTerm / (kappa * tstep)) * (link->flow - link->pastFlow);
            double flows = (link->inertialTerm / (kappa * tstep)) * (link->pastFlow) + link->hGrad * link->flow;
            double dq = -(dh - link->hLoss + flows - pastTerms) / (link->hGrad + (link->inertialTerm / (kappa * tstep))) +link->flow;

            // ... special case to prevent negative flow in constant HP pumps

            if (link->isHpPump() &&
                link->status == Link::LINK_OPEN &&
                dq > link->flow) dq = link->flow / 2.0;

            // ... save flow change

            dQ[i] = -dq;
        }
    }
}

//-----------------------------------------------------------------------------

//  Find how much of the head and flow changes to apply to a new solution.

double CGGASolver::findRWCStepSize(int trials, double currentTime, double tstep, double rwcTheta)
{
    // ... find the new error norm at full step size

    double lamda = 1.0;
    errorNorm = findErrorNorm(lamda, currentTime, tstep, currentMode);

    if (stepSizing == RELAXATION && oldErrorNorm < ErrorThreshold)
    {
        lamda = 0.5;
        double errorNorm2 = findErrorNorm(lamda, currentTime, tstep, currentMode);
        if (errorNorm2 < errorNorm) errorNorm = errorNorm2;
        else
        {
            lamda = 1.0;
            errorNorm = findErrorNorm(lamda, currentTime, tstep, currentMode);
        }
    }

    // ... if called for, implement a lamda search procedure
    //     to find the best step size lamda to take

    if (stepSizing == ARF || stepSizing == BRF)
    {
        {
            minErrorNorm = 0;
            double testError = 1000000;
            lambdaNumber = 1 / dl;
            Lambda.resize(lambdaNumber, 0);

            memset(&Lambda[0], 0, lambdaNumber * sizeof(double));

            for (int i = 0; i < lambdaNumber; i++)
            {
                Lambda[i] += (i + 1) * dl;

                int errorCode = findRWCHeadChanges(currentTime, tstep, rwcTheta);
                if (errorCode >= 0)
                {
                    Node* node = network->node(errorCode);
                    network->msgLog << endl << s_IllConditioned << node->name;
                    return HydSolver::FAILED_ILL_CONDITIONED;
                }
                findRWCFlowChanges(tstep, rwcTheta);                  //*/
                errorNorm = findErrorNorm(Lambda[i], currentTime, tstep, currentMode);

                if (errorNorm < testError)
                {
                    testError = errorNorm;
                    lamda = Lambda[i];
                    updateSolution(Lambda[i]);
                }
            }
            minErrorNorm = testError;
        }
        return lamda;
    }
    return errorNorm;
}

//-----------------------------------------------------------------------------

//  Compute the error norm associated with a given step size.

double CGGASolver::findRWCErrorNorm(double lamda, double currentTime, double tstep)
{

    hLossEvalCount++;
    return hydBalance.evaluate(lamda, (double*)&dH[0], (double*)&dQ[0],
        (double*)&xQ[0], network, currentTime, tstep, currentMode);
}

//-----------------------------------------------------------------------------

//  Update heads and flows for a given step size.

void CGGASolver::updateRWCSolution(double lamda)
{
    for (int i = 0; i < nodeCount; i++)
    {
        network->node(i)->head += lamda * dH[i];
    }
    for (int i = 0; i < linkCount; i++)
    {
        network->link(i)->flow += lamda * dQ[i];

		Pipe* pipe = dynamic_cast<Pipe*>(network->link(i));

        if (pipe)
        {
            pipe->startFlow += lamda * dQ[i];
            pipe->endFlow += lamda * dQ[i];
        }

    }
}


int CGGASolver::executeWaterHammerSolve(double dt, int& trials, double currentTime) {

    network->msgLog << "\n--- Water Hammer Solve (t=" << currentTime << " s) ---";

    // Configure pipes for this time step
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe) {
            configurePipeForTimestep(pipe, dt, network);
        }
    }

    bool converged = false;
    bool statusChanged = true;
    double lambda = 1.0;
    double prevErrorNorm = HUGE_VAL;
    double dl = 1.0;

    const int maxIterations = 33;
    int iterations = 0;

    setConvergenceLimits();

    while (!converged && iterations < maxIterations) {

        prevErrorNorm = errorNorm;
        setFixedGradeNodes();

        if (statusChanged) {
            prevErrorNorm = findUnsteadyErrorNorm(1.0, currentTime, dt);
            lambda = 1.0;
        }
        statusChanged = false;
        dl = 1.0;

        // Store state for line search
        storeIterationState();

        int errorCode = findUnsteadyHeadChanges(currentTime);
        if (errorCode >= 0) {
            network->msgLog << "\n  WH ill-conditioned at position " << errorCode;
            return HydSolver::FAILED_ILL_CONDITIONED;
        }

        findUnsteadyFlowChanges(currentMode, currentTime);
    
        // Step sizing with line search
        if (stepSizing == ARF) {
            lambda = 1.0;
            double trialError = calculateTrialError(lambda, currentTime, dt);

            if (trialError < prevErrorNorm * 0.99) {
                errorNorm = trialError;
                restoreIterationState();
                updateWHSolution(lambda);
            }
            else {
                while (minErrorNorm >= prevErrorNorm * 0.99 && dl > 0.001) {
                    dl *= 0.25;
                    lambda = findUnsteadyStepSize(iterations, currentTime, dt, prevErrorNorm, dl);
                }
                errorNorm = minErrorNorm;
            }
        }
        else if (stepSizing == BRF) {
            while (minErrorNorm >= prevErrorNorm && dl > 0.001) {
                dl *= 0.25;
                lambda = findUnsteadyStepSize(iterations, currentTime, dt, prevErrorNorm, dl);
            }
            errorNorm = minErrorNorm;
        }
        else {
            errorNorm = findUnsteadyErrorNorm(1.0, currentTime, dt);
            updateWHSolution(1.0);
        }

        converged = hasConverged();
        if (converged) {
            statusChanged = linksChangedStatus();
        }

        if (converged && !statusChanged) break;
        iterations++;
    
    }

    trials = iterations;

    network->msgLog << "\n  WH: " << iterations << " iterations, error=" << errorNorm;

    if (converged) {
        return HydSolver::SUCCESSFUL;
    }
    else if (errorNorm < 0.1) {
        network->msgLog << "\n  WH accepting solution (error < 0.1)";
        return HydSolver::SUCCESSFUL;
    }
    else {
        return HydSolver::FAILED_NO_CONVERGENCE;
    }
}

void CGGASolver::storeIterationState() {
    iterHeads.resize(nodeCount);
    iterFlows.resize(linkCount);
    iterStartFlows.resize(linkCount);
    iterEndFlows.resize(linkCount);

    for (int i = 0; i < nodeCount; i++) {
        Node* node = network->node(i);
        if (node) iterHeads[i] = node->head;
    }

    for (int i = 0; i < linkCount; i++) {
        Link* link = network->link(i);
        if (link) {
            iterFlows[i] = link->flow;
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                iterStartFlows[i] = pipe->startFlow;
                iterEndFlows[i] = pipe->endFlow;
            }
        }
    }
}

void CGGASolver::restoreIterationState() {
    for (int i = 0; i < nodeCount && i < iterHeads.size(); i++) {
        Node* node = network->node(i);
        if (node) node->head = iterHeads[i];
    }

    for (int i = 0; i < linkCount && i < iterFlows.size(); i++) {
        Link* link = network->link(i);
        if (link) {
            link->flow = iterFlows[i];
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                pipe->startFlow = iterStartFlows[i];
                pipe->endFlow = iterEndFlows[i];
            }
        }
    }
}

std::string CGGASolver::getModeString(SolverMode mode) const {
    switch (mode) {
    case QUASI_STEADY: return "QUASI_STEADY";
    case RWC: return "RWC";
    case WATER_HAMMER: return "WATER_HAMMER";
    default: return "UNKNOWN";
    }
}

double CGGASolver::findUnsteadyStepSize(int trials, double currentTime, double tstep, double prevErrorNorm, double dl) {
    // This function searches for the optimal relaxation factor
    // It's called when we need to find a lambda that reduces the error
    
    // First, always check the full step
    double lambda = 1.0;
    double fullStepError = findUnsteadyErrorNorm(lambda, currentTime, tstep);
    
    // For the relaxation method, try a half step near convergence
    if (stepSizing == RELAXATION && prevErrorNorm < 1e-5) {
        lambda = 0.5;
        double halfStepError = findUnsteadyErrorNorm(lambda, currentTime, tstep);
        
        if (halfStepError < fullStepError) {
            errorNorm = halfStepError;
            network->msgLog << "\n    Half-step relaxation applied";
        } else {
            lambda = 1.0;
            errorNorm = fullStepError;
        }
        
        // Apply the chosen lambda
        restoreState();
        updateWHSolution(lambda);
        return lambda;
    }
    
    // For ARF and BRF methods, perform a comprehensive search
    if (stepSizing == ARF || stepSizing == BRF) {
        // Initialize search parameters
        minErrorNorm = std::numeric_limits<double>::max();
        double bestLambda = 1.0;
        double bestError = fullStepError;
        
        // Calculate number of lambda values to test based on dl
        int lambdaCount = static_cast<int>(1.0 / dl);
        if (lambdaCount < 2) lambdaCount = 2;  // At least test 0.5 and 1.0
        if (lambdaCount > 100) lambdaCount = 100;  // Cap for efficiency
        
        std::vector<double> lambdaValues(lambdaCount);
        
        // Generate lambda values to test
        for (int i = 0; i < lambdaCount; i++) {
            lambdaValues[i] = (i + 1) * dl;
            if (lambdaValues[i] > 1.0) lambdaValues[i] = 1.0;
        }
        
        network->msgLog << "\n    Testing " << lambdaCount << " lambda values with dl=" << dl;
        
        // Test each lambda value
        for (int i = 0; i < lambdaCount; i++) {
            // Restore to the state before any updates
            restoreState();
            
            // Apply this lambda and calculate error
            double testLambda = lambdaValues[i];
            updateWHSolution(testLambda);
            double testError = findUnsteadyErrorNorm(1.0, currentTime, tstep);
            
            // Track the best lambda
            if (testError < bestError) {
                bestError = testError;
                bestLambda = testLambda;
                network->msgLog << "\n      λ=" << testLambda << " gives error=" << testError << " (best so far)";
            }
        }
        
        // Set the minimum error found for the calling function
        minErrorNorm = bestError;
        
        // Apply the best lambda found
        restoreState();
        updateWHSolution(bestLambda);
        errorNorm = bestError;
        
        network->msgLog << "\n    Selected λ=" << bestLambda << " with error=" << bestError;
        
        return bestLambda;
    }
    
    // Default: return full step
    restoreState();
    updateWHSolution(1.0);
    errorNorm = fullStepError;
    return 1.0;
}

// Helper function to calculate error for a trial λ without permanently updating
double CGGASolver::calculateTrialError(double lambda, double currentTime, double tstep) {
    // This function tests what the error would be with a given lambda
    // without permanently changing the state
    
    // First, restore to the state before any updates
    restoreState();
    
    // Apply the trial relaxation factor
    updateWHSolution(lambda);
    
    // Calculate and return the resulting error norm
    // This evaluates how well the equations are satisfied with this lambda
    double trialError = findUnsteadyErrorNorm(1.0, currentTime, tstep);
    
    return trialError;
}

void CGGASolver::restoreState() {
    // Restore node heads - heads are stored in the node objects
    for (int i = 0; i < nodeCount && i < saved_nodeHeads.size(); i++) {
        Node* node = network->nodes[i];
        if (node) {
            node->head = saved_nodeHeads[i];  // Restore the head to the node object
        }
    }
    
    // Restore link flows - flows are stored in the link objects
    for (int i = 0; i < linkCount && i < saved_linkFlows.size(); i++) {
        Link* link = network->links[i];
        if (!link) continue;
        
        link->flow = saved_linkFlows[i];  // Restore the flow to the link object
        
        // For water hammer pipes, also restore start and end flows
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe) {
            pipe->startFlow = saved_startFlows[i];
            pipe->endFlow = saved_endFlows[i];
        } 
    }
    
    // Restore the computed changes (dH, dQ, etc.)
    // These ARE simple arrays in your solver, so we can restore them directly
    for (int i = 0; i < saved_dH.size() && i < dH.size(); i++) {
        dH[i] = saved_dH[i];
    }
    
    for (int i = 0; i < saved_dQ.size() && i < dQ.size(); i++) {
        dQ[i] = saved_dQ[i];
    }
    
    for (int i = 0; i < saved_dQ_start.size() && i < dQ_start.size(); i++) {
        dQ_start[i] = saved_dQ_start[i];
    }
    
    for (int i = 0; i < saved_dQ_end.size() && i < dQ_end.size(); i++) {
        dQ_end[i] = saved_dQ_end[i];
    }
}

//  Compute the error norm associated with a given step size for unsteady analysis.
double CGGASolver::findUnsteadyErrorNorm(double lamda, double currentTime, double tstep)
{
    hLossEvalCount++;

    // Pass transition info to HydBalance
    return hydBalance.evaluateUnsteady(
        lamda, &dH[0], &dQ[0], &dQ_start[0], &dQ_end[0], &xQ[0], network, currentTime, tstep);
}

void CGGASolver::setConvergenceLimits()
{
    // ... maximum trials
    trialsLimit = network->option(Options::MAX_TRIALS);

    // ... limit on relative flow change ratio (old EPANET2 criterion)
    flowRatioLimit = network->option(Options::RELATIVE_ACCURACY);

    // ... tolerance in satisfying hydraulic balances
    headErrLimit = network->option(Options::HEAD_TOLERANCE) /
        network->ucf(Units::LENGTH);
    flowErrLimit = network->option(Options::FLOW_TOLERANCE) / network->ucf(Units::FLOW);

    // ... limit on largest flow change
    flowChangeLimit = network->option(Options::FLOW_CHANGE_LIMIT) /
        network->ucf(Units::FLOW);

    // ... use a default head error limit if need be
    if (flowRatioLimit == 0.0 && headErrLimit == 0.0 &&
        flowErrLimit == 0.0 && flowChangeLimit == 0.0)
    {
        headErrLimit = 0.005;
    }

    // ... convert missing limits to a huge number
    if (flowRatioLimit == 0.0) flowRatioLimit = Huge;
    if (headErrLimit == 0.0) headErrLimit = Huge;
    if (flowErrLimit == 0.0) flowErrLimit = Huge;
    if (flowChangeLimit == 0.0) flowChangeLimit = Huge;
}

int CGGASolver::findUnsteadyHeadChanges(double currentTime)
{
    // ... setup the coeff. matrix of the GGA linearized system

    setUnsteadyMatrixCoeffs(currentTime);

    // ... temporarily use the head change array dH[] to store new heads

    double* h = &dH[0];

    // ... solve the linearized GGA system for new nodal heads
    //     (matrixSolver returns a negative integer if it runs successfully;
    //      otherwise it returns the index of the row that caused it to fail.)

    int errorCode = matrixSolver->solve(nodeCount, h);
    if (errorCode >= 0) return errorCode;

    // ... save new heads as head changes

    for (int i = 0; i < nodeCount; i++)
    {
        dH[i] = h[i] - network->nodes[i]->head;
    }

    // ... return a negative number indicating that
    //     the matrix solver ran successfully

    return -1;
}

void CGGASolver::findUnsteadyFlowChanges(int currentMode, double currentTime) {
    // Main loop for all links in the network
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;
        
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;
        
        // Reset flow changes
        dQ_start[i] = 0.0;
        dQ_end[i] = 0.0;
        dQ[i] = 0.0;
        
        // Get nodal heads (including changes from linear system solution)
        double H_A = link->fromNode->head + dH[n1];
        double H_B = link->toNode->head + dH[n2];

        // COMPRESSIBLE FLOW MODEL (WATER HAMMER)
        if (link->type() == Link::PIPE) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue; // Skip if not a pipe
            
            // === CALCULATE CHARACTERISTIC PARAMETERS ===
            
            // Physical constants
            double c = pipe->getWaveSpeed();
            double A = pipe->getArea();
            double L = pipe->length;
            double g = GRAVITY;
            double Bj = c / (g * A);  // Wave impedance term (B')
            double D = pipe->diameter;
            double epsilon = 0.5;  // Friction integration parameter
            
            // Calculate wave travel time and reach-back time
            double physicalWaveTime = L / c;
            double reachBackTime;
            
            // Determine reach-back index based on number of reaches
            int reachBackSteps = pipe->numReaches;
            double waveTravel = reachBackSteps * tstep;  // Wave travel distance

            // Retrieve reach-back flows and heads from history
            double Qa = pipe->pastStartFlow;  // Default fallback
            double Qb = pipe->pastEndFlow;    // Default fallback
            double Ha = pipe->fromNode->pastHead;  // Default fallback
            double Hb = pipe->toNode->pastHead;    // Default fallback
            
            // Use the FlowHistoryManager to get reach-back values
            FlowHistoryResult historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, network);
            //FlowHistoryResult historyResult = getReachBackValues(pipe, currentTime, physicalWaveTime, network);

            // If history was found, use the values from the result
            if (historyResult.found) {
                Qa = historyResult.startFlow;
                Qb = historyResult.endFlow;
                Ha = historyResult.startHead;
                Hb = historyResult.endHead;
            }
    
            // === FRICTION AND RESISTANCE CALCULATIONS (CORRECTED) ===
            
            double KA, KB;

            if (network->option(Options::HEADLOSS_MODEL) == "D-W")
            {
                double frictionFactor = pipe->getFrictionFactor(pipe->getRe(pipe->flow, network->option(Options::KIN_VISCOSITY)));
                double K = frictionFactor * pipe->length / (2 * g * D * A * A);

                double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, network->option(Options::KIN_VISCOSITY)));
                double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, network->option(Options::KIN_VISCOSITY)));
                // Split K for each end as needed (KA for upstream, KB for downstream)
                KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
                KB = frictionFactorB * pipe->length / (2 * g * D * A * A);

            }
            else
            {
                KA = pipe->resistance;
                KB = pipe->resistance;
            }
            
            const double MIN_RESISTANCE = 1e-10;
            KA = std::max(KA, MIN_RESISTANCE);
            KB = std::max(KB, MIN_RESISTANCE);
            
            // Unsteady friction terms
            double kUFa = pipe->computeShearDecayCoeff(Qa);
            double kUFb = pipe->computeShearDecayCoeff(Qb);
            double Ua = 0.0;
            double Ub = 0.0;
            
            if (kUFa > 0.0 && kUFb > 0.0) {
                auto sign = [](double val) -> double { 
                    return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0); 
                };
                
                // CORRECTED: Unsteady friction terms
                Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
                Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
            }
            
            // CORRECTED: Characteristic equation coefficients
            double BA, BB, RA, RB;
            
            if (network->option(Options::HEADLOSS_MODEL) == "D-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
            }
            else if (network->option(Options::HEADLOSS_MODEL) == "H-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;
            }
            
            pipe->charB_plus = BA;
            pipe->charB_minus = BB;
            pipe->charC_plus = RA;
            pipe->charC_minus = RB;
            
            // === CALCULATE NEW FLOWS ===
            
            double newStartFlow = (H_A + RA) / BA;
            double newEndFlow = (-H_B + RB) / BB;

            double avgFlow = (newStartFlow + newEndFlow) / 2;
            
            // Calculate and store flow changes
            dQ_start[i] = newStartFlow - pipe->startFlow;
            dQ_end[i] = newEndFlow - pipe->endFlow;
            dQ[i] = 0.5 * (dQ_start[i] + dQ_end[i]);
            
        }
        // QUASI-STEADY FLOW MODEL
        else {
            handleQuasiSteadyFlowChange(link, H_A, H_B, i);
        }
    }
}

// Helper function for quasi-steady flow calculations
void CGGASolver::handleQuasiSteadyFlowChange(Link* link, double H_A, double H_B, int linkIndex) {
    int n1 = link->fromNode->index;
    int n2 = link->toNode->index;

    
    // Handle pressure regulating valves
    if (link->hGrad == 0.0) {
        if (link->isPRV()) {
            dQ[linkIndex] = -xQ[n2] - link->flow;
        }
        else if (link->isPSV()) {
            dQ[linkIndex] = xQ[n1] - link->flow;
        }
        return;
    }

    // Regular flow change based on head difference
    double dh = (link->fromNode->head + dH[n1]) - (link->toNode->head + dH[n2]);
    double dq = -(dh - link->hLoss) / (link->hGrad); // */

    /*double dh = (link->fromNode->head + dH[n1]) - (link->toNode->head + dH[n2]);
    double dhpast = (link->fromNode->pastHead) - (link->toNode->pastHead);
    double pastTerms = ((1 - kappa) / kappa) * (link->pastHloss - dhpast);
    //double flows = (link->inertialTerm / (kappa * tstep)) * (link->flow - link->pastFlow);
    double flows = (link->inertialTerm / (kappa * tstep)) * (link->pastFlow) + link->hGrad * link->flow;
    double dq = -(dh - link->hLoss + flows - pastTerms) / (link->hGrad + (link->inertialTerm / (kappa * tstep))) + link->flow; // */

    // Special case for constant HP pumps
    if (link->isHpPump() && link->status == Link::LINK_OPEN && dq > link->flow) {
        dq = link->flow / 2.0;
    }
        
    dQ[linkIndex] = -dq;

    // Save flow changes (same for start and end in quasi-steady flow)
    dQ_start[linkIndex] = -dq;
    dQ_end[linkIndex] = -dq; // */
}

//  Compute the coefficient matrix of the linearized set of equations for heads.

void CGGASolver::setUnsteadyMatrixCoeffs(double currentTime)
{
    // Ensure xQ is properly sized
    if (xQ.size() < nodeCount) {
        xQ.resize(nodeCount, 0.0);
    } else {
        // Zero out the vector using vector operations instead of memset
        std::fill(xQ.begin(), xQ.begin() + nodeCount, 0.0);
    }
    
    matrixSolver->reset();
    setUnsteadyLinkCoeffs(currentTime);
    setUnsteadyNodeCoeffs(currentTime);
    setUnsteadyValveCoeffs(currentTime);
}

//  Compute matrix coefficients for links in unsteady flow conditions.
void CGGASolver::setUnsteadyLinkCoeffs(double currentTime)
{
    // Clear node flow balances
    for (int i = 0; i < nodeCount; i++) {
        if (i < xQ.size()) {
            xQ[i] = 0.0;
        }
    }

    // Create and initialize the vector tracking nodes connected to incompressible links
    std::vector<bool> nodeHasIncompressibleLink(nodeCount, false);

    std::vector<int> nodeDomain(nodeCount, -1);
    
    // Process all links for flow contributions
    for (int j = 0; j < linkCount; j++)
    {
        Link* link = network->links[j];
        if (!link) continue;

        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        if (!node1 || !node2) continue;
        
        int n1 = node1->index;
        int n2 = node2->index;

        // Skip links with zero head gradient
        if (link->hGrad == 0.0) {
            continue;
        }

        // Update node flow balances based on flow model and node connectivity
        if (link->type() == Link::PIPE) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                // Check if either node is an "intermediate node" connected to incompressible links
                bool node1IsIntermediate = nodeHasIncompressibleLink[n1];
                bool node2IsIntermediate = nodeHasIncompressibleLink[n2];

                double avgFlow = (pipe->startFlow + pipe->endFlow) / 2.0;
            }
        } else {
            // Incompressible flow links always use single flow value
            xQ[n1] -= link->flow;
            xQ[n2] += link->flow;
        }

        // Handle hydraulic coefficients based on flow model
        if (link->type() == Link::VALVE || link->type() == Link::PUMP) {
            handleQuasiSteadyFlowCoeffs(link, node1, node2, n1, n2);
        }
        else if (link->type() == Link::PIPE) {
            handleCompressibleFlowCoeffs(link, node1, node2, currentTime, n1, n2);
        }
    }
}

//  Compute matrix coefficients for dynamic tanks and external node outflows.
void CGGASolver::setUnsteadyNodeCoeffs(double currentTime)
{
    for (int i = 0; i < nodeCount; i++)
    {
        // ... if node's head not fixed

        Node* node = network->node(i);
        if ( !node->fixedGrade )
        {
            // ... for dynamic tanks, add area terms to row i
            //     of the head solution matrix & r.h.s. vector

            if ( node->type() == Node::TANK && theta != 0.0 )
            {
                Tank* tank = static_cast<Tank*>(node);

                /*double a = tank->area / (theta * tstep);
				matrixSolver->addToDiag(i, a);

				a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
				matrixSolver->addToRhs(i, a); // */

				if (tank->head == tank->pastHead)
				{

					double a = tank->area / (theta * tstep);
					matrixSolver->addToDiag(i, a);

					a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
					matrixSolver->addToRhs(i, a); // 
				}
				
				else
				{
                    double a = (tank->area + ((tank->area - tank->pastArea) / (tank->head - tank->pastHead)) * (tank->head)) / (theta * tstep);
                    matrixSolver->addToDiag(i, a);

                    double b = (tank->pastArea) * tank->pastHead / (theta * tstep) + (1.0 - theta) * tank->pastOutflow / theta;
                    double c = b + ((tank->area - tank->pastArea) / (tank->head - tank->pastHead)) * (tank->head) / (theta * tstep);
                    matrixSolver->addToRhs(i, c); //

                    //xQ[i] -= tank->outflow;
				} // */
            }

            // ... for junctions, add effect of external outflows

            else if ( node->type() == Node::JUNCTION )
            {
                // ... update junction's net inflow
                xQ[i] -= node->outflow;
                double pastOutflow = -(1.0 - theta) * node->pastOutflow / theta;
                matrixSolver->addToDiag(i, node->qGrad);
                matrixSolver->addToRhs(i, node->qGrad * node->head);
                matrixSolver->addToRhs(i, pastOutflow);
            }

            // ... add node's net inflow to r.h.s. row
            matrixSolver->addToRhs(i, (double)xQ[i]);
        }

        // ... if node has fixed head, force solution to produce it

        else
        {
            matrixSolver->setDiag(i, 1.0);
            matrixSolver->setRhs(i, node->head);
        }
    }
}

//  Compute matrix coefficients for pressure regulating valves.
void  CGGASolver::setUnsteadyValveCoeffs(double currentTime)
{
    for (Link* link : network->links)
    {
        // ... skip links that are not active pressure regulating valves

        if (link->hGrad > 0.0) continue;

        // ... determine end node indexes of link

        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... add net inflow of downstream node of a PRV to the
        //     r.h.s. row of its upstream node

        if (link->isPRV())
        {
            matrixSolver->addToRhs(n1, (double)xQ[n2]);
        }

        // ... add net inflow of upstream node of a PSV to the
        //     r.h.s. row of its downstream node

        if (link->isPSV())
        {
            matrixSolver->addToRhs(n2, (double)xQ[n1]);
        }
    }
}

//  Find changes in nodal heads by solving a linearized system of equations.
int CGGASolver::findSteadyHeadChanges()
{
    // ... setup the coeff. matrix of the GGA linearized system

    setSteadyMatrixCoeffs();

    // ... temporarily use the head change array dH[] to store new heads

    double* h = &dH[0];

    // ... solve the linearized GGA system for new nodal heads
    //     (matrixSolver returns a negative integer if it runs successfully;
    //      otherwise it returns the index of the row that caused it to fail.)

    int errorCode = matrixSolver->solve(nodeCount, h);
    if (errorCode >= 0) return errorCode;

    // ... save new heads as head changes

    for (int i = 0; i < nodeCount; i++)
    {
        dH[i] = h[i] - network->nodes[i]->head;
    }

    // ... return a negative number indicating that
    //     the matrix solver ran successfully

    return -1;
}

//  Find the changes in link flows resulting from a set of nodal head changes.
void CGGASolver::findSteadyFlowChanges()
{
    for (int i = 0; i < linkCount; i++)
    {
        // ... get link object and its end node indexes

        dQ[i] = 0.0;
        Link* link = network->link(i);
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... flow change for pressure regulating valves

        if (link->hGrad == 0.0)
        {
            if (link->isPRV()) dQ[i] = -xQ[n2] - link->flow;
            if (link->isPSV()) dQ[i] = xQ[n1] - link->flow;
            continue;
        }

        // ... apply GGA flow change formula:

        double dh = (link->fromNode->head + dH[n1]) -
            (link->toNode->head + dH[n2]);
        double dq = (link->hLoss - dh) / link->hGrad;

        // ... special case to prevent negative flow in constant HP pumps

        if (link->isHpPump() &&
            link->status == Link::LINK_OPEN &&
            dq > link->flow) dq = link->flow / 2.0;

        // ... save flow change

        dQ[i] = -dq;
    }
}

//  Compute the coeffciient matrix of the linearized set of equations for heads.

void CGGASolver::setSteadyMatrixCoeffs()
{
    memset(&xQ[0], 0, nodeCount * sizeof(double));
    matrixSolver->reset();
    setSteadyLinkCoeffs();
    setSteadyNodeCoeffs();
    setSteadyValveCoeffs();
}

//  Compute matrix coefficients for link head loss gradients.

void CGGASolver::setSteadyLinkCoeffs()
{
    // CRITICAL: Double-check that our linkCount matches reality
    if (linkCount != network->links.size()) {
        network->msgLog << "\nWARNING: linkCount mismatch in setSteadyLinkCoeffs() - " 
                       << linkCount << " vs " << network->links.size();
        // Update to correct value
        linkCount = network->links.size();
    }
    
    // Clear node flow balances
    for (int i = 0; i < nodeCount; i++) {
        if (i < xQ.size()) {
            xQ[i] = 0.0;
        }
    }

    for (int j = 0; j < linkCount; j++)
    {
        Link* link = network->links[j];

        // Skip links with zero head gradient.
        if (link->hGrad == 0.0) continue;

        // Identify end nodes and their indices.
        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        int n1 = node1->index;
        int n2 = node2->index;

        // Update node flow balances.
        xQ[n1] -= link->flow;
        xQ[n2] += link->flow;

        // Compute coefficients: a is the link contribution factor, and b is head loss.
        double a = 1.0 / link->hGrad;
        double b = a * link->hLoss;

        //double b = a * (link->hLoss - link->hGrad * link->flow);

        // Update off-diagonal coefficient if both nodes are not fixed.
        if (!node1->fixedGrade && !node2->fixedGrade) {
            matrixSolver->addToOffDiag(j, -a);
        }

        // For the start node:
        if (node1->fixedGrade) {
            matrixSolver->addToRhs(n2, a * node1->head);
        }
        else {
            matrixSolver->addToDiag(n1, a);
            matrixSolver->addToRhs(n1, b);
        }

        // For the end node:
        if (node2->fixedGrade) {
            matrixSolver->addToRhs(n1, a * node2->head);
        }
        else {
            matrixSolver->addToDiag(n2, a);
            matrixSolver->addToRhs(n2, -b);
        }
    }
}

//  Compute matrix coefficients for dynamic tanks and external node outflows.
void  CGGASolver::setSteadyNodeCoeffs()
{
    for (int i = 0; i < nodeCount; i++)
    {
        // ... if node's head not fixed

        Node* node = network->nodes[i];
        if (!node->fixedGrade)
        {
            // ... for dynamic tanks, add area terms to row i
            //     of the head solution matrix & r.h.s. vector

            if (node->type() == Node::TANK && theta != 0.0)
            {
                Tank* tank = static_cast<Tank*>(node);
                double a = tank->area / (theta * tstep);
                matrixSolver->addToDiag(i, a);

                a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
                matrixSolver->addToRhs(i, a);
            }

            // ... for junctions, add effect of external outflows

            else if (node->type() == Node::JUNCTION)
            {
                // ... update junction's net inflow
                xQ[i] -= node->outflow;
                matrixSolver->addToDiag(i, node->qGrad);
                matrixSolver->addToRhs(i, node->qGrad * node->head);
            }

            // ... add node's net inflow to r.h.s. row
            matrixSolver->addToRhs(i, (double)xQ[i]);
        }

        // ... if node has fixed head, force solution to produce it

        else
        {
            matrixSolver->setDiag(i, 1.0);
            matrixSolver->setRhs(i, node->head);
        }
    }
}

//  Compute matrix coefficients for pressure regulating valves.
void  CGGASolver::setSteadyValveCoeffs()
{
    for (Link* link : network->links)
    {
        // ... skip links that are not active pressure regulating valves

        if (link->hGrad > 0.0) continue;

        // ... determine end node indexes of link

        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... add net inflow of downstream node of a PRV to the
        //     r.h.s. row of its upstream node

        if (link->isPRV())
        {
            matrixSolver->addToRhs(n1, (double)xQ[n2]);
        }

        // ... add net inflow of upstream node of a PSV to the
        //     r.h.s. row of its downstream node

        if (link->isPSV())
        {
            matrixSolver->addToRhs(n2, (double)xQ[n1]);
        }
    }
}

// Check convergence of the solution
bool CGGASolver::hasConverged() {
    return
        (hydBalance.maxHeadErr < headErrLimit) &&
        (hydBalance.maxFlowErr < flowErrLimit) &&
        (hydBalance.maxFlowChange < flowChangeLimit) &&
        (hydBalance.totalFlowChange < flowRatioLimit);
}

void CGGASolver::handleCompressibleFlowCoeffs(Link* link, Node* node1, Node* node2, double currentTime, int n1, int n2) {
    // Skip if link isn't a pipe
    Pipe* pipe = dynamic_cast<Pipe*>(link);
    if (!pipe) return;        
    
    // === CALCULATE CHARACTERISTIC PARAMETERS ===
    
    // Physical constants
    double c = pipe->getWaveSpeed();
    double A = pipe->getArea();
    double L = pipe->length;
    double g = GRAVITY;
    double Bj = c / (g * A);  // Wave impedance term (B')
    double D = pipe->diameter;
    double epsilon = 0.5;  // Friction integration parameter
    
    // CORRECTED: Calculate proper reach-back time based on discretization
    double reachBackTime;
    
   // Calculate wave travel time and reach-back time
    double physicalWaveTime = L / c;
            
    // Determine reach-back index based on number of reaches
    int reachBackSteps = pipe->numReaches;
    double waveTravel = reachBackSteps * tstep;  // Wave travel distance

    // Retrieve reach-back flows and heads from history
    double Qa = pipe->pastStartFlow;  // Default fallback
    double Qb = pipe->pastEndFlow;    // Default fallback
    double Ha = pipe->fromNode->pastHead;  // Default fallback
    double Hb = pipe->toNode->pastHead;    // Default fallback
            
    // Use the FlowHistoryManager to get reach-back values
    FlowHistoryResult historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, network);
    //FlowHistoryResult historyResult = getReachBackValues(pipe, currentTime, physicalWaveTime);

    // If history was found, use the values from the result
    if (historyResult.found) {
        Qa = historyResult.startFlow;
        Qb = historyResult.endFlow;
        Ha = historyResult.startHead;
        Hb = historyResult.endHead;
    } 

    // Convert gradients to resistance coefficients for characteristic equations
    double KA, KB;

    if (network->option(Options::HEADLOSS_MODEL) == "D-W")
    {
        double frictionFactor = pipe->getFrictionFactor(pipe->getRe(pipe->flow, network->option(Options::KIN_VISCOSITY)));
        double K = frictionFactor * pipe->length / (2 * g * D * A * A);

        double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, network->option(Options::KIN_VISCOSITY)));
        double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, network->option(Options::KIN_VISCOSITY)));
                // Split K for each end as needed (KA for upstream, KB for downstream)
        KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
        KB = frictionFactorB * pipe->length / (2 * g * D * A * A);

    }
    else
    {
        KA = pipe->resistance;
        KB = pipe->resistance;
    }

    // === ENSURE NUMERICAL STABILITY ===
    const double MIN_RESISTANCE = 1e-10;
    KA = std::max(KA, MIN_RESISTANCE);
    KB = std::max(KB, MIN_RESISTANCE);
    
    // === UNSTEADY FRICTION TERMS ===
    double kUFa = pipe->computeShearDecayCoeff(Qa);
    double kUFb = pipe->computeShearDecayCoeff(Qb);
    double Ua = 0.0;
    double Ub = 0.0;
    
    if (kUFa > 0.0 && kUFb > 0.0) {
        auto sign = [](double val) -> double { 
            return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0); 
        };
        
        // CORRECTED: Unsteady friction terms
        Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
        Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
    }
    
    // CORRECTED: Characteristic equation coefficients
    double BA, BB, RA, RB;

    if (network->option(Options::HEADLOSS_MODEL) == "D-W") {
        // For boundary A: using C+ from point b
        BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
        RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;

        // For boundary B: using C- from point a
        BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
        RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
    }
    else if (network->option(Options::HEADLOSS_MODEL) == "H-W") {
        // For boundary A: using C+ from point b
        BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
        RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
        
        // For boundary B: using C- from point a
        BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
        RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;

    }

    // === APPLY WAVE ATTENUATION SCHEME IF ENABLED ===
    if (pipe->impingingWaveCount > 0) {
        //updateWaveAttenuationValues(pipe, BA, BB, RA, RB);
    }
    
    // Store for reference and debugging
    pipe->charB_plus = BA;  // Coefficient for boundary A (using C+)
    pipe->charB_minus = BB; // Coefficient for boundary B (using C-)
    pipe->charC_plus = RA;  // Reach-back term for boundary A
    pipe->charC_minus = RB; // Reach-back term for boundary B

    // === APPLY BOUNDARY CONDITIONS TO MATRIX ===

    if (!node1->fixedGrade) 
    {
        double diagContribA = 1.0 / BA;
        
        matrixSolver->addToRhs(n1, -RA / BA);
        matrixSolver->addToDiag(n1, 1.0 / BA);
    }
    else {
        double rhsValue = node1->head / BA;
        //matrixSolver->addToRhs(n2, rhsValue);
    }

    if (!node2->fixedGrade) 
    {
        double diagContribB = 1.0 / BB;
        
        matrixSolver->addToRhs(n2, RB / BB);
        matrixSolver->addToDiag(n2, 1.0 / BB);
    }
    else {
        double rhsValue = node2->head / BB;
        //matrixSolver->addToRhs(n1, rhsValue);
    }

    if (RB == 0.00 || RA == 0.00) {
        network->msgLog << "\nWARNING: Zero coefficient detected in CGGASolver::handleCompressibleFlowCoeffs for pipe " 
                       << pipe->name << " at time " << currentTime;
    }
}

void CGGASolver::configurePipeForTimestep(Pipe* pipe, double timestep, Network* network) {
    // This function configures a pipe's discretization to work correctly
    // with the given timestep while satisfying the Courant condition
    
    double c = pipe->getWaveSpeed();
    double L = pipe->length;
    
    // Calculate minimum number of reaches needed for stability
    // From Courant condition: Δt ≤ Δx/c, where Δx = L/N
    // Therefore: N ≥ L/(c·Δt)
    
    double minReachesFloat = L / (c * timestep);
    int minReaches = static_cast<int>(std::ceil(minReachesFloat));
    
    // Ensure at least 1 reach
    if (minReaches < 1) minReaches = 1;
    
    // Limit maximum reaches for computational efficiency
    const int MAX_REACHES = 100;
    
    // Set the number of reaches
    pipe->numReaches = minReaches;
    
    // Calculate and verify the configuration
    double dx = L / minReaches;
    double actualCourant = c * timestep / dx;
}

// Handle quasi-steady flow
void CGGASolver::handleQuasiSteadyFlowCoeffs(Link* link, Node* node1, Node* node2, int n1, int n2)
{
    double a = 1.0 / link->hGrad;
    double b = a * link->hLoss;

    // Calculate coefficient a = 1/hGrad

    if (link->hLoss == 0.0) {
        b = MIN_GRADIENT; // Avoid division by zero if hLoss is zero
    }
    // Process upstream node (node1)
    if (!node1->fixedGrade) {
        // Add to diagonal
        matrixSolver->addToDiag(n1, a);
        
        // Add to RHS
        matrixSolver->addToRhs(n1, b);
    } else {
        matrixSolver->addToRhs(n2, a * node1->head);
        //matrixSolver->addToRhs(n2, b); // Only for trial
    }
    
    // Process downstream node (node2)
    if (!node2->fixedGrade) {
        // Add to diagonal
        matrixSolver->addToDiag(n2, a);
        
        // Add to RHS
        matrixSolver->addToRhs(n2, -b);
    } else {
        matrixSolver->addToRhs(n1, a * node2->head);
        //matrixSolver->addToRhs(n1, -b); // Only for trial
    }
    
    // Add off-diagonal terms if both nodes are not fixed
    if (!node1->fixedGrade && !node2->fixedGrade) {
        matrixSolver->addToOffDiag(link->index, -a);
    }

    if (a == 0.0 || b == 0.0) {
        network->msgLog << "\nWARNING: Zero coefficient detected in CGGASolver::handleQuasiSteadyFlowCoeffs for link " 
                       << link->name;
    }
}

void CGGASolver::updateSolution(double lamda, bool useShadowMode) {
    // First, ensure vectors are properly sized
    if (currentHeads.size() != nodeCount) {
        currentHeads.resize(nodeCount, 0.0);
    }

    if (currentFlows.size() != linkCount) {
        currentFlows.resize(linkCount, 0.0);
    }
    
    // Initialize shadow vectors if needed
    if (useShadowMode) {
        if (shadowHeads.size() != nodeCount) {
            shadowHeads.resize(nodeCount);
            // Copy current heads to shadow heads if first time
            for (int i = 0; i < nodeCount && i < network->nodes.size(); i++) {
                if (network->nodes[i]) {
                    shadowHeads[i] = network->nodes[i]->head;
                }
            }
        }
        
        if (shadowFlows.size() != linkCount) {
            shadowFlows.resize(linkCount);
            // Copy current flows to shadow flows if first time
            for (int i = 0; i < linkCount && i < network->links.size(); i++) {
                if (network->links[i]) {
                    shadowFlows[i] = network->links[i]->flow;
                }
            }
        }
    }

    // Update nodes
    for (int i = 0; i < nodeCount && i < dH.size(); i++) {
        if (network->nodes[i]) {
            if (useShadowMode) {
                // Update shadow values only
                shadowHeads[i] += lamda * dH[i];
                currentHeads[i] = shadowHeads[i]; // For flow indicators
            } else {
                // Normal update
                network->nodes[i]->head += lamda * dH[i];
                currentHeads[i] = network->nodes[i]->head;
            }
        }
    }

    // Update links
    for (int i = 0; i < linkCount && i < dQ.size(); i++) {
        if (network->links[i]) {
            Link* link = network->links[i];
            
            if (useShadowMode) {
                // Update shadow flows only
                shadowFlows[i] += lamda * dQ[i];
                currentFlows[i] = shadowFlows[i];
            } else {
                // Normal updates based on flow model
                if (link->type() == Link::PIPE) {
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    if (pipe) {
                        pipe->startFlow += lamda * dQ[i];
                        pipe->endFlow += lamda * dQ[i];
                        pipe->flow += lamda * dQ[i];
                        currentFlows[i] = 0.5 * (pipe->startFlow + pipe->endFlow);
                    } else {
                        link->flow += lamda * dQ[i];
                        currentFlows[i] = link->flow;
                    }
                } else {
                    link->flow += lamda * dQ[i];
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    if (pipe) {
                        pipe->startFlow = link->flow;
                        pipe->endFlow = link->flow;
                    }
                    currentFlows[i] = link->flow;
                }
            }
        }
    }
}

void CGGASolver::updateWHSolution(double lamda) {
    // Track maximum solution changes for monitoring convergence
    double maxHeadChange = 0.0;
    double maxFlowChange = 0.0;
    
    // Apply progressive relaxation in early time steps
    double relaxationFactor = 1.0;
    const int initialStepsCount = 10;
    if (tstep < initialStepsCount) {
        // Start with heavier relaxation, gradually increasing to full strength
        relaxationFactor = 0.3 + 0.7 * (tstep / (double)initialStepsCount);
        network->msgLog << "\nApplying relaxation factor: " << relaxationFactor;
        lamda *= 1.0; // relaxationFactor;
    }
    
    // Resize solution tracking vectors if needed
    if (currentHeads.size() != nodeCount) currentHeads.resize(nodeCount);
    if (currentFlows.size() != linkCount) currentFlows.resize(linkCount);

    // ==== STEP 1: UPDATE NODE HEADS ====
    for (int i = 0; i < nodeCount && i < dH.size(); i++) {
        Node* node = network->nodes[i];
        if (!node) continue;  // Skip null nodes
        
        // Skip fixed-grade nodes (reservoirs) - their heads don't change
        if (node->fixedGrade) continue;

        // Store previous head for tracking changes
        double prevHead = node->head;
        
        // Apply change with lambda damping factor
        node->head += lamda * dH[i];
        
        // Track maximum head change
        double headChange = std::abs(node->head - prevHead);
        if (headChange > maxHeadChange) {
            maxHeadChange = headChange;
        }
        
        // Store updated value
        currentHeads[i] = node->head;
    }

    // ==== STEP 2: UPDATE LINK FLOWS ====
    // First update compressible links (water hammer pipes)
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;  // Skip null links
        
        // ==== GAWH HANDLING FOR PIPES ====
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe) {
            // Store previous values for stability control and change tracking
            double prevStartFlow = pipe->startFlow;
            double prevEndFlow = pipe->endFlow;
            double prevFlow = pipe->flow;
            
            // Normal pipes - apply flow changes with lambda damping
            if (i < dQ_start.size()) pipe->startFlow += lamda * dQ_start[i];
            if (i < dQ_end.size()) pipe->endFlow += lamda * dQ_end[i];

            // Set reference flow (average for reporting)
            pipe->flow = 0.5 * (pipe->startFlow + pipe->endFlow);
            currentFlows[i] = pipe->flow;
            
            // Track maximum flow change
            double flowChange = std::max(
                std::abs(pipe->startFlow - prevStartFlow),
                std::abs(pipe->endFlow - prevEndFlow)
            );
            if (flowChange > maxFlowChange) {
                maxFlowChange = flowChange;
            }
        }
        // Handle all other links after water hammer pipes
    }
    
    // Now update incompressible flow links
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;
        
        // Skip water hammer pipes (already updated)
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe) continue;
        
        // Process all other links (incompressible)
        // Store previous value for stability control
        double prevFlow = link->flow;
            
        // Apply changes with lambda factor
        link->flow += lamda * dQ[i];
           
        // For pipes with incompressible flow, set start/end flows equal
        if (pipe) {
            pipe->startFlow = link->flow;
            pipe->endFlow = link->flow;
            pipe->flow = link->flow ;  // Use average flow for reporting
        }
        
        currentFlows[i] = link->flow;
        
        // Track maximum flow change
        double flowChange = std::abs(link->flow - prevFlow);
        if (flowChange > maxFlowChange) {
            maxFlowChange = flowChange;
        }
    }
    
    // Log solution update metrics
    network->msgLog << "\nSolution update metrics:"
                   << "\n  Max head change: " << maxHeadChange
                   << "\n  Max flow change: " << maxFlowChange;
}

// Dynamically check for link status changes
bool CGGASolver::linksChangedStatus() {
    bool statusChanged = false;

    for (int i = 0; i < linkCount; i++) {
        Link* link = network->link(i);
        if (!link) continue;

        double q = link->flow;                  // Current flow
        double h1 = link->fromNode->head;       // Upstream head
        double h2 = link->toNode->head;         // Downstream head

        // Update link status based on flow and head conditions
        link->updateStatus(q, h1, h2);

        // ... check for flow into full or out of empty tanks

        if ( link->status > Link::LINK_CLOSED )
        {
            if ( link->fromNode->isClosed(q) || link->toNode->isClosed(-q) )
            {
                link->status = Link::TEMP_CLOSED;
                link->flow = ZERO_FLOW;
                
            }
        }

        // Check if the link status changed
        if (link->status != link->previousStatus) {
            statusChanged = true;
            link->previousStatus = link->status;
            if (link->status == Link::TEMP_CLOSED || link->status == Link::LINK_CLOSED) {
                link->flow = ZERO_FLOW;
                dQ[i] = 0.0;
            }
        }
    }

    return statusChanged;
}

// Report the results of each trial
void CGGASolver::reportTrial(int trials, double stepSize) {
	std::cout << "Trial " << trials << ": Step size = " << stepSize << std::endl;
	std::cout << "Total error norm = " << errorNorm << std::endl;
}


// Find the error norm for a given step size
double CGGASolver::findErrorNorm(double lamda, double currentTime, double tstep, double currentMode) {
	
    // Safety check for array sizes
    if (xQ.size() < nodeCount || dH.size() < nodeCount || dQ.size() < linkCount) {
        std::cerr << "Error: Vector sizes too small in findErrorNorm" << std::endl;
        return std::numeric_limits<double>::max();
    }

    hLossEvalCount++;
    return hydBalance.evaluate(lamda, (double*)&dH[0], (double*)&dQ[0],
        (double*)&xQ[0], network, currentTime, tstep, currentMode);
}

// Find the optimal step size for the next trial
double CGGASolver::findStepSize(int trials, double currentTime) {
	double lamda = 1.0;
	errorNorm = findErrorNorm(lamda, currentTime, tstep, currentMode);
	return errorNorm;
}

double CGGASolver::findWHStepSize(int trials, double currentTime, double dl, double minErrorNorm) {
    // Initialize variables for your experimental approach
    double bestLambda = 1.0;
    double testError = 1000000.0;  // Start with a very large error
    minErrorNorm = 0.0;
    
    // Calculate how many lambda values to test based on dl
    lambdaNumber = static_cast<int>(1.0 / dl);
    Lambda.resize(lambdaNumber, 0.0);
    
    // Clear the lambda array
    memset(&Lambda[0], 0, lambdaNumber * sizeof(double));
    
    network->msgLog << "\nTesting experimental step sizing approach with " << lambdaNumber << " lambda values:";
    
    // Your experimental approach: recompute search direction for each lambda
    for (int i = 0; i < lambdaNumber; i++) {
        // Calculate this test lambda value
        Lambda[i] = (i + 1) * dl;
        
        network->msgLog << "\n  Testing lambda = " << Lambda[i];
        
        // KEY DIFFERENCE: Recompute search direction for this specific lambda
        // This is what you want to test - does this improve water hammer convergence?
        int errorCode = findUnsteadyHeadChanges(currentTime);
        if (errorCode >= 0) {
            Node* node = network->node(errorCode);
            network->msgLog << "\n  Matrix singular at lambda = " << Lambda[i] << " for node " << node->name;
            continue; // Skip this lambda and try the next one
        }
        
        // Recompute flow changes for this lambda
        findUnsteadyFlowChanges(currentMode, currentTime);
        
        // Calculate error norm with this lambda and these recomputed changes
        errorNorm = findUnsteadyErrorNorm(Lambda[i], currentTime, tstep);
        
        network->msgLog << ", Error = " << errorNorm;
        
        // Check if this is the best lambda so far
        if (errorNorm < testError) {
            testError = errorNorm;
            bestLambda = Lambda[i];
            
            // Apply this solution since it's the best so far
            updateWHSolution(Lambda[i]);
            
            network->msgLog << " <- New best";
        }
        
        // Early termination if we find a very good solution
        if (errorNorm < 1e-10) {
            network->msgLog << "\n  Excellent convergence found - stopping search";
            break;
        }
    }
    
    // Store the best error found for the main solver loop
    minErrorNorm = testError;
    
    network->msgLog << "\nExperimental method selected lambda = " << bestLambda 
                   << " with final error = " << minErrorNorm;
    
    return minErrorNorm;
}

void CGGASolver::setFixedGradeNodes()
{
    // Loop through each active link to check for PRV/PSV valves.
    for (Link* link : network->links)
    {
        Node* node = nullptr;  // Declare a local pointer for this iteration.

        // For a PRV, the control node is the toNode;
        // for a PSV, the control node is the fromNode.
        if (link->isPRV())
            node = link->toNode;
        else if (link->isPSV())
            node = link->fromNode;
        else
            continue;  // Skip other types.

        // Set the fixed grade status based on link status.
        if (link->status == Link::VALVE_ACTIVE)
        {
            node->fixedGrade = true;
            node->head = link->setting + node->elev;
        }
        else {
            node->fixedGrade = false;
        }
    }

    // After time 0, if theta and tstep are nonzero, ensure that tank nodes are not fixed.
    if (theta > 0.0 && tstep > 0.0)
    {
        for (Node* tankNode : network->nodes)
        {
            if (tankNode->type() == Node::TANK)
                tankNode->fixedGrade = false;
        }
    }
}

// Sign function helper
int CGGASolver::sign(double val) {
    return (val > 0.0) ? 1 : ((val < 0.0) ? -1 : 0);
}    