/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION */


// dualflowhistory.h
#ifndef DUALFLOWHISTORY_H
#define DUALFLOWHISTORY_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Core/network.h"
#include <set>
#include <tuple>
#include <deque>

// Forward declarations
class Pipe;
class Pump;
class Link;
class Node;
class Network;
class Junction;

/**
 * @struct DualFlowHistory
 * @brief Stores historical flow and head data for water hammer analysis
 *
 * This structure maintains time-series data for pipe flows and heads at both
 * ends of a pipe, allowing for reach-back operations needed in algebraic
 * water hammer (AWH) calculations.
 */


 /**
  * @struct FlowHistoryResult
  * @brief Simple structure to return reach-back flow and head values
  *
  * This structure encapsulates the results of a reach-back operation,
  * providing a clean interface for retrieving historical values.
  */

struct PipeHistoryEntry {
    double time;           // Simulation time for this entry
    double startFlow;      // Flow at pipe start (upstream)
    double endFlow;        // Flow at pipe end (downstream)
    double startHead;      // Head at pipe start node
    double endHead;        // Head at pipe end node
    double timestep;       // Timestep used when this was recorded
    int    solverMode;     // 0=QS, 1=RWC, 2=WH

    PipeHistoryEntry()
        : time(0), startFlow(0), endFlow(0), startHead(0), endHead(0),
        timestep(0), solverMode(0) {
    }

    PipeHistoryEntry(double t, double qStart, double qEnd, double hStart, double hEnd,
        double dt, int mode)
        : time(t), startFlow(qStart), endFlow(qEnd), startHead(hStart), endHead(hEnd),
        timestep(dt), solverMode(mode) {
    }
};

struct FlowHistoryResult {
    bool   found;       // Whether valid history was found
    double startFlow;   // Interpolated flow at pipe start
    double endFlow;     // Interpolated flow at pipe end
    double startHead;   // Interpolated head at pipe start
    double endHead;     // Interpolated head at pipe end
    double time;        // The actual time of the retrieved/interpolated data
    bool   wasInterpolated;  // True if values were interpolated

    FlowHistoryResult()
        : found(false), startFlow(0), endFlow(0), startHead(0), endHead(0),
        time(0), wasInterpolated(false) {
    }
};

struct WaveAttenuationPoint {
    double flow = 0.0;
    double head = 0.0;
    double position = 0.0;
	double amplitude = 0.0;
}; // */

class DualFlowHistory {
private:
    // Fixed-size storage vectors (pre-allocated)
    std::vector<double> timestamps;
    std::vector<std::vector<double>> startFlowHistory;
    std::vector<std::vector<double>> endFlowHistory;
    std::vector<std::vector<double>> startHeadHistory;
    std::vector<std::vector<double>> endHeadHistory;
    std::vector<std::vector<double>> junctionHeadHistory;
    std::vector<std::vector<WaveAttenuationPoint>> waveAttenuationHistory;

    // Circular buffer management
    size_t capacity;      // Maximum number of entries
    size_t currentSize;   // Current number of valid entries

public:
    /**
     * Constructor with specified capacity
     * @param maxCapacity Maximum number of time points to store (default 1000)
     */
    explicit DualFlowHistory(size_t maxCapacity = 1000) :
        capacity(maxCapacity),
        currentSize(0)
    {

        // Pre-allocate all storage to avoid reallocations
        timestamps.resize(capacity);
        startFlowHistory.resize(capacity);
        endFlowHistory.resize(capacity);
        startHeadHistory.resize(capacity);
        endHeadHistory.resize(capacity);
        junctionHeadHistory.resize(capacity);
        waveAttenuationHistory.resize(capacity);
    }

    // Copy constructor for compatibility
    DualFlowHistory(const DualFlowHistory& other) = default;

    // Assignment operator
    DualFlowHistory& operator=(const DualFlowHistory& other) = default;

    /**
     * Clear all history data and reset indices
     */
    void clear();

    /**
     * Get the current number of valid entries
     */
    size_t size() const { return currentSize; }

    /**
     * Check if the buffer is empty
     */
    bool empty() const { return currentSize == 0; }

    /**
     * Get the buffer capacity
     */
    size_t getCapacity() const { return capacity; }

    // Friend class for any special access needs
    friend class FlowHistoryManager;
};


class FlowHistoryManager {
public:
    static FlowHistoryManager& getInstance() {
        static FlowHistoryManager instance;
        return instance;
    }

    // Delete copy constructor and assignment operator
    FlowHistoryManager(const FlowHistoryManager&) = delete;
    FlowHistoryManager& operator=(const FlowHistoryManager&) = delete;

    // Core functionality
    void initializeHistory(Network* network, double currentTime, double timeStep, double historyStartTime);
    void addHistory(Pipe* pipe, const DualFlowHistory& history);
    FlowHistoryResult getReachBackValues(Pipe* pipe, double currentTime,
        double reachBackTime, Network* nw);

    void FlowHistoryManager::updateHistoryForPipe(
        Pipe* pipe,
        double currentTime,
        double actualTimeStep,
        const std::vector<double>& startFlows,
        const std::vector<double>& endFlows,
        const std::vector<double>& startHeads,
        const std::vector<double>& endHeads,
        const std::vector<double>& junctionHeads,
        const std::vector<WaveAttenuationPoint>& wavePoints);

    // Utility functions
    void clear();
    bool hasHistoryForPipe(Pipe* pipe) const;
    double FlowHistoryManager::linearInterpolate(double v1, double v2, double alpha);

private:
    FlowHistoryManager()
        : currentSolverMode(0),
        lastModeChangeTime(0.0) {}
    std::map<Pipe*, std::unique_ptr<DualFlowHistory>> historyData;

    int currentSolverMode;

    double lastModeChangeTime;

    // Singleton instance
    static FlowHistoryManager* instance;

    double transitionTime;

};
#endif // DUALFLOWHISTORY_H