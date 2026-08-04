//! \file test_flowhistory.cpp
//! \brief Unit tests for the flow history subsystem of EPANET-CGGA.
//!
//! The history subsystem stores committed pipe flows and heads and serves
//! reach-back queries during water hammer steps. Its correctness is a
//! prerequisite for the characteristic equations at every mode transition,
//! so it is verified here in isolation from the solver.
//!
//! No external test framework is required: the file builds as a standalone
//! executable and returns a non-zero exit code if any check fails.

#include "Core/dualflowhistory.h"
#include "Elements/pipe.h"

#include <cmath>
#include <cstdio>
#include <string>

namespace {

int  failures = 0;
int  checks   = 0;

void check(bool condition, const std::string& name)
{
    ++checks;
    if (!condition) { ++failures; std::printf("  FAIL  %s\n", name.c_str()); }
    else            {             std::printf("  ok    %s\n", name.c_str()); }
}

void checkClose(double actual, double expected, double tol, const std::string& name)
{
    check(std::fabs(actual - expected) <= tol, name);
}

// ---------------------------------------------------------------------------
// DualFlowHistory: circular buffer state
// ---------------------------------------------------------------------------
void testBufferConstruction()
{
    std::printf("DualFlowHistory construction\n");

    DualFlowHistory h;
    check(h.empty(),                       "default buffer is empty");
    check(h.size() == 0,                   "default buffer reports zero entries");
    check(h.getCapacity() == 1000,         "default capacity is 1000 entries");

    DualFlowHistory small(16);
    check(small.getCapacity() == 16,       "explicit capacity is honoured");
    check(small.empty(),                   "explicit-capacity buffer starts empty");

    DualFlowHistory copy(small);
    check(copy.getCapacity() == 16,        "copy preserves capacity");
    check(copy.empty(),                    "copy preserves empty state");
}

void testBufferClear()
{
    std::printf("DualFlowHistory clear\n");

    DualFlowHistory h(8);
    h.clear();
    check(h.empty(),                       "clear on an empty buffer leaves it empty");
    check(h.getCapacity() == 8,            "clear does not change capacity");
}

// ---------------------------------------------------------------------------
// Interpolation: the arithmetic the reach-back path depends on
// ---------------------------------------------------------------------------
void testLinearInterpolation()
{
    std::printf("FlowHistoryManager::linearInterpolate\n");

    FlowHistoryManager& m = FlowHistoryManager::getInstance();

    checkClose(m.linearInterpolate(10.0, 20.0, 0.0), 10.0, 1e-12,
               "alpha = 0 returns the first value exactly");
    checkClose(m.linearInterpolate(10.0, 20.0, 1.0), 20.0, 1e-12,
               "alpha = 1 returns the second value exactly");
    checkClose(m.linearInterpolate(10.0, 20.0, 0.5), 15.0, 1e-12,
               "alpha = 0.5 returns the midpoint");
    checkClose(m.linearInterpolate(10.0, 20.0, 0.25), 12.5, 1e-12,
               "interpolation is linear in alpha");

    // Sign handling: reach-back values routinely straddle zero during a transient.
    checkClose(m.linearInterpolate(-5.0, 5.0, 0.5), 0.0, 1e-12,
               "interpolation across a sign change returns zero at the midpoint");
    checkClose(m.linearInterpolate(-8.0, -2.0, 0.5), -5.0, 1e-12,
               "interpolation between two negative values is correct");

    // Degenerate interval: both endpoints equal.
    checkClose(m.linearInterpolate(3.7, 3.7, 0.37), 3.7, 1e-12,
               "equal endpoints return that value for any alpha");
}

// ---------------------------------------------------------------------------
// Manager identity and lifetime
// ---------------------------------------------------------------------------
void testManagerSingleton()
{
    std::printf("FlowHistoryManager singleton\n");

    FlowHistoryManager& a = FlowHistoryManager::getInstance();
    FlowHistoryManager& b = FlowHistoryManager::getInstance();
    check(&a == &b, "getInstance returns the same object on repeated calls");
}

// ---------------------------------------------------------------------------
// Registration and retrieval for a single pipe
// ---------------------------------------------------------------------------
void testRegistrationAndClear()
{
    std::printf("FlowHistoryManager registration\n");

    FlowHistoryManager& m = FlowHistoryManager::getInstance();
    m.clear();

    Pipe pipe("TEST_PIPE_1");
    check(!m.hasHistoryForPipe(&pipe),
          "a pipe with no registered history reports none");

    DualFlowHistory h(32);
    m.addHistory(&pipe, h);
    check(m.hasHistoryForPipe(&pipe),
          "a pipe reports history once it has been registered");

    Pipe other("TEST_PIPE_2");
    check(!m.hasHistoryForPipe(&other),
          "registration is per pipe and does not leak to other pipes");

    m.clear();
    check(!m.hasHistoryForPipe(&pipe),
          "clear removes registered history");
}

// ---------------------------------------------------------------------------
// Reach-back outside the retained window must be reported, not extrapolated
// ---------------------------------------------------------------------------
void testReachBackOutsideWindow()
{
    std::printf("FlowHistoryManager reach-back guards\n");

    FlowHistoryManager& m = FlowHistoryManager::getInstance();
    m.clear();

    Pipe pipe("TEST_PIPE_3");
    FlowHistoryResult r = m.getReachBackValues(&pipe, 100.0, 99.0, nullptr);
    check(!r.found,
          "reach-back on an unregistered pipe reports found == false");

    DualFlowHistory h(8);
    m.addHistory(&pipe, h);
    FlowHistoryResult r2 = m.getReachBackValues(&pipe, 100.0, 99.0, nullptr);
    check(!r2.found,
          "reach-back into an empty history reports found == false rather than "
          "returning silently extrapolated values");

    m.clear();
}

} // namespace

int main()
{
    std::printf("EPANET-CGGA flow history unit tests\n");
    std::printf("-----------------------------------\n");

    testBufferConstruction();
    testBufferClear();
    testLinearInterpolation();
    testManagerSingleton();
    testRegistrationAndClear();
    testReachBackOutsideWindow();

    std::printf("-----------------------------------\n");
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
