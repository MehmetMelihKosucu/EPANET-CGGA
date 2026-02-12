/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "pattern.h"
#include "Utilities/mempool.h"

#include <limits>
using namespace std;

//-----------------------------------------------------------------------------

// Pattern Factory

Pattern* Pattern::factory(int type_, string name_, MemPool* memPool)
{
    switch ( type_ )
    {
    case FIXED_PATTERN:
        return new(memPool->alloc(sizeof(FixedPattern))) FixedPattern(name_);
        break;
    case VARIABLE_PATTERN:
        return new(memPool->alloc(sizeof(VariablePattern))) VariablePattern(name_);
        break;
    default:
        return nullptr;
    }
}

//-----------------------------------------------------------------------------

// Abstract Pattern Constructor

Pattern::Pattern(string name_, int type_) :
    Element(name_),
    type(type_),
    currentIndex(0),
    interval(0.0),
    interpType(LINEAR)  // Default to LINEAR for water hammer
{}

Pattern::~Pattern()
{
    factors.clear();
}

//-----------------------------------------------------------------------------

//  Returns a Pattern's factor value at the current time period.

double Pattern::currentFactor()
{
    if ( factors.size() == 0 ) return 1.0;
    return factors[currentIndex];
}

double Pattern::factorAt(double t)
{
    if (interpType == DISCRETE) {
        return currentFactor();
    }
    else {
        return interpolatedFactor(t);
    }
}

// Fixed Pattern Constructor/Destructor

FixedPattern::FixedPattern(string name_) :
    Pattern(name_, FIXED_PATTERN),
    startTime(0.0)
{}

FixedPattern::~FixedPattern() {}

//-----------------------------------------------------------------------------

//  Initializes the state of a Fixed Pattern.

void FixedPattern::init(double intrvl, double tStart)
{
    startTime = tStart;
    if (interval <= 0.0) interval = intrvl;
    if (factors.size() == 0) factors.push_back(1.0);

    int nPeriods = static_cast<int>(factors.size());
    if (interval > 0.0)
    {
        currentIndex = static_cast<int>(startTime / interval) % nPeriods;
    }
    currentTime = 0.0;
}

//-----------------------------------------------------------------------------

//  Finds the time (sec) until the next change in a FixedPattern.

double FixedPattern::nextTime(double t)
{
    if (interval <= 0.0) return numeric_limits<double>::max();

    int nPeriods = static_cast<int>((startTime + t) / interval);
    return (nPeriods + 1) * interval;
}

//-----------------------------------------------------------------------------

//  Advances a FixedPattern to the period associated with time t (sec).

void FixedPattern::advance(double t)
{
    currentTime = t;  // Store current time for interpolation

    if (interval <= 0.0 || factors.size() == 0) return;

    int nPeriods = static_cast<int>((startTime + t) / interval);
    currentIndex = nPeriods % static_cast<int>(factors.size());
}

double FixedPattern::interpolatedFactor(double t)
{
    if (factors.size() == 0) return 1.0;
    if (factors.size() == 1) return factors[0];
    if (interval <= 0.0) return factors[0];

    // Adjust time relative to pattern start
    double adjustedTime = t - startTime;
    if (adjustedTime < 0.0) adjustedTime = 0.0;

    // Calculate which interval we're in
    double position = adjustedTime / interval;
    int index1 = static_cast<int>(floor(position));

    // Handle boundary conditions
    int numFactors = static_cast<int>(factors.size());

    // Check if we're past the last defined factor
    if (index1 >= numFactors - 1)
    {
        // Option 1: Return last factor (no wrapping) - typical for valve closure
        return factors[numFactors - 1];

        // Option 2: Wrap around (cyclic pattern) - uncomment if needed
        // index1 = index1 % numFactors;
    }

    int index2 = index1 + 1;

    // Ensure index2 is valid
    if (index2 >= numFactors)
    {
        return factors[numFactors - 1];
    }

    // Calculate interpolation fraction (0.0 to 1.0 within the interval)
    double fraction = position - static_cast<double>(index1);
    fraction = max(0.0, min(fraction, 1.0));

    // Linear interpolation
    double f1 = factors[index1];
    double f2 = factors[index2];

    return f1 + (f2 - f1) * fraction;
}

//  Variable Pattern Constructor/Destructor

VariablePattern::VariablePattern(string name_) :
    Pattern(name_, VARIABLE_PATTERN)
{}

VariablePattern::~VariablePattern()
{
    times.clear();
}

//-----------------------------------------------------------------------------

//  Initializes the state of a VariablePattern.
//  (Variable patterns have no initial offset time.)

void VariablePattern::init(double intrvl, double tstart)
{
    if (factors.size() == 0)
    {
        factors.push_back(1.0);
        times.push_back(0.0);
    }
    currentIndex = 0;
    currentTime = 0.0;
}

//-----------------------------------------------------------------------------

//  Finds the time (sec) until the next change in a VariablePattern.

double VariablePattern::nextTime(double t)
{
    if (currentIndex >= static_cast<int>(times.size()) - 1)
    {
        return numeric_limits<double>::max();
    }
    return times[currentIndex + 1];
}

//-----------------------------------------------------------------------------

//  Advances a VariablePattern to the period associated with time t (sec).

void VariablePattern::advance(double t)
{
    currentTime = t;  // Store current time for interpolation

    for (size_t i = currentIndex + 1; i < times.size(); i++)
    {
        if (t < times[i])
        {
            currentIndex = static_cast<int>(i) - 1;
            return;
        }
    }
    currentIndex = static_cast<int>(times.size()) - 1;
}

double VariablePattern::interpolatedFactor(double t)
{
    if (factors.size() == 0) return 1.0;
    if (factors.size() == 1) return factors[0];
    if (times.size() != factors.size()) return factors[0];

    // Before first time point
    if (t <= times[0]) return factors[0];

    // After last time point
    if (t >= times.back()) return factors.back();

    // Find the interval containing time t
    for (size_t i = 0; i < times.size() - 1; i++)
    {
        if (t >= times[i] && t < times[i + 1])
        {
            // Calculate interpolation fraction
            double t1 = times[i];
            double t2 = times[i + 1];
            double fraction = (t - t1) / (t2 - t1);

            // Linear interpolation
            double f1 = factors[i];
            double f2 = factors[i + 1];

            return f1 + (f2 - f1) * fraction;
        }
    }

    // Fallback (shouldn't reach here)
    return factors.back();
}