/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "pipe.h"
#include "Core/network.h"
#include "Core/constants.h"
#include "Models/headlossmodel.h"
#include "Models/leakagemodel.h"
#include "junction.h"
#include <cstdlib>

#include <iostream>
#include <cmath>
using namespace std;

//-----------------------------------------------------------------------------

Pipe::Pipe(std::string name)
      : Link(name), hasCheckValve(false), length(0), roughness(0),
        resistance(0), lossFactor(0), leakCoeff1(0), leakCoeff2(0),
        bulkCoeff(0), wallCoeff(0), massTransCoeff(0), waveSpeed(0),
        numSegments(1), isSegmented(false), _isSegment(false)
{ 
    // Initialize area based on diameter
    diameter = 0.0; // Default initialization
    area = 0.0;     // Default initialization
} //*/


Pipe::~Pipe() {}

//-----------------------------------------------------------------------------

// Copy constructor
Pipe::Pipe(const Pipe& other)
    : Link(other)  // calls the Link copy constructor
{
    // Copy all Pipe-specific fields
    this->hasCheckValve = other.hasCheckValve;
    this->length = other.length;
    this->roughness = other.roughness;
    this->resistance = other.resistance;
    this->lossFactor = other.lossFactor;
    this->leakCoeff1 = other.leakCoeff1;
    this->leakCoeff2 = other.leakCoeff2;
    this->bulkCoeff = other.bulkCoeff;
    this->wallCoeff = other.wallCoeff;
    this->massTransCoeff = other.massTransCoeff;
    this->waveSpeed = other.waveSpeed;
    this->numSegments = other.numSegments;
    this->segment = other.segment;  // shallow copy
    this->isSegmented = other.isSegmented;
    this->_isSegment = other._isSegment;

    // Copy vector members (shallow copy of pointers within them)
    this->segments = other.segments;
    this->interiorJunctions = other.interiorJunctions;
    this->segmentHeads = other.segmentHeads;
    this->segmentFlows = other.segmentFlows;
    this->pastSegmentHeads = other.pastSegmentHeads;
    this->pastSegmentFlows = other.pastSegmentFlows;
}

Pipe* Pipe::clone() const {
    // Create a new Pipe using the copy constructor.
    Pipe* newPipe = new Pipe(name);
    // Copy base class members:
    newPipe->diameter = diameter;
    newPipe->phiA = phiA;
    newPipe->phiR = phiR;
    newPipe->fromNode = fromNode;
    newPipe->toNode = toNode;
    newPipe->fromNode->index = fromNode->index;
    newPipe->toNode->index = toNode->index;
    newPipe->flow = flow;
    newPipe->index = index;
    newPipe->status = status;
    // Copy dual-flow members:
    newPipe->startFlow = startFlow;
    newPipe->endFlow = endFlow;
    newPipe->pastStartFlow = pastStartFlow;
    newPipe->pastEndFlow = pastEndFlow;
    // Copy Pipe-specific members:
    newPipe->length = length;
    newPipe->roughness = roughness;
    newPipe->isSegmented = isSegmented;
    newPipe->_isSegment = _isSegment;
    return newPipe;
} // */


void Pipe::convertUnits(Network* nw)
{
    diameter /= nw->ucf(Units::DIAMETER);
    length   /= nw->ucf(Units::LENGTH);

    // ... convert minor loss coeff. from V^2/2g basis to Q^2 basis

    lossFactor = 0.02517 * lossCoeff / pow(diameter, 4);

    // ... convert roughness length units of Darcy-Weisbach headloss model
    //     (millifeet or millimeters to feet)

    if ( nw->option(Options::HEADLOSS_MODEL ) == "D-W")
    {
        roughness = roughness / nw->ucf(Units::LENGTH) / 1000.0;
    }

    // ... apply global default leakage coeffs.

    if ( leakCoeff1 == MISSING ) leakCoeff1 = nw->option(Options::LEAKAGE_COEFF1);
    if ( leakCoeff2 == MISSING ) leakCoeff2 = nw->option(Options::LEAKAGE_COEFF2);

    // ... apply global default reaction coeffs.

    if ( bulkCoeff == MISSING ) bulkCoeff = nw->option(Options::BULK_COEFF);
    if ( wallCoeff == MISSING ) wallCoeff = nw->option(Options::WALL_COEFF);
}

//-----------------------------------------------------------------------------

bool Pipe::isReactive()
{
    if ( bulkCoeff != 0.0 ) return true;
    if ( wallCoeff != 0.0 ) return true;
    return false;
}

//-----------------------------------------------------------------------------

void Pipe::setResistance(Network* nw)
{
    nw->headLossModel->setResistance(this);
}

//-----------------------------------------------------------------------------

void Pipe::setInitStatus(int s)
{
    initStatus = s;
}

//-----------------------------------------------------------------------------

void Pipe::setInitSetting(double s)
{
    if ( s == 0.0 ) initStatus = LINK_CLOSED;
    else            initStatus = LINK_OPEN;
}

//-----------------------------------------------------------------------------

void Pipe::setInitFlow()
{
    // ... flow at velocity of 1 ft/s
    double initialFlow = PI * diameter * diameter / 4.0;
    startFlow = PI * diameter * diameter / 4.0;
    endFlow = PI * diameter * diameter / 4.0;
    pastStartFlow = PI * diameter * diameter / 4.0;
    pastEndFlow = PI * diameter * diameter / 4.0;
    // Set overall flow as average of start and end flows
    flow = PI * diameter * diameter / 4.0;
    pastFlow = 0;
    pastHloss = 0;
}

//-----------------------------------------------------------------------------

void Pipe::setLossFactor()
{
    lossFactor = 0.02517 * lossCoeff / pow(diameter, 4);
}

double Pipe::getVelocity()
{
    // Use average of startFlow and endFlow for velocity calculation
    double avgFlow = (startFlow + endFlow) / 2.0;
    double area = PI * diameter * diameter / 4.0;
    return abs(flow) / area;
}

//-----------------------------------------------------------------------------

double Pipe::getUnitHeadLoss()
{
    if ( length > 0.0 ) return abs(hLoss) * 1000.0 / length;
    return 0.0;
}

//-----------------------------------------------------------------------------

double Pipe::getRe(const double q, const double viscos)
{
    return abs(q) / (PI * diameter * diameter / 4.0) * diameter / viscos;
}

//-----------------------------------------------------------------------------

void Pipe::findHeadLoss(Network* nw, double q) {
    // Use average of start and end flows for head loss calculation
    //double avgFlow = (startFlow + endFlow) / 2.0;
    
    // Handle closed links
    if (status == LINK_CLOSED || status == TEMP_CLOSED) {
        HeadLossModel::findClosedHeadLoss(q, hLoss, hGrad);
        inertialTerm = MIN_GRADIENT;
        return;
    }
    else {
        // Standard head loss calculation (no segmentation)
        nw->headLossModel->findHeadLoss(this, q, hLoss, hGrad);
        if (hasCheckValve) HeadLossModel::addCVHeadLoss(q, hLoss, hGrad);
        inertialTerm = (length * 4) / (32.174 * PI * diameter * diameter);
    }
    
    // Update overall flow as average of start and end flows
    //flow = avgFlow;
}


//-----------------------------------------------------------------------------

double Pipe::findLeakage(Network* nw, double h, double& dqdh)
{
    return nw->leakageModel->findFlow(leakCoeff1, leakCoeff2, length, h, dqdh);
}


//-----------------------------------------------------------------------------

bool Pipe::changeStatus(int s, bool makeChange, const string reason, ostream& msgLog)
{
    if ( status != s )
    {
        if ( makeChange )
        {
            msgLog << "\n    " << reason;
            status = s;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------

// For debugging only

void Pipe::validateStatus(Network* nw, double qTol)
{
    if ( hasCheckValve && flow < -qTol )
    {
        nw->msgLog << "\nCV " << name << " flow = " << flow*nw->ucf(Units::FLOW);
    }
}

double Pipe::computeShearDecayCoeff(double flow) const {
    const double D = diameter;
    const double A = PI * D * D / 4.0;
    const double nu = 1.08e-5;  // Kinematic viscosity (ft²/s)

    double Re = fabs(flow) * D / (A * nu);

    // Minimum effective value to ensure some damping always occurs
    const double k_min = 0.00000001; // 0.06;

    if (Re < 2300) {
        // Laminar flow: theoretical value is reasonable
        double C_star = 4.76e-3;
        return std::max(sqrt(C_star), k_min);  // ≈ 0.069
    }
    else {
        // Turbulent flow: use empirical multiplier on theoretical value
        // The factor of 3.0 accounts for real-world effects not in theory
        double C_star = 7.41 * pow(Re, -log10(14.3 * pow(Re, -0.05)));
        double k_theoretical = sqrt(C_star);
        return std::max(3.0 * k_theoretical, k_min);
    }
}

/*double Pipe::getWaveSpeed() const {
    // Dalga hizi (m/s). CGGA_WAVE_SPEED ortam degiskeni tanimliysa onu,
    // degilse varsayilan degeri kullanir. Her kosu ayri bir surec oldugu
    // icin deger surecin basinda bir kez okunur ve kullanilan deger
    // dogrulama icin stderr'e bir kez yazilir.
    static const double a_ms = [] {
        const double VARSAYILAN = 319.0;   // m/s
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
        const char* s = std::getenv("CGGA_WAVE_SPEED");
        double v = s ? std::atof(s) : 0.0;
        if (v <= 0.0) v = VARSAYILAN;
        std::cerr << "[CGGA] wave speed = " << v << " m/s\n";
        return v;
        }();
    return a_ms / 0.3048;   // motor ft/s ile calisiyor
} // */

double Pipe::getWaveSpeed() const {
    // Constants in FPS (Imperial) units
    double E = 3.05e8;    // Young's modulus of ductile iron (psi)
    double rho = 1.94;    // Water density (slugs/ft^3)
    double K = 3.19e8;    // Bulk modulus of water (psi)
    double D = diameter;  // Outer diameter in feet

    // Wall thickness for 6-inch ductile iron pipe (Pressure Class 150)
    double e = 0.007 * diameter; // Wall thickness in feet

    // Safety check to prevent division by zero
    if (e <= 0.0) {
        throw std::invalid_argument("Wall thickness must be greater than zero.");
    }

    // Type coefficient (C), accounting for Poisson's ratio (? = 0.25)
    double typeCoefficient = 1;

    double waveSpeed = 0;

    /*if (index == 4) {
        typeCoefficient = 12;
    }
    else if (index == 7 || index == 8) {
        typeCoefficient = 1;
    }
    else{
        typeCoefficient = 1;
    } // */

    /*if (index == 4) {
        waveSpeed = 1000 / 0.3048;
    }
    else if (index == 7 || index == 8) {
        waveSpeed = 1000 / 0.3048;
    }
    else{
        waveSpeed = 1000 / 0.3048;
    } // */

    /*if (diameter <= 0.4 / 0.3048) {
        waveSpeed = 350 / 0.3048;
    }
    else{
        waveSpeed = 1050 / 0.3048;
    } // */

    // Compute wave speed (a) during water hammer
    waveSpeed = 1000 / 0.3048; // 319 / 0.3048; // sqrt((K / rho) / (1 + (K / E) * (D / e) * typeCoefficient));

    return waveSpeed;  // Returns wave speed in ft/s // */
} // */

/*double Pipe::getWaveSpeed() const { // ASCE
    // --- Properties in SI units ---
    // Ductile iron pipe
    const double E = 1.65e11;   // Young's modulus, Pa  (~165 GPa)
    // Water at ~15-20 C
    const double K = 2.19e9;    // Bulk modulus, Pa     (~2.19 GPa)
    const double rho = 998.0;     // Density, kg/m^3

    // Wall-thickness-to-diameter ratio for ductile iron (AWWA C150/C151
    // gives roughly 0.015-0.025 over this size range; 0.02 is representative).
    // D/e is dimensionless, so the diameter's unit does not affect the result.
    // Replace with the actual pressure-class thickness per pipe if known.
    const double eOverD = 0.02;
    if (eOverD <= 0.0)
        throw std::invalid_argument("Wall-thickness ratio must be > 0.");
    const double DoverE = 1.0 / eOverD;   // = 50

    // Pipe restraint coefficient (Korteweg c1):
    //   1.0            -> thin wall, frequent expansion joints (no axial restraint)
    //   1 - nu^2 (~0.92, nu=0.28 for ductile iron) -> anchored against axial motion
    const double c1 = 1.0;

    // Korteweg wave speed in SI, m/s
    const double a_ms = std::sqrt((K / rho) / (1.0 + (K / E) * DoverE * c1));

    // Engine expects wave speed in ft/s (matches the original 1000/0.3048 usage).
    return a_ms / 0.3048;   // m/s -> ft/s   (~1148 m/s -> ~3768 ft/s)
} // */

/*double Pipe::getWaveSpeed() const { // Turkish Journal of Hydraulic
    // --- HDPE (PE100) pressure pipe, SI units ---
    // Short-term (instantaneous) dynamic modulus for a rapid transient.
    // NOT the long-term creep modulus (~0.2-0.3 GPa); a fast event sees the
    // stiff, short-term response. Representative range 0.8-1.5 GPa.
    const double E = 1.1e9;     // Young's modulus, Pa  (~1.1 GPa, short-term PE100)
    const double K = 2.19e9;    // Water bulk modulus, Pa (~2.19 GPa, 15-20 C)
    const double rho = 998.0;     // Water density, kg/m^3

    // HDPE is specified by SDR = outside diameter / wall thickness.
    // SDR 11 (PN16) and SDR 17 (PN10) are typical. e/D = 1/SDR, so D/e = SDR.
    const double SDR = 11.0;      // set to the actual pressure class used
    if (SDR <= 0.0)
        throw std::invalid_argument("SDR must be > 0.");
    const double DoverE = SDR;    // D/e = SDR

    // Restraint coefficient (Korteweg c1). For thin-wall thermoplastic with
    // expansion joints, 1.0; if axially restrained, 1 - nu^2 (nu ~ 0.45 for HDPE).
    const double c1 = 1.0;

    // Korteweg wave speed, m/s
    const double a_ms = std::sqrt((K / rho) / (1.0 + (K / E) * DoverE * c1));

    // Engine expects ft/s.
    return a_ms / 0.3048;   // ~310 m/s for these values
} // */

double Pipe::getFrictionFactor(double Re) const {
    // For zero or extremely low Reynolds numbers, avoid division by zero
    if (Re < 1.0) {
        return 64.0;  // Limiting case as Re approaches zero
    }

    // For laminar flow (Re < 2000)
    if (Re < 2000.0) {
        return 64.0 / Re;
    }

    // For transitional flow (2000 = Re < 4000)
    else if (Re < 4000.0) {
        // Smooth transition between laminar and turbulent regimes
        double fLaminar = 64.0 / 2000.0;
        double fTurbulent = calcSwameeJain(4000.0);
        double t = (Re - 2000.0) / 2000.0;  // Transition factor
        return fLaminar * (1.0 - t) + fTurbulent * t;
    }

    // For turbulent flow (Re = 4000)
    else {
        return calcSwameeJain(Re);
    }
}

// Helper method using the Swamee-Jain approximation formula
double Pipe::calcSwameeJain(double Re) const {
    // Relative roughness (e/D)
    double relativeRoughness = roughness / diameter;

    // Swamee-Jain approximation formula
    return 0.25 / pow(log10(relativeRoughness / 3.7 + 5.74 / pow(Re, 0.9)), 2);
}

// Getter for pipe length
double Pipe::getLength() const {
    return length;
}






