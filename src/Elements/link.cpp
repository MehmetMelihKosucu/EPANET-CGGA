/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "link.h"
#include "pipe.h"
#include "pump.h"
#include "valve.h"
#include "Core/constants.h"
#include "Core/network.h"
#include "Utilities/mempool.h"
#include "Core/datamanager.h"
#include "epanet3.h"
#include "node.h"
#include "junction.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using namespace std;

//-----------------------------------------------------------------------------

static const string s_From = " status changed from ";
static const string s_To =   " to ";
static const string linkStatusWords[] = {"CLOSED", "OPEN", "ACTIVE", "TEMP_CLOSED"};

//-----------------------------------------------------------------------------

/// Constructor

Link::Link(std::string name_) :
    Element(name_),
    rptFlag(false),
    fromNode(nullptr),
    toNode(nullptr),
    initStatus(LINK_OPEN),
    diameter(0.0),
    lossCoeff(0.0),
    initSetting(1.0),
    status(0),
    previousStatus(0),
    flow(0.0),
    pastFlow(0.0),
    leakage(0.0),
    hLoss(0.0),
    pastHloss(0.0),
    hGrad(0.0),
    setting(0.0),
    quality(0.0),
    inertialTerm(0.0)
    //category(Node::NodeCategory::OTHER) // Adjusted to include proper namespace
{
}

/// Destructor

Link::~Link() = default;

Link::Link(const Link& other)
    : Element(other) // If Element has a copy constructor we can call
{
    // Copy all Link data members
    this->rptFlag = other.rptFlag;
    this->initStatus = other.initStatus;
    this->diameter = other.diameter;
    this->lossCoeff = other.lossCoeff;
    this->initSetting = other.initSetting;
    this->status = other.status;
    this->previousStatus = other.previousStatus;
    this->flow = other.flow;
    this->pastFlow = other.pastFlow;
    this->leakage = other.leakage;
    this->hLoss = other.hLoss;
    this->pastHloss = other.pastHloss;
    this->hGrad = other.hGrad;
    this->setting = other.setting;
    this->quality = other.quality;
    this->inertialTerm = other.inertialTerm;
    this->phiA = other.phiA;
    this->phiR = other.phiR;
    this->area = other.area;
    this->wallThickness = other.wallThickness;
    this->typeCoefficient = other.typeCoefficient;

    // Pointers are shallow-copied; usually the network owns them
    this->fromNode = other.fromNode;
    this->toNode = other.toNode;
    this->network = other.network;
    this->memPool = other.memPool;

    // Also copy index, if that’s used
    this->index = other.index;
}
//-----------------------------------------------------------------------------

/// Factory Method

Link* Link::factory(int type, std::string name_, MemPool* memPool, Network* network)
{
    Link* link = nullptr;
    switch (type)
    {
    case Link::PIPE:
        link = new Pipe(name_);
        break;
    case Link::PUMP:
        link = new Pump(name_);
        break;
    case Link::VALVE:
        link = new Valve(name_);
        break;
    default:
        throw std::runtime_error("Unknown link type in Link::factory().");
    }
    // Assign pointers if needed
    link->memPool = memPool;
    link->network = network;
    return link;
}


Link* Link::clone() const
{
    // Because this is a const method, we can only call
    // 'this->type()' if it's also const. If not, you can 
    // cast away const or change the design slightly.
    switch (type())
    {
    case Link::PIPE:
        // Downcast to Pipe and call Pipe::clone()
        return dynamic_cast<const Pipe*>(this)->clone();

    case Link::PUMP:
        // If you later implement a Pump class, do:
        return dynamic_cast<const Pump*>(this)->clone();
        throw std::runtime_error("Pump clone() not implemented.");

    case Link::VALVE:
        // If you later implement a Valve class, do:
        return dynamic_cast<const Valve*>(this)->clone();
        throw std::runtime_error("Valve clone() not implemented.");

    default:
        throw std::runtime_error("Unknown link type in Link::clone().");
    }
}

//-----------------------------------------------------------------------------

void Link::initialize(bool reInitFlow)
{
    status = initStatus;
    setting = initSetting;
    if ( reInitFlow )
    {
        if ( status == LINK_CLOSED ) flow = ZERO_FLOW;
        else setInitFlow();
    }
    leakage = 0.0;
	
    inertialTerm = 0; // length / (GRAVITY * area);

    // Assign a default wave speed for the link (ft/s)
    // You can use a typical value for water in pipes (~4,000 ft/s) or assign it dynamically
    //waveSpeed = 4000.0;
}

//-----------------------------------------------------------------------------

double Link::getUnitHeadLoss()
{
    return hLoss;
}

//-----------------------------------------------------------------------------

string Link::writeStatusChange(int oldStatus)
{
    stringstream ss;
    ss << "    " << typeStr() << " " <<
	    name << s_From << linkStatusWords[oldStatus] << s_To <<
	    linkStatusWords[status];
    return ss.str();
}

int Link::getIndex() const {
    return index;  // Ensure `index` is properly initialized for each Link object.
}

Network* Link::getNetwork() const {
    return network;  // Ensure `network` is set during the Link object's initialization.
}











