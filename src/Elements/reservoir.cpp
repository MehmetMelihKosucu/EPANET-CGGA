/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "reservoir.h"
#include "pattern.h"
#include "Core/network.h"

using namespace std;

//-----------------------------------------------------------------------------

//    Constructor

Reservoir::Reservoir(std::string name_, NodeCategory category_) :
    Node(name_, Node::NodeCategory::IT1),       // Set category to IT1
    headPattern(0)
{
    fullDemand = 0.0;
    fixedGrade = true;
    pastHead = head;
    ph = head;
    category = Node::NodeCategory::IT1;  // Reservoirs are always Independent Type 1
}

Reservoir::~Reservoir() {}

Reservoir* Reservoir::clone() const
{
	Reservoir* newReservoir = new Reservoir(*this);
	newReservoir->headPattern = headPattern;
	newReservoir->qualSource = qualSource;
	newReservoir->actualDemand = actualDemand;
	newReservoir->outflow = outflow;
	newReservoir->fixedGrade = fixedGrade;
	newReservoir->fullDemand = fullDemand;
	newReservoir->pastHead = pastHead;
	newReservoir->ph = ph;
	return newReservoir;
}

Reservoir::Reservoir(const Reservoir& other)
	: Node(other.name, other.category) // Manually call the Node constructor
{
	// Now copy the Reservoir-specific fields
	this->headPattern = other.headPattern;
	this->qualSource = other.qualSource;  // watch out for pointer ownership
	this->actualDemand = other.actualDemand;
	this->outflow = other.outflow;
	this->fixedGrade = other.fixedGrade;
	this->fullDemand = other.fullDemand;
	this->pastHead = other.pastHead;
	this->ph = other.ph;
}

//-----------------------------------------------------------------------------

void Reservoir::convertUnits(Network* nw)
{
    elev /= nw->ucf(Units::LENGTH);
    initQual /= nw->ucf(Units::CONCEN);
}

//-----------------------------------------------------------------------------

void Reservoir::setFixedGrade()
{
    double f = 1.0;
    if ( headPattern )
    {
        f = headPattern->currentFactor();
    }
    head = elev * f;
    fixedGrade = true;
}
