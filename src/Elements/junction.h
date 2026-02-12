/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

//! \file junction.h
//! \brief Describes the Junction class.

#ifndef JUNCTION_H_
#define JUNCTION_H_

#include "Elements/node.h"
#include "Elements/demand.h"
#include "Elements/pipe.h"
#include "Elements/link.h"

#include <string>
#include <list>
#include <unordered_map>

class Network;
class Emitter;
class Pipe;

//! \class Junction
//! \brief A variable head Node with no storage volume.

class Junction: public Node
{
  public:
    int id;  // Add this line to define the 'id' member variable
	int junctionData;
    Junction(std::string name_, NodeCategory category_ = NodeCategory::OTHER);

    ~Junction();

    Junction(const Junction& other);

    Node* clone() const override;

    int type() const override { return JUNCTION; }
    void   convertUnits(Network* nw);
    void   initialize(Network* nw);
    void   findFullDemand(double multiplier, double patternFactor);
    double findActualDemand(Network* nw, double h, double& dqdh);
    double findEmitterFlow(double h, double& dqdh);
    bool   isPressureDeficient(Network* nw);
    bool   hasEmitter() const override { return emitter != nullptr; }
    bool hasConsumerDemand() const;

    Demand            primaryDemand;   //!< primary demand
    std::list<Demand> demands;         //!< collection of additional demands
    double            pMin;            //!< minimum pressure head to have demand (ft)
    double            pFull;           //!< pressure head required for full demand (ft)
    Emitter*          emitter;         //!< emitter object
	double            pastHead;        //!< Head on the previous time step
	double            ph;             //!< synonym of past head
    bool headCalculated = false;

private:
	bool virtual_ = false;  // Flag to indicate if this is a virtual junction
    //static std::unordered_map<std::string, Junction*> virtualJunctions;
};

#endif
