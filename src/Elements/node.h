/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

 //! \file node.h
 //! \brief Describes the Node class.

#ifndef NODE_H_
#define NODE_H_

#include "Elements/element.h"
#include <string>

// Forward declarations
class Network;
class Emitter;
class QualSource;
class MemPool;

//! \class Node
//! \brief A connection point between links in a network.
//!
//! A Node is an abstract class that represents the end connection points
//! between links in a pipe network. Specific concrete sub-classes of Node
//! include Junction, Tank, and Reservoir.

class Node : public Element {
public:
    // Enums
    enum NodeType { JUNCTION, TANK, RESERVOIR };
    enum class NodeCategory { IT1, DT1, OTHER }; // New Node Categories

    NodeCategory category; //!< Node category (IT1, DT1, OTHER)

	void setCategory(NodeCategory cat) { category = cat; } // Setter for category
	NodeCategory getCategory() const { return category; } // Getter for category

    // Constructor and Destructor
    Node(const std::string& name_, NodeCategory cat);  // Just declaration

    virtual ~Node(); // {}// = default;

    // Virtual clone method for deep-copying the node.
    virtual Node* clone() const = 0;

    // Copy constructor
    Node(const Node& other);  // Delete the copy constructor

    virtual int type() const { return -1; } // or make it pure virtual

	std::string name; //!< Node name
    // Factory Method
    static Node* factory(int type_, std::string name_, MemPool* memPool);

    // Pure virtual methods to be overridden by subclasses
    //virtual int type() const = 0;
    virtual void convertUnits(Network* nw) = 0;
    virtual void initialize(Network* nw);

    // Methods for Junction nodes
    virtual void findFullDemand(double multiplier, double patternFactor) {}
    virtual double findActualDemand(Network* nw, double h, double& dqdh) { return 0.0; }
    virtual double findEmitterFlow(double h, double& dqdh) { return 0.0; }
    virtual void setFixedGrade() { fixedGrade = false; }
    virtual bool isPressureDeficient(Network* nw) { return false; }
    virtual bool hasEmitter() const { return false; } // Default for base class

    // Methods for Tank nodes
    virtual void validate(Network* nw) {}
    virtual bool isReactive() { return false; }
    virtual bool isFull() { return false; }
    virtual bool isEmpty() { return false; }
    virtual bool isClosed(double flow) { return false; }
    virtual double getVolume() { return 0.0; }

    // Public utility methods
  
    int getIndex() const;
    bool isIndependentType1() const;
    bool isDependentType1() const;
    void determineCategory(Network* nw);
    bool hasNonlinearCharacteristics() const;
    
    // Input Parameters
    bool rptFlag;       //!< True if results are reported
    double elev;        //!< Elevation (ft)
    double xCoord;      //!< X-coordinate
    double yCoord;      //!< Y-coordinate
    double initQual;    //!< Initial water quality concentration
    QualSource* qualSource; //!< Water quality source information

    // Computed Variables
    bool fixedGrade;    //!< Fixed grade status
    double head;        //!< Hydraulic head (ft)
    double h1ini;       //!< Hydraulic head of upstream node in the initial moment of iteration
    double h2ini;       //!< Hydraulic head of downstream node in the initial moment of iteration
    double pastHead;    //!< Hydraulic head in previous timestep
    double ph;          //!< Synonym for past head
    double qGrad;       //!< Gradient of outflow w.r.t. head (cfs/ft)
    double fullDemand;  //!< Full demand required (cfs)
    double actualDemand; //!< Actual demand delivered (cfs)
    double outflow;     //!< Demand + emitter + leakage flow (cfs)
    double pastOutflow; //!< Outflow in previous timestep (cfs)
    double quality;     //!< Water quality concentration (mass/ft³)
    double flow;    // Flow value
	double pastFlow;    // Past flow value // */
    Network* getNetwork() const { return network; }
    
protected :
	//int index;          //!< Index in the network's node list
	Network* network;   //!< Pointer to the network object
    
};

#endif
