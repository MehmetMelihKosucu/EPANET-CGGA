#include "node.h"
#include "junction.h"
#include "reservoir.h"
#include "tank.h"
#include "qualsource.h"
#include "Utilities/mempool.h"
#include "Core/network.h"
#include <iostream>
#include "link.h"
using namespace std;

//-----------------------------------------------------------------------------

// Constructor

Node::Node(const std::string& name_, NodeCategory cat) :
    Element(name_),
    name(name_),             // Initialize the name
    category(cat)            // Initialize the category
{
    // Initialize other members in the constructor body
    rptFlag = false;
    elev = 0.0;
    xCoord = -1e20;
    yCoord = -1e20;
    initQual = 0.0;
    qualSource = nullptr;
    fixedGrade = false;
    head = 0.0;
    h1ini = 0.0;
    h2ini = 0.0;
    pastHead = 0.0;
    ph = 0.0;
    qGrad = 0.0;
    fullDemand = 0.0;
    actualDemand = 0.0;
    outflow = 0.0;
    quality = 0.0;
    pastOutflow = 0.0;
}

// Destructor

Node::~Node()
{
    delete qualSource;
}

//-----------------------------------------------------------------------------

// Factory Method

Node* Node::factory(int type_, string name_, MemPool* memPool)
{
    Node* newNode = nullptr;

    switch (type_) {
    case JUNCTION:
        newNode = new(memPool->alloc(sizeof(Junction))) Junction(name_, Node::NodeCategory::OTHER);
        break;
    case TANK:
        newNode = new(memPool->alloc(sizeof(Tank))) Tank(name_, NodeCategory::IT1); // Tanks are IT1
        break;
    case RESERVOIR:
        newNode = new(memPool->alloc(sizeof(Reservoir))) Reservoir(name_, NodeCategory::IT1); // Reservoirs are IT1
        break;
    default:
        return nullptr;
    }

    if (!newNode) {
        throw std::runtime_error("Node::factory() failed: Memory allocation returned nullptr");
    }

    return newNode;
}

// Copy constructor
Node::Node(const Node& other)
    : Element(other.name),
    name(other.name),
    category(other.category),
    rptFlag(other.rptFlag),
    elev(other.elev),
    xCoord(other.xCoord),
    yCoord(other.yCoord),
    initQual(other.initQual),
    qualSource(nullptr), // or copy if necessary
    fixedGrade(other.fixedGrade),
    head(other.head),
    h1ini(other.h1ini),
    h2ini(other.h2ini),
    pastHead(other.pastHead),
    ph(other.ph),
    qGrad(other.qGrad),
    fullDemand(other.fullDemand),
    actualDemand(other.actualDemand),
    outflow(other.outflow),
    quality(other.quality)
{
    // If qualSource needs a deep copy, do it here
    // e.g. if (other.qualSource) qualSource = new QualSource(*other.qualSource);
}



//-----------------------------------------------------------------------------

void Node::initialize(Network* nw)
{
    head = elev;
    pastHead = elev;   // head on the previous timestep
    ph = elev;         // synonym of the past head
    h1ini = elev;
    h2ini = elev;
    quality = initQual;

    if (qualSource)
        qualSource->quality = quality;

    actualDemand = 0.0;
    outflow = 0.0;

    if (type() == JUNCTION)
        fixedGrade = false;
    else
        fixedGrade = true;
}


//-----------------------------------------------------------------------------

int Node::getIndex() const
{
    return index; // Assuming `index` is a member variable
}

//-----------------------------------------------------------------------------

// Check if Node is IT1
bool Node::isIndependentType1() const {
    if (this == nullptr) {
        std::cerr << "Error: Attempted to call isIndependentType1() on a null Node pointer.\n";
        return false;  // Default return value if the node is null
    }
    return category == NodeCategory::IT1;
}

//-----------------------------------------------------------------------------

// Check if Node is DT1
bool Node::isDependentType1() const
{
    return category == NodeCategory::DT1;
}

//-----------------------------------------------------------------------------

// Dynamically determine and set node category
void Node::determineCategory(Network* nw) {
    // Tanks and Reservoirs are always Independent Type 1 (IT1)
    if (type() == RESERVOIR || type() == TANK) {
        category = NodeCategory::IT1;
        return;
    }

    // Junctions
    if (type() == JUNCTION) {
        Junction* junction = dynamic_cast<Junction*>(this);
        if (junction) {
            // Check if the demand model is pressure-dependent
            if (nw->isDemandModelPressureDependent() && junction->hasConsumerDemand()) {
                category = NodeCategory::DT1;
                return;
            }
        }

        // Default to IT1 if no pressure-dependent demands
        category = NodeCategory::IT1;
        return;
    }

    // Default to OTHER for any unhandled cases
    category = NodeCategory::OTHER;
}


bool Node::hasNonlinearCharacteristics() const {
    // Ensure hasEmitter is implemented for relevant node types
    return (this->type() == JUNCTION && hasEmitter());
}

