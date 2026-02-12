/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file pipe.h
//! \brief Describes the Pipe class.

#ifndef PIPE_H_
#define PIPE_H_

#include "Elements/link.h"
#include "Elements/junction.h"

class Network;
class Node;
class Junction;

struct TimestepInitialState {
    double flow;            ///< Pipe flow at timestep start
    double startFlow;       ///< Pipe start flow at timestep start (WH mode)
    double endFlow;         ///< Pipe end flow at timestep start (WH mode)
    double fromHead;        ///< Upstream node head at timestep start
    double toHead;          ///< Downstream node head at timestep start
    double timestamp;       ///< Time when this state was stored
    bool isValid;           ///< Flag indicating if state is valid

    TimestepInitialState()
        : flow(0.0), startFlow(0.0), endFlow(0.0)
        , fromHead(0.0), toHead(0.0)
        , timestamp(0.0), isValid(false) {
    }

    void clear() {
        flow = startFlow = endFlow = fromHead = toHead = timestamp = 0.0;
        isValid = false;
    }
};

//! \class Pipe
//! \brief A circular conduit Link through which water flows.

class Pipe: public Link
{
  public:

    // Constructor/Destructor

    Pipe(std::string name);

    ~Pipe();

    // Override clone() to handle Pipe-specific properties.
    Pipe(const Pipe& other);
    Pipe* clone() const override;
    

    // Methods
    int         type() const override { return Link::PIPE; }
    std::string typeStr() override { return "Pipe"; }
    void        convertUnits(Network* nw);
    bool        isReactive();
    void        setInitFlow();
    void        setInitStatus(int s);
    void        setInitSetting(double s);
    void        setResistance(Network* nw);
    void		setLossFactor();
    double      getRe(const double q, const double viscos);
    double      getResistance() {return resistance;}
    double      getVelocity();
    double      getUnitHeadLoss();
    double      getSetting(Network* nw) { return roughness; }
    double      getVolume() { return 0.785398 * length * diameter * diameter; }
    double inertiaCoeff = 0.0;
    double getArea() const 
    {
        // π × diameter²/4
        if (diameter <= 0.0) {
            // Return a small but non-zero area for invalid diameters
            return 0.001;  // Some minimum value to prevent division by zero
        }
        return 0.7853981634 * diameter * diameter; // π/4 = 0.7853981634
    }

    double      getWaveSpeed() const;
    void        findHeadLoss(Network* nw, double q);
   
    bool        canLeak() { return leakCoeff1 > 0.0; }
    double      findLeakage(Network* nw, double h, double& dqdh);
    bool        changeStatus(int s, bool makeChange,
                            const std::string reason,
                            std::ostream& msgLog);
    void        validateStatus(Network* nw, double qTol);

    // Properties
    bool   hasCheckValve;    //!< true if pipe has a check valve
    double length;           //!< pipe length (ft)
    double roughness;        //!< roughness parameter (units depend on head loss model)
    double resistance;       //!< resistance factor (units depend head loss model)
    double lossFactor;       //!< minor loss factor (ft/cfs^2)
    double leakCoeff1;       //!< leakage coefficient (user units)
    double leakCoeff2;       //!< leakage coefficient (user units)
    double bulkCoeff;        //!< bulk reaction coefficient (mass^n/sec)
    double wallCoeff;        //!< wall reaction coefficient (mass^n/sec)
    double massTransCoeff;   //!< mass transfer coefficient (mass^n/sec)
	double diameter;         //!< pipe diameter (ft)
	double waveSpeed;       //!< Wave propagation speed
    
    double startFlow;      // Flow at upstream end (NodeA)
    double pastStartFlow;  // Previous startFlow
    double endFlow;        // Flow at downstream end (NodeB)
    double pastEndFlow;    // Previous endFlow
   
    double getFrictionFactor(double Re) const;

    int numSegments;                       // Number of segments.
    // Boolean function to indicate if this pipe is a segment.
    bool isSegmentPipe() const { return _isSegment; }
    void setIsSegment(bool seg) { _isSegment = seg; }
    bool isSegmented = false;  // Marks if the pipe has already been subdivided
    // --- Segmentation Data Members ---
    struct Segment {
        double length;
        double diameter;
        double roughness;
        // Replace single flow with start and end flows
        double startFlow;      // Flow at upstream end (NodeA)
        double endFlow;        // Flow at downstream end (NodeB)
        double pastStartFlow;  // Previous startFlow
        double pastEndFlow;    // Previous endFlow
        double headLoss;
        double resistance;
        double dx;
        double area;
        double waveSpeed;

        double startFlowChange; //<-- Added missing type for startFlowChange
        double endFlowChange;   //<-- Added missing type for endFlowChange
        double startDistance;
        double endDistance;
        Node* NodeA;
        Node* NodeB;
        double junctionFlow;  // Flow at the junction (Qp)
        double startHead;     // Head at segment start
        double endHead;       // Head at segment end

    };
    Pipe* segment; 

    double requiredSegmentLength; 
    double courantNumber; 

    // Segment storage
    std::vector<Junction*> interiorJunctions;  // Store interior junction pointers
    std::vector<Segment> segments;             // Store segment properties

    std::vector<double> segmentHeads;      // Heads at segment boundaries.
    std::vector<double> segmentFlows;      // Flows at segment boundaries.
    std::vector<double> pastSegmentHeads;  // Past heads.
    std::vector<double> pastSegmentFlows;  // Past flows.

    const std::vector<Segment>& getSegments() const { return segments; }
    double computeShearDecayCoeff(double flow) const;

    double getLength() const;
    bool _isSegment; // Boolean to indicate if this pipe is a segment.

    double Pipe::calcSwameeJain(double Re) const;

    // Characteristic method (MOC) properties
    double charC_plus;     // Positive characteristic value (C+)
    double charC_minus;    // Negative characteristic value (C-)
    double charB_plus;     // Positive characteristic coefficient (B+)
    double charB_minus;    // Negative characteristic coefficient (B-)
    int _numDiscretizedReaches;  // Number of theoretical reaches for Courant condition

    int numReaches; // Number of reaches
    int reachBackSteps; // Number of reach-back steps
    int impingingWaveCount = 0; // Number of impinging waves
    double waveAttenuation; // Wave attenuation factor
    std::vector<std::vector<double>>    wavePathStorage;
    //std::vector<WaveAttenuationPoint> wavePathStorage;
    double characteristicB;
    double resistanceCoeff;

    double startFlowCoeff;     // Coefficient of head term for startFlow calculation
    double startFlowConstant;  // Constant term for startFlow calculation  
    double endFlowCoeff;       // Coefficient of head term for endFlow calculation
    double endFlowConstant;    // Constant term for endFlow calculation

    TimestepInitialState timestepInitial;

  private:
    
    
 };

 


#endif
