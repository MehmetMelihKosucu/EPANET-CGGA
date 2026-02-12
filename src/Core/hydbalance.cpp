/* EPANET 3.1.1 Pressure Management Extension
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

//////////////////////////////////////////////
// Implementation of the HydBalance class.  //
//////////////////////////////////////////////

// TO DO:
// - compute and report system wide cumulative flow balance

#include "hydbalance.h"
#include "network.h"
#include "Elements/node.h"
#include "Elements/junction.h"
#include "Elements/link.h"
#include "Elements/valve.h"
#include "Elements/pipe.h"
#include "Elements/tank.h"
#include "Core/hydengine.h"
#include "Solvers/rwcggasolver.h"
#include "Models/headlossmodel.h"
#include "constants.h"
#include <utility>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include "dualflowhistory.h"
#include <queue>
//#include "Solvers/cggasolver.h" // Include the CGGASolver header file
using std::vector;
using std::map;
using std::set;
using std::less;
using std::allocator;
using std::pair;
using std::make_pair;

void   findNodeOutflows(double lamda, double dH[], double xQ[], Network* nw);
void   findLeakageFlows(double lamda, double dH[], double xQ[], Network* nw);
double findTotalFlowChange(double lamda, double dQ[], Network* nw);

//-----------------------------------------------------------------------------

//  Evaluate the error in satisfying the conservation of flow and energy
//  equations by an updated set of network heads and flows.
//

double HydBalance::evaluate(
    double lamda,      // step size
    double dH[],       // change in nodal heads
    double dQ[],       // change in link flows
    double xQ[],       // nodal inflow minus outflow
    Network* nw,       // network being analyzed
    double currentTime, 
    double tstep,
    int currentMode)   
{
    // Initialize error tracking metrics
    maxFlowErr = 0.0;
    maxHeadErr = 0.0;
    maxFlowChange = 0.0;
    maxHeadErrLink = -1;
    maxFlowErrNode = -1;
    maxFlowChangeLink = -1;
    
    // Calculate the total number of nodes including all interior junctions
    int regularNodeCount = nw->nodes.size();
    int interiorJunctionCount = 0;
    int totalSegmentCount = 0;
    
    // Count interior junctions from segmented pipes if we're in Water Hammer mode
    if (currentMode == 2) {  // Water Hammer mode
        for (Link* link : nw->links) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe && pipe->isSegmented && pipe->numSegments > 1) {
                // Count interior junctions directly from the pipe
                interiorJunctionCount += pipe->interiorJunctions.size();
                
                // Count segments in this pipe
                totalSegmentCount += pipe->numSegments;
            } else {
                // Regular link counts as one segment
                totalSegmentCount += 1;
            }
        }
        
        int totalNodeCount = regularNodeCount + interiorJunctionCount;
        
        nw->msgLog << "\nHydBalance evaluating with " << totalNodeCount 
                  << " total nodes (" << regularNodeCount << " regular + " 
                  << interiorJunctionCount << " interior) and "
                  << totalSegmentCount << " segments";
                  
        // Ensure the xQ array is large enough for all nodes
        // This should already be handled by the calling function, but we check here as well
        if (totalNodeCount > regularNodeCount) {
            // Make sure we're not accessing beyond array bounds
            // Note: This is just a safety check - actual resizing should happen in the caller
            for (int i = 0; i < totalNodeCount; i++) {
                if (i < regularNodeCount) {
                    // Regular nodes - initialize to zero
                    xQ[i] = 0.0;
                }
            }
        }
    } else {
        // For non-Water Hammer modes, just initialize regular nodes
        for (int i = 0; i < regularNodeCount; i++) {
            xQ[i] = 0.0;
        }
    }
    
    // Different error evaluation based on flow regime
    double norm = 0.0;
    
    // Find head loss errors and update xQ with internal link flows
    norm = findHeadErrorNorm(lamda, dH, dQ, xQ, nw, currentTime, tstep, currentMode);

    // Update xQ with external outflows (demands, emitters, etc.)
    findNodeOutflows(lamda, dH, xQ, nw);

    // Add the flow balance error norm
    norm += findFlowErrorNorm(xQ, nw);
    
    // Calculate total relative flow change across all links
    totalFlowChange = findTotalFlowChange(lamda, dQ, nw);
    
    // Return the root mean square error
    return sqrt(norm);
}

//-----------------------------------------------------------------------------

//  Find the error norm in satisfying the head loss equation across each link.


double HydBalance::findHeadErrorNorm(
        double lamda, double dH[], double dQ[], double xQ[], Network* nw, double currentTime, double tstep, double currentMode)
{
    double norm = 0.0;
    double count = 0.0;
    maxHeadErr = 0.0;
    maxFlowChange = 0.0;
    maxFlowChangeLink = 0;

    int linkCount = nw->links.size();
    for (int i = 0; i < linkCount; i++)
    {

        // ... identify link's end nodes

        Link* link = nw->links[i];
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... apply updated flow to end node flow balances

        double flowChange = lamda * dQ[i];
        double flow = link->flow + flowChange;
        xQ[n1] -= flow;
        xQ[n2] += flow;

        // ... update network's max. flow change

		previousMaxFlowChange = maxFlowChange;

        double err = abs(flowChange);
        if ( err > maxFlowChange )
        {
            maxFlowChange = err;
            maxFlowChangeLink = i;
        } 
        
        //currentMode = determineSolverMode(currentTime, tstep, nw);
		
        // ... compute head loss and its gradient (head loss is saved
        // ... to link->hLoss and its gradient to link->hGrad)
//*******************************************************************
		//if ((currentMode == 0 || currentMode == 1) || currentTime == 0 || link->type() != Link::PIPE || nw->option(Options::HYD_SOLVER) != "CGGA")
        link->findHeadLoss(nw, flow);
		//else
            //pipe->findSegmentHeadLoss(nw, flow);
//*******************************************************************
	    
		double unsteadyTerm = 0;

        // ... evaluate head loss error according to Steady and Unsteady Flow Conditions
		

		if (currentTime == 0 || nw->option(Options::HYD_SOLVER) == "GGA" || currentMode == 0 )
		{
			unsteadyTerm = 0;
		}

		else if (currentMode == 1)
		{
			unsteadyTerm = (link->inertialTerm) * (flow - link->pastFlow) / tstep;
		}
		else if (currentMode == 2)
		{
			unsteadyTerm = (link->inertialTerm) * (flow - link->pastFlow) / tstep;
		}
        h1 = link->fromNode->head + lamda * dH[n1];
        h2 = link->toNode->head + lamda * dH[n2];
        if ( link->hGrad == 0.0 ) link->hLoss = h1 - h2;
        //err = h1 - h2 - link->hLoss;


		err = unsteadyTerm - h1 + h2 + link->hLoss;
       
		if ( abs(err) > maxHeadErr )
        {
            maxHeadErr = abs(err);
            maxHeadErrLink = i;
        }

        // ... update sum of squared errors

        norm += err * err;
        count += 1.0;
    }

    // ... return sum of squared errors normalized by link count

    if ( count == 0.0 ) return 0;
    else return norm / count;

}

//-----------------------------------------------------------------------------

//  Find net external outflow at each network node.

void findNodeOutflows(double lamda, double dH[], double xQ[], Network* nw)
{
    // ... initialize node outflows and their gradients w.r.t. head

    // First validate network pointer
    if (!nw) {
        throw std::runtime_error("Network pointer is null in findNodeOutflows");
    }

    // Initialize outflows and gradients with proper null checks
    for (Node* node : nw->nodes) {
        if (!node) {
            std::cerr << "Warning: Encountered null node during initialization" << std::endl;
            continue;
        }
        node->outflow = 0.0;
        node->qGrad = 0.0;
    }

    // ... find pipe leakage flows & assign them to node outflows

    if ( nw->leakageModel ) findLeakageFlows(lamda, dH, xQ, nw);

    // ... add emitter flows and demands to node outflows

    int nodeCount = nw->nodes.size();
    for (int i = 0; i < nodeCount; i++)
    {
        Node* node = nw->node(i);

        // First, perform thorough null checking
        if (!node) {
            std::cerr << "Error: Null node encountered in findNodeOutflows" << std::endl;
            continue;
        }

        // Safely get the node type once
        
        try {
			int mmk = node->type();
        }
        catch (const std::exception& e) {
            std::cerr << "Error: Failed to get node type at index " << i << ": " << e.what() << std::endl;
            continue;
        }

        // Now perform the type checking
        if (Node::NodeType() != Node::JUNCTION && Node::NodeType() != Node::TANK && Node::NodeType() != Node::RESERVOIR) {
            std::cerr << "Error: Invalid node type at index " << i << std::endl;
            continue;
        }

        double h = node->head + lamda * dH[i];
        double q = 0.0;
        double dqdh = 0.0;

        // ... for junctions, outflow depends on head

        if (node->type() == Node::JUNCTION) {
            // First check if this is an interior node (segment point)
            bool isInteriorNode = false;
            for (Link* link : nw->links) {
                if (link->type() == Link::PIPE) {
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    // Check if this node belongs to pipe segments
                    if (node->name.find(pipe->name + "_seg") != std::string::npos) {
                        isInteriorNode = true;
                        break;
                    }
                }
            }

            if (isInteriorNode) {
                // Special handling for interior nodes
                // Interior nodes don't have emitter or demand flows
                node->qGrad = 0.0;
                node->outflow = 0.0;
                node->actualDemand = 0.0;
                // Interior nodes just pass flow through
                q = 0.0;
                dqdh = 0.0;
            }
            else {
                // Original handling for main junctions
                // Contribution from emitter flow
                q = node->findEmitterFlow(h, dqdh);
                node->qGrad += dqdh;
                node->outflow += q;
                xQ[i] -= q;

                // Contribution from demand flow
                if (node->fixedGrade) {
                    q = xQ[i];
                    xQ[i] -= q;
                }
                else {
                    q = node->findActualDemand(nw, h, dqdh);
                    node->qGrad += dqdh;
                    xQ[i] -= q;
                }
                node->actualDemand = q;
                node->outflow += q;
            }
        }


        // ... for tanks and reservoirs all flow excess becomes outflow

        else
        {
            node->outflow = xQ[i];
            xQ[i] = 0.0;
        }
    }
}

//-----------------------------------------------------------------------------

//  Find the error norm in satisfying flow continuity at each node.

double HydBalance::findFlowErrorNorm(double xQ[], Network* nw)
{
// Note: contributions to the nodal flow imbalance array xQ[] were
//       found previously from findHeadErrorNorm() and findNodeOutflows())

    double norm = 0.0;
    maxFlowErr = 0.0;

    int nodeCount = nw->nodes.size();
    for (int i = 0; i < nodeCount; i++)
    {
        // ... update network's max. flow error

        if ( abs(xQ[i]) > maxFlowErr )
        {
            maxFlowErr = abs(xQ[i]);
            maxFlowErrNode = i;
        }

        // ... update sum of squared errors (flow imbalances)

        norm += xQ[i] * xQ[i];
    }

    // ... return sum of squared errors normalized by number of nodes

    return norm / nodeCount;
}

//-----------------------------------------------------------------------------

//  Assign the leakage flow along each network pipe to its end nodes.

void findLeakageFlows(double lamda, double dH[], double xQ[], Network* nw)
{
    double dqdh = 0.0;  // gradient of leakage outflow w.r.t. pressure head

    for (Link* link : nw->links)
    {
        // ... skip links that don't leak

        link->leakage = 0.0;
        dqdh = 0.0;
        if ( !link->canLeak() ) continue;

        // ... identify link's end nodes and their indexes

        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        int n1 = node1->index;
        int n2 = node2->index;

        // ... no leakage if neither end node is not a junction

        bool canLeak1 = (node1->type() == Node::JUNCTION);
        bool canLeak2 = (node2->type() == Node::JUNCTION);
        if ( !canLeak1 && !canLeak2 ) continue;

        // ... find link's average pressure head

        double h1 = node1->head + lamda * dH[n1] - node1->elev;
        double h2 = node2->head + lamda * dH[n2] - node2->elev;
        double h = (h1 + h2) / 2.0;
        if ( h <= 0.0 ) continue;

        // ... find leakage and its gradient

        link->leakage = link->findLeakage(nw, h, dqdh);

        // ... split leakage flow between end nodes, unless one cannot
        //     support leakage or has negative pressure head

        double q = link->leakage / 2.0;
        if ( h1 * h2 <= 0.0 || canLeak1 * canLeak2 == 0 ) q = 2.0 * q;

        // ... add leakage to each node's outflow

        if ( h1 > 0.0 && canLeak1 )
        {
            node1->outflow += q;
            node1->qGrad += dqdh;
            xQ[n1] -= q;
        }
        if ( h2 > 0.0 && canLeak2 )
        {
            node2->outflow += q;
            node2->qGrad += dqdh;
            xQ[n2] -= q;
        }
    }
}

//-----------------------------------------------------------------------------

//  Find the sum of all link flow changes relative to the sum of all link flows.

double findTotalFlowChange(double lamda, double dQ[], Network* nw)
{
    double qSum = 0.0;
    double dqSum = 0.0;
    double dq;

    for ( int i = 0; i < nw->count(Element::LINK); i++ )
    {

        Link* link = nw->links[i];

		if (link->index == 8)
		{
			int mmk = 36;
		}

        dq = lamda * dQ[i];
        dqSum += abs(dq);
        qSum += abs(link->flow + dq);
    }
    if ( qSum > 0.0 ) return dqSum / qSum;
    else return dqSum;
}

// Add this helper function to check for NaN or Inf
bool HydBalance::isValidValue(double value, const char* name) {
    if (std::isnan(value)) {
        std::cout << "WARNING: " << name << " is NaN!" << std::endl;
        return false;
    }
    if (std::isinf(value)) {
        std::cout << "WARNING: " << name << " is Infinity!" << std::endl;
        return false;
    }
    return true;
}

//  Evaluate the error in satisfying the conservation of flow and energy
//  equations for unsteady (water hammer) analysis.
double HydBalance::evaluateUnsteady(
    double lamda,
    double dH[],
    double dQ[],
    double dQ_start[],
    double dQ_end[],
    double xQ[],
    Network* nw,
    double currentTime,
    double tstep)
{
    // ... initialize which elements have the maximum errors
    maxFlowErr = 0.0;
    maxHeadErr = 0.0;
    maxCharErr = 0.0;
    maxFlowChange = 0.0;
    maxHeadErrLink = -1;
    maxFlowErrNode = -1;
    maxCharErrLink = -1;
    maxFlowChangeLink = -1;

    // ... initialize nodal flow imbalances to 0
    int nodeCount = nw->count(Element::NODE);
    memset(&xQ[0], 0, nodeCount * sizeof(double));

    // ... find the error norm with transition-aware reach-back
    double norm = findUnsteadyHeadErrorNorm(
        lamda, dH, dQ, dQ_start, dQ_end, xQ, nw, currentTime, tstep);

    // ... update xQ with external outflows
    findUnsteadyNodeOutflows(lamda, dH, xQ, nw, tstep);

    // ... add the error norm in satisfying conservation of flow
    norm += findUnsteadyFlowErrorNorm(xQ, nw);

    // ... evaluate the total relative flow change
    totalFlowChange = findUnsteadyTotalFlowChange(lamda, dQ, dQ_start, dQ_end, nw);

    // ... return the root mean square error
    return sqrt(norm);
}

//  Find the error norm in satisfying the characteristic equations and momentum
//  conservation for unsteady flow analysis
double HydBalance::findUnsteadyHeadErrorNorm(
    double lamda, double dH[], double dQ[], double dQ_start[], double dQ_end[],
    double xQ[], Network* nw, double currentTime, double tstep)
{
    double norm = 0.0;
    double count = 0.0;
    maxHeadErr = 0.0;
    maxCharErr = 0.0;
    maxFlowChange = 0.0;
    maxFlowChangeLink = 0;

    int nodeCount = nw->nodes.size();
    int linkCount = nw->count(Element::LINK);

    std::vector<bool> nodeHasIncompressibleLink(nodeCount, false);
    std::vector<int> nodeDomain(nodeCount, -1);

    // Identify nodes connected to incompressible flow links
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];
        if (!node) continue;

        for (Link* link : nw->getLinksConnectedToNode(node)) {
            if (link->fromNode == node || link->toNode == node) {
                if (link->type() == Link::PIPE) {
                    nodeHasIncompressibleLink[i] = true;
                }
            }
        }
    }

    for (int i = 0; i < linkCount; i++)
    {
        Link* link = nw->link(i);
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        double h1 = link->fromNode->head + lamda * dH[n1];
        double h2 = link->toNode->head + lamda * dH[n2];

        if (link->type() == Link::PIPE) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                bool node1IsIntermediate = nodeHasIncompressibleLink[n1];
                bool node2IsIntermediate = nodeHasIncompressibleLink[n2];

                if (node1IsIntermediate) {
                    xQ[n1] -= pipe->startFlow;
                }
                else {
                    xQ[n1] -= pipe->startFlow;
                }

                if (node2IsIntermediate) {
                    xQ[n2] += pipe->endFlow;
                }
                else {
                    xQ[n2] += pipe->endFlow;
                }
            }
            else {
                xQ[n1] -= link->flow;
                xQ[n2] += link->flow;
            }
        }
        else {
            xQ[n1] -= link->flow;
            xQ[n2] += link->flow;
        }

        if (link->type() == Link::PIPE)
        {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue;

            double startFlow = pipe->startFlow + lamda * dQ_start[i];
            double endFlow = pipe->endFlow + lamda * dQ_end[i];

            double c = pipe->getWaveSpeed();
            double A = pipe->getArea();
            double L = pipe->length;
            double g = GRAVITY;
            double Bj = c / (g * A);
            double D = pipe->diameter;
            double epsilon = 0.5;

            double physicalWaveTime = L / c;

            // Default fallback values
            double Qa = pipe->pastStartFlow;
            double Qb = pipe->pastEndFlow;
            double Ha = pipe->fromNode->pastHead;
            double Hb = pipe->toNode->pastHead;

            // *** TRANSITION-AWARE REACH-BACK ***
            FlowHistoryResult historyResult;

            historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, nw);

            if (historyResult.found) {
                Qa = historyResult.startFlow;
                Qb = historyResult.endFlow;
                Ha = historyResult.startHead;
                Hb = historyResult.endHead;
            }

            double KA, KB;

            if (nw->option(Options::HEADLOSS_MODEL) == "D-W") {
                double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, nw->option(Options::KIN_VISCOSITY)));
                double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, nw->option(Options::KIN_VISCOSITY)));
                KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
                KB = frictionFactorB * pipe->length / (2 * g * D * A * A);
            }
            else {
                KA = pipe->resistance;
                KB = pipe->resistance;
            }

            const double MIN_RESISTANCE = 1e-10;
            KA = std::max(KA, MIN_RESISTANCE);
            KB = std::max(KB, MIN_RESISTANCE);

            double kUFa = pipe->computeShearDecayCoeff(Qa);
            double kUFb = pipe->computeShearDecayCoeff(Qb);
            double Ua = 0.0;
            double Ub = 0.0;

            if (kUFa > 0.0 && kUFb > 0.0) {
                auto sign = [](double val) -> double {
                    return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0);
                    };
                Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
                Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
            }

            double BA, BB, RA, RB;

            if (nw->option(Options::HEADLOSS_MODEL) == "D-W") {
                BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;
                BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
            }
            else if (nw->option(Options::HEADLOSS_MODEL) == "H-W") {
                BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
                BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;
            }

            double charErrStart = 0.0;
            double charErrEnd = 0.0;

            if (link->fromNode->type() == Node::RESERVOIR) {
                charErrStart = BA * startFlow - h1 - RA;
            }
            else if (link->fromNode->type() == Node::TANK) {
                charErrStart = BA * startFlow - h1 - RA;
            }
            else {
                charErrStart = BA * startFlow - h1 - RA;
            }

            if (link->toNode->type() == Node::RESERVOIR) {
                charErrEnd = BB * endFlow + h2 - RB;
            }
            else if (link->toNode->type() == Node::TANK) {
                charErrEnd = BB * endFlow + h2 - RB;
            }
            else {
                charErrEnd = BB * endFlow + h2 - RB;
            }

            double flowChangeStart = lamda * dQ_start[i];
            if (abs(flowChangeStart) > maxFlowChange) {
                maxFlowChange = abs(flowChangeStart);
                maxFlowChangeLink = i;
            }

            double flowChangeEnd = lamda * dQ_end[i];
            if (abs(flowChangeEnd) > maxFlowChange) {
                maxFlowChange = abs(flowChangeEnd);
                maxFlowChangeLink = i;
            }

            double maxCharErrLink = std::max(abs(charErrStart), abs(charErrEnd));
            if (maxCharErrLink > maxCharErr) {
                maxCharErr = maxCharErrLink;
            }

			if (charErrStart > 2 || charErrEnd > 2) {
				int mmk = 1409;  // Debugging breakpoint
            }

            norm += charErrStart * charErrStart + charErrEnd * charErrEnd;
            count += 2.0;
        }
        else {
            double flow = link->flow + lamda * dQ[i];

            double flowChange = lamda * dQ[i];
            if (abs(flowChange) > maxFlowChange) {
                maxFlowChange = abs(flowChange);
                maxFlowChangeLink = i;
            }

            link->findHeadLoss(nw, flow);

            if (link->hGrad == 0.0) link->hLoss = h1 - h2;

            double err = -h1 + h2 + link->hLoss;

            if (abs(err) > maxHeadErr) {
                maxHeadErr = abs(err);
                maxHeadErrLink = i;
            }

            norm += err * err;
            count += 1.0;
        }
    }

    if (count == 0.0) return 0;
    else return norm / count;
}


//  Find net external outflow at each network node for unsteady analysis
void HydBalance::findUnsteadyNodeOutflows(double lamda, double dH[], double xQ[], Network* nw, double tstep)
{
    // ... initialize node outflows and their gradients w.r.t. head
    for (Node* node : nw->nodes)
    {
        node->outflow = 0.0;
        node->qGrad = 0.0;
    }

    // ... find pipe leakage flows & assign them to node outflows
    if (nw->leakageModel) findLeakageFlows(lamda, dH, xQ, nw);

    // ... add emitter flows and demands to node outflows
    int nodeCount = nw->count(Element::NODE);
    for (int i = 0; i < nodeCount; i++)
    {
        Node* node = nw->node(i);
        double h = node->head + lamda * dH[i];
        double q = 0.0;
        double dqdh = 0.0;

        if (i == 1292) {
            int mmk = 1409;  // Debugging breakpoint
        }

        // ... for junctions, outflow depends on head
        if (node->type() == Node::JUNCTION)
        {
            // ... contribution from emitter flow
            q = node->findEmitterFlow(h, dqdh);
            node->qGrad += dqdh;
            node->outflow += q;
            xQ[i] -= q;
            double HASTIRLAN;
            // ... contribution from demand flow
            // ... for fixed grade junction, demand is remaining flow excess
            if (node->fixedGrade)
            {
                q = xQ[i];
                xQ[i] -= q;
            }
            // ... otherwise junction has pressure-dependent demand
            else
            {
                
                q = node->findActualDemand(nw, h, dqdh);
                node->qGrad += dqdh;
                xQ[i] -= q;
                HASTIRLAN = xQ[i];
            }
            node->actualDemand = q;
            node->outflow += q;
        }
        // ... for tanks and reservoirs all flow excess becomes outflow
        // (same handling as in the original EPANET function)
        else
        {
            node->outflow = xQ[i];
            xQ[i] = 0.0;

            
        }
    }
}

double HydBalance::findUnsteadyFlowErrorNorm(double xQ[], Network* nw) {
    double norm = 0.0;
    maxFlowErr = 0.0;
    maxFlowErrNode = -1;
    
    // Get basic network information
    int nodeCount = nw->nodes.size();
    
    // Step 1: Identify hydraulic domains created by closed valves
    std::vector<int> nodeDomain(nodeCount, -1);
    
    // Step 2: Identify nodes that are directly connected to closed valves
    // These nodes require special treatment in error assessment
    std::vector<bool> nodeConnectedToClosedValve(nodeCount, false);
    
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];

        if (i == 1308) {
            int mmk = 1409;  // Debugging breakpoint
        }

        if (!node || node->fixedGrade) continue;
        
        // Check all links connected to this node
        for (Link* link : nw->getLinksConnectedToNode(node)) {
            if (link->status == Link::LINK_CLOSED) // || abs(link->flow) < 1e-3) 
            {
                // This node is directly connected to a closed valve
                nodeConnectedToClosedValve[i] = true;
                break; // No need to check further links for this node
            } 
        }
    } // */
    
    // Step 3: Calculate domain-specific errors (excluding closed valve nodes)
    std::map<int, double> domainTotalError;
    std::map<int, int> domainNodeCount;
    std::map<int, int> domainActiveNodeCount; // Nodes not connected to closed valves
    
    int totalActiveNodes = 0; // Counter for nodes included in error calculation
    
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];
        if (!node || node->fixedGrade) continue;
        //if (!node) continue; // Skip null nodes

        int domain = nodeDomain[i];

        if (i == 1308) {
            int mmk = 1409;  // Debugging breakpoint
        }
        
        // Always track basic domain statistics
        domainTotalError[domain] += xQ[i] * xQ[i];
        domainNodeCount[domain]++;
        
        // Skip nodes connected to closed valves for active error calculation
        if (nodeConnectedToClosedValve[i]) {
            // Log this exclusion for transparency
            if (abs(xQ[i]) > 0.1) {
                nw->msgLog << "\nExcluding node " << node->name 
                          << " from error calculation (connected to closed valve, error=" 
                          << abs(xQ[i]) << ")";
            }
            continue; // Skip this node in active error calculations
        } // */
        
        // Count active nodes for this domain
        domainActiveNodeCount[domain]++;
        totalActiveNodes++;
        
        // Track global maximum error (only among active nodes)
        if (abs(xQ[i]) > maxFlowErr) {
            maxFlowErr = abs(xQ[i]);
            maxFlowErrNode = i;
        }
        
        // Add to overall error norm (only active nodes)
        norm += xQ[i] * xQ[i];
    }
    
    // Step 4: Provide comprehensive diagnostic information
    if (maxFlowErr > 0.1) {
        nw->msgLog << "\n--- Flow Error Analysis by Hydraulic Domain ---";
        
        for (auto& pair : domainTotalError) {
            int domain = pair.first;
            int totalNodes = domainNodeCount[domain];
            int activeNodes = domainActiveNodeCount[domain];
            
            if (activeNodes > 0) {
                double avgDomainErrorActive = pair.second / activeNodes;
                
                nw->msgLog << "\nDomain " << domain << ": " << activeNodes 
                          << " active nodes (out of " << totalNodes << " total)";
                
                // Only report significant errors in active nodes
                if (avgDomainErrorActive > 0.01) {
                    nw->msgLog << "  Average continuity error: " << sqrt(avgDomainErrorActive);
                }
            } else {
                nw->msgLog << "\nDomain " << domain << ": All " << totalNodes 
                          << " nodes connected to closed valves (excluded from analysis)";
            }
        }
        
        nw->msgLog << "\nTotal active nodes in error calculation: " << totalActiveNodes 
                  << " out of " << nodeCount;
    }
    
    // Step 5: Return normalized error based only on hydraulically active nodes
    // This gives a more accurate picture of simulation health
    if (totalActiveNodes > 0) {
        return norm / totalActiveNodes;
    } else {
        // Special case: if all nodes are connected to closed valves
        nw->msgLog << "\nWarning: All nodes connected to closed valves - cannot assess flow error";
        return 0.0;
    }
}

//  Find the sum of all link flow changes relative to the sum of all link flows
//  for unsteady analysis, considering both start and end flows
double HydBalance::findUnsteadyTotalFlowChange(double lamda, double dQ[], double dQ_start[], 
                                             double dQ_end[], Network* nw)
{
    double qSum = 0.0;
    double dqSum = 0.0;
    double dq;

    for (int i = 0; i < nw->count(Element::LINK); i++)
    {
        Link* link = nw->link(i);
        
        if (link->type() == Link::PIPE)
        {
            // For compressible flow, consider both start and end flows
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue;
            
            double dq_s = lamda * dQ_start[i];
            double dq_e = lamda * dQ_end[i];
            
            dqSum += abs(dq_s) + abs(dq_e);
            qSum += abs(pipe->startFlow + dq_s) + abs(pipe->endFlow + dq_e);
        }
        else
        {
            // For incompressible flow
            dq = lamda * dQ[i];
            dqSum += abs(dq);
            qSum += abs(link->flow + dq);
        }
    }
    
    if (qSum > 0.0) return dqSum / qSum;
    else return dqSum;


} 



