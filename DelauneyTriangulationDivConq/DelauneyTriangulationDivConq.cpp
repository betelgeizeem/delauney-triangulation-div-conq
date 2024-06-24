#include "DelauneyTriangulationDivConq.h"

#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <stdexcept>
#include <algorithm>
#include <fstream>

#include "double2.h"

namespace DelauneyTriangulationDivConq
{
    const double angularEps = 0;//0.0348995; // arcsin(2 grad)

    Box::Box(const std::vector<double2> points)
    {
        minCoord = std::numeric_limits<double>::max();
        maxCoord = std::numeric_limits<double>::min();

        for (auto node : points)
        {
            minCoord = min(node, minCoord);
            maxCoord = max(node, maxCoord);
        }

        size = maxCoord - minCoord;
    }

    Triangulation::EdgePtr Triangulation::findTwin(const Edge& edge) const
    {
        auto node2 = edge->getNode2();
        return node2.lock()->getEdgeToNode(edge->getNode1());
    }

    Triangulation::EdgePtr Triangulation::findConnectingEdge(NodePtr node1, NodePtr node2) const
    {
        return node1.lock()->getEdgeToNode(node2);
    }

    Triangulation::EdgePtr Triangulation::findNext(const Edge& edge) const
    {
        auto triangle = edge->getTriangle().lock();
        for (auto i = 0; i < 3; ++i)
            if (triangle->edges[i].lock() == edge)
                return triangle->edges[(i + 1) % 3];

        return EdgePtr();
    }

    Triangulation::EdgePtr Triangulation::findNextForFreeEdgeTwin(const Edge& freeEdgeTwin) const
    {
        EdgePtr twin, prev, thisEdge = freeEdgeTwin;

        do
        {
            prev = findPrev(thisEdge.lock());
            twin = findTwin(prev.lock());
            thisEdge = twin;
        } while (!(twin.expired() || twin.lock() == freeEdgeTwin));

        return twin.lock() == freeEdgeTwin ? EdgePtr() : prev;
    }

    Triangulation::EdgePtr Triangulation::findPrev(const Edge& edge) const
    {
        auto triangleWPtr = edge->getTriangle();
        for (auto i = 0; i < 3; ++i)
            if (auto triangle = triangleWPtr.lock(); triangle->edges[i].lock() == edge)
                return triangle->edges[(i + 3 - 1) % 3];

        return EdgePtr();
    }

    void Triangulation::flip(const Edge& edge)
    {
        // edge fippes to twinNextNode2 -> edgeNextNode2
        auto twinW = findTwin(edge);
        if (twinW.expired())
        {
            return;
        }
        auto twin = twinW.lock();
        auto twinPrev = findPrev(twin).lock();
        auto twinNext = findNext(twin).lock();

        auto edgePrev = findPrev(edge).lock();
        auto edgeNext = findNext(edge).lock();

        // remove nodes linkage:
        edge->getNode1().lock()->removeEdge(edge);
        twin->getNode1().lock()->removeEdge(twin);

        // change triangles:
        auto triangle1W = edge->getTriangle(); // will become edgePrev, twinNext, edge
        auto triangle2W = twin->getTriangle(); // will become twin, twinPrev, edgeNext

        auto triangle1 = triangle1W.lock(); // will become edgePrev, twinNext, edge
        auto triangle2 = triangle2W.lock(); // will become twin, twinPrev, edgeNext

        triangle1->edges = { edgePrev, twinNext, edge };
        triangle2->edges = { twin, twinPrev, edgeNext };

        // link edges to new triangles:
        for (auto t : { triangle1 , triangle2 })
            for (auto e : t->edges)
                e.lock()->setTriangle(t);

        // change edges:

        (*edge) = EdgeData(twinNext->getNode2(), edgeNext->getNode2(), triangle1W);
        (*twin) = EdgeData(edgeNext->getNode2(), twinNext->getNode2(), triangle2W);

        // relink edges to nodes:
        twinNext->getNode2().lock()->addEdge(edge);
        edgeNext->getNode2().lock()->addEdge(twin);
    }

    bool counterClockWise(double2 p0, double2 p1, double2 p2)
    {
        auto zvec = crossProductZ(p1 - p0, p2 - p0) / length(p1 - p0) / length(p2 - p0);
        return zvec > angularEps;
    }

    void Triangulation::makeConvex()
    {
        // find first outer edge
        Triangulation::Edge startEdgeTwin;
        for (auto edge : getEdges())
        {
            if (findTwin(edge).expired())
            {
                startEdgeTwin = edge;
                break;
            }
        }

        // loop to build new triangles outside of existing triangulation
        // finishes when traversal through all the free edges does not find any new non-convex angles
        auto currEdgeTwin = startEdgeTwin;

        bool createdTriangle = false;
        do
        {
            createdTriangle = false;
            auto nextEdgeTwin = findNextForFreeEdgeTwin(currEdgeTwin).lock();
            auto node0 = currEdgeTwin->getNode1().lock();
            auto node1 = nextEdgeTwin->getNode1().lock();
            auto node2 = currEdgeTwin->getNode2().lock();

            if (counterClockWise(node0->coords, node1->coords, node2->coords))
            {
                createTriangle(node0, node1, node2);
                currEdgeTwin = startEdgeTwin = findConnectingEdge(node1, node2).lock();
                createdTriangle = true;
            }
            else
                currEdgeTwin = nextEdgeTwin;
        } while (createdTriangle || startEdgeTwin != currEdgeTwin);
    }

    Triangulation::EdgePtr Triangulation::findFreeOutcomingEdge(const Node& node) const
    {
        auto outcomingEdges = node->getEdges(); // all outcoming edges
        // need to pick up free one (its twin is expired)
        auto freeIt = std::find_if(outcomingEdges.begin(), outcomingEdges.end(), [this](EdgePtr& edgePtr)
            { return this->isFreeEdge(edgePtr.lock()); });
        if (freeIt == outcomingEdges.end())
            return EdgePtr();

        return *freeIt;
    }

    Triangulation::TrianglePtr Triangulation::createTriangle(Node node1, Node node2, Node node3)
    {
        // create edges
        auto edge1 = std::make_shared<EdgeData>(node1, node2);
        auto edge2 = std::make_shared<EdgeData>(node2, node3);
        auto edge3 = std::make_shared<EdgeData>(node3, node1);

        // put them into triangulation
        edges.insert(edges.end(), { edge1, edge2, edge3 });

        // link points to edges
        node1->addEdge(edge1);
        node2->addEdge(edge2);
        node3->addEdge(edge3);

        // make triangle
        auto triangle = std::make_shared<TriangleData>(edge1, edge2, edge3);

        // remember triangle
        triangles.push_back(triangle);

        // link edges to triangle
        for (auto& edge : { edge1, edge2, edge3 })
            edge->setTriangle(triangle);

        return triangle;
    }

    void Triangulation::removeTriangle(TrianglePtr triangleWPtr)
    {
        auto triangle = triangleWPtr.lock();
        auto it = std::find(triangles.begin(), triangles.end(), triangle);
        triangles.erase(it);

        for (auto edge : triangle->edges)
            removeEdge(edge);
    }

    void Triangulation::removeEdge(EdgePtr edgeWPtr)
    {
        auto edge = edgeWPtr.lock();
        auto it = std::find(edges.begin(), edges.end(), edge);
        edges.erase(it);

        for (auto node : { edge->getNode1(), edge->getNode2() })
            node.lock()->removeEdge(edge);
    }

    bool isOKCounterClockwiseTriangle(double2 p0, double2 p1, double2 p2)
    {
        auto vec1 = p2 - p1;
        auto vec2 = p0 - p1;
        return crossProductZ(vec1, vec2) > 0;
    }

    double findAngle(double2 vec1, double2 vec2)
    {
        auto dot = dotProduct(vec1, vec2);
        auto cos = dot / length(vec1) / length(vec2);

        return std::acos(cos);
    }

    double findAngle(double2 p0, double2 p1, double2 p2)
    {
        return findAngle(p0 - p1, p2 - p1);
    }

    double getSmallestAngleOfTriangle(double2 A, double2 B, double2 C)
    {
        std::array<double2, 3> vecs = { B - A, C - B, A - C };

        auto smallestAngle = std::numeric_limits<double>::max();
        for (auto i = 0; i < 3; ++i)
        {
            auto angle = findAngle(vecs[(i + 1) % 3], -vecs[i]);
            smallestAngle = std::min(smallestAngle, angle);
        }

        return smallestAngle;
    }

    bool isInsideTriangle(double2 p, double2 A, double2 B, double2 C)
    {
        if (crossProductZ(B - A, p - A) < angularEps)
            return false;
        if (crossProductZ(C - B, p - B) < angularEps)
            return false;
        if (crossProductZ(A - C, p - C) < angularEps)
            return false;
        return true;
    }

    Triangulation::Node chooseNewTriangleNodeInMergingZone(
        Triangulation::Node nodeL,
        Triangulation::Node nodeR,
        Triangulation::Node nodeNewL,
        Triangulation::Node nodeNewR)
    {
        auto A = nodeL->coords;
        auto C = nodeR->coords;

        auto nodeNewChosen = Triangulation::Node();

        std::array<double, 2> minAngles;

        for (auto i = 0; i < 2; ++i)
        {
            auto& nodeNew = i == 0 ? nodeNewL : nodeNewR;

            auto B = nodeNew->coords;

            minAngles[i] = isOKCounterClockwiseTriangle(A, B, C) ? getSmallestAngleOfTriangle(A, B, C) : -1;
        }

        if (minAngles[0] < 0 && minAngles[1] < 0)
            return nullptr;

        if (minAngles[0] < 0)
            return nodeNewR;
        if (minAngles[1] < 0)
            return nodeNewL;

        if (isInsideTriangle(nodeNewR->coords, A, nodeNewL->coords, C))
            return nodeNewR;
        if (isInsideTriangle(nodeNewL->coords, A, nodeNewR->coords, C))
            return nodeNewL;

        return minAngles[0] > minAngles[1] ? nodeNewL : nodeNewR;
    }

    void Triangulation::mergeAsConvexStripes(Triangulation& stripe)
    {
        auto thisUpperRightNode = nodes[0];
        auto thisLowerRightNode = nodes[0];

        for (auto& node : nodes)
        {
            if (node->coords.y > thisUpperRightNode->coords.y)
                thisUpperRightNode = node;
            if (node->coords.y == thisUpperRightNode->coords.y &&
                node->coords.x > thisUpperRightNode->coords.x)
                thisUpperRightNode = node;

            if (node->coords.y < thisLowerRightNode->coords.y)
                thisLowerRightNode = node;
            if (node->coords.y == thisLowerRightNode->coords.y &&
                node->coords.x > thisLowerRightNode->coords.x)
                thisLowerRightNode = node;
        }

        auto otherUpperLeftNode = stripe.nodes[0];
        auto otherLowerLeftNode = stripe.nodes[0];

        for (auto& node : stripe.nodes)
        {
            if (node->coords.y > otherUpperLeftNode->coords.y)
                otherUpperLeftNode = node;
            if (node->coords.y == otherUpperLeftNode->coords.y &&
                node->coords.x < otherUpperLeftNode->coords.x)
                otherUpperLeftNode = node;

            if (node->coords.y < otherLowerLeftNode->coords.y)
                otherLowerLeftNode = node;
            if (node->coords.y == otherLowerLeftNode->coords.y &&
                node->coords.x < otherLowerLeftNode->coords.x)
                otherLowerLeftNode = node;
        }

        // merge data
        merge(std::move(stripe));

        auto firstIteration = true;
        auto counter = 0;
        while (thisUpperRightNode != thisLowerRightNode || otherUpperLeftNode != otherLowerLeftNode)
        {
            auto leftNode = thisUpperRightNode;

            // find node candidate on this triangulation
            while (!isOKCounterClockwiseTriangle(
                otherUpperLeftNode->coords,
                thisUpperRightNode->coords,
                leftNode->coords))
            {
                thisUpperRightNode = leftNode;

                if (thisUpperRightNode != thisLowerRightNode)
                {
                    auto outcomingFreeEdgePtr = findFreeOutcomingEdge(leftNode);
                    if (outcomingFreeEdgePtr.expired())
                        throw std::runtime_error("Delauney triangulation: Could not find free outcoming edge for node.");

                    auto incomingFreeEdge = findNextForFreeEdgeTwin(outcomingFreeEdgePtr.lock());
                    if (incomingFreeEdge.expired())
                        throw std::runtime_error("Delauney triangulation: Could not find next free edge.");
                    leftNode = incomingFreeEdge.lock()->getNode1().lock();
                }
                else 
                    break;

                if (!firstIteration)
                    break;
            }

            firstIteration = false;

            // find node candidate on other triangulation
            auto rightNode = otherUpperLeftNode;
            if (otherUpperLeftNode != otherLowerLeftNode)
            {
                auto outcomingFreeEdgePtr = findFreeOutcomingEdge(rightNode);
                if (outcomingFreeEdgePtr.expired())
                    throw std::runtime_error("Delauney triangulation: Could not find free outcoming edge for node.");

                rightNode = outcomingFreeEdgePtr.lock()->getNode2().lock();
            }

            // choose triangle (thisUpperRightNode, thisEdgeCandidate->node1, otherUpperLeftNode) or
            // (thisUpperRightNode, otherEdgeCandidate->node2, otherUpperLeftNode)

            auto chosenNode = chooseNewTriangleNodeInMergingZone(
                thisUpperRightNode, otherUpperLeftNode,
                leftNode, rightNode);

            if (counter++ == 0 && !chosenNode)
            {
                throw std::runtime_error("Delauney triangulation: Could not merge two convex Delauney triangulation stripes.");
            }
            else if (!chosenNode)
            {
                //print();
                break;
            }

            createTriangle(thisUpperRightNode, chosenNode, otherUpperLeftNode);
            (chosenNode == leftNode ? thisUpperRightNode : otherUpperLeftNode) = chosenNode;
        }

        makeConvex();
    }

    void Triangulation::merge(Triangulation&& other)
    {
        nodes.reserve(nodes.size() + other.nodes.size());
        edges.reserve(edges.size() + other.edges.size());
        triangles.reserve(triangles.size() + other.triangles.size());

        std::move(other.nodes.begin(), other.nodes.end(), std::back_inserter(nodes));
        std::move(other.edges.begin(), other.edges.end(), std::back_inserter(edges));
        std::move(other.triangles.begin(), other.triangles.end(), std::back_inserter(triangles));
    }

    bool Triangulation::nodesAreOnSameLine(const std::vector<Triangulation::Node>& nodes)
    {
        auto n = nodes.size();
        auto vector = nodes[0]->coords - nodes[std::max(2, int(n / 2))]->coords;
        const auto& stripe0Coords = nodes[0]->coords;
        bool allNodesAreOnTheSameLine = true;
        for (auto i = 1; i < nodes.size() && allNodesAreOnTheSameLine; ++i)
        {
            allNodesAreOnTheSameLine = allNodesAreOnTheSameLine
                && std::abs(crossProductZ(vector, nodes[i]->coords - stripe0Coords)) < angularEps;
        }

        return allNodesAreOnTheSameLine;
    }

    void Triangulation::print(const std::string& filename)
    {
        if (!OUTPUT)
            return;

        std::ofstream file;
        file.open(filename + ".graph");

        if (!file.is_open())
            throw;

        file << "{ \"x0\": 0, \"y0\" : 0, \"vertices\" : [";

        for (auto i = 0; i < nodes.size(); ++i)
        {
            file << "{";
            file << "\"x\": " << nodes[i]->coords.x * 100 << ",";
            file << "\"y\": " << nodes[i]->coords.y * 100 << ",";
            file << "\"name\" : \"" << nodes[i]->id << "\",";
            file << "\"radius\" : 20,";
            file << "\"background\" : \"#ffffff\",";
            file << "\"fontSize\" : 18,";
            file << "\"color\" : \"#000000\",";
            file << "\"border\" : \"#000000\"";
            file << "}";
            if (i != nodes.size() - 1)
                file << ",";
        }
        file << "],\"edges\": [";

        for (auto i = 0; i < getEdges().size(); ++i)
        {
            auto edge = getEdges()[i];
            auto node1 = std::distance(nodes.begin(), std::find_if(nodes.begin(), nodes.end(),
                [edge](Node node) { return node->id == edge->getNode1().lock()->id; }));
            auto node2 = std::distance(nodes.begin(), std::find_if(nodes.begin(), nodes.end(),
                [edge](Node node) { return node->id == edge->getNode2().lock()->id; }));
            file << "{";
            file << "\"vertex1\": " << node1 << ",";
            file << "\"vertex2\" : " << node2 << ",";
            file << "\"weight\" : \"\",";
            file << "\"isDirected\" : true,";
            file << "\"controlStep\" : 0,";
            file << "\"fontSize\" : 18,";
            file << "\"lineWidth\" : 2,";
            file << "\"background\" : \"#000000\",";
            file << "\"color\" : \"#000000\"";
            file << "}";
            if (i != getEdges().size() - 1)
                file << ",";
        }

        file << "],\"texts\": [] }";

        file.close();
    }

    //void Triangulation::printLinks()
    //{
    //    if (!OUTPUT)
    //        return;
    //    std::ofstream file;
    //    file.open("graph_links.graph");

    //    if (!file.is_open())
    //        throw;

    //    file << "{ \"x0\": 0, \"y0\" : 0, \"vertices\" : [";

    //    for (auto i = 0; i < nodes.size(); ++i)
    //    {
    //        file << "{";
    //        file << "\"x\": " << nodes[i]->coords.x * 100 << ",";
    //        file << "\"y\": " << nodes[i]->coords.y * 100 << ",";
    //        file << "\"name\" : \"" << nodes[i]->id << "\",";
    //        file << "\"radius\" : 20,";
    //        file << "\"background\" : \"#ffffff\",";
    //        file << "\"fontSize\" : 18,";
    //        file << "\"color\" : \"#000000\",";
    //        file << "\"border\" : \"#000000\"";
    //        file << "}";
    //        if (i != nodes.size() - 1)
    //            file << ",";
    //    }
    //    file << "],\"edges\": [";

    //    //Then traverse all its edges to get the longest one:
    //    //auto lowestD = std::numeric_limits<double>::min();
    //    for (auto i = 0; i < nodes.size(); ++i)
    //        for (auto edgeW : nodes[i]->getEdges())
    //        {
    //            auto edge = edgeW.lock();
    //            auto node1 = std::distance(nodes.begin(), std::find_if(nodes.begin(), nodes.end(),
    //                [edge](Node node) { return node->id == edge->getNode1().lock()->id; }));
    //            auto node2 = std::distance(nodes.begin(), std::find_if(nodes.begin(), nodes.end(),
    //                [edge](Node node) { return node->id == edge->getNode2().lock()->id; }));
    //            file << "{";
    //            file << "\"vertex1\": " << node1 << ",";
    //            file << "\"vertex2\" : " << node2 << ",";
    //            file << "\"weight\" : \"\",";
    //            file << "\"isDirected\" : true,";
    //            file << "\"controlStep\" : 0,";
    //            file << "\"fontSize\" : 18,";
    //            file << "\"lineWidth\" : 2,";
    //            file << "\"background\" : \"#000000\",";
    //            file << "\"color\" : \"#000000\"";
    //            file << "}";
    //            if (i != nodes.size() - 1)
    //                file << ",";
    //        }

    //    file << "],\"texts\": [] }";

    //    file.close();
    //}

    Triangulation TriangulationDelauneyBuilder::build(const std::vector<double2>& points)
    {
        // make Box
        box = Box(points);

        //eps = std::sqrt(box.size.x * box.size.y / points.size()) / 100.;

        // all nodes:
        std::vector<Triangulation::Node> nodes(points.size());
        for (auto i = 0; i < points.size(); ++i)
            nodes[i] = std::make_shared<Triangulation::NodeData>(i, points[i]);

        // divide nodes into vertical stripes [SkvortsovAV2002]
        std::vector<Triangulation> triangulationsOnStripes;
        auto m = (int)std::sqrt(0.13f * box.size.x / box.size.y * nodes.size());
        m = std::max(m, 1);
        bool success = false;

        for (auto attempts = 0; attempts < 3 && !success; ++attempts)
        {
            success = false;
            auto stripes = divideIntoStripes(nodes, m * (1 + attempts));
            if (stripes.size() == 1 && isBadStripe(stripes.front()))
            {
                throw std::runtime_error("Delauney triangulation: All points are on a single line, but it was not preprocessed.");
            }

            // run simple build of Delauney triangulation for stripes, then make triangulations convex
            triangulationsOnStripes.resize(stripes.size());
            auto stripeSuccess = true;
            for (auto i = 0; i < stripes.size() && stripeSuccess; ++i)
            {
                triangulationsOnStripes[i] = buildForStripe(stripes[i]);
                stripeSuccess = stripeSuccess && !triangulationsOnStripes[i].empty();

                if (!stripeSuccess)
                    break;

                triangulationsOnStripes[i].makeConvex();

                stripeSuccess = stripeSuccess && flipTrianglesForDelauneyResult(triangulationsOnStripes[i]);
            }

            if (!stripeSuccess)
            {
                for (auto i = 0; i < stripes.size(); ++i)
                    triangulationsOnStripes[i].clear();
                for (auto& node : nodes)
                    node->clear();

                continue;
            }


            // merge Delauney triangulation for stripes
            for (auto i = 1; i < stripes.size(); ++i)
            {
                triangulationsOnStripes[0].mergeAsConvexStripes(triangulationsOnStripes[i]);
                //triangulationsOnStripes[0].print("graph");
            }

            success = true;
        }

        if (success)
        {
            auto& mergedTriangulation = triangulationsOnStripes[0];
            flipTrianglesForDelauneyResult(mergedTriangulation);;
            return mergedTriangulation;
        }
        else
            throw std::runtime_error("Delauney triangulation: Failed to build Delauney triangulation for points.");
    }

    std::vector<std::vector<Triangulation::Node>> TriangulationDelauneyBuilder::divideIntoStripes(const std::vector<Triangulation::Node>& nodes, int m)
    {
        std::vector<std::vector<Triangulation::Node>> stripes(m);

        auto nodesCopy = nodes;

        // sort points by X and Y (if X1 == X2)
        std::sort(nodesCopy.begin(), nodesCopy.end(), [](const Triangulation::Node& a, const Triangulation::Node& b)
            {
                return (a->coords.x != b->coords.x) ? (a->coords.x < b->coords.x) : (a->coords.y < b->coords.y);
            });

        // fill stripes with nodes
        std::size_t n = 0;
        auto widthStripe = box.size.x / m + box.size.x / 1000.;

        for (std::size_t i = 0; i < m - 1; ++i)
        {
            auto b = box.minCoord.x + (i + 1) * widthStripe;
            while (n < nodesCopy.size() && nodesCopy[n]->coords.x <= b)
                stripes[i].push_back(nodesCopy[n++]);
        };

        while (n < nodesCopy.size())
            stripes[m - 1].push_back(nodesCopy[n++]);

        // bad stripe - a stripe with less than 3 points or with points on one line
        auto foundBadStripe = true;
        while (foundBadStripe && m > 1)
        {
            foundBadStripe = false;
            int lastOkStripe = 0;
            for (auto i = 0; i < m; ++i)
            {
                if (isBadStripe(stripes[i]))
                {
                    if (i == 0)
                    {
                        stripes[0].insert(stripes[0].end(), stripes[1].begin(), stripes[1].end());
                        stripes[1].clear();
                        lastOkStripe = 0;
                    }
                    else
                    {
                        stripes[lastOkStripe].insert(stripes[lastOkStripe].end(), stripes[i].begin(), stripes[i].end());
                        stripes[i].clear();
                    }
                    foundBadStripe = true;
                }
                else if (!stripes[i].empty())
                    lastOkStripe = i;
            }

            auto newEnd = std::remove_if(stripes.begin(), stripes.end(), [](const auto& stripe) { return stripe.empty(); });
            stripes.erase(newEnd, stripes.end());
            m = (int)stripes.size();
        }

        return stripes;
    }

    std::array<Triangulation::Node, 3> TriangulationDelauneyBuilder::
        chooseNewTriangle(
            const Triangulation& triangulation,
            std::array<Triangulation::EdgePtr, 3> edges,
            Triangulation::Node newNode) const
    {
        auto bestTriangleSmallestAngle = -1.;
        std::array<Triangulation::Node, 3> newTriangleNodes;
        for (auto i = 0; i < edges.size(); ++i)
        {
            if (edges[i].expired())
                continue;

            auto thisEdge = edges[i].lock();
            auto nextEdge = triangulation.findNext(thisEdge).lock();//edges[(i + 1) % edges.size()].lock();

            auto oppNode = nextEdge->getNode2().lock(); // oppositeNode
            auto baseNode = thisEdge->getNode1().lock();
            auto endNode = thisEdge->getNode2().lock();

            auto vec1 = newNode->coords - baseNode->coords;
            auto vec2 = oppNode->coords - baseNode->coords;
            auto thisEdgeVec = thisEdge->getVector();
            if (crossProductZ(vec1, thisEdgeVec) * crossProductZ(thisEdgeVec, vec2) <= 0) // would be bad triangle
                continue;

            auto A = endNode->coords;
            auto B = baseNode->coords;
            auto C = newNode->coords;

            auto leastAngle = getSmallestAngleOfTriangle(A, B, C);

            if (leastAngle > bestTriangleSmallestAngle)
            {
                newTriangleNodes = { endNode, baseNode, newNode };
                bestTriangleSmallestAngle = leastAngle;
            }
        }

        return newTriangleNodes;
    }

    bool TriangulationDelauneyBuilder::isBadStripe(const std::vector<Triangulation::Node>& stripe)
    {
        // less than 3 points
        if (stripe.size() < 3)
            return true;

        // all nodes in a stripe are on the same line
        return Triangulation::nodesAreOnSameLine(stripe);
    }

    bool onSameLine(double2 p0, double2 p1, double2 p2)
    {
        return std::abs(crossProductZ(p1 - p0, p2 - p0)) < angularEps;
    }

    Triangulation TriangulationDelauneyBuilder::buildForStripe(std::vector<Triangulation::Node> nodes)
    {
        // sort points by Y
        std::sort(nodes.begin(), nodes.end(), [](const Triangulation::Node& a, const Triangulation::Node& b)
            {
                return a->coords.y < b->coords.y;
            });
        Triangulation triangulation(nodes);

        std::unordered_set<int> includedNodesIndices;
        // create initial triangle:
        auto lastFreeEdges = std::array<Triangulation::EdgePtr, 3>();
        auto k = 2;
        while (lastFreeEdges[0].expired() && k < nodes.size())
        {
            if (!onSameLine(nodes[0]->coords, nodes[1]->coords, nodes[k]->coords))
            {
                auto lastCreatedTriangle = Triangulation::TrianglePtr();
                if (counterClockWise(nodes[0]->coords, nodes[1]->coords, nodes[k]->coords))
                    lastCreatedTriangle = triangulation.createTriangle(nodes[0], nodes[1], nodes[k]);
                else
                    lastCreatedTriangle = triangulation.createTriangle(nodes[k], nodes[1], nodes[0]);

                includedNodesIndices.insert({ 0, 1, k });
                lastFreeEdges = lastCreatedTriangle.lock()->edges;
            }
            else ++k;
        }

        // cuold not create initial triangle
        if (lastFreeEdges[0].expired())
            return Triangulation();

        auto startNodeIndex = k == 2 ? 3 : 2;
        // now build triangulation by adding points one by one
        std::vector<Triangulation::TrianglePtr> lastCreatedTriangles(1); //todo
        auto& lastCreatedTriangle = lastCreatedTriangles[0];
        while (includedNodesIndices.size() < nodes.size())
        {
            for (auto i = startNodeIndex; i < nodes.size(); ++i)
            {
                if (includedNodesIndices.count(i) != 0)
                    continue;

                auto triangleNodes = chooseNewTriangle(triangulation, lastFreeEdges, nodes[i]);
                if (!triangleNodes[0] || !triangleNodes[1] || !triangleNodes[2])
                {
                    return Triangulation();
                }
                lastCreatedTriangle = triangulation.createTriangle(triangleNodes[0], triangleNodes[1], triangleNodes[2]);
                auto newEdges = lastCreatedTriangle.lock()->edges;
                for (auto e = 0; e < 3; ++e)
                {
                    if (triangulation.isFreeEdge(newEdges[e].lock()))
                        lastFreeEdges[e] = newEdges[e];
                    else
                        lastFreeEdges[e] = Triangulation::EdgePtr();
                }

                includedNodesIndices.insert(i);
            }
        }

        //correct triangulation by flipping triangles
        if (!flipTrianglesForDelauneyResult(triangulation))
            return Triangulation();

        for (auto t : triangulation.getTriangles())
        {
            auto edges = t->edges;

            for (auto e : t->edges)
            {
                auto twin = triangulation.findTwin(e.lock());

                if (twin.expired())
                    continue;

                auto node = triangulation.findNext(twin.lock()).lock()->getNode2();
                //if (!checked && !isDelauneyCriteriumFor4Points(fourPoints[1], fourPoints[2], fourPoints[3], fourPoints[0]))
                if (/*!checked && */!t->fulfilsDelauneyCriteriumWith(node))
                {

                    triangulation.flip(e.lock());

                    //checked = true;
                    return triangulation;
                }
            }
        }

        return triangulation;
    }

    bool TriangulationDelauneyBuilder::flipTrianglesForDelauneyResult(
        Triangulation& triangulation)
    {
        auto flipped = false;
        auto counter = 0;
        auto const MAX_NUMBER_OF_FLIPPING = triangulation.getTriangles().size() * 2;
        do
        {
            flipped = false;
            for (auto triangle : triangulation.getTriangles())
                for (auto edge : triangle->edges)
                {
                    auto twin = triangulation.findTwin(edge.lock());

                    if (twin.expired())
                        continue;

                    auto notMyNodeInNeibTriangle = triangulation.findNext(twin.lock()).lock()->getNode2();
                    if (!triangle->fulfilsDelauneyCriteriumWith(notMyNodeInNeibTriangle))
                    {
                        triangulation.flip(edge.lock());

                        flipped = true;
                        break;
                    }
                }
            counter++;
        } while (flipped && counter < MAX_NUMBER_OF_FLIPPING);

        if (counter == MAX_NUMBER_OF_FLIPPING)
        {
            return false;
        }

        return true;
    }

    double2 Triangulation::EdgeData::getVector()
    {
        return getNode2().lock()->coords - getNode1().lock()->coords;
    }

    void Triangulation::NodeData::addEdge(EdgePtr edge)
    {
        bool hasEdgeToSameNode = hasEdgeToNode(edge.lock()->getNode2());
        if (!hasEdgeToSameNode)
            edges.push_back(edge);
    }

    void Triangulation::NodeData::removeEdge(EdgePtr edge)
    {
        auto itCandidate = std::find_if(edges.begin(), edges.end(), [edge](auto edgeStored)
            {
                return edgeStored.lock() == edge.lock();
            });
        if (itCandidate != std::end(edges))
            edges.erase(itCandidate);
    }

    Triangulation::EdgePtr Triangulation::NodeData::getEdgeToNode(NodePtr node) const
    {
        auto itCandidate = std::find_if(edges.begin(), edges.end(), [node](EdgePtr candidate)
            {
                return candidate.lock()->getNode2().lock() == node.lock();
            });
        if (itCandidate == std::end(edges))
            return EdgePtr();
        else
            return (*itCandidate);
    }

    bool Triangulation::TriangleData::fulfilsDelauneyCriteriumWith(NodePtr node)
    {
        auto nodesOfTriangle = getNodes();

        std::array<double2, 3> points;
        for (auto i = 0; i < 3; ++i)
            points[i] = nodesOfTriangle[i].lock()->coords;

        double x1 = points[0].x;
        double x2 = points[1].x;
        double x3 = points[2].x;

        double y1 = points[0].y;
        double y2 = points[1].y;
        double y3 = points[2].y;

        auto det2 = [](
            double a11, double a12,
            double a21, double a22)
            {
                return (a11 * a22 - a12 * a21);
            };

        auto det3 = [&det2](
            double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33)
            {
                return
                    a11 * det2(a22, a23, a32, a33) -
                    a12 * det2(a21, a23, a31, a33) +
                    a13 * det2(a21, a22, a31, a32);
            };

        double l1 = x1 * x1 + y1 * y1;
        double l2 = x2 * x2 + y2 * y2;
        double l3 = x3 * x3 + y3 * y3;

        auto a = det3(
            x1, y1, 1,
            x2, y2, 1,
            x3, y3, 1);

        auto b = det3(
            l1, y1, 1,
            l2, y2, 1,
            l3, y3, 1);

        auto c = det3(
            l1, x1, 1,
            l2, x2, 1,
            l3, x3, 1);

        auto d = det3(
            l1, x1, y1,
            l2, x2, y2,
            l3, x3, y3);

        auto x0 = node.lock()->coords.x;
        auto y0 = node.lock()->coords.y;

        auto l0 = x0 * x0 + y0 * y0;
        return a * l0 - b * x0 + c * y0 - d >= 0.;
    }
}