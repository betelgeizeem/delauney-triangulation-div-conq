#pragma once

#include "double2.h"

#include <memory>
#include <array>
#include <algorithm>
#include <vector>

static auto OUTPUT = true;

namespace DelauneyTriangulationDivConq
{
    struct Box
    {
        Box(const std::vector<double2> points = {});

        mutable double2 minCoord;
        mutable double2 maxCoord;

        mutable double2 size;
    };

    class Triangulation
    {
    public:
        using IdVertex = std::size_t;

        class NodeData;
        class EdgeData;
        class TriangleData;

        using Node = std::shared_ptr<Triangulation::NodeData>;
        using Edge = std::shared_ptr<Triangulation::EdgeData>;
        using Triangle = std::shared_ptr<Triangulation::TriangleData>;

        using NodePtr = std::weak_ptr<Triangulation::NodeData>;
        using EdgePtr = std::weak_ptr<Triangulation::EdgeData>;
        using TrianglePtr = std::weak_ptr<Triangulation::TriangleData>;

        Triangulation(std::vector<Node> nodes) : nodes(nodes) {};
        Triangulation() {}

        const std::vector<Node>& getNodes() const { return nodes; }
        const std::vector<Edge>& getEdges() const { return edges; }
        const std::vector<Triangle>& getTriangles() const { return triangles; }

        void clear() { nodes.clear(); edges.clear(); triangles.clear(); }
        bool empty() { return edges.size() == 0; }

        EdgePtr findTwin(const Edge& edge) const;

        EdgePtr findConnectingEdge(NodePtr node1, NodePtr node2) const;
        EdgePtr findNext(const Edge& edge) const;
        EdgePtr findNextForFreeEdgeTwin(const Edge& freeEdgeTwin) const;
        EdgePtr findPrev(const Edge& edge) const;

        void flip(const Edge& edge);

        EdgePtr findFreeOutcomingEdge(const Node& node) const;
        bool isFreeEdge(const Edge& edge) const { return findTwin(edge).expired(); }

        TrianglePtr createTriangle(Node node1, Node node2, Node node3);
        void removeTriangle(TrianglePtr triangleWPtr);
        void removeEdge(EdgePtr edgeWPtr);

        void makeConvex();
        // this stripe is left, other is situated to the right of this one. For non-convex triangulations result is undefined
        void mergeAsConvexStripes(Triangulation& stripe);
        void merge(Triangulation&& other);

        static bool nodesAreOnSameLine(const std::vector<Triangulation::Node>& nodes);

        void print(const std::string& filename);
        // void printLinks();
    private:
        std::vector<Node> nodes;
        std::vector<Edge> edges;
        std::vector<Triangle> triangles;

    };

    // This class creates Delauney Triangulation for set of points. 
    class TriangulationDelauneyBuilder
    {
    public:
        Triangulation build(const std::vector<double2>& points);
    private:
        // builds triangulation, adding triangles vertically point-by-point
        Triangulation buildForStripe(std::vector<Triangulation::Node> stripe);
        std::array<Triangulation::Node, 3> chooseNewTriangle(
            const Triangulation& triangulation,
            std::array<Triangulation::EdgePtr, 3> edges,
            Triangulation::Node newNode) const;
        // divide area into m stripes [SkvortsovAV2002]
        std::vector<std::vector<Triangulation::Node>>
            divideIntoStripes(const std::vector<Triangulation::Node>& nodes, int m);
        bool isBadStripe(const std::vector<Triangulation::Node>& stripe);
        bool flipTrianglesForDelauneyResult(Triangulation& triangulation);
    private:
        Box box;
        //float eps = 0.;
    };


    class Triangulation::EdgeData
    {
    public:
        EdgeData()
            : node1(NodePtr())
            , node2(NodePtr())
            , triangle{ TrianglePtr() }
        {}

        EdgeData(NodePtr node1, NodePtr node2, TrianglePtr triangle = TrianglePtr())
            : node1{ node1 }
            , node2{ node2 }
            , triangle{ triangle }
        {}

        NodePtr getNode1or2(int i) const { return i == 1 ? getNode1() : getNode2(); }
        NodePtr getNode1() const { return node1; }
        NodePtr getNode2() const { return node2; }
        NodePtr getOtherNode(NodePtr node) const { return getNode1().lock() == node.lock() ? getNode2() : getNode1(); }
        bool doesEndOn(NodePtr node) const { return getNode1().lock() == node.lock() || getNode2().lock() == node.lock(); }

        double2 getVector();

        void setTriangle(TrianglePtr tr) { triangle = tr; };
        TrianglePtr getTriangle() { return triangle; }
        virtual ~EdgeData() {}
    private:
        Node node1;
        Node node2;

        TrianglePtr triangle;
    };

    class Triangulation::NodeData
    {
    public:
        NodeData(IdVertex id, double2 coords) : id(id), coords(coords) {}

        const IdVertex id;
        const double2 coords;

        const std::vector<EdgePtr>& getEdges() const { return edges; };

        void addEdge(EdgePtr edge);

        void removeEdge(EdgePtr edge);
        bool hasEdgeToNode(NodePtr node) { return !getEdgeToNode(node).expired(); };
        EdgePtr getEdgeToNode(NodePtr node) const;
        void clear() { edges.clear(); }
    private:
        std::vector<EdgePtr> edges;
    };

    class Triangulation::TriangleData
    {
    public:
        TriangleData(EdgePtr edge1, EdgePtr edge2, EdgePtr edge3) : edges({ edge1, edge2, edge3 }) {}
        std::array<EdgePtr, 3> edges;
        std::array<NodePtr, 3> getNodes() { return { edges[0].lock()->getNode1(), edges[0].lock()->getNode2(), edges[1].lock()->getNode2() }; }

        bool fulfilsDelauneyCriteriumWith(NodePtr node);
    };
}