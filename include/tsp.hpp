#ifndef TSP_HPP
#define TSP_HPP

#include "cost_matrix.hpp"
#include "tsp_waypoint_struct.hpp"
#include "viewpoint_struct.hpp"

#include <vector>

class TSP {
    public:

        TSP();
        TSP(CostMatrix cm);
        void loadCM(int vgd, Viewpoint start);
        void greedyInit();
        void twoOpt();
        float pathCost();
        float pathCost(std::vector<TSPWaypoint>& path);
        void getPath(std::vector<std::vector<float>>& path); // return path with view directions
        void reassignModuleMembership(const std::vector<size_t>& module_membership);

    private:

        size_t n_vp;
        CostMatrix cm;
        // std::vector<size_t> path;
        // std::vector<size_t> nodes;
        std::vector<TSPWaypoint> path;
        std::vector<TSPWaypoint> nodes;

        float insertionCost(size_t idx, size_t insert_idx);
        size_t nearestNeighbor(TSPWaypoint wpt, float& best_cost);
        size_t nearest(float& best_cost, size_t& node_idx);

};

#endif // TSP_HPP