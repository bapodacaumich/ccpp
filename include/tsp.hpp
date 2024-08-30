#ifndef TSP_HPP
#define TSP_HPP

#include "cost_matrix.hpp"
#include <vector>

class TSP {
    public:

        TSP();
        TSP(CostMatrix cm);
        void loadCM(int vgd, vec3 start);
        void greedyInit();
        void twoOpt();
        float pathCost();
        float pathCost(std::vector<size_t>& path);
        void getPath(std::vector<std::vector<float>>& path); // return path with view directions

    private:

        size_t n_vp;
        CostMatrix cm;
        std::vector<size_t> path;
        std::vector<size_t> nodes;

        float insertionCost(size_t idx, size_t insert_idx);
        size_t nearestNeighbor(size_t idx, float& best_cost);
        size_t nearest(float& best_cost, size_t& node_idx);

};

#endif // TSP_HPP