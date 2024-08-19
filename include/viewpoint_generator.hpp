#ifndef VIEWPOINT_GENERATOR_HPP
#define VIEWPOINT_GENERATOR_HPP

#include "cone_camera.hpp"
#include "obs.hpp"
#include "viewpoint_struct.hpp"

#include <vector>

class ViewpointGenerator {
    public:

        // generate all viewopints when initialized
        ViewpointGenerator(); // default constructor -- no structure, no viewpoint generation
        ViewpointGenerator( std::vector<OBS> structure);
        ViewpointGenerator(
            std::vector<OBS> structure,
            float vgd
        );
        ViewpointGenerator(
            std::vector<OBS> structure,
            float vgd,
            ConeCamera cam
        );

        ~ViewpointGenerator(); // destructor for coverage map cleanup

        std::vector<Viewpoint> getCoverageViewpoints();

    private:

        float vgd;
        ConeCamera cam;
        size_t num_mesh_faces;
        std::vector<OBS> structure;
        std::vector<Triangle*> all_faces;
        std::vector<Viewpoint> unfiltered_viewpoints;
        std::vector<Viewpoint> coverage_viewpoints;

        // coverage map
        bool** coverage_map; // viewpoints x num_mesh_faces

        // initialize for constructor
        void initialize();
        bool collision(Viewpoint vp);
        bool collision(vec3 pose, const std::vector<vec3*>& points);

        bool populateViewpoints();
        void countMeshFaces();
        void updateMarginalGainPtrs(
            const bool* filtered_coverage,
            std::vector<int*> sorted_marginal_gain
            );

        void greedy(); // greedy algorithm to select viewpoints and put in coverage_viewpoints

        // need to create map of coverage for each viewpoint
        void populateCoverage();
};

#endif // VIEWPOINT_GENERATOR_HPP