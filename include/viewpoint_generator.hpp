#ifndef VIEWPOINT_GENERATOR_HPP
#define VIEWPOINT_GENERATOR_HPP

#include "obs.hpp"
#include "triangle_coverage_struct.hpp"
#include "viewpoint_coverage_gain_struct.hpp"
#include "viewpoint_struct.hpp"

#include <vector>

class ViewpointGenerator {
    public:

        // generate all viewopints when initialized
        ViewpointGenerator(); // default constructor -- no structure, no viewpoint generation
        // ViewpointGenerator( std::vector<OBS> structure);
        // ViewpointGenerator(
        //     std::vector<OBS> structure,
        //     float vgd
        // );
        ViewpointGenerator(
            std::vector<OBS> structure,
            float vgd=2.0f,
            float inc_angle_max=90.0f * M_PI / 180.0f,
            float inc_improvement_minimum=70.0f * M_PI / 180.0f,
            float inc_improvement_threshold=10.0f * M_PI / 180.0f
        );

        std::vector<Viewpoint> getCoverageViewpoints(bool local);
        void printCoverageMap();

        // save coverage map so we don't have to regenerate
        void saveCoverageMap(std::string& filename);
        void loadCoverageMap(std::string& filename);

        // need to create map of coverage for each viewpoint
        void populateCoverage();
        void printIncidenceAngles();
        void getFilteredCoverage(std::vector<bool>& filtered_coverage_data);
        void missedCoverage();
        void assignModuleMembership(std::vector<Viewpoint>& viewpoints);
        void remapModuleMembership();
        void saveUnfilteredViewpoints(std::string& filename);

    private:

        float vgd;
        float inc_angle_max;
        float inc_improvement_minimum;
        float inc_improvement_threshold;
        size_t num_mesh_faces;

        std::vector<OBS> structure;
        std::vector<Triangle*> all_faces;
        std::vector<Viewpoint> unfiltered_viewpoints;
        std::vector<Viewpoint> coverage_viewpoints;
        std::vector<bool> filtered_coverage; // coverage of 'coverage_viewpoints'
        std::vector<float> filtered_inc_angles; // best inc_angles of 'coverage_viewpoints'
        std::vector<VPCoverageGain> vpcg_unfiltered;
        std::vector<VPCoverageGain> vpcg_filtered;
        std::vector<TriangleCoverage> triangle_coverage;

        // coverage map
        std::vector<std::vector<bool>> coverage_map; // viewpoints x num_mesh_faces
        std::vector<std::vector<float>> inc_angle_map; // viewpoints x num_mesh_faces

        // initialize for constructor
        void initialize();
        bool collision(Viewpoint vp);
        bool collision(vec3 pose, const std::vector<vec3*>& points);

        bool populateViewpoints();
        void rotateViewpoint(const Viewpoint& vp, std::vector<Viewpoint>& rotated_viewpoints, float angle);
        void countMeshFaces();
        void sortUpdateMarginalGain(size_t module_idx=-1);

        void greedy(); // greedy algorithm to select viewpoints and put in coverage_viewpoints
        void greedyModule(size_t module_idx); // greedy algorithm to select viewpoints and put in coverage_viewpoints based on coverage per module  (local constraint)

        // use coverage map and update filtered_coverage
        void updateBestIncAngles();
        void updateCoverage(float inc_angle_threshold);
        void setUpCoverageGain();

        // reassign module membership to four modules
        void reassignModuleMembership();
};

#endif // VIEWPOINT_GENERATOR_HPP