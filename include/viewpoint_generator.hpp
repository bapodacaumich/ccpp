#ifndef VIEWPOINT_GENERATOR_HPP
#define VIEWPOINT_GENERATOR_HPP

#include "cone_camera.hpp"
#include "obs.hpp"
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
            ConeCamera cam=ConeCamera(),
            float inc_angle_max=90.0f * M_PI / 180.0f,
            float inc_improvement_minimum=60.0f * M_PI / 180.0f,
            float inc_improvement_threshold=5.0f * M_PI / 180.0f
        );

        std::vector<Viewpoint> getCoverageViewpoints();
        void printCoverageMap();

        // save coverage map so we don't have to regenerate
        void saveCoverageMap(std::string& filename);
        void loadCoverageMap(std::string& filename);

        // need to create map of coverage for each viewpoint
        void populateCoverage();
        void printIncidenceAngles();
        void getFilteredCoverage(std::vector<bool>& filtered_coverage_data);
        void missedCoverage();

    private:

        float vgd;
        float inc_angle_max;
        float inc_improvement_minimum;
        float inc_improvement_threshold;
        ConeCamera cam;
        size_t num_mesh_faces;

        std::vector<OBS> structure;
        std::vector<Triangle*> all_faces;
        std::vector<Viewpoint> unfiltered_viewpoints;
        std::vector<Viewpoint> coverage_viewpoints;
        std::vector<bool> filtered_coverage; // coverage of 'coverage_viewpoints'
        std::vector<float> filtered_inc_angles; // best inc_angles of 'coverage_viewpoints'
        std::vector<VP_Coverage_Gain> vpcg_unfiltered;
        std::vector<VP_Coverage_Gain> vpcg_filtered;

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
        void sortUpdateMarginalGain();

        void greedy(); // greedy algorithm to select viewpoints and put in coverage_viewpoints


        // use coverage map and update filtered_coverage
        void updateBestIncAngles();
        void updateCoverage(float inc_angle_threshold);
        void setUpCoverageGain();
};

#endif // VIEWPOINT_GENERATOR_HPP