#include "utils.hpp"

#include <string>
#include <vector>

int main() {
    // get all ordered viewpoint files:
    std::vector<std::string> files = {
        "2m_local.csv",
        "2m_global.csv",
        "4m_local.csv",
        "4m_global.csv",
        "8m_local.csv",
        "8m_global.csv",
        "16m_local.csv",
        "16m_global.csv"
    };

    std::string dir = "../data/ordered_viewpoints/";

    for (auto file : files) {
        // compute the cost of the ordered viewpoint file
    }
}