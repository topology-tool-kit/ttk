#pragma once

#include <vector>
#include <set>

// ttk common includes
#include <Debug.h>

namespace ttk {

    class BranchMappingDistance : virtual public Debug {

        private:
        int baseMetric = 0;
        inline float editCost_Wasserstein(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2);
        inline float editCost_Persistence(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2);

        public:
        BranchMappingDistance() {
            this->setDebugMsgPrefix(
                "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                                    // the beginning of every msg
        }
        ~BranchMappingDistance() override = default;

        void setBaseMetric(int m) {
            baseMetric = m;
        }

        float editDistance_branch(  std::vector<float> &nodes1,
                                    std::vector<std::vector<int>> &topo1,
                                    int rootID1,
                                    std::vector<float> &nodes2,
                                    std::vector<std::vector<int>> &topo2,
                                    int rootID2);
    };

}