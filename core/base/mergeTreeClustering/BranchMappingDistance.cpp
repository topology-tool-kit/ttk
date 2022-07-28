#include "BranchMappingDistance.h"
#include "munkres.hpp"

#include <iostream>
#include <vector>
#include <stack>
#include <tuple>
#include <limits>
#include <cmath>
#include <float.h>
#include <set>
#include <algorithm>
#include <chrono>

inline float ttk::BranchMappingDistance::editCost_Wasserstein1(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2){
    if(n1<0){
        float b1 = nodes2[n2];
        float d1 = nodes2[p2];
        float b2 = (b1+d1)*0.5;
        float d2 = (b1+d1)*0.5;
        return std::abs(b1-b2)+std::abs(d1-d2);
    }
    if(n2<0){
        float b1 = nodes1[n1];
        float d1 = nodes1[p1];
        float b2 = (b1+d1)*0.5;
        float d2 = (b1+d1)*0.5;
        return std::abs(b1-b2)+std::abs(d1-d2);
    }
    float b1 = nodes1[n1];
    float d1 = nodes1[p1];
    float b2 = nodes2[n2];
    float d2 = nodes2[p2];
    float d = std::abs(b1-b2)+std::abs(d1-d2);
    return d*d;
}

inline float ttk::BranchMappingDistance::editCost_Wasserstein2(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2){
    if(n1<0){
        float b1 = nodes2[n2];
        float d1 = nodes2[p2];
        float b2 = (b1+d1)*0.5;
        float d2 = (b1+d1)*0.5;
        return std::sqrt(std::abs(b1-b2)*std::abs(b1-b2)+std::abs(d1-d2)*std::abs(d1-d2));
    }
    if(n2<0){
        float b1 = nodes1[n1];
        float d1 = nodes1[p1];
        float b2 = (b1+d1)*0.5;
        float d2 = (b1+d1)*0.5;
        return std::sqrt(std::abs(b1-b2)*std::abs(b1-b2)+std::abs(d1-d2)*std::abs(d1-d2));
    }
    float b1 = nodes1[n1];
    float d1 = nodes1[p1];
    float b2 = nodes2[n2];
    float d2 = nodes2[p2];
    float d = std::sqrt(std::abs(b1-b2)*std::abs(b1-b2)+std::abs(d1-d2)*std::abs(d1-d2));
    return d*d;
}

inline float ttk::BranchMappingDistance::editCost_Persistence(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2){
    if(n1<0){
        float b1 = nodes2[n2];
        float d1 = nodes2[p2];
        return std::abs(d1-b1);
    }
    if(n2<0){
        float b1 = nodes1[n1];
        float d1 = nodes1[p1];
        return std::abs(d1-b1);
    }
    float b1 = nodes1[n1];
    float d1 = nodes1[p1];
    float b2 = nodes2[n2];
    float d2 = nodes2[p2];
    float d = std::abs(std::abs(b1-d1)-std::abs(b2-d2));
    return d;
}

inline float ttk::BranchMappingDistance::editCost_Shifting(int n1, int p1, int n2, int p2, std::vector<float> &nodes1, std::vector<float> &nodes2){
    if(n1<0){
        float b1 = nodes2[n2];
        float d1 = nodes2[p2];
        return std::abs(d1-b1);
    }
    if(n2<0){
        float b1 = nodes1[n1];
        float d1 = nodes1[p1];
        return std::abs(d1-b1);
    }
    float b1 = nodes1[n1];
    float d1 = nodes1[p1];
    float b2 = nodes2[n2];
    float d2 = nodes2[p2];
    float pers1 = std::abs(b1-d1);
    float pers2 = std::abs(b2-d2);
    float d = std::abs(b1-b2)+std::abs(pers1-pers2);
    return d;
}

float ttk::BranchMappingDistance::editDistance_branch(   std::vector<float> &nodes1,
                                                        std::vector<std::vector<int>> &topo1,
                                                        int rootID1,
                                                        std::vector<float> &nodes2,
                                                        std::vector<std::vector<int>> &topo2,
                                                        int rootID2){

    // initialize memoization tables

    std::vector<std::vector<int>> predecessors1(nodes1.size());
    std::vector<std::vector<int>> predecessors2(nodes2.size());

    int depth1=0;
    int depth2=0;
    std::stack<int> stack;
    stack.push(rootID1);
    while(!stack.empty()){
        int nIdx = stack.top();
        stack.pop();
        depth1 = std::max((int)predecessors1[nIdx].size(),depth1);
        for(int cIdx : topo1[nIdx]){
            stack.push(cIdx);
            predecessors1[cIdx].reserve(predecessors1[nIdx].size()+1);
            predecessors1[cIdx].insert(predecessors1[cIdx].end(),predecessors1[nIdx].begin(),predecessors1[nIdx].end());
            predecessors1[cIdx].push_back(nIdx);
        }
    }
    stack.push(rootID2);
    while(!stack.empty()){
        int nIdx = stack.top();
        stack.pop();
        depth2 = std::max((int)predecessors2[nIdx].size(),depth2);
        for(int cIdx : topo2[nIdx]){
            stack.push(cIdx);
            predecessors2[cIdx].reserve(predecessors2[nIdx].size()+1);
            predecessors2[cIdx].insert(predecessors2[cIdx].end(),predecessors2[nIdx].begin(),predecessors2[nIdx].end());
            predecessors2[cIdx].push_back(nIdx);
        }
    }

    size_t nn1 = nodes1.size();
    size_t nn2 = nodes2.size();
    size_t dim1 = 1;
    size_t dim2 = (nn1 + 1) * dim1;
    size_t dim3 = (depth1 + 1) * dim2;
    size_t dim4 = (nn2 + 1) * dim3;
    // size_t dim1_ = 1;
    // size_t dim2_ = (nn1 + 1) * dim1;
    // size_t dim3_ = (nn1 + 1) * dim2;
    // size_t dim4_ = (nn2 + 1) * dim3;

    std::unique_ptr<float[]> memT(new float[(nn1+1)*(depth1+1)*(nn2+1)*(depth2+1)]);
    //std::vector<float> memT((nn1+1)*(depth1+1)*(nn2+1)*(depth2+1));

    memT[nn1+0*dim2+nn2*dim3+0*dim4] = 0;
    for(size_t i=0; i<nn1; i++){
        for(size_t l=1; l<=predecessors1[i].size(); l++){

            int curr1 = i;
            int parent1 = predecessors1[i][predecessors1[i].size()-l];

            //-----------------------------------------------------------------------
            // If first subtree has only one branch, return deletion cost of this branch
            if(topo1[curr1].size()==0){
                memT[curr1+l*dim2+nn2*dim3+0*dim4] =
                      this->baseMetric == 0 ?   editCost_Wasserstein1(curr1,parent1,-1,-1,nodes1,nodes2)
                    : this->baseMetric == 1 ?   editCost_Wasserstein2(curr1,parent1,-1,-1,nodes1,nodes2)
                    : this->baseMetric == 2 ?   editCost_Persistence(curr1,parent1,-1,-1,nodes1,nodes2)
                    :                           editCost_Shifting(curr1,parent1,-1,-1,nodes1,nodes2);
            }
            //-----------------------------------------------------------------------
            // If first subtree has more than one branch, try all decompositions
            else{
                float c = FLT_MAX;
                for(auto child1_mb : topo1[curr1]){
                    float c_ = memT[child1_mb+(l+1)*dim2+nn2*dim3+0*dim4];
                    for(auto child1 : topo1[curr1]){
                        if(child1==child1_mb){
                            continue;
                        }
                        c_ += memT[child1+1*dim2+nn2*dim3+0*dim4];
                    }
                    c = std::min(c,c_);
                }
                memT[curr1+l*dim2+nn2*dim3+0*dim4] = c;
            }
            
        }
    }
    for(size_t j=0; j<nn2; j++){
        for(size_t l=1; l<=predecessors2[j].size(); l++){

            int curr2 = j;
            int parent2 = predecessors2[j][predecessors2[j].size()-l];

            //-----------------------------------------------------------------------
            // If first subtree has only one branch, return deletion cost of this branch
            if(topo2[curr2].size()==0){
                memT[nn1+0*dim2+curr2*dim3+l*dim4] =
                  this->baseMetric == 0 ?   editCost_Wasserstein1(-1,-1,curr2,parent2,nodes1,nodes2)
                : this->baseMetric == 1 ?   editCost_Wasserstein2(-1,-1,curr2,parent2,nodes1,nodes2)
                : this->baseMetric == 2 ?   editCost_Persistence(-1,-1,curr2,parent2,nodes1,nodes2)
                :                           editCost_Shifting(-1,-1,curr2,parent2,nodes1,nodes2);
            }
            //-----------------------------------------------------------------------
            // If first subtree has more than one branch, try all decompositions
            else{
                float c = FLT_MAX;
                for(auto child2_mb : topo2[curr2]){
                    float c_ = memT[nn1+0*dim2+child2_mb*dim3+(l+1)*dim4];
                    for(auto child2 : topo2[curr2]){
                        if(child2==child2_mb){
                            continue;
                        }
                        c_ += memT[nn1+0*dim2+child2*dim3+1*dim4];
                    }
                    c = std::min(c,c_);
                }
                memT[nn1+0*dim2+curr2*dim3+l*dim4] = c;
            }

        }
    }

    for(size_t i=0; i<nn1; i++){
        for(size_t j=0; j<nn2; j++){
            for(size_t l1=1; l1<=predecessors1[i].size(); l1++){
                for(size_t l2=1; l2<=predecessors2[j].size(); l2++){

                    int curr1 = i;
                    int parent1 = predecessors1[i][predecessors1[i].size()-l1];

                    int curr2 = j;
                    int parent2 = predecessors2[j][predecessors2[j].size()-l2];

                    //===============================================================================
                    // If both trees not empty, find optimal edit operation
                    
                    //---------------------------------------------------------------------------
                    // If both trees only have one branch, return edit cost between the two branches
                    if(topo1[curr1].size()==0 and topo2[curr2].size()==0){
                        memT[curr1+l1*dim2+curr2*dim3+l2*dim4] =
                              this->baseMetric == 0 ?   editCost_Wasserstein1(curr1,parent1,curr2,parent2,nodes1,nodes2)
                            : this->baseMetric == 1 ?   editCost_Wasserstein2(curr1,parent1,curr2,parent2,nodes1,nodes2)
                            : this->baseMetric == 2 ?   editCost_Persistence(curr1,parent1,curr2,parent2,nodes1,nodes2)
                            :                           editCost_Shifting(curr1,parent1,curr2,parent2,nodes1,nodes2);
                    }
                    //---------------------------------------------------------------------------
                    // If first tree only has one branch, try all decompositions of second tree
                    else if(topo1[curr1].size()==0){
                        float d = FLT_MAX;
                        for(auto child2_mb : topo2[curr2]){
                            float d_ = memT[curr1+l1*dim2+child2_mb*dim3+(l2+1)*dim4];
                            for(auto child2 : topo2[curr2]){
                                if(child2==child2_mb){
                                    continue;
                                }
                                d_ += memT[nn1+0*dim2+child2*dim3+1*dim4];
                            }
                            d = std::min(d,d_);
                        }
                        memT[curr1+l1*dim2+curr2*dim3+l2*dim4] = d;
                    }
                    //---------------------------------------------------------------------------
                    // If second tree only has one branch, try all decompositions of first tree
                    else if(topo2[curr2].size()==0){
                        float d = FLT_MAX;
                        for(auto child1_mb : topo1[curr1]){
                            float d_ = memT[child1_mb+(l1+1)*dim2+curr2*dim3+l2*dim4];
                            for(auto child1 : topo1[curr1]){
                                if(child1==child1_mb){
                                    continue;
                                }
                                d_ += memT[child1+1*dim2+nn2*dim3+0*dim4];
                            }
                            d = std::min(d,d_);
                        }
                        memT[curr1+l1*dim2+curr2*dim3+l2*dim4] = d;
                    }
                    //---------------------------------------------------------------------------
                    // If both trees have more than one branch, try all decompositions of both trees
                    else{
                        float d = FLT_MAX;
                        //-----------------------------------------------------------------------
                        // Try all possible main branches of first tree (child1_mb) and all possible main branches of second tree (child2_mb)
                        // Then try all possible matchings of subtrees
                        if(topo1[curr1].size()==2 && topo2[curr2].size()==2){
                            int child11 = topo1[curr1][0];
                            int child12 = topo1[curr1][1];
                            int child21 = topo2[curr2][0];
                            int child22 = topo2[curr2][1];
                            d = std::min(d,memT[child11+(l1+1)*dim2+child21*dim3+(l2+1)*dim4] + memT[child12+1*dim2+child22*dim3+1*dim4]);
                            d = std::min(d,memT[child12+(l1+1)*dim2+child22*dim3+(l2+1)*dim4] + memT[child11+1*dim2+child21*dim3+1*dim4]);
                            d = std::min(d,memT[child11+(l1+1)*dim2+child22*dim3+(l2+1)*dim4] + memT[child12+1*dim2+child21*dim3+1*dim4]);
                            d = std::min(d,memT[child12+(l1+1)*dim2+child21*dim3+(l2+1)*dim4] + memT[child11+1*dim2+child22*dim3+1*dim4]);
                        }
                        else{
                            for(auto child1_mb : topo1[curr1]){
                                auto topo1_ = topo1[curr1];
                                topo1_.erase(std::remove(topo1_.begin(),topo1_.end(),child1_mb),topo1_.end());
                                for(auto child2_mb : topo2[curr2]){
                                    float d_ = memT[child1_mb+(l1+1)*dim2+child2_mb*dim3+(l2+1)*dim4];
                                    auto topo2_ = topo2[curr2];
                                    topo2_.erase(std::remove(topo2_.begin(),topo2_.end(),child2_mb),topo2_.end());
                                    auto f = [&] (unsigned r, unsigned c) {
                                        int c1 = r<topo1_.size() ? topo1_[r] : -1;
                                        int c2 = c<topo2_.size() ? topo2_[c] : -1;
                                        return memT[c1+1*dim2+c2*dim3+1*dim4];
                                    };
                                    int size = std::max(topo1_.size(),topo2_.size()) + 1;
                                    auto matching = munkres_algorithm<float>(size, size, f);
                                    for(auto m : matching) d_ += f(m.first,m.second);
                                    d = std::min(d,d_);
                                }
                            }
                        }
                        //-----------------------------------------------------------------------
                        // Try to continue main branch on one child of first tree and delete all other subtrees
                        // Then match continued branch to current branch in second tree
                        for(auto child1_mb : topo1[curr1]){
                            float d_ = memT[child1_mb+(l1+1)*dim2+curr2*dim3+l2*dim4];
                            for(auto child1 : topo1[curr1]){
                                if(child1 == child1_mb){
                                    continue;
                                }
                                d_ += memT[child1+1*dim2+nn2*dim3+0*dim4];
                            }
                            d = std::min(d,d_);
                        }
                        //-----------------------------------------------------------------------
                        // Try to continue main branch on one child of second tree and delete all other subtrees
                        // Then match continued branch to current branch in first tree
                        for(auto child2_mb : topo2[curr2]){
                            float d_ = memT[curr1+l1*dim2+child2_mb*dim3+(l2+1)*dim4];
                            for(auto child2 : topo2[curr2]){
                                if(child2 == child2_mb){
                                    continue;
                                }
                                d_ += memT[nn1+0*dim2+child2*dim3+1*dim4];
                            }
                            d = std::min(d,d_);
                        }
                        memT[curr1+l1*dim2+curr2*dim3+l2*dim4] = d;
                    }

                }
            }
        }
    }

    float res = memT[topo1[rootID1][0]+1*dim2+topo2[rootID2][0]*dim3+1*dim4];
    //delete[] memT;
    return res;
}