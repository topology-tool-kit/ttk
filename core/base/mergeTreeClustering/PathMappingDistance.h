/// \ingroup base
/// \class PathMappingDistance
/// \author Florian Wetzels (wetzels@cs.uni-kl.de)
/// \date 2022.
///
/// This module defines the %PathMappingDistance class that computes distances
/// between two merge trees.
///
/// \b Related \b publication \n
/// "A Deformation-based Edit Distance for Merge Trees" \n
/// Florian Wetzels, Christoph Garth.\n
/// TBP.

#pragma once

#include <vector>
#include <set>

#include <iostream>
#include <vector>
#include <stack>
#include <tuple>
#include <limits>
#include <cmath>
#include <cfloat>
#include <set>
#include <algorithm>
#include <chrono>

// ttk common includes
#include <Debug.h>
#include <AssignmentMunkres.h>
#include <AssignmentAuction.h>
#include <AssignmentExhaustive.h>

namespace ttk {

    class PathMappingDistance : virtual public Debug {

        private:
        int baseMetric = 0;
        int assignmentSolverID_ = 0;
        bool squared = false;

        template <class dataType>
        inline dataType editCost_Persistence(int n1, int p1, int n2, int p2, std::vector<dataType> &nodes1, std::vector<dataType> &nodes2){
            if(n1<0){
                dataType b1 = nodes2[n2];
                dataType d1 = nodes2[p2];
                return squared ? std::abs(d1-b1)*std::abs(d1-b1) : std::abs(d1-b1);
            }
            if(n2<0){
                dataType b1 = nodes1[n1];
                dataType d1 = nodes1[p1];
                return squared ? std::abs(d1-b1)*std::abs(d1-b1) : std::abs(d1-b1);
            }
            dataType b1 = nodes1[n1];
            dataType d1 = nodes1[p1];
            dataType b2 = nodes2[n2];
            dataType d2 = nodes2[p2];
            dataType d = std::abs(std::abs(b1-d1)-std::abs(b2-d2));
            return squared ? d*d : d;
        }

        public:
        PathMappingDistance() {
            this->setDebugMsgPrefix(
                "MergeTreeDistance"); // inherited from Debug: prefix will be printed at
                                    // the beginning of every msg
        }
        ~PathMappingDistance() override = default;

        void setBaseMetric(int m) {
            baseMetric = m;
        }

        void setAssignmentSolver(int assignmentSolver) {
            assignmentSolverID_ = assignmentSolver;
        }

        void setSquared(bool s){
            squared = s;
        }

        template <class dataType>
        dataType editDistance_path(std::vector<dataType> &nodes1,
                                std::vector<std::vector<int>> &topo1,
                                int rootID1,
                                std::vector<dataType> &nodes2,
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

            std::unique_ptr<dataType[]> memT(new dataType[(nn1+1)*(depth1+1)*(nn2+1)*(depth2+1)]);
            //std::vector<dataType> memT((nn1+1)*(depth1+1)*(nn2+1)*(depth2+1));

            memT[nn1+0*dim2+nn2*dim3+0*dim4] = 0;
            for(size_t i=0; i<nn1; i++){
                for(size_t l=1; l<=predecessors1[i].size(); l++){

                    int curr1 = i;
                    int parent1 = predecessors1[i][predecessors1[i].size()-l];

                    //-----------------------------------------------------------------------
                    // Delete curr path and full subtree rooted in path
                    memT[curr1+l*dim2+nn2*dim3+0*dim4] = editCost_Persistence(curr1,parent1,-1,-1,nodes1,nodes2);
                    for(auto child1 : topo1[curr1]){
                        memT[curr1+l*dim2+nn2*dim3+0*dim4] += memT[child1+1*dim2+nn2*dim3+0*dim4];
                    }
                    
                }
            }
            for(size_t j=0; j<nn2; j++){
                for(size_t l=1; l<=predecessors2[j].size(); l++){

                    int curr2 = j;
                    int parent2 = predecessors2[j][predecessors2[j].size()-l];

                    //-----------------------------------------------------------------------
                    // Delete curr path and full subtree rooted in path
                    memT[nn1+0*dim2+curr2*dim3+l*dim4] = editCost_Persistence(-1,-1,curr2,parent2,nodes1,nodes2);
                    for(auto child2 : topo2[curr2]){
                        memT[nn1+0*dim2+curr2*dim3+l*dim4] += memT[nn1+0*dim2+child2*dim3+1*dim4];
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
                                memT[curr1+l1*dim2+curr2*dim3+l2*dim4] = editCost_Persistence(curr1,parent1,curr2,parent2,nodes1,nodes2);
                            }
                            //---------------------------------------------------------------------------
                            // If first tree only has one branch, try all decompositions of second tree
                            else if(topo1[curr1].size()==0){
                                dataType d = FLT_MAX;
                                for(auto child2_mb : topo2[curr2]){
                                    dataType d_ = memT[curr1+l1*dim2+child2_mb*dim3+(l2+1)*dim4];
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
                                dataType d = FLT_MAX;
                                for(auto child1_mb : topo1[curr1]){
                                    dataType d_ = memT[child1_mb+(l1+1)*dim2+curr2*dim3+l2*dim4];
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
                                dataType d = FLT_MAX;
                                //-----------------------------------------------------------------------
                                // Try all possible main branches of first tree (child1_mb) and all possible main branches of second tree (child2_mb)
                                // Then try all possible matchings of subtrees
                                if(topo1[curr1].size()==2 && topo2[curr2].size()==2){
                                    int child11 = topo1[curr1][0];
                                    int child12 = topo1[curr1][1];
                                    int child21 = topo2[curr2][0];
                                    int child22 = topo2[curr2][1];
                                    d = std::min(d,memT[child11+1*dim2+child21*dim3+1*dim4] + memT[child12+1*dim2+child22*dim3+1*dim4] + editCost_Persistence(curr1,parent1,curr2,parent2,nodes1,nodes2));
                                    d = std::min(d,memT[child11+1*dim2+child22*dim3+1*dim4] + memT[child12+1*dim2+child21*dim3+1*dim4] + editCost_Persistence(curr1,parent1,curr2,parent2,nodes1,nodes2));
                                }
                                else{
                                    auto f = [&] (unsigned r, unsigned c) {
                                        size_t c1 = r<topo1[curr1].size() ? topo1[curr1][r] : nn1;
                                        size_t c2 = c<topo2[curr2].size() ? topo2[curr2][c] : nn2;
                                        int l1_ = c1==nn1?0:1;
                                        int l2_ = c2==nn2?0:1;
                                        return memT[c1+l1_*dim2+c2*dim3+l2_*dim4];
                                    };
                                    int size = std::max(topo1[curr1].size(),topo2[curr2].size()) + 1;
                                    auto costMatrix = std::vector<std::vector<dataType>>(size,std::vector<dataType>(size,0));
                                    std::vector<MatchingType> matching;
                                    for(int r=0; r<size; r++){
                                        for(int c=0; c<size; c++){
                                            costMatrix[r][c] = f(r,c);
                                        }
                                    }

                                    AssignmentSolver<dataType> *assignmentSolver;
                                    AssignmentExhaustive<dataType> solverExhaustive;
                                    AssignmentMunkres<dataType> solverMunkres;
                                    AssignmentAuction<dataType> solverAuction;
                                    switch(assignmentSolverID_) {
                                        case 1:
                                        solverExhaustive = AssignmentExhaustive<dataType>();
                                        assignmentSolver = &solverExhaustive;
                                        break;
                                        case 2:
                                        solverMunkres = AssignmentMunkres<dataType>();
                                        assignmentSolver = &solverMunkres;
                                        break;
                                        case 0:
                                        default:
                                        solverAuction = AssignmentAuction<dataType>();
                                        //solverAuction.setEpsilon(auctionEpsilon_);
                                        //solverAuction.setEpsilonDiviserMultiplier(auctionEpsilonDiviser_);
                                        //solverAuction.setNumberOfRounds(auctionRound_);
                                        assignmentSolver = &solverAuction;
                                    }
                                    assignmentSolver->setInput(costMatrix);
                                    assignmentSolver->setBalanced(true);
                                    assignmentSolver->run(matching);
                                    dataType d_ = editCost_Persistence(curr1,parent1,curr2,parent2,nodes1,nodes2);
                                    for(auto m : matching) d_ += std::get<2>(m);
                                    d = std::min(d,d_);
                                }
                                //-----------------------------------------------------------------------
                                // Try to continue main branch on one child of first tree and delete all other subtrees
                                // Then match continued branch to current branch in second tree
                                for(auto child1_mb : topo1[curr1]){
                                    dataType d_ = memT[child1_mb+(l1+1)*dim2+curr2*dim3+l2*dim4];
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
                                    dataType d_ = memT[curr1+l1*dim2+child2_mb*dim3+(l2+1)*dim4];
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

            dataType res = memT[topo1[rootID1][0]+1*dim2+topo2[rootID2][0]*dim3+1*dim4];
            //delete[] memT;
            return squared ? std::sqrt(res) : res;
        }
    };
}
