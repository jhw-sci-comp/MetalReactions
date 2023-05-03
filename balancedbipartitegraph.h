#ifndef BALANCEDBIPARTITEGRAPH_H
#define BALANCEDBIPARTITEGRAPH_H


#include <iostream>


#include <utility>
#include <vector>
#include <algorithm>


/* Balanced bipartite Graph providing methods to find a perfect matching */
class BalancedBipartiteGraph {
public:
    BalancedBipartiteGraph(int);
    BalancedBipartiteGraph(int, std::vector<std::pair<int, int>>);
    bool findMatch(int, std::vector<bool>&, std::vector<int>&, std::vector<int>&);
    std::vector<int> findMaxMatch();

private:
    int n;
    std::vector<std::pair<int, int>> edges;

};

#endif