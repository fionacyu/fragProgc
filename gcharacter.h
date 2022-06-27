#ifndef GCHARACTER_H
#define GCHARACTER_H
#include "mgraph.h"
#include<vector>
#include<iostream>
#include<string>
#include<map>
#include<algorithm>
#include<Eigen/Dense>

Eigen::MatrixXd edm(Eigen::MatrixXd &xyzMat);
double get_bo(std::vector<std::string> &atompairs, double &dist);
void graph_edges(mgraph &graph, Eigen::MatrixXd &edm);
void check_hybrid(mgraph &graph);

#endif