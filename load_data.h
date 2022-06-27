#ifndef LOADDATA_H
#define LOADDATA_H
#include<Eigen/Dense>
#include "input.h"
#include "mgraph.h"
#include "atomic_data.h"

void process_files(input_f &input_object, int argc, char** argv);
void extract_coordinates(std::string coordPath, mgraph &graph, Eigen::MatrixXd &xyzMat);
int get_valence(std::string elem);
std::map<std::string, aType> parameter_definition();
// void parameter_definition();
#endif