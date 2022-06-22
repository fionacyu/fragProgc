#include "input.h"
#include "load_data.h"
#include "mgraph.h"
#include "gcharacter.h"
#include<iostream>
#include<string>
#include<stdio.h>
#include<Eigen/Dense>

int main(int argc, char** argv){
    // std::cout << "cmdline args count= " << argc << std::endl;

    input_f input;
    
    if (argc != 2)
    {
        std::cout << "Usage: " << "./main filename.input" << std::endl; 
        exit(1);
    }

    process_files(input, argc, argv);

    mgraph graph;
    Eigen::MatrixXd xyzMat;
    Eigen::MatrixXd edMat;
    extract_coordinates(input.get_coordfile(), graph, xyzMat);

    edMat = edm(xyzMat);
    // std::cout << edMat << std::endl;
    graph_edges(graph, edMat);
    
    /*for (int i = 0; i < graph.get_nodes().size(); i++)
    {
        int nodelabel = graph.get_nodes()[i];
        std::cout << "node label: " << nodelabel << std::endl;
        std::cout << "\t" << "element: " << graph.nodes(nodelabel).get_element() << std::endl;
        std::cout << "\t" << "coord: " << graph.nodes(nodelabel).get_coordinates()[0] << " " << graph.nodes(nodelabel).get_coordinates()[1] << " " << graph.nodes(nodelabel).get_coordinates()[2] <<std::endl;
        std::cout << "\t" << "ed: " << graph.nodes(nodelabel).get_coordno() << std::endl;
        std::cout << "\t" << "at: " << graph.nodes(nodelabel).get_at() << std::endl;
    }*/


    /*for (int i = 0; i < graph.get_edges().size(); i++)
    {
        std::array<int,2> elabel;
        elabel = graph.get_edges()[i];
        std::cout << "edge: " << elabel[0] << ", " << elabel[1] << std::endl;
        double bo, bl;
        edge egraph;
        egraph = graph.edges(elabel);
        std::cout << "\t" << "BO: " << egraph.get_bo() << std::endl;
        std::cout << "\t" << "BL: " << egraph.get_bl() << std::endl;
    }*/

    return 0;
}