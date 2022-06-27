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

    std::map<std::string, aType> prm_map;
    prm_map = parameter_definition();

    mgraph graph;
    Eigen::MatrixXd xyzMat;
    Eigen::MatrixXd edMat;
    extract_coordinates(input.get_coordfile(), graph, xyzMat);

    edMat = edm(xyzMat);
    // std::cout << edMat << std::endl;
    graph_edges(graph, edMat);
    check_hybrid(graph);
    
    std::vector<std::vector<std::array<int,2>>> conjugated_edges;
    conjugated_edges = conjugate_region(graph);

    std::vector<std::vector<int>> cc;
    cc = graph.connected_components();
    for (int i = 0; i < cc.size(); i++)
    {
        std::vector<int> _cc;
        _cc = cc[i];
        std::cout << "component" << i << std::endl;
        for (int j = 0; j < _cc.size(); j++)
        {
            std::cout << "\t" << _cc[j] << std::endl;
        }
    }

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

    /*for (auto iter = prm_map.begin(); iter != prm_map.end(); iter++)
    {
        aType atomtype;
        atomtype = iter->second;
        std::cout << iter->first << " " << std::endl;
        std::cout << "\t" << "r1: " << atomtype.get_r1() << std::endl;
        std::cout << "\t" << "theta0: " << atomtype.get_theta0() << std::endl;
        std::cout << "\t" << "x1: " << atomtype.get_x1() << std::endl;
        std::cout << "\t" << "D1: " << atomtype.get_D1() << std::endl;
        std::cout << "\t" << "zeta: " << atomtype.get_zeta() << std::endl;
        std::cout << "\t" << "Z1: " << atomtype.get_Z1() << std::endl;
        std::cout << "\t" << "Vi: " << atomtype.get_Vi() << std::endl;
        std::cout << "\t" << "Uj: " << atomtype.get_Uj() << std::endl;
        std::cout << "\t" << "Xi: " << atomtype.get_Xi() << std::endl;
    }*/
    
    return 0;
}