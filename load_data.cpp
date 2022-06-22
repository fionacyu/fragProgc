#include<string>
#include<regex>
#include<sstream>
#include<iostream>
#include<Eigen/Dense>
#include "load_data.h"
#include "input.h"
#include "mgraph.h"
#include "atomic_data.h"

void process_files(input_f &input_object, int argc, char** argv)
{
    input_object = process_input(argv[1]);
    // std::cout << "fragSize " << input_object.get_fragSize() << std::endl;
    // std::cout << "coordPath " << input_object.get_coordfile() << std::endl;
}

void extract_coordinates(std::string coordPath, mgraph &graph, Eigen::MatrixXd &xyzMat)
{
    std::string coord_contents;
    std::string coord_string; 
    std::regex coord_pattern("-?[[:digit:]]{1,2}.[[:digit:]]");
    std::regex element_pattern("[[:alpha:]]{1,2}");
    std::regex int_pattern("^0$|^[1-9][0-9]*$");
    coord_contents = readFile(coordPath);

    int nodecount = 0;
    int coord_count = 0; 
    std::stringstream coordstream(coord_contents);
    std::vector<double>xyz_vec; 


    while (coordstream >> coord_string)
    {
        int atomno;
        if (regex_search(coord_string, int_pattern))
        {
            atomno = stoi(coord_string);
        }
        xyzMat.resize(atomno, 3);
        if (regex_search(coord_string, element_pattern)) // element will be encountered first
        {
            std::string element;
            element += coord_string[coord_string.length() - 1];
            coord_string = regex_replace(coord_string, element_pattern, "");

            node graphn;
            graphn.set_element(element);
            graph.add_node(nodecount, graphn);
            nodecount += 1; 
        }


        if (regex_search(coord_string, coord_pattern)) // xyz coordinates
        {
            coord_count += 1;
            if (coord_count < 4)
            {
                xyz_vec.push_back(stod(coord_string));
                if (xyz_vec.size() == 3)
                {
                    int correspond_node = nodecount - 1;
                    node cnode;
                    cnode = graph.nodes(correspond_node);
                    cnode.set_coordinates(xyz_vec);
                    graph.update_node(correspond_node, cnode);

                    Eigen::Vector3d v1;
                    v1<< xyz_vec[0], xyz_vec[1], xyz_vec[2];
                    // std::cout << v1 << std::endl;
                    xyzMat.row(correspond_node) = v1;
                    xyz_vec.clear(); 
                    coord_count = 0;
                    
                }
            }
        }
    }

    
    // std::cout << xyzMat << std::endl;
}

int get_valence(std::string elem)
{
    int val;
    if (valence.count(elem))
    {
        val = valence.find(elem)->second;
    }

    else
    {
        std::cout << "valence for " << elem << " not found." << std::endl;
        exit(1);
    }
    return val;
}