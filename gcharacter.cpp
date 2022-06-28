#include "mgraph.h"
#include "atomic_data.h"
#include "load_data.h"
#include<vector>
#include<iostream>
#include<string>
#include<map>
#include<algorithm>
#include<set>
#include<Eigen/Dense>

Eigen::MatrixXd edm(Eigen::MatrixXd &xyzMat)
{
    
    Eigen::MatrixXd mat2, mat2T; // one is a column vector, another is a row vector
    mat2 = xyzMat.array().square().rowwise().sum(); 
    mat2T = mat2.transpose();

    int size;
    size = xyzMat.rows(); // no. of rows and columns are identical

    Eigen::MatrixXd sum2;
    sum2.resize(size, size);
    sum2 = mat2.replicate(1, size) + mat2T.replicate(size,1); // matrix containing the sum of squares

    Eigen::MatrixXd abdot;
    abdot.resize(size, size);
    abdot = -2 * xyzMat * xyzMat.transpose();

    Eigen::MatrixXd edm;
    edm.resize(size, size);
    edm = abdot + sum2;
    edm = edm.array().sqrt();
    edm.diagonal().setZero();
    return edm;
}

double get_bo(std::vector<std::string> &atompairs, double &dist)
{
    double tol = 0.003;
    double bo;
    std::string atom1, atom2;
    atom1 = atompairs[0];
    atom2 = atompairs[1];

    std::vector<std::string> atompairsr;
    atompairsr = {atom2, atom1};

    std::map<double, std::vector<double>> pairBO;
    std::vector<double> bondorders;
    if (bondr.count(atompairs) || bondr.count(atompairsr))
    {
        
        if (bondr.count(atompairs))
        {
            pairBO = bondr.find(atompairs)->second;
        }
        else if (bondr.count(atompairsr))
        {
            pairBO = bondr.find(atompairsr)->second;
        }

        for (auto iter = pairBO.begin(); iter != pairBO.end(); iter++)
        {
            std::vector<double> boundary;
            boundary = iter->second;
            double lower, upper;
            lower = boundary[0] - tol;
            upper = boundary[1] + tol;
            // std::cout << "lower and upper: " << lower << " " << upper << std::endl;
            
            if (dist >= lower && dist <= upper)
            {
                // std::cout << "bo: " << iter->first << std::endl;
                bondorders.push_back(iter->first);
            }
            
        }
    }
    else
    {
        // std::cout << atom1 << " " << atom2 << std::endl;
        std::cout << "this program only handles the elements H, B, C, N, O, F, P, S, Cl, Br and I" << std::endl;
        exit(1);
    }

    if (bondorders.empty())
    {
        bo = 0.0;
    }
    else
    {
        std::vector<double>::iterator result;
        result = std::max_element(bondorders.begin(), bondorders.end());
        bo = *result;
    }

    return bo;
    
}


void graph_edges(mgraph &graph, Eigen::MatrixXd &edm)
{
    int size;
    size = edm.rows();
    std::vector<std::array<int,2>> edge_vector;
    for (int i = 0; i < size; i++) // iterating over the nodes, i = node label
    {   
        if (graph.nodes(i).get_element() != "H") // non hydrogen atoms
        {
            int bondED = 0;
            int bondElec = 0;
            // std::cout << "ha node: " << i  << " " << graph.nodes(i).get_element() << std::endl;
            Eigen::RowVectorXd bondVec;
            // bondVec.resize(1, size); // row vector, double check this
            
            bondVec = edm.row(i);

            std::vector<int> idxVec; // contains the indices where the bonding is less than 3 angstroms
            for (int j = 0; j < size; j++)
            {
                if (bondVec[j] < 3.0 && bondVec[j] > 0.0)
                {
                    // std::cout << "bondVec[j]: " << bondVec[j] << std::endl;
                    idxVec.push_back(j);
                }
            }
            
            for (int k = 0; k < idxVec.size(); k++)
            {
                // std::cout << "\t" << idxVec[k] << std::endl;
                std::string atom1, atom2;
                atom1 = graph.nodes(i).get_element();
                atom2 = graph.nodes(idxVec[k]).get_element(); // second node label is idxVec[k]
                // std::cout << atom1 << " " << atom2 << std::endl;
                std::vector<std::string> atompairs;
                atompairs = {atom1, atom2};
                double dist = bondVec[idxVec[k]];
                // std::cout << "\t" << "dist: " << dist << std::endl;
                double bo = get_bo(atompairs, dist);

                // std::cout << "bond order: " << bo << std::endl;
                if (bo != 0.0)
                {
                    bondElec += bo * 2;
                    bondED += 1;
                    std::array<int,2> elabel;
                    if (i < idxVec[k])
                    {
                        elabel[0] = i;
                        elabel[1] = idxVec[k];
                        // elabel = {i, idxVec[k]};
                        // elabel.fill({ {i, idxVec[k]} });
                    }
                    else
                    {
                        elabel[0] = idxVec[k];
                        elabel[1] = i;
                        // elabel = {idxVec[k], i}
                        // elabel.fill({ {idxVec[k], i} });
                    }

                    if (std::find(edge_vector.begin(), edge_vector.end(),elabel)!=edge_vector.end())
                    {
                        continue;
                    }
                    else
                    {
                        edge egraph;
                        egraph.set_bo(bo);
                        egraph.set_bl(dist);
                        graph.add_edge(elabel, egraph);
                        edge_vector.push_back(elabel);
                    }
                    
                }
                
            }
            int valElec;
            int elecDom;
            valElec = get_valence(graph.nodes(i).get_element());
            elecDom = (int)ceil(bondED + 0.5 * (valElec - 0.5 * bondElec ));
            // std::cout << "ed: " << elecDom << std::endl;

            std::string element;
            std::string atomtype;
            element = graph.nodes(i).get_element();

            // getting preliminary atom types
            // atom types are updated later to account for atoms in groups participating in resonance
            if (element == "F" || element == "I")
            {
                atomtype = element + "_";
            }
            else if (element == "Cl" || element == "Br")
            {
                atomtype = element;
            }
            else if (element == "O" || element == "N" || element == "C")
            {
                atomtype = element + "_" + std::to_string(elecDom - 1);
            }
            else if (element == "S")
            {
                if (elecDom == 3)
                {
                    atomtype = "S_2";
                }
                else if (elecDom == 4)
                {
                    std::vector<int> neighNodes;
                    std::vector<std::string> neighElem;
                    neighNodes = graph.neighbors(i);
                    for (int v = 0; v < neighNodes.size(); v++)
                    {
                        std::string neigh;
                        neigh = graph.nodes(neighNodes[v]).get_element();
                        neighElem.push_back(neigh);
                    }
                    int ocount;
                    ocount = std::count(neighElem.begin(), neighElem.end(), "O");
                    if (ocount == 2)
                    {
                        atomtype = "S_3+4";
                    }
                    
                    else if (ocount == 3)
                    {
                        atomtype = "S_3+6";
                    }

                    else 
                    {
                        atomtype = "S_3+2";
                    }
                }
            }
            node ngraph;
            ngraph = graph.nodes(i);
            ngraph.set_coordno(elecDom);
            ngraph.set_at(atomtype);
            graph.update_node(i, ngraph);
            }

            
        
        else
        {
            std::string atomtype;
            node ngraph;
            ngraph = graph.nodes(i);
            atomtype = "H_";
            ngraph.set_at(atomtype);
            graph.update_node(i, ngraph);
        }
        }
}

void check_hybrid(mgraph &graph) // check for hybridisation of oxygen, sulfur and nitrogen (cases like furan, pyrrole and thiophene where the coordination number of the heteroatom is initially determined to be 4)
{
    std::vector<int> nos_nodes;
    std::vector<int> node_vec;
    node_vec = graph.get_nodes();

    for (int i = 0; i < node_vec.size(); i++)
    {
        node ngraph;
        ngraph = graph.nodes(i);
        std::string element;
        element = ngraph.get_element();
        int coordno;
        coordno = ngraph.get_coordno();

        if ((element == "N" && coordno == 4) || element == "O" && coordno == 4|| element == "S" && coordno == 4)
        {
            // std::cout << "node i: " << i << std::endl;
            std::vector<int> neigh_vec;
            neigh_vec = graph.neighbors(i);
            std::vector<int> neigh_coordno;
            for (int j = 0; j < neigh_vec.size(); j++)
            {
                int nneigh_label;
                nneigh_label = neigh_vec[j];
                int coordno;
                node nneigh;
                nneigh = graph.nodes(nneigh_label);
                coordno = nneigh.get_coordno();
                neigh_coordno.push_back(coordno);
            }
            int sp2neigh;
            sp2neigh = std::count(neigh_coordno.begin(), neigh_coordno.end(), 3);

            if (sp2neigh > 1)
            {
                int elecDom = 3;
                std::string atomtype;
                atomtype = element + "_2";
                ngraph.set_coordno(elecDom);
                ngraph.set_at(atomtype);
                graph.update_node(i, ngraph);
            }
        }

        if (element == "N" && coordno == 4) // cases like formamide
        {
            std::vector<int> neigh_vec;
            neigh_vec = graph.neighbors(i);
            int status = 0;
            for (int j = 0; j < neigh_vec.size(); j++)
            {
                int neigh_node;
                neigh_node = neigh_vec[j];
                std::vector<int> nneigh_vec;
                nneigh_vec = graph.neighbors(neigh_node);
                nneigh_vec.erase(std::remove(nneigh_vec.begin(), nneigh_vec.end(), i), nneigh_vec.end());
                for (int k = 0; k < nneigh_vec.size(); k++)
                {
                    int pair1, pair2;
                    pair1 = neigh_node;
                    pair2 = nneigh_vec[k];

                    int coord1, coord2;
                    coord1 = graph.nodes(pair1).get_coordno();
                    coord2 = graph.nodes(pair2).get_coordno();

                    if (coord1 == 3 && coord2 == 3)
                    {
                        status = 1;
                        break;
                    }
                }

                if (status == 1)
                {
                    int elecDom = 3;
                    std::string atomtype;
                    atomtype = "N_2";
                    ngraph.set_coordno(elecDom);
                    ngraph.set_at(atomtype);
                    graph.update_node(i, ngraph);
                    break;
                }
            }
        }
    }
}

std::vector<std::vector<std::array<int,2>>> conjugate_region(mgraph &graph)
{
    std::vector<std::vector<std::array<int,2>>> conjugated_edges;
    std::vector<std::array<int,2>> edge_vec;
    edge_vec = graph.get_edges();
    std::vector<std::array<int,2>> unsat_edges;
    std::set<int> unsat_nodes;

    // std::cout << "unsaturated edges" << std::endl;
    for (int i = 0; i < edge_vec.size(); i++)
    {
        std::array<int,2> egraph;
        egraph = edge_vec[i];
        node node1, node2;
        int node1_label, node2_label;
        node1_label = egraph[0];
        node2_label = egraph[1];
        node1 = graph.nodes(node1_label);
        node2 = graph.nodes(node2_label);
        int cn1, cn2; // coordination number of node1 and node2

        cn1 = node1.get_coordno();
        cn2 = node2.get_coordno();

        if (cn1 <= 3 && cn1 > 1 && cn2 <= 3 && cn2 > 1)
        {
            unsat_edges.push_back(egraph);
            unsat_nodes.insert(node1_label);
            unsat_nodes.insert(node2_label);
            // std::cout << node1_label << " " << node2_label << std::endl;
        }
    }   

    mgraph cgraph;
    for (auto it = unsat_nodes.begin(); it != unsat_nodes.end(); it++)
    {
        node cnode;
        int node_label = *it;
        cgraph.add_node(node_label, cnode);
    }

    for (int j =0; j < unsat_edges.size(); j++)
    {
        edge cedge;
        cgraph.add_edge(unsat_edges[j], cedge);
    }

    std::vector<std::vector<int>> connected_comp;
    connected_comp = cgraph.connected_components();

    for (int k = 0; k < connected_comp.size(); k++)
    {
        std::vector<int> comp_nodes;
        std::vector<int> cnodes;
        comp_nodes = connected_comp[k];

        if (comp_nodes.size() > 2)
        {
            std::vector<std::array<int,2>> comp_edges;
            for (int d = 0; d < unsat_edges.size(); d++)
            {
                std::array<int,2> ue;
                ue = unsat_edges[d];
                int uen1, uen2; // nodes in the unsaturated edge
                uen1 = ue[0];
                uen2 = ue[1];
                if ((std::find(comp_nodes.begin(), comp_nodes.end(), uen1) != comp_nodes.end()) && (std::find(comp_nodes.begin(), comp_nodes.end(), uen2) != comp_nodes.end()))
                {
                    comp_edges.push_back(ue);
                    edge egraph;
                    egraph = graph.edges(ue);
                    bool conjstat = true;
                    egraph.set_conjstatus(conjstat);
                    graph.update_edge(ue, egraph);
                }
            }
            conjugated_edges.push_back(comp_edges);
        }
        /*for (int v = 0; v < comp_nodes.size(); v++)
        {
            std::string element;
            element = graph.nodes(comp_nodes[v]).get_element();
            if (element == "C")
            {
                cnodes.push_back(comp_nodes[v]);
            }
        }

        std::vector<int> prob_nodes; // problematic nodes e.g. case of allenes

        for (int cn = 0; cn < cnodes.size(); cn++)
        {
            int carbon_node;
            carbon_node = cnodes[cn];
            std::vector<int> cneigh;
            cneigh = cgraph.neighbors(carbon_node);
            std::vector<double> bo_vec;

            for (int n = 0; n < cneigh.size(); n++)
            {
                double bond_order;
                int neighbor;
                neighbor = cneigh[n];
                std::array<int,2> elabel;
                
                if (neighbor < carbon_node)
                {
                    elabel[0] = neighbor;
                    elabel[1] = carbon_node;
                }

                else
                {
                    elabel[0] = carbon_node;
                    elabel[1] = neighbor;
                }

                bond_order = graph.edges(elabel).get_bo();
                bo_vec.push_back(bond_order);

                
            }

            dtcount = std::count(bo_vec.begin(), bo_vec.end(), 2.00);

            if (dtcount >= 2)
            {
                prob_nodes.push_back(carbon_node);
            }


        }*/
    }

    
    return conjugated_edges;
}