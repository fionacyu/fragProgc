#include "mgraph.h"
#include<string>
#include<vector>
#include<array>
#include<algorithm>

// node functions
void node::set_element(std::string &ele)
{
    element = ele;
}

void node::set_coordno(int &cn)
{
    coordno = cn;
}

void node::set_coordinates(std::vector<double> &coordinates)
{
    coord = coordinates;
}

void node::set_box(std::vector<double> &bcoord)
{
    box = bcoord;
}

void node::set_at(std::string &atomtype)
{
    at = atomtype;
}

void node::set_pi(int &piElec)
{
    pi = piElec;
}

std::string node::get_element()
{
    return element;
}

int node::get_coordno()
{
    return coordno;
}

std::vector<double> node::get_coordinates()
{
    return coord;
}

std::vector<double> node::get_box()
{
    return box;
}

std::string node::get_at()
{
    return at;
}

int node::get_pi()
{
    return pi;
}

// edge functions
void edge::set_bo(double &bondorder)
{
    bo = bondorder;
}

void edge::set_conjstatus(bool &status)
{
    conjugated = status;
}

void edge::set_bl(double &length)
{
    bl = length;
}

double edge::get_bo()
{
    return bo;
}

bool edge::get_conjstatus()
{
    return conjugated;
}

double edge::get_bl()
{
    return bl;
}

// mgraph functions
void mgraph::add_node(int &label, node &n)
{
    mnodes.insert(std::pair<int, node>(label, n));
    nodeVec.push_back(label);
}

void mgraph::add_edge(std::array<int,2> &elabel, edge &e)
{
    medges.insert(std::pair<std::array<int,2>, edge>(elabel, e));
    edgeVec.push_back(elabel);
}

void mgraph::update_node(int &nodelabel, node &n)
{
    mnodes.find(nodelabel)->second = n;
}

void mgraph::update_edge(std::array<int,2> & edgelabel, edge &e)
{
    medges.find(edgelabel)->second = e;
}

std::vector<int> mgraph::get_nodes()
{
    return nodeVec;
}

std::vector<std::array<int,2>> mgraph::get_edges()
{
    return edgeVec;
}

node mgraph::nodes(int &nodelabel)
{
    return mnodes.find(nodelabel)->second; // returns the key value which is of data type node
}

edge mgraph::edges(std::array<int,2> &edgelabel)
{
    return medges.find(edgelabel)->second;
}

std::vector<int> mgraph::neighbors(int &nodelabel)
{
    std::vector<int> neighVec;
    std::vector<std::array<int,2>> edgeVec;
    edgeVec = get_edges();

    for (int i =0; i < edgeVec.size(); i++)
    {
        std::array<int,2> elabel;
        elabel = edgeVec[i];
        if (std::find(elabel.begin(), elabel.end(),nodelabel)!=elabel.end())
        {
            if (elabel[0] == nodelabel)
            {
                neighVec.push_back(elabel[1]);
            }
            else if (elabel[1] == nodelabel)
            {
                neighVec.push_back(elabel[0]);
            }
        }
    }
    return neighVec;
}

std::vector<int> mgraph::_plain_bfs(int &node_label)
{
    // std::map<int, std::vector<int>> graphmap;
    // std::vector<int> nodes_vec;
    // nodes_vec = get_nodes();
    // for (int i = 0; i < nodes_vec.size(); i++)
    // {
    //     std::vector<int> nneigh;
    //     nneigh = neighbors(node_label);
    //     graphmap.insert(std::pair<int, std::vector<int>>(node_label, nneigh));
    // }

    std::vector<int> seen;
    std::vector<int> nextlevel;
    nextlevel.push_back(node_label);
    while (nextlevel.size() > 0)
    {
        std::vector<int> thislevel;
        thislevel = nextlevel;
        nextlevel.clear();
        for (int j = 0; j < thislevel.size(); j++)
        {
            int v;
            v = thislevel[j];
            if (std::find(seen.begin(), seen.end(), v) == seen.end())
            {
                seen.push_back(v);
                std::vector<int> vneigh;
                vneigh = neighbors(v);
                nextlevel.insert(nextlevel.end(), vneigh.begin(), vneigh.end());
            }
        }
    }
    return seen;

}

std::vector<std::vector<int>> mgraph::connected_components() // algorithm the same as that in networkx python
{
    std::vector<int> seen;
    std::vector<int> nodes_vec;
    nodes_vec = get_nodes();
    std::vector<std::vector<int>> cc;
    for (int i = 0; i< nodes_vec.size(); i++)
    {
        int node_label = nodes_vec[i];
        if (std::find(seen.begin(), seen.end(), node_label) == seen.end()) // node_label not in seen
        {
            std::vector<int> bfs_comp;
            bfs_comp = _plain_bfs(node_label);
            seen.insert(seen.end(), bfs_comp.begin(), bfs_comp.end());
            cc.push_back(bfs_comp);
        }
    }
    return cc;
}