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