#ifndef MGRAPH_H
#define MGRAPH_H
#include<vector>
#include<string>
#include<map>
#include<array>

class node
{
    std::string element;
    int coordno;
    std::vector<double> box;
    std::vector<double> coord;
    std::string at;
    int pi; // number of pi electrons contributing to conjugated system

    public:
    node() //defaults
    {
        coordno = 1; // set to 1 for the hydrogens
        coord = {0.0, 0.0, 0.0};
        box = {0.0, 0.0, 0.0};
        at = "at";
        pi = 0;
    }

    void set_element(std::string &element);
    void set_coordno(int &coordno);
    void set_coordinates(std::vector<double> &coord);
    void set_box(std::vector<double> &box);
    void set_at(std::string &at);
    void set_pi(int &pi);

    // getter functions
    std::string get_element();
    int get_coordno();
    std::vector<double> get_coordinates();
    std::vector<double> get_box();
    std::string get_at();
    int get_pi();
};

class edge
{
    double bo;
    bool conjugated;
    double bl;

    public:
    edge() //defaults
    {
        conjugated = false;
    }
    void set_bo(double &bo);
    void set_conjstatus(bool &conjugated);
    void set_bl(double &bl);

    // getter functions
    double get_bo();
    bool get_conjstatus();
    double get_bl();
};

class mgraph // molecular graph
{
    std::vector<int> nodeVec; // vector containing the node labels
    std::vector<std::array<int,2>> edgeVec; // vector containing the edges
    std::map<int, node> mnodes; //nodes represented by integer
    std::map<std::array<int,2>, edge> medges; // edges represented by array of size 2

    public:
    void add_node(int &nodelabel, node &n); // insert nodes individually 
    void add_edge(std::array<int,2> &edgelabel, edge &e); // insert edge individually 
    void update_node(int &nodelabel, node &n);
    void update_edge(std::array<int,2> &edgelabel, edge &e);

    std::vector<int> get_nodes();
    std::vector<std::array<int,2>> get_edges();
    node nodes(int &nlabel);
    edge edges(std::array<int,2> &elabel);
    std::vector<int> neighbors(int &nodelabel);
    std::vector<int> _plain_bfs(int &node_label);
    std::vector<std::vector<int>> connected_components();

};

#endif