#ifndef INPUT_H
#define INPUT_H
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<filesystem>

class input_f
{
    int fragSize;
    std::string coordfile;

    public:
    void set_fragSize(const int fragSize);
    int get_fragSize();
    void set_coordfile(std::string coordfile);
    std::string get_coordfile();
    
};

// input_f process_input(std::string inputfile);
std::string readFile(std::string &input);
input_f get_inputs(std::string &input_contents);
input_f process_input(std::string inputfile);

class aType
{
    double r1;
    double theta0;
    double x1;
    double D1;
    double zeta;
    double Z1;
    double Vi;
    double Uj;
    double Xi;

    public:
    void set_r1(double &r1);
    void set_theta0(double &theta0);
    void set_x1(double &x1);
    void set_D1(double &D1);
    void set_zeta(double &zeta);
    void set_Z1(double &Z1);
    void set_Vi(double &Vi);
    void set_Uj(double &Uj);
    void set_Xi(double &Xi);

    // getter functions
    double get_r1();
    double get_theta0();
    double get_x1();
    double get_D1();
    double get_zeta();
    double get_Z1();
    double get_Vi();
    double get_Uj();
    double get_Xi();
};

#endif