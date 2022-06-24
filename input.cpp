#include "input.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<filesystem>

void input_f::set_fragSize(int size)
{
    fragSize = size;
}

int input_f::get_fragSize()
{
    return fragSize;
}

void input_f::set_coordfile(std::string coordPath)
{
    coordfile = coordPath;
}

std::string input_f::get_coordfile()
{
    return coordfile;
}

std::string readFile(std::string &input)
{
    std::ifstream inputfile;
    std::string inputstring;
    std::string line;

    inputfile.open(input);
    if (inputfile.fail())
    {
        std::cout << "error with " << input << " file" << std::endl;
        exit (1);
    }
    while (getline(inputfile, line))
    {
        std::string newline;
        newline = " " + line;
        inputstring += newline;         
    }
    return inputstring; 
}

input_f get_inputs(std::string &input_contents)
{
    input_f input_obj;
    std::string input_string; 
    std::stringstream inputstream(input_contents);
    while (inputstream >> input_string)
    {
        if (input_string == "coordinates")
        {
            std::string filename;
            inputstream >> filename;
            input_obj.set_coordfile(filename);
        }

        if (input_string == "fragSize")
        {
            std::string fragSize;
            inputstream >> fragSize;
            input_obj.set_fragSize(stoi(fragSize));
        }

        
    }
    return input_obj;
}

// void extract_coordinates(std::string coordpath, )

input_f process_input(std::string inputfile)
{
    std::string input_contents;
    input_f input_obj;
    input_contents = readFile(inputfile); 
    input_obj = get_inputs(input_contents);
return input_obj;
}

void aType::set_r1(double &radius)
{
    r1 = radius;
}

void aType::set_theta0(double &theta)
{
    theta0 = theta;
}

void aType::set_x1(double &x)
{
    x1 = x;
}

void aType::set_D1(double &D)
{
    D1 = D;
}

void aType::set_zeta(double &z)
{
    zeta = z;
}

void aType::set_Z1(double &Z)
{
    Z1 = Z;
}

void aType::set_Vi(double &V)
{
    Vi = V;
}

void aType::set_Uj(double &U)
{
    Uj = U;
}

void aType::set_Xi(double &chi)
{
    Xi = chi;
}

double aType::get_r1()
{
    return r1;
}

double aType::get_theta0()
{
    return theta0;
}

double aType::get_x1()
{
    return x1;
}

double aType::get_D1()
{
    return D1;
}

double aType::get_zeta()
{
    return zeta;
}

double aType::get_Z1()
{
    return Z1;
}

double aType::get_Vi()
{
    return Vi;
}

double aType::get_Uj()
{
    return Uj;
}

double aType::get_Xi()
{
    return Xi;
}
