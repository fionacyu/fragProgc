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