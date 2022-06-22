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

#endif