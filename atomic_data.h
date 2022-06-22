#ifndef ATOMIC_H
#define ATOMIC_H
#include<map>
#include<cmath>
#include<string>
#include<vector>

const std::map<std::string, int>valence=
{
    { "H", 1 },
    { "B", 3 },
    { "C", 4 },
    { "N", 5 },
    { "O", 6 },
    { "F", 7 },
    { "P", 5 },
    { "S", 6 },
    { "Cl", 7 },
    { "Br", 7 },
    { "I", 7 }
};

// vdW radii provided by Bondi
const std::map<std::string, double>vdwradii=
{
    { "H" , 1.20 },
    { "He", 1.40 },
    { "C" , 1.70 },
    { "N" , 1.55 },
    { "O" , 1.52 },
    { "F" , 1.47 },
    { "Ne", 1.54 },
    { "Si", 2.10 },
    { "P" , 1.80 },
    { "S" , 1.80 },
    { "Cl", 1.75 },
    { "Ar", 1.88 },
    { "As", 1.85 },
    { "Se", 1.90 },
    { "Br", 1.85 },
    { "Kr", 2.02 },
    { "Te", 2.06 },
    { "I" , 1.98 },
    { "Xe", 2.16 }
};

// cov radii taken from Cordero, Beatriz, et al. "Covalent radii revisited." Dalton Transactions 21 (2008): 2832-2838
const std::map<std::string, double>covradii=
{
    { "H", 0.31 },
    { "B", 0.84 },
    { "C4", 0.76 }, 
    { "C3", 0.73 },
    { "C2", 0.69 },
    { "N", 0.71 },
    { "O", 0.66 },
    { "F", 0.57 },
    { "P", 1.07 },
    { "S", 1.05 },
    { "Cl", 1.02 },
    { "Br", 1.20 },
    { "I", 1.39 }
};

const std::map<std::vector<std::string>, std::map<double, std::vector<double>>>bondr=
{
    { { "C", "C" }, { {1.00, {1.370, 1.596}}, {1.5, {1.370, 1.432}}, {2.00, {1.243, 1.382}}, {3.00, {1.187, 1.268}} } },
    { { "C", "Br" }, { {1.00, {1.789, 1.950}} } },
    { { "C", "Cl" }, { {1.00, {1.612, 1.813}} } },
    { { "C", "F" }, { {1.00, {1.262, 1.401}} } },
    { { "C", "I" }, { {1.00, {1.992, 2.157}} } },
    { { "C", "N" }, { {1.00, {1.347, 1.492}}, {1.5, {1.328, 1.350}}, {2.00, {1.207, 1.338}}, {3.00, {1.14, 1.177}} } },
    { { "C", "O" }, { {1.00, {1.273, 1.448}}, {2.00, {1.135, 1.272}}, {3.00, {1.115, 1.145}} } },
    { { "C", "P" }, { {1.00, {1.858, 1.858}}, {2.00, {1.673, 1.673}}, {3.00, {1.542, 1.562}} } },
    { { "C", "S" }, { {1.00, {1.714, 1.849}}, {2.00, {1.553, 1.647}}, {3.00, {1.478, 1.535}} } },
    { { "O", "O" }, { {1.00, {1.116, 1.516}}, {2.00, {1.2, 1.208}} } },
    { { "N", "O" }, { {1.00, {1.184, 1.507}}, {2.00, {1.066, 1.258}} } },
    { { "O", "S" }, { {2.00, {1.405, 1.5}} } },
    { { "C", "H" }, { {1.00, {0.931, 1.14}} } },
    { { "N", "H" }, { {1.00, {0.836, 1.09}} } },
    { { "S", "H" }, { {1.00, {1.322, 1.4}} } },
    { { "O", "H" }, { {1.00, {0.912, 1.033}} } },
    { { "F", "H" }, { {1.00, {0.917, 1.014}} } },
    { { "N", "N" }, { {1.00, {1.181, 1.864}}, {1.5, {1.332, 1.332}}, {2.00, {1.139, 1.252}}, {3.00, {1.098, 1.133}} } },
    { { "S", "S" }, { {1.00, {1.89, 2.155}}, {2.00, {1.825, 1.898}} } },
    { { "H", "H" }, { {1.00, {0.741, 0.741}} } },
    { { "F", "O" }, { {1.00, {1.421, 1.421}} } },
    { { "F", "F" }, { {1.00, {1.322, 1.412}} } },
    { { "Cl", "H" }, { {1.00, {1.275, 1.321}} } },
    { { "Cl", "O" }, { {1.00, {1.641, 1.704}}, {2.00, {1.404, 1.414}} } },
    { { "Cl", "Cl" }, { {1.00, {1.9879, 1.9879}} } },
    { { "N", "F" }, { {1.00, {1.317, 1.512}} } }

};

const double pi = 4.0*atan(1.0);
const double KCAL_TO_KJ = 4.184;
const double DEG_TO_RADIAN = pi/180.00;
const double RAD_TO_DEG = 180.0/pi;


#endif