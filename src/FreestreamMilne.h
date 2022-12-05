#pragma once
//This file contains a wrapper class for freestream-milne

#include <vector>

using namespace std;

class FREESTREAMMILNE {
 private:

 public:
    FREESTREAMMILNE();
    ~FREESTREAMMILNE();

    int run_freestream_milne(const char * inputPath = nullptr);

    // IS THIS VARIABLE NECESSARY
    int gridSize; //the total number of grid points in x, y, and eta : used for vector memory allocation

    //the total freestreaming time, useful for passing to JETSCAPE
    float tau_LandauMatch;

    //support to initilialize the energy density from a vector - useful for JETSCAPE
    //note units of argument should be GeV / fm^3
    //then we convert to fm^(-4)
    void initialize_from_vector(std::vector<float>);
    std::vector<float> init_energy_density;

    //support to write final hydro variables to vectors - useful for JETSCAPE
    //note we need to convert back to GeV / fm^3 units here
    void output_to_vectors(std::vector<double>&, //e
                            std::vector<double>&, //p
                            std::vector<double>&, //ut
                            std::vector<double>&, //ux
                            std::vector<double>&, //uy
                            std::vector<double>&, //un
                            std::vector<double>&, //pitt
                            std::vector<double>&, //pitx
                            std::vector<double>&, //pity
                            std::vector<double>&, //pitn
                            std::vector<double>&, //pixx
                            std::vector<double>&, //pixy
                            std::vector<double>&, //pixn
                            std::vector<double>&, //piyy
                            std::vector<double>&, //piyn
                            std::vector<double>&, //pinn
                            std::vector<double>&); //Pi

    std::vector<double> final_energy_density;
    std::vector<double> final_pressure;
    std::vector<double> final_ut;
    std::vector<double> final_ux;
    std::vector<double> final_uy;
    std::vector<double> final_un;
    std::vector<double> final_pitt;
    std::vector<double> final_pitx;
    std::vector<double> final_pity;
    std::vector<double> final_pitn;
    std::vector<double> final_pixx;
    std::vector<double> final_pixy;
    std::vector<double> final_pixn;
    std::vector<double> final_piyy;
    std::vector<double> final_piyn;
    std::vector<double> final_pinn;
    std::vector<double> final_Pi;

};
