#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "neighborlist.h"
#include "rectboundary.h"

// types of interactions allowed 
enum class InteractionType {
    Hertzian, // Pair coeffs: [kn, kt, gamman, gammat]
    Hookean, // TODO
    Coulomb, // Pair coeffs: [cutoff]
    Gravity, // TODO
    Polarizable // Pair coeffs: [cutoff]
};

// System takes care of all large scale processes
class System {
protected:
    int timestep = 0;
    int n;
    double dt;
    Boundary* bounds; // contact boundary conditions
    std::vector<Particle> particles;
    bool init = false; //system must be initialized after creating (using one of the "generate" methods)
    std::vector<InteractionType> interactions;
    ///////////////////
    // pair_coeff format:
    // pair_coeff stores the pair coefficients in the order listed in "InteractionType"
    // IE if hertzian and coulomb are chosen, then the pair coeffs will be [kn, kt, gamman, gammat, coulomb/cutoff]
    ///////////////////
    std::vector<double> pair_coeff;
    NeighborList nlist;
    std::vector<Eigen::Vector3d> x_at_last_build; // previous x values for computing neighbor_list displacement

    // energy
    double U = 0;
    double K = 0;
    double H = 0;

public:
    System(std::vector<InteractionType> interacts, double deltat, double cutoff, double skin) : interactions(interacts), bounds(nullptr), dt(deltat), nlist(cutoff, skin) {}
    System(std::vector<InteractionType> interacts, double deltat, double cutoff, double skin, Boundary* b) : interactions(interacts), bounds(b), dt(deltat), nlist(cutoff, skin) {}
    bool GenerateRandom(int n_particles, int maxtry, std::vector<double> pairc); //  TODO: move pair_coeff generation to the constructor
    bool GenerateLattice(int n_particles, double spacing, std::vector<double> pairc); //TODO
    void CheckAndRebuildNeighbors(); // For checking every timestep
    void ComputeForces();
    void ComputeK();
    void BuildNeighborList(); // Only for the first time building 
    void WriteTimestep(std::string fname); 
    void WriteLogHeader(std::string fname);
    void WriteLog(std::string fname);
    void Step(); // Integration, checking neighbor list, etc
};