#include<vector>
#include "particle.h"

// abstract class that all boundaries inherit from
class Boundary {
public:
    virtual Eigen::Vector3d GenerateRandomCoord() = 0;
    virtual double Overlap(Particle p) = 0; // return the overlap of a given particle with the wall
    virtual void ApplyBoundaryForces(std::vector<Particle>& particles, const std::vector<double>& pair_coeff) {}
};