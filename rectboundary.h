#include <vector>
#include "boundary.h"

class RectBoundary : public Boundary {
protected:
    std::vector<double> bounds;
public:
    RectBoundary(double minusx, double plusx, double minusy, double plusy, double minusz, double plusz);
    bool Contains(Particle p);
    Eigen::Vector3d GenerateRandomCoord() override;
    double Overlap(Particle p) override;
    void ApplyBoundaryForces(std::vector<Particle>& particles, const std::vector<double>& pair_coeff) override;
};