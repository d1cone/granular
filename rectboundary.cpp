#include <random>
#include "rectboundary.h"

RectBoundary::RectBoundary(double minusx, double plusx, double minusy, double plusy, double minusz, double plusz) {
    bounds = {minusx, plusx, minusy, plusy, minusz, plusz};
}

bool RectBoundary::Contains(Particle p) {
    if ((p.x[0] >= bounds[0] && p.x[0] <= bounds[1]) && (p.x[1] >= bounds[2] && p.x[1] <= bounds[3]) && (p.x[2] >= bounds[4] && p.x[2] <= bounds[5])) {
        return true;
    }
    return false;
}

double RectBoundary::Overlap(Particle p) {
    Eigen::Vector3d closest_pt = {std::max(bounds[0], std::min(p.x[0], bounds[1])),
    std::max(bounds[2], std::min(p.x[1], bounds[3])),
    std::max(bounds[4], std::min(p.x[2], bounds[5]))
    };

    Eigen::Vector3d disp_vector = p.x - closest_pt;
    return disp_vector.norm();
}

Eigen::Vector3d RectBoundary::GenerateRandomCoord() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    
    double buffer = 0.5; 
    
    std::uniform_real_distribution<double> randx(bounds[0] + buffer, bounds[1] - buffer);
    std::uniform_real_distribution<double> randy(bounds[2] + buffer, bounds[3] - buffer);
    std::uniform_real_distribution<double> randz(bounds[4] + buffer, bounds[5] - buffer);

    return {randx(gen), randy(gen), randz(gen)};
}

void RectBoundary::ApplyBoundaryForces(std::vector<Particle>& particles, const std::vector<double>& pair_coeff) {
    if (pair_coeff.size() < 4) return; // Need kn, gamman, gammat
    
    double kn = pair_coeff[0];
    double gamman = pair_coeff[2];
    double gammat = pair_coeff[3];

    for (auto& p : particles) {
        // Normal vectors pointing INWARD from the 6 walls
        Eigen::Vector3d normals[6] = {
            {1, 0, 0}, {-1, 0, 0},   // minusx, plusx
            {0, 1, 0}, {0, -1, 0},   // minusy, plusy
            {0, 0, 1}, {0, 0, -1}    // minusz, plusz
        };
        
        // Orthogonal distances from particle center to the 6 walls
        double distances[6] = {
            p.x[0] - bounds[0], // dist to minusx
            bounds[1] - p.x[0], // dist to plusx
            p.x[1] - bounds[2], // dist to minusy
            bounds[3] - p.x[1], // dist to plusy
            p.x[2] - bounds[4], // dist to minusz
            bounds[5] - p.x[2]  // dist to plusz
        };

        for (int w = 0; w < 6; w++) {
            double overlap = p.r - distances[w];
            
            // If overlap > 0, the particle is colliding with the wall
            if (overlap > 0.0) {
                Eigen::Vector3d n_vec = normals[w];
                
                // Relative velocity (particle velocity minus stationary wall velocity of 0)
                Eigen::Vector3d v_rel = p.v; 
                double vn = v_rel.dot(n_vec);
                
                // Normal Force (Hertzian Spring + Dashpot)
                double fn_spring = kn * std::pow(overlap, 1.5);
                double fn_damping = -gamman * vn;
                double fn_total = std::max(0.0, fn_spring + fn_damping); // Cannot be attractive
                
                Eigen::Vector3d force_wall = fn_total * n_vec;
                
                // Tangential Force (Damping only)
                Eigen::Vector3d vt_vec = v_rel - (vn * n_vec);
                force_wall += -gammat * vt_vec;
                
                // Apply the force to the particle
                p.f += force_wall;
            }
        }
    }
}