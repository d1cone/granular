#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include "system.h"

// System methods

void System::ComputeForces() {
    const auto& head = nlist.getHead();
    const auto& neighbors = nlist.getNeighbors();

    double k = 4.3;
    double f_coeff = (k - 1) * (k - 1)/9;
    
    bool use_hertz = std::find(interactions.begin(), interactions.end(), InteractionType::Hertzian) != interactions.end();
    bool use_coulomb = std::find(interactions.begin(), interactions.end(), InteractionType::Coulomb) != interactions.end();
    bool use_polarized = std::find(interactions.begin(), interactions.end(), InteractionType::Polarizable) != interactions.end();

    // get pair_coefficients
    int pair_idx = 0;

    // hertzian
    double kn = use_hertz ? pair_coeff[pair_idx++] : 0.0;
    double kt = use_hertz ? pair_coeff[pair_idx++] : 0.0;
    double gamman = use_hertz ? pair_coeff[pair_idx++] : 0.0;
    double gammat = use_hertz ? pair_coeff[pair_idx++] : 0.0;

    // coulomb
    double k_coulomb = pair_coeff.size() > 4 ? pair_coeff[pair_idx++] : 0.0; 

    // reset potential energy calculation
    U = 0;

    // iterate through particles
    for (int i = 0; i < particles.size(); i++) {
        int start = head[i];
        int end = head[i + 1];
        
        for (int idx = start; idx < end; idx++) {
            int j = neighbors[idx];
            
            Eigen::Vector3d r_vec = particles[i].x - particles[j].x;
            double r_sq = r_vec.squaredNorm();
            
            if (r_sq <= nlist.getCutoff() * nlist.getCutoff() && r_sq > 0.0) {
                double r = std::sqrt(r_sq);
                Eigen::Vector3d n_vec = r_vec / r; 
                Eigen::Vector3d force_ij = Eigen::Vector3d::Zero();
                
                // Hertzian Interaction
                if (use_hertz) {
                    double overlap = (particles[i].r + particles[j].r) - r;
                    
                    if (overlap > 0.0) {

                        Eigen::Vector3d v_ij = particles[i].v - particles[j].v;
                        double vn = v_ij.dot(n_vec); 
                        
                        // Normal Force
                        double fn_spring = kn * std::pow(overlap, 1.5);
                        double fn_damping = -gamman * vn;
                        double fn_total = fn_spring + fn_damping;
                        
                        fn_total = std::max(0.0, fn_total); 
                        
                        force_ij += fn_total * n_vec;
                        
                        // Tangential Force: Damping only
                        // Eigen::Vector3d vt_vec = v_ij - (vn * n_vec);
                        // force_ij += -gammat * vt_vec;

                        U += (2.0 / 5.0) * kn * std::pow(overlap, 2.5);
                    }
                }
                
                // Coulombic Interaction
                if (use_coulomb) {
                    double f_mag = k_coulomb * (particles[i].q * particles[j].q) / r_sq;
                    force_ij += f_mag * n_vec;
                    U += k_coulomb * (particles[i].q * particles[j].q) / r;
                }
                
                // Polarizable interaction
                if (use_polarized) {
                    double y = f_coeff * std::pow(particles[i].r*particles[j].r/r_sq, 3);
                    double f_mag = k_coulomb / std::pow((1 - 4*y)*r, 2) * (particles[i].q*particles[j].q*(1 + 2*y*(3 + 4*y)) - 2*((k - 1)/3 * std::pow(particles[i].r/r, 3)*particles[j].q*particles[j].q + (k - 1)/3 * std::pow(particles[j].r/r, 3)*particles[i].q*particles[i].q)*(1 + 2*y));
                    force_ij += f_mag * n_vec;
                    U += k_coulomb * (particles[i].q * particles[j].q) / r;
                }

                // Update forces
                particles[i].f += force_ij;
                particles[j].f -= force_ij; 
            }
        }
    }
    // Enforce boundary conditions, if they exist
    if (bounds != nullptr) {
        bounds->ApplyBoundaryForces(particles, pair_coeff);
    }
}

// Checks the neighbor list and rebuilds every timestep
void System::CheckAndRebuildNeighbors() {
    if (!init) throw std::logic_error("Neighbors cannot be built before the system is initialized");

    double max_disp_sq = 0.0;

    // Find the max displacement
    for (int i = 0; i < particles.size(); i++) {
        double disp_sq = (particles[i].x - x_at_last_build[i]).squaredNorm();
        if (disp_sq > max_disp_sq) {
            max_disp_sq = disp_sq;
        }
    }

    double trigger_dist = nlist.getSkinCutoff() / 2.0;
    
    // If the max diplacement is greater than the trigger distance squared, rebuild the neighbor list
    if (max_disp_sq > (trigger_dist * trigger_dist)) {
        BuildNeighborList();
    }
}

bool System::GenerateRandom(int n_particles, int maxtry, std::vector<double> pairc) {
    if (init) { 
        std::cout << "System already initialized\n";
        return false;
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> rand(-9.0, 9.0);

    for (int i = 1; i <= n_particles; i++) {
        int t = 0;
        bool placed = false;
        
        while (t < maxtry) {
            Eigen::Vector3d coords;
            if (!bounds) {
                coords = {rand(gen), rand(gen), rand(gen)};
            } else {
                coords = bounds->GenerateRandomCoord();
            }
            
            bool failed = false;
            for (const auto& p : particles) {
                if ((std::abs(coords[0] - p.x[0]) < p.r * 2.0) && 
                    (std::abs(coords[1] - p.x[1]) < p.r * 2.0) && 
                    (std::abs(coords[2] - p.x[2]) < p.r * 2.0)) { 
                    failed = true;
                    break;
                }
            }
            
            if (!failed) {
                Particle p;
                p.id = i;
                if (i % 2 == 1) p.type = 1;
                else p.type = 2;
                p.x = coords;
                p.v = Eigen::Vector3d::Zero(); 
                p.f = Eigen::Vector3d::Zero();
                p.r = 0.5;  
                p.m = 1.0;
                if (i % 2 == 1) p.q = 1;
                else p.q = -1;
                
                particles.push_back(p);
                placed = true;
                break; 
            }
            t++;
        }
        if (!placed) std::cout << "Warning: Could not place particle " << i << "\n";
    }
    
    pair_coeff = pairc; 
    init = true;
    return true;
}

void System::Step() {
    if (!init) throw std::logic_error("Simulation cannot be run before generation");


    for (auto& p : particles) {
        p.v += 0.5 * dt * (p.f / p.m);
        
        p.x += p.v * dt;
        
        p.f.setZero(); 
    }

    CheckAndRebuildNeighbors();

    ComputeForces();

    ComputeK();

    for (auto& p : particles) {
        p.v += 0.5 * dt * (p.f / p.m);
    }
    
    timestep++;
}

// OVITO friendly dump writer (followed dump format from lammps)
void System::WriteTimestep(std::string fname) {
    std::ofstream dump(fname, std::ios_base::app); 
    
    dump << "ITEM: TIMESTEP\n" << timestep << "\n";
    dump << "ITEM: NUMBER OF ATOMS\n" << particles.size() << "\n";
    
    dump << "ITEM: BOX BOUNDS pp pp pp\n";
    dump << "-10.0 10.0\n-10.0 10.0\n-10.0 10.0\n"; 
    
    dump << "ITEM: ATOMS id type x y z radius q vx vy vz\n";
    for (const auto& p : particles) {
        dump << p.id << " " << p.type << " " 
             << p.x[0] << " " << p.x[1] << " " << p.x[2] << " " 
             << p.r << " " << p.q << " "
             << p.v[0] << " " << p.v[1] << " " << p.v[2] << "\n";
    }
}

void System::WriteLogHeader(std::string fname) {
    std::ofstream log(fname);

    log << "timestep | Total Energy | Kinetic Energy | Potential Energy" << std::endl;
}

void System::WriteLog(std::string fname) {
    std::ofstream log(fname, std::ios_base::app);

    log << timestep << "     " << H << "     " << K << "     " << U << "     " << "\n";
}

void System::BuildNeighborList() {
    if (!init) throw std::logic_error("Neighbors cannot be built before the system is initialized");

    nlist.Build(particles);

    x_at_last_build.resize(particles.size());
    for (int i = 0; i < particles.size(); i++) {
        x_at_last_build[i] = particles[i].x;
    }
}

void System::ComputeK() {
    K = 0;
    for (const auto& p : particles) {
        K += 0.5 * p.m * p.v.squaredNorm();
    }

    H = K + U;
}