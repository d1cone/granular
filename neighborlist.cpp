#include "neighborlist.h"

// 
void NeighborList::Build(const std::vector<Particle>& particles) {
    int n = particles.size();
    head.assign(n + 1, 0); 
    neighbors.clear();    
    
    double cut_sq = ghost_cutoff * ghost_cutoff; 
    
    for (int i = 0; i < n; i++) {
        head[i] = neighbors.size(); 
        
        for (int j = i + 1; j < n; j++) {
            Eigen::Vector3d r_vec = particles[i].x - particles[j].x;
            
            if (r_vec.squaredNorm() <= cut_sq) {
                neighbors.push_back(j);
            }
        }
    }
    head[n] = neighbors.size(); 
}