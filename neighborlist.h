#include <vector>
#include "particle.h"

class NeighborList {
protected:
    double cutoff; 
    double skin_cutoff;
    double ghost_cutoff; // cutoff + skin
    std::vector<int> head; // each index in the head represents a particle, and the number
    // in the index gives the index of where the neighborlist of neighbors starts 
    std::vector<int> neighbors; // implemented using a flat array to save computing power

public:
    NeighborList(int c, int cskin) : cutoff(c), skin_cutoff(cskin) {
        ghost_cutoff = c + cskin;
    }
    double getCutoff() { return cutoff; }
    const std::vector<int>& getHead() const { return head; }
    const std::vector<int>& getNeighbors() const { return neighbors; }
    double getSkinCutoff() const { return skin_cutoff; } 
    void Build(const std::vector<Particle>& particles);
    void ClearList();
};