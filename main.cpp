#include <iostream>
#include <vector>
#include <fstream>
#include "system.h"

int main() {
    std::string run_name = "polarizable_q=pm1";
    std::string dump_name = "./dumps/" + run_name;
    std::string log_name = "./logs/" + run_name;

    std::vector<InteractionType> interactions = {InteractionType::Hertzian, InteractionType::Polarizable};
    
    std::vector<double> pair_coeffs = {1000.0, 0.0, 50.0, 0.0, 1.0}; 

    RectBoundary bounds = {-5,5,-5,5,-5,5};
    RectBoundary* b = &bounds;

    double dt = 0.001;
    double cutoff = 5; 
    double skin = 6;

    System sys(interactions, dt, cutoff, skin, b);
    
    sys.GenerateRandom(250, 10000, pair_coeffs); 
    
    sys.BuildNeighborList(); 
    
    sys.WriteLogHeader(log_name);

    std::ofstream clear_dump(dump_name);
    clear_dump.close();

    int num_steps = 1000000;
    for (int step = 0; step < num_steps; step++) {
        sys.Step();
        
        if (step % 10 == 0) {
            sys.WriteLog(log_name);
        }

        if (step % 100 == 0) {
            sys.WriteTimestep(dump_name);
            std::cout << "Step: " << step << " completed." << std::endl;
        }
    }

    return 0;
}