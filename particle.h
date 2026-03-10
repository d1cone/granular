#pragma once
#include <Eigen/Dense>

// values stores by each particle
struct Particle {
    int id;
    int type;
    Eigen::Vector3d x = Eigen::Vector3d::Zero();
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    Eigen::Vector3d f = Eigen::Vector3d::Zero();
    double r;
    double q;
    double m;
};