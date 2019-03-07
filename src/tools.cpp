#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;
    
    for(unsigned i = 0; i < estimations.size(); i++){
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array()*residual.array(); //element-wise multiplication
        rmse += residual;
    }
    rmse = rmse/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /**
     * Jacobian used to calculate predicted radar measurements - first-order Taylor expansion
     * required to map state x from Cartesian coordinates px, py, vx, vy to Polar coordinates rho, phi, rho_dot
     *
     * NOTE: For mathematical derivation of range rate rho_dot, see https://en.wikipedia.org/wiki/Range_rate
     * NOTE: To arrive at partial derivatives of rho_dot with respect to px, py,
     * use the product rule and expand the resulting first summand/fraction by multiplying
     * numerator and denominator by (px^2+py^2)
     */
    
    MatrixXd Hj(3, 4);
    
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    if(px == 0 && py == 0) {
        std::cout << "Error: Division by 0 in Tools::CalculateJacobian" << std::endl;
        return Hj;
    }
    float rho = std::sqrt(px*px+py*py);
    float rho_squared = rho*rho;
    float rho_cubed = rho_squared*rho;
    Hj << px/rho, py/rho, 0, 0,
        -py/rho_squared, px/rho_squared,0, 0,
        py*(vx*py-vy*px)/rho_cubed, px*(vy*px-vx*py)/rho_cubed, px/rho, py/rho;
    return Hj;
    
}
