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
        residual = residual.array()*residual.array() //element-wise multiplication
        rmse += residual;
    }
    rmse = residual/estimations.size();
    rmse = rmse.array().sqrt();
    return rmse
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
}
