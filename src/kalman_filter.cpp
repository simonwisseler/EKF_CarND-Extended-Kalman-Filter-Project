#include <math.h>

#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &Q_in, MatrixXd &H_in, MatrixXd &R_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    Q_ = Q_in;
    H_ = H_in;
    R_ = R_in;
}

void KalmanFilter::Predict() {
    /**
     * NOTE: LINEAR motion model
     */
    x_ = F_*x_;
    MatrixXd F_t = F_.transpose();
    P_ = F_*P_*F_t + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    MatrixXd H_t = H_.transpose();
    MatrixXd S = H_*P_*H_t + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_*H_t*Si;
    
    VectorXd z_pred = H_*x_;
    VectorXd innovation = z - z_pred;
    
    x_ = x_ + (K*innovation);
    MatrixXd I = MatrixXd::Identity(4,4);
    P_ = (I - K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    MatrixXd H_t = H_.transpose();
    MatrixXd S = H_*P_*H_t + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_*H_t*Si;
    
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    
    VectorXd z_pred(3);
    float rho = sqrt(px*px + py*py);
    float phi = atan2(py, px);
    float rho_dot;
    
    if (px == 0 && py == 0) {
        rho_dot = 0;
    } else {
        rho_dot = (px*vx + py*vy)/rho;
    }
    z_pred << rho, phi, rho_dot;
    
    VectorXd innovation = z - z_pred;
    
    x_ = x_ + (K*innovation);
    MatrixXd I = MatrixXd::Identity(4,4);
    P_ = (I - K*H_)*P_;
}
