#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    previous_timestamp_ = 0;

    R_laser_ = MatrixXd(2,2);
    R_laser_ << 0.0225, 0,
                0, 0.0225;
    
    R_radar_ = MatrixXd(3,3);
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;
    
    H_laser_ = MatrixXd(2,4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;
    
    Hj_ = MatrixXd(3,4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
     * Initialization
     */
    if (!is_initialized_) {
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ <<  1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1000, 0,
                    0, 0, 0, 1000;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x_(0) = measurement_pack.raw_measurements_(0);
            ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            float rho     = measurement_pack.raw_measurements_(0);
            float phi    = measurement_pack.raw_measurements_(1);
            float rho_dot = measurement_pack.raw_measurements_(2);
            ekf_.x_(0) = rho     * cos(phi);
            ekf_.x_(1) = rho     * sin(phi);
            // best estimate for vx, vy as component of velocity normal to rho not captured by radar measurement
            ekf_.x_(2) = rho_dot * cos(phi);
            ekf_.x_(3) = rho_dot * sin(phi);
        }
        
        previous_timestamp_ = measurement_pack.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    
    /**
    * Prediction
    */
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, dt, 0,
                0, 1, 0, dt,
                0, 0, 1, 0,
                0, 0, 0, 1;
    
    
    double noise_ax = 9.0; // random noise (acceleration in x) provided by task
    double noise_ay = 9.0; // random noise (acceleration in y) provided by task
    
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    double dt_4_4 = dt_4 / 4;
    double dt_3_2 = dt_3 / 2;
    
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
                0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
                dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
                0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;
    
    ekf_.Predict();

    
    /**
    * Update
    */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    
    
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
