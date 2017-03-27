#include "kalman_filter.h"
#include <iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();

  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */

    VectorXd z_pred(3);
    double& px = x_[0];
    double& py = x_[1];
    double pxy_sqrt = sqrt(px*px + py*py);
    double& vx = x_[2];
    double& vy = x_[3];
    z_pred << pxy_sqrt, 0, 0;
    // deal with py = 0 or px = 0 for phi
    if(py == 0 && px == 0){ // at origin
        z_pred[1] = 0;
    }
    else if(py == 0){ // at x axis
        if(px > 0)
            z_pred[1] = 0;
        else
            z_pred[1] = M_PI;
    }
    else if(px == 0){ // at y axis
        if(py > 0)
            z_pred[1] = M_PI / 2;
        else
            z_pred[1] = - M_PI / 2;
    }
    else{ // normal case: px > 0, py > 0
        z_pred[1] = atan2(py, px);
    }


    // ensure not divide by zero
    if(pxy_sqrt != 0)
        z_pred[2] = (px*vx + py*vy)/pxy_sqrt;

    VectorXd y = z - z_pred;
    //set phi of y in range of -pi to pi
    if(fabs(y[1]) > M_PI){
        double phi = 2*M_PI - fabs(y[1]);
        y[1] = y[1] > 0 ? -phi : phi;
    }


    MatrixXd Ht = H_.transpose();
    MatrixXd PHt = P_ * Ht;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = PHt * Si;


    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
