#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
        std::cout<<"estimations size error"<<std::endl;
        return rmse;
    }
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd res = 	estimations[i] - ground_truth[i];
        res = res.array()*res.array();
        rmse = rmse + res;
    }

    //calculate the mean
    rmse = rmse/ground_truth.size();
    //calculate the squared root
    rmse = rmse.array().sqrt();
    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
    MatrixXd Hj(3,4);
    Hj << 0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0;
    //recover state parameters
    const double& px = x_state(0);
    const double& py = x_state(1);
    const double& vx = x_state(2);
    const double& vy = x_state(3);

    double pxy = px*px + py*py;
    double pxysqrt = sqrt(pxy);
    double pxy32 = pxy*pxysqrt;
    //check division by zero
    if(fabs(pxy) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
        return Hj;
    }



    //compute the Jacobian matrix
    Hj <<   px/pxysqrt, py/pxysqrt, 0, 0,
            -py/pxy,    px/pxy,     0, 0,
            py*(vx*py-vy*px)/pxy32, px*(vy*px-vx*py)/pxy32, px/pxysqrt, py/pxysqrt;

    return Hj;
}
