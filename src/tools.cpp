#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    /**
    TODO:
    * Calculate the RMSE here.
    */
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean
    rmse = rmse/estimations.size();

    //calculate the squared root
    rmse = rmse.array().sqrt();

    //return the result
    return rmse;
}

VectorXd Tools::CalculateStates(const MeasurementPackage &measurement_pack) {
    float rho = measurement_pack.raw_measurements_[0];
    float phi = measurement_pack.raw_measurements_[1];
    float rhodot = measurement_pack.raw_measurements_[2];

    VectorXd x(4);
    x << rho*cos(phi),rho*sin(phi),rhodot*cos(phi),rhodot*sin(phi);

    return x;
}

MatrixXd Tools::CalculateObservation(const Eigen::VectorXd &x_state) {

    VectorXd hx(3);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //TODO: YOUR CODE HERE

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3;
    if(fabs(px) < 0.0001){
        c3 = 3.1415926/2;
    }
    else{
        c3 = atan(py/px);
    }
    float c4 = px*vx+py*vy;
    float c5 = vx*vx + vy*vy;
    float c6;
    if(fabs(c1) < 0.0001){
        c6 = sqrt(c5);
    }
    else{
        c6 = c4/c2;
    }

    hx << c2,c3,c6;

    return hx;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    /**
    TODO:
    * Calculate a Jacobian here.
    */
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    //TODO: YOUR CODE HERE

    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        //return the null Jacobian matrix
        Hj << 0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0;
    }
    else{
        //compute the Jacobian matrix
        Hj << (px/c2), (py/c2), 0, 0,
                -(py/c1), (px/c1), 0, 0,
                py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    }

    return Hj;
}

MatrixXd Tools::CalculateF(float dt) {
    MatrixXd F(4,4);
    F << 1, 0, dt, 0,
         0, 1, 0, dt,
         0, 0, 1, 0,
         0, 0, 0, 1;

    return F;
}

MatrixXd Tools::CalculateQ(float dt, float noise_ax, float noise_ay) {
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    //set the process covariance matrix Q
    MatrixXd Q(4, 4);
    Q << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

    return Q;
}
