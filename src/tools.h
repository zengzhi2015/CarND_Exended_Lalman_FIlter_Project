#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
    * A helper method to calculate RMSE.
    */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

    /**
    * A helper method to transform the polar coord. meas. to Cardissain coord. meas.
    */
    Eigen::VectorXd CalculateStates(const MeasurementPackage &measurement_pack);

    /**
    * A helper method to calculate non-linear observation.
    */
    Eigen::MatrixXd CalculateObservation(const Eigen::VectorXd& x_state);

    /**
    * A helper method to calculate Jacobians.
    */
    Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

    /**
    * A helper method to calculate process matrix.
    */
    Eigen::MatrixXd CalculateF(float dt);

    /**
    * A helper method to calculate process covarience matrix.
    */
    Eigen::MatrixXd CalculateQ(float dt, float noise_ax, float noise_ay);
};

#endif /* TOOLS_H_ */
