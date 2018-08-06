#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.001 // Just a small number

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5,5);

  //set state dimension
  n_x = 5;

  //set augmented dimension
  n_aug = 7;

  n_sig = 2 * n_aug + 1;

  weights = VectorXd(n_sig);

  lambda_ = 3 - n_aug;

  weights(0) = lambda_ / (lambda_ + n_aug);
  for (int i = 1; i < weights.size(); i++) {
        weights(i) = 0.5 / (n_aug + lambda_);
  }

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.57;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;


 ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // Measurement noise covariance matrices initialization
  R_radar = MatrixXd(3, 3);
  R_radar << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0,std_radrd_*std_radrd_;
  R_lidar = MatrixXd(2, 2);
  R_lidar << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

}

UKF::~UKF() {}

void UKF::Norm(double *ang) {
    while (*ang > M_PI) *ang -= 2. * M_PI;
    while (*ang < -M_PI) *ang += 2. * M_PI;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
   if (!is_initialized_) {

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */

              double rho = measurement_pack.raw_measurements_(0);
              double phi = measurement_pack.raw_measurements_(1);
              x_(0) = rho*cos(phi);
              x_(1) = rho*sin(phi);
              x_(2) = 0 ;
              x_(3) = 0 ;
              x_(4) = 0;

          }
          else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_(0) = measurement_pack.raw_measurements_(0);
            x_(1) = measurement_pack.raw_measurements_(1);
            x_(2) = 0 ;
            x_(3) = 0 ;
            x_(4) = 0;

            if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
		            x_(0) = EPS;
		              x_(1) = EPS;
	          }
          }

          // done initializing, no need to predict or update
          is_initialized_ = true;
          previous_timestamp_ = measurement_pack.timestamp_;

          return;
        }

     double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
     previous_timestamp_ = measurement_pack.timestamp_;

     Prediction(dt);

     if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
       UpdateRadar(measurement_pack);
     }else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
       UpdateLidar(measurement_pack);
     }

}




/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  AugmentedSigmaPoints(&augPoints);
  SigmaPointPrediction(&Xsig_pred,augPoints,delta_t);
  PredictMeanAndCovariance(&x_,&P_,Xsig_pred);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;
  // Create matrix for sigma points in measurement space
  // Transform sigma points into measurement space
  MatrixXd Zsig = Xsig_pred.block(0, 0, n_z, n_sig);
  UpdateUKF(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  int n_z = 3;
  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig);
  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig; i++) {
    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    // Measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);          //r
    Zsig(1,i) = atan2(p_y,p_x);                   //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / Zsig(0,i);   //r_dot
  }
  UpdateUKF(meas_package, Zsig, n_z);
}



void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {



  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  double sqrt_lambda_n_aug = sqrt(lambda_+n_aug);

  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt_lambda_n_aug*L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt_lambda_n_aug*L.col(i);
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;
}



void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd Xsig_aug,double delta_t) {

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > EPS) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }
  *Xsig_out = Xsig_pred;
}



void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out,MatrixXd Xsig_pred) {

  //create vector for predicted state
  VectorXd x = *x_out;

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);


/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predicted state mean
  x = Xsig_pred * weights;

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    Norm(&(x_diff(3)));
    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}



void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred  = Zsig * weights;
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig; i++) {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // Angle normalization
    Norm(&(z_diff(1)));
    S = S + weights(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    R = R_radar;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
    R = R_lidar;
  }
  S = S + R;

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
      // Angle normalization
      Norm(&(z_diff(1)));
    }
    // State difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    // Angle normalization
    Norm(&(x_diff(3)));
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  // Measurements
  VectorXd z = meas_package.raw_measurements_;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // Residual
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    // Angle normalization
    Norm(&(z_diff(1)));
  }
  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate NIS
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
	NIS_radar = z.transpose() * S.inverse() * z;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){ // Lidar
	NIS_laser = z.transpose() * S.inverse() * z;
  }
}
