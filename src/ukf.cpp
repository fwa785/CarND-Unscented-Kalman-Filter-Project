#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  n_x_ = 5;

  n_aug_ = n_x_ + 2;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  x_ << 1, 1, 0, 0, 0;

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;


  is_initialized_ = false;

  // open file to write NIS
  NIS_laser_f.open("nis_laser.dat", std::ofstream::out);
  NIS_radar_f.open("nis_radar.dat", std::ofstream::out);

  if (!NIS_laser_f.is_open()) {
    cout << "Error Open nis_laser.dat file" << endl;
  }
  else {
    NIS_laser_f.close();
  }

  if (!NIS_radar_f.is_open()) {
    cout << "Error Open nis_radar.dat file" << endl;
  }
  else {
    NIS_radar_f.close();
  }
}

UKF::~UKF() {
  if (!NIS_laser_f.is_open()) {
    NIS_laser_f.close();
  }
  if (!NIS_radar_f.is_open()) {
    NIS_radar_f.close();
  }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && !use_laser_) {
    return;
  }

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && !use_radar_) {
    return;
  }

  if (!is_initialized_) {
    cout << "UKF:" << endl;
    // initialize x_
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      x_(0) = px;
      x_(1) = py;
      x_(2) = 1;
      x_(3) = atan2(py, px);
      x_(4) = 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
      x_(2) = rho_dot;
      x_(3) = phi;
      x_(4) = 0;
    }

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  /* prediction */
  Prediction(dt);

  /* update */
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_aug_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0.0;
  x_aug(n_x_ + 1) = 0.0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  *Xsig_aug_out = Xsig_aug;
}

void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t) {
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x = Xsig_aug.col(i);
    double px = x(0);
    double py = x(1);
    double v = x(2);
    double yaw = x(3);
    double yawd = x(4);
    double va = x(5);
    double vyawdd = x(6);

    double delta_t_2 = delta_t * delta_t;

    if (fabs(yawd) <= 0.001) {
      Xsig_pred_(0, i) = px + v * cos(yaw) * delta_t + cos(yaw)*(delta_t_2)*va / 2.;
      Xsig_pred_(1, i) = py + v * sin(yaw) * delta_t + sin(yaw)*(delta_t_2)*va / 2.;
      Xsig_pred_(2, i) = v + va * delta_t;
      Xsig_pred_(3, i) = yaw + yawd * delta_t + vyawdd * (delta_t_2) / 2.;
      Xsig_pred_(4, i) = yawd + vyawdd * delta_t;
    }
    else {
      Xsig_pred_(0, i) = px + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) + cos(yaw)*(delta_t_2)*va / 2.;
      Xsig_pred_(1, i) = py + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t)) + sin(yaw)*(delta_t_2)*va / 2.;
      Xsig_pred_(2, i) = v + va * delta_t;
      Xsig_pred_(3, i) = yaw + yawd * delta_t + vyawdd * (delta_t_2) / 2.;
      Xsig_pred_(4, i) = yawd + vyawdd * delta_t;
    }
  }
}

void UKF::PredictMeanCoverance() {
  //predicted state mean
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Augument the Sigma Points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  GenerateSigmaPoints(&Xsig_aug);

  // Predict the Sigma Points
  PredictSigmaPoints(Xsig_aug, delta_t);

  // Predict Mean and Covariance
  PredictMeanCoverance();
}

void UKF::PredictLaserMeasurement(MatrixXd *Zsig_out, VectorXd *z_pred_out, MatrixXd *S_out) {
  // Laser measurement has px and py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd xsig_pred = Xsig_pred_.col(i);
    double px = xsig_pred(0);
    double py = xsig_pred(1);
    double v = xsig_pred(2);
    double yaw = xsig_pred(3);
    double yaw_rate = xsig_pred(4);

    VectorXd zsig = VectorXd(n_z);
    zsig(0) = px;
    zsig(1) = py;
    Zsig.col(i) = zsig;
  }

  //mean predicted measurement
  z_pred = Zsig * weights_;

  //innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
    0, std_laspy_*std_laspy_;
  S = S + R;

  *Zsig_out = Zsig;
  *z_pred_out = z_pred;
  *S_out = S;
}

void UKF::UpdateLaserState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, MeasurementPackage meas_package) {
  int n_z = 2;

  VectorXd z = meas_package.raw_measurements_;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // predict the laser measurement
  PredictLaserMeasurement(&Zsig, &z_pred, &S);

  // Update the state
  UpdateLaserState(Zsig, z_pred, S, meas_package);

  // Calculate NIS
  VectorXd z = meas_package.raw_measurements_;
  double nis = (z - z_pred).transpose() * S.inverse() * (z - z_pred);
  NIS_laser_f.open("nis_laser.dat", std::ofstream::out | std::ofstream::app);
  if (NIS_laser_f.is_open()) {
    NIS_laser_f << nis << endl;
    NIS_laser_f.close();
  }
}

void UKF::PredictRadarMeasurement(MatrixXd *Zsig_out, VectorXd *z_pred_out, MatrixXd *S_out) {
  // Laser measurement has px and py
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd xsig_pred = Xsig_pred_.col(i);
    double px = xsig_pred(0);
    double py = xsig_pred(1);
    double v = xsig_pred(2);
    double yaw = xsig_pred(3);
    double yaw_rate = xsig_pred(4);

    VectorXd zsig = VectorXd(n_z);
    zsig(0) = sqrt(px * px + py * py);
    zsig(1) = atan2(py, px);
    zsig(2) = (px * cos(yaw) * v + py * sin(yaw) * v) / zsig(0);
    Zsig.col(i) = zsig;
 }

  //mean predicted measurement
  z_pred = Zsig * weights_;

  //innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
    0, std_radphi_* std_radphi_, 0,
    0, 0, std_radrd_ *std_radrd_;
  S = S + R;

  *Zsig_out = Zsig;
  *z_pred_out = z_pred;
  *S_out = S;
}

void UKF::UpdateRadarState(MatrixXd Zsig, VectorXd z_pred, MatrixXd S, MeasurementPackage meas_package) {
  int n_z = 3;

  VectorXd z = meas_package.raw_measurements_;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  // predict the laser measurement
  PredictRadarMeasurement(&Zsig, &z_pred, &S);

  // Update the state
  UpdateRadarState(Zsig, z_pred, S, meas_package);

  // Calculate NIS
  VectorXd z = meas_package.raw_measurements_;
  double nis = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  NIS_radar_f.open("nis_radar.dat", std::ofstream::out | std::ofstream::app);
  if (NIS_radar_f.is_open()) {
    NIS_radar_f << nis << endl;
    NIS_radar_f.close();
  }
}

