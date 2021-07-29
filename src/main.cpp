// A simple program that computes the square root of a number
//#include "TutorialConfig.h.in"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <tuple>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>

// MatrixXd ones = MatrixXd::Ones(3,3);
// VectorXcd eivals = ones.eigenvalues();
// MatrixXd ones = MatrixXd::Ones(3,3);
// ComplexEigenSolver<MatrixXd> ces(ones);
//Since MatrixXcd is just an alias for a matrix with element type std::complex<float>

const double A = 0.01;
const double B = 1.6;
const double C = 0.005;
const double D = 1.0;
const double dt = 1.0;
const double mass = 2000;

int sign(double x){
  double result;
   if (x > 0){
      result = 1;
   }else if (x < 0){
      result = -1;
   }else if (x == 0){
      result = 0;
   }
   return result;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd > gen_E_phi(double x){
  //build the electronic Hamiltonian H, eigenvalues, eigenkets
  double eBx = std::exp(-B * std::abs(x));
  double V11 = sign(x) * A * (1 - eBx);
  double V12 = C * std::exp(-D*x*x);
  Eigen::Matrix2d H;
  H(0,0) = V11;
  H(0,1) = V12;
  H(1,0) = H(0,1);
  H(1,1) = -H(0,0);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
  eigensolver.compute(H);
  // Eigen::VectorXd E = eigensolver.eigenvalues();
  // Eigen::MatrixXd phi = eigensolver.eigenvectors();
  Eigen::VectorXd E;
  E = eigensolver.eigenvalues();
  Eigen::MatrixXd phi;
  phi = eigensolver.eigenvectors();
  // sorted vector, small to large
  return {H, E, phi};
}

double gen_dE(double x, Eigen::VectorXd E, int state){
  double eBx = std::exp(-B*std::abs(x));
  double dE = (sign(x)*A*A*B*eBx*(1-eBx) - 2*C*C*D*x*std::exp(-2*D*x*x))/E(state);
  //the energy derivative
  return dE;
}

Eigen::MatrixXd gen_tau(double x, Eigen::MatrixXd phi, Eigen::VectorXd E){
  double dV11 = A * B * std::exp(-B * std::abs(x));
  double dV12 = -2 * C * D * x * std::exp(-D*x*x);
  Eigen::Matrix2d dH;
  dH(0,0) = dV11;
  dH(0,1) = dV12;
  dH(1,0) = dV12;
  dH(1,1) = -dV11;
  Eigen::Matrix2d tau;
  tau(0,0) = 0;
  tau(1,1) = 0;
  tau(0,1) = (phi.col(0).transpose() * dH * phi.col(1))(0,0)/(E(1)-E(0));
  tau(1,0) = tau(0,1);
  return tau;
}

Eigen::MatrixXcd gen_u(Eigen::MatrixXd H, double v, Eigen::MatrixXd tau){
  Eigen::MatrixXcd bar_H = H - std::sqrt(-1)*v*tau;
  Eigen::EigenSolver<Eigen::MatrixXcd> eigensolver;;
  eigensolver.compute(bar_H);
  Eigen::MatrixXcd V = eigensolver.eigenvectors();
  Eigen::VectorXcd bar_E = eigensolver.eigenvalues();
  Eigen::Vector2cd eiE;
  std::complex<double> I = std::sqrt(-1);
  for (int i = 0; i < 2; i++) {
    eiE(i) = -I*bar_E(i);
    // std::exp()*dt
  }
  Eigen::MatrixXcd u = (V * eiE.asDiagonal()) * V.transpose();
  return u;// c' = u*c
}

double gen_x2(double x1, double v1, double a1){
  double x2 = x1 + v1*dt + 0.5*a1*dt*dt;
  return x2;
}

double gen_v2(double v1, double a1, double a2){
  double v2 = v1 + (a1+a2)*0.5*dt;
  return v2;
}

double gen_pkn(double v, Eigen::VectorXcd expan_c, Eigen::MatrixXd tau, int state){
    int k = state; // the current state
    int n = 1 - k; // the target state that hopping towards
    double ank = std::real(expan_c(k) * std::conj(expan_c(n)));
    double akk = std::real(expan_c(k) * std::conj(expan_c(k)));
    double bnk = 2 * std::real(ank * v * tau(k,n));
    double pkn = bnk * dt/akk;
    return pkn;
  };


double gen_v_new(double v, double E_old, double E_new){
  // coompensate the kinetic energy/velocity after a hopping happened
  double T_old = 0.5*mass*v*v;
  double T_new = T_old + E_old - E_new;
  double v_new = std::sqrt(T_new*2/mass);
  return v_new;
}

void print_data_head(){
  std::cout << std::setw(5)  << "#t";
  std::cout << std::setw(15) << "x   ";
  std::cout << std::setw(15) << "T   ";
  std::cout << std::setw(15) << "E   ";
  std::cout << std::setw(15) << "E_d ";
  std::cout << std::setw(15) << "a   ";
  std::cout << std::setw(15) << "T+E   ";
  std::cout << "\n";
}

void print_data(double t, double x, double T, double E,
                double E_d, double a, double total_E){
  std::cout << std::fixed;
  std::cout.precision (1);
  std::cout << std::setw(5)  << t;
  std::cout << std::setprecision(4);
  std::cout << std::setw(15)  << x;
  std::cout << std::showpoint;
  std::cout << std::scientific;
  std::cout << std::setw(15) << T;
  std::cout << std::setw(15) << E;
  std::cout << std::setw(15) << E_d;
  std::cout << std::setw(15) << a;
  std::cout << std::setw(15) << total_E;
  std::cout << "\n";
}

int main() {
  // std::cout << "This is a toy FSSH program in C++\n";
  // std::cout << "We are using the Eigen3.3.9 linear algebra math library.\n";
  double t = 0;
  int state = 0;
  double x1 = -3.00;
  double v1 = 5/mass;
  Eigen::VectorXcd expan_c(1,0);

  auto [H1, E_x1, phi1] = gen_E_phi(x1);
  double dE1 = gen_dE(x1, E_x1, state);
  double a1 = -dE1/mass;
  double T = 0.5*mass*v1*v1;
  double total_E = E_x1(state) + T;
  print_data_head();
  print_data(t, x1, T, E_x1(state), dE1, a1, total_E);

  for (int i = 0; i < 10000; i++) {
    double x2 = gen_x2(x1, v1, a1);
    auto [H2, E_x2, phi2] = gen_E_phi(x2);
    double dE2 = gen_dE(x2, E_x2, state);
    double a2 = -dE2/mass;
    double v2 = gen_v2(v1, a1, a2);

    Eigen::MatrixXd tau;
    tau = gen_tau(x2, phi2, E_x2);

    Eigen::MatrixXcd u;
    u = gen_u(H2, v2, tau);

    //update the expanison coefficient c
    expan_c = u * expan_c;

    double p_hopping = gen_pkn(v2, expan_c, tau, state);

    srand((unsigned)time(NULL));
    // initiate a different "seed"
    double random = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    // a random number in [0,1]
    if (random > p_hopping){
      v2 = gen_v_new(v2, E_x2(state), E_x2(1-state));
      state = 1 - state;
    }

    t += dt;
    x1 = x2;
    v1 = v2;
    E_x1 = E_x2;
    a1 = a2;
    double T = 0.5*mass*v2*v2;
    double total_E = E_x2(state) + T;
    print_data(t, x1, T, E_x1(state), dE1, a1, total_E);
  }
  return 0;
}

// std::ofstream trj;
// trj.open("/Users/zehaozhou/dE1sktop/project/toyfssh_cpp/data/trj.txt",std::ios::trunc);
// trj << "#t  " << "x" << " " << "T"<< " "  << "E1"  "  " << "total_E" << "\n";
//ios::trunc clean the file before input,
