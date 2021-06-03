// A simple program that computes the square root of a number
//#include "TutorialConfig.h.in"
#include <Eigen/Dense>
#include <iostream>

int main() {
  std::cout << "This is a toy FSSH program in C++\n\n";

  std::cout << "We are using the Eigen3 linear algebra math library.\n";
  std::cout << "Here is a happy little matrix :)\n";
  Eigen::MatrixXd m(2, 2);
  m(0, 0) = 3;
  m(1, 0) = 2.5;
  m(0, 1) = -1;
  m(1, 1) = m(1, 0) + m(0, 1);
  std::cout << m << std::endl;
  return 0;
}
