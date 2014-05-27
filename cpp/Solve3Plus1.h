#include <Eigen/Core>
#include <string>

inline Eigen::Matrix3d skewMatrix(const Eigen::Vector3d& vec);
inline Eigen::Matrix3d skewMatrix(double x, double y, double z);
Eigen::MatrixXd SolveConstrained3ptClosed(Eigen::Matrix3d& q, Eigen::Matrix3d& qp);
Eigen::MatrixXd Solve3Plus1(Eigen::Matrix3d& Q, Eigen::Matrix3d& Qp);
Eigen::MatrixXd Solve3Plus1EMatrix(Eigen::Matrix3d& K, Eigen::Matrix3d& Kp, 
        Eigen::Vector3d& g, Eigen::Vector3d& gp);
void testInit();
void testSolveConstrained3ptClosed();
void testSolve3Plus1();
std::string testSolve3Plus1EMatrix();
