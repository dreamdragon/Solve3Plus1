#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unsupported/Eigen/Polynomials>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

#include "Solve3Plus1.h"

using namespace Eigen;

typedef Matrix<double,5,1> Vector5d; 
typedef Array<double,5,1> Array5d; 
typedef Matrix<double, Dynamic, 1, 0, 4, 1> Vector4di;
typedef Matrix<double, 4, Dynamic, 0, 4, 4> Matrix4di;

inline Matrix3d skewMatrix(const Vector3d& vec)
{
    return (Matrix3d()  << 
            0, -vec(2,0), vec(1,0),
            vec(2,0), 0, -vec(0,0),
            -vec(1,0), vec(0,0), 0).finished();
}
inline Matrix3d skewMatrix(double x, double y, double z)
{
    return (Matrix3d()  << 
            0, -z, y,
            z, 0, -x,
            -y, x, 0).finished();
}

// input, 2 3x3 matrix
MatrixXd SolveConstrained3ptClosed(Matrix3d& q, Matrix3d& qp)
{
    double a11 = -1-qp(0,0)*q(0,0);
    double a12 = qp(1,0)*q(0,0);
    double a13 = -q(0,0)+qp(0,0);
    double a14 = -qp(1,0);
    double a15 = q(1,0);
    double a16 = -qp(0,0)*q(1,0);

    double a21 = -1-qp(0,1)*q(0,1);
    double a22 = qp(1,1)*q(0,1);
    double a23 = -q(0,1)+qp(0,1);
    double a24 = -qp(1,1);
    double a25 = q(1,1);
    double a26 = -qp(0,1)*q(1,1);

    double a31 = -1-qp(0,2)*q(0,2);
    double a32 = qp(1,2)*q(0,2);
    double a33 = -q(0,2)+qp(0,2);
    double a34 = -qp(1,2);
    double a35 = q(1,2);
    double a36 = -qp(0,2)*q(1,2);

    double g1 = a31*a14*a22+a34*a12*a21-a34*a22*a11-a32*a14*a21-a31*a12*a24+a32*a24*a11;
    double g2 = a33*a14*a22+a34*a12*a23+a32*a24*a13-a34*a22*a13-a33*a12*a24-a32*a14*a23;
    double g3 = -a34*a11*a25+a31*a14*a25+a34*a15*a21-a31*a24*a15+a36*a22*a11+a35*a24*a11-a36*a12*a21-a31*a22*a16+a31*a12*a26-a35*a14*a21-a32*a11*a26+a32*a16*a21;
    double g4 = -a33*a22*a16+a34*a15*a23-a34*a11*a26+a33*a12*a26+a34*a16*a21+a31*a14*a26-a31*a24*a16-a32*a13*a26-a32*a15*a21+a32*a16*a23+a32*a11*a25+a31*a22*a15-a36*a14*a21-a36*a12*a23+a36*a22*a13+a36*a24*a11-a34*a13*a25-a35*a14*a23+a35*a12*a21-a35*a22*a11+a35*a24*a13+a33*a14*a25-a33*a24*a15-a31*a12*a25;
    double g5 = -a36*a14*a23+a36*a24*a13-a34*a13*a26+a35*a12*a23-a33*a24*a16-a35*a22*a13+a33*a22*a15+a33*a14*a26-a32*a15*a23-a33*a12*a25+a34*a16*a23+a32*a13*a25;
    double g6 = a36*a11*a25-a35*a11*a26-a36*a21*a15+a35*a16*a21-a31*a16*a25+a31*a15*a26;
    double g7 = -a36*a23*a15+a33*a15*a26+a35*a16*a23+a36*a13*a25-a33*a16*a25-a35*a13*a26;

    // the order in reverse in matlab and eigen
    Vector5d H;
    H << (-2*g2*g7-g7*g7+g5*g5-g2*g2),
      (-2*g7*g4+2*g1*g5+2*g6*g5-2*g2*g4),
      (2*g2*g7+2*g3*g5+g1*g1+g2*g2+g7*g7+g6*g6+2*g6*g1-2*g5*g5-g4*g4),
      (2*g7*g4-2*g6*g5+2*g2*g4+2*g6*g3-2*g1*g5+2*g1*g3), 
      (-2*g3*g5+g3*g3+g5*g5+g4*g4);

    // solve the 4th order poly:
    PolynomialSolver<double,4> psolve( H );

    //y = roots(H);

    // kill off imaginary roots:
    std::vector<double> realRoots;
    psolve.realRoots( realRoots );
    // reverse the order of y, can be disabled---------------------------------------
    std::reverse(realRoots.begin(),realRoots.end());
    // ------------------------------------------------------------------------------
    Map<ArrayXd> y( &realRoots[0], realRoots.size() );
    std::cout << "Real roots: " << y.transpose() << std::endl;

    // or we can use Array directly
    ArrayXd z = -(g6*y+g1*y+g3*y*y-g5*y*y+g5)/(g2+g7+g4*y);
    ArrayXd d = (a22*a11-a12*a21)*y*y+((-a14*a21+a22*a13+a24*a11-a12*a23)*z+a11*a25-a15*a21)*y+(a24*a13-a14*a23)*z*z+(a13*a25-a15*a23)*z;
    // kill off zero denominators:
    //nz = d~=0 & (g2+g7+g4*y) ~= 0;
    //y = y( nz );
    //z = z( nz );
    //d = d( nz );

    ArrayXd t = g2+g7+g4*y;
    std::vector<double> vy, vz, vd;
    vy.reserve(t.size()*sizeof(double));
    vz.reserve(t.size()*sizeof(double));
    vd.reserve(t.size()*sizeof(double));
    for ( int i = 0; i < t.size(); ++i )
    {
        if ((d(i) != 0) && (t(i) !=0))
        {
            vy.push_back(y(i));
            vz.push_back(z(i));
            vd.push_back(d(i));
        }
    }
    ArrayXd dy( vy.size() ); 
    ArrayXd dz( vz.size() ); 
    ArrayXd dd( vd.size() ); 
    for ( int i = 0; i < vy.size(); ++i )
    {
        dy(i) = vy[i]; 
        dz(i) = vz[i]; 
        dd(i) = vd[i]; 
    }

    std::cout << "dy: " << dy.transpose() << std::endl;
    std::cout << "dz: " << dz.transpose() << std::endl;
    std::cout << "dd: " << dd.transpose() << std::endl;

    ArrayXd dw = -((-a14*a22+a12*a24)*dy*dy+(a15*a24+a16*a22-a14*a25-a12*a26)*dy+(-a14*a22+a12*a24)*dz*dz+(-a14*a26+a12*a25-a15*a22+a16*a24)*dz-a15*a26+a16*a25)/d;
    ArrayXd dx = -((a21*a14-a11*a24)*dy*dy+((-a21*a12+a23*a14+a11*a22-a13*a24)*dz+a11*a26-a21*a16)*y+(a13*a22-a23*a12)*z*z+(a13*a26-a23*a16)*z)/d;

    std::cout << "dw: " << dw.transpose() << std::endl;
    std::cout << "dx: " << dx.transpose() << std::endl;

    //E = [dw';dx';dy';dz'];
    Matrix4di E (vy.size(), 4);
    E.row(0) = dw.transpose();
    E.row(1) = dx.transpose();
    E.row(2) = dy.transpose();
    E.row(3) = dz.transpose();
    return E;    
}

MatrixXd Solve3Plus1(Matrix3d& Q, Matrix3d& Qp)
{
    Matrix4d R = SolveConstrained3ptClosed(Q,Qp);
    MatrixXd E;

    if (R.size() > 0)
    {
        for ( int i = 0; i < R.rows(); ++i )
        {
            R(3,i) = atan2(R(2,i),R(3,i));
        }

        std::cout << "R: \n" << R << std::endl;

        MatrixXd S = MatrixXd::Ones(5,R.cols());
        S.row(3) = R.row(2);
        S.row(4) = R.row(3);

        S.row(0) = R.row(1);
        S.row(1) = R.row(0);

        MatrixXd subS = S.topRows(3).cwiseProduct(S.topRows(3));
        VectorXd N = (subS.row(0)+subS.row(1)+subS.row(2)).cwiseSqrt();
        S.topRows(3) = S.topRows(3).cwiseProduct((N.transpose().replicate(3,1).cwiseInverse()));

        std::cout << "S: \n" << S << std::endl;

        // construct essential matrix
        int n = 2;
        E = MatrixXd::Zero(9,S.cols()*n);
        Matrix3d u;
        u <<  0, 0, 1, 
              0, 0, 0,
              -1, 0, 0;
        MatrixXd Tx1, Tx2;
        MatrixXd E1,E2,R;
        for ( int i = 0; i < S.cols(); ++i )
        {
            Tx1 = skewMatrix(S(0,i),S(1,i),S(2,i));
            Tx2 = skewMatrix(-S(0,i),-S(1,i),-S(2,i));
            //std::cout << "sk: \n" << sk << std::endl;
            R = (u*S(4,i)).exp();
            E1 = Tx1*R;
            E2 = Tx2*R;
            //std::cout << "E1: \n" << E1 << std::endl;
            //std::cout << "E2: \n" << E2 << std::endl;
            E1.resize(9,1);
            E2.resize(9,1);
            E.col(i*n+0) = E1;
            E.col(i*n+1) = E2;
            //E.col(i*n+2) = -t1;
            //E.col(i*n+3) = -t2;
        }
    }
    
    return E;    
}

MatrixXd Solve3Plus1EMatrix(Matrix3d& K, Matrix3d& Kp, Vector3d& g, Vector3d& gp)
{
    Vector3d c;
    c << 0,1,0;

    Matrix3d R, Q;
    Vector3d a = g.cross(c);
    if (a.norm() > 0.000001)
    {
        a = a/a.norm();
        double theta = acos(g.dot(c));
        R = (skewMatrix(a)*theta).exp();
        Q = R*K;
        for ( int i = 0; i < 3; ++i )
        {
            Q.row(i) = Q.row(i).array()/Q.row(2).array();
        }
    }
    else
    {
        R = Matrix3d::Identity(3,3);
        Q = K;
    }

    Matrix3d Rp, Qp;
    Vector3d ap = gp.cross(c);
    if (ap.norm() > 0.000001)
    {
        ap = ap/ap.norm();
        double thetap = acos(gp.dot(c));
        Rp = (skewMatrix(ap)*thetap).exp();
        Qp = Rp*Kp;
        for ( int i = 0; i < 3; ++i )
        {
            Qp.row(i) = Qp.row(i).array()/Qp.row(2).array();
        }
    }
    else
    {
        Rp = Matrix3d::Identity(3,3);
        Qp = Kp;
    }

    return Solve3Plus1(Q,Qp);
}



Matrix3d q,qp;
Matrix4d E;

void testInit()
{
    q << -0.1964,   -0.3617,   -0.1811,
      0.5821,    0.3890 ,   0.8287,
      1.0000,    1.0000,    1.0000; 

    qp <<   -0.1650,   -0.3438,   -0.1576,
       0.6253,    0.4380,    0.8586,
       1.0000,    1.0000,    1.0000;

    std::cout << "q:\n" << q << std::endl;
    std::cout << "qp:\n" << qp << std::endl;

    E <<   1.4593,    1.1779,   -5.0968,    0.1040,
      -4.4852,    1.4482,   -0.0304,   -0.3527,
      0.2896,   -0.0653,   -0.0580,    0.0147,
      0.9572,    0.9979,   -0.9983,    0.9999;
}


void testSolveConstrained3ptClosed()
{
    testInit();
    MatrixXd result = SolveConstrained3ptClosed(q,qp);

    std::cout << "ground trouth:\n" << E << std::endl;
    std::cout << "result:\n" << result << std::endl;
}

void testSolve3Plus1()
{
    testInit();
    MatrixXd result = Solve3Plus1(q,qp);
    std::cout << "result:\n" << result << std::endl;
}

std::string testSolve3Plus1EMatrix()
{
    Matrix3d K, Kp;
    Vector3d g, gp;
    K << 
    0.0253,   -0.1879,    0.1826,
   -0.1937,   -0.1300,   -0.1317,
    1.0000,    1.0000,    1.0000;
    
    Kp <<
        0.4116,    0.1779,    0.6017,
    0.3273,    0.4019,    0.4124,
    1.0000,    1.0000,    1.0000;

    g << 0.7984,
        0.2711,
        0.5377;

    gp << 0.8919,
        0.4450,
        0.0804;

    MatrixXd E = Solve3Plus1EMatrix(K,Kp,g,gp);
    std::stringstream ss;
    ss << "Essential matrices:\n " << E;
    std::cout << ss.str() << std::endl;
    return ss.str();
}

int main()
{
    testSolve3Plus1EMatrix();
    //testSolve3Plus1();
    //testSolveConstrained3ptClosed();
}
