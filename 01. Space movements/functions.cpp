#include <iostream>
#include <vector>
#include </usr/local/include/eigen3/Eigen/Dense>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace Eigen;


MatrixXd Euler2A(double fi, double teta, double psi){
    MatrixXd M1(3, 3), M2(3, 3), M3(3, 3);
    M1 << 1.0, 0.0,       0.0,
          0.0, cos(fi), -sin(fi),
          0.0, sin(fi), cos(fi);

    M2 << cos(teta),  0.0, sin(teta),
          0.0      ,  1.0, 0.0,
          -sin(teta), 0.0, cos(teta);

    M3 << cos(psi), -sin(psi), 0.0,
          sin(psi), cos(psi),  0.0,
          0.0,        0.0,     1.0;

    return M3 * M2 * M1;
}

pair<VectorXd, double> AxisAngle(Matrix3d A){
    MatrixXd E(3, 3);
    E.setIdentity();

    float dA = round(A.determinant());  
    Matrix3d AAtransp = A * A.transpose();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            AAtransp(i, j) = round(abs(AAtransp(i, j)));
        }
    }

    if (dA != 1 || AAtransp != E) {
        cout << "Matrica je neispravna!" << endl;
        exit(1);
    }

    Matrix3d Ae;
    Ae = A - E;

    Vector3d array = Ae.row(0).cross(Ae.row(1));
    array.normalize();

    Vector3d v = {0, -(array(2) / array(1)), 1};
    Vector3d v1 = A * v;

    double fi = acos(v.dot(v1) / (v.norm() * v1.norm()));

    return pair<VectorXd, double>(array, fi);
}

MatrixXd Rodrigez(Vector3d array, double fi){
    MatrixXd E(3, 3);
    E.setIdentity();

    MatrixXd M(3, 3);
    M << 0, -array(2), array(1),
          array(2), 0, -array(0),
         -array(1), array(0), 0;

    MatrixXd Rod = (array * array.transpose()) + (cos(fi) * (E - array * array.transpose()) + sin(fi) * M);

    return Rod;
}

tuple<double, double, double> Auler2E(MatrixXd A){
    double x, y, z;

    if(A(2, 0) < 1){
        if(A(2, 0) > -1){
            z = atan2(A(1, 0), A(0, 0));
            y = asin(-A(2, 0));
            x = atan2(A(2, 1), A(2, 2));
        } else {
            x = atan2(-A(0, 1), A(1, 1));
            y = M_PI / 2;
            x = 0;
        }
    } else {
        x = atan2(-A(0, 1), A(1, 1));
        y = -M_PI / 2;
        x = 0;
    }

    return tuple<double, double, double> (x, y, z);
}

Vector4d AxisAngle2Q(Vector3d array, double fi){
    double o = cos(fi/2);

    array.normalize();
    Vector3d m = sin(fi/2) * array;

    Vector4d sol;
    sol << m, o;

    return sol;
}

pair<VectorXd, double> Q2AxisAngle(Vector4d v){
    v.normalize();
    if(v(3) < 0) v = -v;

    double fi = 2 * acos(v(3));

    Vector3d array, q3 = {v(0), v(1), v(2)};
    q3.normalize();

    if(abs(v(3)) == 1) array = {1, 0, 0};
    else array = q3;

    return pair<VectorXd, double>(array, fi);
}

int main() {
    cout << fixed;
    cout << setprecision(8);
    double fi = -atan(0.142857143);     // -arctg(1/7)
    double teta = -asin(0.333333333);          // -arcsin(1/3)
    double psi = atan(2);               // arctg(2)

    cout << "fi: " << fi << endl;
    cout << "teta: " << teta << endl;
    cout << "psi: " << psi << endl;
    cout << endl;

    MatrixXd A = Euler2A(fi, teta, psi);
    cout << "--- Euler2A ---" << endl;
    cout << A << endl;
    cout << endl;

    pair<VectorXd, double> pfi = AxisAngle(A);
    cout << "--- AxisAngle ---" << endl;
    cout << "[" << pfi.first << "]" << endl;
    cout << pfi.second << endl;
    cout << endl;

    MatrixXd Rod = Rodrigez(pfi.first, pfi.second);
    cout << "--- Rodrigez ---" << endl;
    cout << Rod << endl;
    cout << endl;

    tuple<double, double, double> auler2e = Auler2E(A);
    cout << "--- A2Euler ---" << endl;
    cout << "fi: " << get<0>(auler2e) << endl;
    cout << "teta: " << get<1>(auler2e) << endl;
    cout << "psi: " << get<2>(auler2e) << endl;
    cout << endl;

    Vector4d aa2q = AxisAngle2Q(pfi.first, pfi.second);
    cout << "--- AxisAngle2Q ---" << endl;
    cout << aa2q << endl;
    cout << endl;

    pair<VectorXd, double> q2fi = Q2AxisAngle(aa2q);
    cout << "--- Q2AxisAngle ---" << endl;
    cout << "array: " << endl << "[" << q2fi.first << "]" << endl;
    cout << "fi: " << q2fi.second << endl;
    cout << endl;

    return 0;
}
