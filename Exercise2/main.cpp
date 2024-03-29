#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


Vector2d sol_palu(Matrix2d& A,Vector2d& b)
{
    Vector2d x = A.fullPivLu().solve(b);
    return x;
}

Vector2d sol_qr(Matrix2d& A,Vector2d& b)
{
    Vector2d x = A.fullPivHouseholderQr().solve(b);
    return x;
}

double err(Vector2d& x,Vector2d& exact)
{
    double err= (x-exact).norm()/exact.norm();
    return err;
}

void resolution(Matrix2d& A,Vector2d& b, Vector2d& exact)
{
    cout << "A= \n" << A <<"\n"<<"b= \n" << b << endl;
    Vector2d x_palu = sol_palu(A,b) ;
    double err_palu = err(x_palu ,exact) ;
    cout << "L'errore relativo commesso risolvendo Ax=b con la fattorizzazione PALU e' err= " << err_palu << endl;
    Vector2d x_qr = sol_qr(A,b) ;
    double err_qr = err(x_qr ,exact) ;
    cout << "L'errore relativo commesso risolvendo Ax=b con la fattorizzazione QR e' err= " << err_qr << "\n" <<endl;
}

int main()
{
    Vector2d exactSol(-1.0,-1.0) ;

    // Case 1
    Matrix2d A;
    A << 5.547001962252291e-01, -3.770900990025203e-02, 8.320502943378437e-01, 9.992887623566787e-01;
    Vector2d b1(-5.169911863249772e-01, 1.672384680188350e-01);
    resolution(A,b1,exactSol);

    // Case 2
    Matrix2d B;
    B << 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
    Vector2d b2(-6.394645785530173e-04, 4.259549612877223e-04);
    resolution(B,b2,exactSol);

    // Case 3
    Matrix2d C;
    C << 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01, -8.320502947645361e-01;
    Vector2d b3(-6.400391328043042e-10, 4.266924591433963e-10);
    resolution(C,b3,exactSol);

  return 0;
}
