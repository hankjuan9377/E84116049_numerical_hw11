#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double PI = acos(-1.0);
const int N = 3;         // basis functions phi_1, phi_2, phi_3
const int M = 10000;      // integration steps

// Right-hand side function R(x)
double R(double x) {
    return (1 - x * x) * exp(-x) + (x + 1);
}

// Basis function phi_i(x) = sin((i+1) * pi * x)
double phi(int i, double x) {
    return sin((i + 1) * PI * x);
}

// phi_i'(x)
double dphi(int i, double x) {
    return (i + 1) * PI * cos((i + 1) * PI * x);
}

// phi_i''(x)
double ddphi(int i, double x) {
    return -pow((i + 1) * PI, 2) * sin((i + 1) * PI * x);
}

int main() {
    double A[N][N] = {0};  // A*c=b 
    double b[N] = {0};     
    double c[N];           
    int i, j, m;

    double h = 1.0 / M;
    for (m = 0; m <= M; ++m) {
        double x = m * h;
        double w = (m == 0 || m == M) ? 0.5 : 1.0;//±è§Îªk 

        for (j = 0; j < N; ++j) {
            for (i = 0; i < N; ++i) {
                double term = ddphi(i, x) * phi(j, x)
                            + (x + 1) * dphi(i, x) * phi(j, x)
                            - 2 * phi(i, x) * phi(j, x);
                A[j][i] += w * term * h;
            }
            b[j] += w * R(x) * phi(j, x) * h;
        }
    }

    // Gaussian
    for (i = 0; i < N; ++i) {
        double pivot = A[i][i];
        for (j = i; j < N; ++j) A[i][j] /= pivot;
        b[i] /= pivot;

        for (int k = i + 1; k < N; ++k) {
            double factor = A[k][i];
            for (j = i; j < N; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    for (i = N - 1; i >= 0; --i) {
        c[i] = b[i];
        for (j = i + 1; j < N; ++j)
            c[i] -= A[i][j] * c[j];
    }

    // Output coefficients
    cout << fixed << setprecision(6);
    cout << "(n = 3):" << endl;
    for (i = 0; i < N; ++i) {
        cout << "c[" << i + 1 << "] = " << c[i] << endl;
    }

    // Print approximate solution at 11 points
    cout << "\nApproximate y(x) = 1 + x + sum(c_i * sin(i*pi*x)):" << endl;
    for (m = 0; m <= 10; ++m) {
        double x = m / 10.0;
        double y = 1 + x;
        for (i = 0; i < N; ++i)
            y += c[i] * phi(i, x);
        cout << "x = " << setw(4) << x << ",\ty(x) = " << y << endl;
    }

    return 0;
}

