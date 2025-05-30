#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

const double h = 0.1;
const int N = 10;


// y'' = -(x+1)y' + 2y + (1 - x^2)e^{-x}
double f_1(double x, double y, double dy) {
    return -(x + 1) * dy + 2 * y + (1 - x * x) * exp(-x);
}

// y'' = -(x+1)y' + 2y
double f_2(double x, double y, double dy) {
    return -(x + 1) * dy + 2 * y;
}

// Runge-Kutta method for y1 or y2
void RK4(vector<double>& y_out, vector<double>& x_out,
         double y0, double dy0,
         double (*f)(double, double, double)) {
    double x = 0.0;		//x0 
    double y = y0;		//y0
    double z = dy0;		//dy0

    y_out.clear();
    x_out.clear();

    for (int i = 0; i <= N; ++i) {
        y_out.push_back(y);
        x_out.push_back(x);

        double k1 = h * z;
        double l1 = h * f(x, y, z);

        double k2 = h * (z + l1 / 2);
        double l2 = h * f(x + h / 2, y + k1 / 2, z + l1 / 2);

        double k3 = h * (z + l2 / 2);
        double l3 = h * f(x + h / 2, y + k2 / 2, z + l2 / 2);

        double k4 = h * (z + l3);
        double l4 = h * f(x + h, y + k3, z + l3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        z += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        x += h;
    }
}

int main() {
    vector<double> y1, y2, x;

    // Step 1: Solve y1'' = f_nonhom, y1(0)=1, y1'(0)=0
    RK4(y1, x, 1.0, 0.0, f_1);

    // Step 2: Solve y2'' = f_hom, y2(0)=0, y2'(0)=1
    RK4(y2, x, 0.0, 1.0, f_2);

    // Step 3: Compute c = (2 - y1(1)) / y2(1)
    double c = (2.0 - y1.back()) / y2.back();

    cout << fixed << setprecision(6);
    cout << "c = " << c << endl;
    cout << "y(x) = y1(x) + c * y2(x):\n";
    cout << "x\t\ty(x)\n";

    for (int i = 0; i <= N; ++i) {
        double y = y1[i] + c * y2[i];
        cout << x[i] << "\t" << y << endl;
    }

    return 0;
}

