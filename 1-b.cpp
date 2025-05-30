#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

const int N = 10;             // 區間數
const double h = 0.1;         // 步長
const double Y_LEFT = 1.0;    // y(0)
const double Y_RIGHT = 2.0;   // y(1)


double rhs(double x) {
    return (1 - x * x) * exp(-x);
}

int main() {
    vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = i * h;

    int M = N - 1; // y1~y9
    vector<vector<double> > A(M);
    for (int i = 0; i < M; ++i)
        A[i].resize(M, 0.0); // 初始化為 0.0

    vector<double> d(M), y(M);

    // 建立 A 矩陣與 d 向量
    for (int i = 0; i < M; ++i) {
        double xi = x[i + 1]; // 對應 x1 ~ x9

        double ai = 1.0 - (xi + 1) * h / 2.0;      // A_i
        double bi = -2.0 - 2.0 * h * h;            // B_i
        double ci = 1.0 + (xi + 1) * h / 2.0;      // C_i

        d[i] = h * h * rhs(xi); // D_i

        if (i == 0) {
            A[i][i] = bi;
            A[i][i + 1] = ci;
            d[i] -= ai * Y_LEFT;
        } else if (i == M - 1) {
            A[i][i - 1] = ai;
            A[i][i] = bi;
            d[i] -= ci * Y_RIGHT;
        } else {
            A[i][i - 1] = ai;
            A[i][i] = bi;
            A[i][i + 1] = ci;
        }
    }

    // Gaussian elimination - 前向消去
    for (int i = 0; i < M - 1; ++i) {
        double factor = A[i + 1][i] / A[i][i];
        for (int j = i; j < M; ++j)
            A[i + 1][j] -= factor * A[i][j];
        d[i + 1] -= factor * d[i];
    }

    // 回代求解
    y[M - 1] = d[M - 1] / A[M - 1][M - 1];
    for (int i = M - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < M; ++j)
            sum += A[i][j] * y[j];
        y[i] = (d[i] - sum) / A[i][i];
    }

    // 輸出結果
    cout << fixed << setprecision(6);
    cout << "x\t\ty(x)" << endl;
    cout << x[0] << "\t" << Y_LEFT << endl;
    for (int i = 1; i < N; ++i)
        cout << x[i] << "\t" << y[i - 1] << endl;
    cout << x[N] << "\t" << Y_RIGHT << endl;

    return 0;
}

