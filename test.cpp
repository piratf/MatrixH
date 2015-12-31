#include <iostream>
#include <complex>
#include <time.h>
#include "Matrix.h"
using namespace std;

// 矩阵大小
const unsigned D = 10;

int main() {
    ios::sync_with_stdio(false);
    freopen("test.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    srand((unsigned)time(NULL));

    Matrix<double> m(D, D), n(D, D);
    for (unsigned i = 0; i < D; ++i) {
        for (unsigned j = 0; j < D; ++j) {
            m[i][j] = rand() % 10;
            // cin >> m[i][j];
        }
    }
    for (unsigned i = 0; i < D; ++i) {
        for (unsigned j = 0; j < D; ++j) {
            n[i][j] = rand() % 10;
            // cin >> n[i][j];
        }
    }
    cout << Matrix<double>::eye(D, D) << endl;
    cout << m << endl;
    cout << n << endl;
    cout << m + n << endl;
    cout << m - n << endl;
    cout << m * n << endl;
    cout << m / n << endl;
    cout << m % n << endl;
    cout << m.cut(0, 1, 4, 7);
    cout << n / m << endl;
    auto bitch = n.eig();
    cout << bitch << endl;
    cout << (m ^ 100) << endl;
    cout << m.avg() << endl;
    cout << m.eig() << endl;
    m.cov();
    m.inv();
    m.hess();
    cout << m.cond2() << endl;
    return 0;
}