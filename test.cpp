#include <iostream>
#include <complex>
#include <time.h>
#include "Matrix.h"
using namespace std;

// 矩阵大小
const unsigned D = 100;

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
    cerr << "print" << endl;
    cout << m << endl;
    cout << n << endl;
    cerr << "swap" << endl;
    m.swap(n);
    cerr << "print" << endl;
    cout << m << endl;
    cout << n << endl;
    cerr << "get contribute" << endl;
    cout << m.row() << ' ' << m.col() << endl;
    cout << m.max() << ' ' << m.min() << ' ' << m.avg() << endl;
    cerr << "operator +" << endl;
    cout << m + n << endl;
    cerr << "operator -" << endl;
    cout << m - n << endl;
    cerr << "operator *" << endl;
    cout << m * n << endl;
    cerr << "operator /" << endl;
    cout << m / n << endl;
    cerr << "operator %" << endl;
    cout << m % n << endl;
    cerr << "cut" << endl;
    cout << m.cut(0, m.row() / 2, 0, m.col() / 2);
    cerr << "eig" << endl;
    cout << n.eig() << endl;
    cerr << "operator ^" << endl;
    cout << (m ^ 100) << endl;
    cerr << "avg" << endl;
    cout << m.avg() << endl;
    cerr << "copy" << endl;
    n = m;
    cerr << "cov" << endl;
    cout << m.cov() << endl;
    cerr << "inv" << endl;
    cout << m.inv() << endl;
    cerr << "hess" << endl;
    cout << m.hess() << endl;
    cerr << "cond2" << endl;
    cout << m.cond2() << endl;
    cerr << "eye" << endl;
    cout << Matrix<double>::eye(D, D) << endl;
    return 0;
}