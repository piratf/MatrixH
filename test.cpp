#include <iostream>
#include <complex>
#include <time.h>
#ifdef __WIN32
#include <windows.h>
#endif
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

#ifdef __WIN32
    Matrix<double> c(m.row(), n.row());
    Matrix<double> d(m.row(), n.row());
    LARGE_INTEGER freq;
    LARGE_INTEGER start_t, stop_t;
    LARGE_INTEGER freq1;
    LARGE_INTEGER start_t1, stop_t1;
    double exe_time;
    QueryPerformanceFrequency(&freq);
    m.inv();
    // fprintf(stdout, "The frequency of your pc is %d.\n", freq.QuadPart);
    cerr << "The frequency of your pc is  " << freq.QuadPart << endl;
    QueryPerformanceCounter(&start_t);
    QueryPerformanceCounter(&stop_t);
    exe_time = 1e3 * (stop_t.QuadPart - start_t.QuadPart) / freq.QuadPart;
    cerr << "inv use Time: " << exe_time << endl;
    // fprintf(stdout, "Your program executed time is %fms.\n", exe_time);


    QueryPerformanceFrequency(&freq1);
    // fprintf(stdout, "The frequency of your pc is %d.\n", freq.QuadPart);
    cerr << "The frequency of your pc is  " << freq1.QuadPart << endl;
    QueryPerformanceCounter(&start_t1);
    m.smul(d, n);
    QueryPerformanceCounter(&stop_t1);
    exe_time = 1e3 * (stop_t1.QuadPart - start_t1.QuadPart) / freq1.QuadPart;
    cerr << "Use Time: " << exe_time << endl;
    // fprintf(stdout, "Your program executed time is %fms.\n", exe_time);
#endif


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