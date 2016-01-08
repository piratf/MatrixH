#include <iostream>
#include <complex>
#include <map>
#include <time.h>
#ifdef __WIN32
#include <windows.h>
#endif
#include "Matrix.h"
using namespace std;

// 矩阵大小
// 超过 100 以后在普通电脑上求逆和除法会比较慢
const unsigned D = 30;

void AccurateTest() {

    Matrix<double> m = Matrix<double>::rand(D, D);
    Matrix<double> n = Matrix<double>::rand(D, D);

#ifdef __WIN32
    LARGE_INTEGER freq;
    LARGE_INTEGER start_t, stop_t;
    LARGE_INTEGER freq1;
    LARGE_INTEGER start_t1, stop_t1;
    double exe_time;
    QueryPerformanceFrequency(&freq);
    // fprintf(stdout, "The frequency of your pc is %d.\n", freq.QuadPart);
    cout << "The frequency of your pc is  " << freq.QuadPart << endl;
    QueryPerformanceCounter(&start_t);
    m.mul(m, n);
    QueryPerformanceCounter(&stop_t);
    exe_time = 1e3 * (stop_t.QuadPart - start_t.QuadPart) / freq.QuadPart;
    cout << "use Time of mul: " << exe_time << "ms" << endl;
    // fprintf(stdout, "Your program executed time is %fms.\n", exe_time);


    QueryPerformanceFrequency(&freq1);
    // fprintf(stdout, "The frequency of your pc is %d.\n", freq.QuadPart);
    cout << "The frequency of your pc is  " << freq1.QuadPart << endl;
    QueryPerformanceCounter(&start_t1);
    m.smul(m, n);
    QueryPerformanceCounter(&stop_t1);
    exe_time = 1e3 * (stop_t1.QuadPart - start_t1.QuadPart) / freq1.QuadPart;
    cout << "use Time of smul: " << exe_time << "ms" << endl;
    // fprintf(stdout, "Your program executed time is %fms.\n", exe_time);
#endif
}

void testWithPrint() {
    srand((unsigned)time(NULL));

    Matrix<double> m = Matrix<double>::rand(D, D);
    Matrix<double> n = Matrix<double>::rand(D, D);
    // for (unsigned i = 0; i < D; ++i) {
    //     for (unsigned j = 0; j < D; ++j) {
    //         // cin >> m[i][j];
    //     }
    // }
    // for (unsigned i = 0; i < D; ++i) {
    //     for (unsigned j = 0; j < D; ++j) {
    //         // cin >> n[i][j];
    //     }
    // }
    cout << "print" << endl;
    cout << m << endl;
    cout << n << endl;
    cout << "swap" << endl;
    m.swap(n);
    cout << "print" << endl;
    cout << m << endl;
    cout << n << endl;
    cout << "get contribute" << endl;
    cout << m.row() << ' ' << m.col() << endl;
    cout << m.max() << ' ' << m.min() << ' ' << m.avg() << endl;
    cout << "operator +" << endl;
    cout << m + n << endl;
    cout << "operator -" << endl;
    cout << m - n << endl;

    
    cout << "operator *" << endl;
    cout << m * n << endl;
    cout << "operator /" << endl;
    cout << m / n << endl;
    cout << "operator %" << endl;
    cout << m % n << endl;
    cout << "cut" << endl;
    cout << m.cut(0, m.row() / 2, 0, m.col() / 2);
    cout << "eig" << endl;
    cout << n.eig() << endl;
    // cout << "cond2" << endl;
    // cout << n.cond2() << endl;
    cout << "operator ^" << endl;
    cout << (m ^ 15) << endl;
    cout << "avg" << endl;
    cout << m.avg() << endl;
    cout << "copy" << endl;
    n = m;
    cout << "cov" << endl;
    cout << m.cov() << endl;
    cout << "inv" << endl;
    cout << m.inv() << endl;
    // cout << "hess" << endl;
    // cout << m.hess() << endl;
    cout << "cond2" << endl;
    cout << m.cond2() << endl;
    cout << "eye" << endl;
    cout << Matrix<double>::eye(D, D) << endl;
}

void testRandnWithMap() {
    std::map<int, int> hist;
    for (int n = 0; n < 10000; ++n) {
        // cout << Matrix<double>::randn() << endl;
        ++hist[std::round(Matrix<double>::randn())];
    }
    for (auto p : hist) {
        std::cout << std::fixed << std::setprecision(1) << std::setw(2)
                  << p.first << ' ' 
                  << std::string(p.second / 200, '*') << '\n';
    }
}

void test() {
    cout << Matrix<double>::randn(5, 10) << endl;
}

int main() {
    ios::sync_with_stdio(false);
    freopen("test.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    testWithPrint();
    // testRandnWithMap();
    // test();
    // AccurateTest();

    return 0;
}