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
const unsigned D = 5;

void AccurateTest() {
#ifdef __WIN32
    LARGE_INTEGER freq;
    LARGE_INTEGER start_t, stop_t;
    LARGE_INTEGER freq1;
    LARGE_INTEGER start_t1, stop_t1;
    double exe_time;
    QueryPerformanceFrequency(&freq);
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
    QueryPerformanceCounter(&stop_t1);
    exe_time = 1e3 * (stop_t1.QuadPart - start_t1.QuadPart) / freq1.QuadPart;
    cerr << "Use Time: " << exe_time << endl;
    // fprintf(stdout, "Your program executed time is %fms.\n", exe_time);
#endif
}

void testWithPrint() {
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
    cout << Matrix<int>::rand(5, 10) << endl;
}

int main() {
    ios::sync_with_stdio(false);
    freopen("test.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    // testWithPrint();
    // testRandnWithMap();
    test();

    return 0;
}