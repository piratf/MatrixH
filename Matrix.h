#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <cassert>
#include <algorithm>
#include <complex>
#include <iomanip>
#include <random>
#include <omp.h>
#include <functional>
using std::cout;
using std::endl;

namespace piratfMatrixH {


    using vecSizeT = size_t;

    // 单例随机数生成器
    static std::mt19937 *_genPtr;

    static std::normal_distribution<double> normDis(0, 1);
    static std::uniform_real_distribution<double> unifDoubleDis(0, 1);
    static std::uniform_real_distribution<float> unifFloatDis(0, 1);
    static std::uniform_int_distribution<int> unifIntDis(0, 65535);

    // 辅助代码
    // =========================================
    // 用来判断类别是否相同
    // =========================================
    template <typename T, typename U>
    struct SameType {
        static const bool isSame = false;
    };

    template <typename T>
    struct SameType<T, T> {
        static const bool isSame = true;
    };

    // 用来判断类别是否是复数
    // =========================================
    template <typename T>
    struct ComplexType {
        static const bool isComplex = false;
    };

    template <typename T>
    struct ComplexType <std::complex<T> > {
        static const bool isComplex = true;
    };

    // 重载判断复数的大小
    // =========================================
    template <typename T>
    bool operator < (const std::complex<T> &lhs, const std::complex<T> &rhs) {
        return lhs.real() < rhs.real();
    }

    template <typename T>
    bool operator < (const std::complex<T> &lhs, const T &rhs) {
        return lhs.real() < rhs.real();
    }


    class Exception {

      public:
        explicit Exception(const std::string &_m): message(_m) {
        }

        void printMessage() const {
            std::cout << message << std::endl;
        }

      private:
        std::string message;
    };

    // 矩阵的逆不存在 - 异常类
    class NoInverseException : public Exception {

    };

    // 复数矩阵无法进行这个运算 - 异常类
    class ComplexTypeCantDo : public Exception {

    };

    // 辅助代码结束
    // =========================================

    template <typename T>
    class Matrix {
      public:
        Matrix() = default;
        Matrix(Matrix<T> &&other);
        Matrix(const Matrix<T> &other);
        Matrix(vecSizeT _x);
        Matrix(vecSizeT _x, vecSizeT _y);
        Matrix(std::vector<std::vector<T> > dvec);

        static void inline setRand();
        static void inline setRand(std::mt19937 *_p);
        static T inline randn();
        static Matrix<T> inline randn(vecSizeT n);
        static Matrix<T> inline randn(vecSizeT r, vecSizeT c);

        static T inline rand();
        static Matrix<T> inline rand(vecSizeT n);
        static Matrix<T> inline rand(vecSizeT r, vecSizeT c);

        // static T inline random();

        void inline clear();
        void inline swap(Matrix<T> &rhs) {
            data.swap(rhs.getData());
        }
        void inline setZero();
        std::vector<std::vector<T> > &getData();
        const std::vector<std::vector<T> > &getData() const ;

        vecSizeT inline col() const;
        vecSizeT inline row() const;

        Matrix<T> inline cut(vecSizeT rs, vecSizeT re, vecSizeT cs, vecSizeT ce) const;
        Matrix<T> inline row(vecSizeT index) const;

        Matrix<double> toDouble() const;
        Matrix<double> inv() const;
        Matrix<T> cov() const;
        Matrix<double> hess() const;
        Matrix<std::complex<double>> eig(double eps = 1e-12, unsigned LOOP = 100000) const;

        void mul(Matrix<T> &ret, const Matrix<T> &other) const ;
        void smul(Matrix<T> &ret, const Matrix<T> &other) const ;
        void StrassenMul(vecSizeT rs, vecSizeT re, vecSizeT cs, vecSizeT ce, const Matrix<T> &other, Matrix<T> &ret) const ;
        Matrix<double> div(const Matrix<T> &other) const ;
        Matrix<double> rdiv(const Matrix<T> &other) const ;

        T inline max() const;
        T inline min() const;
        T inline avg() const;
        std::complex<double> cond2() const;

        static T vectorDotProduct(const std::vector<T> lhs, std::vector<T> rhs);
        static Matrix<T> cov(const Matrix<T> &input);
        static Matrix<T> eye(vecSizeT _x, vecSizeT _y);
        static T inline avg(const std::vector<T> &vec);

        bool inline empty() const;
        bool inline isSquare() const;
        bool inline isSingular() const;

        std::complex<double> det() const;
        Matrix<T> getTransposition() const;

        std::vector<T> &operator [](vecSizeT index);
        const std::vector<T> &operator [](vecSizeT index) const;
        Matrix<T> operator = (const Matrix<T> &other);
        Matrix<T> operator = (Matrix<T> &&other);
        Matrix<T> inline operator + (const Matrix<T> &other) const ;
        Matrix<T> inline operator - (const Matrix<T> &other) const ;
        Matrix<T> inline operator * (const Matrix<T> &other) const ;
        Matrix<double> inline operator / (const Matrix<T> &other) const ;
        Matrix<double> inline operator % (const Matrix<T> &other) const ;
        Matrix<T> inline operator *= (const Matrix<T> &other);
        Matrix<T> inline operator += (const Matrix<T> &other);
        Matrix<T> inline operator -= (const Matrix<T> &other);
        Matrix<double> inline operator /= (const Matrix<T> &other);
        /**
        * 矩阵快速幂
        * 参数需要一个当前类类型成员，需要使用 friend 实现
        * 友元模板的原因需要在声明处实现
        */
        friend Matrix<T> inline operator ^ (const Matrix<T> &mat, unsigned exponent) {
            Matrix<T> ans = eye(mat.row(), mat.col());
            Matrix<T> src = mat;

            for (; exponent; exponent >>= 1) {
                if (exponent & 1) {
                    ans *= src;
                }

                src *= src;
            }

            return ans;
        }

        /**
        * 友元模板的原因需要放在声明处实现
        */
        friend std::ostream &operator << (std::ostream &o, const Matrix<T> &mat) {
            // 判断要输出的数是否是 double 类型
            bool dFlag = SameType<T, double>::isSame;

            if (dFlag) {
                using std::ios;
                using std::setprecision;
                o << setiosflags(ios::right) << setiosflags(ios::scientific) << setprecision(4);
            }

            for (vecSizeT i = 0; i != mat.data.size(); ++i) {
                for (vecSizeT j = 0; j != mat.data[i].size(); ++j) {
                    if (dFlag) {
                        o << std::setw(12) << mat.data[i][j] << ' ';
                    } else {
                        o << mat.data[i][j] << ' ';
                    }
                }

                if (i < mat.data.size() - 1) {
                    o << '\n';
                }
            }

            o << std::endl;
            return o;
        }

      private:
        std::vector< std::vector<T> > data;
        static std::mt19937 gen;
    };

    /**
     * 移动构造
     */
    template <typename T>
    Matrix<T>::Matrix(Matrix<T> &&other) {
        data.swap(other.data);
    }

    /**
     * 拷贝构造
     */
    template <typename T>
    Matrix<T>::Matrix(const Matrix<T> &other) {
        data = other.getData();
    }

    /**
     * 建立一个n行的空矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(vecSizeT _x) {
        std::vector< std::vector<T> > temp(_x);
        data = temp;
    }

    /**
     * 通过矩阵的行列数构造空矩阵
     */
    template <typename T>
    Matrix<T>::Matrix(vecSizeT _x, vecSizeT _y) {
        std::vector< std::vector<T> > temp(_x, std::vector<T>(_y));
        data = temp;
    }

    /**
     * 矩阵深拷贝构造函数，调用了vector的拷贝方法
     */
    template <typename T>
    Matrix<T>::Matrix(std::vector<std::vector<T> > dvec) {
        data = dvec;
    }

    /**
     * 生成一个指定大小的单位阵
     * @author piratf
     * @param  _x 行数
     * @param  _y  列数
     * @return    一个单位矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::eye(vecSizeT _x, vecSizeT _y) {
        Matrix<T> mat(_x, _y);

        for (vecSizeT i = 0; i < _x; ++i) {
            for (vecSizeT j = 0; j < _y; ++j) {
                if (i == j) {
                    mat.data[i][j] = 1;
                }
            }
        }

        return mat;
    }

    /**
     * 获得矩阵的转置
     * @author piratf
     * @return 新的矩阵，内容为原矩阵的转置
     */
    template <typename T>
    Matrix<T> Matrix<T>::getTransposition() const {
        decltype(data.size()) sizeRow = data.size();

        if (sizeRow == 0) {
            std::cerr << "error** Matrix<T>::getTransposition -> empty Matrix!" << std::endl;
        }

        using vecSizeT = decltype(data.size());
        vecSizeT sizeCol = data[0].size();

        Matrix tran(sizeCol, sizeRow);

        for (vecSizeT i = 0; i < sizeRow; ++i) {
            for (vecSizeT j = 0; j < sizeCol; ++j) {
                tran.data[j][i] = data[i][j];
            }
        }

        return tran;
    }

    /**
     * 静态函数：获取两个向量的点乘结果
     * @author piratf
     * @param  lhs 向量1
     * @param  rhs 向量2
     * @return     double:点乘的结果
     */
    template <typename T>
    T Matrix<T>::vectorDotProduct(const std::vector<T> lhs, std::vector<T> rhs) {
        T ans = 0;

        for (decltype(lhs.size()) i = 0; i != lhs.size(); ++i) {
            ans += lhs[i] * rhs[i];
        }

        return ans;
    }

    /**
     * 一个向量的均值
     */
    template <typename T>
    T Matrix<T>::avg(const std::vector<T> &vec) {
        T sum = 0;

        for (T var : vec) {
            sum += var;
        }

        return sum / vec.size();
    }

    /**
     * 矩阵的协方差矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::cov() const {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!empty());
#endif
        const vecSizeT sizeRow = row();
        const vecSizeT sizeCol = col();
        auto mat = this ->getTransposition();
        std::vector<T> avgVec;

        // 对于每一行求其均值
        for (auto &row : mat.data) {
            avgVec.push_back(avg(row));
        }

        // 获得协方差矩阵参数
        Matrix temp(sizeRow, sizeCol);

        for (vecSizeT i = 0; i != sizeRow; ++i) {
            for (vecSizeT j = 0; j != sizeCol; ++j) {
                temp.data[i][j] = mat.data[i][j] - avgVec[i];
            }
        }

        // 获得协方差矩阵
        Matrix cov(sizeRow, sizeRow);

        for (vecSizeT i = 0; i != sizeRow; ++i) {
            for (vecSizeT j = 0; j != sizeRow; ++j) {
                cov.data[i][j] = Matrix<T>::vectorDotProduct(temp.data[i], temp.data[j]) / (sizeCol - 1);
            }
        }

        return cov;
    }

    /**
     * 获得样本矩阵的协方差矩阵
     * 参数矩阵中的各样本需按行排列
     * @author piratf
     * @return 一个新的协方差矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::cov(const Matrix<T> &input) {
        return input.cov();
    }

    /**
     * 判断是否是空矩阵
     * @author piratf
     * @return 1 : 0 -> 空矩阵 : 不是空矩阵
     */
    template <typename T>
    bool Matrix<T>::empty() const {
        return !data.size();
    }

    /**
     * 判断矩阵是否是方阵
     * @author piratf
     * @return 1 : 0 -> 方阵 : 不是方阵
     */
    template <typename T>
    bool Matrix<T>::isSquare() const {
        return empty() ? 0 : data.size() == data[0].size();
    }

    /**
     * 求矩阵行列式
     * @author piratf
     * @return double: 行列式的值
     */
    template <typename T>
    std::complex<double> Matrix<T>::det() const {
        // 所有特征根的乘积
        auto e = eig();
        std::complex<double> ret = e[0];

        for (vecSizeT i = 1; i < data.size(); ++i) {
            ret *= e[i];
        }

        return ret;
    }

    /**
     * 判断当前矩阵是否是奇异矩阵
     * @author piratf
     * @return 1: 是奇异矩阵 0: 不是奇异矩阵
     */
    template <typename T>
    bool Matrix<T>::isSingular() const {
        return det() ? false : true;
    }

    /**
     * 部分主元的高斯消去法求逆
     */
    template <typename T>
    Matrix<double> Matrix<T>::inv() const {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!empty());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(isSquare());
#endif
        vecSizeT i, j, k, len = row();
        double maxVal, temp;
        //将A矩阵存放在临时矩阵中
        Matrix<double> TMat;

        if (SameType<T, double>::isSame) {
            TMat = *this;
        } else {
            TMat = this ->toDouble();
        }

        //初始化ans矩阵为单位阵
        Matrix<double> ans = Matrix<double>::eye(row(), col());

        for (i = 0; i < len; i++) {
            //寻找主元
            maxVal = TMat[i][i];
            k = i;

            for (j = i + 1; j < len; j++) {
                if (std::abs(TMat[j][i]) > std::abs(maxVal)) {
                    maxVal = TMat[j][i];
                    k = j;
                }
            }

            //如果主元所在行不是第i行，进行行交换
            if (k != i) {
                TMat[i].swap(TMat[k]);
                ans[i].swap(ans[k]);
            }

            assert(cond2().real() < 1e11);
            //判断主元是否为0, 若是, 则矩阵A不是满秩矩阵,不存在逆矩阵
            // if (cond2().real() > 1e10) {
            //     throw (Exception("ERROR** Matrix::inv -> there is no inverse matrix!"));
            // }
            //消去A的第i列除去i行以外的各行元素
            temp = TMat[i][i];

            for (j = 0; j < len; j++) {
                TMat[i][j] = TMat[i][j] / temp;     //主对角线上的元素变为1
                ans[i][j] = ans[i][j] / temp;       //伴随计算
            }

            // 遍历行
            for (j = 0; j < len; j++) {
                // 不是第i行
                if (j != i) {
                    temp = TMat[j][i];

                    // 第j行元素 - i行元素 * j列i行元素
                    for (k = 0; k < len; k++) {
                        TMat[j][k] -= TMat[i][k] * temp;
                        ans[j][k] -= ans[i][k] * temp;
                    }
                }
            }
        }

        return ans;
    }

    /**
     * 矩阵的行数
     */
    template <typename T>
    vecSizeT Matrix<T>::row() const {
        return data.size();
    }

    /**
     * 矩阵的列数
     */
    template <typename T>
    vecSizeT Matrix<T>::col() const {
        if (data.size()) {
            return data[0].size();
        } else {
            return 0;
        }
    }

    /**
     * 按行号取某一行，返回新的矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::row(vecSizeT index) const {
        return Matrix<T>(data[index]);
    }

    /**
     * 取矩阵中的某一行
     */
    template <typename T>
    std::vector<T> &Matrix<T>::operator [](vecSizeT index) {
        assert(index >= 0 && index < data.size());
        return data[index];
    }

    /**
     * 取矩阵中的某一行
     */
    template <typename T>
    const std::vector<T> &Matrix<T>::operator [](vecSizeT index) const {
        assert(index >= 0 && index < data.size());
        return data[index];
    }

    /**
     * 生成当前矩阵的 Hessenberg 形式，以新矩阵返回
     */
    template <typename T>
    Matrix<double> Matrix<T>::hess() const {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!empty());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(isSquare());
#endif
        Matrix<double> A = toDouble();

        vecSizeT n = data.size();
        vecSizeT i, j, k;
        Matrix<double> ret(n, n);
        T temp = 0;
        vecSizeT max;

        for (k = 1; k < n - 1; ++k) {
            i = k - 1;
            max = k;
            temp = std::abs(A[k][i]);

            for (j = k + 1; j < n; ++j) {
                if (temp < std::abs(A[j][i])) {
                    temp = std::abs(A[j][i]);
                    max = j;
                }
            }

            ret[0][0] = A[max][i];
            i = max;

            if (ret[0][0]) {
                if (i != k) {
                    for (j = k - 1; j < n; ++j) {
                        temp = A[i][j];
                        A[i][j] = A[k][j];
                        A[k][j] = temp;
                    }

                    for (j = 0; j < n; ++j) {
                        temp = A[j][i];
                        A[j][i] = A[j][k];
                        A[j][k] = temp;
                    }
                }

                for (i = k + 1; i < n; ++i) {
                    temp = A[i][k - 1] / ret[0][0];
                    A[i][k - 1] = 0;

                    for (vecSizeT j = k; j < n; ++j) {
                        A[i][j] -= temp * A[k][j];
                    }

                    for (j = 0; j < n; ++j) {
                        A[j][k] += temp * A[j][i];
                    }
                }
            }
        }

        return A;
    }

    /**
     * 返回矩阵的全部特征根，以复数表示
     * QR 分解法
     */
    template<typename T>
    Matrix<std::complex<double>> Matrix<T>::eig(double eps, unsigned LOOP) const {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!empty());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(isSquare());
#endif
        unsigned loop = LOOP;
        const vecSizeT n = data.size();
        vecSizeT m = n;
        Matrix<double> A = hess();
        Matrix<double> ret(n, 2);
        vecSizeT i, j, k, t;
        double tempVar, indexVar, sign, p, q;
        double r, x, s, e, f, z, y, temp;
        double num;

        while (m != 0) {
            t = m - 1;

            while (t > 0) {
                temp = std::abs(A[t - 1][t - 1]);
                temp += std::abs(A[t][t]);
                temp *= eps;

                if (std::abs(A[t][t - 1]) > temp) {
                    --t;
                } else {
                    break;
                }
            }

            if (t == m - 1) {
                ret[m - 1][0] = A[m - 1][m - 1];
                ret[m - 1][1] = 0;
                m -= 1;
                loop = LOOP;
            } else if (t == m - 2) {
                tempVar = -A[m - 1][m - 1] - A[m - 2][m - 2];
                indexVar = A[m - 1][m - 1] * A[m - 2][m - 2]
                           - A[m - 1][m - 2] * A[m - 2][m - 1];
                num = tempVar * tempVar - 4 * indexVar;
                y = std::sqrt(std::abs(num));

                if (num > 0) {
                    sign = 1;

                    if (tempVar < 0) {
                        sign = -1;
                    }

                    ret[m - 1][0] = -(tempVar + sign * y) / 2;
                    ret[m - 1][1] = 0;
                    ret[m - 2][0] = indexVar / ret[m - 1][0];
                    ret[m - 2][1] = 0;
                } else {
                    ret[m - 1][0] = -tempVar / 2;
                    ret[m - 2][0] = -tempVar / 2;
                    ret[m - 1][1] = y / 2;
                    ret[m - 2][1] = -y / 2;
                }

                m -= 2;
                loop = LOOP;
            } else {
                if (loop < 1) {
                    return Matrix<std::complex<double> >();
                }

                --loop;
                j = t + 2;

                while (j < m) {
                    A[j][j - 2] = 0;
                    ++j;
                }

                j = t + 3;

                while (j < m) {
                    A[j][j - 3] = 0;
                    ++j;
                }

                k = t;

                while (k < m - 1) {
                    if (k != t) {
                        p = A[k][k - 1];
                        q = A[k + 1][k - 1];

                        if (k != m - 2) {
                            r = A[k + 2][k - 1];
                        } else {
                            r = 0;
                        }
                    } else {
                        tempVar = A[m - 1][m - 1];
                        indexVar = A[m - 2][m - 2];
                        x = tempVar + indexVar;
                        y = tempVar * indexVar - A[m - 2][m - 1] * A[m - 1][m - 2];
                        p = A[t][t] * (A[t][t] - x) + A[t][t + 1] * A[t + 1][t] + y;
                        q = A[t + 1][t] * (A[t][t] + A[t + 1][t + 1] - x);
                        r = A[t + 1][t] * A[t + 2][t + 1];
                    }

                    if (p != 0 || q != 0 || r != 0) {
                        if (p < 0) {
                            sign = -1;
                        } else {
                            sign = 1;
                        }

                        s = sign * std::sqrt(p * p + q * q + r * r);

                        if (k != t) {
                            A[k][k - 1] = -s;
                        }

                        e = -q / s;
                        f = -r / s;
                        x = -p / s;
                        y = -x - f * r / (p + s);
                        num = e * r / (p + s);
                        z = -x - e * q / (p + s);

                        for (j = k; j < m; ++j) {
                            tempVar = A[k][j];
                            indexVar = A[k + 1][j];
                            p = x * tempVar + e * indexVar;
                            q = e * tempVar + y * indexVar;
                            r = f * tempVar + num * indexVar;

                            if (k != m - 2) {
                                tempVar = A[k + 2][j];
                                p += f * tempVar;
                                q += num * tempVar;
                                r += z * tempVar;
                                A[k + 2][j] = r;
                            }

                            A[k + 1][j] = q;
                            A[k][j] = p;
                        }

                        j = k + 3;

                        if (j > m - 2) {
                            j = m - 1;
                        }

                        for (i = t; i < j + 1; ++i) {
                            tempVar = A[i][k];
                            indexVar = A[i][k + 1];
                            p = x * tempVar + e * indexVar;
                            q = e * tempVar + y * indexVar;
                            r = f * tempVar + num * indexVar;

                            if (k != m - 2) {
                                tempVar = A[i][k + 2];
                                p += f * tempVar;
                                q += num * tempVar;
                                r += z * tempVar;
                                A[i][k + 2] = r;
                            }

                            A[i][k + 1] = q;
                            A[i][k] = p;
                        }
                    }

                    ++k;
                }
            }
        }

        // 返回一个复数
        Matrix<std::complex<double> >res(n, 1);

        for (vecSizeT i = 0; i < ret.row(); ++i) {
            // 判断是否是复数类型
            bool flag = ComplexType<T>::isComplex;

            if (flag) {
                res[i][0] = std::complex<double>(static_cast<std::complex<T> >(ret[i][0]).real(), static_cast<std::complex<T> >(ret[i][1]).real());
            } else {
                res[i][0] = std::complex<double>(ret[i][0], ret[i][1]);
            }
        }

        return res;
    }

    /**
     * 返回矩阵中最大的元素
     */
    template <typename T>
    T Matrix<T>::max() const {
        if (!data.size()) {
            return static_cast<T>(0);
        }

        T maxv = data[0][0];

        for (vecSizeT i = 0; i < data.size(); ++i) {
            for (vecSizeT j = 0; j < data[0].size(); ++j) {
                maxv = data[i][j] < maxv ? maxv : data[i][j];
            }
        }

        return maxv;
    }

    /**
     * 返回矩阵中最小的元素
     */
    template <typename T>
    T Matrix<T>::min() const {
        if (!data.size()) {
            return static_cast<T>(0);
        }

        T minv = data[0][0];

        for (vecSizeT i = 0; i < data.size(); ++i) {
            for (vecSizeT j = 0; j < data[0].size(); ++j) {
                minv = data[i][j] < minv ? data[i][j] : minv;
            }
        }

        return minv;
    }

    /**
     * 矩阵在 二范式 下的条件数
     */
    template <typename T>
    std::complex<double> Matrix<T>::cond2() const {
        // 获取特征值
        auto e = eig();
        return e.max() / e.min();
    }

    /**
     * 矩阵的均值
     */
    template <typename T>
    T Matrix<T>::avg() const {
        if (!empty()) {
            return static_cast<T>(0);
        }

        T sum;

        for (vecSizeT i = 0; i < data.size(); ++i) {
            for (vecSizeT j = 0; j < data[0].size(); ++j) {
                sum += data[i][j];
            }
        }

        // for (const auto row : data) {
        //     for (T var : row) {
        //         sum += var;
        //     }
        // }
        return sum / (row() * col());
    }

    /**
     * 矩阵相加赋值
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator += (const Matrix<T> &other) {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(row() == other.row() && col() == other.col());
#endif

        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT j = 0; j < col(); ++j) {
                data[i][j] += other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵相加
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator + (const Matrix<T> &other) const {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(row() == other.row() && col() == other.col());
#endif

        Matrix<T> ret(row(), col());

        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT j = 0; j < col(); ++j) {
                ret[i][j] = data[i][j] + other[i][j];
            }
        }

        return ret;
    }


    /**
     * 矩阵相减
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator -= (const Matrix<T> &other) {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(row() == other.row() && col() == other.col());
#endif

        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT j = 0; j < col(); ++j) {
                data[i][j] -= other[i][j];
            }
        }

        return *this;
    }

    /**
     * 矩阵相减
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator - (const Matrix<T> &other) const {
#ifndef EQUALSIZE
#define EQUALSIZE
        assert(data.size());
        assert(row() == other.row() && col() == other.col());
#endif

        Matrix<T> ret(row(), col());

        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT j = 0; j < col(); ++j) {
                ret[i][j] = data[i][j] - other[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵相乘
     * 普通算法
     */
    template <typename T>
    void Matrix<T>::mul(Matrix<T> &ret, const Matrix<T> &other) const {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(col() == other.row());
#endif

        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT k = 0; k < col(); ++k) {
                for (vecSizeT j = 0; j < col(); ++j) {
                    ret[i][j] += data[i][k] * other[k][j];
                }
            }
        }

        return;
    }

    /**
     * 矩阵相乘
     * 普通算法
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator * (const Matrix<T> &other) const {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(col() == other.row());
#endif
        Matrix<T> ret(col(), other.row());

        // 对称矩阵使用 斯特拉森乘法
        if (row() == col()) {
            smul(ret, other);
        } else {
            // 普通乘法
            mul(ret, other);
        }

        return ret;
    }

    /**
     * 矩阵 *=
     * 普通算法
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator *= (const Matrix<T> &other) {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(col() == other.row());
#endif
        Matrix<T> ret(col(), other.row());

        // 大的对称矩阵使用 斯特拉森乘法
        if (row() > 50 && row() == col()) {
            smul(ret, other);
        } else {
            // 普通乘法
            mul(ret, other);
        }

        return ret;
    }

    /**
     * 矩阵左除
     */
    template <typename T>
    Matrix<double> Matrix<T>::div(const Matrix<T> &other) const {
#ifndef DIVCHECK
#define DIVCHECK
        assert(data.size());
        assert(other.row() == other.col());
#endif

        return this ->toDouble() * other.inv();
    }
    /**
     * 矩阵右除
     */
    template <typename T>
    Matrix<double> Matrix<T>::rdiv(const Matrix<T> &other) const {
#ifndef DIVCHECK
#define DIVCHECK
        assert(data.size());
        assert(other.row() == other.col());
#endif

        return this ->inv() * other.toDouble();
    }

    /**
     * 矩阵相除
     * 左除
     */
    template <typename T>
    Matrix<double> Matrix<T>::operator / (const Matrix<T> &other) const {
#ifndef DIVCHECK
#define DIVCHECK
        assert(data.size());
        assert(other.row() == other.col());
#endif

        return div(other);
    }

    /**
     * 矩阵相除
     * 左除
     */
    template <typename T>
    Matrix<double> Matrix<T>::operator /= (const Matrix<T> &other) {
#ifndef DIVCHECK
#define DIVCHECK
        assert(data.size());
        assert(other.row() == other.col());
#endif

        return *this = div(other);
    }

    /**
     * 矩阵相除
     * 右除
     */
    template <typename T>
    Matrix<double> Matrix<T>::operator % (const Matrix<T> &other) const {
#ifndef DIVCHECK
#define DIVCHECK
        assert(data.size());
        assert(other.row() == other.col());
#endif

        return rdiv(other);
    }

    /**
     * 转换为 double 矩阵
     */
    template <typename T>
    Matrix<double> Matrix<T>::toDouble() const {
        Matrix<double> mat(row(), col());

        // 如果矩阵是复数，则只取实数部分，忽略虚数部分
        // 未实现
        for (vecSizeT i = 0; i < row(); ++i) {
            for (vecSizeT j = 0; j < col(); ++j) {
                mat[i][j] = static_cast<double>(data[i][j]);
            }
        }

        return mat;
    }

    /**
     * 斯特拉森乘法主函数，需要两个 n * n 矩阵
     */
    template <typename T>
    void Matrix<T>::smul(Matrix<T> &ret, const Matrix<T> &other) const {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(col() == other.row());
#endif
#ifndef ISSQUARE
#define ISSQUARE
        assert(isSquare());
#endif
#ifndef OTHERISSQUARE
#define OTHERISSQUARE
        assert(other.isSquare());
#endif
        vecSizeT n = row();
        StrassenMul(0, n, 0, n, other, ret);
    }

    /**
     * 斯特拉森乘法，复杂度 n ^ 2.80
     */
    template <typename T>
    void Matrix<T>::StrassenMul(vecSizeT rs, vecSizeT re, vecSizeT cs, vecSizeT ce, const Matrix<T> &other, Matrix<T> &ret) const {
#ifndef MULCHECK
#define MULCHECK
        assert(data.size());
        assert(col() == other.row());
#endif

        if (re - rs == 2 && ce - cs == 2) {
            vecSizeT rs1 = -~rs;
            vecSizeT cs1 = -~cs;
            T P1 = data[rs][cs] * (other[rs][cs1] - other[rs1][cs1]);
            T P2 = (data[rs][cs] + data[rs][cs1]) * other[rs1][cs1];
            T P3 = (data[rs1][cs] + data[rs1][cs1]) * other[rs][cs];
            T P4 = data[rs1][cs1] * (other[rs1][cs] - other[rs][cs]);
            T P5 = (data[rs][cs] + data[rs1][cs1]) * (other[rs][cs] + other[rs1][cs1]);
            T P6 = (data[rs][cs1] - data[rs1][cs1]) * (other[rs1][cs] + other[rs1][cs1]);
            T P7 = (data[rs][cs] - data[rs1][cs]) * (other[rs][cs] + other[rs][cs1]);
            ret[rs][cs] = P5 + P4 - P2 + P6;
            ret[rs][cs1] = P1 + P2;
            ret[rs1][cs] = P3 + P4;
            ret[rs1][cs1] = P1 + P5 - P3 - P7;
        } else if (re - rs < 2 || rs - rs < 2) {
            for (vecSizeT i = rs; i < re; ++i) {
                for (vecSizeT k = cs; k < ce; ++k) {
                    for (vecSizeT j = cs; j < ce; ++j) {
                        ret[i][j] += data[i][k] * other[k][j];
                    }
                }
            }
        } else {
            vecSizeT rm = rs + ((re - rs) / 2);
            vecSizeT cm = cs + ((ce - cs) / 2);
            StrassenMul(rs, rm, cs, cm, other, ret);
            StrassenMul(rm, re, cs, cm, other, ret);
            StrassenMul(rs, rm, cm, ce, other, ret);
            StrassenMul(rm, re, cm, ce, other, ret);
        }
    }

    /**
     * 从矩阵中取一部分
     * 从 0 开始，到 re 结束
     * 返回新的另一个矩阵，不影响原矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::cut(vecSizeT rs, vecSizeT re, vecSizeT cs, vecSizeT ce) const {
#ifndef NOTEMPTY
#define NOTEMPTY
        assert(!empty());
#endif
        assert(re < row() && rs >= 0 && ce < col() && cs >= 0);
        Matrix<T> ret(-~re - rs, -~ce - cs);

        for (vecSizeT i = rs, ri = 0; i <= re; ++i, ++ri) {
            for (vecSizeT j = cs, rj = 0; j <= ce; ++j, ++rj) {
                ret[ri][rj] = data[i][j];
            }
        }

        return ret;
    }

    /**
     * 矩阵清空
     */
    template <typename T>
    void Matrix<T>::clear() {
        data.clear();
    }

    /**
     * 矩阵清零
     */
    template <typename T>
    void Matrix<T>::setZero() {
        vecSizeT n = row();
        vecSizeT m = col();

        for (vecSizeT i = 0; i < n; ++i) {
            std::memset(&data[i][0], 0, sizeof(T) * m);
        }
    }

    /**
     * 拷贝矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator = (const Matrix<T> &other) {
        data = other.getData();
        return *this;
    }

    /**
     * 移动矩阵
     */
    template <typename T>
    Matrix<T> Matrix<T>::operator = (Matrix<T> &&other) {
        data.swap(other.getData());
        return *this;
    }

    /**
     * 获取 data 成员，可用于整块更新
     */
    template <typename T>
    std::vector<std::vector<T> > &Matrix<T>::getData() {
        return data;
    }

    /**
     * 以 const 方式获取成员，可用于安全读
     */
    template <typename T>
    const std::vector<std::vector<T> > &Matrix<T>::getData() const {
        return data;
    }

    /**
     * 生成 0 - 1 之间的一个正态随机数
     */
    template <typename T>
    T Matrix<T>::randn() {
        setRand();
        return normDis(*_genPtr);
    }

    /**
     * 生成一个 n * n 正态随机矩阵
     * 输入参数为矩阵的边长
     */
    template <typename T>
    Matrix<T> Matrix<T>::randn(vecSizeT n) {
        setRand();
        Matrix<T> mat(n, n);

        for (vecSizeT i = 0; i < n; ++i) {
            for (vecSizeT j = 0; j < n; ++j) {
                mat[i][j] = normDis(*_genPtr);
            }
        }

        return mat;
    }

    /**
     * 生成一个正态随机矩阵
     * 输入参数为矩阵的行数和列数
     */
    template <typename T>
    Matrix<T> Matrix<T>::randn(vecSizeT r, vecSizeT c) {
        setRand();
        Matrix<T> mat(r, c);
        auto fun = std::bind(normDis, *_genPtr);

        for (vecSizeT i = 0; i < r; ++i) {
            for (vecSizeT j = 0; j < c; ++j) {
                mat[i][j] = normDis(*_genPtr);
            }
        }

        return mat;
    }

    /**
     * 初始化单例随机数生成器
     */
    template <typename T>
    void Matrix<T>::setRand() {
        if (!_genPtr) {
#ifdef linux
            std::random_device rd;
            _genPtr = new std::mt19937(rd());
#else
            _genPtr = new std::mt19937(std::time(NULL));
#endif
        }
    }

    /**
     * 传入初始化单例随机数生成器
     */
    template <typename T>
    void Matrix<T>::setRand(std::mt19937 *_p) {
        if (_genPtr) {
            delete _genPtr;
        }

        _genPtr = _p;
    }

    /**
     * 如果矩阵是浮点数，生成 0 - 1 之间的一个均匀分布随机数
     * 如果矩阵是整形，生成 0 - 65525 之间的一个均匀分布随机数
     */
    template <typename T>
    T Matrix<T>::rand() {
        setRand();

        if (SameType<T, double>::isSame) {
            return unifDoubleDis(*_genPtr);
        } else if (SameType<T, float>::isSame) {
            return unifFloatDis(*_genPtr);
        } else {
            return unifIntDis(*_genPtr);
        }
    }

    /**
     * 生成一个 n * n 均匀分布随机矩阵
     * 输入参数为矩阵的边长
     */
    template <typename T>
    Matrix<T> Matrix<T>::rand(vecSizeT n) {
        setRand();
        Matrix<T> mat(n, n);

        for (vecSizeT i = 0; i < n; ++i) {
            for (vecSizeT j = 0; j < n; ++j) {
                mat[i][j] = rand();
            }
        }

        return mat;
    }

    /**
     * 生成一个均匀分布随机矩阵
     * 输入参数为矩阵的行数和列数
     */
    template <typename T>
    Matrix<T> Matrix<T>::rand(vecSizeT r, vecSizeT c) {
        setRand();
        Matrix<T> mat(r, c);

        for (vecSizeT i = 0; i < r; ++i) {
            for (vecSizeT j = 0; j < c; ++j) {
                mat[i][j] = rand();
            }
        }

        return mat;
    }

}

namespace pmh = piratfMatrixH;

#endif