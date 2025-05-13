#ifndef MATRIXOPERATION_H
#define MATRIXOPERATION_H

#include "matrix.h"
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

// 矩阵高级操作类
template <typename T>
class MatrixOperation {
public:
    // 计算矩阵的秩
    static int rank(const Matrix<T>& mat) {
        // 使用高斯消元法计算秩
        Matrix<T> copy = mat; // 复制矩阵，避免修改原矩阵
        return gaussianElimination(copy);
    }
    
    // 计算特征值（仅对方阵有效）
    static std::vector<std::complex<double>> eigenvalues(const Matrix<T>& mat) {
        if (!mat.isSquare()) {
            throw std::invalid_argument("只有方阵才能计算特征值");
        }
        
        // 对于小矩阵，我们可以使用简单的方法
        std::vector<std::complex<double>> result;
        
        if (mat.getRows() == 1) {
            result.push_back(std::complex<double>(mat(0, 0), 0));
            return result;
        }
        
        if (mat.getRows() == 2) {
            // 对于2x2矩阵，我们可以直接使用二次公式
            double a = 1.0;
            double b = -(mat(0, 0) + mat(1, 1));
            double c = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
            
            double discriminant = b * b - 4 * a * c;
            if (discriminant >= 0) {
                double sqrtDiscriminant = std::sqrt(discriminant);
                result.push_back(std::complex<double>((-b + sqrtDiscriminant) / (2 * a), 0));
                result.push_back(std::complex<double>((-b - sqrtDiscriminant) / (2 * a), 0));
            } else {
                double realPart = -b / (2 * a);
                double imagPart = std::sqrt(-discriminant) / (2 * a);
                result.push_back(std::complex<double>(realPart, imagPart));
                result.push_back(std::complex<double>(realPart, -imagPart));
            }
            return result;
        }
        
        // 对于更大的矩阵，我们可以使用幂法来估计最大特征值
        // 注意：这只是一个简化的实现，不能处理所有情况
        result.push_back(powerMethod(mat));
        
        // 提示用户这是一个简化版本
        throw std::runtime_error("计算完整特征值需要更复杂的算法，建议使用Eigen库");
        
        return result;
    }
    
    // LU分解（返回L和U矩阵）
    static std::pair<Matrix<T>, Matrix<T>> luDecomposition(const Matrix<T>& mat) {
        if (!mat.isSquare()) {
            throw std::invalid_argument("只有方阵才能进行LU分解");
        }
        
        size_t n = mat.getRows();
        Matrix<T> L(n, n, 0);
        Matrix<T> U(n, n, 0);
        
        // 初始化L的对角线为1
        for (size_t i = 0; i < n; ++i) {
            L(i, i) = 1;
        }
        
        // 执行LU分解
        for (size_t i = 0; i < n; ++i) {
            // 计算U的第i行
            for (size_t j = i; j < n; ++j) {
                T sum = 0;
                for (size_t k = 0; k < i; ++k) {
                    sum += L(i, k) * U(k, j);
                }
                U(i, j) = mat(i, j) - sum;
            }
            
            // 计算L的第i列
            for (size_t j = i + 1; j < n; ++j) {
                T sum = 0;
                for (size_t k = 0; k < i; ++k) {
                    sum += L(j, k) * U(k, i);
                }
                if (std::abs(U(i, i)) < 1e-10) {
                    throw std::runtime_error("LU分解中遇到除以零的情况");
                }
                L(j, i) = (mat(j, i) - sum) / U(i, i);
            }
        }
        
        return {L, U};
    }
    
    // QR分解（返回Q和R矩阵）
    static std::pair<Matrix<T>, Matrix<T>> qrDecomposition(const Matrix<T>& mat) {
        size_t m = mat.getRows();
        size_t n = mat.getCols();
        
        // 创建一个与A相同的矩阵
        Matrix<T> A = mat;
        
        // 创建Q和R矩阵
        Matrix<T> Q(m, n, 0);
        Matrix<T> R(n, n, 0);
        
        // 使用Gram-Schmidt正交化
        for (size_t j = 0; j < n; ++j) {
            // 获取当前列向量
            std::vector<T> vj(m);
            for (size_t i = 0; i < m; ++i) {
                vj[i] = A(i, j);
            }
            
            // 正交化
            for (size_t i = 0; i < j; ++i) {
                // 计算点积
                T dot = 0;
                for (size_t k = 0; k < m; ++k) {
                    dot += A(k, j) * Q(k, i);
                }
                R(i, j) = dot;
                
                // 减去投影
                for (size_t k = 0; k < m; ++k) {
                    vj[k] -= R(i, j) * Q(k, i);
                }
            }
            
            // 计算范数
            T norm = 0;
            for (size_t i = 0; i < m; ++i) {
                norm += vj[i] * vj[i];
            }
            norm = std::sqrt(norm);
            
            if (norm < 1e-10) {
                throw std::runtime_error("QR分解中遇到线性相关的列");
            }
            
            // 归一化
            for (size_t i = 0; i < m; ++i) {
                Q(i, j) = vj[i] / norm;
            }
            
            R(j, j) = norm;
        }
        
        return {Q, R};
    }
    
    // 求解线性方程组 Ax = b
    static Matrix<T> solveLinearSystem(const Matrix<T>& A, const Matrix<T>& b) {
        if (!A.isSquare()) {
            throw std::invalid_argument("系数矩阵必须是方阵");
        }
        
        if (A.getRows() != b.getRows() || b.getCols() != 1) {
            throw std::invalid_argument("方程右侧必须是与系数矩阵行数相同的列向量");
        }
        
        // 使用LU分解求解线性方程组
        auto [L, U] = luDecomposition(A);
        
        size_t n = A.getRows();
        
        // 求解 Ly = b
        Matrix<T> y(n, 1, 0);
        for (size_t i = 0; i < n; ++i) {
            T sum = 0;
            for (size_t j = 0; j < i; ++j) {
                sum += L(i, j) * y(j, 0);
            }
            y(i, 0) = b(i, 0) - sum;
        }
        
        // 求解 Ux = y
        Matrix<T> x(n, 1, 0);
        for (int i = n - 1; i >= 0; --i) {
            T sum = 0;
            for (size_t j = i + 1; j < n; ++j) {
                sum += U(i, j) * x(j, 0);
            }
            
            if (std::abs(U(i, i)) < 1e-10) {
                throw std::runtime_error("方程组求解中遇到除以零的情况");
            }
            
            x(i, 0) = (y(i, 0) - sum) / U(i, i);
        }
        
        return x;
    }
    
private:
    // 使用高斯消元法计算矩阵的秩
    static int gaussianElimination(Matrix<T>& matrix) {
        size_t rows = matrix.getRows();
        size_t cols = matrix.getCols();
        
        size_t rank = 0;
        std::vector<bool> rowUsed(rows, false);
        
        for (size_t col = 0; col < cols; ++col) {
            size_t row;
            for (row = 0; row < rows; ++row) {
                if (!rowUsed[row] && std::abs(matrix(row, col)) > 1e-10) {
                    break;
                }
            }
            
            if (row != rows) {
                ++rank;
                rowUsed[row] = true;
                
                // 将主元归一化
                T pivot = matrix(row, col);
                for (size_t j = col; j < cols; ++j) {
                    matrix(row, j) /= pivot;
                }
                
                // 消去其他行
                for (size_t i = 0; i < rows; ++i) {
                    if (i != row && std::abs(matrix(i, col)) > 1e-10) {
                        T factor = matrix(i, col);
                        for (size_t j = col; j < cols; ++j) {
                            matrix(i, j) -= factor * matrix(row, j);
                        }
                    }
                }
            }
        }
        
        return rank;
    }
    
    // 使用幂法估计最大特征值
    static std::complex<double> powerMethod(const Matrix<T>& matrix) {
        size_t n = matrix.getRows();
        
        // 初始化向量
        std::vector<double> x(n, 1.0);
        double prevLambda = 0.0;
        
        // 迭代计算
        for (int iter = 0; iter < 100; ++iter) {
            // 计算 y = Ax
            std::vector<double> y(n, 0.0);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    y[i] += matrix(i, j) * x[j];
                }
            }
            
            // 找到最大分量
            double maxComp = 0.0;
            for (size_t i = 0; i < n; ++i) {
                if (std::abs(y[i]) > std::abs(maxComp)) {
                    maxComp = y[i];
                }
            }
            
            // 归一化
            for (size_t i = 0; i < n; ++i) {
                x[i] = y[i] / maxComp;
            }
            
            // 检查收敛性
            if (std::abs(maxComp - prevLambda) < 1e-10) {
                return std::complex<double>(maxComp, 0);
            }
            
            prevLambda = maxComp;
        }
        
        return std::complex<double>(prevLambda, 0);
    }
};

#endif // MATRIXOPERATION_H
