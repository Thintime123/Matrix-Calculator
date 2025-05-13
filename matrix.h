#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <complex>
#include <algorithm>
#include <QString>

template <typename T>
class Matrix {
private:
    std::vector<std::vector<T>> data;
    size_t rows;
    size_t cols;

public:
    // 构造函数
    Matrix() : rows(0), cols(0) {}
    
    Matrix(size_t r, size_t c, const T& val = T()) : rows(r), cols(c) {
        data.resize(r, std::vector<T>(c, val));
    }
    
    Matrix(const std::vector<std::vector<T>>& mat) : data(mat) {
        if (mat.empty()) {
            rows = 0;
            cols = 0;
        } else {
            rows = mat.size();
            cols = mat[0].size();
            // 检查是否每行的列数相同
            for (const auto& row : mat) {
                if (row.size() != cols) {
                    throw std::invalid_argument("每行的列数必须相同");
                }
            }
        }
    }
    
    // 拷贝构造函数
    Matrix(const Matrix<T>& other) : data(other.data), rows(other.rows), cols(other.cols) {}
    
    // 移动构造函数
    Matrix(Matrix<T>&& other) noexcept : data(std::move(other.data)), rows(other.rows), cols(other.cols) {
        other.rows = 0;
        other.cols = 0;
    }
    
    // 拷贝赋值运算符
    Matrix<T>& operator=(const Matrix<T>& other) {
        if (this != &other) {
            data = other.data;
            rows = other.rows;
            cols = other.cols;
        }
        return *this;
    }
    
    // 移动赋值运算符
    Matrix<T>& operator=(Matrix<T>&& other) noexcept {
        if (this != &other) {
            data = std::move(other.data);
            rows = other.rows;
            cols = other.cols;
            other.rows = 0;
            other.cols = 0;
        }
        return *this;
    }
    
    // 获取矩阵大小
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    
    // 访问元素
    T& operator()(size_t i, size_t j) {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("矩阵索引越界");
        }
        return data[i][j];
    }
    
    const T& operator()(size_t i, size_t j) const {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("矩阵索引越界");
        }
        return data[i][j];
    }
    
    // 基本矩阵运算
    Matrix<T> operator+(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("矩阵维度不匹配，无法进行加法运算");
        }
        
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = (*this)(i, j) + other(i, j);
            }
        }
        return result;
    }
    
    Matrix<T> operator-(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw std::invalid_argument("矩阵维度不匹配，无法进行减法运算");
        }
        
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = (*this)(i, j) - other(i, j);
            }
        }
        return result;
    }
    
    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("矩阵维度不匹配，无法进行乘法运算");
        }
        
        Matrix<T> result(rows, other.cols, T());
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                T sum = T();
                for (size_t k = 0; k < cols; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }
    
    Matrix<T> operator*(const T& scalar) const {
        Matrix<T> result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = (*this)(i, j) * scalar;
            }
        }
        return result;
    }
    
    // 转置
    Matrix<T> transpose() const {
        Matrix<T> result(cols, rows);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }
    
    // 判断是否为方阵
    bool isSquare() const {
        return rows == cols;
    }
    
    // 计算行列式(仅对方阵有效)
    T determinant() const {
        if (!isSquare()) {
            throw std::invalid_argument("只有方阵才能计算行列式");
        }
        
        if (rows == 1) {
            return (*this)(0, 0);
        }
        
        if (rows == 2) {
            return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
        }
        
        T det = T();
        int sign = 1;
        
        for (size_t j = 0; j < cols; ++j) {
            Matrix<T> subMatrix(rows - 1, cols - 1);
            
            for (size_t i = 1; i < rows; ++i) {
                for (size_t k = 0, l = 0; k < cols; ++k) {
                    if (k != j) {
                        subMatrix(i - 1, l++) = (*this)(i, k);
                    }
                }
            }
            
            det += sign * (*this)(0, j) * subMatrix.determinant();
            sign = -sign;
        }
        
        return det;
    }
    
    // 计算矩阵的逆(仅对可逆方阵有效)
    Matrix<T> inverse() const {
        if (!isSquare()) {
            throw std::invalid_argument("只有方阵才能计算逆矩阵");
        }
        
        T det = determinant();
        if (std::abs(det) < 1e-10) {
            throw std::invalid_argument("矩阵不可逆");
        }
        
        if (rows == 1) {
            Matrix<T> result(1, 1);
            result(0, 0) = T(1) / (*this)(0, 0);
            return result;
        }
        
        Matrix<T> result(rows, cols);
        
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                Matrix<T> subMatrix(rows - 1, cols - 1);
                
                for (size_t k = 0, r = 0; k < rows; ++k) {
                    if (k != i) {
                        for (size_t l = 0, c = 0; l < cols; ++l) {
                            if (l != j) {
                                subMatrix(r, c++) = (*this)(k, l);
                            }
                        }
                        r++;
                    }
                }
                
                int sign = ((i + j) % 2 == 0) ? 1 : -1;
                result(j, i) = sign * subMatrix.determinant() / det;
            }
        }
        
        return result;
    }
    
    // 转换为字符串
    QString toString() const {
        QString str;
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                str += QString::number(static_cast<double>(data[i][j])) + " ";
            }
            str += "\n";
        }
        return str;
    }
    
    // 从字符串加载矩阵
    static Matrix<T> fromString(const QString& str) {
        auto lines = str.split('\n', Qt::SkipEmptyParts);
        if (lines.isEmpty()) {
            return Matrix<T>();
        }
        
        auto firstLine = lines[0].split(' ', Qt::SkipEmptyParts);
        size_t cols = firstLine.size();
        size_t rows = lines.size();
        
        Matrix<T> result(rows, cols);
        
        for (size_t i = 0; i < rows; ++i) {
            auto values = lines[i].split(' ', Qt::SkipEmptyParts);
            if (values.size() != cols) {
                throw std::invalid_argument("矩阵格式错误，每行的列数必须相同");
            }
            
            for (size_t j = 0; j < cols; ++j) {
                result(i, j) = static_cast<T>(values[j].toDouble());
            }
        }
        
        return result;
    }
};

// 额外的运算符重载，用于标量乘法
template <typename T>
Matrix<T> operator*(const T& scalar, const Matrix<T>& matrix) {
    return matrix * scalar;
}

#endif // MATRIX_H
