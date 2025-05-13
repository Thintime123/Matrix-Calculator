#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <stdexcept>
#include <QTextStream>
#include <QStringList>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    
    // 设置窗口标题
    setWindowTitle("矩阵计算器");
    
    // 初始化操作下拉框
    ui->cmbOperation->addItem("矩阵加法 (A + B)", 0);
    ui->cmbOperation->addItem("矩阵减法 (A - B)", 1);
    ui->cmbOperation->addItem("矩阵乘法 (A * B)", 2);
    ui->cmbOperation->addItem("矩阵转置 (A^T)", 3);
    ui->cmbOperation->addItem("矩阵行列式 |A|", 4);
    ui->cmbOperation->addItem("矩阵求逆 (A^-1)", 5);
    ui->cmbOperation->addItem("矩阵的秩", 6);
    ui->cmbOperation->addItem("矩阵特征值", 7);
    ui->cmbOperation->addItem("LU分解", 8);
    ui->cmbOperation->addItem("QR分解", 9);
    
    // 默认选择矩阵A
    ui->rbMatrixA->setChecked(true);
    
    // 初始更新操作界面状态
    enableSecondMatrixInput(true);
    
    // 连接信号和槽
    connect(ui->cmbOperation, SIGNAL(currentIndexChanged(int)), this, SLOT(on_cmbOperation_currentIndexChanged(int)));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_btnCreateMatrix_clicked()
{
    try {
        int rows = ui->spnRows->value();
        int cols = ui->spnCols->value();
        
        if (rows <= 0 || cols <= 0) {
            showErrorMessage("行数和列数必须大于0");
            return;
        }
        
        // 创建矩阵
        Matrix<double> matrix(rows, cols, 0.0);
        
        // 更新当前选中的矩阵
        if (ui->rbMatrixA->isChecked()) {
            matrixA = matrix;
            updateMatrixDisplay(matrixA, ui->txtMatrixA);
        } else if (ui->rbMatrixB->isChecked()) {
            matrixB = matrix;
            updateMatrixDisplay(matrixB, ui->txtMatrixB);
        }
    } catch (const std::exception& e) {
        showErrorMessage(QString("创建矩阵失败: %1").arg(e.what()));
    }
}

void MainWindow::on_btnCalculate_clicked()
{
    try {
        // 从文本编辑器获取矩阵A
        if (!parseMatrixInput(ui->txtMatrixA, matrixA)) {
            showErrorMessage("矩阵A格式错误");
            return;
        }
        
        // 根据选择的操作执行计算
        int operation = ui->cmbOperation->currentData().toInt();
        
        switch (operation) {
            case 0: { // 矩阵加法
                if (!parseMatrixInput(ui->txtMatrixB, matrixB)) {
                    showErrorMessage("矩阵B格式错误");
                    return;
                }
                resultMatrix = matrixA + matrixB;
                break;
            }
            case 1: { // 矩阵减法
                if (!parseMatrixInput(ui->txtMatrixB, matrixB)) {
                    showErrorMessage("矩阵B格式错误");
                    return;
                }
                resultMatrix = matrixA - matrixB;
                break;
            }
            case 2: { // 矩阵乘法
                if (!parseMatrixInput(ui->txtMatrixB, matrixB)) {
                    showErrorMessage("矩阵B格式错误");
                    return;
                }
                resultMatrix = matrixA * matrixB;
                break;
            }
            case 3: { // 矩阵转置
                resultMatrix = matrixA.transpose();
                break;
            }
            case 4: { // 行列式
                if (!matrixA.isSquare()) {
                    showErrorMessage("只有方阵才能计算行列式");
                    return;
                }
                double det = matrixA.determinant();
                ui->txtResult->setPlainText(QString::number(det));
                return;
            }
            case 5: { // 矩阵求逆
                if (!matrixA.isSquare()) {
                    showErrorMessage("只有方阵才能计算逆矩阵");
                    return;
                }
                resultMatrix = matrixA.inverse();
                break;
            }
            case 6: { // 矩阵的秩
                int rank = MatrixOperation<double>::rank(matrixA);
                ui->txtResult->setPlainText(QString("矩阵的秩: %1").arg(rank));
                return;
            }
            case 7: { // 特征值
                if (!matrixA.isSquare()) {
                    showErrorMessage("只有方阵才能计算特征值");
                    return;
                }
                auto eigenvalues = MatrixOperation<double>::eigenvalues(matrixA);
                QString result = "特征值:\n";
                for (const auto& val : eigenvalues) {
                    if (std::abs(val.imag()) < 1e-10) {
                        result += QString::number(val.real()) + "\n";
                    } else {
                        result += QString::number(val.real()) + " + " + 
                                  QString::number(val.imag()) + "i\n";
                    }
                }
                ui->txtResult->setPlainText(result);
                return;
            }
            case 8: { // LU分解
                if (!matrixA.isSquare()) {
                    showErrorMessage("只有方阵才能进行LU分解");
                    return;
                }
                auto [L, U] = MatrixOperation<double>::luDecomposition(matrixA);
                QString result = "L 矩阵:\n" + L.toString() + "\nU 矩阵:\n" + U.toString();
                ui->txtResult->setPlainText(result);
                return;
            }
            case 9: { // QR分解
                auto [Q, R] = MatrixOperation<double>::qrDecomposition(matrixA);
                QString result = "Q 矩阵 (正交矩阵):\n" + Q.toString() + "\nR 矩阵 (上三角矩阵):\n" + R.toString();
                ui->txtResult->setPlainText(result);
                return;
            }
            default:
                showErrorMessage("不支持的操作");
                return;
        }
        
        // 显示结果矩阵
        displayMatrix(resultMatrix);
        
    } catch (const std::exception& e) {
        showErrorMessage(QString("计算错误: %1").arg(e.what()));
    }
}

void MainWindow::on_btnClear_clicked()
{
    ui->txtMatrixA->clear();
    ui->txtMatrixB->clear();
    ui->txtResult->clear();
}

void MainWindow::on_rbMatrixA_toggled(bool checked)
{
    if (checked) {
        // 显示矩阵A的尺寸
        ui->spnRows->setValue(matrixA.getRows());
        ui->spnCols->setValue(matrixA.getCols());
    }
}

void MainWindow::on_rbMatrixB_toggled(bool checked)
{
    if (checked) {
        // 显示矩阵B的尺寸
        ui->spnRows->setValue(matrixB.getRows());
        ui->spnCols->setValue(matrixB.getCols());
    }
}

void MainWindow::on_btnSaveMatrix_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, "保存矩阵", "", "文本文件 (*.txt)");
    if (fileName.isEmpty()) {
        return;
    }
    
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        showErrorMessage("无法打开文件进行写入");
        return;
    }
    
    QTextStream out(&file);
    
    // 保存当前显示的矩阵
    if (ui->rbMatrixA->isChecked()) {
        out << matrixA.toString();
    } else if (ui->rbMatrixB->isChecked()) {
        out << matrixB.toString();
    } else {
        out << resultMatrix.toString();
    }
    
    file.close();
}

void MainWindow::on_btnLoadMatrix_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, "加载矩阵", "", "文本文件 (*.txt)");
    if (fileName.isEmpty()) {
        return;
    }
    
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        showErrorMessage("无法打开文件进行读取");
        return;
    }
    
    QTextStream in(&file);
    QString content = in.readAll();
    file.close();
    
    try {
        Matrix<double> loadedMatrix = Matrix<double>::fromString(content);
        
        // 更新当前选中的矩阵
        if (ui->rbMatrixA->isChecked()) {
            matrixA = loadedMatrix;
            updateMatrixDisplay(matrixA, ui->txtMatrixA);
            ui->spnRows->setValue(matrixA.getRows());
            ui->spnCols->setValue(matrixA.getCols());
        } else if (ui->rbMatrixB->isChecked()) {
            matrixB = loadedMatrix;
            updateMatrixDisplay(matrixB, ui->txtMatrixB);
            ui->spnRows->setValue(matrixB.getRows());
            ui->spnCols->setValue(matrixB.getCols());
        }
    } catch (const std::exception& e) {
        showErrorMessage(QString("加载矩阵失败: %1").arg(e.what()));
    }
}

void MainWindow::on_cmbOperation_currentIndexChanged(int index)
{
    // 根据操作类型启用或禁用矩阵B输入
    int operation = ui->cmbOperation->itemData(index).toInt();
    
    // 矩阵加法、减法、乘法需要矩阵B
    bool needMatrixB = (operation == 0 || operation == 1 || operation == 2);
    enableSecondMatrixInput(needMatrixB);
}

void MainWindow::updateMatrixDisplay(const Matrix<double>& matrix, QPlainTextEdit* textEdit)
{
    textEdit->setPlainText(matrix.toString());
}

void MainWindow::displayMatrix(const Matrix<double>& matrix)
{
    ui->txtResult->setPlainText(matrix.toString());
}

bool MainWindow::parseMatrixInput(QPlainTextEdit* textEdit, Matrix<double>& matrix)
{
    QString text = textEdit->toPlainText().trimmed();
    if (text.isEmpty()) {
        return false;
    }
    
    try {
        matrix = Matrix<double>::fromString(text);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

void MainWindow::showErrorMessage(const QString& message)
{
    QMessageBox::critical(this, "错误", message);
}

void MainWindow::enableSecondMatrixInput(bool enable)
{
    ui->lblMatrixB->setEnabled(enable);
    ui->txtMatrixB->setEnabled(enable);
}
