#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QButtonGroup>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QString>
#include <QPlainTextEdit>
#include <QFile>
#include "matrix.h"
#include "matrixoperation.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_btnCreateMatrix_clicked();
    void on_btnCalculate_clicked();
    void on_btnClear_clicked();
    void on_rbMatrixA_toggled(bool checked);
    void on_rbMatrixB_toggled(bool checked);
    void on_btnSaveMatrix_clicked();
    void on_btnLoadMatrix_clicked();
    void on_cmbOperation_currentIndexChanged(int index);

private:
    Ui::MainWindow *ui;
    Matrix<double> matrixA;
    Matrix<double> matrixB;
    Matrix<double> resultMatrix;
    
    void updateMatrixDisplay(const Matrix<double>& matrix, QPlainTextEdit* textEdit);
    void displayMatrix(const Matrix<double>& matrix);
    bool parseMatrixInput(QPlainTextEdit* textEdit, Matrix<double>& matrix);
    void showErrorMessage(const QString& message);
    void enableSecondMatrixInput(bool enable);
};
#endif // MAINWINDOW_H
