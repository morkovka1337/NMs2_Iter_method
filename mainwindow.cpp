#define _USE_MATH_DEFINES
#include <QMessageBox>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <string>


double f_right(double x, double y)
{
    return -exp(-x*y*y);
}
double mu1(double y)
{
    return (y-2)*(y-3);
}
double mu2(double y)
{
    return y*(y-2)*(y-3);
}
double mu3(double x)
{
    return (x-1)*(x-2);
}
double mu4(double x)
{
    return x* (x-1)*(x-2);
}

///////
double myu1_test(double y)
{
    return sin(M_PI *1*y);
}
double myu2_test(double y)
{
    return sin(M_PI *2*y);
}
double myu3_test(double x)
{
    return sin(M_PI *x*2);
}
double myu4_test(double x)
{
    return sin(M_PI *x*3);
}

double acc_u(double x, double y)
{
    return sin(M_PI *x*y);
}
double f_test(double x, double y)
{
    return - M_PI * M_PI * sin(M_PI * x * y) * (x * x + y * y);
}

/*std::vector<std::vector<double>> solveMinDisrapency(std::vector<std::vector<double>> &startSolution,
                                                    std::vector<std::vector<double>> &f, unsigned int n, unsigned int m,
                                                    unsigned long N, double eps)
{
    return
}*/

std::vector<std::vector<double>> solve_Zeidel(std::vector<std::vector<double>> &startSolution,
                                 std::vector<std::vector<double>> &f, unsigned int n, unsigned int m,
                                 unsigned long N, double eps)
{
    unsigned int S = 0;
    double epsMax = 0.0;
    double epsCur = 0.0;
    double a2, k2, h2;
    double a = 1, b = 2, c = 2, d = 3;
    int i, j;
    double vOld;
    double vNew;
    h2 = (n/(b-a))*(n/(b-a));
    k2 = (m/(d-c))*(m/(d-c));
    a2 = 2*(h2+k2);
    auto v = startSolution;
    while (true)
    {
        epsMax = 0;
        for (j = 1; j < m; j++)
        {
            for (i = 1; i < n; i++)
            {
                vOld = v[i][j];
                vNew = h2*(v[i+1][j] + v[i-1][j]) + k2*(v[i][j-1] + v[i][j+1]);
                vNew += f[i - 1][j - 1];
                vNew /= a2;
                epsCur = fabs(vOld - vNew);
                if (epsCur > epsMax)
                {
                    epsMax = epsCur;
                }
                v[i][j] = vNew;

            }
        }
        S += 1;
        if ((epsMax < eps) || (S >= N))
        {
            break;
        }
    }
    double r = 0.0;
    for (j = 1; j < m; j++)
    {
        for (i = 1; i < n; i++)
        {
            double rCur = a2 * v[i][j] - (v[i+1][j] + v[i-1][j]) * h2  - (v[i][j+1] + v[i][j-1]) * k2 - f[i-1][j-1];
            if (abs(rCur) > r)
            {
                r = abs(rCur);
            }
        }
    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("При решении Р.С. с помощью метода Зейделя с параметрами NMax = ") +
                              QString::number(N) + QString("\nи eps = ") +
                              QString::number(eps) + QString(" за S = ") +
                              QString::number(S) + QString(" итераций было получено решение \nс точностью eps max = ") +
                              QString::number(epsMax)+ QString(" и максимальной невязкой ||r|| =  ") +
                              QString::number(r));
    msgBox.exec();
    return v;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    // тестовая задача
    // a = 1, b = 2, c = 2, d = 3

    // Utest = sin(pi*x*y)

    unsigned int n = ui->spinBox->value();
    unsigned int m = ui->spinBox_2->value();
    unsigned int limit = ui->spinBox_3->value();
    double eps = ui->textEdit->toPlainText().toDouble();
    double h = double(1) / double(n);
    double k = double(1) / double(m);
    auto startSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    auto f = std::vector<std::vector<double>>(n - 1, std::vector<double> (m - 1, 0));
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m -1; ++j)
        {
            f[i][j] = - f_test(1 + (i + 1) * h, 2 + (j + 1) * k);
        }
    }
    for (int j = 0; j <= m; ++j)
    {
        startSolution[0][j] = myu1_test(2 + j * k);
    }
    for (int j = 0; j <= m; ++j)
    {
        startSolution[n][j] = myu2_test(2 + j * k);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[i][0] = myu3_test(1 + i * h);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[i][m] = myu4_test(1 + i * h);
    }

    auto solution = solve_Zeidel(startSolution, f, n, m, limit, eps);
    auto accurateSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            accurateSolution[i][j] = acc_u(1 + i * h, 2 + j * k);
        }
    }

    ui->tableWidget->setRowCount(m + 1);
    ui->tableWidget->setColumnCount(n + 1);

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
        }
    }
    auto maxDiffSol = 0.0;
    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            double currentElem = fabs(solution[i][j] - accurateSolution[i][j]);
            if (currentElem > maxDiffSol){
                maxDiffSol = currentElem;
            }
        }
    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("Тестовая задача решена с точностью ") + QString::number(maxDiffSol));
    msgBox.exec();
}

void MainWindow::on_pushButton_2_clicked()
{
    //основная задача
     // f = -e^(-xy^2)
    unsigned int n = ui->spinBox->value();
    unsigned int m = ui->spinBox_2->value();
    unsigned int limit = ui->spinBox_3->value();
    double eps = ui->textEdit->toPlainText().toDouble();
    double h = double(1) / double(n);
    double k = double(1) / double(m);
    auto startSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    auto f = std::vector<std::vector<double>>(n - 1, std::vector<double> (m - 1, 0));
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m -1; ++j)
        {
            f[i][j] = - f_right(1 + (i + 1) * h, 2 + (j + 1) * k);
        }
    }
    for (int j = 0; j <= m; ++j)
    {
        startSolution[0][j] = mu1(2 + j * k);
    }
    for (int j = 0; j <= m; ++j)
    {
        startSolution[n][j] = mu2(2 + j * k);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[i][0] = mu3(1 + i * h);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[i][m] = mu4(1 + i * h);
    }

    auto solution = solve_Zeidel(startSolution, f, n, m, limit, eps);

    ui->tableWidget->setRowCount(m + 1);
    ui->tableWidget->setColumnCount(n + 1);

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
        }
    }

    auto n2 = n * 2;
    auto m2 = m * 2;
    auto limit2 = limit;

    h = double(1) / double(n2);
    k = double(1) / double(m2);
    startSolution = std::vector<std::vector<double>>(n2 + 1, std::vector<double> (m2 + 1, 0));
    f = std::vector<std::vector<double>>(n2 - 1, std::vector<double> (m2 - 1, 0));
    for (int i = 0; i < n2 - 1; ++i)
    {
        for (int j = 0; j < m2 -1; ++j)
        {
            f[i][j] = - f_right(1 + (i + 1) * h, 2 + (j + 1) * k);
        }
    }
    for (int j = 0; j <= m2; ++j)
    {
        startSolution[0][j] = mu1(2 + j * k);
    }
    for (int j = 0; j <= m2; ++j)
    {
        startSolution[n2][j] = mu2(2 + j * k);
    }
    for (int i = 0; i <= n2; ++i)
    {
        startSolution[i][0] = mu3(1 + i * h);
    }
    for (int i = 0; i <= n2; ++i)
    {
        startSolution[i][m2] = mu4(1 + i * h);
    }

    auto solutionDouble = solve_Zeidel(startSolution, f, n2, m2, limit2, eps);
    auto maxXDiff = 0, maxYDiff = 0;
    auto maxDiffSol = 0.0;
    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            double currentElem = fabs(solution[i][j] - solutionDouble[2*i][2*j]);
            if (currentElem > maxDiffSol){
                maxDiffSol = currentElem;
                maxXDiff = i;
                maxYDiff = j;
            }
        }
    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("Основная задача решена с погрешностью ") + QString::number(maxDiffSol) +
                              QString("\nМаксимальное отклонение численного от точного наблюдается в точке: x = ") +
                              QString::number(1 + double(maxXDiff) * 2* h) +
                              QString(", y = ") +
                              QString::number(2 + double(maxYDiff) * 2* k) );
    msgBox.exec();

}
