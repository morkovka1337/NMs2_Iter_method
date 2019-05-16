#define _USE_MATH_DEFINES
#include <QMessageBox>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <string>

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
double myu5_test(double x)
{
    return sin(M_PI *x* 2.5);
}
double myu6_test(double y)
{
    return sin(M_PI * y * 1.5);
}

std::vector<std::vector<double>> computeResidual(std::vector<std::vector<double>> &v,
                                                 std::vector<std::vector<double>> &F,
                                                 bool test)
{
    double a = 1, b = 2, c = 2, d = 3;
    auto n = v.size() - 1;
    auto m = v[0].size() - 1;
    double h2 = -(n/(b-a))*(n/(b-a));
    double k2 = -(m/(d-c))*(m/(d-c));
    double a2 = -2*(h2+k2);
    auto r = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    if (test)
    {
        for (unsigned i = 0; i < r.size(); i++)
        {
            for (unsigned j = 0; j < r[i].size(); j++)
            {
                if (i >= 1 && i <= n-1 && j >= 1 && j <= m - 1)
                {
                    r[i][j] = F[i-1][j-1] - a2 * v[i][j] - h2 * (v[i+1][j] + v[i-1][j]) -
                            k2 * (v[i][j-1] + v[i][j+1]);
                }
                else if (i == 0)
                {
                    r[i][j] = -v[i][j] + myu1_test(c + j * (d - c) / m);
                }
                else if (i == n)
                {
                    r[i][j] = -v[i][j] + myu2_test(c + j * (d - c) / m);
                }
                else if (j == 0)
                {
                    r[i][j] = -v[i][j] + myu3_test(a + i * (b - a) / n);
                }
                else if (j == m)
                {
                    r[i][j] = -v[i][j] + myu4_test(a + i * (b - a) / n);
                }
            }
        }
    }
    else
    {
        for (unsigned i = 0; i < r.size(); i++)
        {
            for (unsigned j = 0; j < r[i].size(); j++)
            {
                if (i >= 1 && i <= n-1 && j >= 1 && j <= m - 1)
                {
                    r[i][j] = F[i-1][j-1] - a2 * v[i][j] - h2 * (v[i+1][j] + v[i-1][j]) -
                            k2 * (v[i][j-1] + v[i][j+1]);
                }
                else if (i == 0)
                {
                    r[i][j] = -v[i][j] + mu1(c + j * (d - c) / m);
                }
                else if (i == n)
                {
                    r[i][j] = -v[i][j] + mu2(c + j * (d - c) / m);
                }
                else if (j == 0)
                {
                    r[i][j] = -v[i][j] + mu3(a + i * (b - a) / n);
                }
                else if (j == m)
                {
                    r[i][j] = -v[i][j] + mu4(a + i * (b - a) / n);
                }
            }
        }
    }
    return r;
}

std::vector<std::vector<double>> computeResidualCut(std::vector<std::vector<double>> &v,
                                                 std::vector<std::vector<double>> &F)
{
    double a = 1, b = 2, c = 2, d = 3;
    auto n = v.size() - 1;
    auto m = v[0].size() - 1;
    double h2 = -(n/(b-a))*(n/(b-a));
    double k2 = -(m/(d-c))*(m/(d-c));
    double a2 = -2*(h2+k2);
    auto r = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
    for (int i = 0; i < n/2; ++i) {
        r.push_back(std::vector<double> (m/2 + 1, 0));
    }
    for (unsigned i = 0; i < r.size(); i++)
    {
        for (unsigned j = 0; j < r[i].size(); j++)
        {
            if ((i >= 1 && i <= n-1 && j >= 1 && j <= m/2 - 1) ||
                  (i >= 1 && i <= n/2-1 && j >= m/2-1 && j <= m - 1)  )
            {
                r[i][j] = F[i-1][j-1] - a2 * v[i][j] - h2 * (v[i+1][j] + v[i-1][j]) -
                        k2 * (v[i][j-1] + v[i][j+1]);
            }
            else if (i == 0)
            {
                r[i][j] = -v[i][j] + myu1_test(c + j * (d - c) / m);
            }
            else if (i == n)
            {
                r[i][j] = -v[i][j] + myu2_test(c + j * (d - c) / m);
            }
            else if (j == 0)
            {
                r[i][j] = -v[i][j] + myu3_test(a + i * (b - a) / n);
            }
            else if (j == m)
            {
                r[i][j] = -v[i][j] + myu4_test(a + i * (b - a) / n);
            }
            else if (j == m/2)
            {
                r[i][j] = -v[i][j] + myu5_test(a + i * (b - a) / n);
            }
            else if (i == n/2)
            {
                r[i][j] = -v[i][j] + myu6_test(c + j * (d - c) / m);
            }
        }
    }

    return r;
}


std::vector<std::vector<double>> matrixOnVecCut(std::vector<std::vector<double>> &h)
{
    double a = 1, b = 2, c = 2, d = 3;
    auto n = h.size() - 1;
    auto m = h[0].size() - 1;
    double h2 = -(n/(b-a))*(n/(b-a));
    double k2 = -(m/(d-c))*(m/(d-c));
    double a2 = -2*(h2+k2);
    auto Ah = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
    for (int i = 0; i < n/2; ++i) {
        Ah.push_back(std::vector<double> (m/2 + 1, 0));
    }for (unsigned j = 0; j < h.size(); j++)
    {
        for (unsigned i = 0; i < h[j].size(); i++)
        {
            if ((i >= 1 && i <= n-1 && j >= 1 && j <= m/2 - 1) ||
                  (i >= 1 && i <= n/2-1 && j >= m/2-1 && j <= m - 1)  )
            {
                Ah[i][j] = a2 * h[i][j] + h2 * (h[i+1][j] + h[i-1][j]) +
                        k2 * (h[i][j-1] + h[i][j+1]);
            }
            else
            {
                Ah[i][j] = h[i][j];
            }
        }
    }
    return Ah;
}

std::vector<std::vector<double>> matrixOnVec(std::vector<std::vector<double>> &h)
{
    double a = 1, b = 2, c = 2, d = 3;
    auto n = h.size() - 1;
    auto m = h[0].size() - 1;
    double h2 = -(n/(b-a))*(n/(b-a));
    double k2 = -(m/(d-c))*(m/(d-c));
    double a2 = -2*(h2+k2);
    auto Ah = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    for (unsigned j = 0; j < h.size(); j++)
    {
        for (unsigned i = 0; i < h[j].size(); i++)
        {
            if (i >= 1 && i <= n-1 && j >= 1 && j <= m - 1)
            {
                Ah[i][j] = a2 * h[i][j] + h2 * (h[i+1][j] + h[i-1][j]) +
                        k2 * (h[i][j-1] + h[i][j+1]);
            }
            else
            {
                Ah[i][j] = h[i][j];
            }
        }
    }
    return Ah;
}

double vectorOnVector(std::vector<std::vector<double>> &h, std::vector<std::vector<double>> &r)
{
    double res = 0.0;
    for (unsigned j = 0; j < r.size(); j++)
    {
        for (unsigned i = 0; i < r[j].size(); i++)
        {
            res += r[i][j] * h[i][j];
        }
    }
    return res;
}


double f_right(double x, double y)
{
    return -exp(-x*y*y);
}


double acc_u(double x, double y)
{
    return sin(M_PI *x*y);
}
double f_test(double x, double y)
{
    return - M_PI * M_PI * sin(M_PI * x * y) * (x * x + y * y);
}

std::vector<std::vector<double>> solveConjugateGradient(std::vector<std::vector<double>> &startSolution,
                                                    std::vector<std::vector<double>> &F, unsigned int n, unsigned int m,
                                                    unsigned long N, double eps, bool test)
{    
    auto v = startSolution;
    double residual = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double epsMax = 0;
    double epsCur = 0.0;
    unsigned int it = 0;
    auto r = computeResidual(v, F, test);
    auto Ah = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    auto h = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    for (unsigned j = 0; j < h.size(); j++)
    {
        for (unsigned i = 0; i < h[j].size(); i++)
        {
            h[i][j] = r[i][j];
        }
    }

    while (true)
    {

        Ah = matrixOnVec(h);
        epsMax = 0.0;
        alpha =  vectorOnVector(r, r) / vectorOnVector(h, Ah);
        /* считаем x1 */
        for (unsigned j = 0; j < m+1; j++)
        {
            for (unsigned i = 0; i < n+1; i++)
            {
                double vOld = v[i][j];
                v[i][j] = v[i][j] + alpha * h[i][j];
                epsCur = abs(vOld - v[i][j]);
                if (epsCur > epsMax)
                {
                    epsMax = epsCur;
                }
            }
        }
        residual = 0.0;
        for (unsigned j = 0; j < r.size(); j++)
        {
            for (unsigned i = 0; i < r[j].size(); i++)
            {
               if (abs(r[i][j]) > residual)
               {
                   residual = abs(r[i][j]);
               }
            }
        }
        double betaDenominator = vectorOnVector(r, r);
        for (unsigned j = 0; j < r.size(); j++)
        {
            for (unsigned i = 0; i < r[j].size(); i++)
            {
                r[i][j] = r[i][j] - alpha * Ah[i][j];
            }
        }
        beta = vectorOnVector(r, r) / betaDenominator;

        for (unsigned j = 0; j < h.size(); j++)
        {
            for (unsigned i = 0; i < h[j].size(); i++)
            {
                h[i][j] = r[i][j] + beta * h[i][j];
            }
        }

        it += 1;
        if (it >= N || abs(beta) < 1e-14 || epsMax < eps)
        {
            break;
        }


    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("При решении Р.С. с помощью метода сопряженных градиентов с параметрами NMax = ") +
                              QString::number(N) + QString("\nи eps = ") +
                              QString::number(eps) + QString(" за S = ") +
                              QString::number(it) + QString(" итераций было получено решение \nс точностью eps max = ") +
                              QString::number(epsMax) + QString(" Невязка составила: ") + QString::number(residual)

                              );
    msgBox.exec();
    return v;
}

std::vector<std::vector<double>> solveConjugateGradientCut(std::vector<std::vector<double>> &startSolution,
                                                    std::vector<std::vector<double>> &F, unsigned int n, unsigned int m,
                                                    unsigned long N, double eps)
{
    auto v = startSolution;
    double residual = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    double epsMax = 0;
    double epsCur = 0.0;
    unsigned int it = 0;
    auto r = computeResidualCut(v, F);
    auto Ah = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
    for (int i = 0; i < n/2; ++i) {
        Ah.push_back(std::vector<double> (m/2 + 1, 0));
    }
    auto h = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
    for (int i = 0; i < n/2; ++i) {
        h.push_back(std::vector<double> (m/2 + 1, 0));
    }
    for (unsigned j = 0; j < h.size(); j++)
    {
        for (unsigned i = 0; i < h[j].size(); i++)
        {
            h[i][j] = r[i][j];
        }
    }

    while (true)
    {

        Ah = matrixOnVecCut(h);
        epsMax = 0.0;
        alpha =  vectorOnVector(r, r) / vectorOnVector(h, Ah);
        /* считаем x1 */
        for (unsigned j = 0; j < m/2+1; j++)
        {
            for (unsigned i = 0; i < n+1; i++)
            {
                double vOld = v[i][j];
                v[i][j] = v[i][j] + alpha * h[i][j];
                epsCur = abs(vOld - v[i][j]);
                if (epsCur > epsMax)
                {
                    epsMax = epsCur;
                }
            }
        }
        for (unsigned j = m / 2 + 1; j < m + 1; j++)
        {
            for (unsigned i = 0; i < n / 2 +1; i++)
            {
                double vOld = v[i][j];
                v[i][j] = v[i][j] + alpha * h[i][j];
                epsCur = abs(vOld - v[i][j]);
                if (epsCur > epsMax)
                {
                    epsMax = epsCur;
                }
            }
        }
        residual = 0.0;
        for (unsigned j = 0; j < r.size(); j++)
        {
            for (unsigned i = 0; i < r[j].size(); i++)
            {
               if (abs(r[i][j]) > residual)
               {
                   residual = abs(r[i][j]);
               }
            }
        }
        double betaDenominator = vectorOnVector(r, r);
        for (unsigned j = 0; j < r.size(); j++)
        {
            for (unsigned i = 0; i < r[j].size(); i++)
            {
                r[i][j] = r[i][j] - alpha * Ah[i][j];
            }
        }
        beta = vectorOnVector(r, r) / betaDenominator;

        for (unsigned j = 0; j < h.size(); j++)
        {
            for (unsigned i = 0; i < h[j].size(); i++)
            {
                h[i][j] = r[i][j] + beta * h[i][j];
            }
        }

        it += 1;
        if (it >= N || abs(beta) < 1e-14 || epsMax < eps)
        {
            break;
        }


    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("При решении Р.С. с помощью метода сопряженных градиентов с параметрами NMax = ") +
                              QString::number(N) + QString("\nи eps = ") +
                              QString::number(eps) + QString(" за S = ") +
                              QString::number(it) + QString(" итераций было получено решение \nс точностью eps max = ") +
                              QString::number(epsMax) + QString(" Невязка составила: ") + QString::number(residual)

                              );
    msgBox.exec();
    return v;
}

std::vector<std::vector<double>> solve_Zeidel(std::vector<std::vector<double>> &startSolution,
                                 std::vector<std::vector<double>> &f, unsigned int n, unsigned int m,
                                 unsigned long N, double eps)
{
    unsigned int S = 0;
    double epsMax = 0.0;
    double epsCur = 0.0;
    double a2, k2, h2;
    double a = 1, b = 2, c = 2, d = 3;
    unsigned i, j;
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
                epsCur = abs(vOld - vNew);
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
                              QString::number(r)
                              );
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
    double a = 1, b = 2, c = 2, d = 3;

    // Utest = sin(pi*x*y)

    unsigned int n = ui->spinBox->value();
    unsigned int m = ui->spinBox_2->value();
    unsigned int limit = ui->spinBox_3->value();
    double eps = ui->textEdit->toPlainText().toDouble();
    double h = double(b-a) / double(n);
    double k = double(d-c) / double(m);
    std::vector<std::vector<double>> startSolution, f, accurateSolution;
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        startSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        startSolution = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
        for (int i = 0; i < n/2; ++i) {
            startSolution.push_back(std::vector<double> (m/2 + 1, 0));
        }
    }
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        f = std::vector<std::vector<double>>(n - 1, std::vector<double> (m - 1, 0));
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        f = std::vector<std::vector<double>>(n/2, std::vector<double> (m - 1, 0));
        for (int i = 0; i < n/2-1; ++i) {
            f.push_back(std::vector<double> (m/2 - 1, 0));
        }
    }
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        for (unsigned int i = 0; i < n - 1; ++i)
        {
            for (unsigned int j = 0; j < m - 1; ++j)
            {
                f[i][j] = - f_test(a + (i + 1) * h, c + (j + 1) * k);
            }
        }

    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        for (unsigned int i = 0; i < n - 1 ; ++i)
        {
            for (unsigned int j = 0; j < m/2 - 1 ; ++j)
            {
                 f[i][j] = - f_test(a + (i + 1) * h, c + (j + 1) * k);
            }
        }
        for (unsigned int i = 0; i < n/2 ; ++i)
        {
            for (unsigned int j = m/2-1; j < m - 1 ; ++j)
            {
                 f[i][j] = - f_test(a + (i + 1) * h, c + (j + 1) * k);
            }
        }
    }
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
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
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        for (int j = 0; j <= m; ++j)
        {
            startSolution[0][j] = myu1_test(2 + j * k);
        }
        for (int j = 0; j <= m/2; ++j)
        {
            startSolution[n][j] = myu2_test(2 + j * k);
        }
        for (int i = 0; i <= n; ++i)
        {
            startSolution[i][0] = myu3_test(1 + i * h);
        }
        for (int i = 0; i <= n/2; ++i)
        {
            startSolution[i][m] = myu4_test(1 + i * h);
        }
        for (int i = n/2; i <= n; ++i)
        {
            startSolution[i][m/2] = myu5_test(1 + i * h);
        }
        for (int j = m/2; j <= m; ++j)
        {
            startSolution[n/2][j] = myu6_test(2 + j * k);
        }
    }
    auto solution = startSolution;
    if (ui->comboBox->currentText() == "Метод Зейделя")
    {
        solution = solve_Zeidel(startSolution, f, n, m, limit, eps);
    }

    else if (ui->comboBox->currentText() == "Метод сопряженных градиентов")
    {
        if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
        {
            solution = solveConjugateGradient(startSolution, f, n, m, limit, eps, true);
        }
        else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
        {
            solution = solveConjugateGradientCut(startSolution, f, n, m, limit, eps);
        }
    }
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        accurateSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        accurateSolution = std::vector<std::vector<double>>(n/2 + 1, std::vector<double> (m + 1, 0));
        for (int i = 0; i < n/2; ++i) {
            accurateSolution.push_back(std::vector<double> (m/2 + 1, 0));
        }
    }
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        for (unsigned int i = 0; i < n + 1 ; ++i)
        {
            for (unsigned int j = 0; j < m + 1 ; ++j)
            {
                accurateSolution[i][j] = acc_u(a + i * h, c + j * k);
            }
        }
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        for (unsigned int i = 0; i < n + 1 ; ++i)
        {
            for (unsigned int j = 0; j < m/2 + 1 ; ++j)
            {
                accurateSolution[i][j] = acc_u(a + i * h, c + j * k);
            }
        }
        for (unsigned int i = 0; i < n/2 + 1 ; ++i)
        {
            for (unsigned int j = m/2+1; j < m + 1 ; ++j)
            {
                accurateSolution[i][j] = acc_u(a + i * h, c + j * k);
            }
        }
    }
    ui->tableWidget->clear();
    ui->tableWidget->setRowCount(m + 1);
    ui->tableWidget->setColumnCount(n + 1);
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        for (unsigned int i = 0; i < n+1; ++i)
        {
            for (unsigned int j = 0; j < m+1; ++j)
            {
                ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
            }
        }
    }
    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        for (unsigned int i = 0; i < n + 1 ; ++i)
        {
            for (unsigned int j = 0; j < m/2 + 1 ; ++j)
            {
                ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
            }
        }
        for (unsigned int i = 0; i < n/2 + 1 ; ++i)
        {
            for (unsigned int j = m/2+1; j < m + 1 ; ++j)
            {
                ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
            }
        }
    }
    auto maxDiffSol = 0.0;
    auto maxXDiff = 0, maxYDiff = 0;
    if (ui->comboBox_2->currentText() == "Задача на прямоугольнике")
    {
        for (unsigned int i = 0; i < n + 1 ; ++i)
        {
            for (unsigned int j = 0; j < m + 1; ++j)
            {
                double currentElem = abs(solution[i][j] - accurateSolution[i][j]);
                if (currentElem > maxDiffSol){
                    maxDiffSol = currentElem;
                    maxXDiff = i;
                    maxYDiff = j;
                }
            }
        }
    }

    else if (ui->comboBox_2->currentText() == "Задача на вырезанном прямоугольнике")
    {
        for (unsigned int i = 0; i < n + 1 ; ++i)
        {
            for (unsigned int j = 0; j < m/2 + 1 ; ++j)
            {
                double currentElem = abs(solution[i][j] - accurateSolution[i][j]);
                if (currentElem > maxDiffSol){
                    maxDiffSol = currentElem;
                    maxXDiff = i;
                    maxYDiff = j;
                }
            }
        }
        for (unsigned int i = 0; i < n/2 + 1 ; ++i)
        {
            for (unsigned int j = m/2+1; j < m + 1 ; ++j)
            {
                double currentElem = abs(solution[i][j] - accurateSolution[i][j]);
                if (currentElem > maxDiffSol){
                    maxDiffSol = currentElem;
                    maxXDiff = i;
                    maxYDiff = j;
                }
            }
        }
    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("Тестовая задача решена с точностью ") +
                              QString::number(maxDiffSol)+
                              QString(" x =  ") + QString::number(a + maxXDiff * h) +
                              QString(", y = ") + QString::number(c + maxYDiff * k)
                              );
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

    auto solution = startSolution;
    if (ui->comboBox->currentText() == "Метод Зейделя")
    {
        solution = solve_Zeidel(startSolution, f, n, m, limit, eps);
    }

    if (ui->comboBox->currentText() == "Метод сопряженных градиентов")
    {
        solution = solveConjugateGradient(startSolution, f, n, m, limit, eps, false);
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

    auto solutionDouble = startSolution;
    if (ui->comboBox->currentText() == "Метод Зейделя")
    {
        solutionDouble = solve_Zeidel(startSolution, f, n2, m2, limit, eps);
    }

    if (ui->comboBox->currentText() == "Метод сопряженных градиентов")
    {
        solutionDouble = solveConjugateGradient(startSolution, f, n2, m2, limit, eps, false);
    }
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

                              QString::number(1 + double(maxXDiff) * 2 * h) +
                              QString(", y = ") +
                              QString::number(2 + double(maxYDiff) * 2 * k) );
    msgBox.exec();

}
