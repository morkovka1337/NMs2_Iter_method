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

std::vector<std::vector<double>> solve_Conjugate_Gradient(std::vector<std::vector<double>> &startSolution,
                                                    std::vector<std::vector<double>> &f, unsigned int n, unsigned int m,
                                                    unsigned long N, double eps)
{    
    auto v = startSolution;
    double a = 1, b = 2, c = 2, d = 3;
    double residual = 0.0;
    double h2 = (n/(b-a))*(n/(b-a));
    double k2 = (m/(d-c))*(m/(d-c));
    double a2 = -2*(h2+k2);
    double alpha = 0.0;
    double beta = 0.0;
    double epsMax = 0, epsCur = 0;
    unsigned int it = 0;
    /* невязку будем считать по
     * внутренней сетке, аналогично и с p, A*p */
    auto r = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    /* считаем r0, h0 */
    for (auto i = 1; i < r.size()-1; i++)
    {
        for (auto j = 1; j < r[i].size()-1; j++)
        {
            r[i][j] = -f[i-1][j-1];
        }
    }
    auto h = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    for (auto j = 0; j < h.size(); j++)
    {
        for (auto i = 0; i < h[j].size(); i++)
        {
            h[i][j] = -r[i][j];
        }
    }

    /* считаем скалярное невязки на h0 */
    double temp = 0.0;
    for (auto j = 0; j < r.size(); j++)
    {
        for (auto i = 0; i < r[j].size(); i++)
        {
            temp += r[i][j] * h[i][j];
        }
    }
    /*считаем A * h0 */
    auto Ah = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    for (auto j = 0; j < h.size(); j++)
    {
        for (auto i = 0; i < h[j].size(); i++)
        {
            if (i >= 1 && i <= n-1 && j >= 1 && j <= m-1)
            {
                Ah[i][j] = a2 * h[i][j] + h2*(h[i+1][j] + h[i-1][j]) +
                        k2*(h[i][j-1] + h[i][j+1]);
            }

            else
            {
                Ah[i][j] = h[i][j];
            }
        }
    }
    /*скалярное A*h0, h0 */
    double temp2 = 0.0;
    for (auto j = 0; j < h.size(); j++)
    {
        for (auto i = 0; i < h[j].size(); i++)
        {
            temp2 += Ah[i][j] * h[i][j];
        }
    }
    alpha = - temp / temp2;
    /* считаем x1 */
    for (auto j = 1; j < m; j++)
    {
        for (auto i = 1; i < n; i++)
        {
            v[i][j] = v[i][j] + alpha * h[i][j];
        }
    }
    while (true)
    {
        /*
         * 1) пересчитать r1
         * 2) Скалярно Ah0, r1
         * 3) beta = поделить Ah0, r1 на Ah0, h0
         * 4) пересчитать h1
         * 5) посчитать Ah1
         * 6) скалярное Ah1, h1
         * 7) скалярное r1, h1
         * 8) посчитать alpha1
         * 9) посчитать x2
         */
        residual = 0.0;
        for (auto j = 0; j < v.size(); j++)
        {
            for (auto i = 0; i < v[j].size(); i++)
            {
                if (i >= 1 && i <= n-1 && j >= 1 && j <= m-1)
                {
                    r[i][j] = a2 * v[i][j] + h2*(v[i+1][j] + v[i-1][j]) +
                            k2*(v[i][j-1] + v[i][j+1]) - f[i-1][j-1];
                    if (abs(r[i][j]) > residual)
                    {
                        residual = abs(r[i][j]);
                    }
                }

                else
                {
                    r[i][j] = 0;
                }
            }
        }
        std::cout << residual << std::endl;
        double Ah0r1 = 0.0;
        for (auto j = 0; j < r.size(); j++)
        {
            for (auto i = 0; i < r[j].size(); i++)
            {
                Ah0r1 += Ah[i][j] * r[i][j];
            }
        }
        double Ah0h0 = 0.0;
        for (auto j = 0; j < h.size(); j++)
        {
            for (auto i = 0; i < h[j].size(); i++)
            {
                Ah0h0 += Ah[i][j] * h[i][j];
            }
        }

        beta = Ah0r1 / Ah0h0;

        for (auto j = 0; j < h.size(); j++)
        {
            for (auto i = 0; i < h[j].size(); i++)
            {
                h[i][j] = -r[i][j] + beta * h[i][j];
            }
        }

        double r1h1 = 0.0;
        for (auto j = 0; j < h.size(); j++)
        {
            for (auto i = 0; i < h[j].size(); i++)
            {
                r1h1 += r[i][j]* h[i][j];
            }
        }


        for (auto j = 0; j < h.size(); j++)
        {
            for (auto i = 0; i < h[j].size(); i++)
            {
                if (i >= 1 && i <= n-1 && j >= 1 && j <= m-1)
                {
                    Ah[i][j] = a2 * h[i][j] + h2*(h[i+1][j] + h[i-1][j]) +
                            k2*(h[i][j-1] + h[i][j+1]);
                }

                else
                {
                    Ah[i][j] = h[i][j];
                }
            }
        }

        double Ah1h1 = 0.0;
        for (auto j = 0; j < r.size(); j++)
        {
            for (auto i = 0; i < r[j].size(); i++)
            {
                Ah1h1 += Ah[i][j] * h[i][j];
            }
        }
        alpha = -r1h1 / Ah1h1;


        for (auto j = 1; j < m; j++)
        {
            for (auto i = 1; i < n; i++)
            {
                v[i][j] = v[i][j] + alpha * h[i][j];
            }
        }
        it += 1;
        if (it >= N || residual < eps)
        {
            break;
        }
//        epsMax = 0.0;
//        scalarMultipleResidual = scalarMultipleResidualNext =
//                scalarMultipleMatrixOnParameter = 0.0;
//        for (auto i = 0; i < r.size(); i++)
//        {
//            for (auto j = 0; j < r[i].size(); j++)
//            {
//            scalarMultipleResidual += r[i][j] * r[i][j];
//            }
//        }
//        for (auto j = 0; j < Apk.size(); j++)
//        {
//            for (auto i = 0; i < Apk[j].size(); i++)
//            {
//                auto temp = 0.0;
//                if (i >= 1 && i <= n-1 && j >= 1 && j <= m-1)
//                {
//                    temp = a2 * p[i][j] + h2*(p[i+1][j] + p[i-1][j]) +
//                            k2*(p[i][j-1] + p[i][j+1]);
//                }

//                else
//                {
//                    temp = p[i][j];
//                }
//                Apk[i][j] = temp;
//            }
//        }
//        for (auto i = 0; i < p.size(); i++)
//        {
//            for (auto j = 0; j < p[i].size(); j++)
//            {
//                scalarMultipleMatrixOnParameter += Apk[i][j] * p[i][j];
//            }
//        }
//        alpha = scalarMultipleResidual / scalarMultipleMatrixOnParameter;
//        for (auto j = 1; j < m; j++)
//        {
//            for (auto i = 1; i < n; i++)
//            {
//                auto vOld = v[i][j];
//                auto vNew = v[i][j] + alpha * p[i][j];
//                epsCur = abs(vOld - vNew);
//                if (epsCur > epsMax)
//                {
//                    epsMax = epsCur;
//                }
//                v[i][j] = vNew;
//            }
//        }
//        for (auto i = 0; i < r.size(); i++)
//        {
//            for (auto j = 0; j < r[i].size(); j++)
//            {
//                rNext[i][j] = r[i][j] - alpha * Apk[i][j];
//            }
//        }
//        for (auto i = 0; i < r.size(); i++)
//        {
//            for (auto j = 0; j < r[i].size(); j++)
//            {
//                scalarMultipleResidualNext += rNext[i][j] * rNext[i][j];
//            }
//        }
//        beta = scalarMultipleResidualNext / scalarMultipleResidual;

//        for (auto i = 0; i < p.size(); i++)
//        {
//            for (auto j = 0; j < p[i].size(); j++)
//            {

//                p[i][j] = rNext[i][j] + beta * p[i][j];
//            }


    }
    QMessageBox msgBox;
    msgBox.setInformativeText(QString("При решении Р.С. с помощью метода Зейделя с параметрами NMax = ") +
                              QString::number(N) + QString("\nи eps = ") +
                              QString::number(eps) + QString(" за S = ") +
                              QString::number(it) + QString(" итераций было получено решение \nс точностью eps max = ") +
                              QString::number(epsMax) + QString(" Невязка составила: ") + QString::number(residual));
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
    double a = 1, b = 2, c = 2, d = 3;

    // Utest = sin(pi*x*y)

    unsigned int n = ui->spinBox->value();
    unsigned int m = ui->spinBox_2->value();
    unsigned int limit = ui->spinBox_3->value();
    double eps = ui->textEdit->toPlainText().toDouble();
    double h = double(b-a) / double(n);
    double k = double(d-c) / double(m);
    auto startSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    auto f = std::vector<std::vector<double>>(n - 1, std::vector<double> (m - 1, 0));
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m -1; ++j)
        {
            f[i][j] = - f_test(a + (i + 1) * h, c + (j + 1) * k);
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
    auto solution = startSolution;
    if (ui->comboBox->currentText() == "Метод Зейделя")
    {
        solution = solve_Zeidel(startSolution, f, n, m, limit, eps);
    }
    else if (ui->comboBox->currentText() == "Метод сопряженных градиентов")
    {
        solution = solve_Conjugate_Gradient(startSolution, f, n, m, limit, eps);
    }
    auto accurateSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            accurateSolution[i][j] = acc_u(a + i * h, c + j * k);
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
            double currentElem = abs(solution[i][j] - accurateSolution[i][j]);
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

                              QString::number(1 + double(maxXDiff) * 2 * h) +
                              QString(", y = ") +
                              QString::number(2 + double(maxYDiff) * 2 * k) );
    msgBox.exec();

}
