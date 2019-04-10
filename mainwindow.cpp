#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

double norm(std::vector<double> one, std::vector<double> two)
{
    double res = 0;
    for (auto c1 = one.cbegin(); c1 != one.cend(); c1++)
    {
        for (auto c2 = two.cbegin(); c1 != two.cend(); c1++)
        {
            res += pow(((*c1) * (*c1) - (*c2) * (*c2)), 2);
        }
    }
    return sqrt(res);
}

std::vector<double> solve_Zeidel(std::vector<double> &startSolution,
                                 std::vector<double> &f, unsigned int n, unsigned int m,
                                 unsigned long N, double eps)
{
    auto xNew = startSolution;
    auto xOld = startSolution;
    unsigned long k = 0;
    while ((k < N) || (norm(xNew, xOld) < eps))
    {
        for (unsigned int i = 0; i < ((n-1)*(m-1)); i++)
        {
            double val = 0;
            for(unsigned int j = 0; j < i; j++)
            {
                if (j == )
            }
        }
    }


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

}

void MainWindow::on_pushButton_2_clicked()
{
    //основная задача
     // f = -e^(-xy^2)

}
