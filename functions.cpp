#include <iostream>
#include "functions.h"
#include "eigen.h"

void test()
{
    std::cout << "test ~~" << std::endl;
}

void cov(double** doubleMatrix, int row, int col)
{
    double meanValue, total;
    double *mean = new double [col];
    double **tmp = newDoubleMatrix(row, col);
    double **covMatrix = newDoubleMatrix(col, col);
    copyMatrix(doubleMatrix, tmp, row, col);

    for(int i = 0; i < col; i++)
    {
        meanValue = 0;
        for (int j = 0; j < row; j++)
        {
            meanValue += tmp[j][i];
        }
        meanValue /= row;
        mean[i] = meanValue;
    }
    for(int i = 0; i < col; i++)
    {
        for (int j = 0; j < row; j++)
        {
            tmp[j][i] -= mean[i];
        }
    }

    for(int i = 0; i < col; i++)
    {
        for (int j = 0; j < col; j++)
        {
            total = 0;
            for (int k = 0; k < row; k++)
            {
                total += tmp[k][i] * tmp[k][j];
            }
            total /= (row - 1);
            covMatrix[i][j] = total;
        }
    }

    std::cout << "covariance matrix:\n";
    for(int i = 0; i < col; i++)
    {
        for (int j = 0; j < col; j++)
        {
            std::cout << covMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    

    eigTest(covMatrix);



    delete [] mean;
    deleteDoubleMatrix(tmp, row);
}

void copyMatrix(double** source, double** destination, int row, int col)
{
    for(int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            destination[i][j] = source[i][j];
        }
    }
}
void deleteDoubleMatrix(double** matrix, int row)
{
    for(int i = 0; i < row; i++)
    {
        delete [] matrix[i];
    }
    delete [] matrix;
}

double** newDoubleMatrix(int row, int col)
{
    double** tmp = new double*[row];
    for (int i = 0; i < row; i++)
    {
        tmp[i] = new double[col];
        for (int j = 0; j < col; j++)
        {
            tmp[i][j] = 0;
        }
    }
    return tmp;
}
