using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Polymer_brush
{
    static class Matrix
    {
        public static double determinantGauss(int n, double[,] M, int currentColumn, bool changeSign)
        {
            double[,] A = CopyMatrix(n, M);
            //Calculate determinant if there are no rows
            if (currentColumn == n)
            {
                double det = 1;
                for (int i = 0; i < n; i++)
                    det *= A[i,i];
                if (changeSign)
                    det *= -1;
                return det;
            }
            //Find first non-zero element in a current column
            int nonZeroRow = -1;
            for (int i = currentColumn; i < n; i++)
            {
                if (A[i,currentColumn] != 0)
                {
                    nonZeroRow = i;
                    break;
                }
            }
            if (nonZeroRow == -1)
                return 0;
            //Change strings, change sign
            if (nonZeroRow != currentColumn)
            {
                double nonZeroRowCopy = 0;
                for (int i = 0; i < n; i++)
                {
                    nonZeroRowCopy = A[nonZeroRow,i];
                    A[nonZeroRow,i] = A[currentColumn,i];
                    A[currentColumn,i] = nonZeroRowCopy;
                }
                //
                if (changeSign)
                    changeSign = false;
                else
                    changeSign = true;
            }
            //Substract first column from other columns
            for (int i = currentColumn + 1; i < n; i++)
            {
                if (A[i,currentColumn] != 0)
                {
                    //Substract rows
                    double coeffitient = A[i,currentColumn] / A[currentColumn,currentColumn];
                    for (int j = 0; j < n; j++)
                    {
                        A[i,j] -= A[currentColumn,j] * coeffitient;
                    }
                }
            }
            //Recursion
            currentColumn++;
            return determinantGauss(n, A, currentColumn, changeSign);
        }
        public static double[,] CopyMatrix(int n, double[,] A)
        {
            double[,] copy = new double[n,n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    copy[i,j] = A[i,j];
            }
            return copy;
        }

        //Функции матрицчной алгебры
        public static double[,] unitMatrix(int n)
        {
            double[,] matrix = new double[n,n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i,j] = 0;
                    if (i == j)
                        matrix[i,j] = 1;
                }
            }
            return matrix;
        }
        public static double[,] multiplyMatrixes(int ai, int bj, int m, double[,] A, double[,] B)
        {
            double[,] C = new double[ai, bj];
            for (int i = 0; i < ai; i++)
            {
                for (int j = 0; j < bj; j++)
                {
                    double c = 0;
                    for (int k = 0; k < m; k++)
                        c += A[i,k] * B[k,j];
                    C[i,j] = c;
                }
            }
            return C;
        }
        public static double determinant(int n, double[,] A, int[] expeledColumns, int[] expeledRows, int expeledRowsAndColumns)
        {
            if (expeledColumns == null)
            {
                expeledColumns = new int[n];
                for (int i = 0; i < n; i++)
                    expeledColumns[i] = -1;
            }
            if (expeledRows == null)
            {
                expeledRows = new int[n];
                for (int i = 0; i < n; i++)
                    expeledRows[i] = -1;
            }
            //for (int i = 0; i < n; i++)
            //    std::cout << expeledColumns[i] << "\n";
            //expeledRowsAndColumns++;

            if (expeledRowsAndColumns == n - 1)
            {
                int leftColumn = -1;
                int leftRow = -1;
                for (int i = 0; i < n; i++)
                    if (arrayContains(expeledColumns, i, n) == -1)
                        leftColumn = i;
                for (int i = 0; i < n; i++)
                    if (arrayContains(expeledRows, i, n) == -1)
                        leftRow = i;
                return A[leftRow,leftColumn];
            }

            //Find not expeled columns for next determinant
            int[] notExpeledColumns = new int[n - expeledRowsAndColumns];
            int j = 0;
            for (int i = 0; i < n; i++)
                if (arrayContains(expeledColumns, i, n) == -1)
                {
                    notExpeledColumns[j] = i;
                    j++;
                }
            int[] notExpeledRows = new int[n - expeledRowsAndColumns];
            j = 0;
            for (int i = 0; i < n; i++)
                if (arrayContains(expeledRows, i, n) == -1)
                {
                    notExpeledRows[j] = i;
                    j++;
                }

            double output = 0;


            for (int i = 0; i < n - expeledRowsAndColumns; i++)
            {
                int indexOfExpelingRow = notExpeledRows[0];
                int indexOfExpelingColumn = notExpeledColumns[i];
                expeledColumns[expeledRowsAndColumns] = indexOfExpelingColumn;
                expeledRows[expeledRowsAndColumns] = indexOfExpelingRow;

                expeledRowsAndColumns++;
                double dt = determinant(n, A, expeledColumns, expeledRows, expeledRowsAndColumns);
                expeledRowsAndColumns--;

                if (i % 2 == 0)
                    output += A[indexOfExpelingRow,indexOfExpelingColumn] * dt;
                else
                    output -= A[indexOfExpelingRow,indexOfExpelingColumn] * dt;

            }
            expeledColumns[expeledRowsAndColumns] = -1;
            expeledRows[expeledRowsAndColumns] = -1;
            return output;
        }
        public static int arrayContains(int[] a, int v, int aSize)
        {
            //int n= sizeof(a) / sizeof(int);
            for (int i = 0; i < aSize; i++)
            {
                if (a[i] == v)
                    return i;
            }
            return -1;
        }
        public static double[,]reverseMatrix(int n, double[,] M)
        {
            double[,] reverseMatrix = zeroMatrix(n);
            double[,] A = CopyMatrix(n, M);
            //double detA = determinant(n, A, NULL, NULL, 0);
            double detA = determinantGauss(n, A, 0, false);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    /* int* expeledColumn = new int[n];
                     int* expeledRow = new int[n];
                     for (int k = 0; k < n; k++)
                         expeledColumn[i] = -1;
                     for (int k = 0; k < n; k++)
                         expeledRow[i] = -1;
                     expeledColumn[0] = j;
                     expeledRow[0] = i;

                     double minor = determinant(n, A, expeledColumn, expeledRow, 1);*/
                    double[,] minorMatrix = new double[n-1, n - 1];
                    //for (int k = 0; k < n - 1; k++)
                    //    minorMatrix[k] = new double[n - 1];
                    for (int a = 0; a < n; a++)
                        for (int b = 0; b < n; b++)
                        {
                            if (a < i)
                            {
                                if (b < j)
                                    minorMatrix[a,b] = A[a,b];
                                if (b > j)
                                    minorMatrix[a,b - 1] = A[a,b];
                            }
                            if (a > i)
                            {
                                if (b < j)
                                    minorMatrix[a - 1,b] = A[a,b];
                                if (b > j)
                                    minorMatrix[a - 1,b - 1] = A[a,b];
                            }
                        }
                    ////////////////////////////////////////////////////////////////
                    //double minor = determinantGauss(n, A, expeledColumn, expeledRow, 1);
                    double minor = determinantGauss(n - 1, minorMatrix, 0, false);
                    //double v = (pow(-1, (i + j) % 2)) * minor / detA;
                    double v = minor / detA;
                    if ((i + j) % 2 == 1)
                        v = -v;
                    reverseMatrix[j,i] = v;
                }
            return reverseMatrix;
        }
        public static double[,] zeroMatrix(int n)
        {
            double[,] matrix = new double[n,n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    matrix[i,j] = 0;
                }
            }
            return matrix;
        }
        public static double[] multiplyMatrixAndVector(int ai, int bj, double[,] A, double[] B)
        {
            double[] C = new double[ai];
            for (int i = 0; i < ai; i++)
            {
                C[i] = 0;
                for (int k = 0; k < bj; k++)
                    C[i] += A[i,k] * B[k];
            }
            return C;
        }
    }
}
