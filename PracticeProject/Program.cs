using System;
using MathNet.Numerics.LinearAlgebra;

namespace Appstarter{
    public class Appstarter{
        static void Main(string[] args)
        {
            Console.WriteLine("씨밤 세팅하다 날 끝나네");

            Matrix<double> A = Matrix<double>.Build.DenseOfArray(new double[,]
        {
            {2.0, 1.0, -1.0},
            {-3.0, -1.0, 2.0},
            {-2.0, 1.0, 2.0}
        });

        // 3x1 벡터 B 생성
        Vector<double> B = Vector<double>.Build.DenseOfArray(new double[] {8.0, -11.0, -3.0});

        GaussJordan(ref A, ref B);

        Console.WriteLine(B);
        }

        public static bool Solve(float[,] M)
        {
            // input checks
            int rowCount = M.GetUpperBound(0) + 1;
            if (M == null || M.Length != rowCount * (rowCount + 1))
            throw new ArgumentException("The algorithm must be provided with a (n x n+1) matrix.");
            if (rowCount < 1)
            throw new ArgumentException("The matrix must at least have one row.");

            // pivoting
            for (int col = 0; col + 1 < rowCount; col++) if (M[col, col] == 0)
            // check for zero coefficients
            {
                // find non-zero coefficient
                int swapRow = col + 1;
                for (;swapRow < rowCount; swapRow++) if (M[swapRow, col] != 0) break;

                if (M[swapRow, col] != 0) // found a non-zero coefficient?
                {
                    // yes, then swap it with the above
                    float[] tmp = new float[rowCount + 1];
                    for (int i = 0; i < rowCount + 1; i++)
                    { tmp[i] = M[swapRow, i]; M[swapRow, i] = M[col, i]; M[col, i] = tmp[i]; }
                }
                else return false; // no, then the matrix has no unique solution
            }

            // elimination
            for (int sourceRow = 0; sourceRow + 1 < rowCount; sourceRow++)
            {
                for (int destRow = sourceRow + 1; destRow < rowCount; destRow++)
                {
                    float df = M[sourceRow, sourceRow];
                    float sf = M[destRow, sourceRow];
                    for (int i = 0; i < rowCount + 1; i++)
                    M[destRow, i] = M[destRow, i] * df - M[sourceRow, i] * sf;
                }
            }

            // back-insertion
            for (int row = rowCount - 1; row >= 0; row--)
            {
                float f = M[row,row];
                if (f == 0) return false;

                for (int i = 0; i < rowCount + 1; i++) M[row, i] /= f;
                for (int destRow = 0; destRow < row; destRow++)
                { M[destRow, rowCount] -= M[destRow, row] * M[row, rowCount]; M[destRow, row] = 0; }
            }
            return true;
        }


        public static void GaussJordan(ref Matrix<double> A, ref Vector<double> B)
        {
            int rowCount = A.RowCount;

            // Input checks
            if (A == null || B == null || A.RowCount != rowCount || A.ColumnCount != rowCount || B.Count != rowCount)
            {
                throw new ArgumentException("Invalid input. A must be a square matrix, and B must be a vector of the same size as A.");
            }

            // Pivoting
            for (int col = 0; col < rowCount; col++)
            {
                if (A[col, col] == 0)
                {
                    // Find a non-zero coefficient in the same column
                    int swapRow = -1;
                    for (int row = col + 1; row < rowCount; row++)
                    {
                        if (A[row, col] != 0)
                        {
                            swapRow = row;
                            break;
                        }
                    }

                    if (swapRow != -1) // Found a non-zero coefficient
                    {
                        for (int i = 0; i < rowCount; i++)
                        {
                            double temp = A[col, i];
                            A[col, i] = A[swapRow, i];
                            A[swapRow, i] = temp;
                        }
                        // Swap elements in vector B
                        double tempB = B[col];
                        B[col] = B[swapRow];
                        B[swapRow] = tempB;
                    }
                    else
                    {
                        throw new ArgumentException("The matrix does not have a unique solution.");
                    }
                }
            }

            // Elimination
            for (int sourceRow = 0; sourceRow < rowCount; sourceRow++)
            {
                double pivot = A[sourceRow, sourceRow];

                for (int i = 0; i < rowCount; i++)
                {
                    A[sourceRow, i] /= pivot;
                }

                B[sourceRow] /= pivot;

                for (int destRow = 0; destRow < rowCount; destRow++)
                {
                    if (destRow != sourceRow)
                    {
                        double factor = A[destRow, sourceRow];
                        for (int i = 0; i < rowCount; i++)
                        {
                            A[destRow, i] -= factor * A[sourceRow, i];
                        }
                        B[destRow] -= factor * B[sourceRow];
                    }
                }
            }
        }


    }
}