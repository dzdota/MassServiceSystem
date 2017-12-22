using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace testgistogr
{
    class UniformityDoubleMultidimensionalData
    {
        InitialAnalysMultidimensionalData X;
        InitialAnalysMultidimensionalData Y;

        public List<Data> Estimation = new List<Data>();

        public double N1, N2;
        public int n;

        public Data SravnSred;

        public UniformityDoubleMultidimensionalData(List<InitialStatisticalAnalys> ISA,List<List<int>> ind) 
        {
            X = new InitialAnalysMultidimensionalData(ISA,ind[0],0);
            Y = new InitialAnalysMultidimensionalData(ISA,ind[1],0);


            n = ind[0].Count;
            N1 = ISA[ind[0][0]].unsortl.Length;
            N2 = ISA[ind[1][0]].unsortl.Length;
            for (int i = 1; i < n;i++) 
                if (N1 != ISA[ind[0][i]].unsortl.Length || N2 != ISA[ind[0][i]].unsortl.Length) 
                {
                    MessageBox.Show("Not correct input data", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }
            Refresh();

        }
        public void  Refresh()
        {
            Estimation.Clear();
            SravnSred = SravnSredF(X, Y);
            Estimation.Add(SravnSred);
        }
        public Data SravnSredF(InitialAnalysMultidimensionalData X, InitialAnalysMultidimensionalData Y)
        {
            Data rez = new Data();

            double[,] S0 = S0F(X, Y);
            double[,] S1 = S1F(X, Y);
            rez.Name = "Рівність багатомірних середніх при рівності DC матриці";
            rez.Q = -(N1 + N2 - 2 - n / 2) * Math.Log(Matrix.Determinant(S1) / Matrix.Determinant(S0));
            rez.QKvant = Hi.HIF(X.ISA[0].alf.Q, n);
            rez.H = rez.Q <= rez.QKvant;
            return rez;
        }
        private double[,] S0F(InitialAnalysMultidimensionalData X, InitialAnalysMultidimensionalData Y)
        {
            int n = X.ISA.Count;
            double[,] rez = new double[n,n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    rez[i, j] = S0ij(X, Y, i, j, X.ISA[0].unsortl.Length, Y.ISA[0].unsortl.Length);
            return rez;
        }

        private double S0ij(InitialAnalysMultidimensionalData X, InitialAnalysMultidimensionalData Y, int i, int j, double N1, double N2) 
        {
            double rez = 0;
            double xij = 0,
                yij = 0,
                xi = 0,
                xj = 0,
                yi = 0,
                yj = 0;
            xi = X.ISA[i].unsortl.Sum();
            yi = Y.ISA[i].unsortl.Sum();
            xj = X.ISA[j].unsortl.Sum();
            yj = Y.ISA[j].unsortl.Sum();
            for (int L = 0; L < N1; L++) 
                xij += X.ISA[i].unsortl[L] * X.ISA[j].unsortl[L];
            for (int L = 0; L < N2; L++)
                yij += Y.ISA[i].unsortl[L] * Y.ISA[j].unsortl[L];
            rez = -(xi + yi) * (xj + yj) / (N1 + N2);
            rez += xij + yij;
            rez /= N1 + N2 - 2;
            return rez;
        }

        private double[,] S1F(InitialAnalysMultidimensionalData X, InitialAnalysMultidimensionalData Y)
        {
            int n = X.ISA.Count;
            double[,] rez = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    rez[i, j] = S1ij(X, Y, i, j, X.ISA[0].unsortl.Length, Y.ISA[0].unsortl.Length);
            return rez;

        }

        private double S1ij(InitialAnalysMultidimensionalData X, InitialAnalysMultidimensionalData Y, int i, int j, double N1, double N2)
        {
            double rez = 0;
            double xij = 0,
                yij = 0,
                xi = 0,
                xj = 0,
                yi = 0,
                yj = 0;
            xi = X.ISA[i].unsortl.Sum();
            yi = Y.ISA[i].unsortl.Sum();
            xj = X.ISA[j].unsortl.Sum();
            yj = Y.ISA[j].unsortl.Sum();
            for (int L = 0; L < N1; L++)
                xij += X.ISA[i].unsortl[L] * X.ISA[j].unsortl[L];
            for (int L = 0; L < N2; L++)
                yij += Y.ISA[i].unsortl[L] * Y.ISA[j].unsortl[L];
            rez = -(xi * xj) / (N1);
            rez += -(yi * yj) / (N2);
            rez += xij + yij;
            rez /= N1 + N2 - 2;
            return rez;
        }

    }
}
