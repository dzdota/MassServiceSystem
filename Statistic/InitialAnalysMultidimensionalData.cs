using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace testgistogr
{
    class InitialAnalysMultidimensionalData
    {
        public List<InitialStatisticalAnalys> ISA;
        public int n;
        public int N;
        public double[] E;
        public double[,] X;
        public double[,] Y;
        public double[,] DC;
        public double[,] K;
        public Data[] CMC;
        public Data CoefDeter;
        public Data ZnachRegres;
        public bool[] KTurnOn;
        public int k;
        public Data[] A;
        public Data[] Astandart;
        public double[,] C;
        public double Szal;
        public double SigmKv;

        public InitialAnalysMultidimensionalData(List<InitialStatisticalAnalys> ISA, List<int> IND, int k)
        {
            this.k = k;
            this.ISA = new List<InitialStatisticalAnalys>();
            E = new double[IND.Count];
            for (int i = 0; i < IND.Count; i++)
            {
                this.ISA.Add(ISA[IND[i]]);
                E[i] = ISA[IND[i]].Mx.Q;
            }
            n = this.ISA.Count;
            N = this.ISA[0].unsortl.Length;
            DC = DCFound(ISA, IND);
            K = KFound(ISA, IND);
            X = XFound(this.ISA, n, N, k);
            Y = XFound(this.ISA, n, N, k);
            KTurnOn = new bool[n];
            for (int i = 0; i < KTurnOn.Length; i++)
                KTurnOn[i] = true;
            if (k>=0)
            {
                CMC = new Data[n];
                for (int i =0;i<n;i++)
                {
                    CMC[i] = new Data();
                    CMC[i].Name = "Коеф. багатов. корел.";
                    CMC[i].Q = CMCFound(ISA, K, i);
                    CMC[i].QKvant = (N - n - 1) * (Math.Pow(CMC[i].Q, 2)) / (n * (1 - Math.Pow(CMC[i].Q, 2)));
                    CMC[i].Q0 = Distributions.FisherQuantile(1 - ISA[0].alf.Q / 2, n, N - n - 1);
                }

                SigmKv = Math.Sqrt(this.ISA[k].Dx.Q * (1 - Math.Pow(CMC[k].Q, 2)) * (N - 1) / (N - n-1));

                CoefDeter = new Data();
                CoefDeter.Name = "Коефіціент Детермінації = ";
                CoefDeter.Q = Math.Pow(CMC[k].Q, 2);


                ZnachRegres = new Data();
                ZnachRegres.Name = "Перевірка значущості від. регресії";
                ZnachRegres.QKvant = (N - n ) * ((1 / (1 - CoefDeter.Q)) - 1) / (n - 1) ;
                ZnachRegres.Q0 = Distributions.FisherQuantile(1 - ISA[0].alf.Q / 2, n-1, N - n - 2);

                C = CFound(n,X);
                A = MuliplRegresFound(this.ISA, k);
                Astandart = AstndF(this.ISA, A, k);
                Szal = SzalFound(this.ISA, k, A);

            }
        }
        Data[] AstndF(List<InitialStatisticalAnalys> ISA, Data[] A,int k)
        {
            Data[] As = new Data[A.Length - 1];
            for (int i = 0; i < As.Length; i++)
                As[i] = new Data();
            for (int i =0;i<ISA.Count-1;i++)
            {
                int ireal = i;
                if (i >= k)
                    ireal += 1;
                As[i].Q = A[i+1].Q * ISA[ireal].Gx.Q / ISA[k].Gx.Q;

            }
            return As;
        }
        double[,] XFound(List<InitialStatisticalAnalys> ISA, int n, int N, int k)
        {

            double[,] X = new double[N, n - 1];
            for (int row = 0; row < N; row++)
                for (int col = 0; col < n - 1; col++)
                {
                    int realcol = col;
                    if (col >= k)
                        realcol++;
                    X[row, col] = ISA[realcol].unsortl[row];
                }
            return X;
        }
        double[,] YFound(List<InitialStatisticalAnalys> ISA, int n, int N, int k)
        {

            double[,] Y = new double[N, 1];
            for (int row = 0; row < N; row++)
                    Y[row, 0] = ISA[k].unsortl[row];
            return Y;
        }
        double[,] CFound(int n,double[,] X)
        {
            double[,] C = new double[n - 1, n - 1];
            C = Matrix.InverseMatrix(Matrix.MultiplicMatrix(Matrix.TranspMatrix(X),X));
            return C;

        }
        double CMCFound(List<InitialStatisticalAnalys> ISA, double[,] K, int k)
        {
            double rez = 0;

            double detK = Matrix.Determinant(K);
            double detSubK = Matrix.Determinant(Matrix.SubMatrix(K, k, k));
            rez = Math.Sqrt(1 - (detK / detSubK));
            /*
            List<InitialStatisticalAnalys> ISAStand = new List<InitialStatisticalAnalys>();
            double[,] SubK = Matrix.SubMatrix(K, k, k);
            double del = Matrix.Determinant(SubK);
            double[,] delS_Matrix = new double[SubK.GetLength(0)+1, SubK.GetLength(1)+1];
            for (int row = 0; row < K.GetLength(0); row++)
                ISAStand.Add(new InitialStatisticalAnalys(
                    InitialStatisticalAnalys.StandData(ISA[row].unsortl, ISA[row].Gx.Q, ISA[row].Mx.Q))
                    );
            for (int row = 0; row < SubK.GetLength(0); row++)
                for (int col = 0; col < SubK.GetLength(1); col++)
                    delS_Matrix[row + 1, col] = SubK[row, col];
            for (int col = 0; col < SubK.GetLength(1); col++)
            {
                int colreal = col;
                if (col >= k)
                    colreal++;
                delS_Matrix[0, col] = Correlation_RegressionAnalysis.KorelationFound(ISAStand[k], ISAStand[colreal]);
                delS_Matrix[col + 1, delS_Matrix.GetLength(1) - 1] = delS_Matrix[0, col];
            }
            double delS = Matrix.Determinant(delS_Matrix);
            rez = Math.Sqrt(delS / del);*/
            return rez;
        }
        double[,] KFound(List<InitialStatisticalAnalys> ISA, List<int> IND)
        {
            double[,] rez = new double[IND.Count, IND.Count];
            for (int i = 0; i < IND.Count; i++)
                rez[i, i] = 1;
            for (int i = 0; i < IND.Count; i++)
            {
                for (int j = i + 1; j < IND.Count; j++)
                {
                    rez[i, j] = Correlation_RegressionAnalysis.KorelationFound(ISA[IND[i]], ISA[IND[j]]);
                    rez[j, i] = rez[i, j];
                }
            }
            return rez;
        }
        double[,] DCFound(List<InitialStatisticalAnalys> ISA, List<int> IND)
        {
            double[,] rez = new double[IND.Count, IND.Count];
            for (int i = 0; i < IND.Count; i++)
                rez[i, i] = ISA[IND[i]].Dx.Q;
            for (int i = 0; i < IND.Count; i++)
            {
                for (int j = i + 1; j < IND.Count; j++)
                {
                    rez[i, j] = cov(ISA[IND[i]], ISA[IND[j]]);
                    rez[j, i] = rez[i, j];
                }
            }
            return rez;
        }
        double cov(InitialStatisticalAnalys ISAe1, InitialStatisticalAnalys ISAe2)
        {
            return ISAe1.Gx.Q * ISAe2.Gx.Q * Correlation_RegressionAnalysis.KorelationFound(ISAe1, ISAe2);
        }
        public static Data[] MuliplRegresFound(List<InitialStatisticalAnalys> ISA, int k)
        {
            int N = ISA[0].unsortl.Length;
            int n = ISA.Count;
            Data[] A = new Data[ISA.Count];
            for (int i =0; i< A.Length;i++)
                A[i] = new Data();
            double[,] X_EX = new double[N,n-1];
            double[,] Y0= new double[N, 1];
            for (int row = 0; row < N; row++)
            {
                Y0[row, 0] = ISA[k].unsortl[row] - ISA[k].Mx.Q;
                for (int col = 0; col < n-1; col++)
                {
                    int colreal = col;
                    if (col >= k)
                        colreal++;
                    X_EX[row, col] = ISA[colreal].unsortl[row] - ISA[colreal].Mx.Q;
                }
            }
            double[,] X_EXmX_ExTr = Matrix.MultiplicMatrix(Matrix.TranspMatrix(X_EX), X_EX);
            double[,] InversX_EXmX_ExTr = Matrix.InverseMatrix(X_EXmX_ExTr);
            double[,] G = Matrix.MultiplicMatrix(Matrix.TranspMatrix(X_EX), Y0);
            double[,] Ai = Matrix.MultiplicMatrix(InversX_EXmX_ExTr, G);
            double AkExk = 0;
            for (int i = 1; i < n; i++)
            {
                A[i].Q = Ai[i - 1, 0];
                int ireal = i;
                if (i-1 >= k)
                    ireal++;
                AkExk += A[i].Q * ISA[ireal - 1].Mx.Q;
            }
            A[0].Q = ISA[k].Mx.Q - AkExk;

            return A;
        }
        public static double SzalFound(List<InitialStatisticalAnalys> ISA, int k,Data[] A)
        {
            int N = ISA[0].unsortl.Length;
            int n = ISA.Count;
            double Szal = 0;
            for (int i =0;i<N;i++)
            {
                double AX = 0;
                for (int j =0;j < n-1;j++)
                {
                    int jreal = j;
                    if (jreal >= k)
                        jreal++;
                    AX += A[j + 1].Q * ISA[jreal].unsortl[i];
                }
                Szal += Math.Pow(ISA[k].unsortl[i] - A[0].Q - AX, 2);
            }
            Szal /= N - n;
            Szal = Math.Sqrt(Szal);
            return Szal;
        }
        public static double[,] SubK(double[,] K, bool[] KTurnOn)
        {
            double[,] SubK = K;
            for (int i = KTurnOn.Length - 1; i >= 0; i--) 
                if(!KTurnOn[i])
                    SubK = SubK1(SubK, i);
            return SubK;

        }
        public static double[,] SubK1(double[,] K, int index)
        {
            double[,] SubK = new double[K.GetLength(0), K.GetLength(1)];
            for (int col = 0;col <SubK.GetLength(1);col++)
            {
                if (col == index)
                    continue;
                for (int row = 0; row < SubK.GetLength(0); row++)
                {
                    if (row == index)
                        continue;
                    SubK[row, col] = (K[row, col] - K[row, index] * K[col, index]) /
                        Math.Sqrt((1 - Math.Pow(K[row, index], 2)) * (1 - Math.Pow(K[col, index], 2)));
                }
            }
            return SubK;
        }
    }
}
