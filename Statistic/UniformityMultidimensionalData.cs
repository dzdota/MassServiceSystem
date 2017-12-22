using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace testgistogr
{
    class UniformityMultidimensionalData
    {
        List<InitialAnalysMultidimensionalData> IAMD;

        public List<Data> Estimation = new List<Data>();

        Data SravSred;
        Data SravnDisper;
        int n;
        double[] N;
        public UniformityMultidimensionalData(List<InitialStatisticalAnalys> ISA,List<List<int>> ind) 
        {
            IAMD = new List<InitialAnalysMultidimensionalData>();
            n = ind[0].Count;
            N = new double[ind.Count];
            for(int i =0;i<ind.Count;i++)
                IAMD.Add(new InitialAnalysMultidimensionalData(ISA,ind[i],0));

            for(int d = 0;d<ind.Count;d++)
            {
                N[d] = IAMD[d].ISA[0].unsortl.Length;
                for (int i = 1; i < IAMD[d].ISA.Count;i++)
                    if (n != IAMD[d].ISA.Count||N[d] != IAMD[d].ISA[i].unsortl.Length) 
                    {
                        MessageBox.Show("Not correct input data", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                        return;
                    }
            }
            Referesh();

        }
        public void Referesh() 
        {
            Estimation.Clear();
            double[,] Ex = ExF(IAMD);
            SravSred = SravSredF(IAMD);
            SravnDisper = SravnDisperF(IAMD);
            Estimation.Add(SravSred);
            Estimation.Add(SravnDisper);
            
        }

        private Data SravnDisperF(List<InitialAnalysMultidimensionalData> IAMD)
        {
            Data rez = new Data();
            rez.Name = "Збіг DC";
            double Nsum = N.Sum();
            double V = 0;
            double[,] S = new double[n, n];
            List<double[,]> Sd = new List<double[,]>();
            for (int d = 0; d < IAMD.Count;d++ )
                Sd.Add(SdF(IAMD[d]));
            for (int d = 0; d < IAMD.Count; d++)
                S = Matrix.Addition(S,Matrix.MultiplicNumber(Sd[d],N[d]-1));
            S = Matrix.MultiplicNumber(S,1.0/(Nsum - IAMD.Count));

            double detS = Matrix.Determinant(S);
            for (int d = 0; d < IAMD.Count; d++) 
            {
                V += ((N[d] - 1) / 2) * Math.Log(detS / Matrix.Determinant(Sd[d]));
            }
            rez.Q = V;
            rez.QKvant = Hi.HIF(IAMD[0].ISA[0].alf.Q, n * (n + 1) * (IAMD.Count - 1) / 2); 
            rez.H = rez.Q <= rez.QKvant;
            return rez;
        }
        private Data SravSredF(List<InitialAnalysMultidimensionalData> IAMD) 
        {
            Data rez = new Data();
            rez.Name = "Збіг середніх при розбіжності DC";
            double[,] Ex = ExF(IAMD);
            double[,] V = new double[1, 1];
            for (int d = 0; d < IAMD.Count; d++) 
            {
                double[,] Exd = xdF(IAMD[d]);
                double[,] Xdx = new double[1, n];
                for (int i = 0; i < n;i++ )
                {
                    Xdx[0, i] = Exd[0, i] - Ex[0, i];
                }
                double[,] SInv = Matrix.InverseMatrix(SdF(IAMD[d]));
                V = Matrix.Addition(V, Matrix.MultiplicMatrix(Xdx, Matrix.MultiplicMatrix(SInv, Matrix.TranspMatrix(Xdx))));
            }
            rez.Q = V[0, 0];
            rez.QKvant = Hi.HIF(IAMD[0].ISA[0].alf.Q, (int)(n*(IAMD.Count - 1)));
            rez.H = rez.Q <= rez.QKvant;
            return rez;
        }
        private double[,] ExF(List<InitialAnalysMultidimensionalData> IAMD) 
        {
            double[,] Ex;
            double[,] NSd = new double[IAMD[0].n, IAMD[0].n];
            double[,] NSObrExTd = new double[IAMD[0].n, 1];
            for (int d = 0; d < IAMD.Count;d++ )
            {
                double[,] SObr = Matrix.InverseMatrix(SdF(IAMD[d]));
                NSd = Matrix.Addition(NSd,Matrix.MultiplicNumber(SObr,IAMD[d].N));


                double[,] ExdT = Matrix.TranspMatrix(xdF(IAMD[d]));
                double[,] SObrExdT = Matrix.MultiplicMatrix(SObr, ExdT);
                NSObrExTd = Matrix.Addition(NSObrExTd,Matrix.MultiplicNumber(SObrExdT,IAMD[d].N));
            }
            Ex = Matrix.TranspMatrix(Matrix.MultiplicMatrix(Matrix.InverseMatrix(NSd),NSObrExTd));
            return Ex;
        }
        private double[,] xdF(InitialAnalysMultidimensionalData IAMD) 
        {
            double[,] rez = new double[1, IAMD.n];
            for (int i = 0; i < IAMD.E.Length; i++)
                rez[0, i] = IAMD.E[i];
            return rez;
        }
        private double[,] XdF(List<InitialAnalysMultidimensionalData> IAMD, int d,int l)
        {
            double[,] rez = new double[1, IAMD[d].E.Length];
            for (int i = 0; i < IAMD[d].E.Length; i++)
                rez[l, i] = IAMD[d].ISA[i].unsortl[l];
            return rez;
        }

        private double[,] SdF(InitialAnalysMultidimensionalData IAMD)
        {
            double[,] rez ;
            double[,] X = new double[IAMD.N, IAMD.n];
            double[,] XT;
            double[,] Ex = xdF(IAMD);
            for (int i = 0; i < IAMD.n; i++)
                for (int l = 0; l < IAMD.N;l++ )
                    X[l, i] = IAMD.ISA[i].unsortl[l] - Ex[0, i];
            XT = Matrix.TranspMatrix(X);
            rez = Matrix.MultiplicMatrix(XT, X);
            rez = Matrix.MultiplicNumber(rez, 1.0 / (IAMD.N - 1));
            return rez;
        }
    }
}
