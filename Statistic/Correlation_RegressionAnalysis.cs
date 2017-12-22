using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace testgistogr
{
    class Correlation_RegressionAnalysis
    {
        public string RegresTypeVib;

        List<InitialStatisticalAnalys> ML;
        List<int> MSelectGR;
        public bool Doubl = false;
        public bool Nezal = true;
        public double N = 0;
        public double Nmax = 0;
                                               
        public double[] Korelation = new double[4];
        public double[] KorelationVidnoh = new double[3];
        public double[] RangKorelation = new double[2];
        public double[] RangKoefKend = new double[3];

        public double IndexFehnera = 0;
        public double[] KoefSpolF = new double[3];
        public double[] KoefSvazYoola = new double[4];

        public double[] NezXY = new double[2];
        public double[] KoefSpolPirs = new double[2];
        public double[] MeraZvazKendallaStatStyard = new double[2];
        public string MeraZvazKendallaStatStyardName = "";
        public double[] KoefRKorSpir = new double[2];

        public double[] KriterBarkleta = new double[2];
        public double[] ProvRegrs = new double[2];

        public double KoefDeterm = 0;

        public double[] X2f = new double[2];

        public double[,] f;
        public double T = 0;

        public double Szal = 0;
        public double Szal2 = 0;
        private List<double> Z = new List<double>();
        private double[,] Nij = new double[3,3];


        //public double[] AB = new double[4];
        public double[] ABTeil = new double[2];

        public List<Data> Q = new List<Data>();

        public Correlation_RegressionAnalysis(List<InitialStatisticalAnalys> ML, List<int> MSelectGR,string RegresTypeVib)
        {
            this.ML = ML;
            this.MSelectGR = MSelectGR;
            this.RegresTypeVib = RegresTypeVib;
            
            N = NFound(ML, MSelectGR);
            {
                double NezCount = ML[MSelectGR[0]].l.Count;
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (ML[MSelectGR[i]].l.Count != NezCount)
                        Nezal = false;
            }
            T = Distributions.StudentQuantile(1-ML[MSelectGR[0]].alf.Q / 2, ML[MSelectGR[0]].l.Count - 2);
            if (MSelectGR.Count == 2)
            {
                Doubl = true;
                if (Nezal == true)
                {
                    f = fFound(ML[MSelectGR[0]], ML[MSelectGR[1]], 15, 7);
                    Korelation[0] = KorelationFound(ML[MSelectGR[0]],ML[MSelectGR[1]]);
                    KorelationVidnoh[0] = KorelationVidnohFound(ML, MSelectGR);
                    RangKorelation[0] = RangKorelationFound(ML, MSelectGR);
                    RangKorelation[1] = RangKorelation[0] * Math.Sqrt((ML[MSelectGR[0]].l.Count - 2) / (1 - Math.Pow(RangKorelation[0], 2)));


                    X2f[0] = X2fFound(ML[MSelectGR[0]], ML[MSelectGR[1]],f);
                    X2f[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, ML[MSelectGR[0]].l.Count - 2);

                    RangKoefKend[0] = RangKoefKendFound(ML, MSelectGR);//slow
                    RangKoefKend[1] = 3 * RangKoefKend[0] * Math.Sqrt((ML[MSelectGR[0]].l.Count * (ML[MSelectGR[0]].l.Count - 1)) /
                        (2 * (2 * ML[MSelectGR[0]].l.Count + 5)));
                    RangKoefKend[2] = Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2);

                    Nij = NijFound(ML, MSelectGR);
                    int ti = (int)(12*ML[MSelectGR[1]].Dx.Q/ML[MSelectGR[0]].Dx.Q);
                    try
                    {
                        if (Math.Abs(ti) > 20)
                            ti = 20;
                    }
                    catch
                    { ti = 20; }
                    TablPerTab(ML, MSelectGR,fFound(ML[MSelectGR[0]],ML[MSelectGR[1]],12,ti));

                    Q = RegresParamFound(ML[MSelectGR[0]], ML[MSelectGR[1]], Korelation[0], RegresTypeVib, ref Szal2);
                    Szal = SzalF(ML[MSelectGR[0]], ML[MSelectGR[1]], Q, RegresTypeVib);/*
                    AB[1] = Korelation[0] * ML[MSelectGR[1]].Gx.Q / ML[MSelectGR[0]].Gx.Q;
                    AB[0] = ML[MSelectGR[1]].Mx.Q - AB[1] * ML[MSelectGR[0]].Mx.Q;*/
                    NachYslovRegAnal(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                   // Szal = SzalFound(AB, ML[MSelectGR[0]], ML[MSelectGR[1]]);
                    /*AB[2] = Szal * Math.Sqrt(1.0 / ML[MSelectGR[0]].l.Count + Math.Pow(ML[MSelectGR[0]].Mx.Q, 2) /
                        (ML[MSelectGR[0]].Dx.Q * (ML[MSelectGR[0]].l.Count - 1)));
                    AB[3] = Szal / (ML[MSelectGR[0]].Gx.Q * Math.Sqrt(ML[MSelectGR[0]].l.Count - 1));*/


                    ABTailFound(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                    double op = Korelation[0] * (1 - Math.Pow(Korelation[0], 2)) / (2 * ML[MSelectGR[0]].l.Count);
                    double oi = (1 - Math.Pow(Korelation[0], 2)) / Math.Sqrt(ML[MSelectGR[0]].l.Count - 1);
                    double rn = Korelation[0] + op - Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2) * oi;
                    double rv = Korelation[0] + op + Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2) * oi;
                    if (rn < -1)
                        rn = 1;
                    if (rv > 1)
                        rv = 1;
                    Korelation[1] = rn;
                    Korelation[2] = rv;
                    Korelation[3] = op;
                    if (RegresTypeVib == RegresTypeName.ParabRegresion)
                        KoefDeterm = (1 - Math.Pow(Szal2,2)/ML[MSelectGR[1]].Dx.Q)*100;
                    else
                        KoefDeterm = Math.Pow(Korelation[0], 2) * 100;
                }
            }
            else
            {
                Doubl = false;
            }
            Round();
        }
        static public List<Data> RegresParamFound(InitialStatisticalAnalys ISAX, InitialStatisticalAnalys ISAY, double Korelation,string TypeRegVib,ref double Szal2) 
        {
            List<Data> Qm = new List<Data>();
            if (TypeRegVib == RegresTypeName.ParabRegresion)
            {
                Data a = new Data(), b = new Data(), c = new Data();
                Qm.Add(a);
                Qm.Add(b);
                Qm.Add(c);
                a.Name = "a";
                b.Name = "b";
                c.Name = "c";
                double n1 = 0;//107 page andan
                double x2 = InitialStatisticalAnalys.StartMoment(ISAX.l, 2);
                double x3 = InitialStatisticalAnalys.StartMoment(ISAX.l, 3);
                double x4 = InitialStatisticalAnalys.StartMoment(ISAX.l, 4);
                for (int i = 0; i < ISAX.unsortl.Length; i++)
                {
                    n1 += (ISAY.unsortl[i] - ISAY.Mx.Q) * (Math.Pow(ISAX.unsortl[i], 2) - x2);

                }
                n1 /= ISAX.unsortl.Length;
                c.Q = ISAX.Dx.Q * n1 - (x3 - x2 * ISAX.Mx.Q) * Korelation * ISAX.Gx.Q * ISAY.Gx.Q;
                c.Q /= ISAX.Dx.Q * (x4 - Math.Pow(x2, 2)) - Math.Pow(x3 - x2 * ISAX.Mx.Q, 2);
                b.Q = (x4 - Math.Pow(x2, 2)) * Korelation * ISAX.Gx.Q * ISAY.Gx.Q - (x3 - x2 * ISAX.Mx.Q) * n1;
                b.Q /= ISAX.Dx.Q * (x4 - Math.Pow(x2, 2)) - Math.Pow(x3 - x2 * ISAX.Mx.Q, 2);
                a.Q = ISAY.Mx.Q - b.Q * ISAX.Mx.Q - c.Q * ISAX.X_2.Q;
                Data a2 = new Data(), b2 = new Data(), c2 = new Data();
                Qm.Add(a2);
                Qm.Add(b2);
                Qm.Add(c2);
                a2.Name = "a2";
                b2.Name = "b2";
                c2.Name = "c2";
                a2.Q = ISAY.Mx.Q;
                double Mfi2scv = 0;
                {
                    double TD = 0;
                    for (int i = 0; i < ISAX.l.Count; i++)
                    {
                        Mfi2scv += Math.Pow(fi2F(ISAX.unsortl[i], ISAX.Dx.Q, ISAX.Mx.Q, x2, x3),2);
                        b2.Q += (ISAX.unsortl[i] - ISAX.Mx.Q) * ISAY.unsortl[i];
                        c2.Q += fi2F(ISAX.unsortl[i], ISAX.Dx.Q, ISAX.Mx.Q, x2, x3) * ISAY.unsortl[i];
                        TD+= Math.Pow(fi2F(ISAX.unsortl[i], ISAX.Dx.Q, ISAX.Mx.Q, x2, x3), 2);
                    }
                    c2.Q /= TD;
                }
                Mfi2scv /= ISAX.l.Count;
                b2.Q /= ISAX.l.Count;
                b2.Q /= ISAX.Dx.Q;
                Szal2 = 0;
                for (int i = 0; i < ISAX.l.Count; i++)
                {
                    Szal2 += Math.Pow(ISAY.unsortl[i] - a2.Q - b2.Q * fi1F(ISAX.unsortl[i], ISAX.Mx.Q) - c2.Q * fi2F(ISAX.unsortl[i], ISAX.Dx.Q, ISAX.Mx.Q, x2, x3), 2);
                }
                Szal2 /= ISAX.l.Count - 3;
                Szal2 = Math.Sqrt(Szal2);
                a2.QSigma = Szal2 / Math.Sqrt(ISAX.unsortl.Length);
                b2.QSigma = Szal2 / (ISAX.Dx.Q * Math.Sqrt(ISAX.unsortl.Length));
                c2.QSigma = Szal2 / Math.Sqrt(ISAX.unsortl.Length * Mfi2scv);
                double Tt = Distributions.StudentQuantile(1 - ISAX.alf.Q / 2, ISAX.unsortl.Length - 3);
                a2.QButton = a2.Q - Tt * a2.QSigma;
                a2.QUpper = a2.Q + Tt * a2.QSigma;
                b2.QButton = b2.Q - Tt * b2.QSigma;
                b2.QUpper = b2.Q + Tt * b2.QSigma;
                c2.QButton = c2.Q - Tt * c2.QSigma;
                c2.QUpper = c2.Q + Tt * c2.QSigma;
                double at = ISAY.Mx.Q - b.Q * ISAX.Mx.Q - c.Q * Math.Pow(ISAX.Mx.Q, 2) - a.Q;
                Data ta = new Data(), tb = new Data(), tc = new Data();

            }
            else
            {
                Data a = new Data(), b = new Data();
                Qm.Add(a);
                Qm.Add(b);
                a.Name = "a";
                b.Name = "b";
                List<double> t = new List<double>();
                List<double> z = new List<double>();
                for(int i =0;i<ISAX.unsortl.Length;i++)
                {
                    t.Add(RegresType.FiX(ISAX.unsortl[i], TypeRegVib));
                    z.Add(RegresType.FiY(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib));
                }
                /*double Mfx = 0,Mfxfx=0, Mfy = 0, Mfxfy = 0,W = 0;
                for (int i = 0; i < ISAX.unsortl.Length; i++)
                {
                    Mfx += RegresType.FiX(ISAX.unsortl[i], TypeRegVib) *
                        RegresType.W(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib);
                    Mfxfx += Math.Pow(RegresType.FiX(ISAX.unsortl[i], TypeRegVib),2) *
                        RegresType.W(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib);
                    Mfy += RegresType.FiY(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib) *
                        RegresType.W(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib);
                    Mfxfy += RegresType.FiX(ISAX.unsortl[i], TypeRegVib) * RegresType.FiY(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib) *
                        RegresType.W(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib);
                    W += RegresType.W(ISAY.unsortl[i], ISAX.unsortl[i], TypeRegVib);

                }
                Mfx /= W;
                Mfxfx /= W;
                Mfxfy /= W;
                Mfy /= W;*/
               // Qm[1].Q = (Mfxfy - Mfx * Mfy) / (Mfxfx - Math.Pow(Mfx, 2));
                //Qm[0].Q = Mfy - Qm[1].Q * Mfx;
                InitialStatisticalAnalys ISAt = new InitialStatisticalAnalys(t);
                InitialStatisticalAnalys ISAz = new InitialStatisticalAnalys(z);
                double Kor_tz = Correlation_RegressionAnalysis.KorelationFound(ISAt, ISAz);
                Qm[1].Q = Kor_tz * ISAz.Gx.Q / ISAt.Gx.Q;
                Qm[0].Q = RegresType.A(ISAz.Mx.Q - Qm[1].Q * ISAt.Mx.Q, TypeRegVib);
                double Szal = SzalF(ISAX, ISAY, Qm, TypeRegVib);
                Qm[0].QSigma = RegresType.A(Szal * Math.Sqrt(1.0 / ISAX.unsortl.Length + Math.Pow(ISAt.Mx.Q, 2) / (ISAt.Dx.Q * (ISAX.unsortl.Length - 1))), TypeRegVib);
                Qm[1].QSigma = Szal / (ISAX.unsortl.Length - 1);
                Qm[0].QButton = RegresType.A(Qm[0].Q - Distributions.StudentQuantile(1 - ISAX.alf.Q / 2, ISAX.unsortl.Length - 2) * Qm[0].QSigma, TypeRegVib);
                Qm[0].QUpper = RegresType.A(Qm[0].Q + Distributions.StudentQuantile(1 - ISAX.alf.Q / 2, ISAX.unsortl.Length - 2) * Qm[0].QSigma, TypeRegVib);
                Qm[1].QButton = Qm[1].Q - Distributions.StudentQuantile(1 - ISAX.alf.Q / 2, ISAX.unsortl.Length - 2) * Qm[1].QSigma;
                Qm[1].QUpper  = Qm[1].Q + Distributions.StudentQuantile(1 - ISAX.alf.Q / 2, ISAX.unsortl.Length - 2) * Qm[1].QSigma;


            }
            return Qm;

        }
        static public string TypeRegresFound(InitialStatisticalAnalys ISAX, InitialStatisticalAnalys ISAY) 
        {
            string[] Type = Regex.Split(RegresTypeName.TypeRegresion, "\n");
            double[] Szalmas = new double[Type.Length];
            double Korel = KorelationFound(ISAX,ISAY);
            for (int i = 0; i < Type.Length;i++ ) 
            {
                double td = 0;
                List<Data> q = RegresParamFound(ISAX, ISAY,Korel, Type[i],ref td);
                Szalmas[i] =  SzalF(ISAX, ISAY, q, Type[i]);
            }
            double min = Szalmas[0];
            int index = 0;
            for (int i = 0; i < Type.Length; i++)
            {
                if (min!=double.NaN&&Szalmas[i] <= min) 
                {
                    min = Szalmas[i];
                    index = i;
                }
            }
            return Type[index];
        }
        static public double SzalF(InitialStatisticalAnalys ISAX,InitialStatisticalAnalys ISAY,List<Data> Q,string RegresTypeN)
        {
            double Szal = 0;
            for (int i =0;i<ISAX.unsortl.Length;i++)
            {
                Szal += Math.Pow(ISAY.unsortl[i] - RegresType.Model(ISAX.unsortl[i], Q, RegresTypeN), 2);
            }
                Szal/=(ISAX.unsortl.Length - Q.Count);
            return Math.Sqrt(Szal);
        }
        private void ABTailFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2) 
        {
            double B = 0;
            double A = 0;
            List<double> b = new List<double>();
            for (int i = 1; i < gr1.l.Count; i++) 
                for (int j = 0; j < i;j++ )
                    b.Add((gr2.unsortl[j] - gr2.unsortl[i]) / (gr1.unsortl[j] - gr1.unsortl[i]));
            b.Sort();
            if (b.Count % 2 == 0)
                B = (b[b.Count / 2 - 1] + b[b.Count / 2]) / 2;
            else
                B = b[(b.Count - 1) / 2];
            List<double> a = new List<double>();
            for (int i = 0; i < gr1.l.Count; i++)
            {
                a.Add(gr2.unsortl[i] - B * gr1.unsortl[i]);
            }
            a.Sort();
            if (a.Count % 2 == 0)
                A = (a[a.Count / 2 - 1] + a[a.Count / 2]) / 2;
            else
                A = a[(a.Count - 1) / 2];
            ABTeil[0] = A;
            ABTeil[1] = B;
        }
        private double X2fFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2,double[,] f) 
        {
            double rez = 0;
            for (int i =0;i<f.GetLength(0) ;i++ ) 
                for (int j = 0; j < f.GetLength(1); j++)
                {
                    double f1 = Distributions.NormalDoublefFound(
                        ((i + 0.5) * gr1.Len.Q / f.GetLength(0) + gr1.Min.Q), gr1.Mx.Q, gr1.Gx.Q, 
                        ((j + 0.5) * gr2.Len.Q / f.GetLength(1) + gr2.Min.Q), gr2.Mx.Q, gr2.Gx.Q, Korelation[0]) *
                        (gr2.Len.Q / f.GetLength(1))*(gr1.Len.Q / f.GetLength(0));
                    if(Math.Round(f1,4)!=0)
                        rez += Math.Pow(f[i, j] - f1, 2) / f1;
                }
            return rez;
        }
        private bool NachYslovRegAnal(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            bool rez = true;
            if (!(gr1.AvtoType == gr2.AvtoType && gr2.AvtoType == "Нормальний"))
                rez = false;
            int M = 18;
            List<double>[] yi = new List<double>[M];
            for (int i = 0; i < M; i++)
                yi[i] = new List<double>();
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
            {
                for (int j = 0; j < M; j++)
                    if (gr1.unsortl[i] < (j + 1.0001) * gr1.Len.Q / M + gr1.Min.Q && gr1.unsortl[i] >= (j * gr1.Len.Q / M + gr1.Min.Q))
                    {
                        yi[j].Add(gr2.unsortl[i]);
                        break;
                    }
            }
            ///ProvRegr
            double t1 = 0, t2 = 0;

            ///


            double[] Sxj = new double[M];
            double S = 0;
            for (int i = 0; i < M; i++) 
            {
                double ysr = yi[i].Sum() / yi[i].Count;
                if(ysr>0)
                {
                     for (int j = 0; j < yi[i].Count;j++ ) 
                     {
                         t2 += Math.Pow(yi[i][j] - ysr, 2);
                         Sxj[i] += Math.Pow(yi[i][j] - ysr, 2);
                     }
                     if (yi[i].Count  >1)
                         Sxj[i] /= yi[i].Count - 1;
                     t1 += yi[i].Count * Math.Pow(ysr - Q[0].Q - Q[1].Q * ((i + 0.5) * gr1.Len.Q / M + gr1.Min.Q),2);
                }
            }
            ProvRegrs[0] = (gr1.l.Count - M)*t1/((M - 1)*t2);
            ProvRegrs[1] = Distributions.FisherQuantile(gr1.alf.Q, M - 1, gr1.l.Count - M);

            double C = 0;
            for (int i = 0; i < M;i++ ) 
            {
                S += (yi[i].Count - 1) * Sxj[i];
                if (yi[i].Count != 0) 
                    C += 1.0 / yi[i].Count-1.0 / gr1.l.Count;
            }
            S /= gr1.l.Count - M;
            //C -= 1 / gr1.l.Count;
            C /= 3.0 * (M - 1);
            C++;
            for (int i = 0; i < M; i++)
            {
                if (Sxj[i] != 0 )
                    KriterBarkleta[0] += yi[i].Count * Math.Log(Sxj[i] / S);
            }
            KriterBarkleta[0] /= -C;
            KriterBarkleta[1] = Hi.HIF(gr1.alf.Q, M - 1);



            return rez;
        }
        void TablPerTab(List<InitialStatisticalAnalys> ML, List<int> MSelectGR,double[,] f)
        { 
            double[] n = new double[f.GetLength(1)];
            double[] m = new double[f.GetLength(0)];
            for (int i = 0; i < f.GetLength(0);i++ )
                for (int j = 0; j < f.GetLength(1); j++)
                    m[i] += f[i, j] * ML[MSelectGR[0]].l.Count;
            for (int i = 0; i < f.GetLength(1); i++)
                for (int j = 0; j < f.GetLength(0); j++)
                    n[i] += f[j, i] * ML[MSelectGR[1]].l.Count;
            for (int i = 0; i < f.GetLength(1); i++)
                for (int j = 0; j < f.GetLength(0); j++)
                    NezXY[0] += Math.Pow(f[j, i] * ML[MSelectGR[0]].l.Count - n[i] * m[j] / ML[MSelectGR[0]].l.Count, 2) / (n[i] * m[j] / ML[MSelectGR[0]].l.Count);
            NezXY[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, (f.GetLength(1) - 1) * (f.GetLength(0) - 1));
            KoefSpolPirs[0] = Math.Sqrt(NezXY[0] / (NezXY[0] + ML[MSelectGR[0]].l.Count));

            double P = 0, Q = 0;
            double Stsig = 0;
            for (int i = 0; i < f.GetLength(1); i++)
                for (int j = 0; j < f.GetLength(0); j++)
                {
                    double t1 = 0, t2 = 0,A1=0,B1=0,A2=0,B2=0;
                    for (int k = i; k < f.GetLength(1); k++)
                    {
                        for (int l = j ; l < f.GetLength(0); l++)
                        {
                            t1 += f[l, k] * ML[MSelectGR[1]].l.Count;
                            A1 += f[l, k] * ML[MSelectGR[1]].l.Count;
                            A2 += f[l - j, k - i] * ML[MSelectGR[1]].l.Count;
                            B1 += f[l, k - i] * ML[MSelectGR[1]].l.Count;
                            B2 += f[l - j , k] * ML[MSelectGR[1]].l.Count;

                        }
                        for (int l = 0; l < j; l++)
                        {
                            t2 += f[l, k] * ML[MSelectGR[1]].l.Count;
                        }
                    }
                    P += f[j, i] * t1 * ML[MSelectGR[1]].l.Count;
                    Q += f[j, i] * t2 * ML[MSelectGR[1]].l.Count;
                    Stsig += f[j, i] * ML[MSelectGR[1]].l.Count * Math.Pow( A1 + A2 - B1 - B2, 2);
                }
            if (f.GetLength(0) == f.GetLength(1))
            {
                MeraZvazKendallaStatStyardName = "Міра зв'язку Кендалла";
                double T1 = 0, T2 = 0;
                for (int j = 0; j < f.GetLength(1); j++)
                {
                    T1 += n[j] * (n[j] - 1);
                }
                for (int j = 0; j < f.GetLength(0); j++)
                {
                    T2 += m[j] * (m[j] - 1);
                }
                T1 /= 2;
                T2 /= 2;
                MeraZvazKendallaStatStyard[0] = (P - Q) / Math.Sqrt((ML[MSelectGR[1]].l.Count * (ML[MSelectGR[1]].l.Count - 1) / 2.0 - T1) *
                    (ML[MSelectGR[1]].l.Count * (ML[MSelectGR[1]].l.Count - 1) / 2.0 - T2));
                MeraZvazKendallaStatStyard[1] = Math.Sqrt((f.GetLength(1) * 4.0 + 10) / (9.0 * (Math.Pow(f.GetLength(1), 2) - f.GetLength(1))));
            }
            else 
            {
                MeraZvazKendallaStatStyardName = "Статистика Стюарда";
                MeraZvazKendallaStatStyard[0] = 2 * (P - Q) * Math.Min(f.GetLength(0), f.GetLength(1)) /
                    (Math.Pow(ML[MSelectGR[1]].l.Count, 2) * (Math.Min(f.GetLength(0), f.GetLength(1)) - 1));
                MeraZvazKendallaStatStyard[1] = 2 * Math.Min(f.GetLength(0), f.GetLength(1)) *
                    Math.Sqrt(Math.Pow(ML[MSelectGR[1]].l.Count, 2) * Stsig - 4.0 * ML[MSelectGR[1]].l.Count * (P - Q)) /
                    (Math.Pow(ML[MSelectGR[1]].l.Count, 3) * (Math.Min(f.GetLength(0), f.GetLength(1)) - 1));

            }
            double RKS =0;
            for (int i = 0; i < f.GetLength(1); i++)
                for (int j = 0; j < f.GetLength(0); j++)
                {
                    double S1 = 0, S2 = 0;
                    for (int k = 0; k < i; k++)
                        S1 += 1.5 * n[k] ;
                    for (int l = 0; l < j; l++)
                        S2 += 1.5 * m[l] ;
                    RKS += f[j, i] * ML[MSelectGR[1]].l.Count * (S1 - ML[MSelectGR[1]].l.Count / 2.0) * (S2 - ML[MSelectGR[1]].l.Count / 2.0);
                }

            double Mi = 0, Ni = 0;
            for (int i = 0; i < f.GetLength(1); i++) 
            {
                Ni += Math.Pow(n[i], 3) - n[i];
            }

            for (int i = 0; i < f.GetLength(0); i++) 
            {
                Mi += Math.Pow(m[i], 3) - m[i];
            }
            KoefRKorSpir[0] = 12.0 * RKS / Math.Sqrt((Math.Pow(ML[MSelectGR[1]].l.Count, 3) - f.GetLength(1) - Ni) *
                (Math.Pow(ML[MSelectGR[1]].l.Count, 3) - f.GetLength(1) - Mi));
            KoefRKorSpir[1] = (1 - Math.Pow(KoefRKorSpir[0], 2)) / (Math.Sqrt(ML[MSelectGR[1]].l.Count - 2));

        }
        double[,] NijFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR) 
        {
            double[,] rez = new double[3,3];
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++) 
            {
                if (ML[MSelectGR[0]].unsortl[i] > ML[MSelectGR[0]].Mx.Q && ML[MSelectGR[1]].unsortl[i] > ML[MSelectGR[1]].Mx.Q)
                    rez[1, 1]++;
                else if (ML[MSelectGR[0]].unsortl[i] > ML[MSelectGR[0]].Mx.Q && ML[MSelectGR[1]].unsortl[i] <= ML[MSelectGR[1]].Mx.Q)
                    rez[0, 1]++;
                else if (ML[MSelectGR[0]].unsortl[i] <= ML[MSelectGR[0]].Mx.Q && ML[MSelectGR[1]].unsortl[i] > ML[MSelectGR[1]].Mx.Q)
                    rez[1, 0]++;
                else if (ML[MSelectGR[0]].unsortl[i] <= ML[MSelectGR[0]].Mx.Q && ML[MSelectGR[1]].unsortl[i] <= ML[MSelectGR[1]].Mx.Q)
                    rez[0, 0]++;
            }
            rez[2, 0] = rez[0, 0] + rez[0, 1];///N0
            rez[2, 1] = rez[1, 0] + rez[1, 1];///N1
            rez[0, 2] = rez[0, 0] + rez[1, 0];///M0
            rez[1, 2] = rez[0, 1] + rez[1, 1];///M1
            rez[2, 2] = rez[0, 0] + rez[0, 1] + rez[1, 0] + rez[1, 1];///N
            double V = rez[0,0]+rez[1,1];
            double W = rez[0,1]+rez[1,0];
            IndexFehnera = (V - W) / (V + W);
            KoefSpolF[0] = (rez[0, 0] * rez[1, 1] - rez[0, 1] * rez[1, 0])/Math.Sqrt(rez[0,2]*rez[1,2] * rez[2,0] * rez[2,1]);
            KoefSpolF[1] = Math.Pow(KoefSpolF[0], 2) * ML[MSelectGR[0]].l.Count;
            KoefSpolF[2] = Hi.HIF(ML[MSelectGR[1]].alf.Q, 1);
            KoefSvazYoola[0] = (rez[0, 0] * rez[1, 1] - rez[0, 1] * rez[1, 0]) / (rez[0, 0] * rez[1, 1] + rez[0, 1] * rez[1, 0]);
            KoefSvazYoola[1] = (Math.Sqrt(rez[0, 0] * rez[1, 1]) - Math.Sqrt(rez[0, 1] * rez[1, 0])) / 
                (Math.Sqrt(rez[0, 0] * rez[1, 1]) + Math.Sqrt(rez[0, 1] * rez[1, 0]));
            KoefSvazYoola[2] = KoefSvazYoola[0]/(Math.Sqrt(1.0 / rez[0, 0] + 1.0 / rez[1, 0] + 1.0 / rez[0, 1] + 1.0 / rez[1, 1]) *
                (1 - Math.Pow(KoefSvazYoola[0],2))/2.0);
            KoefSvazYoola[3] = KoefSvazYoola[1] / (Math.Sqrt(1.0 / rez[0, 0] + 1.0 / rez[1, 0] + 1.0 / rez[0, 1] + 1.0 / rez[1, 1]) *
                (1 - Math.Pow(KoefSvazYoola[1], 2)) / 4.0);
            return rez;
        }
        double KrKohrenaFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR, double BinSim)
        {
            double Q = 0;
            double[] u = new double[ML[MSelectGR[0]].l.Count];
            double[] t = new double[MSelectGR.Count];
            for (int i=0;i<ML[MSelectGR[0]].l.Count;i++)
                    u[i] = MSelectGR.Count;
            for (int i=0;i<MSelectGR.Count;i++)
                    t[i] = ML[MSelectGR[0]].l.Count;
            for (int i=0;i<ML[MSelectGR[0]].l.Count;i++)
                for (int j=0;j<MSelectGR.Count;j++)
                        if (ML[MSelectGR[j]].unsortl[i]<=BinSim)
                        {
                            t[j]--;
                            u[i]--;
                        }
            double mt =0,mu=0;
            for (int i=0;i<ML[MSelectGR[0]].l.Count;i++)
                    mu+=u[i]/ML[MSelectGR[0]].l.Count;
            for (int i=0;i<MSelectGR.Count;i++)
                    mt+=t[i]/MSelectGR.Count;
            double eu=u.Sum(),et=0,euu=0;
            for (int i=0;i<MSelectGR.Count;i++)
                    et+=Math.Pow(t[i] -mt,2);
            for (int i=0;i<ML[MSelectGR[0]].l.Count;i++)
                    euu+=Math.Pow(u[i],2);
            if (MSelectGR.Count * eu - euu==0)
            {
                return -1;
            }
            Q=MSelectGR.Count*(MSelectGR.Count-1)*et/(MSelectGR.Count*eu - euu); 

            return Q;
        }

        public static double KorelationFound(InitialStatisticalAnalys ISAX,InitialStatisticalAnalys ISAY)
        {
            double Mxy = 0;
            double sx = 0, sy = 0;
            for (int i = 0; i < ISAX.l.Count; i++)
            {
                sx += ISAX.unsortl[i];
                sy += ISAY.unsortl[i];
                Mxy += ISAX.unsortl[i] * ISAY.unsortl[i];
            }
            var ty = Mxy / Math.Sqrt(sx * sy);
            Mxy /= ISAX.l.Count;

            return (ISAX.l.Count / (ISAX.l.Count - 1)) *
                ((Mxy - ISAX.Mx.Q * ISAY.Mx.Q) / (ISAX.Gx.Q * ISAY.Gx.Q));
        }
        static public bool KorelationZnach(double Korelation,int N,double alf)
        {
            bool Znach = false;
            double T = Korelation*Math.Sqrt((N-2)/(1 - Math.Pow(Korelation,2)));
            double TQv = Distributions.StudentQuantile(1 - alf / 2, N - 2);
            Znach = Math.Abs(T) <= TQv;
            return !Znach;
        }
        double KorelationVidnohFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double rez = 0;
            int M=10;
            List<double>[] yi= new List<double>[M];
            for (int i = 0; i < M; i++)
                yi[i] = new List<double>();
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
            {
                for (int j = 0; j < M; j++)
                        if (ML[MSelectGR[0]].unsortl[i] <= (j + 1) * ML[MSelectGR[0]].Len.Q / M + ML[MSelectGR[0]].Min.Q &&
                            ML[MSelectGR[0]].unsortl[i] >= (j * ML[MSelectGR[0]].Len.Q / M + ML[MSelectGR[0]].Min.Q))
                        {
                            yi[j].Add(ML[MSelectGR[1]].unsortl[i]);
                            break;
                        }
            }
            double Emyjy = 0;
            for (int i = 0; i < M; i++)
            {
                var t1 = yi[i].Sum();
                t1/= yi[i].Count;
                Emyjy += yi[i].Count * Math.Pow(yi[i].Sum() / yi[i].Count - ML[MSelectGR[1]].Mx.Q, 2);
            }
            rez = Emyjy / (ML[MSelectGR[1]].Dx.Q * (ML[MSelectGR[1]].l.Count - 1));
            rez = Math.Sqrt(rez);
            return rez;
        }
        double RangKorelationFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR) 
        {
            double d_2 = 0;
            double[] rx = RDvomFound(ML[MSelectGR[0]].unsortl);
            double[] ry = RDvomFound(ML[MSelectGR[1]].unsortl);
            for (int i=0;i<ML[MSelectGR[0]].l.Count ;i++ ) 
            {
                d_2 += Math.Pow(rx[i] - ry[i], 2);
            }
            return 1 -  d_2 / (ML[MSelectGR[0]].l.Count * (Math.Pow(ML[MSelectGR[0]].l.Count,2)-1));
        }
        double RangKoefKendFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double s = 0;
            double[] rangY = RDvomFound(ML[MSelectGR[1]].unsortl);
            for (int i = 0; i < ML[MSelectGR[0]].l.Count-1; i++)
            {
                for (int j = i+1; j < ML[MSelectGR[0]].l.Count; j++)
                {
                    if (rangY[i] < rangY[j])
                        s++;
                    else if (rangY[i] < rangY[j])
                        s--;
                }
            }
            return 2*s/ (ML[MSelectGR[1]].l.Count * (ML[MSelectGR[1]].l.Count - 1));
        }
        double NFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR) 
        {
            double hn=0;
            for (int i = 0; i < MSelectGR.Count; i++)
            {
                if (Nmax < ML[MSelectGR[i]].l.Count)
                    Nmax = ML[MSelectGR[i]].l.Count;
                hn += ML[MSelectGR[i]].l.Count;
            }
            return hn;
        }
        double[] RDvomFound(double[] mas)
        {
            double[] rez = new double[mas.Length];
            List<double[]> dat = new List<double[]>();
            for (int i = 0;i<mas.Length;i++ )
                dat.Add(new double[]{mas[i],i,0});
            dat.Sort(delegate(double[] t1,double[] t2)
            {return(t1[0].CompareTo(t2[0]));});
            for (int i = 0;i<dat.Count ;i++ )
                dat[i][2] = i + 1;
            {
                List<double[]> simpldata = new List<double[]>();
                double prev = dat[0][0];
                simpldata.Add(new double[] { dat[0][2] - 1, prev });
                for (int i = 1; i < dat.Count; prev = dat[i][2], i++)
                {
                    if (dat[i][0] != prev)
                    {
                        if (simpldata.Count > 1)
                        {
                            double sum = 0;
                            for (int j = 0; j < simpldata.Count; sum += simpldata[j++][0])
                                ;
                            for (int j = 0; j < simpldata.Count; j++)
                                dat[i - j - 1][2] = sum / simpldata.Count;
                        }
                        simpldata.Clear();
                    }
                    simpldata.Add(new double[] { dat[i][2], dat[i][0] });
                }
            }
            dat.Sort(delegate(double[] t1, double[] t2)
            { return (t1[1].CompareTo(t2[1])); });
            for (int i = 0; i < mas.Length; i++)
                rez[i]=dat[i][2];
            return rez;
        }
#pragma warning disable IDE1006
        double[,] fFound(InitialStatisticalAnalys gr1,InitialStatisticalAnalys gr2,int Ncol,int Nrow)
        {
            double[,] rez= new double[Ncol,Nrow];
            for (int i=0;i<gr1.l.Count ;i++ ) 
            {
                for (int j = 0; j < Ncol; j++)
                    for (int k = 0; k < Nrow; k++)
                        if (gr1.unsortl[i] <= (j + 1) * gr1.Len.Q / Ncol + gr1.Min.Q && gr1.unsortl[i] >= (j * gr1.Len.Q / Ncol + gr1.Min.Q) &&
                            gr2.unsortl[i] <= ((k + 1) * gr2.Len.Q / Nrow + gr2.Min.Q) && gr2.unsortl[i] >= (k * gr2.Len.Q / Nrow + gr2.Min.Q))
                            {
                                rez[j, k] += 1.0 / gr1.l.Count;
                                goto t1;
                            }
            t1:;
            }
            /*for (int j = 0; j < Ncol; j++)
                for (int k = 0; k < Nrow; k++)
                    if (rez[j, k] <= gr1.alf)
                        rez[j, k] = 0;*/
            return rez;
        }
#pragma warning restore IDE1006
        double CFound(double i, double Nd)
        {
            double iF = 1, NF = 1,NiF=1;
            for (int j = 1; j <= i;j++ )
                iF *= j;
            for (int j = 1; j <= Nd; j++)
                NF *= j;
            for (int j = 1; j <= Nd-i; j++)
                NiF *= j;
            return iF / (NiF * NF);
        }
        /*private double SzalFound(double[] AB,InitialStatisticalAnalys X,InitialStatisticalAnalys Y)
        {
            double rez = 0;
            for (int i =0;i<Y.l.Count ;i++ ) 
            {
                rez += Math.Abs(Y.unsortl[i] - AB[0] - AB[1] * X.unsortl[i]);
            }
            rez /= X.l.Count;
            return rez;

        }*/
        private double ZFound( InitialStatisticalAnalys gr1,InitialStatisticalAnalys gr2)
        {
            double D = 0;
            for (int i = 0; i < gr1.l.Count; i++)
            {
                double l2 = 0;
                if (gr2.AvtoType == "Нормальний")
                    l2 = Distributions.NormalfFound(gr1.l[i], gr2.Mx.Q, gr2.Gx.Q);
                else if (gr2.AvtoType == "Еспоненціальний")
                    l2 = Distributions.ExpFound(gr1.l[i], gr2.Mx.Q);
                else if (gr2.AvtoType == "Рівномірний")
                    l2 = (gr1.l[i] - gr2.Min.Q) * gr2.Len.Q;

                D = Math.Max(Math.Abs(gr1.F[i]- l2), D);
            }
            return D;
        }


#pragma warning disable IDE1006
        static public double[] fi1F(InitialStatisticalAnalys ISA)
        {
            double[] rez = new double[ISA.unsortl.Length]; ///page 108
            for (int i = 0; i < rez.Length; i++)
            {
                rez[i] = ISA.unsortl[i] - ISA.Mx.Q;
            }
            return rez;
            //return element - mx;
        }
        static public double[] fi2F(InitialStatisticalAnalys ISA)
        {
            double[] rez = new double[ISA.unsortl.Length]; ///page 108
            double x3 = InitialStatisticalAnalys.StartMoment(ISA.l, 3);
            double x2 = InitialStatisticalAnalys.StartMoment(ISA.l, 2);
            for (int i = 0; i < rez.Length;i++ )
            {
                rez[i] = Math.Pow(ISA.unsortl[i], 2) - (x3 - x2 * ISA.Mx.Q) * (ISA.unsortl[i] - ISA.Mx.Q) / ISA.Dx.Q - x2;
            }
            return rez;
        }
        static public double fi1F(double element,double mx)
        {///page 108
            return element - mx;
        }
        static public double fi2F(double element, double dx, double mx, double mx2, double mx3)
#pragma warning restore IDE1006
        {
            return Math.Pow(element, 2) - ((mx3 - mx2 * mx) * (element - mx) / dx)- mx2;
        }
        private void Round()
        {
            int o = 4;

            Korelation[0] = Math.Round(Korelation[0], o);
            Korelation[1] = Math.Round(Korelation[1], o);
            Korelation[2] = Math.Round(Korelation[2], o);
            Korelation[3] = Math.Round(Korelation[3], o);


            RangKorelation[0] = Math.Round(RangKorelation[0], o);
            RangKorelation[1] = Math.Round(RangKorelation[1], o);

            KorelationVidnoh[0] = Math.Round(KorelationVidnoh[0], o);



            RangKoefKend[0] = Math.Round(RangKoefKend[0], o);
            RangKoefKend[1] = Math.Round(RangKoefKend[1], o);
            RangKoefKend[2] = Math.Round(RangKoefKend[2], o);




            IndexFehnera = Math.Round(IndexFehnera, o);

            KoefSpolF[0] = Math.Round(KoefSpolF[0], o);
            KoefSpolF[1] = Math.Round(KoefSpolF[1], o);
            KoefSpolF[2] = Math.Round(KoefSpolF[2], o);
            KoefSvazYoola[0] = Math.Round(KoefSvazYoola[0], o);
            KoefSvazYoola[1] = Math.Round(KoefSvazYoola[1], o);
            KoefSvazYoola[2] = Math.Round(KoefSvazYoola[2], o);
            KoefSvazYoola[3] = Math.Round(KoefSvazYoola[3], o);

            MeraZvazKendallaStatStyard[0] = Math.Round(MeraZvazKendallaStatStyard[0], o);
            MeraZvazKendallaStatStyard[1] = Math.Round(MeraZvazKendallaStatStyard[1], o);

            KoefRKorSpir[0] = Math.Round(KoefRKorSpir[0], o);
            KoefRKorSpir[1] = Math.Round(KoefRKorSpir[1], o);

            KriterBarkleta[0] = Math.Round(KriterBarkleta[0], o);
            KriterBarkleta[1] = Math.Round(KriterBarkleta[1], o);

            KoefDeterm = Math.Round(KoefDeterm, 0);

            /*AB[0] = Math.Round(AB[0], o);
            AB[1] = Math.Round(AB[1], o);
            AB[2] = Math.Round(AB[2], o);
            AB[3] = Math.Round(AB[3], o);*/

            ABTeil[0] = Math.Round(ABTeil[0], o);
            ABTeil[1] = Math.Round(ABTeil[1], o);

            ProvRegrs[0] = Math.Round(ProvRegrs[0], o);
            ProvRegrs[1] = Math.Round(ProvRegrs[1], o);

            X2f[0] = Math.Round(X2f[0], o);
            X2f[1] = Math.Round(X2f[1], o);
        }
    }
}
