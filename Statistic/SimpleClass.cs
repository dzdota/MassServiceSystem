using System;
using System.Collections.Generic;
using System.Linq;

namespace testgistogr
{
    class SimpleClass
    {
        List<InitialStatisticalAnalys> ML;
        List<int> MSelectGR;
        public bool Doubl = false;
        public bool Nezal = true;
        public double N = 0;
        public double Nmax = 0;
        public double[] SravnSred = new double[2];
        public double[] SravnDisper = new double[2];
        public double[] VarSv2Sm2 = new double[2];
        public double Sv2= 0;
        public double Sm2= 0;
        public double[] KrSmirnKolmag = new double[2];///krt>alf
        public double[] KrsumRangVils = new double[2];///krt<=u alf/2
        public double[] KrUMannaUit = new double[2];///krt<=u alf
        public double[] RizSerRangVib = new double[2];///krt<=u alf
        public double[] KrZnakiv = new double[2];
        public double[] KrKruskalaUolisa = new double[2];
        public double[] KrKohrena = new double[2];
                                               
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
        public double[] RegresParam;
        public double[] AB = new double[4];
        public double[] ABTeil = new double[2];
        public double[] ProvRegrs = new double[2];

        public double KoefDeterm = 0;

        public double[] X2f = new double[2];

        public double[,] f;
        public double T = 0;

        public double Szal = 0;
        private List<double> Z = new List<double>();
        private double[,] r;
        private double[,] Nij = new double[3,3];

        public SimpleClass(List<InitialStatisticalAnalys> ML, List<int> MSelectGR,double BinSim)
        {
            this.ML = ML;
            this.MSelectGR = MSelectGR;
            
            N = NFound(ML, MSelectGR);
            Nezal = Nezalegni();
             
            r = rFound(ML, MSelectGR);
            Sm2 = Sm2Found(ML, MSelectGR);
            Sv2 = Sv2Found(ML, MSelectGR);
            VarSv2Sm2[0] = Sm2 / Sv2 ;
            VarSv2Sm2[1] = Distributions.FisherQuantile(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1, (int)(N - MSelectGR.Count));
            var tyui = Distributions.FisherQuantile(ML[MSelectGR[0]].alf.Q, 10,10);
            if (Nezal == true)
            {
                KrKohrena[0] = KrKohrenaFound(ML,MSelectGR, BinSim);
                KrKohrena[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1);
                T = Distributions.StudentQuantile(1-ML[MSelectGR[0]].alf.Q / 2, ML[MSelectGR[0]].l.Count - 2);
            }
            if (MSelectGR.Count == 2)
            {
                Doubl = true;
                if (Nezal == true)
                {
                    f = fFound(ML[MSelectGR[0]], ML[MSelectGR[1]], 15, 7);
                    Korelation[0] = KorelationFound(ML, MSelectGR);
                    KorelationVidnoh[0] = KorelationVidnohFound(ML, MSelectGR);
                    RangKorelation[0] = RangKorelationFound(ML, MSelectGR);
                    RangKorelation[1] = RangKorelation[0] * Math.Sqrt((ML[MSelectGR[0]].l.Count - 2) / (1 - Math.Pow(RangKorelation[0], 2)));


                    X2f[0] = X2fFound(ML[MSelectGR[0]], ML[MSelectGR[1]],f);
                    X2f[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, ML[MSelectGR[0]].l.Count - 2);

                    RangKoefKend[0] = RangKoefKendFound(ML, MSelectGR);//slow
                    RangKoefKend[1] = 3 * RangKoefKend[0] * Math.Sqrt((ML[MSelectGR[0]].l.Count * (ML[MSelectGR[0]].l.Count - 1)) /
                        (2 * (2 * ML[MSelectGR[0]].l.Count + 5)));
                    RangKoefKend[2] = Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2);

                    Nij = NijFound(ML, MSelectGR, BinSim);
                    int ti = (int)(12*ML[MSelectGR[1]].Dx.Q/ML[MSelectGR[0]].Dx.Q);
                    try
                    {
                        if (Math.Abs(ti) > 20)
                            ti = 20;
                    }
                    catch
                    { ti = 20; }
                    TablPerTab(ML, MSelectGR,fFound(ML[MSelectGR[0]],ML[MSelectGR[1]],12,ti));

                    RegresParam = RegresParamFound();
                    NachYslovRegAnal(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                    KoefDeterm = Math.Pow(Korelation[0], 2)*100;
                    Szal = SzalFound(AB, ML[MSelectGR[0]], ML[MSelectGR[1]]);
                    AB[2] = Szal * Math.Sqrt(1.0 / ML[MSelectGR[0]].l.Count + Math.Pow(ML[MSelectGR[0]].Mx.Q, 2) /
                        (ML[MSelectGR[0]].Dx.Q * (ML[MSelectGR[0]].l.Count - 1)));
                    AB[3] = Szal / (ML[MSelectGR[0]].Gx.Q * Math.Sqrt(ML[MSelectGR[0]].l.Count - 1));

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
                    Korelation[2] = op;
                    SravnSred[1] = Distributions.StudentQuantile(ML[MSelectGR[0]].alf.Q / 2,
                        ML[MSelectGR[0]].l.Count + ML[MSelectGR[1]].l.Count - 2);
                }
                else 
                {
                    SravnSred[1] = Distributions.StudentQuantile(ML[MSelectGR[0]].alf.Q / 2, ML[MSelectGR[0]].l.Count - 2);
                }

                SravnSred[0] = SimpleMx(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                SravnDisper[0] = SimpleS(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                SravnDisper[1] = Distributions.FisherQuantile(ML[MSelectGR[0]].alf.Q, ML[MSelectGR[0]].l.Count - 1, ML[MSelectGR[1]].l.Count - 1);
                
                KrSmirnKolmag[0] = 1.0 - LzFound(ZFound(ML[MSelectGR[0]], ML[MSelectGR[1]])); 
                KrSmirnKolmag[1] = ML[MSelectGR[0]].alf.Q;

                KrsumRangVils[0] = KrsumRangVilsFound(r,ML[MSelectGR[0]], ML[MSelectGR[1]]); 
                KrsumRangVils[1] = Distributions.NormalQuantile(1-ML[MSelectGR[0]].alf.Q / 2.0);

                KrUMannaUit[0] = KrUMannaUitFound(ML[MSelectGR[0]], ML[MSelectGR[1]]);
                KrUMannaUit[1] = Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2.0);

                RizSerRangVib[0] = RizSerRangVibFound(r,ML[MSelectGR[0]], ML[MSelectGR[1]]);
                RizSerRangVib[1] = Distributions.NormalQuantile(1 - ML[MSelectGR[0]].alf.Q / 2.0);

                KrZnakiv[0] = KrZnakivFound(ML[MSelectGR[0]], ML[MSelectGR[1]]);

            }
            else
            {
                Doubl = false;
                SravnSred[0] = SimpleMx(ML, MSelectGR);
                SravnDisper[0] = SimpleS(ML, MSelectGR);
                SravnDisper[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1);

                KrKruskalaUolisa[0] = KrKruskalaUolisaFound(r, ML, MSelectGR);
                KrKruskalaUolisa[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1);


            }
            Round();
        }
        private double[] RegresParamFound() 
        {
            double[] rez = new double[3];
            return rez;
        }

        private bool Nezalegni()
        {
            bool Nezal = true;
            double NezCount = ML[MSelectGR[0]].l.Count;
            for (int i = 0; i < MSelectGR.Count; i++)
                if (ML[MSelectGR[i]].l.Count != NezCount)
                    Nezal = false;
            return Nezal;
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
            AB[1] = Korelation[0] * gr2.Gx.Q / gr1.Gx.Q;
            AB[0] = gr2.Mx.Q - AB[1] * gr1.Mx.Q;
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
                     t1 += yi[i].Count * Math.Pow(ysr - AB[0] - AB[1] * ((i + 0.5) * gr1.Len.Q / M + gr1.Min.Q),2);
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
        double[,] NijFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR,double BinSim) 
        {
            double[,] rez = new double[3,3];
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++) 
            {
                if (ML[MSelectGR[0]].unsortl[i] > BinSim && ML[MSelectGR[1]].unsortl[i] > BinSim)
                    rez[1, 1]++;
                else if (ML[MSelectGR[0]].unsortl[i] > BinSim && ML[MSelectGR[1]].unsortl[i] <= BinSim)
                    rez[0, 1]++;
                else if (ML[MSelectGR[0]].unsortl[i] <= BinSim && ML[MSelectGR[1]].unsortl[i] > BinSim)
                    rez[1, 0]++;
                else if (ML[MSelectGR[0]].unsortl[i] <= BinSim && ML[MSelectGR[1]].unsortl[i] <= BinSim)
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

        double KorelationFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double Mxy = 0;
            double sx = 0, sy = 0;
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
            {
                sx += ML[MSelectGR[0]].unsortl[i];
                sy += ML[MSelectGR[1]].unsortl[i];
                Mxy += ML[MSelectGR[0]].unsortl[i] * ML[MSelectGR[1]].unsortl[i];
            }
            var ty = Mxy / Math.Sqrt(sx * sy);
            Mxy /= ML[MSelectGR[0]].l.Count;

            return (ML[MSelectGR[0]].l.Count / (ML[MSelectGR[0]].l.Count - 1)) *
                ((Mxy - ML[MSelectGR[0]].Mx.Q * ML[MSelectGR[1]].Mx.Q) / (ML[MSelectGR[0]].Gx.Q * ML[MSelectGR[1]].Gx.Q));
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
            for (int i=0;i<ML[MSelectGR[0]].l.Count ;i++ ) 
            {
                d_2 += Math.Pow(ri(r, ML[MSelectGR[0]].l[i]) - ri(r, ML[MSelectGR[1]].l[i]), 2);
            }
            return 1 -  d_2 / (ML[MSelectGR[0]].l.Count * (Math.Pow(ML[MSelectGR[0]].l.Count,2)-1));
        }
        double RangKoefKendFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double s = 0;
            double[] rangY = rDvomFound(ML[MSelectGR[1]].unsortl);
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
            return 2 * s/ (ML[MSelectGR[1]].l.Count * (ML[MSelectGR[1]].l.Count - 1));
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
        double SimpleMx(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            if (gr1.l.Count != gr2.l.Count)
            {
                if (gr1.l.Count + gr2.l.Count <=25)
                {
                    double A = gr1.Mx.Q - gr2.Mx.Q;
                    double B = Math.Sqrt((gr1.l.Count * gr2.l.Count) / (gr1.l.Count + gr2.l.Count));
                    double C = (gr1.l.Count - 1) * gr1.Dx.Q + (gr2.l.Count - 1) * gr2.Dx.Q;
                    C /= gr1.l.Count + gr2.l.Count - 2;
                    C = Math.Sqrt(C);
                    return A * B / C;
                }
                return (gr1.Mx.Q - gr2.Mx.Q) / Math.Sqrt((gr1.Dx.Q / gr1.l.Count) + (gr2.Dx.Q / gr2.l.Count));
            }
            for (int i = 0; i < gr1.l.Count;i++ )
                Z.Add(gr1.unsortl[i] - gr2.unsortl[i]);
            InitialStatisticalAnalys TG = new InitialStatisticalAnalys(Z);
            if (TG.Mx.Q == TG.Gx.Q && TG.Mx.Q == 0)
                return 0;
            return TG.Mx.Q * Math.Sqrt(TG.l.Count) / TG.Gx.Q;
        }
        double SimpleS(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            if (gr1.Dx.Q >= gr2.Dx.Q)
                return gr1.Dx.Q / gr2.Dx.Q;
            else
                return gr2.Dx.Q / gr1.Dx.Q;
        }
        double SimpleMx(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double S = 0, Mx = 0,hx = 0,hs = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
                hx += ML[MSelectGR[i]].l.Count * ML[MSelectGR[i]].Mx.Q;
            Mx =hx/N;
            for (int i = 0; i < MSelectGR.Count; i++)
                hs += ML[MSelectGR[i]].l.Count * Math.Pow(ML[MSelectGR[i]].Mx.Q - Mx, 2);
            S = hs / (MSelectGR.Count - 1);
            return S;
        }
        double SimpleS(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double ha = 0, hb = 0;
            for (int i = 0; i < MSelectGR.Count;i++ )
            {
                ha += (ML[MSelectGR[i]].l.Count - 1) * ML[MSelectGR[i]].Dx.Q;
                hb += ML[MSelectGR[i]].l.Count - 1;
            }
            double S = ha / hb;
            double B = 0;
            double C = 0;
            ha = hb = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
            {
                B -= (ML[MSelectGR[i]].l.Count - 1) * Math.Log(ML[MSelectGR[i]].Dx.Q/S);
                ha += 1 / (ML[MSelectGR[i]].l.Count - 1);
                hb += (ML[MSelectGR[i]].l.Count - 1);
            }
            hb = 1 / hb;
            C = 1 + (ha - hb) / (3 * MSelectGR.Count - 1);
            return B / C;
        }


        private double Sv2Found(List<InitialStatisticalAnalys> ML, List<int> MSelectGR) 
        {
            double hs = 0,Sv2 = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
                hs += (ML[MSelectGR[i]].l.Count - 1) * ML[MSelectGR[i]].Dx.Q;
            Sv2 = hs / (N - MSelectGR.Count);
            return Sv2;
        }
        private double Sm2Found(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double Sm2 = 0, xsr = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
                xsr += ML[MSelectGR[i]].l.Count * ML[MSelectGR[i]].Mx.Q;
            xsr /= N;
            for (int i = 0; i < MSelectGR.Count; i++)
                Sm2 += ML[MSelectGR[i]].l.Count * Math.Pow(ML[MSelectGR[i]].Mx.Q - xsr , 2);
            Sm2 /= (MSelectGR.Count - 1);
            return Sm2;
        }
#pragma warning disable IDE1006
        double[,] rFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
#pragma warning restore IDE1006
        {
            double[,] rez = new double[(int)N, 3];
            bool[] act = new bool[MSelectGR.Count];
            int[] indm = new int[MSelectGR.Count];
            for (int j = 0; j < N; )
            {
                int index = 0;
                double Min = double.MaxValue, gran = double.MaxValue;
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (act[i]==false && Min > ML[MSelectGR[i]].l[indm[i]])
                    {
                        index = i;
                        Min = ML[MSelectGR[i]].l[indm[i]];
                    }
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (i != index && act[i] == false && gran > ML[MSelectGR[i]].l[indm[i]])
                        gran = ML[MSelectGR[i]].l[indm[i]];
                for (; indm[index] < ML[MSelectGR[index]].l.Count&&ML[MSelectGR[index]].l[indm[index]] <= gran; indm[index]++,j++)
                {
                    rez[j,0] = j+1;
                    rez[j,1] = index;
                    rez[j, 2] = ML[MSelectGR[index]].l[indm[index]];
                }
                if (indm[index] >= ML[MSelectGR[index]].l.Count)
                    act[index] = true;
            }
            {
                List<double[]> simpldata = new List<double[]>();
                double prev = rez[0,2];
                simpldata.Add(new double[]{rez[0,0]-1,prev});
                for (int i = 1; i < N; prev = rez[i, 2],i++ )
                {
                    if (rez[i, 2] != prev)
                    {
                        if (simpldata.Count > 1)
                        {
                            double sum = 0;
                            for (int j = 0; j < simpldata.Count; sum += simpldata[j++][0])
                                ;
                            for (int j = 0; j < simpldata.Count; j++)
                                rez[i - j - 1, 0] = sum / simpldata.Count;
                        }
                        simpldata.Clear();
                    }
                    simpldata.Add(new double[] { rez[i, 0], rez[i, 2] });   
                }
            }
            return rez;
        }
#pragma warning disable IDE1006
        double[] rDvomFound(double[] mas)
#pragma warning restore IDE1006
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
        double KrsumRangVilsFound(double[,] rang,InitialStatisticalAnalys gr1,InitialStatisticalAnalys gr2)
        { 
            double KrsumRangVils = 0,Mw=0,Dw = 0;
            for (int i = 0; i < N;i++ )
                if (rang[i,1]==0)
                    KrsumRangVils+=rang[i,0];
            Mw = gr1.l.Count * (N + 1) / 2;
            Dw = gr1.l.Count * gr2.l.Count * (N + 1) / 12;
            return (KrsumRangVils - Mw) / Math.Sqrt(Dw);
        }
        double KrUMannaUitFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            double u = 0,Mu=0,Du=0;
            double[,] z = new double[gr1.l.Count, gr2.l.Count];
            for (int i = 0; i < gr1.l.Count; i++)
                for (int j = 0; j < gr2.l.Count; j++)
                    if (gr1.unsortl[i] > gr2.unsortl[j])
                    {
                        z[i, j] = 1;
                        u++;
                    }
                    else
                        z[i, j] = 0;
            Mu = gr1.l.Count * gr2.l.Count / 2;
            Du = Mu * (N + 1) / 6;
            return (u - Mu) / Math.Sqrt(Du);
        }
        double RizSerRangVibFound(double[,] rang, InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2) 
        {
            double Mrx = 0, Mry = 0;
            for(int i =0;i<N;i++)
                if (rang[i,1]==0)
                    Mrx+=rang[i,0];
                else if (rang[i,1]==1)
                    Mry+=rang[i,0];
            Mrx /= gr1.l.Count;
            Mry /= gr2.l.Count;
            double v = (Mrx - Mry) / (N * Math.Sqrt((N + 1) / (12.0 * gr1.l.Count * gr2.l.Count)));
            return v;
        }
        double KrZnakivFound(InitialStatisticalAnalys gr1,InitialStatisticalAnalys gr2){
            double S = 0;
            double hN=gr1.l.Count;
            List<double> Z = new List<double>();
            for (int i = 0; i < gr1.l.Count && i < gr2.l.Count;i++ ){
                Z.Add(gr1.unsortl[i] - gr2.unsortl[i]);
                if (gr1.unsortl[i] - gr2.unsortl[i] > 0)
                    S++;
                if (gr1.unsortl[i] - gr2.unsortl[i] == 0)
                    hN--;
            }
            if (hN<=15)
            {
                double alf0 = 0;
                for (int i = 0; i <= hN - S; i++)
                    alf0 += CFound(i, hN);
                alf0 *= Math.Pow(2, -hN);
                KrZnakiv[1] = gr1.alf.Q;
                return alf0;
            }
            else
            {
                double Sn = (S - 0.5 - hN / 2) / Math.Sqrt(hN / 4);
                KrZnakiv[1] = Distributions.NormalQuantile(1 - gr1.alf.Q/2);
                return Sn;
            }
        }
        double KrKruskalaUolisaFound(double[,] rang, List<InitialStatisticalAnalys> ML, List<int> MSelectGR) /// <=Hi alf v = MSelectGR.Count -1
        {
            double rez = 0;
            double[] wi = new double[MSelectGR.Count];
            for (int i = 0; i < MSelectGR.Count; i++)
            {
                for (int j = 0; j < N; j++)
                    if (rang[j, 1] == i)
                        wi[i] += rang[j, 0];
                wi[i] /= ML[MSelectGR[i]].l.Count;
                double A = Math.Pow(wi[i]-(N + 1)/2,2);
                double B = (N+1) * ( N - ML[MSelectGR[i]].l.Count)/(ML[MSelectGR[i]].l.Count*12);
                double C = 1 - ML[MSelectGR[i]].l.Count/N;
                rez += A * C / B;
            }
            
            return rez;
        }
#pragma warning disable IDE1006
        double[,] fFound(InitialStatisticalAnalys gr1,InitialStatisticalAnalys gr2,int Ncol,int Nrow)
#pragma warning restore IDE1006
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

#pragma warning disable IDE1006
        double ri(double[,] r,double value)
#pragma warning restore IDE1006
        {
            for (int i = 0; i < N; i++)
                if(r[i,2] == value){
                    return r[i, 0];
                }
            return -1;
        }
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
        double LzFound(double z)
        {
            double A = 1 - 2.0 * z / (3.0 * Math.Sqrt(N));
            double B = 2.0 * Math.Pow(z, 2) * (1 - 2.0 * Math.Pow(z, 2) / 3.0) / (3.0 * N);
            double C = 4.0 * z * (1.0 / 5.0 - 19.0 * Math.Pow(z, 2) / 15.0 + 2.0 * Math.Pow(z, 4) / 3.0) / (9.0 * Math.Pow(N, 1.5));
            return 1.0 - Math.Exp(-2.0 * Math.Pow(z, 2)) * (A + B + C);
        }
        private double SzalFound(double[] AB,InitialStatisticalAnalys X,InitialStatisticalAnalys Y)
        {
            double rez = 0;
            for (int i =0;i<Y.l.Count ;i++ ) 
            {
                rez += Math.Abs(Y.unsortl[i] - AB[0] - AB[1] * X.unsortl[i]);
            }
            rez /= X.l.Count;
            return rez;

        }
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


        private void Round()
        {
            int o=4;
            SravnSred[0] = Math.Round( SravnSred[0],o);
            SravnDisper[0] = Math.Round( SravnDisper[0],o);
            VarSv2Sm2[0] = Math.Round( VarSv2Sm2[0],o);
            KrSmirnKolmag[0] = Math.Round( KrSmirnKolmag[0],o);
            KrsumRangVils[0] = Math.Round( KrsumRangVils[0],o);
            KrUMannaUit[0] = Math.Round( KrUMannaUit[0],o);
            RizSerRangVib[0] = Math.Round( RizSerRangVib[0],o);
            KrZnakiv[0] = Math.Round(KrZnakiv[0], o);
            KrKruskalaUolisa[0] = Math.Round(KrKruskalaUolisa[0], o);

            SravnSred[1] = Math.Round(SravnSred[1], o);
            SravnDisper[1] = Math.Round(SravnDisper[1], o);
            VarSv2Sm2[1] = Math.Round(VarSv2Sm2[1], o);
            KrSmirnKolmag[1] = Math.Round(KrSmirnKolmag[1], o);
            KrsumRangVils[1] = Math.Round(KrsumRangVils[1], o);
            KrUMannaUit[1] = Math.Round(KrUMannaUit[1], o);
            RizSerRangVib[1] = Math.Round(RizSerRangVib[1], o);
            KrZnakiv[1] = Math.Round(KrZnakiv[1], o);
            KrKruskalaUolisa[1] = Math.Round(KrKruskalaUolisa[1], o);


            KrKohrena[0] = Math.Round(KrKohrena[0], o);
            KrKohrena[1] = Math.Round(KrKohrena[1], o);

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

            Sm2 = Math.Round(Sm2, o);
            Sv2 = Math.Round(Sv2, o);
            KriterBarkleta[0] = Math.Round(KriterBarkleta[0], o);
            KriterBarkleta[1] = Math.Round(KriterBarkleta[1], o);

            KoefDeterm = Math.Round(KoefDeterm, 0);

            AB[0] = Math.Round(AB[0], o);
            AB[1] = Math.Round(AB[1], o);
            AB[2] = Math.Round(AB[2], o);
            AB[3] = Math.Round(AB[3], o);

            ABTeil[0] = Math.Round(ABTeil[0], o);
            ABTeil[1] = Math.Round(ABTeil[1], o);

            ProvRegrs[0] = Math.Round(ProvRegrs[0], o);
            ProvRegrs[1] = Math.Round(ProvRegrs[1], o);

            X2f[0] = Math.Round(X2f[0], o);
            X2f[1] = Math.Round(X2f[1], o);
        }
    }
}
