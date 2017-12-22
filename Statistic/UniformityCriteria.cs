using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace testgistogr
{

    class UniformityCriteria
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
        public double Sv2 = 0;
        public double Sm2 = 0;
        public double[] KrSmirnKolmag = new double[2];///krt>alf
        public double[] KrsumRangVils = new double[2];///krt<=u alf/2
        public double[] KrUMannaUit = new double[2];///krt<=u alf
        public double[] RizSerRangVib = new double[2];///krt<=u alf
        public double[] KrZnakiv = new double[2];
        public double[] KrKruskalaUolisa = new double[2];
        public double[] KrKohrena = new double[2];

        public double T = 0;
        private List<double> Z = new List<double>();
        private double[,] r;


        public UniformityCriteria(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            this.ML = ML;
            this.MSelectGR = MSelectGR;
            
            N = NFound(ML, MSelectGR);
            {
                double NezCount = ML[MSelectGR[0]].l.Count;
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (ML[MSelectGR[i]].l.Count != NezCount)
                        Nezal = false;
            }
            rDvomFound(ML[MSelectGR[0]].unsortl);
            r = rFound(ML, MSelectGR);
            Sm2 = Sm2Found(ML, MSelectGR);
            Sv2 = Sv2Found(ML, MSelectGR);
            VarSv2Sm2[0] = Sm2 / Sv2 ;
            VarSv2Sm2[1] = Distributions.FisherQuantile(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1, (int)(N - MSelectGR.Count));
            var tyui = Distributions.FisherQuantile(ML[MSelectGR[0]].alf.Q, 10,10);
            if (Nezal == true)
            {
                KrKohrena[0] = KrKohrenaFound(ML,MSelectGR);
                KrKohrena[1] = Hi.HIF(ML[MSelectGR[0]].alf.Q, MSelectGR.Count - 1);
                T = Distributions.StudentQuantile(1-ML[MSelectGR[0]].alf.Q / 2, ML[MSelectGR[0]].l.Count - 2);
            }
            if (MSelectGR.Count == 2)
            {
                Doubl = true;
                if (Nezal == true)
                {
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
        double KrKohrenaFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double Q = 0;
            double[] u = new double[ML[MSelectGR[0]].l.Count];
            double[] t = new double[MSelectGR.Count];
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
                u[i] = MSelectGR.Count;
            for (int i = 0; i < MSelectGR.Count; i++)
                t[i] = ML[MSelectGR[0]].l.Count;
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
                for (int j = 0; j < MSelectGR.Count; j++)
                    if (ML[MSelectGR[j]].unsortl[i] <= ML[MSelectGR[j]].Mx.Q)
                    {
                        t[j]--;
                        u[i]--;
                    }
            double mt = 0, mu = 0;
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
                mu += u[i] / ML[MSelectGR[0]].l.Count;
            for (int i = 0; i < MSelectGR.Count; i++)
                mt += t[i] / MSelectGR.Count;
            double eu = u.Sum(), et = 0, euu = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
                et += Math.Pow(t[i] - mt, 2);
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
                euu += Math.Pow(u[i], 2);
            if (MSelectGR.Count * eu - euu == 0)
            {
                return -1;
            }
            Q = MSelectGR.Count * (MSelectGR.Count - 1) * et / (MSelectGR.Count * eu - euu);

            return Q;
        }
        double KorelationVidnohFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double rez = 0;
            int M = 10;
            List<double>[] yi = new List<double>[M];
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
                t1 /= yi[i].Count;
                Emyjy += yi[i].Count * Math.Pow(yi[i].Sum() / yi[i].Count - ML[MSelectGR[1]].Mx.Q, 2);
            }
            rez = Emyjy / (ML[MSelectGR[1]].Dx.Q * (ML[MSelectGR[1]].l.Count - 1));
            rez = Math.Sqrt(rez);
            return rez;
        }
        double RangKorelationFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double d_2 = 0;
            for (int i = 0; i < ML[MSelectGR[0]].l.Count; i++)
            {
                d_2 += Math.Pow(ri(r, ML[MSelectGR[0]].l[i]) - ri(r, ML[MSelectGR[1]].l[i]), 2);
            }
            return 1 - d_2 / (ML[MSelectGR[0]].l.Count * (Math.Pow(ML[MSelectGR[0]].l.Count, 2) - 1));
        }
        double RangKoefKendFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double s = 0;
            double[] rangY = rDvomFound(ML[MSelectGR[1]].unsortl);
            for (int i = 0; i < ML[MSelectGR[0]].l.Count - 1; i++)
            {
                for (int j = i + 1; j < ML[MSelectGR[0]].l.Count; j++)
                {
                    if (rangY[i] < rangY[j])
                        s++;
                    else if (rangY[i] < rangY[j])
                        s--;
                }
            }
            return 2 * s / (ML[MSelectGR[1]].l.Count * (ML[MSelectGR[1]].l.Count - 1));
        }
        double NFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double hn = 0;
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
                if (gr1.l.Count + gr2.l.Count <= 25)
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
            for (int i = 0; i < gr1.l.Count; i++)
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
            double S = 0, Mx = 0, hx = 0, hs = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
                hx += ML[MSelectGR[i]].l.Count * ML[MSelectGR[i]].Mx.Q;
            Mx = hx / N;
            for (int i = 0; i < MSelectGR.Count; i++)
                hs += ML[MSelectGR[i]].l.Count * Math.Pow(ML[MSelectGR[i]].Mx.Q - Mx, 2);
            S = hs / (MSelectGR.Count - 1);
            return S;
        }
        double SimpleS(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double ha = 0, hb = 0;
            for (int i = 0; i < MSelectGR.Count; i++)
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
                B -= (ML[MSelectGR[i]].l.Count - 1) * Math.Log(ML[MSelectGR[i]].Dx.Q / S);
                ha += 1 / (ML[MSelectGR[i]].l.Count - 1);
                hb += (ML[MSelectGR[i]].l.Count - 1);
            }
            hb = 1 / hb;
            C = 1 + (ha - hb) / (3 * MSelectGR.Count - 1);
            return B / C;
        }


        private double Sv2Found(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double hs = 0, Sv2 = 0;
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
                Sm2 += ML[MSelectGR[i]].l.Count * Math.Pow(ML[MSelectGR[i]].Mx.Q - xsr, 2);
            Sm2 /= (MSelectGR.Count - 1);
            return Sm2;
        }
        double[,] rFound(List<InitialStatisticalAnalys> ML, List<int> MSelectGR)
        {
            double[,] rez = new double[(int)N, 3];
            bool[] act = new bool[MSelectGR.Count];
            int[] indm = new int[MSelectGR.Count];
            for (int j = 0; j < N; )
            {
                int index = 0;
                double Min = double.MaxValue, gran = double.MaxValue;
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (act[i] == false && Min > ML[MSelectGR[i]].l[indm[i]])
                    {
                        index = i;
                        Min = ML[MSelectGR[i]].l[indm[i]];
                    }
                for (int i = 0; i < MSelectGR.Count; i++)
                    if (i != index && act[i] == false && gran > ML[MSelectGR[i]].l[indm[i]])
                        gran = ML[MSelectGR[i]].l[indm[i]];
                for (; indm[index] < ML[MSelectGR[index]].l.Count && ML[MSelectGR[index]].l[indm[index]] <= gran; indm[index]++, j++)
                {
                    rez[j, 0] = j + 1;
                    rez[j, 1] = index;
                    rez[j, 2] = ML[MSelectGR[index]].l[indm[index]];
                }
                if (indm[index] >= ML[MSelectGR[index]].l.Count)
                    act[index] = true;
            }
            {
                List<double[]> simpldata = new List<double[]>();
                double prev = rez[0, 2];
                simpldata.Add(new double[] { rez[0, 0] - 1, prev });
                for (int i = 1; i < N; prev = rez[i, 2], i++)
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
        double[] rDvomFound(double[] mas)
        {
            double[] rez = new double[mas.Length];
            List<double[]> dat = new List<double[]>();
            for (int i = 0; i < mas.Length; i++)
                dat.Add(new double[] { mas[i], i, 0 });
            dat.Sort(delegate(double[] t1, double[] t2)
            { return (t1[0].CompareTo(t2[0])); });
            for (int i = 0; i < dat.Count; i++)
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
                rez[i] = dat[i][2];
            return rez;
        }
        double KrsumRangVilsFound(double[,] rang, InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            double KrsumRangVils = 0, Mw = 0, Dw = 0;
            for (int i = 0; i < N; i++)
                if (rang[i, 1] == 0)
                    KrsumRangVils += rang[i, 0];
            Mw = gr1.l.Count * (N + 1) / 2;
            Dw = gr1.l.Count * gr2.l.Count * (N + 1) / 12;
            return (KrsumRangVils - Mw) / Math.Sqrt(Dw);
        }
        double KrUMannaUitFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            double u = 0, Mu = 0, Du = 0;
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
            for (int i = 0; i < N; i++)
                if (rang[i, 1] == 0)
                    Mrx += rang[i, 0];
                else if (rang[i, 1] == 1)
                    Mry += rang[i, 0];
            Mrx /= gr1.l.Count;
            Mry /= gr2.l.Count;
            double v = (Mrx - Mry) / (N * Math.Sqrt((N + 1) / (12.0 * gr1.l.Count * gr2.l.Count)));
            return v;
        }
        double KrZnakivFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            double S = 0;
            double hN = gr1.l.Count;
            List<double> Z = new List<double>();
            for (int i = 0; i < gr1.l.Count && i < gr2.l.Count; i++)
            {
                Z.Add(gr1.unsortl[i] - gr2.unsortl[i]);
                if (gr1.unsortl[i] - gr2.unsortl[i] > 0)
                    S++;
                if (gr1.unsortl[i] - gr2.unsortl[i] == 0)
                    hN--;
            }
            if (hN <= 15)
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
                KrZnakiv[1] = Distributions.NormalQuantile(1 - gr1.alf.Q / 2);
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
                double A = Math.Pow(wi[i] - (N + 1) / 2, 2);
                double B = (N + 1) * (N - ML[MSelectGR[i]].l.Count) / (ML[MSelectGR[i]].l.Count * 12);
                double C = 1 - ML[MSelectGR[i]].l.Count / N;
                rez += A * C / B;
            }

            return rez;
        }
        double[,] fFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2, int Ncol, int Nrow)
        {
            double[,] rez = new double[Ncol, Nrow];
            for (int i = 0; i < gr1.l.Count; i++)
            {
                for (int j = 0; j < Ncol; j++)
                    for (int k = 0; k < Nrow; k++)
                        if (gr1.unsortl[i] <= (j + 1) * gr1.Len.Q / Ncol + gr1.Min.Q && gr1.unsortl[i] >= (j * gr1.Len.Q / Ncol + gr1.Min.Q) &&
                            gr2.unsortl[i] <= ((k + 1) * gr2.Len.Q / Nrow + gr2.Min.Q) && gr2.unsortl[i] >= (k * gr2.Len.Q / Nrow + gr2.Min.Q))
                        {
                            rez[j, k] += 1.0 / gr1.l.Count;
                            goto t1;
                        }
            t1: ;
            }
            /*for (int j = 0; j < Ncol; j++)
                for (int k = 0; k < Nrow; k++)
                    if (rez[j, k] <= gr1.alf)
                        rez[j, k] = 0;*/
            return rez;
        }

        double ri(double[,] r, double value)
        {
            for (int i = 0; i < N; i++)
                if (r[i, 2] == value)
                {
                    return r[i, 0];
                }
            return -1;
        }
        double CFound(double i, double Nd)
        {
            double iF = 1, NF = 1, NiF = 1;
            for (int j = 1; j <= i; j++)
                iF *= j;
            for (int j = 1; j <= Nd; j++)
                NF *= j;
            for (int j = 1; j <= Nd - i; j++)
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
        private double SzalFound(double[] AB, InitialStatisticalAnalys X, InitialStatisticalAnalys Y)
        {
            double rez = 0;
            for (int i = 0; i < Y.l.Count; i++)
            {
                rez += Math.Abs(Y.unsortl[i] - AB[0] - AB[1] * X.unsortl[i]);
            }
            rez /= X.l.Count;
            return rez;

        }
        private double ZFound(InitialStatisticalAnalys gr1, InitialStatisticalAnalys gr2)
        {
            double D = 0;
            for (int i = 0; i < gr1.l.Count; i++)
            {
                double l2 = 0;
                if (gr2.AvtoType == Distributions.Normal)
                    l2 = Distributions.NormalfFound(gr1.l[i], gr2.Mx.Q, gr2.Gx.Q);
                else if (gr2.AvtoType == Distributions.Exp)
                    l2 = Distributions.ExpFound(gr1.l[i], gr2.Mx.Q);
                else if (gr2.AvtoType == Distributions.Line)
                    l2 = (gr1.l[i] - gr2.Min.Q) * gr2.Len.Q;

                D = Math.Max(Math.Abs(gr1.F[i] - l2), D);
            }
            return D;
        }


        private void Round()
        {
            int o = 4;
            SravnSred[0] = Math.Round(SravnSred[0], o);
            SravnDisper[0] = Math.Round(SravnDisper[0], o);
            VarSv2Sm2[0] = Math.Round(VarSv2Sm2[0], o);
            KrSmirnKolmag[0] = Math.Round(KrSmirnKolmag[0], o);
            KrsumRangVils[0] = Math.Round(KrsumRangVils[0], o);
            KrUMannaUit[0] = Math.Round(KrUMannaUit[0], o);
            RizSerRangVib[0] = Math.Round(RizSerRangVib[0], o);
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

            Sm2 = Math.Round(Sm2, o);
            Sv2 = Math.Round(Sv2, o);
        }
    }
}
