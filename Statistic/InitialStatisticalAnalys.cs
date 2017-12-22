using System;
using System.Collections.Generic;
using System.Windows.Forms;

namespace testgistogr
{
    public class InitialStatisticalAnalys
    {
        public List<double> l;
        public List<Data> Estimation = new List<Data>();
        public double[] unsortl;
        public Data m ;
        public Data alf;
        public Data Min;
        public Data Max;
        public Data Len;
        public Data Step;
        public Data Mx;
        public Data X_2;
        public Data Mediana;
        public Data MAD;
        public Data Moda;
        public Data Gx;
        public Data Dx;
        public Data Mx_rang;
        public Data AverYolsha;
        public Data CoefAsim;
        public Data CoefVarPirson;
        public Data CoefEcscecness;
        public Data CoefEcscec;
        public Data W;
        public double[] Quantile;
        public List<double> f;
        public List<double> Y2;
        public List<double> F;
        public double a;
        public double b;
        public double A;
        public double B;
        public Data predictioninterval;
        public Data KrAbbe;


        public int Group = 0;

        public double T = 0;

        public double[] TypeR;
        public string AvtoType;


        public InitialStatisticalAnalys(List<double> l)
        {
            DataClassified();
            this.l = new List<double>(l);
            unsortl = new double[l.Count];
            l.CopyTo(unsortl);
            this.l.Sort();
            this.alf.Q = 0.05;
            if (l.Count <= 100)
            {
                var v1 = (int)Math.Sqrt(l.Count);
                if (v1 % 2 == 0)
                    this.m.Q = v1 - 1;
                else
                    this.m.Q = v1;
            }
            else if (l.Count > 100)
            {
                var v1 = (int)Math.Pow(l.Count, 1.0 / 3.0);
                if (v1 % 2 == 0)
                    this.m.Q = v1 - 1;
                else
                    this.m.Q = v1;
            }
            Refresh();
            //TypeR = TypeRFound();
        }
        public InitialStatisticalAnalys(List<double> l,int Group)
        {
            DataClassified();
            this.l = new List<double>(l);
            this.Group = Group;
            unsortl = new double[l.Count];
            l.CopyTo(unsortl);
            this.l.Sort();
            this.alf.Q = 0.05;
            if (l.Count <= 100)
            {
                var v1 = (int)Math.Sqrt(l.Count);
                if (v1 % 2 == 0)
                    this.m.Q = v1 - 1;
                else
                    this.m.Q = v1;
            }
            else if (l.Count > 100)
            {
                var v1 = (int)Math.Pow(l.Count, 1.0 / 3.0);
                if (v1 % 2 == 0)
                    this.m.Q = v1 - 1;
                else
                    this.m.Q = v1;
            }
            Refresh();
            TypeR = TypeRFound();
        }
        public InitialStatisticalAnalys(List<double> l, int m,double alf)
        {
            DataClassified();
            this.l = l;
            unsortl = new double[l.Count];
            l.CopyTo(unsortl);
            this.l.Sort();
            this.m.Q = m;
            this.alf.Q = alf;
            Refresh();
            TypeR = TypeRFound();
        }

        public double this[int index]
        {
            get
            {
                return unsortl[index];
            }
        }

        public double Length()
        {
            return unsortl.Length;
        }

        private void DataClassified() 
        {
            this.m = new Data()
            {
                Name = "Кількість класів"
            };
            Estimation.Add(this.m);

            this.alf = new Data()
            {
                Name = "Aльфа"
            };
            Estimation.Add(this.alf);

            Min = new Data()///minimal value
            {
                Name = "Мінімальне значення"
            };
            Estimation.Add(Min);

            Max = new Data()///maximum value
            {
                Name = "Максимальне значення"
            };
            Estimation.Add(Max);

            Len = new Data()///minimal value
            {
                Name = "Довжина"
            };
            Estimation.Add(Len);

            Step = new Data()///step data
            {
                Name = "Крок"
            };
            Estimation.Add(Step);

            Mx = new Data()///average value
            {
                Name = "Середнье арифметичне"
            };
            Estimation.Add(Mx);

            Mx_rang = new Data()///averag-mi-max /rangovana mx
            {
                Name = "Середнье арифметичне ранжоване"
            };
            Estimation.Add(Mx_rang);

            Gx = new Data()///S
            {
                Name = "Сігма(S)"
            };
            Estimation.Add(Gx);

            Dx = new Data()///S^2
            {
                Name = "Дисперсія(S^2)"
            };
            Estimation.Add(Dx);

            X_2 = new Data()///average value x^2
            {
                Name = "Середнье арифметичне x^2"
            };
            //Estimation.Add(X_2);

            Mediana = new Data()
            {
                Name = "Вибіркова медіана(MEDIANA)"
            };
            Estimation.Add(Mediana);

            MAD = new Data()///median absolute value
            {
                Name = "MAD"
            };
            Estimation.Add(MAD);

            Moda = new Data()
            {
                Name = "Мода(MODA)"
            };
            Estimation.Add(Moda);

            AverYolsha = new Data()
            {
                Name = "Медіана середніх Уолша"
            };
            Estimation.Add(AverYolsha);

            CoefAsim = new Data()
            {
                Name = "Коефіціент асіметріі(А)"
            };
            Estimation.Add(CoefAsim);

            CoefVarPirson = new Data()
            {
                Name = "Коефіціент віріації Пірсона(W)"
            };
            Estimation.Add(CoefVarPirson);

            CoefEcscecness = new Data()
            {
                Name = "Коефіціент ексцесу(E) незсунений"
            };
            Estimation.Add(CoefEcscecness);

            CoefEcscec = new Data()
            {
                Name = "Коефіціент ексцесу(E) зсунений"
            };
            //Estimation.Add(CoefEcscec);

            W = new Data()
            {
                Name = "Коефіціент ексцесу(E) незсунений"
            };
            Estimation.Add(W);

            predictioninterval = new Data()
            {
                Name = "Довірчий інтервал"
            };
            //Estimation.Add(predictioninterval);

            KrAbbe = new Data()
            {
                Name = "Кр.Аббе"
            };
            Estimation.Add(KrAbbe);

        }
        public void Refresh() 
        {
            Min.Q = l[0];
            Max.Q = l[l.Count - 1];
            Len.Q = Max.Q - Min.Q;
            T = Distributions.StudentQuantile(alf.Q / 2, (int)m.Q - 1);
            f = new List<double>();
            Y2 = new List<double>();
            F = new List<double>();
            if (l[l.Count - 1] - l[0]==0)
                Step.Q = 1;
            else
                Step.Q = 1.0001 * ((l[l.Count - 1] - l[0]) / m.Q);
            MxFound(l);
            DxFound();
            X_2.Q = StartMoment(l, 2);
            Mx_rang.Q = Mx.Q - (l[0] + l[l.Count - 1]) / l.Count;
            Mediana.Q = MEDFound(l);
            MADFound();
            MODFound();
            CoefEcscecFound(l);
            CoefAsimFound(l);
            QuantileFound();
            //AverYolshaFound();
            CutData();
            CoefVarPirson.Q = Gx.Q / Mx.Q;
            W.Q = MAD.Q / Mediana.Q;
            PredictionintervalFound(l);
            A = Mx.Q - Math.Sqrt(3 * (X_2.Q - Mx.Q * Mx.Q));
            B = Mx.Q + Math.Sqrt(3 * (X_2.Q - Mx.Q * Mx.Q));
            KrAbbe.Q = KrAbbeFound(l, unsortl);
            KrAbbe.QKvant = Distributions.NormalQuantile(1 - alf.Q / 2);

            DataRound();
        }

        double KrAbbeFound(List<double> l, double[] unsortl)
        {
            double d2 = 0;
            for (int i = 0; i < l.Count - 1; i++)
                d2 += Math.Pow(unsortl[i + 1] - unsortl[i], 2) / (l.Count - 1);
            double q = d2 / (2.0 * Dx.Q);
            double U = (q - 1) * Math.Sqrt((Math.Pow(l.Count, 2) - 1) / (l.Count - 2));
            return U;
        }
        private void MxFound(List<double> l)
        {
            Mx.Q = StartMoment(l, 1);
        }
        public static double MEDFound(List<double> p)
        {
            p.Sort();
            double rez;
            if (p.Count % 2 == 0)
                rez = (p[p.Count / 2 - 1] + p[p.Count / 2]) / 2;
            else
                rez = p[(p.Count - 1) / 2];
            return rez;
        }
        void MADFound()
        {
            List<double> hlplist = new List<double>();
            for (int i = 0; i < l.Count; i++)
            {
                hlplist.Add(Math.Abs(l[i] - Mediana.Q));
            }
            hlplist.Sort();
            MAD.Q = Math.Round(1.483 * MEDFound(hlplist), 5);
        }
        void MODFound()
        {
            double max = 0;
            double j = 0;
            for (double i = l[0]; i <= l[l.Count - 1]; i += Step.Q)
            {
                int N = 0, N2 = 0;
                for (int p = 0; p < l.Count && l[p] <= i + Step.Q; p++, N2++)
                {
                    if (l[p] >= i && l[p] < i + Step.Q)
                    {
                        N++;
                    }
                }
                if (max <= N)
                {
                    max = N;
                    j = i;
                }
                f.Add((double)N / l.Count);
                Y2.Add((double)N2 / l.Count);
            }
            Moda.Q = j + Step.Q / 2;
            for (int i = 0; i < l.Count; i++)
            {
                int i2 = i+1;
                for (; i2<l.Count&&l[i] == l[i2]; i2++)
                    ;
                F.Add((double)i2 / l.Count);
            }
            double Xd = (l[l.Count - 1] - l[0]) / 1000;
        }
        void DxFound()
        {
            for (int i = 0; i < l.Count; i++)
            {
                Dx.Q += Math.Pow((l[i] - Mx.Q), 2);
            }
            Dx.Q = Math.Round((Dx.Q / (l.Count - 1)), 5);
            Gx.Q = Math.Round(Math.Sqrt(Dx.Q), 5);
            Gx.QSigma=Gx.Q / Math.Sqrt(2*l.Count);
            Dx.QSigma=Gx.QSigma * Gx.QSigma;
            Gx.QButton = Gx.Q - T * Gx.QSigma;
            Gx.QUpper = Gx.Q + T * Gx.QSigma;
            Dx.QButton = Gx.QButton * Gx.QButton;
            Dx.QUpper = Gx.QUpper * Gx.QUpper;


            Mx.QSigma = Gx.Q / Math.Sqrt(l.Count);
            Mx.QButton = Mx.Q - T * Mx.QSigma;
            Mx.QUpper = Mx.Q + T * Mx.QSigma;
            //Mx_G = Gx / Math.Sqrt(l.Count);
            //Mx_low = Mx.Q - T * Mx_G;
            //Mx_upper = Mx.Q + T * Mx_G;
        }
        void QuantileFound()
        {
            Quantile = new double[7];
            Quantile[0] = l[(int)(l.Count * 0.05)];
            Quantile[1] = l[(int)(l.Count * 0.1)];
            Quantile[2] = l[(int)(l.Count * 0.25)];
            Quantile[3] = l[(int)(l.Count * 0.5)];
            Quantile[4] = l[(int)(l.Count * 0.75)];
            Quantile[5] = l[(int)(l.Count * 0.9)];
            Quantile[6] = l[(int)(l.Count * 0.95)];
        }
        void AverYolshaFound()
        {
            List<double> lhlp=new List<double>();
            for(int i = 0;i<l.Count;i++)
            {
                for (int j = 0; j <= i;j++ )
                {
                    try
                    {
                        lhlp.Add((l[i] + l[j]) / 2);
                    }
                    catch
                    {
                        lhlp.Clear();
                        MessageBox.Show("Недостатньо пам'яті для знаходження середнього Уолша");
                        AverYolsha.Q = 0;
                        return;
                    }
                }
            }
            AverYolsha.Q = MEDFound(lhlp);
        }
        void CutData()
        {

            double t1 = 2 + 0.21 * Math.Log10(0.04 * l.Count);
            double t2 = Math.Sqrt(19*Math.Sqrt(CoefEcscec.Q+2)+ 1 );
            //double t = 1.2 + 3.6 * (1 - 1 / Math.Sqrt(Math.Abs(CoefEcscec))) * Math.Log10(l.Count / 10);
            //double t = 1.55 + 0.8 * Math.Sqrt(Math.Abs(CoefEcscecness - 1)) * Math.Log10(l.Count / 10);
            //a = Mx - t * Gx;
            //b = Mx + t * Gx;
            if (CoefAsim.Q>0.2)
            {
                a = Mx.Q - t1 * Gx.Q;
                b = Mx.Q + t2 * Gx.Q;
            }
            else if (CoefAsim.Q < -0.2)
            {
                a = Mx.Q - t2 * Gx.Q;
                b = Mx.Q + t1 * Gx.Q;
            }
            else
            {
                a = Mx.Q - t1 * Gx.Q;
                b = Mx.Q + t1 * Gx.Q;
            }
        }
        static public double StartMoment(List<double> p, int N)
        {
            double rez = 0;
            for (int i = 0; i < p.Count;i++ )
            {
                rez += Math.Pow(p[i], N);
            }
            rez/=p.Count;
            return rez;
        }
        public double CentralMoment(List<double> p, int N)
        {
            double rez = 0;
            for (int i = 0; i < p.Count;i++ )
            {
                rez+=Math.Pow((p[i]-Mx.Q),N);
            }
            rez/=p.Count;
            return (double)rez;
        }
        static public List<double> StandData(double[] unsortl, double Gx, double Ex)
        {
            List<double> rez = new List<double>();
            for (int i =0 ;i<unsortl.Length;i++)
                rez.Add((unsortl[i] - Ex) / Gx);
            return rez;
        }
        private void PredictionintervalFound(List<double> l)
        {
            predictioninterval.QSigma = Gx.Q * Math.Sqrt(1 + 1 / l.Count);
            predictioninterval.QButton = Mx.Q - T * predictioninterval.QSigma;
            predictioninterval.QUpper = Mx.Q + T * predictioninterval.QSigma;
        }
        private double Bi(int index)
        {
            double rez = 0;
            if (index % 2 == 0)
            {
                var v1 = CentralMoment(l, index + 2);
                var v2 = Math.Pow(CentralMoment(l, 2), (index / 2) + 1);
                rez = v1 / v2;
            }
            else
            {
                var v1 = CentralMoment(l, 3) * CentralMoment(l, index + 2);
                var v2 = Math.Pow(CentralMoment(l, 2), (index - 1 / 2) + 3);
                rez = v1 / v2;
            }
            return rez;
        }
        private void CoefAsimFound(List<double> l)
        {
            CoefAsim.Q = CentralMoment(l, 3) / Math.Pow(CentralMoment(l, 2), 1.5);
            CoefAsim.Q *= Math.Sqrt(l.Count * (l.Count - 1)) / (l.Count - 2);
            CoefAsim.QSigma = (Math.Sqrt(Math.Abs((4 * Bi(4) - 12 * Bi(3) - 24 * Bi(2) + 9 * Bi(2) * Bi(1) + 35 * Bi(1) - 36) / (4 * l.Count))));
            CoefAsim.QButton = CoefAsim.Q - T * CoefAsim.QSigma;
            CoefAsim.QUpper = CoefAsim.Q + T * CoefAsim.QSigma;
        }
        private void CoefEcscecFound(List<double> l)
        {
            var v1 = CentralMoment(l, 4);
            var v2 = Math.Pow(CentralMoment(l, 2), 2);
            double Es = v1 / v2;
            // Es /= l.Count;
            CoefEcscecness.Q = Es;
            int N = l.Count;
            CoefEcscec.Q = ((N * N - 1) / ((N - 2) * (N - 3))) * ((Es - 3) + 6 / (N + 1));
            CoefEcscec.QSigma = (Math.Sqrt(Math.Abs(Bi(6) - 4 * Bi(4) * Bi(2) - 8 * Bi(3) + 4 * Math.Pow(Bi(2), 3) - Math.Pow(Bi(2), 2) + 16 * Bi(2) * Bi(1) + 16 * Bi(1)) / (l.Count)));
            CoefEcscec.QButton = CoefEcscec.Q - T * CoefEcscec.QSigma;
            CoefEcscec.QUpper = CoefEcscec.Q + T * CoefEcscec.QSigma;
        }
        private void DataRound() 
        {
            int o=4;
            Step.Q=Math.Round(Step.Q,o);
            Mx.Q=Math.Round(Mx.Q,o);
            Mediana.Q=Math.Round(Mediana.Q,o);
            MAD.Q=Math.Round(MAD.Q,o);
            Moda.Q=Math.Round(Moda.Q,o);
            Gx.Q=Math.Round(Gx.Q,o);
            Dx.Q=Math.Round(Dx.Q,o);
            Mx_rang.Q=Math.Round(Mx_rang.Q,o);
            AverYolsha.Q=Math.Round(AverYolsha.Q,o);
            CoefAsim.Q=Math.Round(CoefAsim.Q,o);
            CoefVarPirson.Q=Math.Round(CoefVarPirson.Q,o);
            CoefEcscecness.Q=Math.Round(CoefEcscecness.Q,o);
            CoefEcscec.Q = Math.Round(CoefEcscec.Q, o);
            W.Q = Math.Round(W.Q, o);
            for (int i = 0; i < Quantile.Length;i++ )
            {
                Quantile[i] = Math.Round(Quantile[i], o);
            }
            a = Math.Round(a, o);
            b = Math.Round(b, o);


            Mx.QButton = Math.Round(Mx.QButton, o);
            Mx.QUpper = Math.Round(Mx.QUpper, o);
            Gx.QButton = Math.Round(Gx.QButton, o);
            Gx.QUpper = Math.Round(Gx.QUpper, o);
            Dx.QButton = Math.Round(Dx.QButton, o);
            Dx.QUpper = Math.Round(Dx.QUpper, o);
            CoefAsim.QButton = Math.Round(CoefAsim.QButton, o);
            CoefAsim.QUpper = Math.Round(CoefAsim.QUpper, o);
            CoefEcscec.QButton = Math.Round(CoefEcscec.QButton, o);
            CoefEcscec.QUpper = Math.Round(CoefEcscec.QUpper, o);
            predictioninterval.QButton = Math.Round(predictioninterval.QButton, o);
            predictioninterval.QUpper = Math.Round(predictioninterval.QUpper, o);

            KrAbbe.Q = Math.Round(KrAbbe.Q, o);
            KrAbbe.QKvant = Math.Round(KrAbbe.QKvant, o);
        }
        
        private double[] TypeRFound() 
        {
            double[] TypeR = new double[3];
            double i2 = -0.5;
            for (int i = 0; i < l.Count - 1; i++, i2 += (double)1 / l.Count)
            {
                double v1 = Distributions.NormalFobrFound(i2);
                double Y_X = ((this.l[i] - Min.Q) / Len.Q);
                //double l = Math.Abs((ML[i] - ML[0]) / ((ML[ML.Count - 1] - ML[0])) - ((v1 * gr.Gx + gr.Mx.Q - ML[0]) / (ML[ML.Count - 1] - ML[0])));
                double li = Math.Abs((l[i] - v1 * Gx.Q - Mx.Q) / Len.Q);
                TypeR[0] += Math.Pow(li, 2);
                TypeR[1] += Math.Pow((this.F[i]) - (1 - Math.Pow(2.73, -(1 / this.Mx.Q) * (this.l[i] - this.Min.Q))), 2);
                TypeR[2] += Math.Pow((this.F[i]) - Y_X, 2);
            }
            TypeR[0] /= (this.l.Count - 1);
            TypeR[1] /= (this.l.Count - 1);
            TypeR[2] /= (this.l.Count - 1);
            if (TypeR[0] < TypeR[1] && TypeR[0] < TypeR[2])
                AvtoType = "Нормальний";
            else if (TypeR[1] < TypeR[0] && TypeR[1] < TypeR[2])
                AvtoType = "Експоненціальний";
            else if (TypeR[2] < TypeR[1] && TypeR[2] < TypeR[0])
                AvtoType = "Рівномірний";
            return TypeR;
        }
    }
}
