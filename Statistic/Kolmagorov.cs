using System;
using System.Collections.Generic;

namespace testgistogr
{
    static class Kolmagorov
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ML"></param>
        /// <param name="gr"></param>
        /// <param name="type"></param>
        /// <returns>P</returns>
        static public double KolmagorovFound(List<double> ML, InitialStatisticalAnalys gr, int type)
        {
            return KolmagorovFound(ML, gr, type, gr.Mx.Q, gr.Gx.Q);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="ML"></param>
        /// <param name="gr"></param>
        /// <param name="type"></param>
        /// <returns>P</returns>
        static public double KolmagorovFound(List<double> ML, InitialStatisticalAnalys gr, int type, double Mx, double Gx)
        {
            double D  = DFound(ML, gr, type,Mx,Gx);
            double Z = Math.Sqrt(ML.Count)*D;
            double rez = 0;
            double A1 = (double)1 / (18 * ML.Count);
            for (int k = 1; k < 75; k++)
            {
                double f1 = k * k - 0.5 * (1 - Math.Pow(-1, k));
                double f2 = 5 * k * k + 22 - 7.5 * (1 - Math.Pow(-1, k));
                //double O = Math.Pow(Z, 13) / Math.Pow(ML.Count,2);
                double C1 = (f1 - 4 * (f1 + 3)) * Math.Pow(Z * k, 2) + 8 * Math.Pow(k * Z, 4);
                double C2 = (f2 * f2 / 5 - 4 * (f2 + 45) * Math.Pow(Z * k, 2) / (15) + 8 * Math.Pow(k * Z, 4));
                double B1 = 1-2*k*k*Z/(3*Math.Sqrt(ML.Count));
                double A2 = k * k * Z / (27 * Math.Pow(ML.Count, 1.5));
                double G = Math.Pow(-1, k) * Math.Exp((double)-2 * Math.Pow(k * Z, 2));
                rez += G * (B1 - A1 * C1 + A2 * C2); 
            }
            rez*=2;
            return -rez;
        }
        static private double DFound( List<double> ML,InitialStatisticalAnalys gr,int Type,double Mx,double Gx)
        {
            double[] D = new double[2];
            for (int i = 0; i < gr.l.Count - 1; i++)
            {
                double l=0;
                double lmin=0;
                if (Type==0)
                {
                    double v1 = Distributions.NormalFFound((gr.l[i] - Mx) / Gx);
                    l = Math.Abs(gr.F[i]-v1);
                    if (i>0)
                    {
                        v1 = Distributions.NormalFFound((gr.l[i - 1] - Mx) / Gx);
                        lmin = Math.Abs(gr.F[i] - v1);    
                    }
                }
                else if (Type == 1)
                {
                    l = Math.Abs(gr.F[i] - 1 + Math.Pow(2.73, -(1 / (Mx - gr.Min.Q)) * (gr.l[i] - gr.Min.Q)));
                    if (i > 0)
                        lmin = Math.Abs(gr.F[i] - 1 + Math.Pow(2.73, -(1 / (Mx - gr.Min.Q)) * (gr.l[i - 1] - gr.l[0])));
                }
                else if (Type == 2)
                {
                    l = Math.Abs(gr.F[i] - ((gr.l[i] - gr.l[0]) / gr.Len.Q));
                    if (i > 0)
                        lmin = Math.Abs(gr.F[i] - ((gr.l[i - 1] - ML[0]) / gr.Len.Q));
                }
                D[1] = Math.Max(l, D[1]);
                D[0] = Math.Max(lmin, D[0]);
            }
            return Math.Max(D[0],D[1]);
        }
    }
}
