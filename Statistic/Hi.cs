using System;
using System.Collections.Generic;

namespace testgistogr
{
    static class Hi
    {
        static public double HiSqurdFound(InitialStatisticalAnalys gr,int l)
        {
            double Hi = 0;
            List<double> Y4 = new List<double>();
            double L = (double)1 / gr.Mx.Q;
            for (double i = gr.l[0], v = 0; v < gr.Y2.Count; i += gr.Step.Q, v++)
            {
                if (l == 0)
                {
                    Y4.Add(Distributions.NormalFFound(((i + gr.Step.Q) - gr.Mx.Q) / gr.Gx.Q) - Distributions.NormalFFound((i - gr.Mx.Q) / gr.Gx.Q));
                }
                else if (l == 1)
                {
                    double F2 = 1 - Math.Exp(-(1 / gr.Mx.Q) * (gr.Min.Q + gr.Step.Q * (v + 1) ));
                    double F1 = 1 - Math.Exp(-(1 / gr.Mx.Q) * (gr.Min.Q + gr.Step.Q * (v)));
                    Y4.Add(F2 - F1);
                }
                else if (l == 2)
                {
                    Y4.Add(gr.Step.Q / gr.Len.Q);
                }
            }
            for (int i = 0; i < gr.m.Q; i++)
            {
                Hi += (Math.Pow(Y4[i] - gr.f[i], 2) / Y4[i] ) * gr.l.Count;
            }
            return Math.Round(Hi, 4);
        }
        static public double HIF(double alf,int m)
        {
            double d=0 ;
            if (0.5<=alf&&alf<=0.999)
                d = 2.0637*Math.Pow(Math.Log(1/(1-alf))-0.16,0.4274)-1.5774;
            if (0.001<=alf&&alf<0.5)
                d = -2.0637*Math.Pow(Math.Log(1/(alf))-0.16,0.4274)+1.5774;
            double A = d*Math.Sqrt(2);
            double B = 2 * (Math.Pow(d, 2) - 1) / 3;
            double C = d * (Math.Pow(d, 2) - 7) / (Math.Sqrt(2) * 9);
            double D = (6 * Math.Pow(d, 4)+14*Math.Pow(d,2)-32) / (405);
            double E = d*(9*Math.Pow(d,4)+256*Math.Pow(d,2)-433)/(Math.Sqrt(2)*4860);
            return m + A * Math.Sqrt(m) + B + C / Math.Sqrt(m)+D/m+E/(Math.Sqrt(m)*m);
        }
    }
}
