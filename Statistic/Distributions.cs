using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.IO;
using TestSimpleRNG;

namespace testgistogr
{
    static class Distributions
    {
        public static string Normal = "Нормальний";
        public static string Exp = "Експоненціальний";
        public static string Line = "Рівномірний";
        static public double NormalfFound(double x , double m , double g) 
        {
            return Math.Exp(-Math.Pow((x-m)/g,2)/2)/(g*Math.Sqrt(2*Math.PI));
        }
        static public double NormalDoublefFound(double x, double mx, double gx, double y, double my, double gy, double r)
        {
            double A = 1.0/(2*Math.PI*gx*gy*Math.Sqrt(1 - Math.Pow(r,2)));
            double B = -1.0/(2.0*(1 - Math.Pow(r,2)));
            double C = Math.Pow((x - mx)/gx,2) - 2*r*((x - mx)/gx)*((y - my)/gy) + Math.Pow((x - my)/gy,2); 
            return A*Math.Exp(B*C);
        }

        static public double ExpFound(double x,double Mx)
        {
            double l = 1.0 / Mx;
            return 1 - Math.Exp(-l * x);
        }
        static public double[,] NormalFRead()
        {
            double[,] rez = new double[291, 2];
            FileStream fs = new FileStream("F.txt", FileMode.Open);
            StreamReader streamReader = new StreamReader(fs, Encoding.ASCII);
            String str = streamReader.ReadToEnd();
            str = str.Replace("\t", " ").Replace("\r", " ").Replace("\n", " ").Replace("   ", " ").Replace("  ", " ").Replace(",", ".");
            int u = 0;
            double[] d = new double[5];

            for (int i = 0; i < str.Length; i++, u++)
            {
                String hstr = "";
                for (; i < str.Length && str[i] != ' ' && str[i] != '\t'; i++)
                {
                    hstr += str[i];
                }
                rez[u, 0] = Convert.ToDouble(hstr);
                hstr = "";
                if (i != str.Length - 2)
                {
                    i++;
                }
                for (; i < str.Length && str[i] != ' ' && str[i] != '\t'; i++)
                {
                    hstr += str[i];
                }
                rez[u, 1] = Convert.ToDouble(hstr);
            }
            streamReader.Close();
            fs.Close();
            return rez;
        }
        static public double NormalFFound(double x, double[,] F)
        {
            double rez = -10;
            double X = Math.Abs(Math.Round(x, 2));
            for (int i = 0; i < 291; i++)
            {
                if (X == F[i, 0])
                    rez = F[i, 1];
            }
            if (rez == -10)
                rez = 0.5;
            if (X * x < 0)
            {
                rez *= -1;
            }
            return rez;
        }
        static public double NormalFobrFound(double x, double[,] F)
        {
            double rez = -10;
            double X = Math.Abs(x);
            for (int i = 0; i < 291; i++)
            {
                if (Math.Round(X, 2) == Math.Round(F[i, 1], 2))
                    rez = F[i, 0];
            }
            if (X * x < 0)
            {
                rez *= -1;
            }
            return rez;
        }
        static public double NormalFFound(double x)
        {
            if (x < 0)
            {
                return 1 - NormalFFound(Math.Abs(x));
            }
            double[] b = { 0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429 };
            double e = 7.8 * Math.Pow(10, -8);
            double p = 0.2316419;
            double T = 1 / (1 + p * x);
            double hlpdoubl = Math.Pow(Math.E, -(x * x) / 2) / Math.Sqrt(2 * Math.PI);
            double rez = 1 - hlpdoubl * (b[0] * T + b[1] * Math.Pow(T, 2) + b[2] * Math.Pow(T, 3) + b[3] * Math.Pow(T, 4) + b[4] * Math.Pow(T, 5)) + e;
            return rez;
        }
        static public double NormalFobrFound(double x)
        {
            double X= x;
            double a1 = -5;
            double b1 = 5;
            double c1 = 0;
            double a2 = NormalFFound(a1) - 0.52;
            double b2 = NormalFFound(b1) - 0.48;
            double c2 = NormalFFound(c1) - 0.5;
            for (;Math.Abs(b1-a1)>0.001 ; )
            {
                if (c2<=x&&x<=b2)
                {
                    a2 = c2;
                    a1 = c1;
                }
                if (a2 <= x && x <= c2)
                {
                    b2 = c2;
                    b1 = c1;
                }
                c1 = (a1 + b1) / 2;
                c2 = NormalFFound(c1) - 0.5;
            }
            return c1;
        }
        static public double NormalQuantile(double x)
        {
            if (x==0)
            {
                x = 0.01;
            }
            double p=x;
            double alf = x * 2;
            double t = Math.Sqrt(Math.Log(1 / (x * x)));
            double e = 4.5 * Math.Pow(10, -4);
            double[] c = {2.515517,0.802853,0.010328 };
            double[] d = {1.432788,0.1892659,0.001308 };
            double hlpdou1 =(c[0]+c[1]*t+c[2]*t*t);
            double hlpdou2 =(1+d[0]*t+d[1]*t*t+d[2]*t*t*t);
            double rez = t - ((hlpdou1) / (hlpdou2))+e;
            return -rez;

        }
        static public double StudentQuantile(double alf, int m)
        {
            double T = 0;
            double u_alf = Math.Abs(Distributions.NormalQuantile(alf/*1 - 0.05 / 2*/));
            double g1 = (Math.Pow(u_alf, 3) + u_alf) / 4;
            double g2 = (5 * Math.Pow(u_alf, 5) + 16 * Math.Pow(u_alf, 3) + 3 * u_alf) / 96;
            double g3 = (34 * Math.Pow(u_alf, 7) + 19 * Math.Pow(u_alf, 5) + 17 * Math.Pow(u_alf, 3) - 15 * u_alf) / 384;
            double g4 = (79 * Math.Pow(u_alf, 9) + 779 * Math.Pow(u_alf, 7) + 1482 * Math.Pow(u_alf, 5) - 1920 * Math.Pow(u_alf, 3) - 945 * u_alf) / 92160;
            T = u_alf + g1 / (m ) + g2 / Math.Pow(m , 2) + g3 / Math.Pow(m , 3) + g4 / Math.Pow(m , 4);
            return T;
        }
        static public double FisherQuantile(double alf, int m1, int m2)
        {
            double Normal = NormalQuantile(1-alf/2);
            double g = 1.0 / m1 + 1.0 / m2;
            double s = 1.0 / m1 - 1.0 / m2;
            double A = Normal*Math.Sqrt(g / 2);
            double B = s*(Math.Pow(Normal,2) + 2)/6;
            double C = Math.Sqrt(g / 2) * (g * (Math.Pow(Normal, 2) + 3.0 * Normal) / 24.0 + Math.Pow(s, 2) * (Math.Pow(Normal, 3) + 11.0 * Normal)/(72.0*g));
            double D = g * s * (Math.Pow(Normal, 4) + 9.0 * Math.Pow(Normal, 4) + 8.0) / 120.0;
            double E = Math.Pow(s, 3) * (3.0 * Math.Pow(Normal, 4) + 7.0 * Math.Pow(Normal, 2) - 16.0)/(3240.0 * g);
            double F1 = Math.Pow(g, 2) * (Math.Pow(Normal, 5) + 20.0 * Math.Pow(Normal, 3) + 15.0 * Normal)/1920.0;
            double F2 = Math.Pow(s, 4) * (Math.Pow(Normal, 5) + 44.0 * Math.Pow(Normal, 3) + 183.0 * Normal) / 2880.0;
            double F3 = Math.Pow(s, 4) * (9.0 * Math.Pow(Normal, 5)  - 284.0 * Math.Pow(Normal, 3) - 1513.0 * Normal) / (15520.0*g*g);
            double F = Math.Sqrt(g / 2) * (F1 + F2 + F3);
            double Z = A - B + C - D + E + F;
            return Math.Exp(2*Z);
        }
    }
}
