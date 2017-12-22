using System;
using System.Collections.Generic;

namespace testgistogr
{
    public static class RegresTypeName
    {
        public static string LineRegresion = "a + b*x";
        public static string ParabRegresion = "a + b*x +c*x^2";
        public static string SqrtRegresion = "sqrt(a + b*x)";
        public static string _1_X_Regresion = "a + b/x";
        public static string _1_Y_Regresion = "1/(a + b*x)";
        public static string X_Y_Regresion = "x/(a + b*x)";
        public static string LnxRegresion = "a + b*Ln(x)";
        public static string StepenRegresion = "a*x^b";
        public static string ERegresion = "a*e^(bx)";
        public static string KubRegresion = "a + b*x^3";
        public static string E_1_X_Regresion = "a*e^(b/x)";
        public static string B_Regresion = "b / (a + x)";
        public static string Eo1_YRegresion = "1/(a+b*e^(-x))";
        public static string TypeRegresion = LineRegresion + "\n" + ParabRegresion + "\n" + SqrtRegresion
            + "\n" + _1_X_Regresion + "\n" + _1_Y_Regresion + "\n" + X_Y_Regresion + "\n" + LnxRegresion + "\n"
            + StepenRegresion + "\n" + ERegresion + "\n" + KubRegresion + "\n" + E_1_X_Regresion 
            /*+ B_Regresion + "\n" + Eo1_YRegresion*/;
    }
    static class RegresType
    {
        public static double Model(double x, List<Data> Q, string Name) 
        {
            if (Name == RegresTypeName.LineRegresion)
                return Q[0].Q + Q[1].Q * x;
            else if (Name == RegresTypeName.ParabRegresion)
                return Q[0].Q + Q[1].Q * x + Q[2].Q * x * x;
            else if (Name == RegresTypeName.SqrtRegresion)
                return Math.Sqrt(Q[0].Q + Q[1].Q * x);
            else if (Name == RegresTypeName._1_X_Regresion)
                return Q[0].Q + Q[1].Q / x;
            else if (Name == RegresTypeName._1_Y_Regresion)
                return 1.0 / (Q[0].Q + Q[1].Q * x);
            else if (Name == RegresTypeName.X_Y_Regresion)
                return x / (Q[0].Q + Q[1].Q * x);
            else if (Name == RegresTypeName.LnxRegresion)
                return Q[0].Q + Q[1].Q * Math.Log(x);
            else if (Name == RegresTypeName.StepenRegresion)
                return Q[0].Q * Math.Pow(x, Q[1].Q);
            else if (Name == RegresTypeName.ERegresion)
                return Q[0].Q * Math.Pow(Math.E, Q[1].Q * x);
            else if (Name == RegresTypeName.KubRegresion)
                return Q[0].Q + Q[1].Q * Math.Pow(x, 3);
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return Q[0].Q * Math.Pow(Math.E, Q[1].Q / x);
            else if (Name == RegresTypeName.B_Regresion)
                return Q[1].Q / (x + Q[0].Q);
            else if (Name == RegresTypeName.Eo1_YRegresion)
                return 1.0 / (Q[0].Q + Q[1].Q * Math.Pow(Math.E, -x));
            return double.NaN;
        }
        public static double Model(double x, double[] Q, string Name)
        {
            if (Name == RegresTypeName.LineRegresion)
                return Q[0] + Q[1] * x;
            else if (Name == RegresTypeName.ParabRegresion)
                return Q[0] + Q[1] * x + Q[2] * x * x;
            else if (Name == RegresTypeName.SqrtRegresion)
                return Math.Sqrt(Q[0] + Q[1] * x);
            else if (Name == RegresTypeName._1_X_Regresion)
                return Q[0] + Q[1] / x;
            else if (Name == RegresTypeName._1_Y_Regresion)
                return 1.0 / (Q[0] + Q[1] * x);
            else if (Name == RegresTypeName.X_Y_Regresion)
                return x / (Q[0] + Q[1] * x);
            else if (Name == RegresTypeName.LnxRegresion)
                return Q[0] + Q[1] * Math.Log(x);
            else if (Name == RegresTypeName.StepenRegresion)
                return Q[0] * Math.Pow(x, Q[1]);
            else if (Name == RegresTypeName.ERegresion)
                return Q[0] * Math.Pow(Math.E, Q[1] * x);
            else if (Name == RegresTypeName.KubRegresion)
                return Q[0] + Q[1] * Math.Pow(x, 3);
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return Q[0] * Math.Pow(Math.E, Q[1] / x);
            else if (Name == RegresTypeName.B_Regresion)
                return Q[1] / (x + Q[0]);
            else if (Name == RegresTypeName.Eo1_YRegresion)
            {
                var er = Math.Pow(Math.E, -x);
                var te = 1.0 / (Q[0] + Q[1] * er);
                return 1.0 / (Q[0] + Q[1] * Math.Pow(Math.E, -x));
            }
            return double.NaN;
        }
        public static double FiX(double x, string Name) 
        {
            if (Name == RegresTypeName.LineRegresion)
                return x;
            else if (Name == RegresTypeName.SqrtRegresion)
                return x;
            else if (Name == RegresTypeName._1_X_Regresion)
                return 1.0 / x;
            else if (Name == RegresTypeName._1_Y_Regresion)
                return x;
            else if (Name == RegresTypeName.X_Y_Regresion)
                return x ;
            else if (Name == RegresTypeName.LnxRegresion)
                return Math.Log(x);
            else if (Name == RegresTypeName.StepenRegresion)
                return Math.Log(x);
            else if (Name == RegresTypeName.ERegresion)
                return x;
            else if (Name == RegresTypeName.KubRegresion)
                return Math.Pow(x, 3);
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return 1.0/x;
            else if (Name == RegresTypeName.B_Regresion)
                return 1.0/x;
            else if (Name == RegresTypeName.Eo1_YRegresion)
                return Math.Pow(Math.E,-x);
            return double.NaN;
        }

        public static double FiY(double y, double x, string Name)
        {
            if (Name == RegresTypeName.LineRegresion)
                return y;
            else if (Name == RegresTypeName.SqrtRegresion)
                return Math.Pow(y, 2);
            else if (Name == RegresTypeName._1_X_Regresion)
                return y;
            else if (Name == RegresTypeName._1_Y_Regresion)
                return 1.0 / y;
            else if (Name == RegresTypeName.X_Y_Regresion)
                return x / y;
            else if (Name == RegresTypeName.LnxRegresion)
                return y;
            else if (Name == RegresTypeName.StepenRegresion)
            {
                double rez = Math.Log(y);
                if (double.IsNaN(rez) || double.IsInfinity(rez))
                    return 1;
                return rez;
            }
            else if (Name == RegresTypeName.ERegresion)
                return Math.Log(y);
            else if (Name == RegresTypeName.KubRegresion)
                return y;
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return Math.Log(y);
            else if (Name == RegresTypeName.B_Regresion)
                return 1.0 / y;
            else if (Name == RegresTypeName.Eo1_YRegresion)
                return 1.0 / y;
            return double.NaN;
        }
        public static double W(double y, double x, string Name)
        {
            if (Name == RegresTypeName.LineRegresion)
                return 1;
            else if (Name == RegresTypeName.SqrtRegresion)
                return 1.0/(4*Math.Pow(y, 2));
            else if (Name == RegresTypeName._1_X_Regresion)
                return 1.0/Math.Pow(x,4);
            else if (Name == RegresTypeName._1_Y_Regresion)
                return Math.Pow(y,4);
            else if (Name == RegresTypeName.X_Y_Regresion)
                return Math.Pow(y,4)/Math.Pow(x,2);
            else if (Name == RegresTypeName.LnxRegresion)
                return 1.0/Math.Pow(x,2);
            else if (Name == RegresTypeName.StepenRegresion)
                return Math.Pow(y, 2) / Math.Pow(x, 2);
            else if (Name == RegresTypeName.ERegresion)
                return Math.Pow(y,2);
            else if (Name == RegresTypeName.KubRegresion)
                return 9*Math.Pow(x,4);
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return Math.Pow(y,2)/Math.Pow(x,4);
            else if (Name == RegresTypeName.B_Regresion)
                return Math.Pow(y, 4) / Math.Pow(x, 4);
            else if (Name == RegresTypeName.Eo1_YRegresion)
                return Math.Pow(y, 4)*Math.Pow(Math.E,-2*x);
            return double.NaN;
        }
        public static double A(double a, string Name) 
        {
            if (Name == RegresTypeName.LineRegresion)
                return a;
            else if (Name == RegresTypeName.SqrtRegresion)
                return a;
            else if (Name == RegresTypeName._1_X_Regresion)
                return a;
            else if (Name == RegresTypeName._1_Y_Regresion)
                return a;
            else if (Name == RegresTypeName.X_Y_Regresion)
                return a ;
            else if (Name == RegresTypeName.LnxRegresion)
                return a;
            else if (Name == RegresTypeName.StepenRegresion)
                return Math.Exp(a);
            else if (Name == RegresTypeName.ERegresion)
                return Math.Exp(a);
            else if (Name == RegresTypeName.KubRegresion)
                return a;
            else if (Name == RegresTypeName.E_1_X_Regresion)
                return Math.Exp(a);
            else if (Name == RegresTypeName.B_Regresion)
                return a;
            else if (Name == RegresTypeName.Eo1_YRegresion)
                return a;
            return double.NaN;
        }

    }
}
