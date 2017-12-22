using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using TMO;

namespace IS.Data
{
    public class List : System.Collections.Generic.List<PointF>
    {
        public double QMax;
        
        public new void Add(PointF item)
        {
            base.Add(item);
            if (base.Count == 1)
                QMax = base[0].Y;
            else
                QMax = Math.Max(QMax, item.Y);
        }
    }

    public class PracticalData
    {
        public string Name;
        public List Q;
        public PracticalData(string Name)
        {
            this.Name = Name;
            Q = new List();
            
        }
    }
    class Evaluation
    {
        public List<TMO.Data> Teorical;
        public List<PracticalData> Practical;
        public List<Series> EvauationSeries;

        TMO.Data PFailure;
        TMO.Data Q;
        TMO.Data A;
        TMO.Data AverZ;
        TMO.Data AverR;
        TMO.Data AverK;
        TMO.Data AverTWait;

        PracticalData PFailurePractical;
        PracticalData QPractical;
        PracticalData APractical;
        PracticalData AverZPractical;
        PracticalData AverRPractical;
        PracticalData AverKPractical;
        PracticalData AverTWaitPractical;


        public void SetData(ImitationSystem IS)
        {
            SetDataHeader();
            SetTeoricalData(IS);
            SetPracticalData(IS);
            SetSeries(IS.lambda);
        }


        private void SetDataHeader()
        {
            PFailure = new TMO.Data()
            {
                Name = "Ймовірність відмови(Pf)"
            };
            Q = new TMO.Data()
            {
                Name = "Відносна пропускна здатність(Q)"
            };
            A = new TMO.Data()
            {
                Name = "Абсолютна пропускна здатність(A)"
            };
            AverZ = new TMO.Data()
            {
                Name = "Середнє число зайнятих каналів(AverZ)"
            };
            AverR = new TMO.Data()
            {
                Name = "Середня кількість вимог в черзі(AverR)"
            };
            AverK = new TMO.Data()
            {
                Name = "Середня кількість вимог в системі(AverK)"
            };
            AverTWait = new TMO.Data()
            {
                Name = "Середній час очікування в черзі(AverTWait)"
            };
            Teorical = new List<TMO.Data>()
            {
                PFailure, Q, A, AverZ, AverR, AverK//, AverTWait
            };

            PFailurePractical = new PracticalData(PFailure.Name);
            QPractical = new PracticalData(Q.Name);
            APractical = new PracticalData(A.Name);
            AverZPractical = new PracticalData(AverZ.Name);
            AverRPractical = new PracticalData(AverR.Name);
            AverKPractical = new PracticalData(AverK.Name);
            AverTWaitPractical = new PracticalData(AverTWait.Name);

            Practical = new List<PracticalData>()
            {
                PFailurePractical, QPractical, APractical, AverZPractical,
                AverRPractical, AverKPractical//, AverTWaitPractical
            };
        }

        private void SetTeoricalData(ImitationSystem IS)
        {
            double ro = IS.lambda / IS.mu;

            if (!IS.UnlimitedQueue)
            {
                PFailure.Q = IS.pTeorical[IS.pTeorical.Length - 1];
                A.Q = IS.lambda * (1 - IS.pTeorical[IS.pTeorical.Length - 1]);
                AverZ.Q = A.Q / IS.mu;
                if (IS.m != 0)
                    AverR.Q = IS.pTeorical[IS.n + 1] * (1 - (IS.m + 1) * Math.Pow(ro / IS.n, IS.m) + IS.m * Math.Pow(ro / IS.n, IS.m + 1)) /
                        Math.Pow(1 - (ro / IS.n), 2);
                AverK.Q = AverR.Q + AverZ.Q;
                AverTWait.Q = AverR.Q / IS.lambda;
            }
            else
            {
                PFailure.Q = 0;
                PFailure.Q0 = IS.NFailure / (IS.NFailure + IS.NProcessing);
                A.Q = IS.lambda;
                AverZ.Q = A.Q / IS.mu;
                if (IS.pTeorical.Length > IS.n + 1)
                    AverR.Q = IS.pTeorical[IS.n + 1] /
                        (Math.Pow(1 - ro / IS.n, 2));
                AverK.Q = AverR.Q + AverZ.Q;
                AverTWait.Q = AverR.Q / IS.lambda;
            }
            Q.Q = 1 - PFailure.Q;
        }

        private void SetPracticalData(ImitationSystem IS)
        {

            double ChennelCount = 1;
            double QueueCount = 0;

            double SumChannel = 0;
            double SumQueue = 0;
            double SumTimeWait = 0;

            double QueueExitCount = 0;
            double ChennelCountChange = 1;
            double NumberProccesing = 0;
            double QueueCountChange = 0;
            for (int i = 1; i < IS.TimeLine.Count; i++)
            {
                SumChannel += (IS.TimeLine[i].TimeEvent - IS.TimeLine[i - 1].TimeEvent) * ChennelCount;
                SumQueue += (IS.TimeLine[i].TimeEvent - IS.TimeLine[i - 1].TimeEvent) * QueueCount;
                ChennelCountChange = IS.TimeLine[i].TimeEvent;
                QueueCountChange = IS.TimeLine[i].TimeEvent;
                if (IS.TimeLine[i].InputOutput == EventType.Input)
                    ChennelCount++;
                else if (IS.TimeLine[i].InputOutput == EventType.Output)
                {
                    ChennelCount--;
                    NumberProccesing++;
                }
                else if (IS.TimeLine[i].InputOutput == EventType.Wait && IS.TimeLine[i].TimeProcces != 0)
                    QueueCount++;
                else if (IS.TimeLine[i].InputOutput == EventType.Wait && IS.TimeLine[i].TimeProcces == 0)
                {
                    SumTimeWait += IS.TimeLine[i].Neighbor.TimeProcces * (IS.TimeLine[i].TimeEvent - IS.TimeLine[i - 1].TimeEvent) * 1.001;
                    QueueExitCount += (IS.TimeLine[i].TimeEvent - IS.TimeLine[i - 1].TimeEvent);
                    AverTWaitPractical.Q.Add(new PointF((float)IS.TimeLine[i].TimeEvent, (float)(SumTimeWait / QueueExitCount)));
                    QueueCount--;
                }
                APractical.Q.Add(new PointF((float)IS.TimeLine[i].TimeEvent, 
                    (float)NumberProccesing / (float)IS.TimeLine[i].TimeEvent));
                AverZPractical.Q.Add(new PointF((float)IS.TimeLine[i].TimeEvent, 
                    (float)(SumChannel / ChennelCountChange)));
                AverRPractical.Q.Add(new PointF((float)IS.TimeLine[i].TimeEvent, 
                    (float)(SumQueue / QueueCountChange)));
                AverKPractical.Q.Add(new PointF((float)IS.TimeLine[i].TimeEvent, 
                    (float)(SumQueue / QueueCountChange) + (float)(SumChannel / ChennelCountChange)));

            }
            SetPFailureAndQPractical(IS);
            SetQ0((float)IS.TimeLine[IS.TimeLine.Count - 1].TimeEvent);
        }

        private void SetQ0(float EndTime)
        {
            for (int i = 0; i < Practical.Count; i++)
                if (Practical[i].Q.Count != 0)
                {
                    Practical[i].Q.Add(new PointF(EndTime, Practical[i].Q[Practical[i].Q.Count - 1].Y));
                    Teorical[i].Q0 = Practical[i].Q[Practical[i].Q.Count - 1].Y;
                }
        }

        private void SetPFailureAndQPractical(ImitationSystem IS)
        {
            double NF = 0, NP = 0;
            for (int f = 1, p = 1; f < IS.NFailureElem.Count && p < IS.NProcessingElem.Count;)
            {
                NF = IS.NFailureElem[f].TimeProcces;
                NP = IS.NProcessingElem[p].TimeProcces;
                PFailurePractical.Q.Add(
                    new PointF((float)Math.Max(IS.NFailureElem[f].TimeEvent, IS.NProcessingElem[p].TimeEvent),
                    (float)(NF / (NP + NF))));
                QPractical.Q.Add(new PointF(
                    PFailurePractical.Q[PFailurePractical.Q.Count - 1].X,
                    1 - PFailurePractical.Q[PFailurePractical.Q.Count - 1].Y));
                if (p < IS.NProcessingElem.Count && IS.NFailureElem[f].TimeEvent > IS.NProcessingElem[p].TimeEvent)
                    p++;
                else
                    f++;
            }
        }

        private void SetSeries(double lambda)
        {
            Color[] colors = new Color[]
            {
                Color.Red,
                Color.Green,
                Color.Orange,
                Color.Blue,
                Color.Purple,
                Color.Teal,
                Color.Black,
                Color.Brown
            };
            EvauationSeries = new List<Series>();
            for (int i = 0; i < Practical.Count; i++)
            {
                Series ser = new Series(Practical[i].Name)
                {
                    ChartType = SeriesChartType.Line,
                    ["PointWidth"] = "2",
                    BorderWidth = 2,
                    BorderColor = Color.Black,
                    Color = colors[i],
                    IsVisibleInLegend = false
                };
                EvauationSeries.Add(ser);
                for (int j = 0; j < Practical[i].Q.Count; j ++)
                {
                    if (ser.Points.Count == 0 ||
                        Math.Abs(Practical[i].Q[j].Y - ser.Points[ser.Points.Count - 1].YValues[0] ) 
                        > Practical[i].Q[Practical[i].Q.Count - 1].Y * 1.5 / 100)
                    ser.Points.AddXY(Practical[i].Q[j].X, Practical[i].Q[j].Y);
                    else if (j == Practical[i].Q.Count - 1)
                        ser.Points.AddXY(Practical[i].Q[j].X, Practical[i].Q[j].Y);
                }
            }
            for (int i = 0; i < Practical.Count; i++)
            {
                Series ser = new Series(Practical[i].Name + "Teorical")
                {
                    ChartType = SeriesChartType.Line,
                    ["PointWidth"] = "2",
                    BorderWidth = 2,
                    BorderColor = Color.Black,
                    Color = colors[i],
                    IsVisibleInLegend = false
                };
                EvauationSeries.Add(ser);
                ser.Points.AddXY(0, Teorical[i].Q);
                ser.Points.AddXY(Practical[i].Q[Practical[i].Q.Count - 1].X, Teorical[i].Q);
            }
        }
    }
}
