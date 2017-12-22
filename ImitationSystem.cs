using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows;
using System.Windows.Forms.DataVisualization.Charting;
using TMO;

namespace IS
{
    class ImitationSystem
    {
        public List<Event> TimeLine = new List<Event>();
        public double lambda = 1;
        public double mu = 1;
        public int n;
        public int m;
        public int MaxQueueCount;
        public List<double> LambdaElem;
        public List<double> MuElem;
        public List<Event> NFailureElem;
        public List<Event> NProcessingElem;
        public double NFailure;
        public double NProcessing;
        public Series[] pPracticalSeries;
        public double[] pPractical;
        public double[] pTeorical;
        public int round = 4;
        public Data.Evaluation Evaluation;
        public bool UnlimitedQueue;

        public ImitationSystem(double lambda, double mu, int n, int m)
        {
            this.lambda = lambda;
            this.mu = mu;
            this.n = n;
            this.m = m;
            UnlimitedQueue = false;
        }

        public ImitationSystem(double lambda, double mu, int n)
        {
            this.lambda = lambda;
            this.mu = mu;
            this.n = n;
            this.m = 0;
            UnlimitedQueue = true;
        }
        
        public void Generate(double EndTime)
        {
            Channels channels = new Channels(n);
            Queue<Event> queue = new Queue<Event>();
            
            Colculate(EndTime, ref channels, ref queue);

            NFailureElem.Add(new Event(NProcessingElem[NProcessingElem.Count - 1].TimeEvent, NFailure, EventType.Input));

            TimeLineMerge(channels, queue);
            PTeoricalFound();
            PPracticalFound();

            Evaluation = new Data.Evaluation();
            Evaluation.SetData(this);
        }

        private void Colculate(double EndTime, ref Channels channels, ref Queue<Event> queue)
        {

            TimeLine = new List<Event>();

            NFailure = 0;
            NProcessing = 0;

            NFailureElem = new List<Event>()
            { new Event(0, 0, EventType.Input) };
            NProcessingElem = new List<Event>()
            { new Event(0, 0, EventType.Input) };

            Random rand = new Random();

            LambdaElem = new List<double>();
            MuElem = new List<double>();
            LambdaElem.Add(Math.Log(1.0 / (1 - rand.NextDouble())) / lambda);
            MuElem.Add(Math.Log(1.0 / (1 - rand.NextDouble())) / mu);

            for (double TimeCome = 0, TimeLeave = MuElem[MuElem.Count - 1];
                TimeCome <= EndTime; 
                TimeCome += LambdaElem[LambdaElem.Count - 1],
                LambdaElem.Add(Math.Log(1.0 / (1 - rand.NextDouble())) / lambda),
                MuElem.Add(Math.Log(1.0 / (1 - rand.NextDouble())) / mu))
            {
                TimeLeave = TimeCome + MuElem[MuElem.Count - 1];

                Event EventCome = new Event(TimeCome, TimeLeave - TimeCome, EventType.Input, LambdaElem.Count - 1),
                EventLeav = new Event(TimeLeave, 0, EventType.Output, LambdaElem.Count - 1);
                EventCome.Neighbor = EventLeav;
                EventLeav.Neighbor = EventCome;

                Channel firstFree = channels.FirstFree();
                queue.QueueClear(channels, firstFree, queue, TimeCome, ref NProcessing, ref NProcessingElem);
                if (queue.Count == 0 && firstFree.FreeTime < EventCome.TimeEvent)
                {
                    firstFree.ChangeEvent(EventCome);
                    NProcessing++;
                    NProcessingElem.Add(new Event(EventCome.TimeEvent, NProcessing, EventType.NFailure));
                }
                else
                {
                    if (queue.Count < m || UnlimitedQueue)
                    {
                        queue.Enqueue(EventCome);
                        MaxQueueCount = Math.Max(MaxQueueCount, queue.Count);
                        queue.QueueClear(channels, firstFree, queue, TimeCome, ref NProcessing, ref NProcessingElem);
                    }
                    else
                    {
                        NFailure++;
                        NFailureElem.Add(new Event(TimeCome, NFailure, EventType.NFailure));
                    }
                }

            }
        }

        private void PTeoricalFound()
        {
            pTeorical = new double[n + m + 1];
            double ro = lambda / mu;
            if (UnlimitedQueue)
            {
                if ((ro / n) >= 1)
                    return;
                pTeorical = new double[n + MaxQueueCount + 1];
                int factor = 1;
                for (int j = 0; j <= n; j++, factor = Factorial(j))
                    pTeorical[0] += Math.Pow(ro, j) / factor;
                factor /= n + 1;
                pTeorical[0] += Math.Pow(ro, n + 1) / ((n - ro) * factor);
                pTeorical[0] = 1.0 / pTeorical[0];
                int i = 1;
                for (; i <= n; i++)
                    pTeorical[i] = Math.Pow(ro, i) / (Factorial(i)) * pTeorical[0];
                for (; i < pTeorical.Length; i++)
                    pTeorical[i] = Math.Pow(ro, i) * pTeorical[0] / (Math.Pow(n, i - n) * factor);
            }
            else
            {
                if (n == 1)
                {
                    for (int i = 0; i < pTeorical.Length; i++)
                        pTeorical[0] += Math.Pow(ro, i);
                    pTeorical[0] = 1.0 / pTeorical[0];
                    for (int i = 1; i < pTeorical.Length; i++)
                        pTeorical[i] = pTeorical[0] * Math.Pow(ro, i);
                }
                else
                {
                    for (int i = 0; i < pTeorical.Length; i++)
                        pTeorical[0] += Math.Pow(ro, i) /
                            (Factorial(Math.Min(i, n)) * Math.Pow(n, Math.Max(i - n, 0)));
                    pTeorical[0] = 1.0 / pTeorical[0];
                    for (int i = 1; i <= n; i++)
                        pTeorical[i] = pTeorical[0] * Math.Pow(ro, i) / Factorial(i);
                    for (int i = n + 1, fact = Factorial(n); i < pTeorical.Length; i++)
                        pTeorical[i] = pTeorical[0] * Math.Pow(ro, i) /
                            (fact * Math.Pow(n, i - n));

                }
            }
        }

        private int Factorial(int num)
        {
            if (num <= 1)
                return 1;
            int result = num;
            return Factorial(num - 1) * num;
        }

        private void PPracticalFound()
        {
            pPractical = new double[n + MaxQueueCount + 1];
            pPracticalSeries = new Series[2 * (n + MaxQueueCount + 1)];
            for (int i = 0; i < pPractical.Length; i++)
            {
                pPracticalSeries[i] = new Series()
                {
                    Name = "p" + i.ToString() + "  " + Math.Round(pTeorical[i], round).ToString(),
                    ChartType = SeriesChartType.Line,
                    ["PointWidth"] = "2",
                    BorderWidth = 2,
                    BorderColor = Color.Black,
                    Color = Color.FromArgb(
                        (int)(i * 255.0 / pPractical.Length),
                        0,
                        (int)(255 - i * 255.0 / pPractical.Length)),
                    IsVisibleInLegend = false
                };
                pPracticalSeries[i + n + MaxQueueCount + 1] = new Series()
                {
                    Name = "pTeorical" + i.ToString(),
                    ChartType = SeriesChartType.Line,
                    ["PointWidth"] = "2",
                    BorderWidth = 2,
                    BorderColor = Color.Black,
                    Color = Color.FromArgb(
                        (int)(i * 255.0 / pPractical.Length),
                        0,
                        (int)(255 - i * 255.0 / pPractical.Length)),
                    IsVisibleInLegend = false
                };
                pPracticalSeries[i + n + MaxQueueCount + 1].Points.AddXY(0, pTeorical[i]);
                pPracticalSeries[i + n + MaxQueueCount + 1].Points.AddXY(TimeLine[TimeLine.Count - 1].TimeEvent, pTeorical[i]);
            }
            int z = 1;
            for (int i = 1; i < TimeLine.Count; i++)
            {
                pPractical[z] += TimeLine[i].TimeEvent - TimeLine[i - 1].TimeEvent;
                if (pPracticalSeries[z].Points.Count == 0 ||
                    Math.Abs(pPracticalSeries[z].Points[pPracticalSeries[z].Points.Count - 1].YValues[0]
                    - pPractical[z] / TimeLine[i].TimeEvent) 
                    > 0.0008)
                    pPracticalSeries[z].Points.AddXY(TimeLine[i].TimeEvent, pPractical[z] / TimeLine[i].TimeEvent);
                if (TimeLine[i].InputOutput == EventType.Input ||
                    (TimeLine[i].InputOutput == EventType.Wait && TimeLine[i].TimeProcces != 0))
                    z++;
                else if (TimeLine[i].InputOutput == EventType.Output ||
                    (TimeLine[i].InputOutput == EventType.Wait && TimeLine[i].TimeProcces == 0))
                    z--;
            }
            for (int i = 0; i < pPractical.Length; i++)
                pPractical[i] /= TimeLine[TimeLine.Count - 1].TimeEvent;
            for (int i = 0; i < pPractical.Length; i++)
                pPracticalSeries[i].Points.AddXY(TimeLine[TimeLine.Count - 1].TimeEvent, pPractical[i]);
        }

        private void TimeLineMerge(Channels channels, Queue<Event> queue)
        {
            TimeLine = new List<Event>(channels[0].TimeLine);
            for (int i = 1; i < channels.Length; i++)
                TimeLine.AddRange(channels[i].TimeLine);
            TimeLine.AddRange(queue.TimeLine);
            TimeLine.Sort(Compare);
        }

        private int Compare(Event e1, Event e2)
        {
            if (e1.TimeEvent < e2.TimeEvent)
                return -1;
            if ((e1.TimeEvent - e2.TimeEvent) == 0)
            {
                if (e1.InputOutput == EventType.Input &&
                    e2.InputOutput == EventType.Output)
                    return 1;
                else if (e2.InputOutput == EventType.Input &&
                    e1.InputOutput == EventType.Output)
                    return -1;/*
                if (e1.InputOutput == EventType.Wait && e1.TimeProcces != 0 &&
                    e2.InputOutput == EventType.Output)
                    return 1;
                if (e1.InputOutput == EventType.Input && e2.TimeProcces == 0 &&
                    e2.InputOutput == EventType.Wait)
                    return 1;*/
                return 0;
            }
            else 
                return 1;
        }

    }
}
