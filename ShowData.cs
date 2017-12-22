using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using TMO;
using testgistogr;
using IS;
using System.Windows.Forms;

namespace TMOLaba4_5
{
    static  class ShowData
    {
        static public void ShowDigrams(ImitationSystem IS, Chart chart)
        {
            chart.Series.Add(DrawN("N відмови", IS.NFailureElem, IS));
            chart.Series.Add(DrawN("N обробки", IS.NProcessingElem, IS));
        }

        static public void ShowLamdaMu(ImitationSystem IS, Chart Lambdachart, Chart Muchart)
        {
            Lambdachart.Series.Add(GistogPaint(new InitialStatisticalAnalys(IS.LambdaElem), Color.DarkBlue, "Lambda"));
            Muchart.Series.Add(GistogPaint(new InitialStatisticalAnalys(IS.MuElem), Color.DarkBlue, "Mu"));
        }
        static public Series GistogPaint(InitialStatisticalAnalys ISA, Color col, string SerName)
        {
            Series SerCor = new Series(SerName)
            {
                ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Column,
                BorderWidth = 1,
                Color = col,
                BorderColor = Color.Black,
                ["PointWidth"] = "1"
            };
            for (double v = 0; v < ISA.Y2.Count; v++)
                SerCor.Points.AddXY(Math.Round(ISA.Min.Q + ISA.Step.Q * v + ISA.Step.Q / 2, 3), ISA.f[(int)v]);
            return SerCor;
        }

        static public void ShowPTreeView(Series[] pPracticalSeries, TreeView EvaluationTreeView)
        {
            EvaluationTreeView.Nodes.Clear();
            for (int i = 0; i < pPracticalSeries.Length / 2; i++)
            {
                EvaluationTreeView.Nodes.Add(new TreeNode()
                {
                    Text = pPracticalSeries[i].Name,
                    ForeColor = pPracticalSeries[i].Color,
                    Checked = true
                });
            }
        }

        static public void ShowEvaluationTreeView(IS.Data.Evaluation Evaluation, TreeView EvaluationTreeView)
        {
            EvaluationTreeView.Nodes.Clear();
            for (int i = 0; i < Evaluation.EvauationSeries.Count / 2; i++)
            {
                EvaluationTreeView.Nodes.Add(new TreeNode()
                {
                    Text = Evaluation.EvauationSeries[i].Name,
                    ForeColor = Evaluation.EvauationSeries[i].Color,
                    Checked = true
                });
            }

        }

        static public void ShowP(Chart EvaluationChart, Series[] pPracticalSeries)
        {
            EvaluationChart.ChartAreas[0].AxisY.Maximum = 0;

            for (int i = 0; i < pPracticalSeries.Length / 2; i++)
            {
                if(pPracticalSeries[i].Points.Count != 0)
                { 
                    EvaluationChart.ChartAreas[0].AxisY.Maximum = Math.Max(EvaluationChart.ChartAreas[0].AxisY.Maximum,
                    pPracticalSeries[i].Points[pPracticalSeries[i].Points.Count - 1].YValues[0] * 1.5);
                    EvaluationChart.Series.Add(pPracticalSeries[i]);
                    EvaluationChart.Series.Add(pPracticalSeries[i + pPracticalSeries.Length / 2]);
                }
            }

        }
        static public void ShowISData(Chart EvaluationChart, IS.Data.Evaluation Evaluation)
        {
            EvaluationChart.ChartAreas[0].AxisY.Maximum = 0;
            
            for (int i = 0; i < Evaluation.Practical.Count; i++)
            {
                EvaluationChart.ChartAreas[0].AxisY.Maximum = 
                    Math.Max(EvaluationChart.ChartAreas[0].AxisY.Maximum,
                    Evaluation.Practical[i].Q[Evaluation.Practical[i].Q.Count - 1].Y * 1.5);
                EvaluationChart.Series.Add(Evaluation.EvauationSeries[i]);
                EvaluationChart.Series.Add(Evaluation.EvauationSeries[i + Evaluation.Practical.Count]);
            }

        }

        static public void ShowQueueAndChannel(ImitationSystem IS, Chart chart)
        {
            Series QueueSer = new Series("В черзі")
            {
                ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line,
                ["PointWidth"] = "3",
                BorderWidth = 2,
                BorderColor = Color.Black
            };
            Series ChannelSer= new Series("Активних каналів")
            {
                ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line,
                ["PointWidth"] = "3",
                BorderWidth = 2,
                BorderColor = Color.Black
            };
            chart.Series.Add(ChannelSer);
            chart.Series.Add(QueueSer);
            ChannelSer.Points.AddXY(0, 1);
            QueueSer.Points.AddXY(0, 0);
            int c = 1;
            int q = 0;
            for (int i = 1; i < IS.TimeLine.Count && IS.TimeLine[i].TimeEvent < 1000 / IS.lambda; i++)
            {
                if (IS.TimeLine[i].InputOutput == EventType.Input)
                {
                    ChannelSer.Points.AddXY(IS.TimeLine[i].TimeEvent, c);
                    c++;
                    ChannelSer.Points.AddXY(IS.TimeLine[i].TimeEvent, c);
                }
                else if (IS.TimeLine[i].InputOutput == EventType.Output)
                {
                    ChannelSer.Points.AddXY(IS.TimeLine[i].TimeEvent, c);
                    c--;
                    ChannelSer.Points.AddXY(IS.TimeLine[i].TimeEvent, c);
                }
                if (IS.TimeLine[i].InputOutput == EventType.Wait && IS.TimeLine[i].TimeProcces != 0)
                {
                    QueueSer.Points.AddXY(IS.TimeLine[i].TimeEvent, q);
                    q++;
                    QueueSer.Points.AddXY(IS.TimeLine[i].TimeEvent, q);
                }
                else if (IS.TimeLine[i].InputOutput == EventType.Wait && IS.TimeLine[i].TimeProcces == 0)
                {
                    QueueSer.Points.AddXY(IS.TimeLine[i].TimeEvent, q);
                    q--;
                    QueueSer.Points.AddXY(IS.TimeLine[i].TimeEvent, q);
                }
            }
        }

        static Series DrawN(String Name, List<Event> N, ImitationSystem IS)
        {
            Series result = new Series(Name)
            {
                ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line,
                ["PointWidth"] = "3",
                BorderWidth = 3,
                BorderColor = Color.Black
            };
            for (int i = 0; i < N.Count; i+=Math.Max(N.Count / 1000, 1))
                result.Points.AddXY(N[i].TimeEvent, N[i].TimeProcces);
            return result;
        }

        static public void ClearDiograms(Chart chart)
        {
            chart.Series.Clear();
        }
    }
}
