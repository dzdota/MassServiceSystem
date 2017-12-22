using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using TMO;
using IS;
using testgistogr;

namespace TMOLaba4_5
{
    public partial class Form1 : Form
    {
        ImitationSystem IS;
        public Form1()
        {
            InitializeComponent();

        }

        private void Imitationbutton_Click(object sender, EventArgs e)
        {
            int n;
            double lambda, mu, ImitationTime;
            try
            {
                n = Convert.ToInt32(ntextBox.Text);
                lambda = Convert.ToDouble(lambdatextBox.Text.Replace(",","."));
                mu = Convert.ToDouble(mutextBox.Text.Replace(",", "."));
                ImitationTime = Convert.ToDouble(ImitationTimetextBox.Text);
            }
            catch
            {
                MessageBox.Show("Невірно введенні дані", "Помилка", MessageBoxButtons.OK, MessageBoxIcon.Error);
                IS = null;
                return;
            }
            if (int.TryParse(mtextBox.Text, out int m) && checkBox1.Checked)
                IS = new ImitationSystem(lambda, mu, n, m);
            else
                IS = new ImitationSystem(lambda, mu, n);

            IS.Generate(ImitationTime);

            ShowData.ClearDiograms(LambdaChart);
            ShowData.ClearDiograms(MuChart);
            ShowData.ClearDiograms(QueueChannelchart);
            ShowData.ClearDiograms(Nchart);
            ShowData.ClearDiograms(Evaluationchart);
            ShowData.ClearDiograms(Pchart);

            ShowData.ShowDigrams(IS, Nchart);
            ShowData.ShowQueueAndChannel(IS, QueueChannelchart);
            ShowData.ShowISData(Evaluationchart,  IS.Evaluation);
            ShowData.ShowP(Pchart, IS.pPracticalSeries);
            ShowData.ShowLamdaMu(IS, LambdaChart, MuChart);
            
            if (EvaluationtreeView.Nodes.Count == 0)
                ShowData.ShowEvaluationTreeView(IS.Evaluation, EvaluationtreeView);
            else
                EvaluationtreeView_AfterCheck(IS, new TreeViewEventArgs(EvaluationtreeView.Nodes[0]));
            ShowData.ShowPTreeView(IS.pPracticalSeries, treeView1);

            WraitData.Clear(dataGridView1);
            WraitData.WraitP(IS, dataGridView1);
            WraitData.WraitISData(dataGridView1, IS.Evaluation, IS.round);

            GC.Collect();

        }

        private void checkBox1_CheckedChanged(object sender, EventArgs e)
        {
            mtextBox.Visible = checkBox1.Checked;
        }

        private void EvaluationtreeView_AfterCheck(object sender, TreeViewEventArgs e)
        {
            double max = 0;
            for (int i = 0; i < IS.Evaluation.Practical.Count; i++)
            {
                IS.Evaluation.EvauationSeries[i].Enabled = EvaluationtreeView.Nodes[i].Checked;
                IS.Evaluation.EvauationSeries[IS.Evaluation.Practical.Count + i].Enabled = EvaluationtreeView.Nodes[i].Checked;
                if (EvaluationtreeView.Nodes[i].Checked)
                    max = Math.Max(max, IS.Evaluation.Practical[i].Q[IS.Evaluation.Practical[i].Q.Count - 1].Y * 1.5);
            }
            Evaluationchart.ChartAreas[0].AxisY.Maximum = max;
        }

        private void treeView1_AfterCheck(object sender, TreeViewEventArgs e)
        {
            double max = 0;
            for (int i = 0; i < IS.pPracticalSeries.Length / 2; i++)
            {
                IS.pPracticalSeries[i].Enabled = treeView1.Nodes[i].Checked;
                IS.pPracticalSeries[IS.pPractical.Length + i].Enabled = treeView1.Nodes[i].Checked;
                if (treeView1.Nodes[i].Checked)
                    max = Math.Max(max, IS.pTeorical[i] * 1.5);
            }
            Pchart.ChartAreas[0].AxisY.Maximum = max;
        }
    }
}
