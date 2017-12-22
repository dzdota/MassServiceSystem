using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using IS;
using TMO;


namespace TMOLaba4_5
{
    static class WraitData
    {
        static public void WraitP(ImitationSystem IS, DataGridView datagrid)
        {
            double sumpP = 0, sumpT = 0;
            for (int i = 0; i < IS.pPractical.Length; i++)
            {
                datagrid.Rows.Add("p" + i.ToString(), 
                    Math.Round(IS.pPractical[i], IS.round).ToString(), 
                    Math.Round(IS.pTeorical[i], IS.round).ToString(),
                    Math.Round(IS.pPractical[i] - IS.pTeorical[i], IS.round).ToString());
                sumpP += IS.pPractical[i];
                sumpT += IS.pTeorical[i];
            }

            datagrid.Rows.Add("sum", Math.Round(sumpP, IS.round).ToString(), Math.Round(sumpT, IS.round).ToString(), Math.Round(sumpP - sumpT, IS.round).ToString());
            datagrid.Rows.Add("N відмови", Math.Round(IS.NFailure, IS.round).ToString());
            datagrid.Rows.Add("N оброблених", Math.Round(IS.NProcessing, IS.round).ToString());
            datagrid.Rows.Add("Максимальна кількість запитів в черзі", IS.MaxQueueCount.ToString());
        }

        static public void WraitISData(DataGridView datagrid, IS.Data.Evaluation Information, int round)
        {
            for (int i = 0; i < Information.Teorical.Count; i++)
                datagrid.Rows.Add(
                    Information.Teorical[i].Name,
                    Information.Teorical[i].Q0 == -1 ? "" : Math.Round(Information.Teorical[i].Q0, round).ToString(),
                    Math.Round(Information.Teorical[i].Q, round).ToString(),
                    Information.Teorical[i].Q0 == -1 ? "" : Math.Round(Information.Teorical[i].Q - Information.Teorical[i].Q0, round).ToString());
        }

        static public void Clear(DataGridView datagrid)
        {
            datagrid.Rows.Clear();
        }
    }
}
