using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IS
{
    class Channels
    {
        Channel[] elements;
        public int Length;
        public Channels(int n)
        {
            elements = new Channel[n];
            for (int i = 0; i < n; i++)
            {
                elements[i] = new Channel();
            }
            Length = n;
        }

        public Channel this[int index]
        {
            get { return elements[index]; }
            set { elements[index] = value; }
        }
        public Channel FirstFree()
        {
            double min = double.MaxValue;
            Channel res = null;
            foreach (Channel c in elements)
            {
                if ((min = Math.Min(c.FreeTime, min)) == c.FreeTime)
                    res = c;
            }
            return res;
        }
    }
}
