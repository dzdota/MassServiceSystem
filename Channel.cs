using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IS
{
    class Channel
    {
        public Event WorkEvent = new Event(double.MinValue, 0, EventType.Input);
        public List<Event> TimeLine = new List<Event>();

        public bool IsFree(double time)
        {
            if (WorkEvent == null || WorkEvent.Neighbor.TimeEvent < time)
                return true;
            return false;
        }
        public void ChangeEvent(Event NewEvent)
        {
            WorkEvent = NewEvent;
            TimeLine.Add(WorkEvent);
            TimeLine.Add(WorkEvent.Neighbor);
        }
        public double FreeTime
        {
            get
            {
                if (WorkEvent == null || WorkEvent.Neighbor == null)
                    return -1;
                return WorkEvent.Neighbor.TimeEvent;
            }
        }
    }
}
