using System;
using System.Collections.Generic;
using System.Linq;

namespace IS
{
    class Queue<T> : System.Collections.Generic.Queue<T>
    {
        public List<Event> TimeLine = new List<Event>();
        public List<double> time = new List<double>();


        public void QueueClear(Channels channels, Channel firstFree, Queue<Event> queue, double CorrectTime, ref double NProcessing, ref List<Event> NProcessingEvent)
        {
            while (queue.Count != 0 && firstFree.FreeTime < CorrectTime)
            {
                Event first = queue.Peek();
                firstFree = channels.FirstFree();
                Event EventCome = new Event(Math.Max(firstFree.FreeTime, first.TimeEvent),
                    first.TimeProcces, EventType.Input);
                Event EventLeav = new Event(first.TimeProcces + Math.Max(firstFree.FreeTime, first.TimeEvent),
                    0, EventType.Output);
                EventCome.Neighbor = EventLeav;
                EventLeav.Neighbor = EventCome;

                if (firstFree.FreeTime > first.TimeEvent)
                {
                    time.Add(firstFree.FreeTime - first.TimeEvent);
                    Event WaitCome = new Event(first.TimeEvent + (firstFree.FreeTime - first.TimeEvent) * 0.0001,
                        (firstFree.FreeTime - first.TimeEvent) * 0.9999, EventType.Wait);
                    Event WaitLeave = new Event(firstFree.FreeTime * 0.9999 + first.TimeEvent * 0.0001, 0, EventType.Wait);
                    WaitCome.Neighbor = WaitLeave;
                    WaitLeave.Neighbor = WaitCome;
                    TimeLine.Add(WaitCome);
                    TimeLine.Add(WaitLeave);
                }
                firstFree.ChangeEvent(EventCome);
                queue.Dequeue();
                NProcessing++;
                NProcessingEvent.Add(new Event(Math.Max(firstFree.FreeTime, first.TimeEvent), NProcessing, EventType.NProccesing));
            }
        }
    }
}
