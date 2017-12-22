using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IS
{
    enum EventType
    {
        Input,
        Output,
        Wait,
        NFailure,
        NProccesing
    };
    class Event
    {
        public double TimeEvent;
        public double TimeProcces;
        public EventType InputOutput;
        public Event Neighbor;
        public int Number;
        public Event(double TimeEvent, double TimeProcces, EventType InputOutput)
        {
            this.TimeEvent = TimeEvent;
            this.TimeProcces = TimeProcces;
            this.InputOutput = InputOutput;
        }
        public Event(double TimeEvent, double TimeProcces, EventType InputOutput,int Number)
        {
            this.TimeEvent = TimeEvent;
            this.TimeProcces = TimeProcces;
            this.InputOutput = InputOutput;
            this.Number = Number;
        }
    }
}
