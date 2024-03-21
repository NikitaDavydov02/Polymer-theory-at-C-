using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Polymer_brush.Components;

namespace Polymer_brush
{
    [Serializable]
    class Input
    {
        public List<Component> Components;
        public Dictionary<Component, double> VolumeFractionsInTheBulk;
        public double[,] Chi;

        
    }
}
