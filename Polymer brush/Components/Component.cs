using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Polymer_brush.Components
{
    public class Component
    {
        public ComponentType Type;
        public string Name;
        public int Size;

    }
    public class PolymerComponent : Component
    {
        public double KunLength;
        public double AreaPerChain;
        public double NtotalSegments;
        public double NgroupTypes;
        public List<double> groupsFractions;
    }
    public enum ComponentType
    {
        Polymer,
        Solvent,
        Additive,
    }
}
