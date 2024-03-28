using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.Serialization;

namespace Polymer_brush.Components
{
    [DataContract]
    public class Component
    {
        [DataMember]
        public ComponentType Type;
        [DataMember]
        public string Name;
        [DataMember]
        public int Size;

    }
    [DataContract]
    public class PolymerComponent : Component
    {
        [DataMember]
        public double KunLength;
        [DataMember]
        public double NtotalSegments;
        [DataMember]
        public double NgroupTypes;
        [DataMember]
        public List<double> groupsFractions;
    }
    [DataContract]
    public enum ComponentType
    {
        Polymer,
        Solvent,
        Additive,
    }
}
