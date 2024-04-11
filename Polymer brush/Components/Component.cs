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
        [DataMember]
        public double KuhnLength;
        [DataMember]
        public double c;
        [DataMember]
        public int NtotalSegments;
        [DataMember]
        public int NouterSegments;
        [DataMember]
        public List<double> groupsFractions;

    }
    [DataContract]
    public enum ComponentType
    {
        [EnumMember(Value = "Pol")]
        Polymer,
        [EnumMember(Value = "Sol")]
        Solvent,
        [EnumMember(Value = "Add")]
        Additive,
    }
}
