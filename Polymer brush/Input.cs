using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Polymer_brush.Components;
using System.Runtime.Serialization;
using System.Collections.ObjectModel;
using System.IO;
namespace Polymer_brush
{
    [DataContract]
    class Input
    {
        [DataMember]
        public List<Component> Components;
        [DataMember]
        public List<double> VolumeFractionsInTheBulk;
        [DataMember]
        public List<List<double>> Chi;
        [DataMember]
        public Geometry geometry;
        [DataMember]
        public double Radius;
        [DataMember]
        public double DensityDegree;



    }
    [DataContract]
    enum Geometry
    {
        [EnumMember]
        Sphere,
        Plane,
    }
}
