using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Runtime.Remoting.Messaging;
using System.Reflection;
using Polymer_brush.Components;
using System.Runtime.Serialization;

namespace Polymer_brush
{
    class Program
    {
        public static double R, BA, y_min, y_max, yacc, aA, rNA, rNB, areaPerChain, Lamb_Pol;
        public static double rN, nu, coe, pi, osmbulk, y_edge, y_cur, P_AB;
        public static double[] chemPotInTheBulk, volumeFractionsInTheBulk, Lagrbulk, size;
        public static double[] chemPotOutsideOfTheStep;
        public static double osmoticPressureOutsideOfTheStep;
        public static double c;
        public static int MaxNumberOfComponent = 50;
        //public static double[] segregationPoints;

        public static double[,] chi;
        public static double[,] etas;
        public static double[] fractionsOfGroups;
        public static int NumberOfComponents;
        public static int NumberOfPolymerGroupTypes;
        public static int chiMatrixSize;
        public static double maxSigma;
        public static double areaDensityDegree = 1;
        static double actualSigma;

        static double[] Xbrush, fipolimer;
        static double point_y;
        static StreamWriter sw;
        static StreamWriter newthonWriter;
        static StreamWriter integralLogWriter;
        static StreamWriter outputWriter;
        static double z = 6;
        //static double[] chemPotInTheBulk;
        static double[] chemPotAtTheBorder;
        static MixingPartModule mixingPartModule;
        static MixingPartModule mixingPartModuleCopy=null;
        static CalculationMode calculationMode;
        static string inputPath = "input.xml";
        static Dictionary<Input,string> LookingForInputFiles()
        {
            Dictionary<Input, string> inputs = new Dictionary<Input, string>();
            string currentDirectoryPath = Environment.CurrentDirectory;
            IEnumerable<string> files = Directory.EnumerateFiles(currentDirectoryPath);
            foreach(string file in files)
            {
                if(file.Contains("input")&& file.Contains(".xml"))
                {
                    Input settings = new Input();
                    DataContractSerializer serializer = new DataContractSerializer(typeof(Input));
                    using (FileStream fs = new FileStream(file, FileMode.Open))
                    {
                        try
                        {
                            settings = serializer.ReadObject(fs) as Input;
                        }
                        catch(Exception e)
                        {
                            throw new Exception("Settings error " + e.Message);
                        }
                    }
                    inputs.Add(settings,file);
                }
            }
            return inputs;
        }
        static void Main(string[] args)
        {
            
            if (File.Exists("output.txt"))
                File.Delete("output.txt");
           

            outputWriter = new StreamWriter(File.Create("output.txt"));
            Dictionary<Input, string>tasks = LookingForInputFiles();
            foreach (Input task in tasks.Keys)
            {
                string outputPath = tasks[task];
                RunTask(task);
                try
                {
                    //RunTask(task);
                }
                catch(Exception ex)
                {
                    outputWriter.WriteLine("Task is aborted due to error: ");
                    outputWriter.WriteLine(ex.Message);
                    Console.WriteLine("Task is aborted");
                    Console.WriteLine("ERROR: " + ex.Message);
                    sw.Close();
                    newthonWriter.Close();
                    integralLogWriter.Close();
                  /*  static StreamWriter sw;
                    static StreamWriter newthonWriter;
                    static StreamWriter integralLogWriter;
                    static StreamWriter outputWriter;*/
                }
            }
            Console.WriteLine("All tasks are done");
            //RunTask(null);
            outputWriter.Close();
            Console.ReadLine();
        }
        static void RunTask(Input task)
        {
            
            outputWriter.WriteLine();
            outputWriter.WriteLine();
            outputWriter.WriteLine();
            outputWriter.WriteLine("/////////////////////////////////////////////////////////////");
            mixingPartModuleCopy = null;
           /* double Fmix = 0;
            double u_sol = 0;
            double u_pol = 0;
            double u_bio = 0;*/
            //calculationMode = CalculationMode.InfinitlyDelute;
            //ReadSettings();

            Enter(task);
            if (volumeFractionsInTheBulk[2] <= 0.001)
                calculationMode = CalculationMode.InfinitlyDelute;
            else
                calculationMode = CalculationMode.Usual;
            //calculationMode = CalculationMode.InfinitlyDelute;
            OutputSettings(task);
            //COMMENT IT
            //CreateInputSettings();
            //return;
            //COMMENT IT

            Console.WriteLine("Initialization successful");
            CalculateMixingSurface(0.005, false);
            CalculateMixingSurface(0.005, true);

            
            if (calculationMode == CalculationMode.InfinitlyDelute)
                SwitchToInfinitlyDeluteMode();

            /////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////
            List<KeyValuePair<double, List<double>>> mixingEnergy = CalculateMixingEnergyProfile(1, 0, 20);
            using (StreamWriter sw = new StreamWriter("mixingFenergyOfSolventAndPolymer.txt"))
            {
                sw.WriteLine("x_pol; Fmix; mixChemPot1; osmotic pressure");

                foreach (KeyValuePair<double, List<double>> pair in mixingEnergy)
                {
                    sw.WriteLine(pair.Key + "  ;  " + pair.Value[0] + ";" + pair.Value[1] + ";" + pair.Value[2] + ";");
                }
            }
            /////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////
            //return;
            sw = new StreamWriter("profile.txt");
            integralLogWriter = new StreamWriter("integral_log_writer.txt");


            y_edge = BisectionSolve(y_min, y_max, yacc, NormalizationFunction, new List<double>());
            y_cur = y_min;
            integralLogWriter.Close();

            //!Calculate polymer concentration at the start position(A / B boundary)
            P_AB = NormalizationSubintegralValue(y_cur, new List<double>() { y_edge });
            //P_AB = FindSubintegralValueForNormalization(y_cur, y_edge); // !   Polymer concentration at A/ B boundary(highest)  must work before finding profile to set Lagr multipliers in common block

            sw.WriteLine("y_cur    solvent    polymer    bio    osm_pressure");

            List<KeyValuePair<double, ProfileInfo>> profile = new List<KeyValuePair<double, ProfileInfo>>();
            
            Console.WriteLine("Calculating profile...");
            y_cur = 1;
            int numberOfPoints = 100;
            double stepInRelativeUnits = (y_edge - 1) / numberOfPoints;
            //while (y_cur < y_edge)
            for (y_cur = 1 + stepInRelativeUnits; y_cur < y_edge; y_cur += stepInRelativeUnits)
            {
                if (Math.Abs(y_cur - 1.02927256377978) < 0.000001)
                    ;
                ProfileInfo info = new ProfileInfo();
                Xbrush = new double[NumberOfComponents];
                FindVolumeFractionsInTheBrushForPoint(out Xbrush, y_cur); //  !calculates concentration profile in the brush after the solution is found
                info.polymerEquationError = brushEquationInfo.polymerEquationError;
                info.polymerEquationMixingPart = brushEquationInfo.polymerEquationMixingPart;
                info.polymerEquationStretchingPart = brushEquationInfo.polymerEquationStretchingPart;
                info.mixingEnergyContributionToF = brushEquationInfo.mixingEnergyContributionToF;
                info.entropyContributionToF = brushEquationInfo.entropyContributionToF;

                for (int i = 0; i < NumberOfComponents; i++)
                {
                    fipolimer[i] = Xbrush[i];
                }

                List<double> composition = new List<double>();
                for (int i = 0; i < NumberOfComponents; i++)
                    composition.Add(fipolimer[i]);
                
                info.composition = composition;
                profile.Add(new KeyValuePair<double, ProfileInfo>(y_cur, info));
                
                // y_cur += aA / R;
            }
            y_cur = y_min;
            double additiveNumberOfMollecules=0; 
            double additiveTotalVolume = 0;
            double additiveFraction = 0;
            double additiveNumberOfCellsAccupied;
            if (calculationMode == CalculationMode.InfinitlyDelute)
            {
                Enter(task);
                Console.WriteLine("*************************");
                Console.WriteLine("*************************");
                Console.WriteLine("Finding additive concentration in infinitly delute solution approximation...");
                Console.WriteLine("Raw profile");
                for (int i = 0; i < profile.Count; i++)
                {
                    string line = profile[i].Key.ToString() + "    ";
                    for (int j = 0; j < profile[i].Value.composition.Count; j++)
                        line += profile[i].Value.composition[j] + "    ";
                    Console.WriteLine(line);
                }
                Console.WriteLine("*************************");
                Console.WriteLine("*************************");
                
                for (int i = 0; i < profile.Count; i++)
                {
                    y_cur = profile[i].Key;
                    double[] baseFractions = new double[2];
                    baseFractions[0] = profile[i].Value.composition[0];
                    baseFractions[1] = profile[i].Value.composition[1];
                    //double additiveFraction = FindAdditiveConcentrationForParticularPolymerAndSolventContentInTheBrush(baseFractions);
                    double error = volumeFractionsInTheBulk[2] * 0.01;
                    error = Math.Pow(10, -17);
                    additiveFraction = BisectionSolve(Math.Pow(10, -18), 0.9, error, InfinitlyDeluteAdditivesEquationsForBisectionSolve, new List<double> { baseFractions[0], baseFractions[1] });

                    profile[i].Value.composition.Add(additiveFraction);
                    additiveTotalVolume += (stepInRelativeUnits * R) * 4 * 3.1415 * (y_cur * R) * (y_cur * R) * additiveFraction;
                }
                
            }
            else
            {
                for (int i = 0; i < profile.Count; i++)
                {
                    y_cur = profile[i].Key;
                    additiveFraction = profile[i].Value.composition[2];
                    additiveTotalVolume += (stepInRelativeUnits * R) * 4 * 3.1415 * (y_cur * R) * (y_cur * R) * additiveFraction;
                   

                }
            }
            additiveNumberOfCellsAccupied = additiveTotalVolume / (aA * aA * aA);
            additiveNumberOfMollecules = additiveNumberOfCellsAccupied / size[2];
            Console.WriteLine("");
            Console.WriteLine("Output");
            outputWriter.WriteLine();
            outputWriter.WriteLine();
            outputWriter.WriteLine("///////////////////////////OUTPUT/////////////////////////////");
            outputWriter.WriteLine("y_cur,nm    y_cur    solvent    polymer    bio    polymerEquationError    polymerEquationMixingPart    polymerEquationStretchingPart    mixingEnergyContributionToF    entropyContributionToF");
            for (int i = 0; i < profile.Count; i++)
            {
                string line = ((profile[i].Key-1)*R*Math.Pow(10,9)) + "    ";
                line+= profile[i].Key.ToString() + "    ";
                for (int j = 0; j < profile[i].Value.composition.Count; j++)
                    line += profile[i].Value.composition[j] + "    ";
                line += profile[i].Value.polymerEquationError + "    ";
                line += profile[i].Value.polymerEquationMixingPart + "    ";
                line += profile[i].Value.polymerEquationStretchingPart + "    ";
                line += profile[i].Value.mixingEnergyContributionToF + "    ";
                line += profile[i].Value.entropyContributionToF + "    ";
                sw.WriteLine(line);
                outputWriter.WriteLine(line + (Osmmix(Program.NumberOfComponents, fipolimer) - osmbulk));
            }

            Console.WriteLine("Beta: " + y_edge);
            Console.WriteLine("Number of molecules: " + additiveNumberOfMollecules);
            Console.WriteLine("Calculation is done!");


            sw.Close();
            outputWriter.WriteLine("Adsorbed molecules: " + additiveNumberOfMollecules);
            outputWriter.WriteLine("Adsorbtion (mol/m2): " + (additiveNumberOfMollecules/(6.02*Math.Pow(10,23)))/(4*3.1415*R*R));
            outputWriter.WriteLine("/////////////////////////////////////////////////////////////");
            outputWriter.WriteLine("**************************************************************");
            outputWriter.WriteLine("**************************************************************");
            outputWriter.WriteLine("**************************************************************");
            outputWriter.WriteLine("**************************************************************");

            
        }
        static void ApplySettings(Input settings)
        {
            if (settings == null)
                return;
            NumberOfComponents = settings.Components.Count;
            Component polymer = null;
            foreach (Component comp in settings.Components)
                if (comp.Type == ComponentType.Polymer)
                    polymer = comp;
            c = polymer.c;
            aA = polymer.KuhnLength;
            rN = polymer.NtotalSegments;
            rNA = polymer.NouterSegments;
            rNB = rN - rNA;
            R = settings.Radius;
            //if (R > aA * rNB && nu == 2)
            //    throw new Exception("These chains are too short to fill the core");
            if (settings.geometry == Geometry.Sphere)
                nu = 2;
            areaDensityDegree = settings.DensityDegree;
            NumberOfPolymerGroupTypes = polymer.groupsFractions.Count;

            size = new double[MaxNumberOfComponent];
            for (int i = 0; i < NumberOfComponents; i++)
                size[i] = settings.Components[i].Size;

            fractionsOfGroups = new double[NumberOfPolymerGroupTypes];
            for (int i = 0; i < NumberOfPolymerGroupTypes; i++)
                fractionsOfGroups[i] = polymer.groupsFractions[i];

            chi = new double[MaxNumberOfComponent, MaxNumberOfComponent];
            chiMatrixSize = NumberOfComponents + NumberOfPolymerGroupTypes - 1;
            for (int i = 0; i < chiMatrixSize; i++)
            {
                for (int j = 0; j < chiMatrixSize; j++)
                    chi[i, j] = settings.Chi[i][j];
            }

            volumeFractionsInTheBulk = new double[NumberOfComponents];
            for (int i = 0; i < NumberOfComponents; i++)
                volumeFractionsInTheBulk[i] = settings.VolumeFractionsInTheBulk[i];
        }
        static void CreateInputSettings()
        {
            if (File.Exists(inputPath))
                return;
            Input settings = new Input();
            settings.Components = new List<Component>();
            //<Solvent>
            Component solvent = new Component();
            solvent.Name = "Solvent";
            solvent.Type = ComponentType.Solvent;
            solvent.Size = 1;
            //</Solvent>
            //<Polymer>
            Component polymer = new Component();
            polymer.Name = "Polymer";
            polymer.Type = ComponentType.Polymer;
            polymer.NtotalSegments = 80;
            polymer.NouterSegments = 60;
            polymer.Size = polymer.NouterSegments;
            polymer.groupsFractions = new List<double>();
            polymer.groupsFractions.Add(1);
            polymer.KuhnLength = aA;
            //</Polymer>
            //<Additive>
            Component additive = new Component();
            additive.Name = "Additive";
            additive.Type = ComponentType.Additive;
            additive.Size = 3;
            //</Additive>
            settings.Components.Add(solvent);
            settings.Components.Add(polymer);
            settings.Components.Add(additive);
            settings.VolumeFractionsInTheBulk = new List<double>();
            settings.VolumeFractionsInTheBulk.Add(0.98);
            settings.VolumeFractionsInTheBulk.Add(0);
            settings.VolumeFractionsInTheBulk.Add(0.02);
            settings.Chi = new List<List<double>>();
            for(int i = 0; i < chiMatrixSize; i++)
            {
                settings.Chi.Add(new List<double>());
                for (int j = 0; j < chiMatrixSize; j++)
                    settings.Chi[i].Add(chi[i, j]);
            }
            settings.geometry = Geometry.Sphere;
            settings.DensityDegree = 1;
            settings.Radius= 3.77 * Math.Pow(10, -8);
            DataContractSerializer serializer = new DataContractSerializer(typeof(Input));
            if(File.Exists(inputPath))
                File.Delete(inputPath);
            using (FileStream fs = new FileStream(inputPath, FileMode.Create))
            {
                serializer.WriteObject(fs, settings);
            }

        }
        static void UploadDefaultSettings()
        {
            c = 2.0;
            pi = 3.1415926;
            coe = (3.0 / 8.0) * pi * pi;
            aA = 6.8 * Math.Pow(10, -9);
            rN = 80.0; //total number of polymer segments per chain
            rNA = 60.0; //number of segments in A - subchain
            rNB = rN - rNA;
            nu = 2.0; //spherical micelle
            R = 3.77 * Math.Pow(10, -8);
            if (R > aA * rNB && nu == 2)
                throw new Exception("These chains are too short to fill the core");
            
            areaDensityDegree = 1;

            NumberOfComponents = 3;
            NumberOfPolymerGroupTypes = 2;

            

            size = new double[MaxNumberOfComponent];
            for (int i = 0; i < NumberOfComponents; i++)
                size[i] = 1.0;

            size[0] = 1.0;// ! solvent
            size[1] = rNA;// polymer
            size[2] = 3.0;// bio//

            chi = new double[MaxNumberOfComponent, MaxNumberOfComponent];
            chiMatrixSize = NumberOfComponents + NumberOfPolymerGroupTypes - 1;

            chi[0, 1] = 0.5;//solv-pol
            chi[0, 2] = 1;//solv-bio
            chi[1, 2] = 0;
            for (int i = 0; i < chiMatrixSize; i++)
            {
                chi[i, i] = 0;
                for (int j = i + 1; j < chiMatrixSize; j++)
                    chi[j, i] = chi[i, j];
            }

            fractionsOfGroups = new double[NumberOfPolymerGroupTypes];
            for (int i = 0; i < NumberOfPolymerGroupTypes; i++)
                fractionsOfGroups[i] = 0.1;
            fractionsOfGroups[0] = 1;
            fractionsOfGroups[1] = 0;

            volumeFractionsInTheBulk = new double[NumberOfComponents];
            for (int i = 0; i < NumberOfComponents; i++)
                volumeFractionsInTheBulk[i] = 0.0;
            volumeFractionsInTheBulk[2] = 0.02;//bio
            volumeFractionsInTheBulk[0] = 0.98;//solvent

            /////////////////////////////////////////////////////////////
            
        }
        static void Enter(Input input)
        {
            Console.WriteLine("Initializing...");
            UploadDefaultSettings();
            ApplySettings(input);
            ////////////////////////////////////
            maxSigma = 1 / (aA * aA);
            actualSigma = maxSigma * areaDensityDegree;
            areaPerChain = 1 / actualSigma;
            y_min = 1.0 + aA / R;
            if (y_min > 1.01)
                y_min = 1.01;
            y_max = 1.00001 * (1.0 + 1.0 * aA * rNA / R);
            //y_max /= 2;
            //y_max = 1.0d0 * (1.0d0 + 2.0 * aA * rNA / R)
            yacc = Math.Pow(10, -8);
            BA = coe / ((rNA * aA) * (rNA * aA));
            chemPotInTheBulk = new double[NumberOfComponents];
            chemPotAtTheBorder = new double[NumberOfComponents];
            fipolimer = new double[NumberOfComponents];

            etas = new double[chiMatrixSize, chiMatrixSize];
            for (int i = 0; i < chiMatrixSize; i++)
                for (int j = 0; j < chiMatrixSize; j++)
                    etas[i, j] = Math.Exp(- chi[i, j] / z);
            ////////////////////////////////////

            if (mixingPartModuleCopy == null)
                mixingPartModule = new MixingPartModule();
            else
                mixingPartModule = mixingPartModuleCopy;

            if (NumberOfComponents == 3)
                volumeFractionsInTheBulk[0] = 1 - volumeFractionsInTheBulk[1] - volumeFractionsInTheBulk[2];
            Lagrmix(NumberOfComponents, volumeFractionsInTheBulk, out chemPotInTheBulk);
           
        }
        static void SwitchToInfinitlyDeluteMode()
        {
            Console.WriteLine("Switching to infinitly delute mode");
            NumberOfComponents = 2;
            volumeFractionsInTheBulk = new double[2];
            volumeFractionsInTheBulk[1] = 0;
            volumeFractionsInTheBulk[0] = 1;
            mixingPartModuleCopy = mixingPartModule;
            mixingPartModule = new MixingPartModule();
            Lagrmix(NumberOfComponents, volumeFractionsInTheBulk, out chemPotInTheBulk);
        }
        static void Lagrmix(int numberOfComponents, double[] X, out double[] mixingPartOfExchangeChemicalPotentials)
        {
            mixingPartOfExchangeChemicalPotentials = new double[NumberOfComponents];
            double[] AlternativemixingPartOfExchangeChemicalPotentials = new double[NumberOfComponents];
            for (int i = 0; i < numberOfComponents; i++)
            {
                /*double sum = 0;
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[j] * (chi[i, j] - chi[0, j]);
				mixingPartOfExchangeChemicalPotentials[i] = (Math.Log(X[i]) + 1.0) /  size[i] - (Math.Log(X[0]) + 1.0) /  size[0] + sum;// ! dummy for solvent identically 0
				*/
                if (X[i] == 0)
                    mixingPartOfExchangeChemicalPotentials[i] = 0;
                else
                    mixingPartOfExchangeChemicalPotentials[i] = mixingPartModule.CalculateExchangeChemialPotentialOfComponent(X, i);                                                                                   //AlternativemixingPartOfExchangeChemicalPotentials(i) = CalculateExchangeChemialPotentialOfComponent(3, X, i)
            }
        }
        static void Lagrmix_PolA(int numberOfComponents, double[] X, out double[] mixingPartOfExchangeChemicalPotentials)
        {
            mixingPartOfExchangeChemicalPotentials = new double[NumberOfComponents];
            double[] AlternativemixingPartOfExchangeChemicalPotentials = new double[NumberOfComponents];

            for (int i = 0; i < numberOfComponents; i++)
            {
                /*if (X[i] < 0)
				{
					for (int n = 0; n < numberOfComponents; n++)
						mixingPartOfExchangeChemicalPotentials[n] = 1050;
					return;
				}

				double sum = 0;
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[j] * (chi[i, j] - chi[0, j]);
				if (i == numberOfComponents - 1)
					mixingPartOfExchangeChemicalPotentials[i] = -(Math.Log(X[0]) + 1.0) /  size[0] + sum;// ! dummy for solvent identically 0
				else
					mixingPartOfExchangeChemicalPotentials[i] = (Math.Log(X[i]) + 1.0) /  size[i] - (Math.Log(X[0]) + 1.0) /  size[0] + sum;// ! dummy for solvent  identically 0
				*/
                AlternativemixingPartOfExchangeChemicalPotentials[i] = mixingPartModule.CalculateExchangeChemialPotentialOfComponent(X, i);
                mixingPartOfExchangeChemicalPotentials[i] = AlternativemixingPartOfExchangeChemicalPotentials[i];
            }
        }
        static double Osmmix(int numberOfComponents, double[] X)
        {
            double sum = 0;
            double sum1 = 0;
            for (int i = 0; i < numberOfComponents; i++)
            {
                sum1 = sum1 + X[i] * (1.0 - 1.0 / size[i]);
                for (int j = 0; j < numberOfComponents; j++)
                    sum += X[i] * X[j] * (chi[0, j] - chi[i, j] / 2.0);
            }
            return Math.Log(X[0]) + sum1 + sum;
        }
        static double NormalizationFunction(double y, List<double> parameters)
        {

            double nu = 2.0;
            double norm = rNA / (rNB * (nu + 1.0));
            //norm = actualSigma * rNA * aA * aA * aA / R;
            norm = (aA / R) * rNA * areaDensityDegree;
            double integrationMin = y_min;
            double integrationMax = y;
            y_cur = y;
            double s = 1;


            //Find chemical potentials at border;
            for (int i = 0; i < NumberOfComponents; i++)
                chemPotAtTheBorder[i] = chemPotInTheBulk[i];
            //chemPotAtTheBorder[0] = chemPotInTheBulk[0];//solvent
            //chemPotAtTheBorder[1] = chemPotInTheBulk[1];//bio
            //Finding volume fractions at the border
           
            //double ERREL = Math.Pow(10, -6);

            double[] XBorderGUESS = new double[NumberOfComponents - 1];//Everything except polymer

            int componentInTheBulkNumber = 0;
            for (int i = 0; i < NumberOfComponents; i++)
                if (i != 1)
                {
                    XBorderGUESS[componentInTheBulkNumber] = volumeFractionsInTheBulk[i];
                    componentInTheBulkNumber++;
                }
            //this is the fraction of solvent at the border
            //XBorderGUESS[1]= volumeFractionsInTheBulk[1];//this is the fraction of biocomponnt at the border

            if(c==0.0)
                XBorderGUESS[0] = 0.01;
            else
                XBorderGUESS[0] = 0.99;
            XBorderGUESS[0] = 0.98;
            // XBorderGUESS[0] = 0.01;
            double FNORM;
            double[] _XBorder = new double[NumberOfComponents-1];
            double ERREL = Math.Pow(10, -3);
            ////////////
            //solvent
            //additive
            ////////////
            DNEQNF(BorderEquations, ERREL, NumberOfComponents - 1, 1000, XBorderGUESS, out _XBorder, out FNORM);

            double[] volFractionsAtTheBorder = new double[NumberOfComponents];
            double volumeFractionSum = 0;

            componentInTheBulkNumber = 0;
            for (int i = 0; i < NumberOfComponents; i++)
            {
                if (i != 1)
                {
                    volFractionsAtTheBorder[i] = _XBorder[componentInTheBulkNumber];
                    volumeFractionSum += volFractionsAtTheBorder[i];
                    componentInTheBulkNumber++;
                }
            }
            volFractionsAtTheBorder[1] = 1 - volumeFractionSum;
            
            Lagrmix_PolA(NumberOfComponents, volFractionsAtTheBorder, out chemPotAtTheBorder);

            double mixingChemPotAtTheBorder = chemPotAtTheBorder[1];
            chemPotAtTheBorder[1] += BA * (R * (integrationMax - 1)) * (R * (integrationMax - 1));
            Lamb_Pol = chemPotAtTheBorder[1];

            CalculateIntegral(integrationMin, integrationMax, NormalizationSubintegralValue, new List<double>() { y }, out s);
            return s - norm;

        }
        static void CalculateIntegral(double integrationMin, double integrationMax, Func<double, List<double>, double> func, List<double> parameters, out double s)
        {
            double EPS = 0.1;
            int JMAX = 15;
            double old_s = -1 * Math.Pow(10, -30);
            s = 0;
            integralLogWriter.WriteLine("Calculate integral for y_max=" + integrationMax);
            for (int n = 1; n < JMAX; n++)
            {
                s = CalculateIntegralWithDefeniteNumberOfTrapezoids(integrationMin, integrationMax, s, n, func, parameters);
                //s = CalculateNormalizationIntegralWithDefeniteNumberOfTrapezoids(integrationMin, integrationMax, s, n+1);

                if (n > 4)
                    if (Math.Abs(s - old_s) < EPS * Math.Abs(old_s) || (s == 0 && old_s == 0))
                    {
                        //integralLogWriter.Close();
                        return;
                    }
                old_s = s;
            }
            Console.WriteLine("Too many steps in q trap");
            integralLogWriter.Close();

            throw new Exception();


        }
        static string func_log = "";
        static double CalculateIntegralWithDefeniteNumberOfTrapezoids(double integrationMin, double integrationMax, double s, int n, Func<double, List<double>, double> func, List<double> parameters)
        {
            integralLogWriter.WriteLine("    Sum Trapezoinds " + n);
            //s = 0;
            if (n == 1)
            {
                integralLogWriter.WriteLine("        " + "Iter" + " ; " + "x" + " ; " + "add" + "; " + "sum;" );

                double f_min = func(integrationMin, parameters);
                double f_max = func(integrationMax, parameters);
                s = 0.5 * (integrationMax - integrationMin) * (f_min+f_max);
                integralLogWriter.WriteLine("        " + 1 + " ; " + (integrationMin+ integrationMax)/2 + " ; " + s + "; " + s + "; "+func_log);

            }
            else
            {
                if (Math.Abs(integrationMax - 1.0294357)<0.00001)
                    ;
                int it = (int)Math.Pow(2, n - 2);
                double tnm = it;
                double del = (integrationMax - integrationMin) / tnm;
                double x = integrationMin + 0.5 * del;
                double sum = 0;
                integralLogWriter.WriteLine("        " + "Iter" + " ; " + "x" + " ; " + "add" + "; " + "sum");
                for (int counter = 0; counter < it; counter++)
                {
                    //double add = FindSubintegralValueForNormalization(x, integrationMax);
                    
                    double add = func(x, parameters);
                    if (double.IsNaN(add))
                        ;
                    sum += add;
                    integralLogWriter.WriteLine("        " + counter + " ; " + x + " ; " + add + "; " + sum + "; " + func_log+";");

                    x += del;
                    
                }
                s = 0.5 * (s + (integrationMax - integrationMin) * sum / tnm);
                integralLogWriter.WriteLine("        "+s);
            }
            return s;
        }

        static double NormalizationSubintegralValue(double x, List<double> parameters)
        {
            FindVolumeFractionsInTheBrushForPoint(out Xbrush, x);
            func_log = Xbrush[1].ToString()+"___;";
            
            double fay = Xbrush[1];
            if (fay <= 0)
                fay = fay;
            return fay * x * x;

        }
        static double BisectionSolve(double y_min, double y_max, double yacc, Func<double, List<double>, double> func, List<double> parameters)
        {
            double output = 0;
            int j_max = 100; //!Maximum allowed number of bisections.
                             //!Using bisection, find the root of a function func known to lie between x1 and x2.The
                             //!root, returned as rtbis, will be refined until its accuracy is ï¿½xacc.
            double dx, f, fmid, xmid;
            //!calculating function value at y_max
            fmid = func(y_max, parameters);
            f = func(y_min, parameters); // !R2(beta2, chemPotInTheBulk, itog2)

            if (f * fmid > 0)
            {
                Console.WriteLine("Bisection error");
                integralLogWriter.Close();
                throw new Exception();
                return 0;
            }
            if (f < 0)
            {
                // then!Orient the search so that f> 0 lies at x + dx.
                output = y_min;
                dx = y_max - y_min;
            }
            else
            {
                output = y_max;
                dx = y_min - y_max;

            }
            for (int j = 0; j < j_max; j++)
            {
                dx = dx / 2;
                xmid = output + dx;
                fmid = func(xmid, parameters); ;// !call R3(beta3, chemPotInTheBulk, itog3)
                if (fmid < 0)
                    output = xmid;
                if (Math.Abs(dx) < yacc || fmid == 0)
                    return output;
            }
            Console.WriteLine("Too many iterations");
            integralLogWriter.Close();
            throw new Exception();
        }
        static void FindVolumeFractionsInTheBrushForPoint(out double[] XBrush, double y_cur)
        {
            //!Solving(localy) system of non - linear equations
            double ERREL = Math.Pow(10, -6);
            point_y = y_cur;
            int ITMAX = 600;
            double[] XBrushGUESS = new double[NumberOfComponents - 1];//Everything except solvent
            double[] XBrushReduced = new double[NumberOfComponents - 1];//Everything except solvent
            XBrush = new double[NumberOfComponents];

            //<OLD>
            XBrushGUESS[0] = 0.97;//this is the fraction of polymer in the brush
            for(int i=1;i<XBrushGUESS.Length;i++)
                XBrushGUESS[i] = Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
             //<OLD>

            /*XBrushGUESS[0] = 0.01;
            if (NumberOfComponents > 2)
                XBrushGUESS[1] = volumeFractionsInTheBulk[2];*/

            //XBrushGUESS[0] = 9*Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
            //XBrushGUESS[1] = 0;//this is the fraction of polymer in the brush
            double FNORM;
            if (Math.Abs(point_y - 10) < 0.01)
                ;
            DNEQNF(BrushEquations, ERREL, NumberOfComponents - 1, ITMAX, XBrushGUESS, out XBrushReduced, out FNORM);

            double volumeFractionsSum = 0;
            for (int i = 1; i < NumberOfComponents; i++)
            {
                XBrush[i] = XBrushReduced[i-1];
                volumeFractionsSum += XBrush[i];
            }
            XBrush[0] = 1 - volumeFractionsSum;
        }
        delegate string NonlinearSystem(double[] X, out double[] F, int L);

        
        static void DNEQNF(NonlinearSystem Func, double ERREL, int L, int ITMAX, double[] XGuess, out double[] X, out double FNORM)
        {
            int splitTransitions = 0;
            newthonWriter = new StreamWriter(File.Create("newthon_log.txt"));
            newthonWriter.WriteLine("Iteration; X0; F0; FNORM; J00;");
            int iterations = 0;
            double[] F = new double[L];
           
            X = new double[XGuess.Length];
            double[] oldX = new double[XGuess.Length];
            double[] deltaX = new double[XGuess.Length];
            for (int i = 0; i < XGuess.Length; i++)
            {
                X[i] = XGuess[i];
                deltaX[i] = 0;
                oldX[i] = X[i];
            }
            double[] realComposition = new double[NumberOfComponents];
            FNORM = 0;
            double[,] J = new double[L, XGuess.Length];

            do
            {
                for (int i = 0; i < X.Length; i++)
                    X[i] += deltaX[i];
                //Post-step treatment

                if (oldX.Sum() >= 1 && L == 2)
                {
                    deltaX[1] = -deltaX[0];
                    for (int j = 0; j < X.Length; j++)
                        X[j] = oldX[j] + deltaX[j];
                   
                }
                for (double devisionStepDegree = 1; ContainsOutrangeValues(X); devisionStepDegree++)
                {
                    for (int j = 0; j < X.Length; j++)
                    {
                        deltaX[j] /= 2;
                        X[j] = oldX[j] + deltaX[j];
                    }
                }
                

                //<Split>

                double segregationDelta;
                if (Func.Method.Name == "BorderEquations")
                {
                    double[] composition = new double[3];
                    composition[0] = X[0];
                    if (NumberOfComponents == 3)
                        composition[2] = X[1];
                    composition[1] = 1 - composition[0] - composition[2];
                    bool inside = mixingPartModule.IsCompositionInsideSegregationZone(composition, out segregationDelta);
                    if (inside)
                    {
                        X[0] = mixingPartModule.Nodes[0].secondComposition[0] - 0.0000001;
                        newthonWriter.WriteLine("Split");
                        newthonWriter.Close();
                        return;
                    }
                }
                if (Func.Method.Name == "BrushEquations")
                {
                    double sum = 0;
                    for (int i = 0; i < XGuess.Length; i++)
                    {
                        realComposition[i + 1] = X[i];
                        sum += X[i];
                    }
                    realComposition[0] = 1 - sum;
                }
                double[] firstSegregationPoint;
                double[] secondSegregationPoint;
                if (Func.Method.Name == "BrushEquations" && mixingPartModule.IsCompositionInsideSegregationZone(realComposition, out segregationDelta, out firstSegregationPoint, out secondSegregationPoint))
                {
                    if (NumberOfComponents != 2)
                    {
                        //.X[0] = mixingPartModule.Nodes[0].secondComposition[1] + 0.0000001;
                        if (deltaX[0] > 0)
                            X[0] = secondSegregationPoint[1] - 0.01;
                        if (deltaX[0] < 0)
                            X[0] = firstSegregationPoint[1] + 0.01;
                        //newthonWriter.Close();
                        //return;
                        //throw new NotImplementedException();
                    }
                    else
                    {
                        //Split
                        if (deltaX[0] > 0)
                        {
                            if (mixingPartModule.Nodes[0].secondComposition[1] < 0.99)
                                X[0] = mixingPartModule.Nodes[0].secondComposition[1] + 0.01;
                            else
                                X[0] = mixingPartModule.Nodes[0].secondComposition[1] + (1 - mixingPartModule.Nodes[0].secondComposition[1]) * 0.5;
                        }
                        if (deltaX[0] < 0)
                        {
                            if (mixingPartModule.Nodes[0].firstComposition[1] > 0.1)
                                X[0] = mixingPartModule.Nodes[0].firstComposition[1] - 0.1;
                            else
                                X[0] = mixingPartModule.Nodes[0].firstComposition[1] * 0.5;
                        }

                        splitTransitions++;
                        if (splitTransitions > 10)
                        {
                            //Split
                            newthonWriter.WriteLine("Split");
                            newthonWriter.Close();
                            //if (NumberOfComponents != 2)
                            //    throw new NotImplementedException();
                            //X[0] = firstSegregationPoint[1];
                            X[0] = mixingPartModule.Nodes[0].secondComposition[1];
                            return;
                        }
                    }
                    
                }

                //</Split>
                iterations++;
                if (iterations > ITMAX)
                {
                    newthonWriter.Close();
                    return;
                    integralLogWriter.Close();
                    throw new Exception(" Newton method did not manage to find solution for system of equations");

                }

                string funcLog = Func(X, out F, L);
                FNORM = 0;
                for (int i = 0; i < L; i++)
                    FNORM += F[i] * F[i];
                //////////////////////////////////LOG///////////////////////////////////
                string newthonWriterLine = "";
                for (int i = 0; i < L; i++)
                    newthonWriterLine += X[i] + ";";
                for (int i = 0; i < L; i++)
                    newthonWriterLine += F[i] + ";";
                newthonWriterLine += FNORM + ";";
                for (int i = 0; i < L; i++)
                    for (int j = 0; j < L; j++)
                        newthonWriterLine += J[i, j] + ";";
                newthonWriter.WriteLine(iterations + ";" + newthonWriterLine + funcLog);
                ///////////////////////////////////////////////////////////////////////////

                /*Func(X, out F, L);
                for (int i = 0; i < L; i++)
                    FNORM += F[i] * F[i];*/
                /////////////////////////////////JACOBIAN CALCULATION///////////////////////////////////////
                //Calculate Jacobian
                for (int i = 0; i < L; i++)
                    for (int j = 0; j < X.Length; j++)
                    {
                        Func(X, out F, L);
                        double f_init = F[i];
                        double old_x = X[j];
                        double dx = old_x * 0.01;
                        if (dx < Math.Pow(10, -10))
                            dx = Math.Pow(10, -10);
                        X[j] += dx;


                        for (double devisionStepDegree = 1; X[j] >= 1; devisionStepDegree++)
                        {
                            dx /= 2;
                            X[j] = dx + old_x;
                        }
                        if (X.Sum() > 1)
                        {
                            X[j] = old_x - dx;
                            dx *= -1;
                        }

                        Func(X, out F, L);
                        double f_df = F[i];
                        J[i, j] = (f_df - f_init) / dx;
                        X[j] = old_x;
                        Func(X, out F, L);
                    }
                //J=CalculateJacobian(Func, X);
                double det = Matrix.determinantGauss(L, J, 0, false);
                for (int j = 0; j < X.Length; j++)
                    oldX[j] = X[j];
                if (det == 0)
                {
                    //Node transition
                    newthonWriter.Close();
                    if (X[0] < Math.Pow(10, -18))
                    {
                        return;
                    }
                    integralLogWriter.Close();
                    throw new Exception();
                }


                double[,] reverse = Matrix.reverseMatrix(L, J);
                double[] rightParts = new double[L];
                for (int i = 0; i < L; i++)
                {
                    rightParts[i] = -F[i];
                    for (int j = 0; j < X.Length; j++)
                        rightParts[i] += J[i, j] * X[j];
                }
                if (iterations > 100)
                    ;
                X = Matrix.multiplyMatrixAndVector(L, L, reverse, rightParts);
                for (int i = 0; i < X.Length; i++)
                {
                    deltaX[i] = X[i] - oldX[i];
                    X[i] = oldX[i];
                }
                
            }
            while (!Converged(deltaX,ERREL));

            
            newthonWriter.WriteLine("Converged!");
            newthonWriter.Close();
            
        }
        static bool Converged(double[] deltaX, double dX)
        {
            foreach (double x in deltaX)
                if (Math.Abs(x) > dX)
                    return false;
            return true;
        }
        static bool ContainsOutrangeValues(double[] X)
        {

            foreach (double x in X)
                if (x < 0 || x > 1)
                    return true;
            double sum = 0;
            foreach (double x in X)
                sum += x;
            if (sum > 1)
                return true;

            return false;
        }
        static string BorderEquations(double[] X, out double[] F, int L)
        {
            string logString = "";
            double[] _volumeFractions = new double[NumberOfComponents];
            double volumeFractionSum = 0;

            int numberOfComponentAtTheBorder = 0;
            for (int i = 0; i < NumberOfComponents; i++)
            {
                if (i != 1)
                {
                    _volumeFractions[i] = X[numberOfComponentAtTheBorder];
                    volumeFractionSum += X[numberOfComponentAtTheBorder];
                    numberOfComponentAtTheBorder++;
                }
                
            }
            _volumeFractions[1] = 1 - volumeFractionSum;
            F = new double[L];
            double[] mixingPartOfExchangeChemicalPotentials;
            Lagrmix_PolA(NumberOfComponents, _volumeFractions, out mixingPartOfExchangeChemicalPotentials);
            double osmoticPressure = CalculateOsmoticPressure(_volumeFractions);
            //F[0] = (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]) * (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]);//solvent
            //F[1] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]);//bio
            F[0] = osmoticPressure * osmoticPressure;
            F[0] = osmoticPressure;
            for (int i = 2; i < NumberOfComponents; i++)
                F[i-1] = (mixingPartOfExchangeChemicalPotentials[i] - chemPotAtTheBorder[i]) * (mixingPartOfExchangeChemicalPotentials[i] - chemPotAtTheBorder[i]);//bio
            for (int i = 2; i < NumberOfComponents; i++)
                F[i - 1] = (mixingPartOfExchangeChemicalPotentials[i] - chemPotAtTheBorder[i]);//bio

            return logString;
        }
        static double targetMixingPotential;
        static ProfileInfo brushEquationInfo;
        static string BrushEquations(double[] X, out double[] F, int L)
        {
            string logString = "";
            double nu = 2;
            double y_cur = point_y;
            if (y_cur == 2.7661752970011535)
                ;

            double volumeFractionSum = 0;
            for (int i = 1; i < NumberOfComponents; i++)
            {
                fipolimer[i] = X[i-1];
                volumeFractionSum += fipolimer[i];
            }
            fipolimer[0] = 1 - volumeFractionSum;
            double[] mixingPartOfExchangeChemicalPotentials;
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            Lagrmix_PolA(NumberOfComponents, fipolimer, out mixingPartOfExchangeChemicalPotentials);
            F = new double[L];

            //F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
            targetMixingPotential = -(BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol);
            ;
            for (int i = 1; i < NumberOfComponents; i++)
            {
                if (i == 1)
                    F[i - 1] = Math.Pow((mixingPartOfExchangeChemicalPotentials[1] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol), 2);//!polymer error
                else
                    F[i-1] = (mixingPartOfExchangeChemicalPotentials[i] - chemPotInTheBulk[i]) * (mixingPartOfExchangeChemicalPotentials[i] - chemPotInTheBulk[i]);// !bio contaminant error

                if (i == 1)
                    F[i - 1] = Math.Pow((mixingPartOfExchangeChemicalPotentials[1] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol), 1);//!polymer error
                else
                    F[i - 1] = (mixingPartOfExchangeChemicalPotentials[i] - chemPotInTheBulk[i]);// !bio contaminant error

            }
            //F[NumberOfComponents - 2] = Math.Pow((mixingPartOfExchangeChemicalPotentials[NumberOfComponents - 1] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol), 2);//!polymer error
            brushEquationInfo = new ProfileInfo();
            brushEquationInfo.polymerEquationError = F[0];
            brushEquationInfo.polymerEquationMixingPart = mixingPartOfExchangeChemicalPotentials[1];
            brushEquationInfo.polymerEquationStretchingPart = BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0));
            mixingPartModule.CalculateMixingFreeEnergy(fipolimer);
            brushEquationInfo.entropyContributionToF = mixingPartModule.entropyPart;
            brushEquationInfo.mixingEnergyContributionToF = mixingPartModule.mixingPart;

            string line = "";
            for (int i = 0; i < F.Length; i++)
                line += F[0] + "  ";
            string line2 = "";
            for (int i = 0; i < X.Length; i++)
                line2 += X[0] + "  ";

            Console.WriteLine("BrushEquations values: " + line + " ----------------Volume fractions:   " + line2);
            logString += mixingPartOfExchangeChemicalPotentials[1] + ";";
            logString += targetMixingPotential + ";";
            return logString;
        }
        static double CalculateOsmoticPressure(double[] _volumeFractions)
        {
            double[] mixingPartOfExchangeChemicalPotentials;
            Lagrmix_PolA(NumberOfComponents, _volumeFractions, out mixingPartOfExchangeChemicalPotentials);
            //IS IT CORRECT?
            //mixingPartOfExchangeChemicalPotentials[1]+= BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0));
            //IS IT CORRECT?
            double osmoticPressure = 0;
            osmoticPressure = mixingPartModule.CalculateMixingFreeEnergy(_volumeFractions);
            for (int i = 1; i < NumberOfComponents; i++)
                osmoticPressure -= _volumeFractions[i] * mixingPartOfExchangeChemicalPotentials[i];
            return osmoticPressure;
        }
        static void CalculateMixingSurface(double step, bool linearized)
        {
            Console.WriteLine("Calculating surface");
           // List<KeyValuePair<Vector3, double>> output = new List<KeyValuePair<Vector3, double>>();
            if (NumberOfComponents != 3)
                return;
            double x3;
            double[] fractions = new double[3];
            string fileName = "";
            if (linearized)
                fileName = "surface_linearized.txt";
            else
                fileName = "surface_nonlinearized.txt";
            using(StreamWriter surfaceWriter = new StreamWriter(fileName))
            {
                surfaceWriter.WriteLine("Solvent;Polymer;Bio;F;isSegreagated; osmotic;");
                for (double x1 = 0; x1 <= 1; x1+=step)
                {
                    fractions[0] = x1;
                    for (double x2 = 0; x2 <= 1; x2 += step)
                    {
                        x3 = 1 - x1 - x2;
                        if (x3 <= 0)
                            continue;
                        fractions[1] = x2;
                        fractions[2] = x3;
                        double F = mixingPartModule.CalculateMixingFreeEnergy(fractions,linearized);
                        double osmoticPressure = CalculateOsmoticPressure(fractions);
                        double a;
                        bool segregated = mixingPartModule.IsCompositionInsideSegregationZone(fractions, out a);
                        if(linearized)
                            surfaceWriter.WriteLine(x1 + ";" + x2 + ";" + x3 + ";" + F + ";"+ segregated+";" + osmoticPressure + ";");
                        else
                            surfaceWriter.WriteLine(x1 + ";" + x2 + ";" + x3 + ";" + F + ";");
                        //output.Add(new KeyValuePair<Vector3, double>(new Vector3(x1,x2,x3), F));
                    }
                }
            }
            
        }
        static List<KeyValuePair<double, List<double>>> CalculateMixingEnergyProfile(int AcomponentIndex, int BcomponentIndex, int numberOfPoints)
        {
            List<KeyValuePair<double, List<double>>> output = new List<KeyValuePair<double, List<double>>>();
            double step = 0.99 / (numberOfPoints - 1);
            step = 0.0001;
            double[] volumeFractions = new double[NumberOfComponents];
            //double volumeFractionSum = 0;
            for (int i = 0; i < NumberOfComponents; i++)
                volumeFractions[i] = 0;
            volumeFractions[AcomponentIndex] = 0.005;
            volumeFractions[BcomponentIndex] = 1 - volumeFractions[AcomponentIndex];
            while (volumeFractions[AcomponentIndex] < 1)
            {
                List<double> value = new List<double>();
                if (volumeFractions[AcomponentIndex] > 0.9506)
                    ;
                double Fmix = mixingPartModule.CalculateMixingFreeEnergy(volumeFractions);
                double osmotic = CalculateOsmoticPressure(volumeFractions);
                double[] exchangeChemPotentials;
                Lagrmix(NumberOfComponents, volumeFractions, out exchangeChemPotentials);
                value.Add(Fmix);
                value.Add(exchangeChemPotentials[1]);
                value.Add(osmotic);
                output.Add(new KeyValuePair<double, List<double>>(volumeFractions[AcomponentIndex], value));
                volumeFractions[AcomponentIndex] += step;
                volumeFractions[BcomponentIndex] = 1 - volumeFractions[AcomponentIndex];
            }

            return output;
        }
       // static double[] chemPotAtTheBaseBrush;
        static double[] compositionAtTheBaseBrush;
       /* static double FindAdditiveConcentrationForParticularPolymerAndSolventContentInTheBrush(double[] SolventAndPolymerFractions)
        {
            Console.WriteLine("Finding Additive Concentration...");
            compositionAtTheBaseBrush = new double[3];
            compositionAtTheBaseBrush[0] = SolventAndPolymerFractions[0];
            compositionAtTheBaseBrush[1] = SolventAndPolymerFractions[1];
            compositionAtTheBaseBrush[2] =0;

            //Lagrmix_PolA(NumberOfComponents, compositionAtTheBaseBrush, out chemPotAtTheBaseBrush);

            //double ERREL = Math.Pow(10, -2) * volumeFractionsInTheBulk[2];
            double ERREL = Math.Pow(10, -6);
            int ITMAX = 600;
            double[] tryAdditive = new double[] { volumeFractionsInTheBulk[2] };
            double[] additiveConcentrations;
            double FNORM;
            DNEQNF(InfinitlyDeluteAdditivesEquations, ERREL, 1, ITMAX, tryAdditive, out additiveConcentrations, out FNORM);

            return additiveConcentrations[0];
        }*/
        static string InfinitlyDeluteAdditivesEquations(double[] X, out double[] F, int L)
        {
            Console.WriteLine("Infinitly delute additive equation");
            //X - only additives
            string logString = "";
            double nu = 2;
            double y_cur = point_y;

            double[] tryComposition = new double[NumberOfComponents];
            tryComposition[0] = (1 - X[0]) * compositionAtTheBaseBrush[0];
            tryComposition[1] = (1 - X[0]) * compositionAtTheBaseBrush[1];
            tryComposition[2] = X[0];
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            double[] chemPotInTheComplementaryBrush;
            Lagrmix_PolA(NumberOfComponents, tryComposition, out chemPotInTheComplementaryBrush);
            F = new double[L];
            for (int i = 0; i < L; i++)
                F[i] = (chemPotInTheComplementaryBrush[i + 2] - chemPotInTheBulk[i + 2]);
            return logString;
        }
        static double InfinitlyDeluteAdditivesEquationsForBisectionSolve(double xAdditive, List<double> parameters)
        {
            Console.WriteLine("Infinitly delute additive equation (Bisection solve)");
            //X - only additives
            string logString = "";
            //double nu = 2;
            //double y_cur = point_y;

            double[] tryComposition = new double[NumberOfComponents];
            tryComposition[0] = (1 - xAdditive) * parameters[0];
            tryComposition[1] = (1 - xAdditive) * parameters[1];
            tryComposition[2] = xAdditive;
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            double[] chemPotInTheComplementaryBrush;
            Lagrmix_PolA(NumberOfComponents, tryComposition, out chemPotInTheComplementaryBrush);
            double output = (chemPotInTheComplementaryBrush[2] - chemPotInTheBulk[2]);
            return output;
        }
        static void OutputSettings(Input settings)
        {
            if (settings == null)
                return;
            outputWriter.WriteLine("///////////////////////////INPUT////////////////////////////");
            outputWriter.WriteLine("Calculation mode: " + calculationMode);
            outputWriter.WriteLine("Number of components: " + settings.Components.Count);

            Component polymer = null;
            foreach (Component comp in settings.Components)
                if (comp.Type == ComponentType.Polymer)
                    polymer = comp;

            outputWriter.WriteLine("c: " + polymer.c);
            outputWriter.WriteLine("aA: " + polymer.KuhnLength);
            outputWriter.WriteLine("rN: " + polymer.NtotalSegments);
            outputWriter.WriteLine("rNA: " + polymer.NouterSegments);
            outputWriter.WriteLine("Radius (nm): " + R*Math.Pow(10,9));
            if (settings.geometry == Geometry.Sphere)
                outputWriter.WriteLine("nu: " + 2);
            outputWriter.WriteLine("areaDensityDegree: " + settings.DensityDegree);
            outputWriter.WriteLine("Number of polymer groups: " + polymer.groupsFractions.Count);

            int a = 0;
            foreach(Component comp in settings.Components)
            {
                outputWriter.WriteLine();
                outputWriter.WriteLine("Name: " + comp.Name);
                outputWriter.WriteLine("Size: " + comp.Size);
                outputWriter.WriteLine("Bulk fraction: " + settings.VolumeFractionsInTheBulk[a]);
                a++;
            }
            outputWriter.WriteLine();

            for (int i = 0; i < NumberOfPolymerGroupTypes; i++)
                outputWriter.WriteLine("Polmer group fraaction["+i+"]: " + polymer.groupsFractions[i]);

            outputWriter.WriteLine();
            outputWriter.WriteLine("Chi matrix");
            for (int i = 0; i < chiMatrixSize; i++)
            {
                string line = "";
                for (int j = 0; j < chiMatrixSize; j++)
                    line += settings.Chi[i][j]+" ";
                outputWriter.WriteLine(line);
            }

        }
    }
    public enum CalculationMode
    {
        Usual,
        InfinitlyDelute,
    }
    public struct Vector3
    {
        public double x1;
        public double x2;
        public double x3;
        public Vector3(double x1,double x2, double x3)
        {
            this.x1 = x1;
            this.x2 = x2;
            this.x3 = x3;
        }
        public override string ToString()
        {
            return x1 + ";" + x2 + ";" + x3 + ";";
        }
    }
    public class ProfileInfo
    {
        public List<double> composition = new List<double>();
        public List<double> additiveEquationError = new List<double>();
        public double polymerEquationError;
        public double polymerEquationStretchingPart;
        public double polymerEquationMixingPart;
        public double mixingEnergyContributionToF;
        public double entropyContributionToF;
    }
}
