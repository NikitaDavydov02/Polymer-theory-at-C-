using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Runtime.Remoting.Messaging;
using System.Reflection;
using Polymer_brush.Components;

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

        static double[] Xbrush, fipolimer;
        static double point_y;
        static StreamWriter sw;
        static StreamWriter newthonWriter;
        static StreamWriter integralLogWriter;
        static double z = 6;
        //static double[] chemPotInTheBulk;
        static double[] chemPotAtTheBorder;
        static MixingPartModule mixingPartModule;
        static CalculationMode calculationMode;
        static void Main(string[] args)
        {
            /*double[] guess = new double[2];
			double[] X;
			double norm;
			guess[0] = 2;
			guess[1] = 2;
			DNEQNF(TestSystem, 0.0001, 2, 100, guess, out X, out norm);
			Console.WriteLine("X[0] = " + X[0]);
			Console.WriteLine("X[1] = " + X[1]);
			Console.ReadLine();*/

            double Fmix = 0;
            double u_sol = 0;
            double u_pol = 0;
            double u_bio = 0;
            calculationMode = CalculationMode.InfinitlyDelute;
            Enter();
            CalculateMixingSurface(0.005, false);
            CalculateMixingSurface(0.005,true);
            if (calculationMode == CalculationMode.InfinitlyDelute)
                SwitchToInfinitlyDeluteMode();

            //////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////
            List<KeyValuePair<double, List<double>>> mixingEnergy = CalculateMixingEnergyProfile(1, 0, 20);
            using (StreamWriter sw = new StreamWriter("mixingFenergyOfSolventAndPolymer.txt"))
            {
                foreach (KeyValuePair<double, List<double>> pair in mixingEnergy)
                {
                    sw.WriteLine(pair.Key + "  ;  " + pair.Value[0] + ";" + pair.Value[1] + ";");
                }
            }
            //////////////////////////////////////////////////////////////////
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

            List<KeyValuePair<double, List<double>>> profile = new List<KeyValuePair<double, List<double>>>();
            Console.WriteLine("Calculating profile...");
            while (y_cur < y_edge)
            {
                Xbrush = new double[NumberOfComponents];
                FindVolumeFractionsInTheBrushForPoint(out Xbrush, y_cur); //  !calculates concentration profile in the brush after the solution is found
                

                for (int i = 0; i < NumberOfComponents; i++)
                {
                    fipolimer[i] = Xbrush[i];
                }
                Console.Write("Fmix=" + Fmix);
                Console.Write("u_sol=" + u_sol);
                Console.Write("u_pol=" + u_pol);
                Console.Write("u_bio=" + u_bio);
                Console.Write("/n");

                /*string line = y_cur.ToString() + "    ";
                for (int i = 0; i < NumberOfComponents; i++)
                    line += fipolimer[i] + "    ";
                sw.WriteLine(line + (Osmmix(Program.NumberOfComponents, fipolimer) - osmbulk));*/
                List<double> composition = new List<double>();
                for (int i = 0; i < NumberOfComponents; i++)
                    composition.Add(fipolimer[i]);
                profile.Add(new KeyValuePair<double, List<double>>(y_cur, composition));
                y_cur += aA / R;
            }
            y_cur = y_min;
            if (calculationMode == CalculationMode.InfinitlyDelute)
            {
                Enter();
                Console.WriteLine("*************************");
                Console.WriteLine("*************************");
                Console.WriteLine("Finding additive concentration in infinitly delute solution approximation...");
                for (int i = 0; i < profile.Count; i++)
                {
                    y_cur = profile[i].Key;
                    double[] baseFractions = new double[2];
                    baseFractions[0] = profile[i].Value[0];
                    baseFractions[1] = profile[i].Value[1];
                    double additiveFraction = FindAdditiveConcentrationForParticularPolymerAndSolventContentInTheBrush(baseFractions);
                    profile[i].Value.Add(additiveFraction);
                }
            }
            Console.WriteLine("");
            Console.WriteLine("Output");
            for(int i=0;i<profile.Count;i++)
            {
                string line = profile[i].Key.ToString() + "    ";
                for (int j = 0; j < profile[i].Value.Count; j++)
                    line += profile[i].Value[j] + "    ";
                sw.WriteLine(line + (Osmmix(Program.NumberOfComponents, fipolimer) - osmbulk)); 
            }

            Console.WriteLine("Beta: " + y_edge);
            Console.WriteLine("Calculation is done!");


            sw.Close();
            Console.ReadLine();
        }
        static void SaveInputSettings()
        {
            Input input = new Input();
            input.AreaPerChain = areaPerChain;
            input.Components = new List<Component>();
            Component solvent = new Component();

        }
        static void Enter()
        {
            Console.WriteLine("Initializing...");
            c = 2.0;
            pi = 3.1415926;
            coe = (3.0 / 8.0) * pi * pi;
            aA = 6.8 * Math.Pow(10, -9);
            rN = 80.0; //total number of polymer segments per chain
            rNA = 60.0; //number of segments in A - subchain
            rNB = rN - rNA;
            nu = 2.0; //spherical micelle
            areaPerChain = (7.0 * Math.Pow(10, -9)) * (7.0 * Math.Pow(10, -9)) / 0.12;
            areaPerChain = (7.0 * Math.Pow(10, -9)) * (7.0 * Math.Pow(10, -9)) / 0.6;
            R = rNB * (nu + 1) * aA * aA * aA / areaPerChain; // core radius((nu+1.0)*rNB / areaPerChain)**2    LAGRANGE_STR1 = BA * (  - 1.0) * *2 * ((nu + 1.0) * rNB / areaPerChain) * *2 * aA * *3  CHECK THIS
                                                              //write(*, *) R,areaPerChain
                                                              //
            BA = coe / ((rNA * aA) * (rNA * aA));

            //y_min = 1.0 + 0.1 * aA / R;
            y_min = 1.0 + aA / R;
            y_max = 1.00001 * (1.0 + 1.0 * aA * rNA / R);                 //y_max = 1.0d0 * (1.0d0 + 2.0 * aA * rNA / R)
            yacc = Math.Pow(10, -8);

            //! areaPerChain min for all morphologies:
            //!areaPerChain_MIN = 10.0 * (nu + 1.1) * aA * *2
            //!areaPerChain_MAX = 3.0d0 * aA * *2 * (3.1415926 * 4.0 * rNB * *2 / 3.0d0) **(1.0 / 3.0)
            //				 !write(*, *) BA,R,rNB* aA, areaPerChain,10.0 * (nu + 1.1) * aA * *2,y_max,BA * (R * (y_max - 1.0)) * *2,
            //!stop
            NumberOfComponents = 3;
            NumberOfPolymerGroupTypes = 1;

            chemPotInTheBulk = new double[NumberOfComponents];
            chemPotAtTheBorder = new double[NumberOfComponents];
            fipolimer = new double[NumberOfComponents];

            size = new double[MaxNumberOfComponent];
            for (int i = 0; i < NumberOfComponents; i++)
                size[i] = 1.0;

            size[0] = 1.0;// ! solvent
            size[1] = rNA;// polymer
            size[2] = 3.0;// bio
            //size[2] = 60.0;

            chi = new double[MaxNumberOfComponent, MaxNumberOfComponent];
            chiMatrixSize = NumberOfComponents + NumberOfPolymerGroupTypes - 1;
            //solv
            //bio
            //pol
            fractionsOfGroups = new double[NumberOfPolymerGroupTypes];
            for (int i = 0; i < NumberOfPolymerGroupTypes; i++)
                fractionsOfGroups[i] = 0.1;
            fractionsOfGroups[0] = 1;
            //fractionsOfGroups[1] = 0;

            chi[0, 1] = 0.5;//solv-pol
            chi[0, 2] = 0;//solv-bio
            chi[1, 2] = 0.3;//pol-bio
            /*//Solvent with other
            chi[0, 3] = -1;//! solv - bio
            chi[0, 1] = 0;//! solv - polym first group
            chi[0, 2] = 0;//! solv - polym second group
            //Bio with other
            chi[3, 1] = -1; //bio- polym first group
            chi[3, 2] = -1;  //bio- polym second group
            //Polymer A with other
            chi[1, 2] = 0;*/

            for (int i = 0; i < chiMatrixSize; i++)
            {
                chi[i, i] = 0;
                for (int j = i + 1; j < chiMatrixSize; j++)
                    chi[j, i] = chi[i, j];
            }

            etas = new double[chiMatrixSize, chiMatrixSize];
            for (int i = 0; i < chiMatrixSize; i++)
                for (int j = 0; j < chiMatrixSize; j++)
                    etas[i, j] = Math.Exp(-3 * chi[i, j] / z);

            //!give bulk composition and calculate Lagr.multipliers and osm pressure in the bulk:
            ////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////
            
            volumeFractionsInTheBulk = new double[NumberOfComponents];
            for (int i = 0; i < NumberOfComponents; i++)
                volumeFractionsInTheBulk[i] = 0.0;
            volumeFractionsInTheBulk[2] = 0.02;//bio
            volumeFractionsInTheBulk[0] = 0.98;//solvent
            //volumeFractionsInTheBulk[0] = 1.0;

            ////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////

            
            //Lagrmix(2, volumeFractionsInTheBulk, out Lagrbulk);

            mixingPartModule = new MixingPartModule();

            Lagrmix(NumberOfComponents, volumeFractionsInTheBulk, out chemPotInTheBulk);
            //chemPotInTheBulk[1] = Lagrbulk[1];//this is for biocomponent
            //osmbulk = Osmmix(NumberOfComponents - 1, volumeFractionsInTheBulk);// !  this is for solvent


        }
        static void SwitchToInfinitlyDeluteMode()
        {
            Console.WriteLine("Switching to infinitly delute mode");
            NumberOfComponents = 2;
            volumeFractionsInTheBulk = new double[2];
            volumeFractionsInTheBulk[1] = 0;
            volumeFractionsInTheBulk[0] = 1;
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
            double integrationMin = 1.0;
            double integrationMax = y;
            double s = 1;


            //Find chemical potentials at border;
            for (int i = 0; i < NumberOfComponents; i++)
                chemPotAtTheBorder[i] = chemPotInTheBulk[i];
            //chemPotAtTheBorder[0] = chemPotInTheBulk[0];//solvent
            //chemPotAtTheBorder[1] = chemPotInTheBulk[1];//bio
            //Finding volume fractions at the border
            double ERREL = Math.Pow(10, -4);
            //double ERREL = Math.Pow(10, -4);

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
            //XBorderGUESS[0] = 0.1;

            double FNORM;
            double[] _XBorder = new double[NumberOfComponents-1];
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
                        return;
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
            //Lamb_Pol = chemPotInTheBulk[2] + BA * (R * (  - 1)) * (R * (  - 1));
            
            FindVolumeFractionsInTheBrushForPoint(out Xbrush, x);
            func_log = Xbrush[1].ToString()+"_;";
            
            double fay = Xbrush[1];
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
            if (Math.Abs(y_cur - 1.41210111495579)<0.001)
                ;
            int ITMAX = 600;
            double[] XBrushGUESS = new double[NumberOfComponents - 1];//Everything except solvent
            double[] XBrushReduced = new double[NumberOfComponents - 1];//Everything except solvent
            XBrush = new double[NumberOfComponents];
            
            XBrushGUESS[0] = 0.97;//this is the fraction of polymer in the brush
            for(int i=1;i<XBrushGUESS.Length;i++)
                XBrushGUESS[i] = Math.Pow(10, -8);//this is the fraction of biocomponent in the brush

            //XBrushGUESS[0] = 9*Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
            //XBrushGUESS[1] = 0;//this is the fraction of polymer in the brush
            double FNORM;
            if (Math.Abs(point_y - 1.41210111495579) < 0.001)
                ;
            DNEQNF(BrushEquations, ERREL, NumberOfComponents - 1, ITMAX, XBrushGUESS, out XBrushReduced, out FNORM);

            double volumeFractionsSum = 0;
            for (int i = 1; i < NumberOfComponents; i++)
            {
                XBrush[i] = XBrushReduced[i-1];
                volumeFractionsSum += XBrush[i];
            }
            XBrush[0] = 1 - volumeFractionsSum;
            //TryingToFindStepInTheBrush(out XInside, y_cur, XBrush);
        }
        delegate string NonlinearSystem(double[] X, out double[] F, int L);

        /* static double[,] CalculateJacobian(NonlinearSystem Func,double[] _X)
         {
             double[] X = new double[_X.Length];
             for (int i = 0; i < X.Length; i++)
                 X[i] = _X[i];
             double[] F = new double[X.Length];
             double[,] J = new double[X.Length, X.Length];

             for (int i = 0; i < X.Length; i++)
                 for (int j = 0; j < X.Length; j++)
                 {
                     Func(X, out F, X.Length);
                     double f_init = F[i];
                     double old_x = X[j];
                     double dx = old_x * 0.01;
                     X[j] += dx;


                     for (double devisionStepDegree = 1; X[j] >= 1; devisionStepDegree++)
                     {
                         dx /= 2;
                         X[j] = dx + old_x;
                     }


                     Func(X, out F, X.Length);
                     double f_df = F[i];
                     J[i, j] = (f_df - f_init) / dx;
                     X[j] = old_x;
                     Func(X, out F, X.Length);
                 }
             return J;
         }*/
        static void DNEQNF(NonlinearSystem Func, double ERREL, int L, int ITMAX, double[] XGuess, out double[] X, out double FNORM)
        {
            int splitTransitions = 0;
            newthonWriter = new StreamWriter(File.Create("newthon_log.txt"));
            // newthonWriter.WriteLine("Iteration; X0; X1; F0; F1; FNORM; J00; J01; J10; J11;");
            newthonWriter.WriteLine("Iteration; X0; F0; FNORM; J00;");
            int iterations = 0;
            double[] F = new double[L];
           
            X = new double[XGuess.Length];
            double[] oldX = new double[XGuess.Length];
            double[] deltaX = new double[XGuess.Length];
            for (int i = 0; i < XGuess.Length; i++)
                X[i] = XGuess[i];
            double[] realComposition = new double[NumberOfComponents];
            

            FNORM = 0;
            Func(X, out F, L);
            for (int i = 0; i < L; i++)
                FNORM += F[i] * F[i];

            double[,] J = new double[L, XGuess.Length];
            while (FNORM >= ERREL)
            {
                //Calculate Jacobian
                for (int i = 0; i < L; i++)
                    for (int j = 0; j < X.Length; j++)
                    {
                        Func(X, out F, L);
                        double f_init = F[i];
                        double old_x = X[j];
                        double dx = old_x * 0.01;
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
                double eps = 0.000000000001;
                for (int j = 0; j < X.Length; j++)
                    oldX[j] = X[j];
                if (det == 0)
                {
                    //Node transition
                    integralLogWriter.Close();
                    newthonWriter.Close();
                    throw new Exception();
                    ;
                }


                double[,] reverse = Matrix.reverseMatrix(L, J);
                double[] rightParts = new double[L];
                for (int i = 0; i < L; i++)
                {
                    rightParts[i] = -F[i];
                    for (int j = 0; j < X.Length; j++)
                        rightParts[i] += J[i, j] * X[j];
                }
                //for (int j = 0; j < X.Length; j++)
                //   oldX[j] = X[j];

                X = Matrix.multiplyMatrixAndVector(L, L, reverse, rightParts);
                

                for (int i = 0; i < X.Length; i++)
                    deltaX[i] = X[i] - oldX[i];

                if (double.IsNaN(X[0]))
                    ;

                for (double devisionStepDegree = 1; ContainsOutrangeValues(X); devisionStepDegree++)
                {
                    for (int j = 0; j < X.Length; j++)
                    {
                        deltaX[j] /= 2;
                        X[j] = oldX[j] + deltaX[j];
                    }
                }
                //<Split>
                //Check if in split zone
                double sum = 0;
                for (int i = 0; i < XGuess.Length; i++)
                {
                    realComposition[i + 1] = X[i];
                    sum += X[i];
                }
                realComposition[0] = 1 - sum;
                double segregationDelta;
                double[] firstSegregationPoint;
                double[] secondSegregationPoint;
                if (mixingPartModule.IsCompositionInsideSegregationZone(realComposition, out segregationDelta, out firstSegregationPoint, out secondSegregationPoint))
                {
                    if (NumberOfComponents != 2)
                        throw new NotImplementedException();
                    //Split
                    if (deltaX[0] > 0)
                        X[0] = secondSegregationPoint[1] + 0.01;
                    if (deltaX[0] < 0)
                        X[0] = firstSegregationPoint[1] - 0.1;
                    splitTransitions++;
                }
                /*if(X[0]>=mixingPartModule.segregationPoints[0]&& X[0]<= mixingPartModule.segregationPoints[1])
                {
                    if (NumberOfComponents != 2)
                        throw new Exception();
                    //Split
                    if (deltaX[0] > 0)
                        X[0] = mixingPartModule.segregationPoints[1] + 0.01;
                    if (deltaX[0] < 0)
                        X[0] = mixingPartModule.segregationPoints[0] - 0.1;
                    splitTransitions++;
                }*/
                //</Split>
                iterations++;
                if (iterations > ITMAX)
                {
                    if (splitTransitions > ITMAX * 0.2)
                    {
                        //Split
                        newthonWriter.WriteLine("Split");
                        newthonWriter.Close();
                        if (NumberOfComponents != 2)
                            throw new NotImplementedException();
                        X[0] = firstSegregationPoint[1];

                        return;
                    }
                    else
                    {
                        newthonWriter.Close();
                        integralLogWriter.Close();
                        throw new Exception(" Newton method did not manage to find solution for system of equations");

                    }
                }
                string funcLog = Func(X, out F, L);
                FNORM = 0;
                for (int i = 0; i < L; i++)
                    FNORM += F[i] * F[i];
                ///////////////////////////////////////////////////////////////////////////
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
                
            }
            newthonWriter.WriteLine("Converged!");
            newthonWriter.Close();
            
        }
        static bool ContainsOutrangeValues(double[] X)
        {
            foreach (double x in X)
                if (x < 0 || x > 1)
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
                surfaceWriter.WriteLine("Solvent;Polymer;Bio;F;");
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
                        double a;
                        bool segregated = mixingPartModule.IsCompositionInsideSegregationZone(fractions, out a);
                        if(linearized)
                            surfaceWriter.WriteLine(x1 + ";" + x2 + ";" + x3 + ";" + F + ";"+ segregated+";");
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
                double Fmix = mixingPartModule.CalculateMixingFreeEnergy(volumeFractions);
                double[] exchangeChemPotentials;
                Lagrmix(NumberOfComponents, volumeFractions, out exchangeChemPotentials);
                value.Add(Fmix);
                value.Add(exchangeChemPotentials[1]);
                output.Add(new KeyValuePair<double, List<double>>(volumeFractions[AcomponentIndex], value));
                volumeFractions[AcomponentIndex] += step;
                volumeFractions[BcomponentIndex] = 1 - volumeFractions[AcomponentIndex];
            }

            return output;
        }
       // static double[] chemPotAtTheBaseBrush;
        static double[] compositionAtTheBaseBrush;
        static double FindAdditiveConcentrationForParticularPolymerAndSolventContentInTheBrush(double[] SolventAndPolymerFractions)
        {
            Console.WriteLine("Finding Additive Concentration...");
            compositionAtTheBaseBrush = new double[3];
            compositionAtTheBaseBrush[0] = SolventAndPolymerFractions[0];
            compositionAtTheBaseBrush[1] = SolventAndPolymerFractions[1];
            compositionAtTheBaseBrush[2] =0;
            
            //Lagrmix_PolA(NumberOfComponents, compositionAtTheBaseBrush, out chemPotAtTheBaseBrush);

            double ERREL = Math.Pow(10, -2) * volumeFractionsInTheBulk[2];
            int ITMAX = 600;
            double[] tryAdditive = new double[] { volumeFractionsInTheBulk[2] };
            double[] additiveConcentrations;
            double FNORM;
            DNEQNF(InfinitlyDeluteAdditivesEquations, ERREL, 1, ITMAX, tryAdditive, out additiveConcentrations, out FNORM);

            return additiveConcentrations[0];
        }
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
       /* static double[] FindSegregationPointsBetweenSolventAndPolymer()
        {
            double[] initialComposition = new double[3];

            initialComposition[2] = 0.7; //polymer
            initialComposition[1] = 0;//bio
            initialComposition[0] = 1 - initialComposition[2]; //solvent
            double Finit = mixingPartModule.CalculateMixingFreeEnergy(initialComposition);

            double x1 = initialComposition[2];//polymer molar fraction
            double x2 = initialComposition[2];
            double dx = 0.01;

            double[] X = new double[3];
            for (int i = 0; i < 3; i++)
                X[i] = 0;

            double bestFmix = 10000000;
            double[] best_x = new double[2];

            while (x1 > 0)
            {
                x2 = initialComposition[2] + dx;

                double[] leftComposition = new double[3];
                leftComposition[0] = 1 - x1;
                leftComposition[1] = 0;
                leftComposition[2] = x1;
                double leftF = mixingPartModule.CalculateMixingFreeEnergy(leftComposition);

                while (x2 < 1)
                {
                    double[] rightComposition = new double[3];
                    rightComposition[0] = 1 - x2;
                    rightComposition[1] = 0;
                    rightComposition[2] = x2;
                    double rightF = mixingPartModule.CalculateMixingFreeEnergy(rightComposition);
                    double Fsep = leftF + (initialComposition[2] - x1) * (rightF - leftF) / (x2 - x1);
                    double delta = Fsep - Finit;
                    if (delta < bestFmix)
                    {
                        bestFmix = delta;
                        best_x[0] = x1;
                        best_x[1] = x2;
                    }
                    x2 += dx;
                }
                x1 -= dx;
            }
            return best_x;
        }
        */
        //static void SegregationEquations(double[] X, out double[] F, int L)
        //{
        //    double nu = 2;

        //    double[] mixingPartOfExchangeChemicalPotentials_0;
        //    double[] volumeFractions_0 = new double[3];
        //    volumeFractions_0[0] = X[0];
        //    volumeFractions_0[1] = 0;
        //    volumeFractions_0[2] = 1 - X[0];
        //    //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
        //    Lagrmix(3, volumeFractions_0, out mixingPartOfExchangeChemicalPotentials_0);

        //    double[] mixingPartOfExchangeChemicalPotentials_1;
        //    double[] volumeFractions_1 = new double[3];
        //    volumeFractions_1[0] = X[1];
        //    volumeFractions_1[1] = 0;
        //    volumeFractions_1[2] = 1 - X[1];
        //    //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
        //    Lagrmix(3, volumeFractions_1, out mixingPartOfExchangeChemicalPotentials_1);
        //    F = new double[L];

        //    //F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
        //    F[0] = (mixingPartOfExchangeChemicalPotentials_0[1] - mixingPartOfExchangeChemicalPotentials_1[1]) * (mixingPartOfExchangeChemicalPotentials_0[1] - mixingPartOfExchangeChemicalPotentials_1[1]);
        //    double secondEquation = mixingPartModule.CalculateMixingFreeEnergy(volumeFractions_1) - mixingPartModule.CalculateMixingFreeEnergy(volumeFractions_0) - (volumeFractions_1[1] - volumeFractions_0[1]) * mixingPartOfExchangeChemicalPotentials_0[1];
        //    F[1] = Math.Pow(secondEquation, 2);
        //}

        //static double CalculateMixingFreeEnergy(double[] X)
        //      {
        //	double a = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) /  size[1];
        //	double b = chi[1, 2] * X[1] * X[2];
        //	double c = chi[0, 1] * X[0] * X[1];
        //	double d = chi[0, 2] * X[0] * X[2];
        //	return a + b + c + d;
        //}
        //static double CalculateGugenheimMixingFreeEnergy(double[] X)
        //{
        //	int n = X.Length;
        //	double translationSum = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) /  size[1];
        //	double[] XX = CalculateGugenheimCorrelations(X, etas);
        //	double mixingSum = 0;
        //	for (int i = 0; i < n; i++)
        //		for (int j = 0; j < i; j++)
        //			mixingSum += chi[i, j] * XX[i] * XX[j] * X[i] * X[j] * etas[i, j];

        //	/*double b = chi[1, 2] * X[1] * X[2];
        //	double c = chi[0, 1] * X[0] * X[1];
        //	double d = chi[0, 2] * X[0] * X[2];*/
        //	double output = translationSum + mixingSum;
        //	return output;
        //}
        //static double CalculateExchangeChemialPotentialOfComponent(double[] X, int componenIndex)
        //      {
        //	double f = CalculateMixingFreeEnergy(X);
        //	//double f = CalculateGugenheimMixingFreeEnergy(X);
        //	double x = X[componenIndex];
        //	double max_dx = 1 - x;
        //          if (X[0] < max_dx)
        //		max_dx = X[0];
        //	double dx = 0.01*x;
        //	if (dx == 0)
        //		dx = 0.01;
        //	if (dx > max_dx)
        //		dx = max_dx;

        //	double oldSolventVolumeFraction = X[0];
        //	double x_dx = x + dx;
        //	X[componenIndex] = x_dx;
        //	X[0] -= dx;
        //	double f_df = CalculateMixingFreeEnergy(X);
        //	//double f_df = CalculateGugenheimMixingFreeEnergy(X);
        //	X[0] = oldSolventVolumeFraction;
        //	X[componenIndex] = x;
        //	return (f_df - f) / dx;


        //}
        //static double[] CalculateGugenheimCorrelations(double[] alphas, double[,] etas)
        //      {
        //	int n = alphas.Length;
        //	double[] XX = new double[n];
        //	double[] newXX = new double[n];
        //	double initialGuess = 1;
        //	for(int i = 0; i < n; i++)
        //          {
        //		XX[i] = initialGuess;
        //		initialGuess -= 0.0001;
        //	}
        //	bool converged = false;
        //          while (!converged)
        //          {
        //		for (int i = 0; i < n; i++)
        //		{
        //			double sum = 0;
        //			for (int j = 0; j < n; j++)
        //				sum += alphas[j] * XX[j] * etas[i, j];
        //			newXX[i] = 1 / sum;
        //		}
        //		for (int i = 0; i < n; i++)
        //			XX[i] = (XX[i] + newXX[i]) / 2;
        //		converged = true;
        //		for (int i = 0; i < n; i++)
        //			if (Math.Abs(XX[i] - newXX[i]) > Math.Pow(10, -12))
        //			{
        //				converged = false;
        //				break;
        //			}
        //	}
        //	return XX;

        //}
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

}
