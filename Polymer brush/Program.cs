using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Runtime.Remoting.Messaging;

namespace Polymer_brush
{
	class Program
	{
		public static double R, BA, y_min, y_max, yacc, aA, rNA, rNB, areaPerChain, Lamb_Pol;
		public static double rN, nu, coe, pi, osmbulk, y_edge, y_cur, P_AB;
		public static double[] chemPotInTheBulk, volumeFractionsInTheBulk, Lagrbulk,  size;
		public static double[] chemPotOutsideOfTheStep;
		public static double osmoticPressureOutsideOfTheStep;
		public static double[] segregationPoints;

        public static double[,] chi;
		public static double[,] etas;
		public static double[] fractionsOfGroups;
		public static int NumberOfComponents;
		public static int NumberOfPolymerGroupTypes;
		public static int chiMatrixSize;
		
		static double[] Xbrush, fipolimer;
		static double point_y;
		static StreamWriter sw;
		static StreamWriter logWriter;
		static StreamWriter newthonWriter;
		static double z = 6;
		//static double[] chemPotInTheBulk;
		static double[] chemPotAtTheBorder;
		static MixingPartModule mixingPartModule;

		static double[] XInBrushForPreviousPoint;
		static void Main(string[] args)
		{
			mixingPartModule = new MixingPartModule();
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
			double u_bio=0;
			Enter();
			List<KeyValuePair<double, List<double>>> mixingEnergy = CalculateMixingEnergyProfile(2, 0, 20);
			logWriter = new StreamWriter("log.txt");
			//newthonWriter = new StreamWriter(File.Create("newthonLog.txt"));
            using (StreamWriter sw = new StreamWriter("mixingFenergyOfSolventAndPolymer.txt"))
			{
				foreach (KeyValuePair<double, List<double>> pair in mixingEnergy)
					sw.WriteLine(pair.Key + "  ;  " + pair.Value[0] + ";"+ pair.Value[1] + ";");
			}
            //segregationPoints = FindSegregationPointsBetweenSolventAndPolymer();
			segregationPoints = new double[2];
			segregationPoints[0] = 0.42;
            segregationPoints[1] = 0.95;
            //return;
            sw = new StreamWriter("profile.txt");
            sw.WriteLine("y_cur    solvent    bio    polymer    osm_pressure");

            y_edge = BisectionSolve(y_min, y_max, yacc,NormalizationFunction,new List<double>());
			y_cur = y_min;

			//!Calculate polymer concentration at the start position(A / B boundary)
			//P_AB = NormalizationSubintegralValue(y_cur, new List<double>() { y_edge });
			//P_AB = FindSubintegralValueForNormalization(y_cur, y_edge); // !   Polymer concentration at A/ B boundary(highest)  must work before finding profile to set Lagr multipliers in common block
           
			
			while (y_cur < y_edge)
            {
				//y_cur += aA / R;
				//!write(*, *) 'Phi poly subroutine is called'
				Xbrush = new double[NumberOfComponents];
				//Xbrush - everithing exept solvent
				FindVolumeFractionsInTheBrushForPoint(out Xbrush, y_cur); //  !calculates concentration profile in the brush after the solution is found

				
				for(int i = 0; i < NumberOfComponents; i++)
                {
					fipolimer[i] = Xbrush[i];
				}
				/*fipolimer[1] = Xbrush[0];//biocomponent // ! local compsition of the brush at point y_cur(1)-solvent(2) - biocomponent(3) - polymer
				fipolimer[2] = Xbrush[1];//polymer
				fipolimer[0] = 1.0 - fipolimer[1] - fipolimer[2];//solvent*/

				Console.Write("Fmix="+Fmix);
				Console.Write("u_sol=" +u_sol);
				Console.Write("u_pol=" + u_pol); 
				Console.Write("u_bio=" + u_bio);
				Console.Write("/n");

				string line = y_cur.ToString() + "    ";
				for (int i = 0; i < NumberOfComponents; i++)
					line += fipolimer[i] + "    ";
				sw.WriteLine(line + (Osmmix(3, fipolimer) - osmbulk));
				y_cur += aA / R;
			}

			Console.WriteLine("Beta: " + y_edge);
			Console.WriteLine("Calculation is done!");

			logWriter.Close();
			sw.Close();
			//newthonWriter.Close();
			Console.ReadLine();
		}
		static void Enter()
        {
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
			//y_min = 1.4;

            y_max = 1.00001 * (1.0 + 1.0 * aA * rNA / R);                 //y_max = 1.0d0 * (1.0d0 + 2.0 * aA * rNA / R)
			yacc = Math.Pow(10, -8);

			//! areaPerChain min for all morphologies:
			//!areaPerChain_MIN = 10.0 * (nu + 1.1) * aA * *2
			//!areaPerChain_MAX = 3.0d0 * aA * *2 * (3.1415926 * 4.0 * rNB * *2 / 3.0d0) **(1.0 / 3.0)
			//				 !write(*, *) BA,R,rNB* aA, areaPerChain,10.0 * (nu + 1.1) * aA * *2,y_max,BA * (R * (y_max - 1.0)) * *2,
			//!stop
			NumberOfComponents = 3;
			NumberOfPolymerGroupTypes = 2;

            XInBrushForPreviousPoint= new double[NumberOfComponents-1];
			for (int i = 0; i < NumberOfComponents - 1; i++)
				XInBrushForPreviousPoint[i] = -0.5;

                size = new double[NumberOfComponents];
			for (int i = 0; i < NumberOfComponents; i++)
				size[i] = 1.0;

			 size[0] = 1.0;// ! solvent
			 size[1] = 3.0;// ! bioadditive
             size[2] = rNA;// polymer
            //size[2] = 60.0;

            chi = new double[NumberOfComponents+NumberOfPolymerGroupTypes-1, NumberOfComponents + NumberOfPolymerGroupTypes - 1];
			chiMatrixSize = NumberOfComponents + NumberOfPolymerGroupTypes - 1;
			//solv
			//bio
			//pol
			fractionsOfGroups = new double[NumberOfPolymerGroupTypes];
			for (int i = 0; i < NumberOfPolymerGroupTypes; i++)
				fractionsOfGroups[i] = 0.1;
			fractionsOfGroups[0] = 1;
			fractionsOfGroups[1] = 0;


			//Solvent with other
			chi[0, 1] = 1;//! solv - bio
			chi[0, 2] =0.5;//! solv - polym first group
			chi[0, 3] = 0.5;//! solv - polym second group

			//Bio with other
			chi[1, 2] = 1; //bio- polym first group
			chi[1, 3] = 1;  //bio- polym second group

			//Polymer A with other
			chi[2, 3] = 0;

			for (int i = 0; i < chiMatrixSize; i++)
			{
				chi[i, i] = 0;
				for (int j = i + 1; j < chiMatrixSize; j++)
					chi[j, i] = chi[i, j];
			}

			etas = new double[chiMatrixSize, chiMatrixSize];
			for (int i = 0; i < chiMatrixSize; i++)
				for (int j = 0; j < chiMatrixSize; j++)
					etas[i, j] = Math.Exp(-3*chi[i, j]/z);

			//!give bulk composition and calculate Lagr.multipliers and osm pressure in the bulk:
			volumeFractionsInTheBulk = new double[NumberOfComponents];
			for(int i=0;i<NumberOfComponents;i++)
				volumeFractionsInTheBulk[i] = 0.0;
			volumeFractionsInTheBulk[0] = 1.0;//solvent
			volumeFractionsInTheBulk[1] = 1.0 - volumeFractionsInTheBulk[0];//bio

			chemPotInTheBulk = new double[NumberOfComponents];
			chemPotAtTheBorder = new double[NumberOfComponents];
			fipolimer = new double[NumberOfComponents];
			//Lagrmix(2, volumeFractionsInTheBulk, out Lagrbulk);
			Lagrmix(NumberOfComponents-1, volumeFractionsInTheBulk, out chemPotInTheBulk);
			//chemPotInTheBulk[1] = Lagrbulk[1];//this is for biocomponent
			osmbulk = Osmmix(NumberOfComponents-1, volumeFractionsInTheBulk);// !  this is for solvent

			//chemPotInTheBulk[0] = Lagrbulk[0]; //!  this should be identiacally zero
			//! now calculate  mixing part of Lagrange multiplier at the edge of the brush, where  FiA = 0:
			//volumeFractionsInTheBulk[2] = 0.0;
			//Lagrmix_PolA(3, volumeFractionsInTheBulk, out Lagrbulk);// ! NOTE  that  Lagr multipliers are spoiled in  Lagrbulk   but stored in common
			//chemPotInTheBulk[2] = Lagrbulk[2];

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
				AlternativemixingPartOfExchangeChemicalPotentials[i] = mixingPartModule.CalculateExchangeChemialPotentialOfComponent(X, i);
				mixingPartOfExchangeChemicalPotentials[i] = AlternativemixingPartOfExchangeChemicalPotentials[i];                                                                                     //AlternativemixingPartOfExchangeChemicalPotentials(i) = CalculateExchangeChemialPotentialOfComponent(3, X, i)
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
				sum1 = sum1 + X[i] * (1.0 - 1.0 /  size[i]);
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
			for(int i=0;i<NumberOfComponents-1;i++)
				chemPotAtTheBorder[i] = chemPotInTheBulk[i];
			//chemPotAtTheBorder[0] = chemPotInTheBulk[0];//solvent
			//chemPotAtTheBorder[1] = chemPotInTheBulk[1];//bio
			//Finding volume fractions at the border
			double ERREL = Math.Pow(10, -4);
			double[] XBorderGUESS = new double[NumberOfComponents-1];
			for (int i = 0; i < NumberOfComponents - 1; i++)
				XBorderGUESS[i] = volumeFractionsInTheBulk[i];//this is the fraction of solvent at the border
			//XBorderGUESS[1]= volumeFractionsInTheBulk[1];//this is the fraction of biocomponnt at the border
			//XBorderGUESS[0] = 0.1;

			double FNORM;
			double[] _XBorder = new double[NumberOfComponents];
			DNEQNF(BorderEquations, ERREL, NumberOfComponents-1, 1000, XBorderGUESS,out _XBorder, out FNORM);

			double[] volFractionsAtTheBorder = new double[NumberOfComponents];

			double volumeFractionSum = 0;
			for(int i = 0; i < NumberOfComponents - 1; i++)
            {
				volFractionsAtTheBorder[i] = _XBorder[i];
				volumeFractionSum += volFractionsAtTheBorder[i];
			}
			volFractionsAtTheBorder[NumberOfComponents-1] = 1 - volumeFractionSum;
			Lagrmix_PolA(NumberOfComponents, volFractionsAtTheBorder, out chemPotAtTheBorder);

			chemPotAtTheBorder[NumberOfComponents-1] += BA * (R * (integrationMax - 1)) * (R * (integrationMax - 1));
			Lamb_Pol = chemPotAtTheBorder[NumberOfComponents-1];

			logWriter.WriteLine();
			logWriter.WriteLine("Normalization function at y=" + y);
			CalculateIntegral(integrationMin, integrationMax,NormalizationSubintegralValue,new List<double>() { y ,-1,-1}, out s);
			return s -norm;

		}
		static void CalculateIntegral(double integrationMin, double integrationMax, Func<double,List<double>,double> func,List<double> parameters, out double s)
        {
			double EPS = 0.1;
			int JMAX = 10;
			double old_s = -1 * Math.Pow(10, -30);
			s = 0;
			logWriter.WriteLine("n;s");

            for (int n = 1; n < JMAX;n++)
            {
				s = CalculateIntegralWithDefeniteNumberOfTrapezoids(integrationMin, integrationMax, s, n, func, parameters);
                //s = CalculateNormalizationIntegralWithDefeniteNumberOfTrapezoids(integrationMin, integrationMax, s, n+1);
                logWriter.WriteLine(n + " ; " + s + ";");
                if (n > 4)
					if (Math.Abs(s - old_s) < EPS * Math.Abs(old_s) || (s == 0 && old_s == 0))
					{
                        logWriter.WriteLine("---");
                        return;
                    }
				old_s = s;
            }
            logWriter.Close();
            sw.Close();
            //newthonWriter.Close();
            Console.WriteLine("Too many steps in q trap");
			throw new Exception();


		}
		static double CalculateIntegralWithDefeniteNumberOfTrapezoids(double integrationMin,double integrationMax,double s, int n, Func<double,List<double>,double> func, List<double>parameters)
        {
			string indent = "    ";
			logWriter.WriteLine(indent + "Integrate n=" + n);
			//s = 0;
			if (n == 1)
			{
                s = 0.5 * (integrationMax - integrationMin) * (func(integrationMin, parameters) + func(integrationMax, parameters));
				logWriter.WriteLine(indent + s);
            }
            else
            {
				int numberOfSegments = (int)Math.Pow(2, n - 2);
				//double tnm = numberOfSegments;
				double delta = (integrationMax - integrationMin) / numberOfSegments;
				double x = integrationMin + 0.5 * delta;
				double sum = 0;
				double lastGoodAdd = 0;
				for(int counter = 0; counter < numberOfSegments; counter++)
                {
                    //double add = FindSubintegralValueForNormalization(x, integrationMax);
                    ///////////COMMENT THIS/////////////
                    if (counter == 1)
					{
						parameters[1] = XInBrushForPreviousPoint[0];
						parameters[2] = XInBrushForPreviousPoint[1];
					}
                    //////////////////////////////
                    double add = func(x,parameters);
					logWriter.WriteLine(indent + indent + add);
					///////////COMMENT THIS/////////////
					if (double.IsNaN(add))
					{
                        add = lastGoodAdd;
                    }
					else
						lastGoodAdd = add;
					//////////////////////////////
					sum += add;
					x += delta;
                }
				double newS = (integrationMax - integrationMin) * sum / numberOfSegments;
				if (double.IsNaN(s))
					s = newS;
				else
					s = 0.5 * (s + newS);

			}
			return s;
		}
		
		static double NormalizationSubintegralValue(double x, List<double> parameters)
        {
			//Lamb_Pol = chemPotInTheBulk[2] + BA * (R * (  - 1)) * (R * (  - 1));
			if (parameters.Count > 1 && parameters[1]>=0)
			{
                double[] initialGuess = new double[2];
                initialGuess[0] = parameters[1];
                initialGuess[1] = parameters[2];
                FindVolumeFractionsInTheBrushForPoint(out Xbrush, x, initialGuess);
            }
            else
				FindVolumeFractionsInTheBrushForPoint(out Xbrush, x);
			double fay = Xbrush[NumberOfComponents-1];//polymer volume fraction
			return fay * x * x;

		}
		static double BisectionSolve(double y_min, double y_max, double yacc, Func<double,List<double>, double> func, List<double> parameters)
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
                logWriter.Close();
                sw.Close();
               // newthonWriter.Close();
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
            logWriter.Close();
            sw.Close();
            //newthonWriter.Close();
            throw new Exception();
		}
		static void FindVolumeFractionsInTheBrushForPoint(out double[] XBrush, double y_cur, double[] initialGuess = null)
        {
			
			//!Solving(localy) system of non - linear equations

			double ERREL = Math.Pow(10, -4);
			point_y = y_cur;
			int ITMAX = 600;
			double[] XBrushGUESS = new double[NumberOfComponents-1];
			double[] XBrushReduced = new double[NumberOfComponents-1];
			XBrush = new double[NumberOfComponents];
			if (initialGuess == null)
			{
                XBrushGUESS[0] = Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
                XBrushGUESS[1] = 0.97;//this is the fraction of polymer in the brush
            }
			else
			{
                XBrushGUESS[0] = initialGuess[0];//this is the fraction of biocomponent in the brush
                XBrushGUESS[1] = initialGuess[1];//this is the fraction of polymer in the brush
            }
			double FNORM;
			
			DNEQNF(BrushEquations, ERREL, NumberOfComponents - 1, ITMAX, XBrushGUESS, out XBrushReduced, out FNORM);

			double volumeFractionsSum = 0;
			for(int i = 1; i < NumberOfComponents; i++)
            {
				XBrush[i] = XBrushReduced[i - 1];
				volumeFractionsSum += XBrush[i];
                XInBrushForPreviousPoint[i - 1] = XBrushReduced[i - 1];
            }
			XBrush[0] = 1- volumeFractionsSum;
			double[] XInside;
			//TryingToFindStepInTheBrush(out XInside, y_cur, XBrush);
		}
		/*static void TryingToFindStepInTheBrush(out double[] XInside, double y_cur, double[] volFractionsOutside)
		{
            double ERREL = Math.Pow(10, -4);
            point_y = y_cur;
            int ITMAX = 600;
            double[] XBrushGUESS = new double[NumberOfComponents - 1];
            double[] XBrushReduced = new double[NumberOfComponents - 1];
            XInside = new double[NumberOfComponents];
            XBrushGUESS[0] = Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
            XBrushGUESS[1] = 0.98f;//this is the fraction of polymer in the brush
            double FNORM;

			chemPotOutsideOfTheStep = new double[NumberOfComponents];
            Lagrmix_PolA(NumberOfComponents, volFractionsOutside, out chemPotOutsideOfTheStep);
			osmoticPressureOutsideOfTheStep = CalculateOsmoticPressure(volFractionsOutside);

            DNEQNF(StepEquations, ERREL, NumberOfComponents - 1, ITMAX, XBrushGUESS, out XBrushReduced, out FNORM);

			
            double volumeFractionsSum = 0;
            for (int i = 1; i < NumberOfComponents; i++)
            {
                XInside[i] = XBrushReduced[i - 1];
                volumeFractionsSum += XInside[i];
            }
            XInside[0] = 1 - volumeFractionsSum;
			if (!double.IsNaN(XInside[0]))
				if (Math.Abs(XInside[0] - volFractionsOutside[0]) > 0.01)
					;
        }*/
		delegate void NonlinearSystem(double[] X, out double[] F, int L);
		static void DNEQNF(NonlinearSystem Func, double ERREL, int L, int ITMAX, double[] XGuess, out double[] X, out double FNORM)
		{
			newthonWriter = new StreamWriter(File.Create("newthonLog.txt"));
			newthonWriter.WriteLine("Function: "+Func.Method.Name+";");
			newthonWriter.WriteLine("Iteration; X[0];X[1];F[0];F[1];FNORM;J[00];J[01];J[10];J[11];");
			int iterations = 0;
			double[] F = new double[L];
			X = new double[XGuess.Length];
			double[] oldX = new double[XGuess.Length];
			double[] deltaX = new double[XGuess.Length];
			for (int i = 0; i < XGuess.Length; i++)
				X[i] = XGuess[i];
			FNORM = 0;
			Func(X, out F, L);
			for (int i = 0; i < L; i++)
				FNORM += F[i]*F[i];

			double[,] J = new double[L, XGuess.Length];
            while (FNORM >= ERREL)
            {
				//Calculate Jacobian
				for(int i=0;i<L;i++)
					for(int j=0;j< X.Length; j++)
                    {
                        if (y_cur == 1.706459054209919 && iterations == 22)
                            ;
                        Func(X, out F, L);
						double f_init = F[i];
						double old_x = X[j];
						double dx = old_x * 0.01;
						X[j] += dx;


                        /*double XSum = 0;
						for (int a = 0; a < L; a++)
							XSum += X[a];
;                       for (double devisionStepDegree = 1;XSum >= 1 && devisionStepDegree<10; devisionStepDegree++)
                        {
							dx /= 2;
							X[j] = dx + old_x;

                            XSum = 0;
                            for (int a = 0; a < L; a++)
                                XSum += X[a];
                        }*/
                        for (double devisionStepDegree = 1; X[j] >= 1; devisionStepDegree++)
                        {
                            dx /= 2;
                            X[j] = dx + old_x;
                        }

                        Func(X, out F, L);
						double f_df = F[i];
						J[i, j] = (f_df - f_init) / dx;
						X[j] = old_x;
						Func(X, out F, L);
					}
                for (int i = 0; i < L; i++)
                    FNORM += F[i] * F[i];
				newthonWriter.WriteLine(iterations + ";" + X[0] + ";" + X[1] + ";" + F[0] + ";" + F[1] + ";" + FNORM + ";" + J[0, 0] + ";" + J[0, 1] + ";" + J[1, 0] + ";" + J[1, 1] + ";");

                double det = Matrix.determinantGauss(L, J, 0, false);
				double[,] reverse = Matrix.reverseMatrix(L, J);
				double[] rightParts = new double[L];
				for(int i = 0; i < L; i++)
                {
					rightParts[i] = -F[i];
					for (int j = 0; j < X.Length; j++)
						rightParts[i] += J[i, j] * X[j];
                }

				for (int j = 0; j < X.Length; j++)
					oldX[j] = X[j];
				X = Matrix.multiplyMatrixAndVector(L, L, reverse, rightParts);

				for(int i=0;i<X.Length;i++)
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
				if (X[0] > 1 || X[1] > 1)
					;
				if (double.IsNaN(X[0]))
					;
				FNORM = 0;
				for (int i = 0; i < L; i++)
					FNORM += F[i] * F[i];
				iterations++;
                if (iterations > ITMAX)
                {
                    logWriter.Close();
                    sw.Close();
                    newthonWriter.Close();
                    throw new Exception(" Newton method did not manage to find solution for system of equations");
                }
			}
			newthonWriter.Close();
        }
		static bool ContainsOutrangeValues(double[] X)
        {
			foreach (double x in X)
				if (x < 0 || x > 1)
					return true;
			return false;
		}
        /*static void StepEquations(double[] X, out double[] F, int L)
        {
            double nu = 2;
            double y_cur = point_y;
            double volumeFractionSum = 0;
            for (int i = 0; i < NumberOfComponents - 1; i++)
            {
                fipolimer[i + 1] = X[i];
                volumeFractionSum += X[i];
            }
            fipolimer[0] = 1 - volumeFractionSum;
            double[] mixingPartOfExchangeChemicalPotentials;
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            Lagrmix_PolA(NumberOfComponents, fipolimer, out mixingPartOfExchangeChemicalPotentials);
            F = new double[L];
			double insideOsmoticPressure = CalculateOsmoticPressure(fipolimer);

            //F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
             for (int i = 0; i < NumberOfComponents - 2; i++)
             {
                 F[i] = (mixingPartOfExchangeChemicalPotentials[i + 1] - chemPotOutsideOfTheStep[i + 1]) * (mixingPartOfExchangeChemicalPotentials[i + 1] - chemPotOutsideOfTheStep[i + 1]);// !bio contaminant error

             }
            F[0] = (mixingPartOfExchangeChemicalPotentials[2] - chemPotOutsideOfTheStep[2]) * (mixingPartOfExchangeChemicalPotentials[2] - chemPotOutsideOfTheStep[2]);// !bio contaminant error

			F[1] = (insideOsmoticPressure - osmoticPressureOutsideOfTheStep) * (insideOsmoticPressure - osmoticPressureOutsideOfTheStep);
        }*/
	
        static void BorderEquations(double[] X, out double[] F, int L)
        {
			double[] _volumeFractions = new double[NumberOfComponents];
			double volumeFractionSum = 0;
			for (int i = 0; i < NumberOfComponents-1; i++)
			{
				_volumeFractions[i] = X[i];
				volumeFractionSum += X[i];
			}
			_volumeFractions[NumberOfComponents - 1] = 1 - volumeFractionSum;
			F = new double[L];
			double[] mixingPartOfExchangeChemicalPotentials;
			Lagrmix_PolA(NumberOfComponents, _volumeFractions, out mixingPartOfExchangeChemicalPotentials);
			double osmoticPressure = CalculateOsmoticPressure(_volumeFractions);
			//F[0] = (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]) * (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]);//solvent
			//F[1] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]);//bio
			F[0] = osmoticPressure * osmoticPressure;
			for(int i=1;i<NumberOfComponents-1;i++)
				F[i] = (mixingPartOfExchangeChemicalPotentials[i] - chemPotAtTheBorder[i]) * (mixingPartOfExchangeChemicalPotentials[i] - chemPotAtTheBorder[i]);//bio

		}
		static void BrushEquations(double[] X, out double[] F, int L)
        {
			double nu = 2;
			double y_cur = point_y;


			double volumeFractionSum = 0;
			for (int i = 0; i < NumberOfComponents - 1; i++)
			{
				fipolimer[i+1] = X[i];
				volumeFractionSum += X[i];
			}
			fipolimer[0] = 1 - volumeFractionSum;
			double[] mixingPartOfExchangeChemicalPotentials;
			//!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
			Lagrmix_PolA(NumberOfComponents, fipolimer, out mixingPartOfExchangeChemicalPotentials);
			F = new double[L];

			//F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
			for(int i=0;i< NumberOfComponents - 2; i++)
            {
				F[i] = (mixingPartOfExchangeChemicalPotentials[i+1] - chemPotInTheBulk[i + 1]) * (mixingPartOfExchangeChemicalPotentials[i + 1] - chemPotInTheBulk[i + 1]);// !bio contaminant error

			}
			F[NumberOfComponents-2] = Math.Pow((mixingPartOfExchangeChemicalPotentials[NumberOfComponents-1] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol),2);//!polymer error

			Console.WriteLine("BrushEquations values: " + F[0] + "  " +F[1] + " ----------------Volume fractions:   " + X[0] + "    " + X[1]);
			
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
		static List<KeyValuePair<double,List<double>>> CalculateMixingEnergyProfile(int AcomponentIndex, int BcomponentIndex, int numberOfPoints)
		{
			List<KeyValuePair<double, List<double>>> output = new List<KeyValuePair<double, List<double>>>();
			double step = 0.99 / (numberOfPoints - 1);
			double[] volumeFractions = new double[NumberOfComponents];
			//double volumeFractionSum = 0;
			for(int i = 0; i < NumberOfComponents; i++)
				volumeFractions[i] = 0;
			volumeFractions[AcomponentIndex] = 0.005;
            volumeFractions[BcomponentIndex] = 1- volumeFractions[AcomponentIndex];
			while (volumeFractions[AcomponentIndex] < 1)
			{
				List<double> value = new List<double>();
                double Fmix = mixingPartModule.CalculateMixingFreeEnergy(volumeFractions);
				double[] exchangeChemPotentials;
				Lagrmix(NumberOfComponents, volumeFractions, out exchangeChemPotentials);
				value.Add(Fmix);
				value.Add(exchangeChemPotentials[2]);
				output.Add(new KeyValuePair<double, List<double>>(volumeFractions[AcomponentIndex],value));
				volumeFractions[AcomponentIndex] += step;
				volumeFractions[BcomponentIndex] = 1 - volumeFractions[AcomponentIndex];
			}

			return output;
        }
		static double[] FindSegregationPointsBetweenSolventAndPolymer()
		{
            double ERREL = Math.Pow(10, -6);
            int ITMAX = 600;
            double[] XGUESS = new double[2];
            double[] X = new double[2];
			XGUESS[0] = 0.000001;
			XGUESS[1] = 0.999999;

            //XGUESS[0] = 0.02;
            //XGUESS[1] = 0.99;
            double FNORM;
            DNEQNF(SegregationEquations, ERREL, 2, ITMAX, XGUESS, out X, out FNORM);
			return X;
        }
        static void SegregationEquations(double[] X, out double[] F, int L)
        {
            double nu = 2;

            double[] mixingPartOfExchangeChemicalPotentials_0;
			double[] volumeFractions_0 = new double[3];
			volumeFractions_0[0] = X[0];
            volumeFractions_0[1] = 0;
            volumeFractions_0[2] = 1-X[0];
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            Lagrmix(3, volumeFractions_0, out mixingPartOfExchangeChemicalPotentials_0);

            double[] mixingPartOfExchangeChemicalPotentials_1;
            double[] volumeFractions_1 = new double[3];
            volumeFractions_1[0] = X[1];
            volumeFractions_1[1] = 0;
            volumeFractions_1[2] = 1 - X[1];
            //!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
            Lagrmix(3, volumeFractions_1, out mixingPartOfExchangeChemicalPotentials_1);
            F = new double[L];

			//F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
			F[0] = (mixingPartOfExchangeChemicalPotentials_0[1] - mixingPartOfExchangeChemicalPotentials_1[1]) * (mixingPartOfExchangeChemicalPotentials_0[1] - mixingPartOfExchangeChemicalPotentials_1[1]);
			double secondEquation = mixingPartModule.CalculateMixingFreeEnergy(volumeFractions_1) - mixingPartModule.CalculateMixingFreeEnergy(volumeFractions_0) - (volumeFractions_1[1] - volumeFractions_0[1]) * mixingPartOfExchangeChemicalPotentials_0[1];
			F[1] = Math.Pow(secondEquation,2);
        }

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
}
