using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace Polymer_brush
{
	class Program
	{
		public static double R, BA, y_min, y_max, yacc, aA, rNA, rNB, areaPerChain, Lamb_Pol;
		public static double rN, nu, coe, pi, osmbulk, y_edge, y_cur, P_AB;
		public static double[] chemPotInTheBulk, volumeFractionsInTheBulk, Lagrbulk,  size;
		public static double[,] chi;
		public static double[,] etas;
		static double[] Xbrush, fipolimer;
		static double point_y;
		static StreamWriter sw;
		static double z = 6;
		//static double[] chemPotInTheBulk;
		static double[] chemPotAtTheBorder;
		static MixingPartModule mixingPartModule;
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
			sw = new StreamWriter("profile.txt");
			sw.WriteLine("y_cur    solvent    bio    polymer    osm_pressure");

			Enter();
			y_edge = BisectionSolve(y_min, y_max, yacc);
			y_cur = y_min;

			//!Calculate polymer concentration at the start position(A / B boundary)
			P_AB = FindSubintegralValueForNormalization(y_cur, y_edge); // !   Polymer concentration at A/ B boundary(highest)  must work before finding profile to set Lagr multipliers in common block
           
			
			while (y_cur < y_edge)
            {
				//y_cur += aA / R;
				//!write(*, *) 'Phi poly subroutine is called'
				Xbrush = new double[3];
				FindVolumeFractionsInTheBrushForPoint(out Xbrush, y_cur); //  !calculates concentration profile in the brush after the solution is found

				
				fipolimer[1] = Xbrush[0];//biocomponent // ! local compsition of the brush at point y_cur(1)-solvent(2) - biocomponent(3) - polymer
				fipolimer[2] = Xbrush[1];//polymer
				fipolimer[0] = 1.0 - fipolimer[1] - fipolimer[2];//solvent

				Console.Write("Fmix="+Fmix);
				Console.Write("u_sol=" +u_sol);
				Console.Write("u_pol=" + u_pol); 
				Console.Write("u_bio=" + u_bio);
				Console.Write("/n");

				sw.WriteLine(y_cur.ToString()+"    "+fipolimer[0] + "    " + fipolimer[1] + "    " + fipolimer[2] + "    " + (Osmmix(3, fipolimer) - osmbulk));
				y_cur += aA / R;
			}

			Console.WriteLine("Beta: " + y_edge);
			Console.WriteLine("Calculation is done!");
			

			sw.Close();
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
			y_max = 1.00001 * (1.0 + 1.0 * aA * rNA / R);                 //y_max = 1.0d0 * (1.0d0 + 2.0 * aA * rNA / R)
			yacc = Math.Pow(10, -8);

			//! areaPerChain min for all morphologies:
			//!areaPerChain_MIN = 10.0 * (nu + 1.1) * aA * *2
			//!areaPerChain_MAX = 3.0d0 * aA * *2 * (3.1415926 * 4.0 * rNB * *2 / 3.0d0) **(1.0 / 3.0)
			//				 !write(*, *) BA,R,rNB* aA, areaPerChain,10.0 * (nu + 1.1) * aA * *2,y_max,BA * (R * (y_max - 1.0)) * *2,
			//!stop

			 size = new double[3];
			 size[0] = 1.0;// ! solvent
			 size[1] = 3.0;// ! bioadditive
			 size[2] = rNA;// polymer

			chi = new double[3, 3];

			//solv
			//bio
			//pol

			chi[0, 1] = 0;//! solv - bio
			chi[0, 2] = 5;//! solv - polym
			chi[1, 2] = 0;//d0! bio - polym
			/* does not work
			 * chi[0, 1] = 1;//! solv - bio
			chi[0, 2] = 1;//! solv - polym
			chi[1, 2] = -0.8;//d0! bio - polym*/


			/*	work
			 *	chi[0, 1] = 1;//! solv - bio
				chi[0, 2] = 0.6;//! solv - polym
				chi[1, 2] = -0.8;//d0! bio - polym*/

			for (int i = 0; i < 3; i++)
			{
				chi[i, i] = 0;
				for (int j = i + 1; j < 3; j++)
					chi[j, i] = chi[i, j];
			}

			etas = new double[3, 3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					etas[i, j] = Math.Exp(-3*chi[i, j]/z);

			//!give bulk composition and calculate Lagr.multipliers and osm pressure in the bulk:
			volumeFractionsInTheBulk = new double[3];

			volumeFractionsInTheBulk[0] = 0.999;
			volumeFractionsInTheBulk[1] = 1.0 - volumeFractionsInTheBulk[0];
			volumeFractionsInTheBulk[2] = 0.0;

			chemPotInTheBulk = new double[3];
			chemPotAtTheBorder = new double[3];
			fipolimer = new double[3];
			//Lagrmix(2, volumeFractionsInTheBulk, out Lagrbulk);
			Lagrmix(2, volumeFractionsInTheBulk, out chemPotInTheBulk);
			//chemPotInTheBulk[1] = Lagrbulk[1];//this is for biocomponent
			osmbulk = Osmmix(2, volumeFractionsInTheBulk);// !  this is for solvent

			//chemPotInTheBulk[0] = Lagrbulk[0]; //!  this should be identiacally zero
			//! now calculate  mixing part of Lagrange multiplier at the edge of the brush, where  FiA = 0:
			//volumeFractionsInTheBulk[2] = 0.0;
			//Lagrmix_PolA(3, volumeFractionsInTheBulk, out Lagrbulk);// ! NOTE  that  Lagr multipliers are spoiled in  Lagrbulk   but stored in common
			//chemPotInTheBulk[2] = Lagrbulk[2];

		}
		static void Lagrmix(int numberOfComponents, double[] X, out double[] mixingPartOfExchangeChemicalPotentials)
		{
			mixingPartOfExchangeChemicalPotentials = new double[3];
			double[] AlternativemixingPartOfExchangeChemicalPotentials = new double[3];
			for (int i = 0; i < numberOfComponents; i++)
			{
				double sum = 0;
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[j] * (chi[i, j] - chi[0, j]);
				mixingPartOfExchangeChemicalPotentials[i] = (Math.Log(X[i]) + 1.0) /  size[i] - (Math.Log(X[0]) + 1.0) /  size[0] + sum;// ! dummy for solvent identically 0
				AlternativemixingPartOfExchangeChemicalPotentials[i] = mixingPartModule.CalculateExchangeChemialPotentialOfComponent(X, i);
				mixingPartOfExchangeChemicalPotentials[i] = AlternativemixingPartOfExchangeChemicalPotentials[i];                                                                                     //AlternativemixingPartOfExchangeChemicalPotentials(i) = CalculateExchangeChemialPotentialOfComponent(3, X, i)
			}
		}
		static void Lagrmix_PolA(int numberOfComponents, double[] X, out double[] mixingPartOfExchangeChemicalPotentials)
		{
			mixingPartOfExchangeChemicalPotentials = new double[3];
			double[] AlternativemixingPartOfExchangeChemicalPotentials = new double[3];

			for (int i = 0; i < numberOfComponents; i++)
			{
				if (X[i] < 0)
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
		static double NormalizationFunctionValue(double y)
		{
			
			double nu = 2.0;
			double norm = rNA / (rNB * (nu + 1.0));
			double integrationMin = 1.0;
			double integrationMax = y;
			double s = 1;


			//Find chemical potentials at border;
			chemPotAtTheBorder[0] = chemPotInTheBulk[0];//solvent
			chemPotAtTheBorder[1] = chemPotInTheBulk[1];//bio
			//Finding volume fractions at the border
			double ERREL = Math.Pow(10, -4);
			double[] XBorderGUESS = new double[2];
			XBorderGUESS[0] = volumeFractionsInTheBulk[0];//this is the fraction of solvent at the border
			XBorderGUESS[1]= volumeFractionsInTheBulk[1];//this is the fraction of biocomponnt at the border
			//XBorderGUESS[0] = 0.1;

			double FNORM;
			double[] _XBorder = new double[3];
			DNEQNF(BorderEquations, ERREL, 2, 1000, XBorderGUESS,out _XBorder, out FNORM);

			double[] volFractionsAtTheBorder = new double[3];
			volFractionsAtTheBorder[0] = _XBorder[0];
			volFractionsAtTheBorder[1] = _XBorder[1];
			volFractionsAtTheBorder[2] = 1 - volFractionsAtTheBorder[0] - volFractionsAtTheBorder[1];
			Lagrmix_PolA(3, volFractionsAtTheBorder, out chemPotAtTheBorder);

			chemPotAtTheBorder[2] += BA * (R * (integrationMax - 1)) * (R * (integrationMax - 1));
			Lamb_Pol = chemPotAtTheBorder[2];

			CalculateNormalizationIntegral(integrationMin, integrationMax, out s);
			return s -norm;

		}
		static void CalculateNormalizationIntegral(double integrationMin, double integrationMax, out double s)
        {
			double EPS = 0.1;
			int JMAX = 8;
			double old_s = -1 * Math.Pow(10, -30);
			s = 0;
			for(int n = 0; n < JMAX;n++)
            {
				s = CalculateNormalizationIntegralWithDefeniteNumberOfTrapezoids(integrationMin, integrationMax, s, n+1);
				
                if (n > 4)
					if (Math.Abs(s - old_s) < EPS * Math.Abs(old_s) || (s == 0 && old_s == 0))
						return;
				old_s = s;
            }
			Console.WriteLine("Too many steps in q trap");
			throw new Exception();


		}
		static double CalculateNormalizationIntegralWithDefeniteNumberOfTrapezoids(double integrationMin,double integrationMax,double s, int n)
        {
			//s = 0;
			if (n == 1)
				s = 0.5 * (integrationMax - integrationMin) * (FindSubintegralValueForNormalization(integrationMin, integrationMax) + FindSubintegralValueForNormalization(integrationMax, integrationMax));
            else
            {
				int it = (int)Math.Pow(2, n - 2);
				double tnm = it;
				double del = (integrationMax - integrationMin) / tnm;
				double x = integrationMin + 0.5 * del;
				double sum = 0;
				for(int counter = 0; counter < it; counter++)
                {
					double add = FindSubintegralValueForNormalization(x, integrationMax);
					if (double.IsNaN(add))
						;
					sum += add;
					x += del;
                }
				s = 0.5 * (s + (integrationMax - integrationMin) * sum / tnm);
			}
			return s;
		}
		static double FindSubintegralValueForNormalization(double y_cur, double integrationMax)
        {
			int L = 2;
			double nu = 2;
			//Lamb_Pol = chemPotInTheBulk[2] + BA * (R * (  - 1)) * (R * (  - 1));
			FindVolumeFractionsInTheBrushForPoint(out Xbrush, y_cur);
			double fay = Xbrush[1];
			if (double.IsNaN(fay))
				;
			return fay * y_cur * y_cur;

		}
		static double BisectionSolve(double y_min, double y_max, double yacc)
		{
			double output = 0;
			int j_max = 100; //!Maximum allowed number of bisections.
							 //!Using bisection, find the root of a function func known to lie between x1 and x2.The
							 //!root, returned as rtbis, will be refined until its accuracy is ï¿½xacc.
			double dx, f, fmid, xmid;
			//!calculating function value at y_max
			fmid = NormalizationFunctionValue(y_max);
			f = NormalizationFunctionValue(y_min);// !R2(beta2, chemPotInTheBulk, itog2)

			if (f * fmid > 0)
			{
				Console.WriteLine("Bisection error");
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
				fmid = NormalizationFunctionValue(xmid);// !call R3(beta3, chemPotInTheBulk, itog3)
				if (fmid < 0)
					output = xmid;
				if (Math.Abs(dx) < yacc || fmid == 0)
					return output;
			}
			Console.WriteLine("Too many iterations");
			throw new Exception();
		}
		static void FindVolumeFractionsInTheBrushForPoint(out double[] XBrush, double y_cur)
        {
			//!Solving(localy) system of non - linear equations

			int L = 2;
			double ERREL = Math.Pow(10, -4);
			point_y = y_cur;
			int ITMAX = 600;
			double[] XBrushGUESS = new double[L];
			XBrush = new double[L];
			XBrushGUESS[0] = Math.Pow(10, -8);//this is the fraction of biocomponent in the brush
			XBrushGUESS[1]=0.97;//this is the fraction of polymer in the brush
			double FNORM;

			DNEQNF(BrushEquations, ERREL, L, ITMAX, XBrushGUESS, out XBrush, out FNORM);

			
		}
		delegate void NonlinearSystem(double[] X, out double[] F, int L);
		static void DNEQNF(NonlinearSystem Func, double ERREL, int L, int ITMAX, double[] XGuess, out double[] X, out double FNORM)
		{
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
						Func(X, out F, L);
						double f_init = F[i];
						double old_x = X[j];
						double dx = old_x * 0.01;
						X[j] += dx;


						for (double devisionStepDegree = 1; X[j]>=1; devisionStepDegree++)
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
					throw new Exception(" Newton method did not manage to find solution for system of equations");
                }
			}
			
        }
		static bool ContainsOutrangeValues(double[] X)
        {
			foreach (double x in X)
				if (x < 0 || x > 1)
					return true;
			return false;
		}
		static void BorderEquations(double[] X, out double[] F, int L)
        {
			
			double[] _volumeFractions = new double[3];
			_volumeFractions[0] = X[0];//solv
			_volumeFractions[1] = X[1];//bio
			_volumeFractions[2] = 1 - _volumeFractions[0] - _volumeFractions[1];
			F = new double[L];
			double[] mixingPartOfExchangeChemicalPotentials;
			Lagrmix_PolA(3, _volumeFractions, out mixingPartOfExchangeChemicalPotentials);
			double osmoticPressure = 0;
			osmoticPressure =mixingPartModule.CalculateMixingFreeEnergy(_volumeFractions);
			for (int i = 1; i < 3; i++)
				osmoticPressure -= _volumeFractions[i] * mixingPartOfExchangeChemicalPotentials[i];
			F[0] = (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]) * (mixingPartOfExchangeChemicalPotentials[0] - chemPotAtTheBorder[0]);//solvent
			F[1] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotAtTheBorder[1]);//bio
			F[0] = osmoticPressure * osmoticPressure;
		}
		static void BrushEquations(double[] X, out double[] F, int L)
        {
			double nu = 2;
			double y_cur = point_y;
			fipolimer[1] = X[0];
			fipolimer[2] = X[1];
			fipolimer[0] = 1 - fipolimer[1] - fipolimer[2];
			double[] mixingPartOfExchangeChemicalPotentials;
			//!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
			Lagrmix_PolA(3, fipolimer, out mixingPartOfExchangeChemicalPotentials);
			F = new double[L];

			F[0] = (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]) * (mixingPartOfExchangeChemicalPotentials[1] - chemPotInTheBulk[1]);// !bio contaminant error
			F[1] = Math.Pow((mixingPartOfExchangeChemicalPotentials[2] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol),2);//!polymer error

			Console.WriteLine("BrushEquations values: " + F[0] + "  " +F[1] + " ----------------Volume fractions:   " + X[0] + "    " + X[1]);

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
