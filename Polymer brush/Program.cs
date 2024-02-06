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
		static double R, BA, y_min, y_max, yacc, aA, rNA, rNB, sigma, Lamb_Pol;
		static double rN, nu, coe, pi, osmbulk, y_try, y_edge, y_cur, P_AB;
		static double[] lambda, Fibulk, Lagrbulk, Nal;
		static double[,] chi;
		static double[,] etas;
		static double[] Xbrush, fipolimer;
		static double point_y;
		static StreamWriter sw;
		static double z = 6;

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
			double u_bio=0;
			sw = new StreamWriter("profile.txt");
			sw.WriteLine("y_cur    solvent    bio    polymer    osm_pressure");

			Enter();
			y_try = y_min;
			y_edge = rtbis3(y_min, y_max, yacc);
			y_cur = y_min;

			//!Calculate polymer concentration at the start position(A / B boundary)
			P_AB = fiav(y_cur, y_edge); // !   Polymer concentration at A/ B boundary(highest)  must work before finding profile to set Lagr multipliers in common block
           
			
			while (y_cur < y_edge)
            {
				//y_cur += aA / R;
				//!write(*, *) 'Phi poly subroutine is called'
				Xbrush = new double[3];
				FI_POLI(out Xbrush, y_cur); //  !calculates concentration profile in the brush after the solution is found

				
				fipolimer[1] = Xbrush[0];// ! local compsition of the brush at point y_cur(1)-solvent(2) - biocomponent(3) - polymer
				fipolimer[2] = Xbrush[1];
				fipolimer[0] = 1.0 - fipolimer[1] - fipolimer[2];

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
			sigma = (7.0 * Math.Pow(10, -9)) * (7.0 * Math.Pow(10, -9)) / 0.12;
			sigma = (7.0 * Math.Pow(10, -9)) * (7.0 * Math.Pow(10, -9)) / 0.6;
			R = rNB * (nu + 1) * aA * aA * aA / sigma; // core radius((nu+1.0)*rNB / SIGMA)**2    LAGRANGE_STR1 = BA * (y_try - 1.0) * *2 * ((nu + 1.0) * rNB / SIGMA) * *2 * aA * *3  CHECK THIS
													   //write(*, *) R,sigma
													   //
			BA = coe / ((rNA * aA) * (rNA * aA));

			y_min = 1.0 + 0.1 * aA / R;
			y_max = 1.00001 * (1.0 + 1.0 * aA * rNA / R);                 //y_max = 1.0d0 * (1.0d0 + 2.0 * aA * rNA / R)
			yacc = Math.Pow(10, -8);

			//! sigma min for all morphologies:
			//!sigma_MIN = 10.0 * (nu + 1.1) * aA * *2
			//!sigma_MAX = 3.0d0 * aA * *2 * (3.1415926 * 4.0 * rNB * *2 / 3.0d0) **(1.0 / 3.0)
			//				 !write(*, *) BA,R,rNB* aA, sigma,10.0 * (nu + 1.1) * aA * *2,y_max,BA * (R * (y_max - 1.0)) * *2,
			//!stop

			Nal = new double[3];
			Nal[0] = 1.0;// ! solvent
			Nal[1] = 3.0;// ! bioadditive
			Nal[2] = rNA;// polymer

			chi = new double[3, 3];

			chi[0, 1] = 1.0;//! solv - bio
			chi[0, 2] = 0.0;//! solv - polym
			chi[1, 2] = -0.8;//d0! bio - polym

			for (int i = 0; i < 3; i++)
			{
				chi[i, i] = 0;
				for (int j = i + 1; j < 3; j++)
					chi[j, i] = chi[i, j];
			}

			etas = new double[3, 3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					etas[i, j] = Math.Exp(-chi[i, j] / z);

			//!give bulk composition and calculate Lagr.multipliers and osm pressure in the bulk:
			Fibulk = new double[3];
			Fibulk[0] = 0.98;
			Fibulk[1] = 1.0 - Fibulk[0];
			Fibulk[2] = 0.0;

			lambda = new double[3];
			fipolimer = new double[3];
			Lagrmix(2, Fibulk, out Lagrbulk);
			lambda[1] = Lagrbulk[1];//this is for biocomponent
			osmbulk = Osmmix(2, Fibulk);// !  this is for solvent

			lambda[0] = Lagrbulk[0]; //!  this should be identiacally zero

			//! now calculate  mixing part of Lagrange multiplier at the edge of the brush, where  FiA = 0:

			Fibulk[2] = 0.0;

			Lagrmix_PolA(3, Fibulk, out Lagrbulk);// ! NOTE  that  Lagr multipliers are spoiled in  Lagrbulk   but stored in common

			lambda[2] = Lagrbulk[2];

		}
		static void Lagrmix(int numberOfComponents, double[] X, out double[] DfmixDfi)
		{
			DfmixDfi = new double[3];
			double[] AlternativeDfmixDfi = new double[3];
			for (int i = 0; i < numberOfComponents; i++)
			{
				double sum = 0;
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[j] * (chi[i, j] - chi[0, j]);
				DfmixDfi[i] = (Math.Log(X[i]) + 1.0) / Nal[i] - (Math.Log(X[0]) + 1.0) / Nal[0] + sum;// ! dummy for solvent identically 0
				AlternativeDfmixDfi[i] = CalculateExchangeChemialPotentialOfComponent(X, i);
				DfmixDfi[i] = AlternativeDfmixDfi[i];                                                                                     //AlternativeDfmixDfi(i) = CalculateExchangeChemialPotentialOfComponent(3, X, i)
																									  //DfmixDfi(i) = AlternativeDfmixDfi(i)
			}
		}
		static void Lagrmix_PolA(int numberOfComponents, double[] X, out double[] DfmixDfi)
		{
			DfmixDfi = new double[3];
			double[] AlternativeDfmixDfi = new double[3];

			for (int i = 0; i < numberOfComponents; i++)
			{
				if (X[i] < 0)
				{
					for (int n = 0; n < numberOfComponents; n++)
						DfmixDfi[n] = 1050;
					return;
				}

				double sum = 0;
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[j] * (chi[i, j] - chi[0, j]);
				if (i == numberOfComponents - 1)
					DfmixDfi[i] = -(Math.Log(X[0]) + 1.0) / Nal[0] + sum;// ! dummy for solvent identically 0
				else
					DfmixDfi[i] = (Math.Log(X[i]) + 1.0) / Nal[i] - (Math.Log(X[0]) + 1.0) / Nal[0] + sum;// ! dummy for solvent  identically 0
				
				AlternativeDfmixDfi[i] = CalculateExchangeChemialPotentialOfComponent(X, i);
				DfmixDfi[i] = AlternativeDfmixDfi[i];
			}
		}
		static double Osmmix(int numberOfComponents, double[] X)
		{
			double sum = 0;
			double sum1 = 0;
			for (int i = 0; i < numberOfComponents; i++)
			{
				sum1 = sum1 + X[i] * (1.0 - 1.0 / Nal[i]);
				for (int j = 0; j < numberOfComponents; j++)
					sum += X[i] * X[j] * (chi[0, j] - chi[i, j] / 2.0);
			}
			return Math.Log(X[0]) + sum1 + sum;
		}
		static double normal(double y)
		{
			double nu = 2.0;
			double aINT = 1.0;
			double bINT = y;
			double s = 1;
			qtrap(aINT, bINT, out s);
			return s - rNA / (rNB * (nu + 1.0));

		}
		static void qtrap(double aINT, double bINT, out double s)
        {
			double EPS = 0.1;
			int JMAX = 8;
			double old_s = -1 * Math.Pow(10, -30);
			s = 0;
			for(int n = 0; n < JMAX;n++)
            {
				if (n == 4)
					n = 4;
				s = trapzd(aINT, bINT, s, n+1);
                if (n > 4)
					if (Math.Abs(s - old_s) < EPS * Math.Abs(old_s) || (s == 0 && old_s == 0))
						return;
				old_s = s;
            }
			Console.WriteLine("Too many steps in q trap");
			throw new Exception();


		}
		static double trapzd(double aINT,double bINT,double s, int n)
        {
			//s = 0;
			if (n == 1)
				s = 0.5 * (bINT - aINT) * (fiav(aINT, bINT) + fiav(bINT, bINT));
            else
            {
				int it = (int)Math.Pow(2, n - 2);
				double tnm = it;
				double del = (bINT - aINT) / tnm;
				double x = aINT + 0.5 * del;
				double sum = 0;
				for(int counter = 0; counter < it; counter++)
                {
					sum += fiav(x, bINT);
					x += del;
                }
				s = 0.5 * (s + (bINT - aINT) * sum / tnm);
			}
			return s;
		}
		static double fiav(double y_cur, double bINT)
        {
			int L = 2;
			double nu = 2;
			double y_try = bINT;
			Lamb_Pol = lambda[2] + BA * (R * (y_try - 1)) * (R * (y_try - 1));
			FI_POLI(out Xbrush, y_cur);
			double fay = Xbrush[1];
			return fay * y_cur * y_cur;

		}
		static double rtbis3(double y_min, double y_max, double yacc)
		{
			double output = 0;
			int j_max = 100; //!Maximum allowed number of bisections.
							 //!Using bisection, find the root of a function func known to lie between x1 and x2.The
							 //!root, returned as rtbis, will be refined until its accuracy is ï¿½xacc.
			double dx, f, fmid, xmid;
			//!calculating function value at y_max
			fmid = normal(y_max);
			f = normal(y_min);// !R2(beta2, lambda, itog2)

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
				fmid = normal(xmid);// !call R3(beta3, lambda, itog3)
				if (fmid < 0)
					output = xmid;
				if (Math.Abs(dx) < yacc || fmid == 0)
					return output;
			}
			Console.WriteLine("Too many iterations");
			throw new Exception();
		}
		static void FI_POLI(out double[] XBrush, double y_cur)
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

			DNEQNF(EU_LA, ERREL, L, ITMAX, XBrushGUESS, out XBrush, out FNORM);

			
		}
		delegate void NonlinearSystem(double[] XBrush, out double[] F, int L);
		static void DNEQNF(NonlinearSystem Func, double ERREL, int L, int ITMAX, double[] XGuess, out double[] X, out double FNORM)
		{
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

				if (double.IsNaN(X[0]))
					;
				FNORM = 0;
				for (int i = 0; i < L; i++)
					FNORM += F[i] * F[i];
			}
			
        }
		static bool ContainsOutrangeValues(double[] X)
        {
			foreach (double x in X)
				if (x < 0 || x > 1)
					return true;
			return false;
		}
		static void EU_LA(double[] XBrush, out double[] F, int L)
        {
			double nu = 2;
			double y_cur = point_y;
			fipolimer[1] = Xbrush[0];
			fipolimer[2] = Xbrush[1];
			fipolimer[0] = 1 - fipolimer[1] - fipolimer[2];
			double[] DfmixDfi;
			//!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
			Lagrmix_PolA(3, fipolimer, out DfmixDfi);
			F = new double[L];

			F[0] = (DfmixDfi[1] - lambda[1]) * (DfmixDfi[1] - lambda[1]);// !bio contaminant error
			F[1] = Math.Pow((DfmixDfi[2] + BA * (R * (y_cur - 1.0)) * (R * (y_cur - 1.0)) - Lamb_Pol),2);//!polymer error

			Console.WriteLine("EU_LA " + F + "    " + Xbrush[0] + "    " + Xbrush[1]);

		}
		static void TestSystem(double[] X, out double[] F, int L)
        {
			F = new double[L];
			F[0] = X[0] * X[0] + X[1] * X[1] - 1;
			F[1] = X[1] - X[0] * X[0];
		}
		static double CalculateMixingFreeEnergy(double[] X)
        {
			double a = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Nal[1];
			double b = chi[1, 2] * X[1] * X[2];
			double c = chi[0, 1] * X[0] * X[1];
			double d = chi[0, 2] * X[0] * X[2];
			return a + b + c + d;
		}
		static double CalculateGugenheimMixingFreeEnergy(double[] X)
		{
			int n = X.Length;
			double translationSum = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Nal[1];
			double[] XX = CalculateGugenheimCorrelations(X, etas);
			double mixingSum = 0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j > i; j++)
					mixingSum += chi[i, j] * XX[i] * XX[j] * X[i] * X[j] * etas[i, j];

			/*double b = chi[1, 2] * X[1] * X[2];
			double c = chi[0, 1] * X[0] * X[1];
			double d = chi[0, 2] * X[0] * X[2];*/
			double output = translationSum + mixingSum;
			return output;
		}
		static double CalculateExchangeChemialPotentialOfComponent(double[] X, int componenIndex)
        {
			//double f = CalculateMixingFreeEnergy(X);
			double f = CalculateGugenheimMixingFreeEnergy(X);
			double x = X[componenIndex];
			double max_dx = 1 - x;
            if (X[0] < max_dx)
				max_dx = X[0];
			double dx = 0.01*x;
			if (dx == 0)
				dx = 0.01;
			if (dx > max_dx)
				dx = max_dx;

			double oldSolventVolumeFraction = X[0];
			double x_dx = x + dx;
			X[componenIndex] = x_dx;
			X[0] -= dx;
			//double f_df = CalculateMixingFreeEnergy(X);
			double f_df = CalculateGugenheimMixingFreeEnergy(X);
			X[0] = oldSolventVolumeFraction;
			X[componenIndex] = x;
			return (f_df - f) / dx;


		}
		static double[] CalculateGugenheimCorrelations(double[] alphas, double[,] etas)
        {
			int n = alphas.Length;
			double[] XX = new double[n];
			double[] newXX = new double[n];
			double initialGuess = 1;
			for(int i = 0; i < n; i++)
            {
				XX[i] = initialGuess;
				initialGuess -= 0.0001;
			}
			bool converged = false;
            while (!converged)
            {
				for (int i = 0; i < n; i++)
				{
					double sum = 0;
					for (int j = 0; j < n; j++)
						sum += alphas[j] * XX[j] * etas[i, j];
					newXX[i] = 1 / sum;
				}
				for (int i = 0; i < n; i++)
					XX[i] = (XX[i] + newXX[i]) / 2;
				converged = true;
				for (int i = 0; i < n; i++)
					if (Math.Abs(XX[i] - newXX[i]) > Math.Pow(10, -12))
					{
						converged = false;
						break;
					}
			}
			return XX;

		}
	} 
}
