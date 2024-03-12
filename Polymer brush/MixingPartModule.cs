using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Polymer_brush
{
    class MixingPartModule
    {
		/*public double CalculateDerivativesOfMixingFreeEnergyWithRespectToVolumeFractions(double[] X, int componentIndex)
        {
			double f = CalculateMixingFreeEnergy(X);
			//double f = CalculateGugenheimMixingFreeEnergy(X);
			double x = X[componentIndex];
			double max_dx = 1 - x;
			if (X[0] < max_dx)
				max_dx = X[0];
			double dx = 0.01 * x;
			if (dx == 0)
				dx = 0.01;
			if (dx > max_dx)
				dx = max_dx;

			double oldSolventVolumeFraction = X[0];
			double x_dx = x + dx;
			X[componentIndex] = x_dx;
			X[0] -= dx;
			double f_df = CalculateMixingFreeEnergy(X);
			//double f_df = CalculateGugenheimMixingFreeEnergy(X);
			X[0] = oldSolventVolumeFraction;
			X[componentIndex] = x;
			return (f_df - f) / dx;
		}*/
		public double CalculateMixingFreeEnergy(double[] X)
        {
			return CalculateFloryMixingFreeEnergy(X);
			//return CalculateGugenheimMixingFreeEnergy(X);
		}
		private double[] CalculateFunctionalGroupsMolarFractions(double[] XofMolecules)
        {
			double[] functionalGroupsVolumeFraction = new double[Program.chiMatrixSize];
			for (int i = 0; i < Program.NumberOfComponents - 1; i++)
				functionalGroupsVolumeFraction[i] = XofMolecules[i];
			for (int i = Program.NumberOfComponents - 1; i < Program.chiMatrixSize; i++)
				functionalGroupsVolumeFraction[i] = XofMolecules[Program.NumberOfComponents - 1] * Program.fractionsOfGroups[i - (Program.NumberOfComponents - 1)];
			return functionalGroupsVolumeFraction;
		}
		private double CalculateFloryMixingFreeEnergy(double[] XofMolecules)
		{
			double[] X = CalculateFunctionalGroupsMolarFractions(XofMolecules);
			double a = 0;
			for(int i=0;i<Program.NumberOfComponents;i++)
				if(XofMolecules[i]!=0)
					a += XofMolecules[i] * Math.Log(XofMolecules[i]) / Program.size[i];
			//a=X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Program.size[1];
			/*double b = Program.chi[1, 2] * X[1] * X[2];
			double c = Program.chi[0, 1] * X[0] * X[1];
			double d = Program.chi[0, 2] * X[0] * X[2];*/
			for (int i = 0; i < Program.chiMatrixSize; i++)
				for (int j = 0;j < Program.chiMatrixSize; j++)
					if(j>i && X[i] != 0	&& X[j] != 0)
					{
						if(i==0&&j==2)
                            a += (Program.chi[i, j] + 1.5 * X[2] * X[2]) * X[i] * X[j];
						else
							a += Program.chi[i, j] * X[i] * X[j];
                    }
			//return a + b + c + d;
			return a;
		}
		private double CalculateGugenheimMixingFreeEnergy(double[] XofMolecules)
		{
			double[] X = CalculateFunctionalGroupsMolarFractions(XofMolecules);
			int n = X.Length;
			double translationSum = 0;
			for (int i = 0; i < Program.NumberOfComponents - 1; i++)
				if(X[i]!=0)
					translationSum += XofMolecules[i] * Math.Log(XofMolecules[i]) / Program.size[i]; 
			double[] XX = CalculateGugenheimCorrelations(X, Program.etas);
			double mixingSum = 0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < i; j++)
                    if (X[i] != 0 && X[j] != 0)
                        mixingSum += Program.chi[i, j] * XX[i] * XX[j] * X[i] * X[j] * Program.etas[i, j];

			/*double b = chi[1, 2] * X[1] * X[2];
			double c = chi[0, 1] * X[0] * X[1];
			double d = chi[0, 2] * X[0] * X[2];*/
			double output = translationSum + mixingSum;
			return output;
		}
		private double CalculateExchangeChemialPotentialOfComponentAnaliticaly(double[] X, int componenIndex)
		{
			if (componenIndex == 0 || componenIndex == 1)
				return 0;
			double output = 0;
			output += (Math.Log(X[2]) + 1) / Program.size[2];
            output -= (Math.Log(X[0]) + 1) / Program.size[0];
			output += Program.chi[0, 2] * X[0];
			output -= Program.chi[0, 2] * X[2];
			return output;
        }

        public double CalculateExchangeChemialPotentialOfComponent(double[] X, int componenIndex)
		{
			
			double f = CalculateMixingFreeEnergy(X);
			//double f = CalculateGugenheimMixingFreeEnergy(X);
			double x = X[componenIndex];
			double max_dx = 1 - x;
			if (X[0] < max_dx)
				max_dx = X[0];
			double dx = 0.01 * x;
			if (dx == 0)
				dx = 0.01;
			if (dx > max_dx)
			{
                dx = max_dx / 2;
				if (max_dx == 0)
					dx = -0.01;
            }

			double oldSolventVolumeFraction = X[0];
			double x_dx = x + dx;
			X[componenIndex] = x_dx;
			X[0] -= dx;
			double f_df = CalculateMixingFreeEnergy(X);
			//double f_df = CalculateGugenheimMixingFreeEnergy(X);
			X[0] = oldSolventVolumeFraction;
			X[componenIndex] = x;
			double output = (f_df - f) / dx;
			//double analiticalResult = CalculateExchangeChemialPotentialOfComponentAnaliticaly(X, componenIndex);
            return output;


		}
		public double[] CalculateGugenheimCorrelations(double[] alphas, double[,] etas)
		{
			int n = alphas.Length;
			double[] XX = new double[n];
			double[] newXX = new double[n];
			double initialGuess = 1;
			for (int i = 0; i < n; i++)
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
