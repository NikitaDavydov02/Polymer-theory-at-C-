﻿using System;
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
			//return CalculateFloryMixingFreeEnergy(X);
			return CalculateGugenheimMixingFreeEnergy(X);
		}
		private double CalculateFloryMixingFreeEnergy(double[] X)
		{
		
			double a = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Program.size[1];
			double b = Program.chi[1, 2] * X[1] * X[2];
			double c = Program.chi[0, 1] * X[0] * X[1];
			double d = Program.chi[0, 2] * X[0] * X[2];
			return a + b + c + d;
		}
		private double CalculateGugenheimMixingFreeEnergy(double[] X)
		{
			int n = X.Length;
			double translationSum = X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Program.size[1];
			double[] XX = CalculateGugenheimCorrelations(X, Program.etas);
			double mixingSum = 0;
			for (int i = 0; i < n; i++)
				for (int j = 0; j < i; j++)
					mixingSum += Program.chi[i, j] * XX[i] * XX[j] * X[i] * X[j] * Program.etas[i, j];

			/*double b = chi[1, 2] * X[1] * X[2];
			double c = chi[0, 1] * X[0] * X[1];
			double d = chi[0, 2] * X[0] * X[2];*/
			double output = translationSum + mixingSum;
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
				dx = max_dx/2;

			double oldSolventVolumeFraction = X[0];
			double x_dx = x + dx;
			X[componenIndex] = x_dx;
			X[0] -= dx;
			double f_df = CalculateMixingFreeEnergy(X);
			//double f_df = CalculateGugenheimMixingFreeEnergy(X);
			X[0] = oldSolventVolumeFraction;
			X[componenIndex] = x;
			return (f_df - f) / dx;


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