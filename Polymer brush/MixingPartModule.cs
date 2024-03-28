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
		public double[] segregationPoints { get; private set; }
		public double[] segregationMixingEnergies;
		public MixingPartModule()
		{
			FindSegregationPointsBetweenSolventAndPolymer();
		}
		public double CalculateMixingFreeEnergy(double[] X)
        {
			if (segregationPoints != null && (X[1] >= segregationPoints[0] && X[1] < segregationPoints[1]))
			{
				
				double output = segregationMixingEnergies[0] + (X[1] - segregationPoints[0]) * (segregationMixingEnergies[1] - segregationMixingEnergies[0]) / (segregationPoints[1] - segregationPoints[0]);
				//output += X[1] * Math.Log(X[1]) / Program.size[1];
                //output += Program.chi[0, 1] * X[0] * X[1];
                //output += Program.chi[0, 2] * X[2] * X[1];
                return output;
			}
            return CalculateFloryMixingFreeEnergy(X);
            return CalculateGugenheimMixingFreeEnergy(X);

        }
		private double[] CalculateFunctionalGroupsMolarFractions(double[] XofMolecules)
        {
			double[] functionalGroupsVolumeFraction = new double[Program.chiMatrixSize];
			functionalGroupsVolumeFraction[0]= XofMolecules[0];
			
			//Polymer groups
            for (int i = 0; i < Program.NumberOfPolymerGroupTypes; i++)
                functionalGroupsVolumeFraction[1+i] = XofMolecules[1] * Program.fractionsOfGroups[i];

			//Additives in solution
            for (int i = 2; i < Program.NumberOfComponents; i++)
                functionalGroupsVolumeFraction[i-1+Program.NumberOfPolymerGroupTypes] = XofMolecules[i];
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
						if(i==0&&j==1)
                            a += (Program.chi[i, j] + Program.c * X[1] * X[1]) * X[i] * X[j];
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
			for (int i = 0; i < Program.NumberOfComponents; i++)
				if(X[i]!=0 && i!=1)
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
		/*private double CalculateExchangeChemialPotentialOfComponentAnaliticaly(double[] X, int componenIndex)
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
		*/
        public double CalculateExchangeChemialPotentialOfComponent(double[] X, int componenIndex)
		{
			double relativeDelta = Math.Pow(10, -2);
			double f0 = CalculateMixingFreeEnergy(X);
			double x = X[componenIndex];
			double max_dx = 1 - x;
			if (X[0] < max_dx)
				max_dx = X[0];
			double dx = relativeDelta * x;
			if (dx == 0)
				dx = relativeDelta;
			if (dx > max_dx)
			{
                dx = max_dx / 2;
				if (max_dx == 0)
					dx = -relativeDelta;
            }
			double oldSolventVolumeFraction = X[0];
			//NEW COMPOSITION
			double x_dx = x + dx;
			X[componenIndex] = x_dx;
			X[0] -= dx;
			double f_df = CalculateMixingFreeEnergy(X);

			////////////////RETURN TO OLD//////////////////////
			X[0] = oldSolventVolumeFraction;
			X[componenIndex] = x;
			double output = (f_df - f0) / dx;
			//double analiticalResult = CalculateExchangeChemialPotentialOfComponentAnaliticaly(X, componenIndex);
            return output;


		}
		private double[] CalculateGugenheimCorrelations(double[] alphas, double[,] etas)
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
        private void FindSegregationPointsBetweenSolventAndPolymer()
        {
            double[] initialComposition = new double[Program.NumberOfComponents];

            initialComposition[1] = 0.7; //polymer
            initialComposition[0] = 1 - initialComposition[1];//solvent
			for (int i = 2; i < Program.NumberOfComponents; i++)
				initialComposition[i] = 0;

            double Finit = CalculateMixingFreeEnergy(initialComposition);

            double x1 = initialComposition[1];//polymer molar fraction
            double x2 = initialComposition[1];
            double dx = 0.0001;

            /*double[] X = new double[Program.NumberOfComponents];
            for (int i = 0; i < Program.NumberOfComponents; i++)
                X[i] = 0;*/

            double bestFmix = 10000000;
            double[] best_x = new double[2];
            double[] best_Fedge = new double[2];

            while (x1 > 0)
            {
                x2 = initialComposition[1] + dx;

                double[] leftComposition = new double[3];
                leftComposition[0] = 1 - x1;
                leftComposition[1] = x1;
                for (int i = 2; i < Program.NumberOfComponents; i++)
                    leftComposition[i] = 0;
                double leftF = CalculateMixingFreeEnergy(leftComposition);

                while (x2 < 1)
                {
                    double[] rightComposition = new double[3];
                    rightComposition[0] = 1 - x2;
                    rightComposition[1] = x2;
                    for (int i = 2; i < Program.NumberOfComponents; i++)
                        rightComposition[i] = 0;
                    double rightF = CalculateMixingFreeEnergy(rightComposition);
                    double Fsep = leftF + (initialComposition[1] - x1) * (rightF - leftF) / (x2 - x1);
                    double delta = Fsep - Finit;
                    if (delta < bestFmix)
                    {
                        bestFmix = delta;
                        best_x[0] = x1;
                        best_x[1] = x2;
						best_Fedge[0] = leftF;
						best_Fedge[1] = rightF;
                    }
                    x2 += dx;
                }
                x1 -= dx;
            }
			//return best_x;
			segregationPoints = new double[2];
			segregationMixingEnergies = new double[2];

            for (int i = 0; i < 2; i++)
			{
                segregationPoints[i] = best_x[i];
				segregationMixingEnergies[i] = best_Fedge[i];
            }

        }
    }
}
