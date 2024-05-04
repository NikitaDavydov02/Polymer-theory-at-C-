using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Polymer_brush
{
    class MixingPartModule
    {
		bool correlation = true;
		//public double[] segregationPoints { get; private set; }
		//public double[] segregationMixingEnergies;

		//public List<KeyValuePair<double, double[]>> TernarySegregationPoints { get; private set; }
		//public List<KeyValuePair<double, double[]>> TernarySegregatioMixingEnergies { get; private set; }
		public List<Node> Nodes { get; private set; }
		public double mixingPart;
		public double entropyPart;
		public MixingPartModule()
		{
			//TernarySegregationPoints = new List<KeyValuePair<double, double[]>>();
			//TernarySegregatioMixingEnergies = new List<KeyValuePair<double, double[]>>();
			Nodes = new List<Node>();
			Console.WriteLine("Fining segregation points");
			FindSegregationPointsBetweenSolventAndPolymer();
			Console.WriteLine(Nodes.Count + " nodes are found");
		}
		public bool IsCompositionInsideSegregationZone(double[] X)
		{
			double minPossibleF;
			double[] firstSegregatonPoint;
			double[] secondSegregatonPoint;
			return IsCompositionInsideSegregationZone(X, out minPossibleF, out firstSegregatonPoint, out secondSegregatonPoint);
		}
		public bool IsCompositionInsideSegregationZone(double[] X, out double minPossibleF)
		{
			double[] firstSegregatonPoint;
			double[] secondSegregatonPoint;
			return IsCompositionInsideSegregationZone(X, out minPossibleF, out firstSegregatonPoint, out secondSegregatonPoint);
		}
		public bool IsCompositionInsideSegregationZone(double[] X, out double minPossibleF, out double[]firstSegregatonPoint, out double[] secondSegregatonPoint)
        {
			
			minPossibleF = 0;
			firstSegregatonPoint = new double[Program.NumberOfComponents];
			secondSegregatonPoint = new double[Program.NumberOfComponents];
			if (Nodes.Count == 0)
				return false;
			for (int i = 0; i < Program.NumberOfComponents; i++)
            {
				secondSegregatonPoint[i] = -1;
				firstSegregatonPoint[i] = -1;
			}
			double Xadditive;
			if (Program.NumberOfComponents > 2)
            {
				Xadditive = X[2];
				for (int i = 0; i < Nodes.Count - 1; i++)
				{
					if (Nodes[i].firstComposition[2] <= Xadditive && Nodes[i+1].firstComposition[2] > Xadditive)
					{
						//Aproximation between crossections
						double Xadditive_min = Nodes[i].firstComposition[2];
						double Xadditive_max = Nodes[i+1].firstComposition[2];
						double x1_left_min = Nodes[i].firstComposition[1];
						double x1_left_max = Nodes[i+1].firstComposition[1];
						double x1_right_min = Nodes[i].secondComposition[1];
						double x1_right_max = Nodes[i+1].secondComposition[1];
						/////////////
						double x1_left = x1_left_min + (Xadditive - Xadditive_min) * (x1_left_max - x1_left_min) / (Xadditive_max - Xadditive_min);
						double x1_right = x1_right_min + (Xadditive - Xadditive_min) * (x1_right_max - x1_right_min) / (Xadditive_max - Xadditive_min);
						if (X[1] >= x1_left && X[1] < x1_right)
						{
							double[] leftComposition = new double[3] { 1 - x1_left - Xadditive, x1_left, Xadditive };
							double[] rightComposition = new double[3] { 1 - x1_right - Xadditive, x1_right, Xadditive };
							double leftF = CalculateMixingFreeEnergy(leftComposition, false);
							double rightF = CalculateMixingFreeEnergy(rightComposition, false);
							minPossibleF = leftF + (X[1] - x1_left) * (rightF - leftF) / (x1_right - x1_left);
							firstSegregatonPoint = leftComposition;
							secondSegregatonPoint = rightComposition;
							return true;
						}
						return false;
					}
				}
			}
            else
            {
				Xadditive = 0;
				if (X[1] >= Nodes[0].firstComposition[1] && X[1] <= Nodes[0].secondComposition[1])
                {
					double[] leftComposition = Nodes[0].firstComposition;
					double[] rightComposition = Nodes[0].secondComposition;
					double leftF = CalculateMixingFreeEnergy(leftComposition, false);
					double rightF = CalculateMixingFreeEnergy(rightComposition, false);
					minPossibleF = leftF + (X[1] - Nodes[0].firstComposition[1]) * (rightF - leftF) / (Nodes[0].secondComposition[1] - Nodes[0].firstComposition[1]);
					firstSegregatonPoint = leftComposition;
					secondSegregatonPoint = rightComposition;
					return true;
                }
			}
			
			return false;
        }
		public double CalculateMixingFreeEnergy(double[] X, bool withSegregation = true)
        {
            if (withSegregation)
            {
				double segreagationF;
				if (IsCompositionInsideSegregationZone(X, out segreagationF))
					return segreagationF;
			}
			if(!correlation)
				return CalculateFloryMixingFreeEnergy(X);
			else
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
			entropyPart = a;
			//a=X[0] * Math.Log(X[0]) + X[1] * Math.Log(X[1]) / Program.size[1];
			/*double b = Program.chi[1, 2] * X[1] * X[2];
			double c = Program.chi[0, 1] * X[0] * X[1];
			double d = Program.chi[0, 2] * X[0] * X[2];*/
			double b = 0;
			for (int i = 0; i < Program.chiMatrixSize; i++)
				for (int j = 0;j < Program.chiMatrixSize; j++)
					if(j>i && X[i] != 0	&& X[j] != 0)
					{
						if(i==0&&j==1)
                        {
							//b += (Program.chi[i, j] + Program.c * X[1] * X[1]) * X[i] * X[j];
							double x1 = 0;
							//b += (Program.chi[i, j] + x1* X[1] + Program.c * X[1] * X[1]) *X[i] * X[j];
							b += (Program.chi[i, j] + Program.c * X[1]) * X[i] * X[j];
						}
						else
							b += Program.chi[i, j] * X[i] * X[j];
                    }
			//return a + b + c + d;
			mixingPart = b;
			return (a+b);
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
                    {
						if (i == 1 && j == 0)
							mixingSum += (Program.chi[i, j] + Program.c *  X[1]) * XX[i] * XX[j] * X[i] * X[j] * Program.etas[i, j];
						else
							mixingSum += Program.chi[i, j] * XX[i] * XX[j] * X[i] * X[j] * Program.etas[i, j];

					}

			/*double b = chi[1, 2] * X[1] * X[2];
			double c = chi[0, 1] * X[0] * X[1];
			double d = chi[0, 2] * X[0] * X[2];*/
			double output = translationSum + mixingSum;
			entropyPart = translationSum;
			mixingPart = mixingSum;
			return output;
		}
        public double CalculateExchangeChemialPotentialOfComponent(double[] X, int componenIndex, bool smallConc=true)
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
			if (smallConc && dx < Math.Pow(10, -12))
				dx = Math.Pow(10, -12);
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
		private void FindSegregationPointsBetweenSolventAndPolymerAtPresenceOfAdditive(double Xadditive)
        {
			double[] initialComposition = new double[Program.NumberOfComponents];

			initialComposition[1] = 0.5; //polymer
			if(Program.NumberOfComponents>2)
				initialComposition[2] = Xadditive;
			initialComposition[0] = 1 - initialComposition[1]- Xadditive;//solvent

			double Finit = CalculateMixingFreeEnergy(initialComposition,false);

			double x1 = initialComposition[1];//polymer molar fraction
			double x2 = initialComposition[1];
			double dx = 0.0001;


			double bestFmix = 10000000;
			double[] bestFirstComposition = new double[Program.NumberOfComponents];
			double[] bestSecondComposition = new double[Program.NumberOfComponents];
			double[] best_Fedge = new double[2];

			double leftF;
			double rightF;
			double[] leftComposition;
			double[] rightComposition;

			double Fsep;
			double delta;

			bool search = true;

			do
			{
				/*AcceptDeltaForSegregationSearch(x1, x2, Xadditive, out rightF, out leftF, out leftComposition, out rightComposition);
				if (x2 != x1)
					Fsep = leftF + (initialComposition[1] - x1) * (rightF - leftF) / (x2 - x1);
				else
					Fsep = leftF;
				double deltaBeforeStep = Fsep - Finit;*/

				search = false;
                //step right
                if (x2 + dx < 1)
                {
					x2 += dx;
					AcceptDeltaForSegregationSearch(x1, x2, Xadditive, out rightF, out leftF, out leftComposition, out rightComposition);
					Fsep = leftF + (initialComposition[1] - x1) * (rightF - leftF) / (x2 - x1);
					delta = Fsep - Finit;
					if (delta < bestFmix)
					{
						bestFmix = delta;
						bestFirstComposition = leftComposition;
						bestSecondComposition = rightComposition;
						best_Fedge[0] = leftF;
						best_Fedge[1] = rightF;
						search = true;
					}
					else
						x2 -= dx;
				}

                if (x1 - dx > 0)
                {
					//step left
					x1 -= dx;
					AcceptDeltaForSegregationSearch(x1, x2, Xadditive, out rightF, out leftF, out leftComposition, out rightComposition);
					Fsep = leftF + (initialComposition[1] - x1) * (rightF - leftF) / (x2 - x1);
					delta = Fsep - Finit;
					if (delta < bestFmix)
					{
						bestFmix = delta;
						bestFirstComposition = leftComposition;
						bestSecondComposition = rightComposition;
						best_Fedge[0] = leftF;
						best_Fedge[1] = rightF;
						search = true;
					}
					else
						x1 += dx;
				}

                if (x1 - dx > 0 && x2+dx<1)
                {
					//both steps
					x1 -= dx;
					x2 += dx;
					AcceptDeltaForSegregationSearch(x1, x2, Xadditive, out rightF, out leftF, out leftComposition, out rightComposition);
					Fsep = leftF + (initialComposition[1] - x1) * (rightF - leftF) / (x2 - x1);
					delta = Fsep - Finit;
					if (delta < bestFmix)
					{
						bestFmix = delta;
						bestFirstComposition = leftComposition;
						bestSecondComposition = rightComposition;
						best_Fedge[0] = leftF;
						best_Fedge[1] = rightF;
						search = true;
					}
					else
					{
						x1 += dx;
						x2 -= dx;
					}
				}
			}
			while (search);
			
			if (DistanceBetweenCompositions(bestFirstComposition, bestSecondComposition) < 0.01)
				return;
			Console.WriteLine("Node is found");
			Node node = new Node(bestFirstComposition, bestSecondComposition, best_Fedge[0], best_Fedge[1]);
			Nodes.Add(node);
			//TernarySegregationPoints.Add(new KeyValuePair<double, double[]>(Xadditive,new double[2]{ best_x[0],best_x[1]}));
			//TernarySegregatioMixingEnergies.Add(new KeyValuePair<double, double[]>(Xadditive, new double[2] { best_Fedge[0], best_Fedge[1] }));
			
			/*segregationPoints = new double[2];
			segregationMixingEnergies = new double[2];

			for (int i = 0; i < 2; i++)
			{
				segregationPoints[i] = best_x[i];
				segregationMixingEnergies[i] = best_Fedge[i];
			}*/
		}
		public void AcceptDeltaForSegregationSearch(double x1, double x2, double Xadditive, out double rightF,out double leftF, out double[] leftComposition, out double[] rightComposition)
        {
			leftComposition = new double[Program.NumberOfComponents];
			//left composition
			leftComposition[0] = 1 - Xadditive - x1;
			leftComposition[1] = x1;
			if (Program.NumberOfComponents > 2)
				leftComposition[2] = Xadditive;
			leftF = CalculateMixingFreeEnergy(leftComposition, false);
			//right composition
			rightComposition = new double[Program.NumberOfComponents];
			rightComposition[0] = 1 - Xadditive - x2;
			rightComposition[1] = x2;
			if (Program.NumberOfComponents > 2)
				rightComposition[2] = Xadditive;
			rightF = CalculateMixingFreeEnergy(rightComposition, false);
			
		}
		private double DistanceBetweenCompositions(double[]A, double[] B)
        {
			double output = 0;
			for (int i = 0; i < A.Length; i++)
				output += (A[i] - B[i]) * (A[i] - B[i]);
			return output;
		}
        private void FindSegregationPointsBetweenSolventAndPolymer()
        {
			if(Program.NumberOfComponents==2)
				FindSegregationPointsBetweenSolventAndPolymerAtPresenceOfAdditive(0);
            else
            {
				for (double Xadditive = 0; Xadditive < 0.2; Xadditive += 0.01)
				{
					FindSegregationPointsBetweenSolventAndPolymerAtPresenceOfAdditive(Xadditive);
				}
			}
			PrintSegregationPoints();
			return;
        }
		private void PrintSegregationPoints()
        {
			using (StreamWriter sw = new StreamWriter("segregation_points.txt"))
			{
				for(int i = 0; i < Nodes.Count; i++)
                {
					string first = "";
					string second = "";
					for(int j = 0; j < Nodes[i].firstComposition.Length; j++)
                    {
						first += Nodes[i].firstComposition[j] + ";";
						second += Nodes[i].secondComposition[j] + ";";
					}

					sw.WriteLine(first);
					sw.WriteLine(second);
				}
			}
        }
    }
}
