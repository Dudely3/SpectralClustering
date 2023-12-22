
namespace SpectralClusteringApp
{

	public partial class Spectral(int k, double gamma)
	{
		public int k = k;
		public double gamma = gamma;

		public int[] Cluster(double[][] X)
		{
			double[][] A = MakeAffinityRBF(X);  // RBF
			double[][] L = MakeLaplacian(A);    // normalized
			double[][] E = MakeEmbedding(L);    // eigenvectors
			int[] result = ProcessEmbedding(E); // k-means

			Console.WriteLine("\nAffinity: ");
			SpectralClusteringProgram.MatShow(A, 5, 5, 4, 9);

			Console.WriteLine("\nLaplacian: ");
			SpectralClusteringProgram.MatShow(L, 5, 5, 4, 9);

			return result;
		}


		private double[][] MakeAffinityRBF(double[][] X)
		{
			// 1s on diagonal (x1 == x2), towards 0 dissimilar
			int n = X.Length;
			double[][] result = SpectralHelpers.MatMake(n, n);
			for (int i = 0; i < n; ++i)
			{
				for (int j = i; j < n; ++j) // upper
				{
					double rbf = SpectralHelpers.MyRBF(X[i], X[j], gamma);
					result[i][j] = rbf;
					result[j][i] = rbf;
				}
			}
			return result;
		}

		private double[][] MakeAffinityRNC(double[][] X)
		{
			// radius neighbors connectivity
			// 1 if x1 and x2 are close; 0 if not close
			int n = X.Length;
			double[][] result = SpectralHelpers.MatMake(n, n);
			for (int i = 0; i < n; ++i)
			{
				for (int j = i; j < n; ++j) // upper
				{
					double d = SpectralHelpers.Distance(X[i], X[j]);
					if (d < gamma)
					{
						result[i][j] = 1.0;
						result[j][i] = 1.0;
					}
				}
			}
			return result;
		}

		private double[][] MakeLaplacian(double[][] A)
		{
			// unnormalized
			// clear but not very efficient to construct D
			// L = D - A
			// here A is an affinity-style adjaceny matrix
			int n = A.Length;

			double[][] D = SpectralHelpers.MatMake(n, n);  // degree matrix
			for (int i = 0; i < n; ++i)
			{
				double rowSum = 0.0;
				for (int j = 0; j < n; ++j)
				{
					rowSum += A[i][j];
				}

				D[i][i] = rowSum;
			}
			double[][] result = SpectralHelpers.MatDifference(D, A);  // D-A
			return NormalizeLaplacian(result);

			// more efficient, but less clear
			//int n = A.Length;
			//double[] rowSums = new double[n];
			//for (int i = 0; i < n; ++i)
			//{
			//  double rowSum = 0.0;
			//  for (int j = 0; j < n; ++j)
			//    rowSum += A[i][j];
			//  rowSums[i] = rowSum;
			//}

			//double[][] result = MatMake(n, n);
			//for (int i = 0; i < n; ++i)
			//  result[i][i] = rowSums[i];  // degree
			//for (int i = 0; i < n; ++i)
			//  for (int j = 0; j < n; ++j)
			//    if (i == j)
			//      result[i][j] = rowSums[i] - A[i][j];
			//    else
			//      result[i][j] = -A[i][j];
			//return result;
		}

		
		private double[][] NormalizeLaplacian(double[][] L)
		{
			// scipy library csgraph._laplacian technique
			int n = L.Length;
			double[][] result = SpectralHelpers.MatCopy(L);
			for (int i = 0; i < n; ++i)
			{
				result[i][i] = 0.0;  // zap away diagonal
			}

			// sqrt of col sums: in-degree version
			double[] w = new double[n];
			for (int j = 0; j < n; ++j)
			{
				double colSum = 0.0;
				for (int i = 0; i < n; ++i)
				{
					colSum += Math.Abs(result[i][j]);
				}

				w[j] = Math.Sqrt(colSum);
			}

			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					result[i][j] /= (w[j] * w[i]);
				}
			}

			// restore diagonal
			for (int i = 0; i < n; ++i)
			{
				result[i][i] = 1.0;
			}

			return result;
		}

		private double[][] MakeEmbedding(double[][] L)
		{
			// eigenvectors for k-smallest eigenvalues
			// extremely deep graph theory
			SpectralHelpers.EigenQR(L, out double[] eigVals, out double[][] eigVecs); // QR algorithm
			int[] allIndices = SpectralHelpers.ArgSort(eigVals);
			int[] indices = new int[k]; // small eigenvecs
			for (int i = 0; i < k; ++i)
			{
				indices[i] = allIndices[i];
			}

			double[][] extracted = SpectralHelpers.MatExtractCols(eigVecs, indices);
			return extracted;
		}

		private int[] ProcessEmbedding(double[][] E)
		{
			// cluster a complex transformation of source data
			KMeans km = new(E, k);
			int[] clustering = km.Cluster();
			return clustering;
		}

	} // Spectral

} // ns