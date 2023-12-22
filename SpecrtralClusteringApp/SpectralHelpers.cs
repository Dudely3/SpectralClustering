//namespace SpecrtralClusteringApp;

internal static class SpectralHelpers
{

	// === QR decomposition functions =======================

	public static void MatDecomposeQR(double[][] mat,
	  out double[][] q, out double[][] r,
	  bool standardize)
	{
		// QR decomposition, Householder algorithm.
		// assumes square matrix

		int n = mat.Length;  // assumes mat is nxn
		int nCols = mat[0].Length;
		if (n != nCols)
		{
			Console.WriteLine("M not square ");
		}

		double[][] Q = MatIdentity(n);
		double[][] R = MatCopy(mat);
		for (int i = 0; i < n - 1; ++i)
		{
			double[][] H = MatIdentity(n);
			double[] a = new double[n - i];
			int k = 0;
			for (int ii = i; ii < n; ++ii)  // last part col [i]
			{
				a[k++] = R[ii][i];
			}

			double normA = VecNorm(a);
			if (a[0] < 0.0) { normA = -normA; }
			double[] v = new double[a.Length];
			for (int j = 0; j < v.Length; ++j)
			{
				v[j] = a[j] / (a[0] + normA);
			}

			v[0] = 1.0;

			double[][] h = MatIdentity(a.Length);
			double vvDot = VecDot(v, v);
			double[][] alpha = VecToMat(v, v.Length, 1);
			double[][] beta = VecToMat(v, 1, v.Length);
			double[][] aMultB = MatProduct(alpha, beta);

			for (int ii = 0; ii < h.Length; ++ii)
			{
				for (int jj = 0; jj < h[0].Length; ++jj)
				{
					h[ii][jj] -= (2.0 / vvDot) * aMultB[ii][jj];
				}
			}

			// copy h into lower right of H
			int d = n - h.Length;
			for (int ii = 0; ii < h.Length; ++ii)
			{
				for (int jj = 0; jj < h[0].Length; ++jj)
				{
					H[ii + d][jj + d] = h[ii][jj];
				}
			}

			Q = MatProduct(Q, H);
			R = MatProduct(H, R);
		} // i

		if (standardize == true)
		{
			// standardize so R diagonal is all positive
			double[][] D = MatMake(n, n);
			for (int i = 0; i < n; ++i)
			{
				if (R[i][i] < 0.0)
				{
					D[i][i] = -1.0;
				}
				else
				{
					D[i][i] = 1.0;
				}
			}
			Q = MatProduct(Q, D);
			R = MatProduct(D, R);
		}

		q = Q;
		r = R;
	} // MatDecomposeQR()

	public static double[][] MatProduct(double[][] matA,
	  double[][] matB)
	{
		int aRows = matA.Length;
		int aCols = matA[0].Length;
		int bRows = matB.Length;
		int bCols = matB[0].Length;
		if (aCols != bRows)
		{
			throw new Exception("Non-conformable matrices");
		}

		double[][] result = MatMake(aRows, bCols);

		for (int i = 0; i < aRows; ++i) // each row of A
		{
			for (int j = 0; j < bCols; ++j) // each col of B
			{
				for (int k = 0; k < aCols; ++k)
				{
					result[i][j] += matA[i][k] * matB[k][j];
				}
			}
		}

		return result;
	}

	public static int[] ArgSort(double[] vec)
	{
		int n = vec.Length;
		int[] idxs = new int[n];
		for (int i = 0; i < n; ++i)
		{
			idxs[i] = i;
		}

		Array.Sort(vec, idxs);  // sort idxs based on vec
		return idxs;
	}

	public static double Distance(double[] v1,
	 double[] v2)
	{
		// helper for MakeAffinityRNC()
		int dim = v1.Length;
		double sum = 0.0;
		for (int j = 0; j < dim; ++j)
		{
			sum += (v1[j] - v2[j]) * (v1[j] - v2[j]);
		}

		return Math.Sqrt(sum);
	}

	// === Eigen functions ==================================

	public static void EigenQR(double[][] M,
	  out double[] eigenVals, out double[][] eigenVecs)
	{
		// compute eigenvalues and eigenvectors same time
		// stats.stackexchange.com/questions/20643/finding-
		//   matrix-eigenvectors-using-qr-decomposition

		int n = M.Length;
		double[][] X = MatCopy(M);  // mat must be square
		double[][] pq = MatIdentity(n);
		int maxCt = 10000;

		int ct = 0;
		while (ct < maxCt)
		{
			MatDecomposeQR(X, out double[][] Q, out double[][] R, false);
			pq = MatProduct(pq, Q);
			X = MatProduct(R, Q);  // note order
			++ct;

			if (MatIsUpperTri(X, 1.0e-8) == true)
			{
				break;
			}
		}

		// eigenvalues are diag elements of X
		double[] evals = new double[n];
		for (int i = 0; i < n; ++i)
		{
			evals[i] = X[i][i];
		}

		// eigenvectors are columns of pq
		double[][] evecs = MatCopy(pq);

		eigenVals = evals;
		eigenVecs = evecs;
	}


	public static double[][] MatCopy(double[][] mat)
	{
		int r = mat.Length;
		int c = mat[0].Length;
		double[][] result = MatMake(r, c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				result[i][j] = mat[i][j];
			}
		}

		return result;
	}

	public static double[][] MatDifference(double[][] ma,
	  double[][] mb)
	{
		int r = ma.Length;
		int c = ma[0].Length;
		double[][] result = MatMake(r, c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				result[i][j] = ma[i][j] - mb[i][j];
			}
		}

		return result;
	}

	public static double[][] MatExtractCols(double[][] m,
	  int[] cols)
	{
		int r = m.Length;
		int c = cols.Length;
		double[][] result = MatMake(r, c);

		for (int j = 0; j < cols.Length; ++j)
		{
			for (int i = 0; i < r; ++i)
			{
				result[i][j] = m[i][cols[j]];
			}
		}
		return result;
	}


	private static double[][] MatIdentity(int n)
	{
		double[][] result = MatMake(n, n);
		for (int i = 0; i < n; ++i)
		{
			result[i][i] = 1.0;
		}

		return result;
	}

	private static bool MatIsUpperTri(double[][] mat,
	  double tol)
	{
		int n = mat.Length;
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < i; ++j)
			{  // check lower vals
				if (Math.Abs(mat[i][j]) > tol)
				{
					return false;
				}
			}
		}
		return true;
	}

	// === common ===========================================

	public static double[][] MatMake(int rows, int cols)
	{
		double[][] result = new double[rows][];
		for (int i = 0; i < rows; ++i)
		{
			result[i] = new double[cols];
		}

		return result;
	}

	public static double MyRBF(double[] v1, double[] v2,
	  double gamma)
	{
		// similarity. when v1 == v2, rbf = 1.0
		// less similar returns small values between 0 and 1
		int dim = v1.Length;
		double sum = 0.0;
		for (int i = 0; i < dim; ++i)
		{
			sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
		}

		return Math.Exp(-gamma * sum);
	}

	private static double VecDot(double[] v1,
	  double[] v2)
	{
		double result = 0.0;
		int n = v1.Length;
		for (int i = 0; i < n; ++i)
		{
			result += v1[i] * v2[i];
		}

		return result;
	}


	private static double VecNorm(double[] vec)
	{
		int n = vec.Length;
		double sum = 0.0;
		for (int i = 0; i < n; ++i)
		{
			sum += vec[i] * vec[i];
		}

		return Math.Sqrt(sum);
	}


	private static double[][] VecToMat(double[] vec,
	  int nRows, int nCols)
	{
		double[][] result = MatMake(nRows, nCols);
		int k = 0;
		for (int i = 0; i < nRows; ++i)
		{
			for (int j = 0; j < nCols; ++j)
			{
				result[i][j] = vec[k++];
			}
		}

		return result;
	}
}