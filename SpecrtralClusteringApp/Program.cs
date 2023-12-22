
namespace SpectralClusteringApp
{
	internal class SpectralClusteringProgram
	{
		static void Main(string[] args)
		{
			Console.WriteLine("\nBegin spectral clustering" +
			  " using C# ");

			double[][] X =
			[
				[0.20, 0.20],
				[0.20, 0.70],
				[0.40, 0.90],
				[0.30, 0.50],
				[0.50, 0.60],
				[0.60, 0.50],
				[0.70, 0.90],
				[0.40, 0.10],
				[0.90, 0.70],
				[0.80, 0.20],
				[0.10, 0.50],
				[0.20, 0.40],
				[0.60, 0.10],
				[0.70, 0.10],
				[0.90, 0.40],
				[0.90, 0.50],
				[0.80, 0.80],
				[0.50, 0.90],
			];

			//string fn = "..\\..\\..\\Data\\dummy_data_18.txt";
			//double[][] X = MatLoad(fn, new int[] { 0, 1 },
			//  ',', "#");

			Console.WriteLine("\nData: ");
			MatShow(X, 18, 2, 2, 6);

			int k = 2;
			// double eps = 0.30; // radius neighbors affinity
			double gamma = 50.0;  // RBF affinity

			Console.WriteLine("\nClustering with k = " + k +
			  " RBF gamma = " + gamma.ToString("F2"));

			Spectral sp = new(k, gamma);
			int[] clustering = sp.Cluster(X);

			Console.WriteLine("\nDone");

			Console.WriteLine("\nspectral clustering: ");

			VecShow(clustering, 2);

			Console.WriteLine("\nEnd demo ");

			Console.ReadLine();
		} // Main()


		public static void MatShow(double[][] m, int nRows,
		  int nCols, int dec, int wid)
		{
			double small = 1.0 / Math.Pow(10, dec);
			for (int i = 0; i < nRows; ++i)
			{
				for (int j = 0; j < nCols; ++j)
				{
					double v = m[i][j];
					if (Math.Abs(v) < small)
					{
						v = 0.0;
					}

					Console.Write(v.ToString("F" + dec).
					  PadLeft(wid));
				}
				if (nCols < m[0].Length)
				{
					Console.Write(" . . .");
				}

				Console.WriteLine("");
			}
			if (nRows < m.Length)
			{
				Console.WriteLine(". . . ".PadLeft(wid));
			}
		}


		static void VecShow(int[] vec, int wid)
		{
			for (int i = 0; i < vec.Length; ++i)
			{
				Console.Write(vec[i].ToString().
				  PadLeft(wid));
			}

			Console.WriteLine("");
		}

		static double[][] MatLoad(string fn, int[] usecols,
		  char sep, string comment)
		{
			// self-contained - no dependencies
			int nRows = 0;
			FileStream ifs = new(fn, FileMode.Open);
			StreamReader sr = new(ifs);
			string? line;
			while ((line = sr.ReadLine()) != null)
			{
				if (line.StartsWith(comment) == false)
				{
					++nRows;
				}
			}

			sr.Close(); ifs.Close();  // could reset instead

			int nCols = usecols.Length;
			double[][] result = new double[nRows][];
			for (int r = 0; r < nRows; ++r)
			{
				result[r] = new double[nCols];
			}

			ifs = new FileStream(fn, FileMode.Open);
			sr = new StreamReader(ifs);

			int i = 0;
			while ((line = sr.ReadLine()) != null)
			{
				if (line.StartsWith(comment) == true)
				{
					continue;
				}

				string[] tokens = line.Split(sep);
				for (int j = 0; j < nCols; ++j)
				{
					int k = usecols[j];  // into tokens
					result[i][j] = double.Parse(tokens[k]);
				}
				++i;
			}
			sr.Close();
			ifs.Close();
			return result;
		}

	} // Program

} // ns