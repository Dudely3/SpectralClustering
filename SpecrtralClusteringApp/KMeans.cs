namespace SpectralClusteringApp
{

	public partial class Spectral
	{
		// === nested KMeans ====================================

		private class KMeans
		{
			public double[][] data;
			public int k;
			public int N;
			public int dim;
			public int trials;  // to find best
			public int maxIter; // inner loop
			public int[] clustering;
			public double[][] means;

			public KMeans(double[][] data, int k)
			{
				this.data = data;  // by ref
				this.k = k;
				N = data.Length;
				dim = data[0].Length;
				trials = N;  // for Cluster()
				maxIter = N * 2;  // for ClusterOnce()
				clustering = new int[N];
				means = new double[k][];
				Initialize(0); // seed, means, clustering
			}


			public void Initialize(int seed)
			{
				Random rnd = new(seed);		

				for (int i = 0; i < k; ++i)
				{
					means[i] = new double[dim];
				}

				// Random Partition (not Forgy)
				int[] indices = new int[N];
				for (int i = 0; i < N; ++i)
				{
					indices[i] = i;
				}

				Shuffle(indices, rnd);
				for (int i = 0; i < k; ++i)  // first k items
				{
					clustering[indices[i]] = i;
				}

				for (int i = k; i < N; ++i)
				{					
					clustering[indices[i]] =
					rnd.Next(0, k); // remaining items
				}

				UpdateMeans();
			}


			private void Shuffle(int[] indices, Random rnd)
			{
				int n = indices.Length;
				for (int i = 0; i < n; ++i)
				{
					int r = rnd.Next(i, n);
					(indices[r], indices[i]) = (indices[i], indices[r]);
				}
			}


			private static double SumSquared(double[] v1,
			  double[] v2)
			{
				int dim = v1.Length;
				double sum = 0.0;
				for (int i = 0; i < dim; ++i)
				{
					sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
				}

				return sum;
			}


			private static double Distance(double[] item,
			  double[] mean)
			{
				double ss = SumSquared(item, mean);
				return Math.Sqrt(ss);
			}


			private static int ArgMin(double[] v)
			{
				int dim = v.Length;
				int minIdx = 0;
				double minVal = v[0];
				for (int i = 0; i < v.Length; ++i)
				{
					if (v[i] < minVal)
					{
						minVal = v[i];
						minIdx = i;
					}
				}
				return minIdx;
			}


			private static bool AreEqual(int[] a1, int[] a2)
			{
				int dim = a1.Length;
				for (int i = 0; i < dim; ++i)
				{
					if (a1[i] != a2[i])
					{
						return false;
					}
				}

				return true;
			}


			private static int[] Copy(int[] arr)
			{
				int dim = arr.Length;
				int[] result = new int[dim];
				for (int i = 0; i < dim; ++i)
				{
					result[i] = arr[i];
				}

				return result;
			}


			public bool UpdateMeans()
			{
				// verify no zero-counts
				int[] counts = new int[k];
				for (int i = 0; i < N; ++i)
				{
					int cid = clustering[i];
					++counts[cid];
				}
				for (int kk = 0; kk < k; ++kk)
				{
					if (counts[kk] == 0)
					{
						throw
						  new Exception("0-count in UpdateMeans()");
					}
				}

				// compute proposed new means
				for (int kk = 0; kk < k; ++kk)
				{
					counts[kk] = 0;  // reset
				}

				double[][] newMeans = new double[k][];
				for (int i = 0; i < k; ++i)
				{
					newMeans[i] = new double[dim];
				}

				for (int i = 0; i < N; ++i)
				{
					int cid = clustering[i];
					++counts[cid];
					for (int j = 0; j < dim; ++j)
					{
						newMeans[cid][j] += data[i][j];
					}
				}
				for (int kk = 0; kk < k; ++kk)
				{
					if (counts[kk] == 0)
					{
						return false;  // bad attempt to update
					}
				}

				for (int kk = 0; kk < k; ++kk)
				{
					for (int j = 0; j < dim; ++j)
					{
						newMeans[kk][j] /= counts[kk];
					}
				}

				// copy new means
				for (int kk = 0; kk < k; ++kk)
				{
					for (int j = 0; j < dim; ++j)
					{
						means[kk][j] = newMeans[kk][j];
					}
				}

				return true;
			} // UpdateMeans()


			public bool UpdateClustering()
			{
				// verify no zero-counts
				int[] counts = new int[k];
				for (int i = 0; i < N; ++i)
				{
					int cid = clustering[i];
					++counts[cid];
				}
				for (int kk = 0; kk < k; ++kk)
				{
					if (counts[kk] == 0)
					{
						throw new
						  Exception("0-count in UpdateClustering()");
					}
				}

				// proposed new clustering
				int[] newClustering = new int[N];
				for (int i = 0; i < N; ++i)
				{
					newClustering[i] = clustering[i];
				}

				double[] distances = new double[k];
				for (int i = 0; i < N; ++i)
				{
					for (int kk = 0; kk < k; ++kk)
					{
						distances[kk] =
						  Distance(data[i], means[kk]);
						int newID = ArgMin(distances);
						newClustering[i] = newID;
					}
				}

				if (AreEqual(clustering, newClustering) == true)
				{
					return false;  // no change; short-circuit
				}

				// make sure no count went to 0
				for (int i = 0; i < k; ++i)
				{
					counts[i] = 0;  // reset
				}

				for (int i = 0; i < N; ++i)
				{
					int cid = newClustering[i];
					++counts[cid];
				}
				for (int kk = 0; kk < k; ++kk)
				{
					if (counts[kk] == 0)
					{
						return false;  // bad update attempt
					}
				}

				// no 0 counts so update
				for (int i = 0; i < N; ++i)
				{
					clustering[i] = newClustering[i];
				}

				return true;
			} // UpdateClustering()


			public int[] ClusterOnce()
			{
				int sanityCt = 1;
				while (sanityCt <= maxIter)
				{
					bool ok;
					if ((ok = UpdateClustering() == false))
					{
						break;
					}

					if ((ok = UpdateMeans() == false))
					{
						break;
					}

					++sanityCt;
				}
				return clustering;
			} // ClusterOnce()


			public double WCSS()
			{
				// within-cluster sum of squares
				double sum = 0.0;
				for (int i = 0; i < N; ++i)
				{
					int cid = clustering[i];
					double[] mean = means[cid];
					double ss = SumSquared(data[i], mean);
					sum += ss;
				}
				return sum;
			}


			public int[] Cluster()
			{
				double bestWCSS = WCSS();  // initial clustering
				int[] bestClustering = Copy(clustering);

				for (int i = 0; i < trials; ++i)
				{
					Initialize(i);  // new seed, means, clustering
					int[] clustering = ClusterOnce();
					double wcss = WCSS();
					if (wcss < bestWCSS)
					{
						bestWCSS = wcss;
						bestClustering = Copy(clustering);
					}
				}
				return bestClustering;
			} // Cluster()

		} // class KMeans


	}
} // ns