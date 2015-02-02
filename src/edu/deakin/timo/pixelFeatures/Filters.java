/*Bicubic interpolation copied from imageJ imageProcessor.*/

/*Some filtering functions reproduced here to enable using the code without ImageJ
	2D arrays, first pointer x (i.e. width), second pointer y (i.e. height): data[x][y]

*/

package edu.deakin.timo.pixelFeatures;

import java.util.ArrayList;
import java.text.DecimalFormat;	/*For debugging*/
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;

public class Filters{

	/*
		Get edge image as the resultant of differentials in horizontal, and diagonal directions implemented with 2D convolutions. Similar to Sobel, but with smaller kernel
		@param data the image for edge highlighting
		@return edgeImage the resultant of differentials
	*/
	public static double[][] getEdgeImage(double[][] data){
		/*
		double[][][] convKerns = {
									{{-1,-1},{1,1}},
									{{-1,1},{-1,1}},
									{{-1,0},{0,1}},
									{{0,-1},{1,0}}
									};
		*/
		//~Sobel
		final double[][][] convKerns = {
									{{-1,-2,-1},{0,0,0},{1,2,1}},
									{{-1,0,1},{-2,0,2},{-1,0,1}}
									/*
									{{-2,-1,0},{-1,0,1},{0,1,2}},
									{{0,-1,-2},{1,0,-1},{2,1,0}}
									*/
								};
		double[][][] tempConvs = new double[convKerns.length][data.length][data[0].length];
		for (int i = 0; i< convKerns.length;++i){
			tempConvs[i] =  conv2(data, convKerns[i]);
		}
		double[][] edgeImage = new double[data.length][data[0].length];
		for (int c = 0;c<data.length;++c){
			for (int r = 0;r<data[0].length;++r){			
			   edgeImage[c][r] = 0;
			   	for (int d = 0; d<tempConvs.length;++d){
			   	   edgeImage[c][r]+=tempConvs[d][c][r]*tempConvs[d][c][r];
			   	}
				edgeImage[c][r] = Math.sqrt(edgeImage[c][r]);
			}
		}
		return edgeImage;
	}


	/*
		Calculates the 2D convolution
		@param M 2D image matrix
		@param k 2D convolution kernel to be convolved with the image.
		@return convolution convoluted image of the same size as M

	*/
	public static double[][] conv2(double[][] M, double[][] k){
		int[][] convLims = new int[2][2];
		convLims[0][0] = (k.length%2 != 0) ? -(k.length-1)/2 : -k.length/2;
		convLims[1][0] = (k.length%2 != 0) ? (k.length-1)/2+1 : k.length/2;
		convLims[0][1] = (k[0].length%2 != 0) ? -(k[0].length-1)/2 : -k[0].length/2;
		convLims[1][1] = (k[0].length%2 != 0) ? (k[0].length-1)/2+1 : k[0].length/2;
		
	   	double[][] convolution = new double[M.length][M[0].length];
		for (int i =0; i<M.length;++i){
		 	for (int j = 0; j<M[i].length;++j){
		 	   	//Calc the convolution for a given pixel
		 	   	convolution[i][j] = 0;
				for (int ii =convLims[0][0]; ii<convLims[1][0];++ii){
	 				for (int jj = convLims[0][1]; jj<convLims[1][1];++jj){
	 				   if (i+ii >= 0 && i+ii<M.length &&
	 				   	j+jj >= 0 && j+jj<M[i+ii].length
	 				   	){
	 				    	  convolution[i][j]+=M[i+ii][j+jj]*k[ii-convLims[0][0]][jj-convLims[0][1]];
	 				   }

	 				}
				}
		 	}
		}
		return convolution;
	}

	public static double[][] highlightEdges(double[][] imagePixels){
		int rows = imagePixels.length;
		int columns = imagePixels[0].length;
		double[][] gradientrows = new double[rows][columns];
		double[][] gradientcolumns = new double[rows][columns];
		double[][] gradientr = new double[rows][columns];
		//Using sobel
		//for gx convolutes the following matrix
		//
		//     |-1 0 1|
		//Gx = |-2 0 2|
		//     |-1 0 1|
		for(int i=1;i<rows-1;++i){
			for(int j=1;j<columns-1;++j){
				gradientrows[i][j] =
				-1*(imagePixels[i-1][j-1]) +1*(imagePixels[i+1][j-1])
				-2*(imagePixels[i-1][j]) +2*(imagePixels[i+1][j])
				-1*(imagePixels[i-1][j+1]) +1*(imagePixels[i+1][j+1]);
			}
		}

		//for gy convolutes the following matrix
		//
		//     |-1 -2 -1|
		//Gy = | 0  0  0|
		//     |+1 +2 +1|
		//
		for(int i=1;i<rows-1;++i){
			for(int j=1;j<columns-1;++j){
				gradientcolumns[i][j] =
				-1*(imagePixels[i-1][j-1]) +1*(imagePixels[i-1][j+1])
				-2*(imagePixels[i][j-1]) +2*(imagePixels[i][j+1])
				-1*(imagePixels[i+1][j-1]) +1*(imagePixels[i+1][j+1]);
			}
		}
		for(int i=1;i<rows-1;i++){
			for(int j=1;j<columns-1;j++){
				gradientr[i][j] = Math.sqrt(gradientrows[i][j]*gradientrows[i][j]+gradientcolumns[i][j]*gradientcolumns[i][j]);
			}
		}
		return gradientr;
    }


	public static double[][] getGradientImage(double[][] imagePixels){
		int rows = imagePixels.length;
		int columns = imagePixels[0].length;
		double[][] gradientrows = new double[rows][columns];
		double[][] gradientcolumns = new double[rows][columns];
		double[][] gradientr = new double[rows][columns];
		//Using sobel
		//for gx convolutes the following matrix
		//
		//     |-1 0 1|
		//Gx = |-2 0 2|
		//     |-1 0 1|
		for(int i=1;i<rows-1;++i){
			for(int j=1;j<columns-1;++j){
				gradientrows[i][j] =
				-1*(imagePixels[i-1][j-1]) +1*(imagePixels[i+1][j-1])
				-2*(imagePixels[i-1][j]) +2*(imagePixels[i+1][j])
				-1*(imagePixels[i-1][j+1]) +1*(imagePixels[i+1][j+1]);
			}
		}

		//for gy convolutes the following matrix
		//
		//     |-1 -2 -1|
		//Gy = | 0  0  0|
		//     |+1 +2 +1|
		//
		for(int i=1;i<rows-1;++i){
			for(int j=1;j<columns-1;++j){
				gradientcolumns[i][j] =
				-1*(imagePixels[i-1][j-1]) +1*(imagePixels[i-1][j+1])
				-2*(imagePixels[i][j-1]) +2*(imagePixels[i][j+1])
				-1*(imagePixels[i+1][j-1]) +1*(imagePixels[i+1][j+1]);
			}
		}
		for(int i=1;i<rows-1;i++){
			for(int j=1;j<columns-1;j++){
				gradientr[i][j] = Math.sqrt(gradientrows[i][j]*gradientrows[i][j]+gradientcolumns[i][j]*gradientcolumns[i][j]);
			}
		}
		return gradientr;
    }

	public static double[][][] get3DVarianceStack(ImagePlus imp){
		int width = imp.getWidth();
		int height = imp.getHeight();
		int depth = imp.getStackSize();
		double[][][] data = new double[width][height][depth];
		final Object[] imageArrayPointers = imp.getStack().getImageArray();
		short[] temp;
       // IJ.log("Get data get3dVarianceStack");
		for (int d = 0; d < depth; ++d) {
            temp = (short[]) imageArrayPointers[d];
			for (int r = 0;r<height;++r){
				for (int c = 0;c<width;++c){
					data[c][r][d] = (double) temp[c+r*width];
				}
			}
        }
		//IJ.log("Got data get3dVarianceStack");
		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<VarianceRunnable> runnables = new ArrayList<VarianceRunnable>();
		for (int d = 0; d < depth; ++d) {
			runnables.add(new VarianceRunnable(data,d));
			threads.add(new Thread(runnables.get(runnables.size()-1)));
			threads.get(threads.size()-1).start();
		}
		//IJ.log("Started threads");
		/**Catch threads*/
		double[][][] varianceStack = new double[width][height][depth];
		for (int d = 0; d < depth; ++d) {
			try{
				threads.get(d).join();
				//IJ.log("Joined thread "+d);
			}catch (Exception err){IJ.log("Couldn't catch variance image thread");}
			double[][] vIm = runnables.get(d).varianceIm;
			//double varSum = 0;
			for (int r = 0;r<height;++r){
				for (int c = 0;c<width;++c){
					varianceStack[c][r][d] = vIm[c][r];
					//varSum+=vIm[c][r];
				}
			}
			
		}
		return varianceStack;
	}
	
	/**Helper class to enable multi-threading*/
	public static class VarianceRunnable implements Runnable{
		
		private double[][][] data;
		public double[][] varianceIm;
		public int d;
		public VarianceRunnable(double[][][] data,int d){
			this.data = data;
			this.d = d;
		}
		public void run(){
			int width = data.length;
			int height = data[0].length;
			int depth = data[0][0].length;
			int x;
			int y;
			int[][] neighbourhood = {	{-1,-1},{-1,0},{-1,1},
										{0,-1},{0,0},{0,1},
										{1,-1},{1,0},{1,1}
										};
			
			//IJ.log("Variance "+d+" nx "+ neighbourhood.length+" yx "+ neighbourhood[0].length);
			varianceIm = new double[width][height];
			for (int i = 0; i<width;++i){
				for (int j = 0; j<height;++j){
					ArrayList<Double> tempData = new ArrayList<Double>();
					/**Add the local neighbourhood pixels to tempData*/
					for (int z = d-1;z<=d+1;++z){
						for (int n = 0; n<neighbourhood.length;++n){
							x = i+neighbourhood[n][0];
							y = j+neighbourhood[n][1];
							if (x >=0 && x < width &&
								y >=0 && y < height &&
								z >=0 && z < depth){
								tempData.add(data[x][y][z]);
							}
						}
					}
					/*calculate local variance based on pixels in the neighbourhood*/
					//IJ.log("Variance "+d+" width "+ width+" height "+ height+" i "+i+" j "+j+" tempData "+tempData.size());
					//IJ.log("Variance "+d+" nx "+ neighbourhood.length+" yx "+ neighbourhood[0].length+" tempData "+tempData.size());
					varianceIm[i][j] = Filters.getVariance(tempData.toArray(new Double[tempData.size()]));
				}
			}
			
		}
	}
	
	
	public static double[][] getVarianceImage(double[][] data, int radius){
		int width = data.length;
		int height = data[0].length;
		double[][] varianceImage = new double[width][height];
		double[] coordinates = new double[2];
		for (int i = 0+radius;i<width-radius;++i){
			for (int j = 0+radius;j<height-radius;++j){
				coordinates[0] = i;
				coordinates[1] = j;
				//System.out.println("source x "+coordinates[0]+" y "+coordinates[1]);
				varianceImage[i][j] = getLocalVariance(data,coordinates,radius);
			}
		}
		return varianceImage;
	}

	/*Local variance with circular sampling. Eight samples per integer increment of radius*/
	public static double getLocalVariance(double[][] data,double[] coordinates,int radius){
		/*Init sampling coordinates*/
		double[][] samplingCoordinates = new double[8*radius+1][2];
		samplingCoordinates[8*radius] = coordinates;
		final double sqrt05 = Math.sqrt(0.5);
		final double[][] directions = {{1,0},{sqrt05,sqrt05},{0,1},{-sqrt05,sqrt05},{-1,0},{-sqrt05,-sqrt05},{0,-1},{sqrt05,-sqrt05}};
		for (int r=0;r<radius;++r){
			for (int t = 0;t <8; ++t){
				samplingCoordinates[t+(r*8)][0] = coordinates[0]+directions[t][0]*((double)(r+1));
				samplingCoordinates[t+(r*8)][1] = coordinates[1]+directions[t][1]*((double)(r+1));

			}
		}
		/*Get the values*/
		double[] values = new double[8*radius+1];
		//DecimalFormat f = new DecimalFormat("0.#");
		for (int i = 0; i<samplingCoordinates.length;++i){
			values[i] = getBicubicInterpolatedPixel(samplingCoordinates[i][0],samplingCoordinates[i][1],data);
			//System.out.println("\tsampling x\t"+f.format(samplingCoordinates[i][0])+"\ty\t"+f.format(samplingCoordinates[i][1])+"\tval\t"+f.format(values[i]));
		}
		return getVariance(values);
	}

	public static double getMean(double[] data){
		double sum = 0;
		for (int i = 0; i<data.length; ++i){
			sum+= data[i];
		}
		sum/=((double) data.length);
		return sum;
	}
	
	public static double mean(Double[] data){
		return getMean(data);
	}
	
	public static double mean(Integer[] data){
		double[] temp = new double[data.length];
		for (int i = 0; i<data.length;++i){
			temp[i] = (double) data[i];
		}
		return getMean(temp);
	}
	
	
	public static double getMean(Double[] data){
		double sum = 0;
		for (int i = 0; i<data.length; ++i){
			sum+= data[i];
		}
		sum/=((double) data.length);
		return sum;
	}
	
	public static double getVariance(Double[] data){
		double variance = 0;
		double mean = getMean(data);
		for (int i = 0; i<data.length; ++i){
			variance+= Math.pow(data[i]-mean,2.0);
		}
		variance/=((double) data.length);
		return variance;
	}
	
	public static double getVariance(double[] data){
		double variance = 0;
		double mean = getMean(data);
		for (int i = 0; i<data.length; ++i){
			variance+= Math.pow(data[i]-mean,2.0);
		}
		variance/=((double) data.length);
		return variance;
	}
	
	public static double getBilinearInt(double x,double y, double[][] data){
		double u = Math.floor(x);	//use floor to handle negative coordinates too
		double v = Math.floor(y);
		int width = data.length;
		int height = data[0].length;
		if (u < width-1 && u >= 0 && v < height-1 && v >= 0){
			return	 data[(int)u][(int)v]*(1d-x)*(1d-y)
					+data[(int)u+1][(int)v]*x*(1d-y)
					+data[(int)u][(int)v+1]*(1d-x)*y
					+data[(int)u+1][(int)v+1]*x*y;
		}else{
			/*Clamp the value to within the image*/
			u = min(u,width-1);
			u = max(u,0);
			v = min(v,height-1);
			v = max(v,0);
			return data[(int) u][(int) v];
		}
	}
	
	/** This method is from Chapter 16 of "Digital Image Processing:
		An Algorithmic Introduction Using Java" by Burger and Burge
		(http://www.imagingbook.com/). */
	public static double getBicubicInterpolatedPixel(double x0, double y0, double[][] data) {
		int u0 = (int) Math.floor(x0);	//use floor to handle negative coordinates too
		int v0 = (int) Math.floor(y0);
		int width = data.length;
		int height = data[0].length;
		if (u0<1 || u0>width-3 || v0< 1 || v0>height-3){
			if ((u0 == 0 || u0 < width-1) && (v0 == 0 || v0 < height-1)){ /*Use bilinear interpolation http://en.wikipedia.org/wiki/Bilinear_interpolation*/
				double x = (x0-(double)u0);
				double y = (y0-(double)v0);
				return data[u0][v0]*(1-x)*(1-y) 	/*f(0,0)(1-x)(1-y)*/
						+data[u0+1][v0]*(1-y)*x	/*f(1,0)x(1-y)*/
						+data[u0][v0+1]*(1-x)*y	/*f(0,1)(1-x)y*/
						+data[u0+1][v0+1]*x*y;	/*f(1,1)xy*/
			}
			
			
			
			return getBilinearInt(x0,y0,data); /*Return linear interpolatoin for points outside the interpolable area*/
		}
		double q = 0;
		for (int j = 0; j < 4; ++j) {
			int v = v0 - 1 + j;
			double p = 0;
			for (int i = 0; i < 4; ++i) {
				int u = u0 - 1 + i;
				p = p + data[u][v] * cubic(x0 - u);
			}
			q = q + p * cubic(y0 - v);
		}
		return q;
	}

	/*Min, max and mean*/
	public static int min(int a, int b){return a < b ? a:b;}
	public static int max(int a, int b){return a > b ? a:b;}
	public static double min(double a, double b){return a < b ? a:b;}
	public static double max(double a, double b){return a > b ? a:b;}
	public static float max(float a, float b){return a > b ? a:b;}
	
	public static double max(double[] a){
		double t = Double.NEGATIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			t = max(a[i],t);
		}
		return t;
	}

	public static double max(double[][] a){
		double t = Double.NEGATIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				t = max(a[i][j],t);
			}
		}
		return t;
	}
	
	public static double max(float[] a){
		float t = Float.NEGATIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			t = max(a[i],t);
		}
		return t;
	}
	
	public static double max(float[][] a){
		float t = Float.NEGATIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				t = max(a[i][j],t);
			}
		}
		return t;
	}
	
	public static double max(double[][][] a){
		double t = Double.NEGATIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				for (int k = 0;k<a[i][j].length;++k){
					t = max(a[i][j][k],t);
				}
			}
		}
		return t;
	}
	
	public static double min(double[] a){
		double t = Double.POSITIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			t = min(a[i],t);
		}
		return t;
	}
	
	public static double min(Double[] a){
		double t = Double.POSITIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			t = min(a[i],t);
		}
		return t;
	}
	
	
	public static double min(double[][] a){
		double t = Double.POSITIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				t = min(a[i][j],t);
			}
		}
		return t;
	}
	
	public static double min(double[][][] a){
		double t = Double.POSITIVE_INFINITY;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				for (int k = 0;k<a[i][j].length;++k){
					t = min(a[i][j][k],t);
				}
			}
		}
		return t;
	}
	
	/**Standard deviation
		This method of calculating sd explained in e.g. http://www.johndcook.com/blog/2008/09/26/comparing-three-methods-of-computing-standard-deviation/
		
		@param	arr	array for which the sd is to be calculated
		@return standard deviation
	*/
	public static double std(double[] arr){
		double sum = 0,sum2 = 0;
		for (int i = 0; i<arr.length;++i){
			sum+=arr[i];
			sum2+=arr[i]*arr[i];
		}
		return Math.sqrt((sum2-sum*sum/((double) arr.length))/(((double) arr.length)-1d));		
	}
	

	public static double[][][] getNormalized(double[][][] stack){
		int width = stack.length;
		int height = stack[0].length;
		int depth = stack[0][0].length;
		double[] norm = new double[]{min(stack),max(stack)};
		
		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<NormRunnable> runnables = new ArrayList<NormRunnable>();
		for (int d = 0; d < depth; ++d) {
			runnables.add(new NormRunnable(stack,norm,d));
			threads.add(new Thread(runnables.get(runnables.size()-1)));
			threads.get(threads.size()-1).start();
		}
		/**Catch threads*/
		double[][][] normalized = new double[width][height][depth];
		for (int d = 0; d < depth; ++d) {
			try{
				threads.get(d).join();
			}catch (Exception err){IJ.log("Couldn't catch variance image thread");}
			double[][] nIm = runnables.get(d).normIm;
			for (int r = 0;r<height;++r){
				for (int c = 0;c<width;++c){
					normalized[c][r][d] = nIm[c][r];
				}
			}
			
		}
		return normalized;
	}
	
	/**Helper class for normalizing to enable multi-threading*/
	public static class NormRunnable implements Runnable{
		
		private double[][][] data;
		public double[][] normIm;
		public double[] norm;
		public int d;
		public NormRunnable(double[][][] data,double[] norm,int d){
			this.data = data;
			this.norm = norm;
			this.d = d;
		}
		public void run(){
			normIm = new double[data.length][data[0].length];
			for (int r = 0;r<data[0].length;++r){
				for (int c = 0;c<data.length;++c){
					normIm[c][r] = (data[c][r][d]-norm[0])/(norm[1]-norm[0]);
				}
			}
		}
	}
	
	/**2D min image*/
	public static double[][] getMinImage(double[][] data){
		int width = data.length;
		int height = data[0].length;
		int x;
		int y;
		int[][] neighbourhood = {	{-1,-1},{-1,0},{-1,1},
									{0,-1},{0,0},{0,1},
									{1,-1},{1,0},{1,1}
									};
		double[][] minIm = new double[width][height];
		ArrayList<Double> tempData = new ArrayList<Double>();
		for (int i = 0; i<width;++i){
			for (int j = 0; j<height;++j){
				tempData.clear();
				/**Add the local neighbourhood pixels to tempData*/
				for (int n = 0; n<neighbourhood.length;++n){
					x = i+neighbourhood[n][0];
					y = j+neighbourhood[n][1];
					if (x >=0 && x < width &&
						y >=0 && y < height){
						tempData.add(data[x][y]);
					}
				}

				/*get the min from the neighbourhood*/
				minIm[i][j] = Filters.min(tempData.toArray(new Double[tempData.size()]));
			}
		}
		return minIm;
	}
	
	/**3D min stack*/
	public static double[][][] get3DMinStack(double[][][] stack){
		int width = stack.length;
		int height = stack[0].length;
		int depth = stack[0][0].length;

		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<MinRunnable> runnables = new ArrayList<MinRunnable>();
		for (int d = 0; d < depth; ++d) {
			runnables.add(new MinRunnable(stack,d));
			threads.add(new Thread(runnables.get(runnables.size()-1)));
			threads.get(threads.size()-1).start();
		}
		/**Catch threads*/
		double[][][] minStack = new double[width][height][depth];
		for (int d = 0; d < depth; ++d) {
			try{
				threads.get(d).join();
			}catch (Exception err){IJ.log("Couldn't catch variance image thread");}
			double[][] nIm = runnables.get(d).minIm;
			for (int r = 0;r<height;++r){
				for (int c = 0;c<width;++c){
					minStack[c][r][d] = nIm[c][r];
				}
			}
			
		}
		return minStack;
	}
	
	/**Helper class for normalizing to enable multi-threading*/
	public static class MinRunnable implements Runnable{
		
		private double[][][] data;
		public double[][] minIm;
		public int d;
		public MinRunnable(double[][][] data,int d){
			this.data = data;
			this.d = d;
		}
		public void run(){
			int width = data.length;
			int height = data[0].length;
			int depth = data[0][0].length;
			int x;
			int y;
			int[][] neighbourhood = {	{-1,-1},{-1,0},{-1,1},
										{0,-1},{0,0},{0,1},
										{1,-1},{1,0},{1,1}
										};
			minIm = new double[width][height];
			for (int i = 0; i<width;++i){
				for (int j = 0; j<height;++j){
					ArrayList<Double> tempData = new ArrayList<Double>();
					/**Add the local neighbourhood pixels to tempData*/
					for (int z = d-1;z<=d+1;++z){
						for (int n = 0; n<neighbourhood.length;++n){
							x = i+neighbourhood[n][0];
							y = j+neighbourhood[n][1];
							if (x >=0 && x < width &&
								y >=0 && y < height &&
								z >=0 && z < depth){
								tempData.add(data[x][y][z]);
							}
						}
					}
					/*get the min from the neighbourhood*/
					minIm[i][j] = Filters.min(tempData.toArray(new Double[tempData.size()]));
				}
			}
		}
	}	
	
	
	
	
	public static double[] getNormalized(double[] stack){
		double[] normalized = new double[stack.length];
		double[] norm = new double[]{min(stack),max(stack)};
		for (int d = 0; d < stack.length; ++d) {
			normalized[d] = (stack[d]-norm[0])/(norm[1]-norm[0]);
        }
		return normalized;
	}
	
	
	public static double[][] getNormalized(double[][] stack){
		double[][] normalized = new double[stack.length][stack[0].length];
		double[] norm = new double[]{min(stack),max(stack)};
		for (int d = 0; d < stack.length; ++d) {
			for (int r = 0;r<stack[d].length;++r){
				normalized[d][r] = (stack[d][r]-norm[0])/(norm[1]-norm[0]);
			}
        }
		return normalized;
	}
	
	/*1D mean*/
	public static double mean(int[] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			returnVal+=(double) a[i];
		}
		return returnVal/=(double) a.length;
	}
	public static double mean(double[] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			returnVal+= a[i];
		}
		return returnVal/=(double) a.length;
	}
	/*2D mean*/
	public static double mean(int[][] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				returnVal+= (double) a[i][j];
			}
		}
		return returnVal/=((double) (a.length*a[0].length));
	}
	
	
	public static double mean(double[][] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				returnVal+= a[i][j];
			}
		}
		return returnVal/=((double) (a.length*a[0].length));
	}
	/*3D mean*/
	public static double mean(int[][][] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				for (int k = 0;k<a[i][j].length;++k){
					returnVal+= (double) a[i][j][k];
				}
			}
		}
		return returnVal/=((double) (a.length*a[0].length));
	}
	public static double mean(double[][][] a){
		double returnVal = 0;
		for (int i = 0;i<a.length;++i){
			for (int j = 0;j<a[i].length;++j){
				for (int k = 0;k<a[i][j].length;++k){
					returnVal+= a[i][j][k];
				}
			}
		}
		return returnVal/=((double) (a.length*a[0].length));
	}


	/*Cross-correlation analysis, equations taken from http://paulbourke.net/miscellaneous/correlate/*/
	/*Calculate 1D cross-correlation for two arrays with same length. No wrapping, i.e. correlation length is limited*/
	public static double[] xcorr(double[] series1,double[] series2, int maxDelay){
		double[] xcor = new double[maxDelay*2+1];
		double ms1 =0;
		double ms2 =0;
		int length = min(series1.length,series2.length);
		/*calculate means*/
		ms1 = mean(series1);
		ms2 = mean(series2);
		double mx;
		double my;
		double summxmy;
		double summxSq;
		double summySq;
		double summxmySq;
		for (int i =-maxDelay;i<=maxDelay;i++){//ignore beginning and end of the signal...
			summxmy=0;
			summxSq=0;
			summySq=0;
			for (int j = maxDelay; j< length-maxDelay; j++){
				mx = series1[j]-ms1;
				my = series2[j+i]-ms2;
				summxmy+=mx*my;
				summxSq+=mx*mx;
				summySq+=my*my;
			}
			xcor[i+maxDelay]=summxmy/Math.sqrt(summxSq*summySq);
		}
		return xcor;
	}



	/*
		Calculate 2D cross-correlation for two 2D arrays. matrix2 needs to be smaller in both dimensions. Not calculated for non-overlapping positions.
	*/
	public static double[][] xcorr(double[][] matrix1,double[][] matrix2){
		double[][] xcor = new double[matrix1.length-matrix2.length+1][matrix1[0].length-matrix2[0].length+1];
		double ms1 =0;
		double ms2 =0;
		int width = matrix1.length;
		int height = matrix1[0].length;
		/*calculate means*/
		ms1 = mean(matrix1);
		ms2 = mean(matrix2);
		double mx;
		double my;
		double summxmy;
		double summxSq;
		double summySq;
		double summxmySq;
		for (int i =0;i<=width-matrix2.length;++i){
			for (int j =0;j<=height-matrix2[0].length;++j){//ignore beginning and end of the signal...
				summxmy=0;
				summxSq=0;
				summySq=0;
				for (int i2 = 0; i2< matrix2.length; ++i2){
					for (int j2 = 0; j2< matrix2[i2].length; ++j2){
						mx = matrix1[i+i2][j+j2]-ms1;
						my = matrix2[i2][j2]-ms2;
						summxmy+=mx*my;
						summxSq+=mx*mx;
						summySq+=my*my;
					}
				}
				xcor[i][j]=summxmy/((Math.sqrt(summxSq))*(Math.sqrt(summySq)));
			}
		}
		return xcor;
	}

	/*
		Calculate 3D cross-correlation for two 3D arrays. matrix2 needs to be smaller in all dimensions. Not calculated for non-overlapping positions.
	*/
	public static double[][][] xcorr(double[][][] matrix1,double[][][] matrix2){
		double[][][] xcor = new double[matrix1.length-matrix2.length+1][matrix1[0].length-matrix2[0].length+1][matrix1[0][0].length-matrix2[0][0].length+1];
		double ms1 =0;
		double ms2 =0;
		int width = matrix1.length;
		int height = matrix1[0].length;
		int depth = matrix1[0][0].length;
		/*calculate means*/
		//System.out.println("Calc means");
		ms1 = mean(matrix1);
		ms2 = mean(matrix2);
		double mx;
		double my;
		double summxmy;
		double summxSq;
		double summySq;
		double summxmySq;
		for (int i =0;i<=width-matrix2.length;++i){
			for (int j =0;j<=height-matrix2[0].length;++j){//ignore beginning and end of the signal...
				for (int k =0;k<=depth-matrix2[0][0].length;++k){//ignore beginning and end of the signal...
					summxmy=0;
					summxSq=0;
					summySq=0;
					for (int i2 = 0; i2< matrix2.length; ++i2){
						for (int j2 = 0; j2< matrix2[i2].length; ++j2){
							for (int k2 = 0; k2< matrix2[i2][j2].length; ++k2){
								mx = matrix1[i+i2][j+j2][k+k2]-ms1;
								my = matrix2[i2][j2][k2]-ms2;
								summxmy+=mx*my;
								summxSq+=mx*mx;
								summySq+=my*my;
							}
						}
					}
					xcor[i][j][k]=summxmy/((Math.sqrt(summxSq))*(Math.sqrt(summySq)));
				}
			}
		}
		return xcor;
	}

	public static final double cubic(double x) {
		final double a = 0.5; // Catmull-Rom interpolation
		if (x < 0.0) x = -x;
		double z = 0.0;
		if (x < 1.0)
			z = x*x*(x*(-a+2.0) + (a-3.0)) + 1.0;
		else if (x < 2.0)
			z = -a*x*x*x + 5.0*a*x*x - 8.0*a*x + 4.0*a;
		return z;
	}
	
	

	
	/**Computes the covariance matrix
		copied from The Apache Commons Mathematics Library http://commons.apache.org/proper/commons-math/
		Applies threading for each column pair -> requires a lot of memory. Limited number of simultaneous threads to one i at a time
		@param a matrix for which the covariance matrix is to be computed
		@return covariance matrix
	*/
	public static double[][] getCovM(double[][] a){
		/**Run in parallel threads to make it snappy
			Fire a maximum of 40 threads at a time. Otherwise memory will run out...
			*/
		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<CovRunnable> runnables = new ArrayList<CovRunnable>();	
		double[][] covM = new double[a[0].length][a[0].length];
		
		
		for (int i = 0; i < a[0].length; ++i) {
			double[] temp1 = new double[a.length];
			for (int k = 0; k<a.length;++k){
				temp1[k] = a[k][i];
			}
			runnables.add(new CovRunnable(temp1,temp1,i,i));
			threads.add(new Thread(runnables.get(runnables.size()-1)));
			threads.get(threads.size()-1).start();
            for (int j = 0; j < i; ++j) {
				double[] temp2 = new double[a.length];
				for (int k = 0; k<a.length;++k){
					temp2[k] = a[k][j];
				}
				runnables.add(new CovRunnable(temp1,temp2,i,j));
				threads.add(new Thread(runnables.get(runnables.size()-1)));
				threads.get(threads.size()-1).start();
            }
			/**Catch the threads, and assign values to covariance matrix*/
			for (int j = 0; j<threads.size();++j){
				try{
					threads.get(j).join();
				}catch (Exception err){IJ.log("Cov failed");}
				
				int ti = runnables.get(j).i;
				int tj = runnables.get(j).j;
				double tvar = runnables.get(j).var;
				covM[ti][tj] = tvar;
				covM[tj][ti] = tvar;
			}
			/**Clear threads and runnables*/
			threads.clear();
			runnables.clear();
        }

		return covM;		
	}
	
	/**A helper class to enable multi threading*/
	public static class CovRunnable implements Runnable{
		public double[] a;
		public double[] b;
		public int i;
		public int j;
		public double var;
		public CovRunnable(double[] a, double[] b, int i,int j){
			this.a = a;
			this.b = b;
			this.i = i;
			this.j = j;
		}
		public void run(){
			var = Filters.covariance(a,b);
		}
	}
	
	/**
     * Computes the bias-corrected covariance between two arrays.
	 * copied from The Apache Commons Mathematics Library http://commons.apache.org/proper/commons-math/
	 * @param a data array
     * @param b data array
     * @return returns the covariance between the two arrays
     */
    public static double covariance(double[] a, double[] b){
        double result = 0d;
		double ma = Filters.mean(a);
		double mb = Filters.mean(b);
		for (int i = 0; i < a.length; i++) {
			double xDev = a[i] - ma;
			double yDev = b[i] - mb;
			result += ((a[i] - ma) * (b[i] - mb) - result) / (i + 1);
		}

        return result * ((double) a.length / (double)(a.length-1));
    }
	

}