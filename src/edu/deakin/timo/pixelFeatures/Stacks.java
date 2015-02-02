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

public class Stacks{

	public static double[][][] get3DVarianceStack(double[][][] data){
		int width = data.length;
		int height = data[0].length;
		int depth = data[0][0].length;
		double[][][] varianceStack = new double[width][height][depth];
		ArrayList<Thread> threads = new ArrayList<Thread>();
		ArrayList<Var2DRunnable> runnables = new ArrayList<Var2DRunnable>();
		for (int d = 0; d < depth; ++d) {
			runnables.add(new Var2DRunnable(data,d));
			threads.add(new Thread(runnables.get(d)));
			threads.get(d).start();
		}
		
		/**Collect threads*/
		for (int d = 0; d < threads.size(); ++d){
			try{
				threads.get(d).join();	/**Wait for the thread to finish*/
			}catch (Exception err){IJ.log("TraceEdge threads failed");}
			double[][] var2D = runnables.get(d).var2D;
			for (int i = 0; i<width;++i){
				for (int j = 0; j<height;++j){
					varianceStack[i][j][d] = var2D[i][j];
				}
			}
        }
		return varianceStack;
	}
	
	public static class Var2DRunnable implements Runnable{
		double[][][] data;
		int d;
		public double[][] var2D;
		public Var2DRunnable(double[][][] data, int d){
			this.data = data;
			this.d =d;
		}
		public void run(){
			var2D = get2DvarianceImage(data,d);
		}
		
	}
	
	private static double[][] get2DvarianceImage(double[][][] data, int d){
		int width = data.length;
		int height = data[0].length;
		int depth = data[0][0].length;
		int x,y;
		double[][] varianceImage = new double[width][height];
		final int[][] neighbourhood = {{-1,-1},{-1,0},{-1,1},{0,-1},{0,0},{0,1},{1,-1},{1,0},{1,1}};
		ArrayList<Double> tempData = new ArrayList<Double>();
		for (int i = 0; i<width;++i){
			for (int j = 0; j<height;++j){
				tempData.clear();
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
				varianceImage[i][j] = Filters.getVariance(tempData.toArray(new Double[tempData.size()]));
			}
		}
		return varianceImage;
	}
	
	public static double[][][] getNormalized(double[][][] stack){
		int width = stack.length;
		int height = stack[0].length;
		int depth = stack[0][0].length;
		double[] norm = new double[]{Filters.min(stack),Filters.max(stack)};
		
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

	/**3D min stack*/
	public static double[][][] get3DMinStack(double[][][] stack){
		int width = stack.length;
		int height = stack[0].length;
		int depth = stack[0][0].length;
		double[][][] minStack = new double[width][height][depth];
		int x;
		int y;
		int[][] neighbourhood = {	{-1,-1},{-1,0},{-1,1},
							{0,-1},{0,0},{0,1},
							{1,-1},{1,0},{1,1}
							};
		for (int d = 0; d < depth; ++d) {
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
								tempData.add(stack[x][y][z]);
							}
						}
					}
					/*get the min from the neighbourhood*/
					minStack[i][j][d] = Filters.min(tempData.toArray(new Double[tempData.size()]));
				}
			}
		}
		return minStack;
	}
	
		/**3D min stack*/
	public static double[][][] getEdgeStack(double[][][] stack){
		int width = stack.length;
		int height = stack[0].length;
		int depth = stack[0][0].length;
		double[][][] edgeStack = new double[width][height][depth];
		double[][] slice = new double[width][height];
		double[][] edge;
		for (int d = 0; d < depth; ++d) {
			for (int i = 0; i<width;++i){
				for (int j = 0; j<height;++j){
					slice[i][j] = stack[i][j][d];
				}
			}
			edge = Filters.getEdgeImage(slice);
			for (int i = 0; i<width;++i){
				for (int j = 0; j<height;++j){
					edgeStack[i][j][d] = edge[i][j];
				}
			}
		}
		return edgeStack;
	}
	

}