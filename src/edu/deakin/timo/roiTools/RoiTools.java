package edu.deakin.timo.roiTools;

import ij.*;
import ij.gui.*;
import java.awt.*;
import java.util.*;	
import edu.deakin.timo.detectEdges.*;
import edu.deakin.timo.pixelFeatures.*;

public class RoiTools{
	
	public static double[][][] getFeatureStack(ImagePlus imp){
		int width = imp.getWidth();
		int height = imp.getHeight();
		short[] tempPointer = (short[]) imp.getProcessor().getPixels();
		double[][] tempData = new double[width][height];
		/*Insert data into tempData*/
		for (int i = 0; i< width; ++i){
			for (int j = 0; j< height; ++j){
				tempData[i][j] = tempPointer[i+j*width];
			}
		}
		/*Get edges with convolutions*/
		double[][] edgeImage = getEdgeImage(tempData);
		/**Calculate image features LBP, variance, gradient*/
		//double[][] varianceImage = Filters.getVarianceImage(tempData, 1);
		
		/**Create the feature stack, normalize form 0 to 1*/
		double[][][] featureStack = new double[2][width][height];
		double[] maxs = new double[]{Filters.max(tempData),Filters.max(edgeImage)/*,Filters.max(varianceImage)*/};
		double[] mins = new double[]{Filters.min(tempData),Filters.min(edgeImage)/*,Filters.min(varianceImage)*/};
		
		for (int i = 0;i<featureStack[0].length;++i){
			for (int j = 0;j<featureStack[0][i].length;++j){
				featureStack[0][i][j] = (tempData[i][j]-mins[0])/(maxs[0]-mins[0]);
				featureStack[1][i][j] = (edgeImage[i][j]-mins[1])/(maxs[1]-mins[1]);
				//featureStack[2][i][j] = (varianceImage[i][j]-mins[2])/(maxs[2]-mins[2]);
			}
		}
		return featureStack;
	}

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
		double[][][] convKerns = {
									{{-1,-2,-1},{0,0,0},{1,2,1}},
									{{-1,0,1},{-2,0,2},{-1,0,1}},
									{{-2,-1,0},{-1,0,1},{0,1,2}},
									{{0,-1,-2},{1,0,-1},{2,1,0}}
								};
		double[][][] tempConvs = new double[convKerns.length][data.length][data[0].length];
		for (int i = 0; i< convKerns.length;++i){
			tempConvs[i] =  Filters.conv2(data, convKerns[i]);
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


	/**
		Determine the corners of the roi
		@param roi the PolygonRoi roi to define corners for
		@param imp the image the roi is defined for
		@return corners 4 by 2 array with x-, and y- coordinates of the corners
	
	*/
	public static Roi getCorners(Roi roi,ImagePlus imp){
		byte[][] mask = getMaskFromRoi(imp, roi);
		int[] centreCoords = getMaskCentre(mask);
		int[][] corners = null;
		/*Calculate the roi polygon coordinates in polar coordinates*/
		//IJ.log("InstanceOf");
		if (roi instanceof PolygonRoi){
			//IJ.log("Yes instanceOf");
			Polygon polygon = ((PolygonRoi)roi).getPolygon();
			Double[][] polar = new Double[polygon.npoints][4];	/*0 r, 1 theta, 2 x, 3 y*/
			for (int i=0; i<polygon.npoints;++i){
				double r,theta;
				double a,b;
				a = (double) (polygon.xpoints[i]-centreCoords[0]);
				b = (double) (polygon.ypoints[i]-centreCoords[1]);
				r = Math.sqrt(Math.pow(a,2d)+Math.pow(b,2d));
				theta = Math.atan2(b,a);
				polar[i][0] = r;
				polar[i][1] = theta;
				polar[i][2] = new Double((double) polygon.xpoints[i]);
				polar[i][3] = new Double((double) polygon.ypoints[i]);
			}
			/*sort the coordinates, theta 0 = towards right, theta increases clock-wise in ImageJ image (where y increases towards bottom!*/
			/**Prepare a comparator*/
			Comparator<Double[]> comparator = new Comparator<Double[]>() {
				@Override
				public int compare(Double[] o1, Double[] o2) {
					return o1[1].compareTo(o2[1]);
				}
			};
			Arrays.sort(polar,comparator);	/**Sort the polar coordinates into ascending order according to theta*/
			/**Look for the four corners, take the maximum r in the four sectors -pi to -pi/2; -pi/2 to 0; 0 to pi/2; pi/2 to pi*/
			int ind = 0;
			double[][] polarCorners = new double[4][3];
			for (int i = 0;i<4;++i){
				double maxR = -1;
				double maxTheta = 0;
				double maxInd = -1;
				while (ind < polar.length && polar[ind][1] < ((2*Math.PI*((double)i+1d)/4d)-Math.PI)){
					if (polar[ind][0] > maxR){
						maxR = polar[ind][0];
						maxTheta = polar[ind][1];
						maxInd = ind;
					}
					++ind;
				}
				polarCorners[i] = new double[]{maxR,maxTheta,maxInd};			
			}
			/**Get the points along the edge in between*/
			double stepNumber = 5d;
			double[][] edgePoints = new double[4*((int) stepNumber)][2];
			ind = -1;
			for (int i = 0; i<(polarCorners.length-1); ++i){
				
				++ind;
				/**Go along the border from previous to next corner, sample three points along the way*/
				int initI = (int) polarCorners[i][2];
				/*insert the corner point*/
				edgePoints[ind][0] = polar[initI][2];
				edgePoints[ind][1] = polar[initI][3];
				
				int endI = (int) polarCorners[i+1][2];
				
				double step = ((double)(endI-initI))/stepNumber;
				int multiple = 1;
				IJ.log("Corner "+i+" initI "+initI+" endI "+endI+" r "+edgePoints[ind][0]+" theta "+edgePoints[ind][1]);
				for (int ii = initI+(int)step;ii<(int) (endI-step/2);ii=initI+((int) (step*(double)multiple))){
					++ind;
					++multiple;
					edgePoints[ind][0] = polar[ii][2];
					edgePoints[ind][1] = polar[ii][3];
					//IJ.log("Int Corner "+i+" ind "+ii+" r "+edgePoints[ind][0]+" theta "+edgePoints[ind][1]);
				}
				
			}
			/*The final edge needs to be handled manually due to the discontinuity*/
			/*Create an array with the indices*/
			int initI = (int) polarCorners[polarCorners.length-1][2];
			int tempEndInd = polar.length-1;
			int endI = (int) polarCorners[0][2];
			int[] indiceArray = new int[tempEndInd-initI+1+endI+1];
			int tI = 0;
			for (int i = initI;i<polar.length;++i){
				indiceArray[tI] = i;
				++tI;
			}
			for (int i = 0;i<=endI;++i){
				indiceArray[tI] = i;
				++tI;
			}
			
			/*add the missing points for the last border*/
			initI = 0;
			/*insert the corner point*/
			++ind;
			edgePoints[ind][0] = polar[(int) polarCorners[polarCorners.length-1][2]][2];
			edgePoints[ind][1] = polar[(int) polarCorners[polarCorners.length-1][2]][3];
			endI = indiceArray.length-1;
			double step = ((double)(endI-initI))/stepNumber;
			//IJ.log("Corner "+i+" ind "+initI+" r "+edgePoints[ind][0]+" theta "+edgePoints[ind][1]);
			int multiple = 1;
			for (int ii = initI+(int)step;ii<(int) (endI-step/2);ii=initI+((int) (step*(double)multiple))){
				++ind;
				++multiple;
				edgePoints[ind][0] = polar[indiceArray[ii]][2];
				edgePoints[ind][1] = polar[indiceArray[ii]][3];
				//IJ.log("Int Corner "+i+" ind "+ii+" r "+edgePoints[ind][0]+" theta "+edgePoints[ind][1]);
			}
			
			//IJ.log("I length "+indiceArray.length+" I;s added "+tI);
			
			
			/*Create the roi polygon*/
			int[] xc = new int[edgePoints.length];
			int[] yc = new int[edgePoints.length];
			for (int i = 0; i<edgePoints.length;++i){
				xc[i] = (int) edgePoints[i][0];
				yc[i] = (int) edgePoints[i][1];
				IJ.log("XC "+xc[i]+" YC "+yc[i]);
			}
			//xc[edgePoints.length] = (int) edgePoints[0][0];
			//yc[edgePoints.length] = (int) edgePoints[0][1];
			return new PolygonRoi(new Polygon(xc,yc,xc.length),Roi.POLYGON);
			
		}
		return new PolygonRoi(new Polygon(),Roi.POLYGON);
	}
	
	public static byte[][] getMaskFromRoi(ImagePlus imp, Roi ijROI){
		int width = imp.getWidth();
		int height = imp.getHeight();
		byte[][] segmentationMask = new byte[width][height];
		short[] tempPointer = (short[]) imp.getProcessor().getPixels();
		if (ijROI != null){	/*Set the seedRoi according to manual ROI*/
			/*Check whether pixel is within ROI*/
			for (int j = 0;j< height;j++){
				for (int i = 0; i < width;i++){
					if (ijROI.contains(i,j)){
						segmentationMask[i][j] = 1;
					}else{
						segmentationMask[i][j] = 0;
					}
				}
			}
			/*Check whether a polygon can be acquired and include polygon points too*/
			Polygon polygon = ijROI.getPolygon();
			if (polygon != null){
				for (int j = 0;j< polygon.npoints;j++){
					segmentationMask[polygon.xpoints[j]][polygon.ypoints[j]] = 1;
				}
			}
		}
		return segmentationMask;
	}
	
	public static long getMaskArea(byte[][] segmentationMask){
		long segmentationArea = 0;
		/*Check whether pixel is within ROI*/
		for (int i = 0;i< segmentationMask.length;i++){
			for (int j = 0; j < segmentationMask[i].length;++j){
				if (segmentationMask[i][j]>0){
					++segmentationArea;
				}
			}
		}
		return segmentationArea;
	}
	
	
	public static int[] getMaskCentre(byte[][] mask){
		double[] coords = {0d,0d,0d};
		for (int i = 0; i<mask.length;++i){
			for (int j = 0; j<mask[i].length;++j){
				if (mask[i][j] > 0){
					coords[0]+=i;
					coords[1]+=j;
					++coords[2];
				}
			}
		}
		coords[0]/=coords[2];
		coords[1]/=coords[2];
		int[] coordRet = {(int) coords[0],(int) coords[1]};
		return coordRet;
	}
		
}