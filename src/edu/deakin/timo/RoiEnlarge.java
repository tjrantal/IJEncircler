package edu.deakin.timo;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.measure.Calibration;	/*For obtaining pixel dimensions from original stack...*/
import ij.gui.NewImage;			/*For creating the output stack images*/
import ij.process.ImageProcessor;		/*For setting output stack image properties*/
import ij.io.FileInfo;			/*For setting image voxel dimensions...*/
import java.util.Properties;	/*For getting image properties*/
import java.util.*;				/*For enumeration*/
/*ROI stuff*/
import java.awt.Color;
import java.awt.Polygon;
import java.awt.Rectangle;
import ij.gui.Roi;
import ij.gui.PolygonRoi;
import ij.plugin.frame.RoiManager;
/*Settings*/
import java.awt.Frame;
import java.awt.Window;
/*Own classes*/
import edu.deakin.timo.pixelFeatures.*;
import edu.deakin.timo.roiTools.*;
import edu.deakin.timo.detectEdges.*;

/** class IJEncricler.
*/
public class RoiEnlarge implements PlugIn{
	/*Implement PlugIn*/
	public void run(String arg) {
		//IJ.log("Started Encircler Plugin");
		ImagePlus imp = WindowManager.getCurrentImage();
       	 /*Check that an image was open*/
		if (imp == null) {
	            IJ.noImage();
	            return;
       	 }
		/*Check that the image is 16 bit gray image*/
		if (imp.getType() != ImagePlus.GRAY16){
			IJ.error("IJGrower expects 16-bit greyscale data, e.g. DICOM-image");
			return;
		}
		
		/*Get options*/
		Frame[] frames = Frame.getFrames();
		int fr = 0;
		while (fr < frames.length && frames[fr].getTitle() != "EncirclerOptions"){
			++fr;
		}
		
		EncirclerOptions options;
		if (fr>= frames.length || frames[fr].getTitle() != "EncirclerOptions"){
			options = new EncirclerOptions();
		}else{
			if (frames[fr] instanceof edu.deakin.timo.EncirclerOptions){
				options = (EncirclerOptions) frames[fr];
			}else{
				options = new EncirclerOptions();	//Should never get here!
			}
		}
		String[] settings = options.getSettings();
		options.saveSettings();
		/*
		for (int i =0;i<settings.length;++i){
			System.out.println(options.keys[i]+" "+settings[i]);
		}
		*/
		int enlargeBy = Integer.parseInt(settings[6]);
		int skip = Integer.parseInt(settings[2]) > 0 ? Integer.parseInt(settings[2]) : 1; /*1 is min skip*/
		int width = imp.getWidth();
		int height = imp.getHeight();
		/**Apply ROI mask*/
		byte[][] roiMask = getRoiMask(imp);
		byte[][] binaryMask = new byte[width][height];
		int roiPixels = 0;
		for (int r = 0;r<height;++r){
			for (int c = 0;c<width;++c){
				if (roiMask[c][r] > 0){	/*If within ROI and above threshold, include in binary image*/
					binaryMask[c][r] = (byte) 1;
					++roiPixels;
				} else {
					binaryMask[c][r] = (byte) 0;
				} 
			}
		}
		
		/*Visualize binaryMask*/
		/*
		ImagePlus resultImage = NewImage.createByteImage("Binary Image",width,height,1, NewImage.FILL_BLACK);
		byte[] rPixels = (byte[])resultImage.getProcessor().getPixels();
		for (int r = 0;r<height;++r){
			for (int c = 0;c<width;++c){
				rPixels[c+r*width] = binaryMask[c][r];
			}
		}
		imp.getCalibration().copy();
		resultImage.setDisplayRange(0, 1);
        resultImage.show();
		*/

		EdgeDetector ed = new EdgeDetector(binaryMask,false);	/*Do not allow cleaving ever*/
		Vector<DetectedEdge> edges = ed.edges;
		//IJ.log("Got edges "+edges.size());
		Collections.sort(edges);
		DetectedEdge longestEdge = edges.get(edges.size()-1);
		//IJ.log("Edge length "+longestEdge.iit.size());
		Polygon polygon = new Polygon();
		for (int j = 0;j<longestEdge.iit.size();j = j+skip){
			polygon.addPoint(longestEdge.iit.get(j),longestEdge.jiit.get(j));
		}
		//System.out.println("Get centre "+polygon.npoints);
		//Enlarge ROI from centre of pixel area
		double[] centreCoordinates = new double[2];
		roiPixels = 0;
		for (int j = 0;j< height;++j){
			for (int i = 0; i < width;++i){
				if (binaryMask[i][j] == 1){
					centreCoordinates[0]+=(double) i;
					centreCoordinates[1]+=(double) j;
					++roiPixels;
				}
			}
		}
		centreCoordinates[0]/=(double)roiPixels;
		centreCoordinates[1]/=(double)roiPixels;
		//System.out.println("X "+centreCoordinates[0]+" Y "+centreCoordinates[1]);
		/*Calc polar coordinates, increment r by enlargeBy and calculate the enlarged ROI*/
		double x,y,r,t;
		for (int i =0;i<polygon.npoints;++i){
			x = ((double)polygon.xpoints[i])-centreCoordinates[0];
			y = ((double)polygon.ypoints[i])-centreCoordinates[1];
			t = Math.atan2(y,x);	//Theta
			r = Math.sqrt(x*x+y*y)+enlargeBy;	//r+5
			polygon.xpoints[i] = (int) (centreCoordinates[0]+r*Math.cos(t));
			polygon.ypoints[i] = (int) (centreCoordinates[1]+r*Math.sin(t));
		}
		
		//polygon.addPoint(longestEdge.iit.get(0),longestEdge.jiit.get(0));	/*Add the initial point the 2nd time*/
		//IJ.log("Polygon accumulated "+polygon.npoints);
		PolygonRoi roi = new PolygonRoi(polygon,Roi.POLYGON);
		/**Use colours to highlight ROIs*/
		roi.setStrokeColor(new Color(0f,0f,1f));
		imp.setRoi(roi);
		
    }

	/**
		Get a byte mask for the ROI of the current imp. 
		@param imp the imagePlus of the current ROI
		@return roiMask a byte array mask of the current image ROI; 1 = in ROI, 0 = out
	*/
	private byte[][] getRoiMask(ImagePlus imp){
		/*Determine mask for ROI-specific calculations, set mask to 1, if it belongs to the ROI, 0 otherwise*/
		int width = imp.getWidth();
		int height = imp.getHeight();
		int depth = imp.getStackSize();
		byte[][] roiMask = new byte[width][height];	/*Automatically initialized to zero*/
		Roi ijROI = imp.getRoi();	//Get the current ROI
		Rectangle rect = ijROI.getBoundingRect();
		
		/*Create ROI mask*/
		if (imp.getMask() != null){
			/*irregular roi, use Roi and bounding rectangle*/
			byte[] tempMask = (byte[]) imp.getMask().getPixels();	/*Out of mask = 0*/
			for (int j = rect.y;j< rect.y+rect.height;++j){
				for (int i = rect.x; i < rect.x+rect.width;++i){
					if (tempMask[i-rect.x+(j-rect.y)*rect.width] !=0){
						roiMask[i][j] =1;	/*In ROI = 1, out = 0*/
					}
				}
			}
		}else{
			/*rectangular ROI, use bounding rectangle*/
			for (int j = rect.y;j< rect.y+rect.height;++j){
				for (int i = rect.x; i < rect.x+rect.width;++i){
					roiMask[i][j] =1;	/*In ROI = 1, out = 0*/
				}
			}
		}
		return roiMask;
	}
}
