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
/*Own classes*/
import edu.deakin.timo.pixelFeatures.*;
import edu.deakin.timo.roiTools.*;
import edu.deakin.timo.detectEdges.*;

/** class IJEncricler.
*/
public class IJEncricler implements PlugIn{
	public void run(String arg) {
		IJ.log("Started Encircler Plugin");
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
		
		/**/

		/*Get image data*/
		int width = imp.getWidth();
		int height = imp.getHeight();
		short[] pixels = (short[]) imp.getProcessor().getPixels();
		double[][] data = new double[width][height];
		for (int r = 0;r<height;++r){
			for (int c = 0;c<width;++c){
				data[c][r] = (double) pixels[c+r*width];
			}
		}
		IJ.log("Got data");
		/**Apply ROI mask*/
		byte[][] roiMask = getRoiMask(imp);
		IJ.log("Got ROI");
		/*Threshold binary*/
		double[] minMax = new double[]{Filters.min(data),Filters.max(data)};
		IJ.log("Min "+minMax[0]+" max "+minMax[1]);
		double threshold = (minMax[1]-minMax[0])/10d+minMax[0];
		IJ.log("Threshold "+threshold);
		byte[][] binaryMask = new byte[width][height];
		for (int r = 0;r<height;++r){
			for (int c = 0;c<width;++c){
				if (data[c][r] >= threshold && roiMask[c][r] > 0){	/*If within ROI and above threshold, include in binary image*/
					binaryMask[c][r] = 1;
				}
			}
		}
		EdgeDetector ed = new EdgeDetector(binaryMask);
		Vector<DetectedEdge> edges = ed.edges;
		IJ.log("Got edges "+edges.size());
		Collections.sort(edges);
		DetectedEdge longestEdge = edges.get(0);
		IJ.log("Edge length "+longestEdge.iit.size());
		Polygon polygon = new Polygon();
		for (int j = 0;j<longestEdge.iit.size();j = j+5){
			polygon.addPoint(longestEdge.iit.get(j),longestEdge.jiit.get(j));
		}
		//polygon.addPoint(longestEdge.iit.get(0),longestEdge.jiit.get(0));	/*Add the initial point the 2nd time*/
		IJ.log("Polygon accumulated "+polygon.npoints);
		PolygonRoi roi = new PolygonRoi(polygon,Roi.POLYGON);
		/**Use colours to highlight ROIs*/
		roi.setStrokeColor(new Color(0f,1f,0f));

		RoiManager rMan;
		if (RoiManager.getInstance() == null){
			rMan = new RoiManager();
		}else{
			rMan = RoiManager.getInstance();
		}
		rMan.addRoi(roi);
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