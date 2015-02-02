Encircler ImageJ plugin

takes the ROI in the current image, and traces around the longest border within the ROI, sets the ROI to the border.

The ROI is binarized (pixels assigned a value of 0, and 1), and blob borders are subsequently traced. Borders are identified by starting from image coordinates 0,0, and advancing row by row until finding a 1. Thereafter the border is traced, and all pixels within the border set to 1. The search for the next border then continues until the maxWidth,maxHeight pixel is reached.

Settings panel;
Settings are stored with Java persistent storage java.util.prefs.Preferences;	
threshold				= used to  the image. Anything < will be assigned 0, >= 1
invert					= flips the binarised values (0 to 1, and 1 to 0)
pixels between points	= coarseness of the ROI polygon, the number of border pixels skipped between polygon points