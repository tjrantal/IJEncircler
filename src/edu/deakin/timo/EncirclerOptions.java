package edu.deakin.timo;

/*Saving settings*/
import java.util.prefs.Preferences;		/*Saving the file save path -> no need to rebrowse...*/
/*The window frame*/
import javax.swing.*;	/*JFrame*/
import java.awt.*;		/*Layout*/

/** class EncirclerOptions.
*/
public class EncirclerOptions extends JFrame {
	private Preferences preferences;		/**Saving the default file path*/
	public final String[] keys = {"Threshold","Invert","SkipPixels","EnableCleaving","CleaveMinRatio","CleaveMinLength","EnlargeBy"};
	private final String[] defaults = {"300","false","5","true","10","4","5"};
	private String[] settings;
	private JTextField[] textFields;
	private JCheckBox[] checkBoxes;
	//private JFrame frame;
	public EncirclerOptions(){
		super("EncirclerOptions");
		/*Get Preferences*/
		settings = new String[keys.length];
		preferences = Preferences.userRoot().node(this.getClass().getName());
		//System.out.println("Get prefs from storage");
		for (int i = 0; i<keys.length;++i){
			settings[i] = preferences.get(keys[i],defaults[i]); /*Use current working directory as default*/
			//System.out.println("Storage "+keys[i]+" "+settings[i]);
		}
		
		checkBoxes = new JCheckBox[2];
		textFields = new JTextField[keys.length-checkBoxes.length];
		/*Add Textfields, and tickboxes*/
		//frame = new JFrame("EncriclerOptions");	//Add frame
		JPanel selections = new JPanel();
		selections.setLayout(new GridLayout(keys.length,2,5,5));
		/*Threshold*/
		selections.add(new JLabel("Threshold"));
		textFields[0] = new JTextField(settings[0],6);
		selections.add(textFields[0]);
		/*Invert*/
		selections.add(new JLabel("Invert intensity"));
		checkBoxes[0] = new JCheckBox();
		if (Boolean.parseBoolean(settings[1]) == true){
			checkBoxes[0].setSelected(true);
			//System.out.println("Check selected");
		}else{
			checkBoxes[0].setSelected(false);
		}
		selections.add(checkBoxes[0]);
		/*Skip*/
		selections.add(new JLabel("Skip Pixels"));
		textFields[1] = new JTextField(settings[2],6);
		selections.add(textFields[1]);
		/*Cleaving*/
		selections.add(new JLabel("Enable cleaving"));
		checkBoxes[1] = new JCheckBox();
		if (Boolean.parseBoolean(settings[3]) == true){
			checkBoxes[1].setSelected(true);
			//System.out.println("Check selected");
		}else{
			checkBoxes[1].setSelected(false);
		}
		selections.add(checkBoxes[1]);
		/*Cleaving ratio*/
		selections.add(new JLabel("Cleave min ratio"));
		textFields[2] = new JTextField(settings[4],6);
		selections.add(textFields[2]);
		/*Cleaving length*/
		selections.add(new JLabel("Cleave min length"));
		textFields[3] = new JTextField(settings[5],6);
		selections.add(textFields[3]);
		/*Enlarge by*/
		selections.add(new JLabel("Enlarge by"));
		textFields[4] = new JTextField(settings[6],6);
		selections.add(textFields[4]);
		/*Make frame visible*/
		selections.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		selections.setPreferredSize(new Dimension(200,40*keys.length));
		this.getContentPane().add(selections);
		this.setLocation(20,20);
		this.pack();
		this.setVisible(true);		
	}
	
	public void saveSettings(){
		/*Save Preferences*/
		//System.out.println("Saving prefs");
		for (int i = 0; i<keys.length;++i){
			preferences.put(keys[i],settings[i]);
			//System.out.println(keys[i]+" "+settings[i]);
		}

	}
	
	public String[] getSettings(){
		try{
			settings[0] = textFields[0].getText();
			settings[1] = new Boolean(checkBoxes[0].isSelected()).toString();
			settings[2] = textFields[1].getText();
			settings[3] = new Boolean(checkBoxes[1].isSelected()).toString();
			settings[4] = textFields[2].getText();
			settings[5] = textFields[3].getText();
			settings[6] = textFields[4].getText();
		}catch (Exception err){
			System.out.println(err);
		}
		return settings;
	}
}
