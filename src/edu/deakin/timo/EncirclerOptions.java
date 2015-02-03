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
	private final String[] keys = {"Threshold","Invert","SkipPixels"};
	private final String[] defaults = {"300","false","5"};
	private String[] settings;
	private JTextField[] textFields;
	private JCheckBox checkBox;
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
		
		textFields = new JTextField[2];
		/*Add Textfields, and tickboxes*/
		//frame = new JFrame("EncriclerOptions");	//Add frame
		JPanel selections = new JPanel();
		selections.setLayout(new GridLayout(3,2,5,5));
		/*Threshold*/
		selections.add(new JLabel("Threshold"));
		textFields[0] = new JTextField(settings[0],6);
		selections.add(textFields[0]);
		/*Invert*/
		selections.add(new JLabel("Invert intensity"));
		checkBox = new JCheckBox();
		if (Boolean.parseBoolean(settings[1]) == true){
			checkBox.setSelected(true);
			//System.out.println("Check selected");
		}else{
			checkBox.setSelected(false);
		}
		selections.add(checkBox);
		/*Skip*/
		selections.add(new JLabel("Skip Pixels"));
		textFields[1] = new JTextField(settings[2],6);
		selections.add(textFields[1]);
		selections.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		selections.setPreferredSize(new Dimension(200,100));
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
			settings[1] = new Boolean(checkBox.isSelected()).toString();
			settings[2] = textFields[1].getText();
		}catch (Exception err){
			System.out.println(err);
		}
		return settings;
	}
}
