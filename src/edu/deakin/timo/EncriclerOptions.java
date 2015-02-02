package edu.deakin.timo;

/*Saving settings*/
import java.util.prefs.Preferences;		/*Saving the file save path -> no need to rebrowse...*/


/** class EncriclerOptions.
*/
public class EncriclerOptions {
	private Preferences preferences;		/**Saving the default file path*/
	private final String[] keys = {"Threshold","Invert","SkipPixels"};
	private final String[] defaults = {"300","false","5"};
	private String[] settings;
	public EncirclerOptions(){
		/*Get Preferences*/
		settings = new String[keys.length];
		preferences = Preferences.userRoot().node(this.getClass().getName());
		try{
			for (int i = 0; i<keys.length;++i){
				settings[i] = preferences.get(keys[i],defaults[i]); /*Use current working directory as default*/
			}
		}catch (IOException ex){
			System.out.println(ex);
			for (int i = 0; i<keys.length;++i){
				settings[i] = defaults[i]; /*Use current working directory as default*/
			}
		}
	}
	
	public void savePreferences(){
			/*Save Preferences*/
		try{
			for (int i = 0; i<keys.length;++i){
				preferences.put(keys[i],settings[i]);
			}
		}catch (IOException ex){
			System.out.println(ex);
			savePath = ".";
		}
	}
	
	public void getPreferences(){
	
	}
}