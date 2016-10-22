

public class DB {
	
	public static void main(String[] args){
		SleepTrackerRepository rt = new SleepTrackerRepository();
		rt.readingListFile();
		rt.readNewConsol();
		rt.writeInFile();
	}

}
