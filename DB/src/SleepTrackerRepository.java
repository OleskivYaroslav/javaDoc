import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Locale;
import java.util.Scanner;

public class SleepTrackerRepository{
		public ArrayList<SleepRecord> List = new ArrayList<SleepRecord>();
		
		private void AddElement(String VATIN, Date StartTime, int duration){
			SleepRecord record = new SleepRecord();
			record.VATIN = VATIN;
			record.StartTime = StartTime;
			record.duration = duration;
			List.add(record);
		}
		private void addOrUpdate(String VATIN, Date StartTime, int duration){
			boolean b= true; 
			for(SleepRecord l : List){
				if ((l.VATIN.equals(VATIN)) & (l.StartTime == StartTime)){
					l.duration = duration;
					b=false;
				}
			}
			if(b){
				AddElement(VATIN,StartTime,duration);}
		}
		public void readingListFile(){
			String s2 = "";
			String[] s1 = new String[3];
			String s = "";
			Date d1;
			int i = 0;
			DateFormat format = new SimpleDateFormat("dd/mm/yyyy");
			Scanner scan = new Scanner("/home/syavo/DB/DB/Properties/File.txt");
			try{
			while (scan.hasNext()) 
			{
			    s2=scan.next();
				s1 = s2.split(",");
				s = s1[0];
				d1 = format.parse(s1[1]);
				i = Integer.parseInt(s1[2]);
				addOrUpdate(s, d1, i);
			}
			}catch(Exception exc){}
			scan.close();
		}
		public void readNewConsol(){
			String s = "";
			Date d1;
			DateFormat format = new SimpleDateFormat("E MMM dd HH:mm:ss z yyyy", Locale.US);
			int i = 0;
			Scanner scan1 = new Scanner(System.in);
			System.out.println("str");
			s = scan1.next();
			System.out.println("datetime");
			d1=null;
			try{
			d1 = format.parse(scan1.next());
			}catch(Exception exc){}
			System.out.println("int");
			i = Integer.parseInt(scan1.next());
			scan1.close();
			addOrUpdate(s,d1,i);
		}
		public void writeInFile(){
			try{
				FileWriter f=new FileWriter("File.txt");
				f.append("");
				f.close();
			}catch (IOException exc){}
			
			try{
				FileWriter f=new FileWriter("File.txt",true);
				for(SleepRecord l: List ){
					f.write(l.VATIN+","+l.StartTime+","+l.duration);
				}
				f.close();
			}catch (IOException exc){}
		}		
}
