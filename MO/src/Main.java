
public class Main {
	
	public double func1(double x){
		return Math.pow((x-5), 2);
	}
	
	public double func2(double x){
		return Math.pow((x-0.24),2)*100;
	}
	
	public double func(double x){
		return Math.pow(x,6)+3*Math.pow(x, 3)+(6*x)-1;
	}
	
	
	public static void main(String[] args) {
		Main my = new Main();
		Methods met=new Methods();
		double a=0,b=0;
		double x0=0;
		double h0 = 1;
		if(my.func(x0)<my.func(x0+h0)){h0=-h0;}
		while(my.func(x0)>my.func(x0+h0)){
			x0=x0+h0;
			h0=h0*2;
		}
		double xm=x0+h0;
		double xm_1=x0;
		double xm_2=x0-(h0/2);
		double xm_next=xm-h0/2;
		double[] xmas={my.func(xm_2),my.func(xm_1),my.func(xm),my.func(xm_next)};
		double[] xmas1={xm_2,xm_1,xm,xm_next};
		int j=0;
		double max=0;
		double min=9999999;
		for (int i =0;i<4;i++){
			if(xmas[i]>max){max = xmas[i]; j=i;}
		}
		xmas1[j]=xmas1[3];
		max=0;
		for (int i =0;i<4;i++){
			if(xmas1[i]>max){max = xmas1[i];}
			if(xmas1[i]<min){min = xmas1[i];}
		}
		a=min;
		b=max;
		System.out.println(xm+" "+xm_1+" "+xm_2+" "+xm_next+"  " +a+"  "+b);
		//met.GoldDividing(a, b, 0.00001); 
		met.ParabolMethod(a, b, 0.1);
	}

}
