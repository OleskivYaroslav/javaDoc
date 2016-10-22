
public class Methods {
	
	Main my = new Main();
	
	public void GoldDividing(double a, double b, double eps){
		double a0=a;
		double b0=b;
		double x1=a0+(((3-Math.sqrt(5))/2)*(b0-a));
		double x2=a0+b0-x1;
		int k=1;
		if(my.func(x1)<=my.func(x2)){
			b0=x2;
		}else{
			a0=x1;	
		}
		do{
		if(my.func(x1)<=my.func(x2)){
			b0=x2;
			System.out.println(k+"  "+x1+"  "+
					x2+"  "+my.func(x1)+"  "+
					my.func(x2)+"  "+a0+" "+b0+"  "+Math.abs(a0-b0));
			x2=x1;
			x1=a0+b0-x2;
		}else{
			a0=x1;
			System.out.println(k+"  "+x1+"  "+
					x2+"  "+my.func(x1)+"  "+
					my.func(x2)+"  "+a0+" "+b0+"  "+Math.abs(a0-b0));
			x1=x2;
			x2=a0+b0-x1;
		}
		k++;
		}while(Math.abs(a0-b0)>eps);
		System.out.println(k+"  "+x1+"  "+
						x2+"  "+my.func(x1)+"  "+
						my.func(x2)+"  "+a0+" "+b0+"  "+Math.abs(a0-b0));
		System.out.println("ANSWER = " +((a0+b0)/2));

	}
	
	public void ParabolMethod(double a, double b, double eps){
		
		double x1=a;
		double x2=(a+b)/2;
		double x3=b;
		double delta_m=0;
		double delta_p=0;
		double zm1=0;
		double zm2=0;
		double zm3=0;
		double zm4=0;
		double x_=0;
		double t,k=0;
		
		do{
			 delta_m=(Double)(my.func(x1)-my.func(x2));
			 delta_p=(Double)(my.func(x3)-my.func(x2));
			 zm1=(Double)Math.pow((x3-x2),2)*delta_m;
			 zm2=(Double)Math.pow((x2-x1),2)*delta_p;
			 zm3=(Double)(x3-x2)*delta_m;
			 zm4=(Double)(x2-x1)*delta_p;
			 x_=x2+(0.5*((zm1-zm2)/(zm3+zm4)));
			 if(x_==x2){t=((x2+x1)/2.0);}else{t=x_;}
			 //double[] xmas={my.func(x1),my.func(t),my.func(x2),my.func(x3)};
			 
			 if ((my.func(t) < my.func(x2)) && (my.func(t) < my.func(x1)) && (my.func(t) < my.func(x3)))
             {
                 x3 = x2;
                 x2 = t;
             }
             if ((my.func(x2) < my.func(t)) && (my.func(x2) < my.func(x1)) && (my.func(x2) < my.func(x3)))
             {
                 x1 = t;
             }
			if(k>1000) break;
		k++;
		}while(Math.abs(x1-x3)>eps);
		double vidp=((x1+x3)/2.0);
		System.out.println("ANSWER = " + vidp);
	}

}
