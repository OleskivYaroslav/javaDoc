import java.util.Scanner;


public class Methods {
	int n;
	double macheps,maxoffl, maxadd;
	double[] Sx; 
	double[][] H;
	double[][] L;
	double[] y;
	double[] x;
	double[] b;
	
	public Methods(int n){
		this.n=n;
		this.macheps=1;
		do{
			macheps/=2;
		}while(1+macheps!=1);
		this.macheps*=2;
		this.Sx=new double [n];
		this.b=new double [n];
		this.H=new double [n][n];
		this.L=new double [n][n];
		for(int i=0;i<n;i++){
			this.Sx[i]=1;
			for(int j=0;j<n;j++){
				this.L[i][j]=0;
				
			}
		}
		//this.maxoffl =0;
		//this.maxadd = 0;
		this.y=new double [n];
		this.x=new double[n];
	}
	Scanner scan = new Scanner(System.in);
	public void readMatrix(){
			for (int i=0; i<n;i++){
				for(int j=0; j<n;j++){
					H[i][j] = Double.parseDouble(scan.next());
				}
			}	
	}
	public void readB(){
		for (int i=0; i<n;i++){
			b[i] = Double.parseDouble(scan.next());
		}	
	}
	public void writeMatrix(){
		double[][] a=new double [n][n];
		a[0][0]=(L[0][0]*L[0][0])+(L[0][1]*L[0][1]);
		a[0][1]=(L[0][0]*L[1][0])+(L[0][1]*L[1][1]);
		a[1][0]=(L[1][0]*L[0][0])+(L[1][1]*L[0][1]);
		a[1][1]=(L[1][0]*L[1][0])+(L[1][1]*L[1][1]);
			for (int i=0; i<n;i++){
				for(int j=0; j<n;j++){
					System.out.print(a[i][j]+"  ");
				}
				System.out.println();
			}	
	}
	public void writeMatrix1(){
	
			for (int i=0; i<n;i++){
				for(int j=0; j<n;j++){
					System.out.print(L[i][j]+"  ");
				}
				System.out.println();
			}	
	}

	
	public void MODELHESS1(){
		//1
		for(int i=0;i<n;i++){
			for(int j=i; j<n;j++){
				H[i][j]=H[i][j]/(Sx[i]*Sx[j]);
			}
		}
		//2
		double sqrteps = Math.pow(macheps, 0.5);
		//3,4
		double max = H[0][0];
		double min = H[0][0];
		for(int i=0;i<n;i++){
			if(H[i][i]>max) max=H[i][i];
			if(H[i][i]<min) min=H[i][i];
		}
		double maxdiag = max;
		double mindiag = min;
		//5
		double maxposdiag = Math.max(0, maxdiag);
		//6
		double mu;
		if(mindiag<=(sqrteps*maxposdiag)){
			mu=2.0*(maxposdiag-mindiag)*sqrteps - mindiag;
			maxdiag=maxdiag+mu;		
		}else{
			mu=0;
		}
		//7
		max = Math.abs(H[0][1]);
		for(int i=0;i<n-1;i++){
			for(int j=i+1; j<n;j++){
				if(Math.abs(H[i][j])>max) max= Math.abs(H[i][j]);
			}
		}
		double maxoff = max;
		//8
		if((maxoff*(1+(2*sqrteps)))>maxdiag){
			mu=mu+(maxoff-maxdiag)+(2*sqrteps*maxoff);
			maxdiag = maxoff*(1+(2*sqrteps));
		}
		//9
		if(maxdiag==0){
			mu=1;
			maxdiag=1;
		}
		//10
		if(mu > 0){
			for(int i = 0; i<n;i++){
				H[i][i]=H[i][i]+mu;
			}
		}
		//11
		maxoffl=Math.pow(Math.max(maxdiag,(maxoff/n)),0.5);
		//12
		CHOLDECOMP ();
		//13
		if(maxadd>0){
			double maxev = H[0][0];
			double minev = H[0][0];
			double offrow;
			double sum1,sum2;
			for(int i = 0;i<n;i++){
				sum1=0;
				for(int j=0;j<i;j++){
					sum1+=Math.abs(H[j][i]);
				}
				sum2=0;
				for(int j=i+1;j<n;j++){
					sum2+=Math.abs(H[i][j]);
				}
				offrow = sum1+sum2;
				maxev = Math.max (maxev, H[i][i] + offrow);
				minev = Math.min (minev, (H[i][i] - offrow));
			}
			double sdd = ((maxev-minev)*sqrteps)-minev;
			sdd=Math.max(sdd, 0);
			mu =  Math.min(maxadd, sdd);
			for(int i = 0;i<n;i++){
				H[i][i]=H[i][i]+mu;
				for(int j=0;j<n;j++){
					L[i][j]=0;
				}
			}
			maxoffl=0;
			CHOLDECOMP();
		}
		//14
		for(int i=0;i<n;i++){
			for(int j=i;j<n;j++){
				H[i][j]=H[i][j]*Sx[i]*Sx[j];
			}
		}
		//15
		for(int i=0;i<n;i++){
			for(int j=0;j<i;j++){
				L[i][j]=L[i][j]*Sx[i];
			}
		}

	}
	public void CHOLDECOMP1(){
		//1
		double minl=Math.pow(macheps, 0.25)*maxoffl;
		//2
		double max=Math.abs(H[0][0]);
		if(maxoffl == 0){
			for(int i=0;i<n;i++){
				if(Math.abs(H[i][i])>max) max=Math.abs(H[i][i]);
			}
			maxoffl=Math.sqrt(max);
		}
		double minl2=Math.sqrt(macheps)*maxoffl;
		//3
		maxadd=0;
		//4
		double sum1;
		double minljj;
		for(int j=0;j<n;j++){
			sum1=0;
			for(int i =0;i<j-1;i++){
				sum1+=L[j][i]*L[j][i];
			}
			L[j][j]=H[j][j]-sum1;
			minljj=0;
			for(int i =j+1;i<n;i++){
				sum1=0;
				for(int k =0;k<j;k++){
				sum1+=(L[i][k]*L[j][k]);
				}
				L[i][j]=H[j][i]-sum1;
				minljj = Math.max (Math.abs(L[i][j]), minljj);
			}
			minljj = Math.max((minljj/maxoffl), minl);
			if(L[j][j]>(minljj*minljj)){
				L[j][j]=Math.pow(L[j][j], 0.5);
			}else{
				if(minljj<minl2){
					minljj=minl2;
				}
				maxadd=Math.max(maxadd,((minljj*minljj)-L[j][j]));
				L[j][j]=minljj;
			}
			for(int i =j+1;i<n;i++){
				L[i][j]=L[i][j]/L[j][j];
			}
		}	
	}
	
	public void CHOLSOLVE1(){
		//double s=0;
		LSOLVE ();
		LTSOLVE ();
		//s*=-1;
	}
	public void writeX(){
		for (int i=0; i<n;i++){
			System.out.println("x["+(i+1)+"] = "+x[i]);
		}
	}
	public void LSOLVE1(){	
		double sum;
		y[0]=Sx[0]/L[0][0];
		for(int i=1;i<n;i++){
			sum=0;
			for(int j=0;j<i;j++){
				sum+=L[i][j]*y[j];
			}
			y[i]=(Sx[i]-sum)/L[i][i];
		}
	}
	
	public void LTSOLVE1(){
		double sum;
		x[n-1]=y[n-1]/L[n-1][n-1];
		for(int i=n-2;i>=0;i--){
			sum=0;
			for(int j=i+1;j<n;j++){
				sum+=L[j][i]*x[j];
			}
			x[i]=(y[i]-sum)/L[i][i];
		}
	}


	public void MODELHESS()
    {
        //1
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                H[i][j] /= (Sx[i] * Sx[j]);
            }
        }
        //end 1
        double sqrtEps = Math.pow(macheps, 0.5);

        double maxDiag = H[0][ 0];
        for (int i = 1; i < n; i++)
        {
            if (maxDiag < H[i][i])
            {
                maxDiag = H[i][i];
            }
        }

        double minDiag = H[0][0];
        for (int i = 1; i < n; i++)
        {
            if (minDiag > H[i][i])
            {
                minDiag = H[i][i];
            }
        }
        //Console.WriteLine("sqrtEps=" + sqrtEps);
        double maxPosDiag = Math.max(0, maxDiag);
        double mu = 0;
        if (minDiag <= sqrtEps * maxPosDiag)
        {
            mu = 2 * (maxPosDiag - minDiag) * sqrtEps - minDiag;
            maxDiag += mu;
        }

        double maxOff = Math.abs(H[0][1]);
        //for (int i = 0; i < n; i++)
        //{
        //    for (int j = i + 1; j < n; j++) 
        //    {
        //        if (maxOff < Math.Abs(h[i, j]))
        //        {
        //            maxOff = Math.Abs(h[i, j]);
        //        }
        //    }
        //}
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if (maxOff < Math.abs(H[i][j]))
                {
                    maxOff = Math.abs(H[i][j]);
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (maxOff < Math.abs(H[i][j]))
                {
                    maxOff = Math.abs(H[i][j]);
                }
            }
        }

        if (maxOff * (1 + 2 * sqrtEps) > maxDiag)
        {
            mu += maxOff - maxDiag + 2 * sqrtEps * maxOff;
            maxDiag = maxOff * (1 + 2 * sqrtEps);
        }

        if (maxDiag == 0)
        {
            mu = 1;
            maxDiag = 1;
        }

        if (mu > 0)
        {
            for (int i = 0; i < n; i++)
            {
                H[i][i] += mu;
            }              
        }

        maxoffl = Math.pow((Math.max(maxDiag, maxOff / n)),0.5);

        //double maxAdd;
        CHOLDECOMP();
        //Console.WriteLine("\nmaxAdd1= {0}", maxAdd);
        if (maxadd > 0)
        {
            double maxEv = H[0][0];
            double minEv = H[0][0];
            double offRow = 0;
            double s1 = 0;
            double s2 = 0;
            for (int i = 0; i < n; i++)
            {
                s1 = 0;
                s2 = 0;
                for (int j = 0; j < i; j++)
                {
                    s1 += Math.abs(H[j][i]);
                }
                for (int j = i + 1; j < n; j++)
                {
                    s2 += Math.abs(H[i][j]);
                }
                offRow = s1 + s2;
                maxEv = Math.max(maxEv, H[i][i] + offRow);
                minEv = Math.min(minEv, H[i][i] - offRow);
            }
            double sdd = (maxEv - minEv) * sqrtEps - minEv;
            sdd = Math.max(sdd, 0);
            mu = Math.min(maxadd, sdd);
            for (int i = 0; i < n; i++)
            {
                H[i][i] += mu;
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    L[i][j] = 0;
                }
            }
            maxadd = 0;
            CHOLDECOMP();
            //14
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    H[i][j] *= Sx[i] * Sx[j];
                }
            }
            //end 14

            //15
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    L[i][j] *= Sx[i];
                }
            }
            //end 15
        }
    }

	
	public void CHOLDECOMP ()
    {
        double minL = Math.pow(macheps, 0.25) * maxoffl;
        double s = Math.abs(H[0][0]);
        if (maxoffl == 0)
        {
            for (int i = 1; i < n; i++)
            {
                if (s < Math.abs(H[i][i]))
                {
                    s = Math.abs(H[i][i]);
                }
            }
            maxoffl = Math.sqrt(s);
        }
        double minL2 = Math.sqrt(macheps) * maxoffl;

        maxadd = 0;

        double minLjj;
        for (int j = 0; j < n; j++)
        {
            s = 0;
            for (int i = 0; i < j; i++)
            {
                s += Math.pow(L[j][i], 2);
            }

            L[j][j] = H[j][j] - s;

            minLjj = 0;
            for (int i = j + 1; i < n; i++)
            {
                s = 0;
                for (int k = 0; k < j; k++)
                {
                    s += L[i][k] * L[j][k];
                }
                L[i][j] = H[j][i] - s;
                minLjj = Math.max(Math.abs(L[i][j]), minLjj);
            }

            minLjj = Math.max(minLjj / maxoffl, minL);
            if (L[j][j] > Math.pow(minLjj, 2))
            {
                L[j][j] = Math.sqrt(L[j][j]);
            }
            else
            {
                if (minLjj < minL2)
                {
                    minLjj = minL2;
                }
                maxadd = Math.max(maxadd, Math.pow(minLjj, 2) - L[j][j]);
                L[j][j] = minLjj;
            }
            for (int i = j + 1; i < n; i++)
            {
                L[i][j] /= L[j][j];
            }
        }
    }


	public void CHOLSOLVE()
    {
        
		LSOLVE();
        
        LTSOLVE();
        
    }

    public void LSOLVE()
    {
        double s;
        y[0] = b[0] / L[0][0];
        for (int i = 1; i < n; i++)
        {
            s = 0;
            for (int j = 0; j < i; j++)
            {
                s += L[i][j] * y[j];
            }
            y[i] = (b[i] - s) / L[i][i];
        }
    }

    public void LTSOLVE()
    {
        double s;
        x[n - 1] = y[n - 1] / L[n - 1][ n - 1];
        for (int i = n - 2; i >= 0; i--)
        {
            s = 0;
            for (int j = i + 1; j < n; j++)
            {
                s += L[j][i] * x[j];
            }
            x[i] = (y[i] - s) / L[i][i];
        }
    }

}
