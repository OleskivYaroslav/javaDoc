import java.util.Scanner;
public class QRDECOMP {
	boolean sing = false;
	int n;
	double [][] M;
	double [] M1;
	double [] M2;
	double [] b;
	
	public QRDECOMP(int n){
		this.n = n;
		this.M = new double [n][n];
		this.M1 = new double [n];
		this.M2 = new double [n];
		this.b = new double [n];
	} 
	Scanner scan = new Scanner(System.in);
	public double [][]readMatrixM(){
			for (int i=0; i<n;i++){
				for(int j=0; j<n;j++){
					M[i][j] = Double.parseDouble(scan.next());
				}
			}	
		
		return M;
	}
	public double []readB(){
		for (int i=0; i<n;i++){
			b[i] = Double.parseDouble(scan.next());
		}	
		return b;
	}
	public void writeMatrix(){
		
			for (int i=0; i<n;i++){
				for(int j=0; j<n;j++){
					System.out.print(M[i][j]+"  ");
				}
				System.out.println();
			}	
	}
	public void QRDECOMP1(){
		double eta=0;
		double sigma=0;
		double tau=0;
		for(int k=0; k<n-1;k++){
			for(int j=k;j<n;j++){
				if (eta<Math.abs(M[j][k])) eta=Math.abs(M[j][k]);
			}
			if(eta==0){
				M1[k]=0;
				M2[k]=0;
				sing = true;
			}else{
				for(int i=k;i<n;i++){
					M[i][k]=M[i][k]/eta;
				}
				double sum=0;
				for(int i=k;i<n;i++){
					sum+=M[i][k]*M[i][k];
				}
				sum=Math.sqrt(sum);
				sigma = sum;
				if(M[k][k]<0){
					sigma*=-1;
				}
				M[k][k]+=sigma;
				M1[k]=M[k][k]*sigma;
				M2[k]=-eta*sigma;
				
				for(int j=k+1;j<n;j++){
					double sum1=0;
					//tau=0;
					for(int i=k;i<n;i++){
						sum1+=M[i][k]*M[i][j];
					}
					tau=sum1/M1[k];
					for(int i=k;i<n;i++){
						M[i][j]=M[i][j]-(tau*M[i][k]);
					}
				}
			}
		}
		if(M[n-1][n-1]==0) sing=true;
		M2[n-1]=M[n-1][n-1];
		
	}
	public void QRSOLVE(){
		double sum=0;
		double tau;
		for (int j=0;j<n-1;j++){
			sum = 0;
			for(int k=j; k<n;k++){
				sum+=M[k][j]*b[k];
			}
			tau=sum/M1[j];
			for(int i=j; i<n;i++){
				b[i]-= tau*M[i][j];
			}
		}
		RSOLVE();
	}
	public void RSOLVE(){
		b[n-1]/=M2[n-1];
		double sum;
		for(int i=n-2; i>=0;i--){
			sum=0;
			for(int j=i+1;j<n;j++){
				sum+=M[i][j]*b[j];
			}
			b[i]=(b[i]-sum)/M2[i];
		}
	}
	public void writeX(){
		for (int i=0; i<n;i++){
			System.out.println("x["+(i+1)+"] = "+b[i]);
		}
	}
}
