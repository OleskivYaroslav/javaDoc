
public class ModelHess {

	
	public static void main(String[] args) {
		Methods m =new Methods(2);
		m.readMatrix();
		m.readB();
		m.MODELHESS();
		m.writeMatrix();
		m.writeMatrix1();
		m.CHOLSOLVE();
		m.writeX();

	}

}
