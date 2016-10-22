
public class QR {
	public static void main(String[] args) {
		QRDECOMP qr = new QRDECOMP(2);
		qr.readMatrixM();
		qr.readB();
		qr.QRDECOMP1();
		qr.QRSOLVE();
		qr.writeX();
	}

}
