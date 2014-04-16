package primes;

public class PrimesGUI {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Prime p= new Prime();
		new Thread(p).start();
		GUI g= new GUI();
		new Thread(g).start();

	}

}
