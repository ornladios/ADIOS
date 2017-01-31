import ornl.adios.AdiosFile;
import ornl.adios.AdiosVar;

public class helloadios {

	public static void main(String[] args) {
		System.out.println("HelloAdios");

        System.out.println("user.dir = " + System.getProperty("user.dir"));
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));

		AdiosFile f = new AdiosFile("arrays.bp");
		f.printself();
		
		AdiosVar v = f.getVars().get("double_arr");
		v.printself();
		System.out.println("dim: " + v.getDims()[0] + ", " + v.getDims()[1]);
		double[] out = new double[20];
		int result = v.read(out);
		System.out.println("result: " + result);
		System.out.println("out: " + out[0] + " " + out[1] + " " + out[2] + " ...");
		System.out.println("out: " + out[10] + " " + out[11] + " " + out[12] + " ...");
				
		f.close();
		
		
	}
}
