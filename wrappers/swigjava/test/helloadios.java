import ornl.adios.AdiosFile;
import ornl.adios.AdiosVar;

public class helloadios {

	public static void main(String[] args) {
		System.out.println("HelloAdios");

        System.out.println("user.dir = " + System.getProperty("user.dir"));
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));

		AdiosFile f = new AdiosFile("adios_global.bp");
		f.printself();
		
		AdiosVar v = f.getVars().get("temperature");
		v.printself();
		System.out.println("dim: " + v.getDims()[0] + ", " + v.getDims()[1]);
		double[] out = new double[20];
		int result = v.getValue(out);
		System.out.println("result: " + result);
		System.out.println("out: " + out[0] + " " + out[1] + " " + out[2] + " ...");
		System.out.println("out: " + out[10] + " " + out[11] + " " + out[12] + " ...");
		
		double[][] out2 = new double[2][10];
		//int result = v.getValue(out2);
		//System.out.println("result: " + result);
		//System.out.println("out: " + out2[0][0] + " " + out2[0][1] + " " + out2[0][2] + " ...");
		//System.out.println("out: " + out2[1][0] + " " + out2[1][1] + " " + out2[1][2] + " ...");

		v.read(out2);
		f.close();
	}
}
