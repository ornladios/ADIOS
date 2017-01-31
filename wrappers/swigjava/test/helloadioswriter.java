import ornl.adios.Adios;
import ornl.adios.AdiosWriter;

public class helloadioswriter {

	public static void main(String[] args) {
		System.out.println("HelloAdiosWriter");

        System.out.println("user.dir = " + System.getProperty("user.dir"));
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));

        Adios.init("arrays.xml");
        
        AdiosWriter fw = new AdiosWriter("arrays.bp", "arrays", "w");
		fw.printself();
		
		double[] value = new double[20];
		for (int i = 0; i < value.length; i++) {
			value[i] = i;
		}
		
		fw.write_int("NX", 10);
		fw.write_int("NY", 2);
		fw.write_double_arr("double_arr", value);
		fw.close();
	}

}
