import ornl.adios.AdiosFile;

public class helloadios {

	public static void main(String[] args) {
		System.out.println("HelloAdios");

        System.out.println("user.dir = " + System.getProperty("user.dir"));
        System.out.println("java.library.path = " + System.getProperty("java.library.path"));

		AdiosFile f = new AdiosFile("test.bp");
		f.printself();
		f.close();
	}
}
