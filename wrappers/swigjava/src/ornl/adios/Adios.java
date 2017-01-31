package ornl.adios;

import ornl.adios.ext.*;

public class Adios {

	static {
		try {
			System.loadLibrary("adiosjava");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}
	
	static public int init(String config) {
		return adioslib.adios_init(config, 0);
	}

	static public void fin() {
		adioslib.adios_finalize(0);
	}
}
