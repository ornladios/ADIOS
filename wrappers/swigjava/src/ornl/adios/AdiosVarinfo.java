package ornl.adios;

import ornl.adios.ext.*;

/**
 * Class for variables to be written through Adios write APIs. 
 * 
 * @author jyc
 *
 */
public class AdiosVarinfo {
	static {
		try {
			System.loadLibrary("adiosjava");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}
	
	private long fd;
	private String name;
	private double[] value;
	
	public AdiosVarinfo(long fd, String name) {
		this.fd = fd;
		this.name = name;
	}

	public void printself() {
		System.out.println("=== AdiosVarinfo ===");
		System.out.println("fd: " + getFd());
		System.out.println("getName: " + getName());
	}
	
	public void write() {
		adioslib.adios_write_double_arr(fd, name, value);
	}

	public long getFd() {
		return fd;
	}
	
	public String getName() {
		return name;
	}

	public double[] getValue() {
		return value;
	}

	public void setValue(double[] value) {
		this.value = value;
	}

}
