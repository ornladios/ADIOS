package ornl.adios;

import java.util.HashMap;
import java.util.Map;

import ornl.adios.ext.*;

/**
 * Class for Adios writing.
 * 
 * @author jyc
 *
 */
public class AdiosWriter {
	static {
		try {
			System.loadLibrary("adiosjava");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}
	
	private long fd;
	private String fname;
	private String gname;
	private String mode;
	private int comm;
	
	protected Map<String, AdiosVarinfo> vars = new HashMap<String, AdiosVarinfo>();
	
	public AdiosWriter(String fname, String gname, String mode) {
		this.setFname(fname);
		this.setGname(gname);
		this.mode = mode;
		this.comm = 0;
		
		long[] fd = {0};
		adioslib.adios_open(fd, gname, fname, mode, comm);
		this.fd = fd[0];
	}

	public int write_int_arr(String vname, int[] value) {
		return adioslib.adios_write_int_arr(fd, vname, value);
	}

	public int write_double_arr(String vname, double[] value) {
		return adioslib.adios_write_double_arr(fd, vname, value);
	}
	
	public int write_int(String vname, int value) {
		return adioslib.adios_write_int(fd, vname, value);
	}

	public int write_double(String vname, double value) {
		return adioslib.adios_write_double(fd, vname, value);
	}

	public int close() {
		return adioslib.adios_close(fd);
	}

	public void printself() {
		System.out.println("=== AdiosWriter ===");
		System.out.println("fd: " + fd);
		System.out.println("fname: " + getFname());
		System.out.println("gname: " + getGname());
	}

	public String getFname() {
		return fname;
	}

	public void setFname(String fname) {
		this.fname = fname;
	}

	public String getGname() {
		return gname;
	}

	public void setGname(String gname) {
		this.gname = gname;
	}

	public Map<String, AdiosVarinfo> getVars() {
		return vars;
	}	
}
