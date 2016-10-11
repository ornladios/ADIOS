package ornl.adios;

import java.math.BigInteger;

import java.util.*;

import ornl.adios.adioslib.*;

public class AdiosFile {
	static {
		try {
			System.loadLibrary("adiosjava");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}

	protected ADIOS_FILE f;
	private int comm;
	private int comm_rank;

	protected Map<String, AdiosVar> vars = new HashMap<String, AdiosVar>();

	public AdiosFile(String fname) {
		comm = 1; // MPI_COMM_WORLD
		comm_rank = 0;

		ADIOS_READ_METHOD method = ADIOS_READ_METHOD.ADIOS_READ_METHOD_BP;
		f = adioslib.adios_read_open_file(fname, method, comm);
		for (String vname : getVar_namelist()) {
			AdiosVar v = new AdiosVar(f, vname);
			vars.put(vname, v);
		}
	}

	public int close() {
		return adioslib.adios_read_close(f);
	}

	public void printself() {
		System.out.println("=== AdiosFile ===");
		System.out.println("Fh: " + getFh());
		System.out.println("getNvars: " + getNvars());
		System.out.println("getNattrs: " + getNattrs());
		System.out.println("getCurrent_step: " + getCurrent_step());
		System.out.println("getLast_step: " + getLast_step());
		System.out.println("getIs_streaming: " + getIs_streaming());
		System.out.println("getPath: " + getPath());
		System.out.println("getEndianness: " + getEndianness());
		System.out.println("getFile_size: " + getFile_size());
		/*
		 * System.out.println(adioslib.get_args()); for (String vname:
		 * adioslib.get_args()) { System.out.println(vname); }
		 * adioslib.print_args(adioslib.get_args());
		 */
		System.out.println("getVar_namelist: " + getVar_namelist());
		for (String vname : getVar_namelist()) {
			// Do your stuff here
			System.out.println(vname);
			AdiosVar v = new AdiosVar(f, vname);
			v.printself();
		}

		System.out.println("getVars: " + getVars());
		for (String k : getVars().keySet()) {
			System.out.println("key: " + k + " value: " + vars.get(k));
		}
	}

	public java.math.BigInteger getFh() {
		return f.getFh();
	}

	public int getNvars() {
		return f.getNvars();
	}

	public int getNattrs() {
		return f.getNattrs();
	}

	public int getCurrent_step() {
		return f.getCurrent_step();
	}

	public int getLast_step() {
		return f.getLast_step();
	}

	public int getIs_streaming() {
		return f.getIs_streaming();
	}

	public String getPath() {
		return f.getPath();
	}

	public int getEndianness() {
		return f.getEndianness();
	}

	public int getVersion() {
		return f.getVersion();
	}

	public BigInteger getFile_size() {
		return f.getFile_size();
	}

	public String[] getVar_namelist() {
		return f.getVar_namelist();
	}

	public Map<String, AdiosVar> getVars() {
		return Collections.unmodifiableMap(vars);
	}

}
