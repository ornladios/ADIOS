package ornl.adios;

import java.math.BigInteger;

import ornl.adios.adioslib.*;

public class AdiosVar {
	static {
		try {
			System.loadLibrary("adiosjava");
		} catch (UnsatisfiedLinkError e) {
			System.err.println("Native code library failed to load. \n" + e);
			System.exit(1);
		}
	}

	protected ADIOS_FILE f;
	protected ADIOS_VARINFO v;
	protected String name;
	protected BigInteger[] dims;

	public AdiosVar(ADIOS_FILE f, String vname) {
		this.f = f;
		name = vname;
		v = adioslib.adios_inq_var(f, vname);
		dims = new BigInteger[v.getNdim()];
		for (int i = 0; i < v.getNdim(); i++) {
			dims[i] = v.getDims(i);
		}
	}

	public void printself() {
		System.out.println("=== AdiosVar ===");
		System.out.println("Fh: " + getFh());
		System.out.println("getName: " + getName());
		System.out.println("getVarid: " + getVarid());
		System.out.println("getType: " + getType());
		System.out.println("getNdim: " + getNdim());
		System.out.print("getDims: ");
		for (int i = 0; i < dims.length; i++) {
			System.out.print(dims[i] + " ");
		}
		System.out.println();
		System.out.println("getNsteps: " + getNsteps());
		System.out.println("getGlobal: " + getGlobal());
	}

	public java.math.BigInteger getFh() {
		return f.getFh();
	}

	public String getName() {
		return name;
	}

	public int getVarid() {
		return v.getVarid();
	}

	public ADIOS_DATATYPES getType() {
		return v.getType();
	}

	public int getNdim() {
		return v.getNdim();
	}

	public java.math.BigInteger[] getDims() {
		return dims;
	}

	public int getNsteps() {
		return v.getNsteps();
	}

	public int getGlobal() {
		return v.getGlobal();
	}
}
