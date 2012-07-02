package gov.ornl.ccs;

import java.util.ArrayList;

public class AdiosGroup
{
    // Declaration of the Native (C) function
    private native int adios_gopen (long fp, String grpname);
    private native int adios_gclose ();

    AdiosFile file;
    public long gp;
    public long gh;
    public int grpid;
    public int vars_count;
    public int attrs_count;
    public String[] var_namelist;
    public String[] attr_namelist;

    static
    {
        // The runtime system executes a class's static
        // initializer when it loads the class.
        System.loadLibrary("AdiosJava");
    }

    public AdiosGroup(AdiosFile file)
    {
        this.file = file;
    }

    public int open(String grpname)
    {
        return adios_gopen(file.fp, grpname);
    }

    public int close()
    {
        return adios_gclose();
    }

    public String toString()
    {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = System.getProperty("line.separator");

        result.append(this.getClass().getName() + " Object {" + NEW_LINE);
        result.append(" gp: " + gp + NEW_LINE);
        result.append(" gh: " + gh + NEW_LINE);
        result.append(" vars_count: " + vars_count + NEW_LINE);
        result.append(" attrs_count: " + attrs_count + NEW_LINE);

        if (var_namelist != null) {
            result.append(" var_namelist.length: " + var_namelist.length + NEW_LINE);
            for (int i = 0; i < var_namelist.length; i++) {
                result.append(" var_namelist[" + i + "]: " + var_namelist[i] + NEW_LINE);
            }
        }

        if (attr_namelist != null) {
            result.append(" attr_namelist.length: " + attr_namelist.length + NEW_LINE);
            for (int i = 0; i < attr_namelist.length; i++) {
                result.append(" attr_namelist[" + i + "]: " + attr_namelist[i] + NEW_LINE);
            }
        }

        result.append("}");

        return result.toString();
    }
}