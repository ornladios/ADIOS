// The AdiosJava.java file
package gov.ornl.ccs;

public class AdiosFile
{

    // Declaration of the Native (C) function
    private native int adios_fopen (String path, long comm);
    private native int adios_fclose ();

    public long fp;
    public long fh;
    public int groups_count;
    public int vars_count;
    public int attrs_count;
    public int tidx_start;
    public int ntimesteps;
    public int version;
    public long file_size;
    public int endianness;
    public String[] group_namelist;

    static
    {
        // The runtime system executes a class's static
        // initializer when it loads the class.
        System.loadLibrary("AdiosJava");
    }

    public int open(String path, long comm)
    {
        return adios_fopen(path, comm);
    }

    public int close()
    {
        return adios_fclose();
    }

    public String toString()
    {
        StringBuilder result = new StringBuilder();
        String NEW_LINE = System.getProperty("line.separator");

        result.append(this.getClass().getName() + " Object {" + NEW_LINE);
        result.append(" fp: " + fp + NEW_LINE);
        result.append(" fh: " + fh + NEW_LINE);
        result.append(" groups_count: " + groups_count + NEW_LINE);
        result.append(" vars_count: " + vars_count + NEW_LINE);
        result.append(" attrs_count: " + attrs_count + NEW_LINE);
        result.append(" tidx_start: " + tidx_start + NEW_LINE);
        result.append(" ntimesteps: " + ntimesteps + NEW_LINE);
        result.append(" version: " + version + NEW_LINE);
        result.append(" file_size: " + file_size + NEW_LINE);
        result.append(" endianness: " + endianness + NEW_LINE);
        if (group_namelist != null) {
            for (int i = 0; i < groups_count; i++) {
                result.append(" group_namelist[" + i + "]: " + group_namelist[i] + NEW_LINE);
            }
        }
        result.append("}");

        return result.toString();
    }
}
