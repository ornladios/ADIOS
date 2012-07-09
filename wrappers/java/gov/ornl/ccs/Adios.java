// The AdiosJava.java file
package gov.ornl.ccs;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;

public class Adios
{
    // Declaration of the Native (C) function
    private static native int adios_init (String xml_fname);

    private static native long adios_open (String group_name, String file_name, String mode, long comm);
    private static native long adios_group_size (long fh, long group_size);
    private static native int adios_write (long fh, String var_name, byte val);
    private static native int adios_write (long fh, String var_name, int val);
    private static native int adios_write (long fh, String var_name, long val);
    private static native int adios_write (long fh, String var_name, float val);
    private static native int adios_write (long fh, String var_name, double val);
    private static native int adios_write (long fh, String var_name, byte[] val);
    private static native int adios_write (long fh, String var_name, int[] val);
    private static native int adios_write (long fh, String var_name, long[] val);
    private static native int adios_write (long fh, String var_name, float[] val);
    private static native int adios_write (long fh, String var_name, double[] val);
    private static native int adios_close (long fh);
    private static native int adios_finalize (int id);

    private static native int adios_mpi_init(String[] args);
    private static native int adios_mpi_comm_rank(long comm);
    private static native int adios_mpi_comm_size(long comm);
    private static native int adios_mpi_finalize();
    private static native long adios_mpi_comm_world();

    private static native long adios_open_and_set_group_size (String group_name, String file_name, String mode, long group_size, long comm);

    private static native int adios_init_noxml();
    private static native int adios_allocate_buffer(int when, long size);
    private static native long adios_declare_group(String name, String time_index, int stats);
    private static native int adios_define_var(long group_id, String name, String path, int type, String dimensions, String global_dimensions, String local_offsets);
    private static native int adios_define_attribute(long group_id, String name, String path, int type, String value, String var);
    private static native int adios_select_method(long group_id, String method, String parameters, String base_path);

    static
    {
        // The runtime system executes a class's static
        // initializer when it loads the class.
        System.loadLibrary("AdiosJava");
    }

    public static int Init(String xml_fname)
    {
        return adios_init(xml_fname);
    }
    
    public static long Open(String group_name, String file_name, String mode, long comm)
    {
        return adios_open(group_name, file_name, mode, comm);
    }

    public static long SetGroupSize(long fh, long group_size)
    {
        return adios_group_size (fh, group_size);
    }

    public static long Write (long fh, String var_name, byte value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, int value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, long value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, float value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, double value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, byte[] value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, int[] value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, long[] value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, float[] value)
    {
        return adios_write (fh, var_name, value);
    }

    public static long Write (long fh, String var_name, double[] value)
    {
        return adios_write (fh, var_name, value);
    }

    public static int Close (long fh)
    {
        return adios_close (fh);
    }

    public static int Finalize (int id)
    {
        return adios_finalize(id);
    }

    public static int MPI_Init(String[] args)
    {
        return adios_mpi_init(args);
    }

    public static int MPI_Comm_rank(long comm)
    {
        return adios_mpi_comm_rank(comm);
    }

    public static int MPI_Comm_size(long comm)
    {
        return adios_mpi_comm_size(comm);
    }

    public static int MPI_Finalize()
    {
        return adios_mpi_finalize();
    }

    public static long MPI_COMM_WORLD()
    {
        return adios_mpi_comm_world();
    }

    public static long OpenAndSetGroupSize(String group_name, String file_name, String mode, long group_size, long comm)
    {
        return adios_open_and_set_group_size(group_name, file_name, mode, group_size, comm);
    }

    public static int Init_Noxml()
    {
        return adios_init_noxml();
    }
    
    public static int AllocateBuffer(AdiosBufferAllocWhen when, long size)
    {
        return adios_allocate_buffer(when.getCode(), size);
    }
    

    public static long DeclareGroup(String name, String time_index, AdiosFlag stats)
    {
        return adios_declare_group(name, time_index, stats.getCode());
    }

    public static int DefineVar(long group_id, String name, String path, AdiosDatatype type, String dimensions, String global_dimensions, String local_offsets)
    {
        return adios_define_var(group_id, name, path, type.getCode(), dimensions, global_dimensions, local_offsets);
    }

    public static int DefineVar(long group_id, String name, String path, AdiosDatatype type, String value, String var)
    {
        return adios_define_attribute(group_id, name, path, type.getCode(), value, var);
    }

    public static int SelectMethod(long group_id, String method, String parameters, String base_path)
    {
        return adios_select_method(group_id, method, parameters, base_path);
    }
}
