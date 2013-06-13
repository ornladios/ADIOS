// The AdiosJava.java file
package gov.ornl.ccs;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;

public class Adios
{
    // Declaration of the Native (C) function
    private static native int adios_init (String xml_fname, long comm);

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

    // Note: adios_read is a delayed operation
    // Only non-buffered method will populate the data
    private static native byte adios_read_byte_value (long fh, String var_name);
    private static native int adios_read_int_value (long fh, String var_name);
    private static native long adios_read_long_value (long fh, String var_name);
    private static native float adios_read_float_value (long fh, String var_name);
    private static native double adios_read_double_value (long fh, String var_name);

    private static native int adios_read (long fh, String var_name, byte[] val);
    private static native int adios_read (long fh, String var_name, int[] val);
    private static native int adios_read (long fh, String var_name, long[] val);
    private static native int adios_read (long fh, String var_name, float[] val);
    private static native int adios_read (long fh, String var_name, double[] val);

    private static native int adios_close (long fh);
    private static native int adios_finalize (int id);

    private static native int adios_mpi_init(String[] args);
    private static native int adios_mpi_comm_rank(long comm);
    private static native int adios_mpi_comm_size(long comm);
    private static native int adios_mpi_finalize();
    private static native long adios_mpi_comm_world();

    private static native long adios_open_and_set_group_size (String group_name, String file_name, String mode, long group_size, long comm);

    private static native int adios_init_noxml(long comm);
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

	/* Call adios_init */
    public static int Init(String xml_fname, long comm)
    {
        return adios_init(xml_fname, comm);
    }
    
	/* Call adios_open. Return a group handler */
    public static long Open(String group_name, String file_name, String mode, long comm)
    {
        return adios_open(group_name, file_name, mode, comm);
    }

	/* Call adios_group_size and return the total size */
    public static long SetGroupSize(long fh, long group_size)
    {
        return adios_group_size (fh, group_size);
    }

	/* Call adios_write and return the total size */
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

    public static byte ReadByteValue (long fh, String var_name)
    {
        return adios_read_byte_value(fh, var_name);
    }

    public static int ReadIntValue (long fh, String var_name)
    {
        return adios_read_int_value(fh, var_name);
    }

    public static long ReadLongValue (long fh, String var_name)
    {
        return adios_read_long_value(fh, var_name);
    }

    public static float ReadFloatValue (long fh, String var_name)
    {
        return adios_read_float_value(fh, var_name);
    }

    public static double ReadDoubleValue (long fh, String var_name)
    {
        return adios_read_double_value(fh, var_name);
    }

    public static long Read (long fh, String var_name, byte[] value)
    {
        return adios_read (fh, var_name, value);
    }

    public static long Read (long fh, String var_name, int[] value)
    {
        return adios_read (fh, var_name, value);
    }

    public static long Read (long fh, String var_name, long[] value)
    {
        return adios_read (fh, var_name, value);
    }

    public static long Read (long fh, String var_name, float[] value)
    {
        return adios_read (fh, var_name, value);
    }

    public static long Read (long fh, String var_name, double[] value)
    {
        return adios_read (fh, var_name, value);
    }

	/* Call adios_close */
    public static int Close (long fh)
    {
        return adios_close (fh);
    }

	/* Call adios_finalize */
    public static int Finalize (int id)
    {
        return adios_finalize(id);
    }

	/* Call MPI_Init */
    public static int MPI_Init(String[] args)
    {
        return adios_mpi_init(args);
    }

 	/* Call MPI_Comm_rank */
    public static int MPI_Comm_rank(long comm)
    {
        return adios_mpi_comm_rank(comm);
    }

	/* Call MPI_Comm_size */
    public static int MPI_Comm_size(long comm)
    {
        return adios_mpi_comm_size(comm);
    }

	/* Call MPI_Finalize */
    public static int MPI_Finalize()
    {
        return adios_mpi_finalize();
    }

	/* Get MPI_COMM_WORLD */
    public static long MPI_COMM_WORLD()
    {
        return adios_mpi_comm_world();
    }

    public static long OpenAndSetGroupSize(String group_name, String file_name, String mode, long group_size, long comm)
    {
        return adios_open_and_set_group_size(group_name, file_name, mode, group_size, comm);
    }

	/* Call adios_init_noxml */
    public static int Init_Noxml(long comm)
    {
        return adios_init_noxml(comm);
    }
    
	/* Call adios_allocate_buffer */
    public static int AllocateBuffer(AdiosBufferAllocWhen when, long size)
    {
        return adios_allocate_buffer(when.getCode(), size);
    }
    

	/* Call adios_declare_group */
    public static long DeclareGroup(String name, String time_index, AdiosFlag stats)
    {
        return adios_declare_group(name, time_index, stats.getCode());
    }

	/* Call adios_define_var */
    public static int DefineVar(long group_id, String name, String path, AdiosDatatype type, String dimensions, String global_dimensions, String local_offsets)
    {
        return adios_define_var(group_id, name, path, type.getCode(), dimensions, global_dimensions, local_offsets);
    }

 	/* Call adios_define_attribute */
    public static int DefineAttribute(long group_id, String name, String path, AdiosDatatype type, String value, String var)
    {
        return adios_define_attribute(group_id, name, path, type.getCode(), value, var);
    }

 	/* Call adios_select_method */
    public static int SelectMethod(long group_id, String method, String parameters, String base_path)
    {
        return adios_select_method(group_id, method, parameters, base_path);
    }
}
