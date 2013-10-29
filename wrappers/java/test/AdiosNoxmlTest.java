import gov.ornl.ccs.*;

import java.nio.ByteBuffer;

public class AdiosNoxmlTest
{
    // The main program
    public static void main(String[] args)
    {
        System.out.println(">>> AdiosJava Noxml Test Drive");
        
        System.out.println(">>> ADIOS NOXML API ... ");
        Adios.MPI_Init(new String[0]);
        long comm = Adios.MPI_COMM_WORLD();
        int rank = Adios.MPI_Comm_rank(comm);
        int size = Adios.MPI_Comm_size(comm);
        System.out.println("[DEBUG] MPI rank/size = " + rank + " / " + size);
        System.out.println("[DEBUG] MPI &comm     = " + comm);

        Adios.Init_Noxml(comm);
        Adios.AllocateBuffer(AdiosBufferAllocWhen.NOW, 10);

        long group_id = Adios.DeclareGroup("restart", "iter", AdiosFlag.YES);
        Adios.SelectMethod(group_id, "MPI", "", "");
        Adios.DefineVar(group_id, "NX", "", AdiosDatatype.INTEGER, "", "", "");
        Adios.DefineVar(group_id, "G", "", AdiosDatatype.INTEGER, "", "", "");
        Adios.DefineVar(group_id, "O", "", AdiosDatatype.INTEGER, "", "", "");
        Adios.DefineVar(group_id, "temperature", "", AdiosDatatype.DOUBLE, "NX", "G", "O");

        System.out.println(">>> ADIOS Write API ... ");
        long adios_handle = Adios.Open ("restart", "adios_noxml.bp", "w", comm);

        int NX = 10; 
        int G = NX * size;
        int O = NX * rank;

        double[] t = new double[NX];
        for (int i = 0; i < NX; i++) {
            t[i] = rank * NX + (double) i;
        }

        long groupsize = 4 + 4 + 4 + 8 * (1) * (NX);
        
        long adios_totalsize = Adios.SetGroupSize(adios_handle, groupsize);
        
        Adios.Write (adios_handle, "NX", NX);
        Adios.Write (adios_handle, "G", G);
        Adios.Write (adios_handle, "O", O);
        Adios.Write (adios_handle, "temperature", t);
        Adios.Close (adios_handle);
        
        Adios.Finalize (rank);        
        Adios.MPI_Finalize();

        System.out.println(">>> Done.");
    }
}
