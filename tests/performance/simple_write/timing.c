#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "timing.h"


// allocates and zeros out timing variables
void timing_alloc (int nsteps)
{
    // allocate time arrays
    Tcalc = (double*) calloc (nsteps, sizeof(double));
    Tcomm = (double*) calloc (nsteps, sizeof(double));
    Tio_open = (double*) calloc (nsteps, sizeof(double));
    Tio_group = (double*) calloc (nsteps, sizeof(double));
    Tio_write = (double*) calloc (nsteps, sizeof(double));
    Tio_close = (double*) calloc (nsteps, sizeof(double));
}

void timing_free ()
{
    free (Tcalc);
    free (Tcomm);
    free (Tio_open);
    free (Tio_group);
    free (Tio_write);
    free (Tio_close);
}

typedef  struct { 
    double value; 
    int    rank; 
} DBL_INT;

void timing_report (int nsteps, MPI_Comm comm)
{
    FILE *f = stdout;
    int step;
    double tcalc=0.0, tcomm=0.0, tio=0.0; // total sums */
    int rank; 
    DBL_INT *Tio;
    DBL_INT *Tio_max;

    MPI_Comm_rank (comm, &rank);

    /* Tio = SUM(Tio_*) and rank */
    Tio     = (DBL_INT*) calloc (nsteps, sizeof(DBL_INT));

    /* Tio_max[i] = maximum I/O time (among processors) in step 'i' (computed on rank 0) +
                    the rank of process which provided the maximum */
    Tio_max = (DBL_INT*) calloc (nsteps, sizeof(DBL_INT));

    for (step = 0; step < nsteps; step++) {
        Tio[step].value = Tio_open[step]+Tio_group[step]+Tio_write[step]+Tio_close[step];
        Tio[step].rank = rank;
    }

    MPI_Allreduce (Tio, Tio_max, nsteps, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    if (!rank) {
        fprintf (f, "===================================================================================\n"
                    "                               TIMING   REPORT\n"
                    "\n"
                    "Step     Tcalc     Tcomm     Tio(max)  Tio_open  Tio_group Tio_write Tio_close rank\n"
                );
        fflush(f);
    }

    for (step = 0; step < nsteps; step++)
    {
        MPI_Barrier(comm);
        /* Process with maximum Tio prints this line (lowest rank of equal Tios) */ 
        if (Tio_max[step].rank == rank) {
            fprintf (f, "%6d    %.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f     %d\n",
                    step,
                    Tcalc[step], 
                    Tcomm[step],
                    Tio_max[step].value,
                    Tio_open[step],
                    Tio_group[step],
                    Tio_write[step],
                    Tio_close[step],
                    rank
                    );
            fflush(f);
        }
        tcalc += Tcalc[step];
        tcomm += Tcomm[step];
        tio += Tio_max[step].value;
    }
            
    MPI_Barrier(comm);
    if (!rank) {
        fprintf (f, "===================================================================================\n"
                    "Total    Tcalc     Tcomm     Tio \n"
                );
        fprintf (f, "      %9.3f %9.3f %9.3f\n",
                tcalc, 
                tcomm,
                tio
                );
        fprintf (f, "Calc+comm+io = %.3f\n", tcalc+tcomm+tio);
        fprintf (f, "===================================================================================\n");
        fflush(f);
    }

}


#if 0
void timing_report_rank0 (int nsteps, double *Tio_max)
{
    FILE *f = stdout;
    int step;
    double tcalc=0.0, tcomm=0.0, tio=0.0, tio_open=0.0, tio_group=0.0, tio_write=0.0, tio_close=0.0;
    double tiosum;

    fprintf (f, "==============================================================================\n"
                "                               TIMING   REPORT\n"
                "\n"
                "Step     Tcalc     Tcomm     Tio       Tio_open  Tio_group Tio_write Tio_close\n"
                "                             max        rank 0    rank 0     rank0     rank0  \n"
            );

    for (step = 0; step < nsteps; step++)
    {
        tiosum = Tio_open[step]+Tio_group[step]+Tio_write[step]+Tio_close[step],
        fprintf (f, "%6d    %.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
                step,
                Tcalc[step], 
                Tcomm[step],
                Tio_max[step],
                Tio_open[step],
                Tio_group[step],
                Tio_write[step],
                Tio_close[step]
        );
        tcalc += Tcalc[step];
        tcomm += Tcomm[step];
        tio += tiosum;
        tio_open += Tio_open[step];
        tio_group += Tio_group[step];
        tio_write += Tio_write[step];
        tio_close += Tio_close[step];
    }
            
    fprintf (f, "==============================================================================\n"
                "Total    Tcalc     Tcomm     Tio       Tio_open  Tio_group Tio_write Tio_close\n"
            );
    tiosum = tio_open+tio_group+tio_write+tio_close;
    fprintf (f, "      %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
            tcalc, 
            tcomm,
            tiosum,
            tio_open,
            tio_group,
            tio_write,
            tio_close
            );
    fprintf (f, "Calc+comm+io = %.3f\n", tcalc+tcomm+tiosum);
    fprintf (f, "==============================================================================\n");

}
#endif
