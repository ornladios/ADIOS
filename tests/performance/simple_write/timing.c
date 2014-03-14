#include <stdio.h>
#include <stdlib.h>
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

void timing_report (int nsteps)
{
    FILE *f = stdout;
    int step;
    double tcalc=0.0, tcomm=0.0, tio=0.0, tio_open=0.0, tio_group=0.0, tio_write=0.0, tio_close=0.0;
    double tiosum;

    fprintf (f, "==============================================================================\n"
                "                               TIMING   REPORT\n"
                "\n"
                "Step     Tcalc     Tcomm     Tio       Tio_open  Tio_group Tio_write Tio_close\n"
            );

    for (step = 0; step < nsteps; step++)
    {
        tiosum = Tio_open[step]+Tio_group[step]+Tio_write[step]+Tio_close[step],
        fprintf (f, "%6d    %.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
                step,
                Tcalc[step], 
                Tcomm[step],
                tiosum,
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

