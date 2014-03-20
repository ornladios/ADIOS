#include "mpi.h"

/* Code execution time separated into calculation and communication time 
 * They are filled out in main()
 */
double *Tcalc, *Tcomm;

/* IO timers: open, adios_group_size(), write and close events.
 * They should be filled out by the output_* functions 
 */
double *Tio_open, *Tio_group, *Tio_write, *Tio_close;

// allocates and zeros out timing variables
void timing_alloc (int nsteps); 
void timing_free ();
void timing_report (int nsteps, MPI_Comm comm); //write the report

