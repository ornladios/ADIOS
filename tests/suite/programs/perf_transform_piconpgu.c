/*

   Simple cosine and sine 1D arrays over time from
   a single process.
   Test: Write plain and transformed arrays and compare timings

*/

/* Include other header files */
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <fcntl.h>
#include <stdint.h>

/* ADIOS include files */
#include "adios.h"


static int64_t       m_adios_group;
static int64_t       v_Ex,v_Ey,v_Ez,v_Bx,v_By,v_Bz;
static int64_t       v_pposx,v_pposy,v_pposz;
static int64_t       v_pmomx,v_pmomy,v_pmomz;
static int64_t       v_pcellx, v_pcelly, v_pcellz, v_pweight;


void define_vars ()
{
    adios_define_var (m_adios_group, "ngrid", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "npart", "", adios_integer, 0, 0, 0);
    adios_define_var (m_adios_group, "tmax", "", adios_real, 0, 0, 0);
    adios_define_var (m_adios_group, "dt", "", adios_real, 0, 0, 0);
    adios_define_var (m_adios_group, "t", "", adios_real, 0, 0, 0);

    v_Ex = adios_define_var (m_adios_group, "Ex", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_Ey = adios_define_var (m_adios_group, "Ey", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_Ez = adios_define_var (m_adios_group, "Ez", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_Bx = adios_define_var (m_adios_group, "Bx", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_By = adios_define_var (m_adios_group, "By", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_Bz = adios_define_var (m_adios_group, "Bz", "", adios_real, "ngrid,ngrid,ngrid", "", "");
    v_pposx     = adios_define_var (m_adios_group, "p/PositionX", "", adios_real,    "npart", "", "");
    v_pposy     = adios_define_var (m_adios_group, "p/PositionY", "", adios_real,    "npart", "", "");
    v_pposz     = adios_define_var (m_adios_group, "p/PositionZ", "", adios_real,    "npart", "", "");
    v_pmomx     = adios_define_var (m_adios_group, "p/MomentumX", "", adios_real,    "npart", "", "");
    v_pmomy     = adios_define_var (m_adios_group, "p/MomentumY", "", adios_real,    "npart", "", "");
    v_pmomz     = adios_define_var (m_adios_group, "p/MomentumZ", "", adios_real,    "npart", "", "");
    v_pcellx    = adios_define_var (m_adios_group, "p/CellX",     "", adios_integer, "npart", "", "");
    v_pcelly    = adios_define_var (m_adios_group, "p/CellY",     "", adios_integer, "npart", "", "");
    v_pcellz    = adios_define_var (m_adios_group, "p/CellZ",     "", adios_integer, "npart", "", "");
    v_pweight   = adios_define_var (m_adios_group, "p/Weight",    "", adios_real,    "npart", "", "");
}


/* --------------------------------- Main --------------------------------- */
int main( int argc, char ** argv)
{
    char        filename [256];
    MPI_Comm    comm = MPI_COMM_WORLD;
    int         rank, size;
    /* ADIOS variables declarations for matching gwrite_schema.ch */
    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;
    float       tmax = 0.5;
    float       dt = 0.5;  // run from 0.0 increasing with 'dt' up to 'tmax'
    int         step,i;
    double      tb, te; // timing
    char *transform = NULL;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    if (argc > 1) {
        transform = strndup (argv[1], 128);
        printf ("Selected transformation: %s\n", transform);
    }

    adios_init_noxml (comm);
    strcpy(filename, "perf_transform.bp");

    adios_declare_group (&m_adios_group, "transform", "iter", adios_flag_yes);
    adios_select_method (m_adios_group, "MPI", "", "");
    define_vars();

    // Declare and Initialize essential variables
    const int ngrid = 128;
    const int ngrid3 = ngrid*ngrid*ngrid;
    //const int npart = 67100000;
    const int npart = 30000000;
    float *Ex, *Ey, *Ez, *Bx, *By, *Bz;
    float *pposx, *pposy, *pposz, *pmomx, *pmomy, *pmomz, *pweight;
    int   *pcellx,  *pcelly,  *pcellz;

    Ex       = (float *) malloc (ngrid3 * sizeof(float));
    Ey       = (float *) malloc (ngrid3 * sizeof(float));
    Ez       = (float *) malloc (ngrid3 * sizeof(float));
    Bx       = (float *) malloc (ngrid3 * sizeof(float));
    By       = (float *) malloc (ngrid3 * sizeof(float));
    Bz       = (float *) malloc (ngrid3 * sizeof(float));
    pposx    = (float *) malloc (npart  * sizeof(float));
    pposy    = (float *) malloc (npart  * sizeof(float));
    pposz    = (float *) malloc (npart  * sizeof(float));
    pmomx    = (float *) malloc (npart  * sizeof(float));
    pmomy    = (float *) malloc (npart  * sizeof(float));
    pmomz    = (float *) malloc (npart  * sizeof(float));
    pweight  = (float *) malloc (npart  * sizeof(float));
    pcellx   = (int   *) malloc (npart  * sizeof(int));
    pcelly   = (int   *) malloc (npart  * sizeof(int));
    pcellz   = (int   *) malloc (npart  * sizeof(int));

    if (pcellz == NULL) {
        printf ("ERROR: could not allocate memory for all variables. Exit", step);
        MPI_Abort(comm,  1);
    }

    // Generate data once used for all timesteps
    srand (2); // fix the random sequence

    // Random data for everybody
    for (i=0; i<npart; i++) {
        pposx[i] = (float)rand()/(float)RAND_MAX;    
        pcellx[i] = rand();
    }
    memcpy (pposy, pposx, npart*sizeof(float));
    memcpy (pposz, pposx, npart*sizeof(float));
    memcpy (pmomx, pposx, npart*sizeof(float));
    memcpy (pmomy, pposx, npart*sizeof(float));
    memcpy (pmomz, pposx, npart*sizeof(float));
    memcpy (pweight, pposx, npart*sizeof(float));
    memcpy (pcelly, pcellx, npart*sizeof(int));
    memcpy (pcellz, pcellx, npart*sizeof(int));

    for (i=0; i<ngrid3; i++) {
        Ex[i] = (float)rand()/(float)RAND_MAX;    
    }
    memcpy (Ey, Ex, ngrid3*sizeof(float));
    memcpy (Ez, Ex, ngrid3*sizeof(float));
    memcpy (Bx, Ex, ngrid3*sizeof(float));
    memcpy (By, Ex, ngrid3*sizeof(float));
    memcpy (Bz, Ex, ngrid3*sizeof(float));

    //  Scan over time
    float timestep;
    step = 0;
    for (timestep = dt; timestep <= tmax; timestep = timestep + dt) {

        printf ("Step %d: ", step);
        if (step == 0) {
            adios_open (&adios_handle, "transform", filename, "w", comm);
        } else {
            adios_open (&adios_handle, "transform", filename, "a", comm);
        }

        if (transform) {
            printf (" transformation with %s\n", transform);
            adios_set_transform(v_Ex, transform);
            adios_set_transform(v_Ey, transform);
            adios_set_transform(v_Ez, transform);
            adios_set_transform(v_Bx, transform);
            adios_set_transform(v_By, transform);
            adios_set_transform(v_Bz, transform);
            adios_set_transform(v_pposx, transform);
            adios_set_transform(v_pposy, transform);
            adios_set_transform(v_pposz, transform);
            adios_set_transform(v_pmomx, transform);
            adios_set_transform(v_pmomy, transform);
            adios_set_transform(v_pmomz, transform);
            adios_set_transform(v_pweight, transform);
            adios_set_transform(v_pcellx, transform);
            adios_set_transform(v_pcelly, transform);
            adios_set_transform(v_pcellz, transform);
        } else {
            printf (" no transformation\n");
        }

        adios_groupsize = 4 + 4 + 4 \
                + 6*ngrid3*sizeof(float) \
                + 7*npart*sizeof(float) \
                + 3*npart*sizeof(int);

        if (timestep == 0 && rank == 0) {
            adios_groupsize += 4 + 4;
        }
        adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);

        tb = MPI_Wtime();
        adios_write (adios_handle, "ngrid", &ngrid);
        adios_write (adios_handle, "npart", &npart);
        adios_write (adios_handle, "t", &timestep);
        if (step == 0 && rank == 0) {
            adios_write (adios_handle, "tmax", &tmax);
            adios_write (adios_handle, "dt", &dt);
        }
        adios_write (adios_handle, "Ex", Ex);
        adios_write (adios_handle, "Ey", Ey);
        adios_write (adios_handle, "Ez", Ez);
        adios_write (adios_handle, "Bx", Bx);
        adios_write (adios_handle, "By", By);
        adios_write (adios_handle, "Bz", Bz);
        adios_write (adios_handle, "p/PositionX", pposx);
        adios_write (adios_handle, "p/PositionY", pposy);
        adios_write (adios_handle, "p/PositionZ", pposz);
        adios_write (adios_handle, "p/MomentumX", pmomx);
        adios_write (adios_handle, "p/MomentumY", pmomy);
        adios_write (adios_handle, "p/MomentumZ", pmomz);
        adios_write (adios_handle, "p/CellX", pcellx);
        adios_write (adios_handle, "p/CellY", pcelly);
        adios_write (adios_handle, "p/CellZ", pcellz);
        adios_write (adios_handle, "p/Weight", pweight);
        te = MPI_Wtime();
        printf ("   time to buffer = %.3g\n",te-tb);

        tb = MPI_Wtime();
        adios_close (adios_handle);
        te = MPI_Wtime();
        printf ("   time to close  = %.3g\n",te-tb);

        step++;
    }
   
    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();

    return 0;
}

