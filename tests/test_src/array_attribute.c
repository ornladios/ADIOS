/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS C Example: write/read string array and other array attributes
 *
 * How to run: ./array_attributes
 * Output: adios_global_no_xml.bp
 * ADIOS config file: None
 *
*/

/* This example will write/read string array attributes. */

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_read.h"

#ifdef DMALLOC
#include "dmalloc.h"
#endif

char  filename [] = "array_attribute.bp";
const int someints[5] = {5,4,3,2,1};
const double somedoubles[5] = {5.55555, 4.4444, 3.333, 2.22, 1.1};
char single_string[] = "A single string attribute";
char *three_strings[] = {"X","Yy","ZzZ"};
char *patchnames[] = {"arms", "deflector-bottom", "deflector-edge", "deflector-top",
        "nozzle-inner", "nozzle-outer", "outer", "inlet", "outer-top",
        "procBoundary0to1", "procBoundary0to2", "procBoundary0to5",
        "procBoundary0to6", "procBoundary0to7", "procBoundary0to8",
        "procBoundary0to13", "procBoundary0to14", "procBoundary0to185",
        "procBoundary0to68", "procBoundary0to71", "procBoundary0to72"};
char* patchtypes[] = {"wall", "wall", "wall", "wall", "wall", "wall",
        "patch", "patch", "patch", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor"};
char* U[] = {"fixedValue", "fixedValue", "fixedValue", "fixedValue",
        "fixedValue", "fixedValue", "pressureInletOutletVelocity",
        "flowRateInletVelocity", "pressureInletOutletVelocity",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor"};
char  U_class[] = "volVectorField";
char* alphawater[] = {"zeroGradient", "zeroGradient", "zeroGradient",
        "zeroGradient", "zeroGradient", "zeroGradient", "inletOutlet",
        "fixedValue", "inletOutlet", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor"};
char  alphawater_class[] = "volScalarField";
char* alphaPhi[] = {"calculated", "calculated", "calculated", "calculated",
        "calculated", "calculated", "calculated", "calculated", "calculated",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor"};
char  alphaPhi_class[] = "surfaceScalarField";
char* k[] = {"fixedValue", "fixedValue", "fixedValue", "fixedValue",
        "fixedValue", "fixedValue", "inletOutlet", "fixedValue", "inletOutlet",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor"};
char  k_class[] = "volScalarField";
char* nut[] = {"zeroGradient", "zeroGradient", "zeroGradient", "zeroGradient",
        "zeroGradient", "zeroGradient", "zeroGradient", "zeroGradient",
        "zeroGradient", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor"};
char  nut_class[] = "volScalarField";
char* p[] = {"calculated", "calculated", "calculated", "calculated",
        "calculated", "calculated", "calculated", "calculated", "calculated",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor"};
char  p_class[] = "volScalarField";
char* p_rgh[] = {"fixedFluxPressure", "fixedFluxPressure", "fixedFluxPressure",
        "fixedFluxPressure", "fixedFluxPressure", "fixedFluxPressure",
        "totalPressure", "fixedFluxPressure", "totalPressure", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor", "processor"};
char  p_rgh_class[] = "volScalarField";
char* phi[] = {"calculated", "calculated", "calculated", "calculated",
        "calculated", "calculated", "calculated", "calculated", "calculated",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor", "processor", "processor", "processor",
        "processor", "processor"};
char  phi_class[] = "surfaceScalarField";

MPI_Comm comm;
int rank, size;

int write_attrs();
int read_attrs();

int main (int argc, char ** argv) 
{

    int ret;
    comm = MPI_COMM_WORLD;
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);
    adios_init_noxml (comm);
    adios_read_init_method(ADIOS_READ_METHOD_BP, comm, "");
    adios_set_max_buffer_size (10);

    ret = write_attrs();
    if (!ret)
        ret = read_attrs();

    MPI_Barrier (comm);
    adios_finalize (rank);
    adios_read_finalize_method(ADIOS_READ_METHOD_BP);

    MPI_Finalize ();
    return ret;
}

int write_attrs()
{
    int64_t  m_adios_group;
    int64_t  m_adios_file;
    uint64_t adios_groupsize, adios_totalsize;
    double   var = 3.14159;

    adios_declare_group (&m_adios_group, "attrs", "", adios_stat_default);
    adios_select_method (m_adios_group, "POSIX", "", "");

    adios_define_var (m_adios_group, "v","", adios_double, "", "", "");

    // add some attributes
    adios_define_attribute_byvalue (m_adios_group,
                                    "single_string","", adios_string,  1, single_string);
    adios_define_attribute_byvalue (m_adios_group,
                                    "three_strings","", adios_string_array,  3, three_strings);
    adios_define_attribute_byvalue (m_adios_group,
                                    "single_int",   "", adios_integer, 1, &someints);
    adios_define_attribute_byvalue (m_adios_group,
                                    "single_double","", adios_double,  1, &somedoubles);
    adios_define_attribute_byvalue (m_adios_group,
                                    "five_ints",    "", adios_integer, 5, &someints);
    adios_define_attribute_byvalue (m_adios_group,
                                    "five_double",  "", adios_double,  5, &somedoubles);

    adios_define_attribute_byvalue (m_adios_group, "region0/patch-names","", adios_string_array,  21, patchnames);
    adios_define_attribute_byvalue (m_adios_group, "region0/patch-types","", adios_string_array,  21, patchtypes);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/U/patch-types","", adios_string_array,  21, U);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/U/class","", adios_string,  1, U_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/alpha.water/patch-types","", adios_string_array,  21, alphawater);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/alpha.water/class","", adios_string,  1, alphawater_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/alphaPhi/patch-types","", adios_string_array,  21, alphaPhi);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/alphaPhi/class","", adios_string,  1, alphaPhi_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/k/patch-types","", adios_string_array,  21, k);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/k/class","", adios_string,  1, k_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/nut/patch-types","", adios_string_array,  21, nut);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/nut/class","", adios_string,  1, nut_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/p/patch-types","", adios_string_array,  21, p);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/p/class","", adios_string,  1, p_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/p_rgh/patch-types","", adios_string_array,  21, p_rgh);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/p_rgh/class","", adios_string,  1, p_rgh_class);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/phi/patch-types","", adios_string_array,  21, phi);
    adios_define_attribute_byvalue (m_adios_group, "region0/field/phi/class","", adios_string,  1, phi_class);


    adios_open (&m_adios_file, "attrs", filename, "w", comm);
    adios_groupsize = adios_type_size(adios_double, NULL);
    adios_group_size (m_adios_file, adios_groupsize, &adios_totalsize);
    adios_write (m_adios_file, "v", &var);
    adios_close (m_adios_file);
    MPI_Barrier (comm);
    return 0;
}

void print_attr (enum ADIOS_DATATYPES attr_type, int attr_size, const void * data)
{
    int type_size = adios_type_size (attr_type, data);
    int nelems = attr_size / type_size;
    int k;
    const char *p = (const char*)data;
    for (k=0; k<nelems; k++)
    {
        if (k>0) printf(", ");
        switch (attr_type)
        {
            case adios_integer:
                printf ("%d", *(const int *)p);
                break;
            case adios_double:
                printf ("%e", *(const double *)p);
                break;
            case adios_string:
                printf ("\"%s\"", (const char *)p);
                break;
            case adios_string_array:
                printf ("\"%s\"", *(const char **)p);
                break;
            default:
                printf ("??????\n");
        }
        p=p+type_size;
    }
    printf("\n");
}

int check_attr (ADIOS_FILE *f, const char *aname, const void* orig)
{
    enum ADIOS_DATATYPES attr_type;
    int attr_size, i;
    void * data = NULL;
    int ret = adios_get_attr(f, aname, &attr_type, &attr_size, &data);
    if (ret != err_no_error)
    {
        printf ("Error getting attribute %s: %s\n", aname, adios_errmsg());
        return adios_errno;
    }
    int type_size = adios_type_size (attr_type, data);
    int nelems = attr_size / type_size;
    for (i = 0; i < nelems; ++i) {
        int different;
        if (attr_type == adios_string_array) {
            different = strcmp(  *((char**)data+i),  *((const char**)orig+i)  );
        } else {
            different = (memcmp((char*)data+i*type_size, (const char*)orig+i*type_size, type_size));
        }
        if (different)
        {
            printf ("Attribute %s element %d does not match "
                    "the original written value on %d bytes:\n\t",
                    aname, i, type_size);
            print_attr(attr_type, attr_size, data);
            printf (" != \n\t");
            print_attr(attr_type, attr_size, orig);
            switch (attr_type)
            {
                case adios_integer:
                    printf ("%d != %d\n", *(int *)data, *(const int *)orig);
                    break;
                case adios_double:
                    printf ("%g != %g\n", *(double *)data, *(const double *)orig);
                    break;
                case adios_string:
                    printf ("\"%s\" != \"%s\"\n", (char *)data, (const char *)orig);
                    break;
                case adios_string_array:
                    printf ("\"%s\" != \"%s\"\n", *(char **)data, *(const char **)orig);
                    break;
                default:
                    printf ("??????\n");
            }
            return -999;
        }
    }

    free (data);
    return err_no_error;
}

int read_attrs()
{
    int ret = 0;
    ADIOS_FILE * f = adios_read_open (filename, ADIOS_READ_METHOD_BP,
                                      comm, ADIOS_LOCKMODE_NONE, 0.0);
    if (f == NULL)
    {
        printf ("%s\n", adios_errmsg());
        ret = adios_errno;
        goto finish;
    }

    if (!rank)
    {
        printf("Found %d attributes\n", f->nattrs);
        /*
        // Print all attributes
        int j;
        enum ADIOS_DATATYPES attr_type;
        int attr_size;
        void * data = NULL;
        for (j=0; j < f->nattrs; j++)
        {
            adios_get_attr_byid (f, j, &attr_type, &attr_size, &data);
            printf ("attr: %s %s = ", adios_type_to_string(attr_type), f->attr_namelist[j]);
            print_attr (f, attr_type, attr_size, data);
            free (data);
            data = 0;
        }
        */
        ret = check_attr(f, "single_string", single_string);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "single_int", someints);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "single_double", somedoubles);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "five_ints", someints);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "five_double", somedoubles);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "three_strings", three_strings);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "region0/patch-names", patchnames);
        if (ret != err_no_error)
            goto finish;
        ret = check_attr(f, "region0/field/p_rgh/patch-types", p_rgh);
        if (ret != err_no_error)
            goto finish;

    }

finish:
    MPI_Barrier (comm);
    adios_read_close (f);
    return 0;
}


