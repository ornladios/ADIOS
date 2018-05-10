/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

/* ADIOS test: declare multiple groups and free them. 
 * Test issue that was fixed in 084d601de8fca676361b0ac4b96735e04295fa0f
 *
 * How to run: group_free_test
 * Output: None
 * ADIOS config file: None
 *
 * This is a sequential test.
*/

/* This example will create 3 adios groups for writing, then frees them.
*/
#ifndef _NOMPI
#define _NOMPI
#endif

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "public/adios.h"
#include "public/adios_types.h"
#include "core/adios_internals.h" // adios_get_groups()

const int  NBLOCKS = 3;
const int  NGROUPS = 5;
int        rank, size;
const int  NX = 100;    // local size
MPI_Comm   comm = MPI_COMM_WORLD;

int64_t create_group_byname(const char* groupname)
{
    int64_t groupid;
    char g_str[100], o_str[100], l_str[100];
    int64_t   var_ids[NBLOCKS];
    int        gb, offset;  // global size /offset
    int i;
    gb = NBLOCKS * NX * size;
    sprintf (g_str, "%d", gb);
    sprintf (l_str, "%d", NX);

    adios_declare_group (&groupid, groupname, "", adios_stat_default);
    adios_select_method (groupid, "POSIX", "", "");

    for (i = 0; i < NBLOCKS; i++)
    {
        offset = rank * NBLOCKS * NX + i * NX;
        sprintf (o_str, "%d", offset);
        var_ids[i] = adios_define_var (groupid, "temperature"
                ,"", adios_double
                ,l_str, g_str, o_str
        );
        adios_set_transform (var_ids[i], "identity");
    }

    return groupid;
}

int64_t create_group(int idx)
{
    char groupname[256];
    sprintf (groupname, "group%1.1d", idx);
    return create_group_byname(groupname);
}

int test1()
{
    printf ("Test 1\n");
    int64_t   groupid[NGROUPS];

    int j;

    for (j = 0; j < NGROUPS; j++)
    {
        groupid[j] = create_group(j);
    }

    MPI_Barrier (comm);

    // free groups in an illogical order 1, 0, 2 [3 4 ...]
    if (NGROUPS > 1) adios_free_group (groupid[1]);
    adios_free_group (groupid[0]);
    for (j = 2; j < NGROUPS; j++)
    {
        adios_free_group (groupid[j]);
    }

    // List of groups should be empty
    struct adios_group_list_struct * list = adios_get_groups();
    assert(list == NULL);
    return 0;
}

void print_groups()
{
    struct adios_group_list_struct * list = adios_get_groups();
    struct adios_group_list_struct * p = list;
    printf ("Groups = ");
    while (p)
    {
        //printf ("{addr = %p  name=%s next=%p} -> ", p, p->group->name, p->next);
        printf ("{addr = %p name=%s} -> ", p, p->group->name);
        p = p->next;
    }
    printf(" NULL\n");
}

// Create and free groups in a specific order:
// C-1, C-2, F-1, C-3, F-3, C-4, F-4, F-2
int test2()
{
    printf ("\n=======\nTest 2\n");
    int64_t groupid1 = create_group(1);
    print_groups();

    int64_t groupid2 = create_group(2);
    print_groups();

    adios_free_group (groupid1);
    print_groups();

    int64_t groupid3 = create_group(3);
    print_groups();

    adios_free_group (groupid3);
    print_groups();

    int64_t groupid4 = create_group(4);
    print_groups();

    adios_free_group (groupid4);
    print_groups();

    adios_free_group (groupid2);
    print_groups();

    // List of groups should be empty
    struct adios_group_list_struct * list  = adios_get_groups();
    assert(list == NULL);
    return 0;
}

// Create and free groups in a specific order:
// C-1, C-2, F-1, C-3, F-3, C-4, F-4, F-2
// Same as test2 but also open and close a file for each group
int test3()
{
    printf ("\n=======\nTest 3\n");
    int64_t f1, f2, f3, f4;
    int64_t groupid1 = create_group_byname("group1");
    adios_open(&f1, "group1", "group_free_test1.bp", "w", comm);

    int64_t groupid2 = create_group_byname("group2");
    adios_open(&f2, "group2", "group_free_test2.bp", "w", comm);

    adios_close(f1);
    adios_free_group (groupid1);

    int64_t groupid3 = create_group_byname("group3");
    adios_open(&f3, "group3", "group_free_test3.bp", "w", comm);

    adios_close(f3);
    adios_free_group (groupid3);

    int64_t groupid4 = create_group_byname("group4");
    adios_open(&f4, "group4", "group_free_test4.bp", "w", comm);

    adios_close(f4);
    adios_free_group (groupid4);

    adios_close(f2);
    adios_free_group (groupid2);

    // List of groups should be empty
    struct adios_group_list_struct * list = adios_get_groups();
    assert(list == NULL);

    // remove files created
    remove("group_free_test1.bp.dir/group_free_test1.bp.0");
    remove("group_free_test1.bp.dir");
    remove("group_free_test1.bp");
    remove("group_free_test2.bp.dir/group_free_test2.bp.0");
    remove("group_free_test2.bp.dir");
    remove("group_free_test2.bp");
    remove("group_free_test3.bp.dir/group_free_test3.bp.0");
    remove("group_free_test3.bp.dir");
    remove("group_free_test3.bp");
    remove("group_free_test4.bp.dir/group_free_test4.bp.0");
    remove("group_free_test4.bp.dir");
    remove("group_free_test4.bp");
    return 0;
}

int main (int argc, char ** argv)
{
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &size);

    adios_init_noxml (comm);
    adios_set_max_buffer_size (10);

    int err = test1();
    err = test2();
    err = test3();

    adios_finalize (rank);
    MPI_Finalize ();
    return 0;
}
