//passed the test for 1024 cores and level-3 spatial aggregation

#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// xml parser
#include <mxml.h>

// see if we have MPI or other tools
#include "config.h"

#include "public/adios.h"
#include "public/adios_types.h"
#include "public/adios_error.h"
#include "core/adios_transport_hooks.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/adios_logger.h"
#include "core/common_adios.h"

#include "mpi.h"

extern struct adios_transport_struct * adios_transports;
static int varcnt=0;
static int start_step=1;
static int total_tsteps=0; 
static char io_method[16]; //the IO method for data output
static char io_parameters[256]; //the IO method parameters 
static uint64_t totalsize=0;
static int grpflag=0; //if there's data left in buffer
char *grp_name;
int64_t grp;
static int aggr_level = 2; // currently fixed to 2 level of aggregation the most
static int aggr_chunksize; //default aggregated chunk size = 2MB
static int aggr_cnt[3][2]; //number of clients at each level for 1D, 2D and 3D variables
static int my_aggregator[3][2]; //2 level of aggregators for three dimensions 

enum ADIOS_IO_METHOD transport_method=ADIOS_METHOD_MPI;

static void aggr_chunks(void **output, int *procs, int ndims, uint64_t *ldims_list, uint64_t *gdims, uint64_t *size_list, uint64_t totalsize, int nchunks, int rank, int level, int type_size);
static uint64_t do_spatial_aggr(int level, int *procs, int ndims, uint64_t *ldims, uint64_t *offsets, char *new_ldims, int rank,  void *data, uint64_t varsize, void *output, int type_size, MPI_Comm comm);


//store the info of the client processes for spatial aggregation
struct aggr_client
{
    int rank;
    uint64_t *ldims;
};

struct aggr_var_struct
{
    char * name;
    char * path; 
    enum ADIOS_DATATYPES type;
    enum ADIOS_FLAG multidim;
    char * dimensions;
    char * global_dimensions;
    char * local_offsets;
    void * data;
    int set_aggr; //1D - 3D aggregation flags; 0:1D; 1:2D; 2:3D
//    int decomp[3]; //1D -3D decomposition 0:1D; 1:2D; 2:3D 

    struct aggr_var_struct *prev;
    struct aggr_var_struct *next;
};

struct adios_MPI_data_struct
{
    int64_t fpr;
    MPI_File fh;
    MPI_Request req;
    MPI_Status status;
    MPI_Comm group_comm;
    MPI_Info info;      // set with base
    int rank;
    int size;

    void *comm;
    struct adios_bp_buffer_struct_v1 b;
    struct adios_group_struct * group;
    char * file_mode;

    int vid;
    int64_t group_id; 
 //   int set_aggr; //1D - 3D aggregation flags; 0:1D; 1:2D; 2:3D
    int layout[3]; //1D -3D if the process layout has been determined; 
                   //XXX: we have the assumption here that all the variables with same number of dimensions have the same processes layout 
    int *procs[3]; //the proceess layout, supporting up to 3D right now
                      //e.g. nprocs[0][3]: process layout for 1D variable
    int decomp[3];
};

static struct aggr_var_struct *vars;
static struct aggr_var_struct *header;

struct aggr_client *aggr3d_clients[2]; //fixed to maximum of 2 level of aggregation
struct aggr_client *aggr2d_clients[2]; //fixed to maximum of 2 level of aggregation
struct aggr_client *aggr1d_clients[2]; //fixed to maximum of 2 level of aggregation

static uint64_t cast_var_data_as_uint64 (const char * parent_name
                                        ,enum ADIOS_DATATYPES type
                                        ,void * data
                                        )
{
    if (!data)
    {
        adios_error (err_unspecified,
                     "cannot write var since dim %s not provided\n",
                     parent_name);
        return 0;
    }

    switch (type)
    {
        case adios_byte:
            return (uint64_t) *(int8_t *) data;

        case adios_short:
            return (uint64_t) *(int16_t *) data;

        case adios_integer:
            return (uint64_t) *(int32_t *) data;

        case adios_long:
            return (uint64_t) *(int64_t *) data;

        case adios_unsigned_byte:
            return (uint64_t) *(uint8_t *) data;

        case adios_unsigned_short:
            return (uint64_t) *(uint16_t *) data;

        case adios_unsigned_integer:
            return (uint64_t) *(uint32_t *) data;

        case adios_unsigned_long:
            return (uint64_t) *(uint64_t *) data;

        case adios_real:
            return (uint64_t) *(float *) data;

        case adios_double:
            return (uint64_t) *(double *) data;

        case adios_long_double:
            return (uint64_t) *(long double *) data;

        case adios_string:
        case adios_complex:
        case adios_double_complex:
            adios_error (err_unspecified,
                         "Cannot convert type %s to integer for var %s\n",
                         adios_type_to_string_int (type), parent_name);
            return 0;
    }
    return 0;
}


static uint64_t get_value_for_dim (struct adios_file_struct * fd,
                                   struct adios_dimension_item_struct * dimension)
{
    uint64_t dim = 0;

    if (dimension->id != 0)
    {
        struct adios_var_struct * var = adios_find_var_by_id (fd->group->vars
                                                             ,dimension->id
                                                             );
        if (var)
        {
            if (var->data)
            {
                dim = cast_var_data_as_uint64 (var->name, var->type, var->data);
            }
            else
            {
                adios_error (err_dimension_required, "array dimension data missing\n");
            }
        }
        else
        {
            struct adios_attribute_struct * attr = adios_find_attribute_by_id
                                                        (fd->group->attributes
                                                        ,dimension->id
                                                        );
            if (attr)
            {
                if (attr->var)
                {
                    if (attr->var->data)
                    {
                        dim = cast_var_data_as_uint64 (attr->var->name
                                                      ,attr->var->type
                                                      ,attr->var->data
                                                      );
                    }
                    else
                    {
                        adios_error (err_dimension_required, "array dimension data missing\n");
                    }
                }
                else
                {
                    dim = cast_var_data_as_uint64 (attr->name, attr->type
                                                  ,attr->value
                                                  );
                }
            }
            else
            {
                adios_error (err_invalid_dimension,
                             "invalid dimension member id: %d\n",
                             dimension->id);
            }
        }
    }
    else
    {
        if (dimension->time_index == adios_flag_yes)
            dim = 1;
        else
            dim = dimension->rank;
    }

    return dim;
}

static uint8_t count_dimensions (struct adios_dimension_struct * dimensions)
{   
    uint8_t count = 0;

    while (dimensions)
    {
        count++;
        dimensions = dimensions->next;
    }
    
    return count;
}       

#if 0
static uint8_t check_boundary(int rank, int step)
{
    int offx, offy, offz;

    offx=rank%npx;
    offy=rank/npx%npy;
    offz=rank/(npx*npy);
        
    //it is an aggregator
    if(offx%step==0 && offy%step ==0 && offz%step ==0) 
        return 1;
    else 
        return 0;
}
#endif


//prepare the number of processes on each dimension
static int cal_layout(int *procs, int rank, int nprocs, int ndims, MPI_Comm comm, uint64_t *ldims, uint64_t *gdims, uint64_t *offsets) 
{
    char *sbuf, *recvbuf;
    int slen, recvlen, blen;
    uint64_t *t_ldims, *t_offsets;
    int i,j;
    int decomp=0;

#if 1

    for(i=0;i<ndims;i++)
        procs[i]=1;

    //printf("ndims=%d %llu %llu %llu\n", ndims, gdims[0], gdims[1], gdims[2]);
    slen=0;
    //prepare the local dimensions and offsets into send buffer
    sbuf = (char *)malloc(ndims*2*sizeof(uint64_t));
    memcpy(sbuf, ldims, ndims*sizeof(uint64_t));
    slen+=ndims*sizeof(uint64_t);
    //memcpy(sbuf+offset, gdims, ndims*sizeof(uint64_t));
    //slen+=ndims*sizeof(uint64_t);
    memcpy(sbuf+slen, offsets, ndims*sizeof(uint64_t));
    slen+=ndims*sizeof(uint64_t);
            
    recvbuf=(char *)malloc(nprocs*ndims*2*sizeof(uint64_t));
    recvlen=0;
    //rank 0 calculate the info then send to the rest 
    if(rank==0) {
        //gather all the info to rank 0
        MPI_Gather(MPI_IN_PLACE, slen, MPI_BYTE, recvbuf, slen, MPI_BYTE, 0, comm);

        t_ldims=(uint64_t *)malloc(ndims*sizeof(uint64_t));
        t_offsets=(uint64_t *)malloc(ndims*sizeof(uint64_t));

        blen=2*ndims*sizeof(uint64_t);
        for(i=1;i<nprocs;i++) {
            memcpy(t_ldims, recvbuf+blen, ndims*sizeof(uint64_t));
            blen+=ndims*sizeof(uint64_t);
            memcpy(t_offsets, recvbuf+blen, ndims*sizeof(uint64_t));
            blen+=ndims*sizeof(uint64_t);

#ifdef DEBUG
            printf("receive %d: ldims = %llu  %llu  %llu  offsets=%llu %llu %llu\n", i, t_ldims[0], t_ldims[1], t_ldims[2], t_offsets[0], t_offsets[1], t_offsets[2]);
#endif

            //the last process on the (0,0,k) dimension will be the first
            //edge process on k
            //FIXME: hard coded for 3-D 
            for(j=0;j<ndims;j++) {
                if(t_offsets[j]!=0 && t_ldims[j]+t_offsets[j]==gdims[j]) {
                    if(procs[j]==1) {
                        if(j==0) 
                            procs[j]=i+1;
                        else if(j==1) 
                            procs[j]=i/procs[0]+1; 
                        else if(j==2) { 
                            procs[j]=nprocs/(procs[0]*procs[1]);
//                            printf("ndims = %d  nprocs=%d npz=%d\n", ndims, nprocs, procs[j]);
                        }
                        decomp++;
                    }
                }
            }

#ifdef DEBUG 
    printf("%d: npx= %d  npy=%d  npz=%d\n", rank, procs[0], procs[1], procs[2]); 
#endif
        }
        free(t_ldims);
        free(t_offsets);

        //send out the process info
        slen=0;
        memset(sbuf, 0x00, (ndims+1)*sizeof(int));
        memcpy(sbuf, procs, ndims*sizeof(int));
#if 0
        memcpy(sbuf, &npx, sizeof(int));
        slen+=sizeof(int);
        memcpy(sbuf+slen, &npy, sizeof(int));
        slen+=sizeof(int);
        memcpy(sbuf+slen, &npz, sizeof(int));
        slen+=sizeof(int);
#endif
        slen+=ndims*sizeof(int);
        memcpy(sbuf+slen, &decomp, sizeof(int));
        slen+=sizeof(int);

        MPI_Bcast(sbuf, slen, MPI_BYTE, 0, comm);
    }
    else {
        MPI_Gather(sbuf, slen, MPI_BYTE, recvbuf, slen, MPI_BYTE, 0, comm);

        //receive npx, npy, npz from rank 0
        memset(sbuf, 0x00, (ndims+1)*sizeof(int)); 

        slen=(ndims+1)*sizeof(int);
        MPI_Bcast(sbuf, slen, MPI_BYTE, 0, comm);
#if 0
        //XXX: hard coded
        memcpy(&npx, sbuf, sizeof(int));
        memcpy(&npy, sbuf+sizeof(int), sizeof(int));
        memcpy(&npz, sbuf+2*sizeof(int), sizeof(int));
#endif
        memcpy(procs, sbuf, ndims*sizeof(int));
        memcpy(&decomp, sbuf+ndims*sizeof(int), sizeof(int));
    }

    /*npx=procs[0];
    npy=procs[1];
    npz=procs[2];
    */

#ifdef DEBUG
    //printf("%d: received npx= %d  npy=%d  npz=%d\n", rank, npx, npy, npz);
    printf("%d: received npx= %d  npy=%d  npz=%d\n", rank, procs[0], procs[1], procs[2]); 
#endif

#endif
    free(sbuf);
    free(recvbuf);

    return decomp;
} 


static void cal_offsets(int *procs, int rank, int ndims, int decomp, int *offsets)
{
    if(ndims==1) {
        offsets[0]=rank%procs[0];
        return;
    }
    else if(ndims==2) {
        //k or j
        if(decomp==1) {
            if(procs[1]==1) 
                offsets[0]=rank%procs[0];
            else if(procs[0]==1) 
                offsets[1]=rank%procs[1];
        }
        else if(decomp==2) {
            offsets[0]=rank%procs[0];
            offsets[1]=rank/procs[0]%procs[1];
        }
    }
    else {
        if(decomp==1) {
            //i, j, k
            if(procs[0]==1 && procs[1]==1)
                offsets[2]=rank%procs[2];
            else if(procs[1]==1 && procs[2]==1)
                offsets[0]=rank%procs[0];
            else if(procs[0]==1 && procs[2]==1)
                offsets[1]=rank%procs[1];
        }
        else if(decomp==2) {
            //(j,k), (i,k), (i,j)
            if(procs[0]==1) {
                offsets[1]=rank%procs[1];
                offsets[2]=rank/procs[1]%procs[2];
            }
            else if(procs[1]==1) {
                offsets[0]=rank%procs[0];
                offsets[2]=rank/procs[0]%procs[2];
            }
            else {
                offsets[0]=rank%procs[0];
                offsets[1]=rank/procs[0]%procs[1];
            }
        }
        else {
            offsets[0]=rank%procs[0];
            offsets[1]=rank/procs[0]%procs[1];
            offsets[2]=rank/(procs[0]*procs[1]);
        }
    }

}

//prepare the aggregation group
//XXX: not yet able to handle odd number of process on one dimension
// rank: the rank of MPI process;
// level: the level of aggregation
static void prep_aggr(int *procs, int ndims, int decomp, int rank, int size, int level)
{
    int scale=1, step=2;
    int aggrx, aggry, aggrz;
    int aggr=1;
    int i,j,k;
    int prev_step, hole;
    int *offsets;

#if 1

    offsets=(int *)malloc(ndims*sizeof(int));
    memset(offsets,0x00, ndims*sizeof(int));

    cal_offsets(procs, rank, ndims, decomp, offsets);

    aggr_cnt[ndims-1][0]=aggr_cnt[ndims-1][1]=1;
    prev_step=1;
    for(i=1; i<=level; i++) {
        scale=(int)pow(2, (i-1));
        step=(int)pow(2, i);
        
#ifdef DEBUG 
        if(rank==0) 
            printf("level=%d :  step = %d  scale=%d\n", i-1, step, scale);
#endif

        //detemine the aggregators and clients 
        hole=0;
        for(j=0;j<ndims;j++) {
            if(offsets[j]%step!=0){
                hole=1;
                break;
            }
        }
        if(hole==0) {//I'am aggregator
//            printf("Aggregator %d: %d %d %d\n", rank, offsets[0], offsets[1], offsets[2]); 

//            if(offsets[0]%step==0 && offsets[1]%step ==0 && offsets[2]%step ==0) { //I'am aggregator
            //printf("%d: i am aggregator %d: %d %d %d\n", i-1, rank, offx, offy, offz);
            my_aggregator[ndims-1][i-1]=rank;

            //allocate the space for the list of clients
            int mal_size=(int)pow(2, decomp)-1;
            if(ndims==3){
                //aggr_clients[i-1] = malloc(mal_size*sizeof(int)); 
                aggr3d_clients[i-1] = malloc(mal_size*sizeof(struct aggr_client)); 
                memset(aggr3d_clients[i-1], 0x00, mal_size*sizeof(struct aggr_client));
            }
            else if (ndims==2) {
                aggr2d_clients[i-1] = malloc(mal_size*sizeof(struct aggr_client)); 
                memset(aggr2d_clients[i-1], 0x00, mal_size*sizeof(struct aggr_client));
            }
            else if(ndims==1){
                aggr1d_clients[i-1] = malloc(mal_size*sizeof(struct aggr_client)); 
                memset(aggr1d_clients[i-1], 0x00, mal_size*sizeof(struct aggr_client));
            }

            //3D variable
            if(ndims==3){
            /*FIXME XXX: hard coded for calculating aggr_clients*/
            //if(my_aggregator[i-1] == rank) {
                     aggr_cnt[ndims-1][i-1]=0;
                     //check if the clients'rank is valid on each dimension
                     if((offsets[0]+scale)<procs[0]) {
                         aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+1*scale;
                         aggr_cnt[ndims-1][i-1]++;
                     }

                     if((offsets[1]+scale)<procs[1]) {
                         aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[0]*scale;
                         aggr_cnt[ndims-1][i-1]++;
                         if((offsets[0]+scale)<procs[0]) {
                             aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[0]*scale+1*scale;
                             aggr_cnt[ndims-1][i-1]++;
                         }
                     }
                     
                     if((offsets[2]+scale)<procs[2]) {
                        aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[1]*scale*procs[0];
                        aggr_cnt[ndims-1][i-1]++;
                        if((offsets[0]+scale)<procs[0]) {
                            aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[1]*scale*procs[0]+1*scale;
                            aggr_cnt[ndims-1][i-1]++;
                        }
                        if((offsets[1]+scale)<procs[1]) {
                            aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[1]*scale*procs[0]+procs[0]*scale;
                            aggr_cnt[ndims-1][i-1]++;
                            if((offsets[0]+scale)<procs[0]) {
                                aggr3d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[1]*scale*procs[0]+procs[0]*scale+1*scale;
                                aggr_cnt[ndims-1][i-1]++;
                            }
                        }
                     //}
                 }
            }
            else if(ndims==2){ //2D variable
                aggr_cnt[ndims-1][i-1]=0;
                //check if the clients'rank is valid on each dimension
                if((offsets[0]+scale)<procs[0]) {
                    aggr2d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+1*scale;
                    aggr_cnt[ndims-1][i-1]++;
                }

                if((offsets[1]+scale)<procs[1]) {
                    aggr2d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[0]*scale;
                    aggr_cnt[ndims-1][i-1]++;
                    if((offsets[0]+scale)<procs[0]) {
                        aggr2d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+procs[0]*scale+1*scale;
                        aggr_cnt[ndims-1][i-1]++;
                    }
                }
            }
            else { //1D variable
                aggr_cnt[ndims-1][i-1]=0;
                //check if the clients'rank is valid on each dimension
                if((offsets[0]+scale)<procs[0]) {
                    aggr1d_clients[i-1][aggr_cnt[ndims-1][i-1]].rank=rank+1*scale;
                    aggr_cnt[ndims-1][i-1]++;
                }
            }
#ifdef DEBUG 
            printf("%d: i am aggregator, level %d # of client = %d \n", rank, i, aggr_cnt[ndims-1][i-1]); 
            int m;
            for(m=0;m<aggr_cnt[ndims-1][i-1];m++) {
                if(ndims==1)
                    printf("%d  ", aggr1d_clients[i-1][m].rank);
                if(ndims==2)
                    printf("%d  ", aggr2d_clients[i-1][m].rank);
                if(ndims==3)
                    printf("%d  ", aggr3d_clients[i-1][m].rank);
            }
            printf("\n");
#endif
         }
         else { //I am the clients
             aggrx=aggry=aggrz=0;
             aggrx=offsets[0]-offsets[0]%step;
             if(ndims>=2) 
                aggry=offsets[1]-offsets[1]%step;
             if(ndims>=3) 
                 aggrz=offsets[2]-offsets[2]%step;
             my_aggregator[ndims-1][i-1] = aggrz*procs[0]*procs[1]+aggry*procs[0]+aggrx;
                 
 //            printf("[%d] aggr offset= %d %d %d\n", rank, aggrx, aggry, aggrz);
 
             //check if this process needs to be included within the  
             //communication of this level
             if(i>1) {
                 hole=0;
                 for(k=0;k<ndims;k++) {
                     if(offsets[k]%prev_step!=0) {
                         hole=1;
                         break;
                     }
                 }
                 if(hole==1)
                     my_aggregator[ndims-1][i-1]=-1;
             }
#ifdef DEBUG 
             printf("%d: level=%d  my aggregator = %d\n", rank, i, my_aggregator[ndims-1][i-1]);
#endif
         } //end of if(hole==0)
         prev_step=step; 
    }//end of for()

#endif
#ifdef DEBUG
    if(level>=1) {
        if(my_aggregator[ndims-1][0]==rank)
            printf("%d: i am aggregator, client = %d %d %d\n", rank, aggr3d_clients[0], aggr3d_clients[1], aggr3d_clients[2]);
        else
            printf("level %d:  [%d] = %d %d %d, myaggr=[%d]\n", level, rank, offz, offy, offx, my_aggregator[ndims-1][0]);
    }
#endif
}
 
 

static int do_write (int64_t fd_p, const char * name, void * var)
{
    struct adios_file_struct * fd = (struct adios_file_struct *) fd_p;

    if (!fd)
    {
        adios_error (err_invalid_file_pointer, "Invalid handle passed to adios_write\n");
        return 1;
    }

    struct adios_var_struct * v = fd->group->vars;
    struct adios_method_list_struct * m = fd->group->methods;

    if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL)
    {
        // nothing to do so just return
        return 0;
    }

    v = adios_find_var_by_name (v, name, fd->group->all_unique_var_names);

    if (!v)
    {
        adios_error (err_invalid_varname, "Bad var name (ignored) in adios_write(): '%s'\n", name);

        return 1;
    }

    common_adios_write_byid (fd, v, var);

    return 0;
}
                       

// temporary solution for compiling error
int declare_group (int64_t * id, const char * name
                        ,const char * time_index
                        ,enum ADIOS_FLAG stats
                        )
{
    int ret;
    ret = adios_common_declare_group (id, name, adios_flag_yes
                                      ,""
                                      ,""
                                      ,time_index
                                      ,stats
                                      );
    if (ret == 1) {
        struct adios_group_struct * g = (struct adios_group_struct *) *id;
        g->all_unique_var_names = adios_flag_no;
    }
    return ret;
}

// temporary solution for compiling error
int select_method (int64_t group, const char * method
                        ,const char * parameters
                        ,const char * base_path
                        )
{
    return adios_common_select_method_by_group_id (0, method, parameters, group ,base_path, 0);
}

static int convert_file_mode(enum ADIOS_METHOD_MODE mode, char * file_mode)
{
   if (mode == adios_mode_read)
       strcpy(file_mode,"r");
   else
       if (mode == adios_mode_write)
       strcpy(file_mode,"w");
       else
           if (mode == adios_mode_append)
       strcpy(file_mode,"a");
           else
               if (mode == adios_mode_update) 
       strcpy(file_mode,"u");
               else
               {
                   fprintf (stderr, "adios_open: unknown file mode: %s\n"
                           ,file_mode
                           );

                   return -1;
               }
   return 0;
}

//find the variable within the buffer link list
int var_lookup(const char *varname, char *path, struct aggr_var_struct *list)
{
    int cnt=0;

    for(cnt=0;cnt<varcnt;cnt++) 
    {
        //compare both the variable name and path
        if(strcmp(varname, list->name)==0 && strcmp(path, list->path)==0) {
            return cnt;
        }
        else if(list->next!=NULL){
            list=list->next;
        }
        else
            break;
    }

    //variable is not within the list, return -1
    return -1;
}

void output_vars(struct aggr_var_struct *vars, int varcnt, struct
        adios_MPI_data_struct * md, struct adios_file_struct * fd) 
{
    int i,j;
    char file_mode[1];
    char fname[256];
    uint64_t adios_size, datasize;
    int iocnt=0;
    char *output;

    if(convert_file_mode(fd->mode, file_mode) == -1) //strange file mode
        return;

    common_adios_open (&md->fpr, grp_name, fd->name, file_mode, md->group_comm);
    common_adios_group_size (md->fpr, totalsize, &adios_size);

    //move pointer to the first variable in the list
    vars=header; 
    //write it out
    for(i=0;i<varcnt;i++) {
/*        if(i!=13) {
            int tmp=0;
            memcpy(&tmp, vars->data, 4);
            printf("i = %d,  var=%s   value=%d\n", i, vars->name, tmp);
        }*/

        do_write(md->fpr, vars->name, vars->data);
        //adios_write(md->fpr, vars->name, vars->data);
        vars=vars->next;
    }
    //close the file
    common_adios_close(md->fpr);
}


void define_iogroup(char *group_name)
{
    int len;

    // is it necessary to have different group name? XXX:FIXME 
    len=5+strlen(group_name); //new groupname= tg_groupname
    grp_name=(char *)malloc(len);
    memset(grp_name, 0x00, len);
    sprintf(grp_name, "agg_%s",group_name);
    declare_group (&grp,grp_name, "", adios_flag_no);
    select_method (grp, io_method,io_parameters,"");
    grpflag=1;
}

//initial variable structure
void init_vars(struct aggr_var_struct *var, struct adios_var_struct * v, int ndims)
{
    int i;

    vars->name=(char *)malloc(strlen(v->name)+1);
    strcpy(vars->name, v->name); 
    vars->type=v->type;
    vars->path=(char *)malloc(strlen(v->path)+1);
    strcpy(vars->path, v->path);
    vars->next=NULL;
    vars->dimensions = (char *)malloc(128*sizeof(char));
    vars->global_dimensions= (char *)malloc(128*sizeof(char));
    vars->local_offsets= (char *)malloc(128*sizeof(char));
    memset(vars->dimensions, 0x00, 128*sizeof(char));
    memset(vars->global_dimensions, 0x00, 128*sizeof(char));
    memset(vars->local_offsets, 0x00, 128*sizeof(char));
    /*
    vars->ldims= (uint64_t *)malloc((ndims+1)*sizeof(uint64_t));
    vars->gdims= (uint64_t *)malloc((ndims+1)*sizeof(uint64_t));
    vars->offsets= (uint64_t *)malloc((ndims+1)*sizeof(uint64_t)); 
    vars->ndims=0; */
    //vars->count=1;
    //vars->iosize=0;
    //vars->datasize=0;
    //vars->flag=0;

    for(i=0;i<3;i++) {
        vars->set_aggr=-1;
    //    vars->decomp[i]=0;
    }
}

static void init_layout_flag(struct adios_MPI_data_struct *md)
{
    int i;

    for(i=0;i<3;i++)
        md->layout[i]=-1;
}

void adios_chunk_init (const PairStruct * parameters, 
                     struct adios_method_struct * method)
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

    method->method_data = malloc (sizeof (struct adios_MPI_data_struct));
    md = (struct adios_MPI_data_struct *) method->method_data;

    init_layout_flag(md);
}


static void adios_var_to_comm (const char * comm_name
                              ,enum ADIOS_FLAG host_language_fortran
                              ,void * data
                              ,MPI_Comm * comm
                              )
{
    if (data)
    {
        int t = *(int *) data;

        if (!comm_name)
        {
            if (!t)
            {
                log_warn ("TIAN method: communicator not provided and none "
                          "listed in XML.  Defaulting to MPI_COMM_SELF\n");

                *comm = MPI_COMM_SELF;
            }
            else
            {
                if (host_language_fortran == adios_flag_yes)
                {
                    *comm = MPI_Comm_f2c (t);
                }
                else
                {
                    *comm = *(MPI_Comm *) data;
                }
            }
        }
        else
        {
            if (!strcmp (comm_name, ""))
            {
                if (!t)
                {
                    log_warn ("TIAN method: communicator not provided and none "
                              "listed in XML.  Defaulting to MPI_COMM_SELF\n");

                    *comm = MPI_COMM_SELF;
                }
                 else
                 {
                     if (host_language_fortran == adios_flag_yes)
                     {
                         *comm = MPI_Comm_f2c (t);
                     }
                     else
                     {
                         *comm = *(MPI_Comm *) data;
                     }
                 }
             }
             else
             {
                 if (!t)
                 {
                     log_warn ("TIAN method: communicator not provided but one "
                               "listed in XML.  Defaulting to MPI_COMM_WORLD\n");
 
                     *comm = MPI_COMM_WORLD;
                 }
                 else
                 {
                     if (host_language_fortran == adios_flag_yes)
                     {
                         *comm = MPI_Comm_f2c (t);
                     }
                     else
                     {
                         *comm = *(MPI_Comm *) data;
                     }
                 }
             }
         }
     }
     else
     {
         log_warn ("TIAN method: coordination-communication not provided. "
                   "Using MPI_COMM_WORLD instead\n");
 
         *comm = MPI_COMM_WORLD;
     }
}

void init_output_parameters(const char *parameters)
{
    int len;
    char *temp_string, *p_size;

    // parsing XML to get the number of timesteps 
    temp_string = (char *) malloc (strlen(parameters) + 1);

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);
  
    //set the aggregated chunk size
    strcpy (temp_string, parameters);
    trim_spaces (temp_string);
    if (p_size = strstr (temp_string, "chunk_size"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ",");
        if (!q) 
            aggr_chunksize= atoi(q+1);
        else
            aggr_chunksize= atoi(p+1); 
    }
    else
    {
        aggr_chunksize=1048576*2;
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);
    if (p_size = strstr (temp_string, "io_method"))
    {
        char * p = strchr (p_size, '=');
        char * q = strtok (p, ",");
        memset(io_method, 16, 0x00);
        if (!q) 
            strcpy(io_method, q+1);
        else
            strcpy(io_method, p+1);
    }
    else
    {
        //without specifing in the XML, no time buffering
        memset(io_method, 16, 0x00);
        strcpy(io_method, "MPI"); 
    }

    strcpy (temp_string, parameters);
    trim_spaces (temp_string);
    if (p_size = strstr (temp_string, "parameters"))
    {
        char * p = strchr (p_size, ':');
        char * q = strtok (p, "<");
        memset(io_parameters, 256, 0x00);
        if (!q) 
            strcpy(io_parameters, q+1);
        else
            strcpy(io_parameters, p+1);
    }
    else
    {
        //without specifing in the XML, no time buffering
        memset(io_parameters, 256, 0x00);
    }
}

int adios_chunk_open (struct adios_file_struct * fd
                   ,struct adios_method_struct * method, MPI_Comm comm)
{

    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;

#if 0
 #if COLLECT_METRICS
     gettimeofday (&timing.t0, NULL); // only used on rank == size - 1, but we don't
                               // have the comm yet to get the rank/size
 #endif
#endif
// adios_buffer_struct_clear (&md->b);
 
#if 0
     adios_var_to_comm (fd->group->group_comm
                       ,fd->group->adios_host_language_fortran
                       ,comm
                       ,&md->group_comm
                       );
#endif


    md->group_comm = comm;
    if (md->group_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank (md->group_comm, &md->rank);
        MPI_Comm_size (md->group_comm, &md->size);
    }
    fd->group->process_id = md->rank;

    //need to get the parameters form XML
     init_output_parameters(method->parameters);
//     md->comm=comm;

#if 0
#if COLLECT_METRICS
    timing.write_count = 0;
    timing.write_size = 0;
    if (timing.t24) free (timing.t24);
    timing.t24 = 0;
#endif
#endif

    return 1;
}

enum ADIOS_FLAG adios_chunk_should_buffer (struct adios_file_struct * fd
                                         ,struct adios_method_struct * method)
{

    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    define_iogroup(method->group->name); 

    //this method handles its own buffering
    return adios_flag_no;
}
                                         
static void prepare_data(void **data, uint64_t varsize, int dims)
{
    int i,j,k,l;
    double val;
    
    l=0;
    for(val=0;val<8;val++) {
        for(i=0;i<dims;i++) {
            for(j=0;j<dims;j++) {
                for(k=0;k<dims;k++) {
                    memcpy(*data+l*sizeof(double), &val,sizeof(double));
                    l++;
                }
            }
        }
    }
}

static struct aggr_var_struct *allocate_vars(int varcnt, struct aggr_var_struct *vars)
{
     if(varcnt==0) {
         vars = (struct aggr_var_struct*) malloc (sizeof(struct aggr_var_struct));
         header=vars; //assign the header of the variable list
     }
     else {
         vars->next = (struct aggr_var_struct*) malloc (sizeof(struct aggr_var_struct));
         vars=vars->next;
     }

     return vars;
}

void adios_chunk_write (struct adios_file_struct * fd
                     ,struct adios_var_struct * v
                     ,void * data
                     ,struct adios_method_struct * method
                     )

{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
    struct adios_dimension_struct * d = v->dimensions;
    struct aggr_var_struct *tmp;
    uint8_t dims_count = 0;
    uint64_t total_size, varsize, alloc_size;
    int64_t afd;
    int i, ndims,type_size=0;
    uint64_t *ldims, *offsets, *gdims;
    char *new_ldims, *dimensions, *local_offsets, *global_dimensions;
    int chunk_cnt, decomp=0;

    ndims=0;

//    tmp=(struct adios_MPI_data_struct *) malloc (sizeof(struct adios_MPI_data_struct));
    if(varcnt==0) {
#if 0
//XXX: write attributes
        if(md->rank==0) {
            struct adios_attribute_struct * a = fd->group->attributes;
                adios_common_define_attribute(grp, a->name, a->path, a->type, (char *)a->value, ""); 
                a=a->next;
            }
        }
#endif
        vars = (struct aggr_var_struct *) malloc (sizeof(struct aggr_var_struct));
        vars->prev=NULL;
        header=vars; //assign the header of the variable list
    }
    else {
        tmp=vars;
        vars->next = (struct aggr_var_struct*) malloc (sizeof(struct aggr_var_struct));
        vars=vars->next;
        vars->prev=tmp;
    }

    //initial the variable structure
    init_vars(vars, v, ndims);

    //retrieve the chunk size
    varsize=adios_get_var_size(v, method->group, data); 

    //number of the dimensions of this variable
    ndims=count_dimensions(v->dimensions);
    if(ndims) //multidimensional data
    {
        vars->multidim=adios_flag_yes;
           
        ldims=(uint64_t *)malloc(ndims*sizeof(uint64_t));
        offsets=(uint64_t *)malloc(ndims*sizeof(uint64_t));
        gdims=(uint64_t *)malloc(ndims*sizeof(uint64_t));
        while (d) {
            uint64_t dim = 0;
            //local dimension
            dim = get_value_for_dim (fd, &d->dimension);
            ldims[dims_count]=dim;
            if(dims_count==0)
                sprintf(vars->dimensions,"%llu", dim);
            else
                sprintf(vars->dimensions,"%s,%llu",vars->dimensions, dim);

            //global dimension
            dim = get_value_for_dim (fd, &d->global_dimension);
            gdims[dims_count]=dim;
            if(dims_count==0)
                sprintf(vars->global_dimensions,"%llu", dim);
            else
                sprintf(vars->global_dimensions,"%s,%llu",vars->global_dimensions, dim);

            //local offsets
            dim = get_value_for_dim (fd, &d->local_offset);
            offsets[dims_count]=dim;
            if(dims_count==0)
                sprintf(vars->local_offsets,"%llu", dim);
            else
                sprintf(vars->local_offsets,"%s,%llu",vars->local_offsets, dim);
            
            dims_count++;
            d=d->next;
        } //end of while (d) 


        //determine if we need to apply spatial aggregation
        //first find out the process layout and domain decomposition
        if(ndims<=3 && varsize<=aggr_chunksize && md->size>1 && vars->set_aggr==-1) {
            if(md->layout[ndims-1]==-1) {
                md->procs[ndims-1]=(int *)malloc(ndims*sizeof(int)); //FIXME
                decomp=cal_layout(md->procs[ndims-1], md->rank, md->size, ndims, md->group_comm, ldims, gdims, offsets);
                md->layout[ndims-1]=1;
                md->decomp[ndims-1]=decomp;
            }
            else
                decomp=md->decomp[ndims-1];
#ifdef DEBUG 
            if(md->rank==0) 
                printf("dims=%d decomp = %d npx=%d  npy=%d  npz=%d\n", ndims, decomp, md->procs[ndims-1][0], md->procs[ndims-1][1], md->procs[ndims-1][2]);
#endif
            if(decomp==0) {
                //FIXME: need to fix the error message
                adios_error(err_corrupted_variable, "Wrong decomposition.");
                exit(-1);
            }
            else {
                //the number of chunks to aggregate in one level of spatial aggregation 
                //is decided by the domain decomposition
                chunk_cnt=(int)pow(2, decomp);

                //too few process or chunk size is large enough, we do not need to aggregate
                //XXX: default the aggregated chunksize=2MB
                if (md->size < chunk_cnt || varsize > aggr_chunksize/chunk_cnt)  
                    vars->set_aggr=0;
                else {
                    vars->set_aggr=1;

                    aggr_level=varsize/chunk_cnt;
                    //XXX: currently the maximum level is fixed to 2 considering the overhead
                    if(aggr_level>2) {
                        // we need at least twice the number of chunks for
                        // the first level aggregation in order to do higher
                        // level aggregation
                        if(md->size<chunk_cnt*2) 
                            aggr_level=1;
                        else
                            aggr_level=2; 
                    }

                    //calculating the aggregator and client processes
                    prep_aggr(md->procs[ndims-1], ndims, decomp, md->rank, md->size, aggr_level);
                }
                //printf("aggr chunksize=%llu,  chunk_cnt=%d varsize=%llu\n", aggr_chunksize, chunk_cnt, varsize);
            }
        }

        
        //no spatial aggregation, just copy data
        if(vars->set_aggr!=1) {
            vars->data=malloc(varsize);
            memcpy(vars->data, data, varsize);
        }
        else { //if we need to do spatial aggregation 
            //only the highest level aggregators need to allocate space for output 
            if(my_aggregator[ndims-1][aggr_level-1]==md->rank) {
                //allocate the total buffer for the output data
                alloc_size=varsize;
                for(i=0;i<aggr_level;i++)
                    alloc_size*=(aggr_cnt[ndims-1][i]+1);
                vars->data=malloc(alloc_size);
            } 

            new_ldims=(char *)malloc(128*sizeof(char));
            memset(new_ldims, 0x00, 128);
            type_size=adios_get_type_size(v->type,data);
            varsize=do_spatial_aggr(aggr_level, md->procs[ndims-1], ndims, ldims, offsets, new_ldims, md->rank, data, varsize, vars->data, type_size, md->group_comm);

            //only the highest level aggregators need to output
            if(my_aggregator[ndims-1][aggr_level-1]==md->rank) {
            //    printf("%d: return var size = %llu\n", md->rank, varsize, new_ldims);
                strcpy(vars->dimensions, new_ldims);
                varcnt++;
            }
            else //clients and lower level aggregators skip the variable
              varsize=0;
        } //end of if(do_spatial_aggr)
        //free(ldims);
    } //end of if(ndims)
    else //scalar 
    {
     //   vars=allocate_vars(varcnt, vars);
        vars->multidim=adios_flag_no;
        varsize=adios_get_type_size(v->type,data);
        vars->data=malloc(varsize);
        memcpy(vars->data, data, varsize); 
        varcnt++;
    }

    totalsize+=varsize;
    if(varsize>0) {
        adios_common_define_var(grp, vars->name, vars->path, vars->type, vars->dimensions, vars->global_dimensions, vars->local_offsets);
    }
    else { //move back the pointer, and release the memory
        vars=vars->prev;
        free(vars->next);
    }
}

void adios_chunk_read (struct adios_file_struct * fd
                    ,struct adios_var_struct * v, void * buffer
                    ,uint64_t buffer_size
                    ,struct adios_method_struct * method
                    )

{
}

void release_resource()
{ 
    int cnt;
    for(cnt=0;cnt<varcnt;cnt++)
    {
        free(vars->data);
        free(vars->dimensions);
        free(vars->global_dimensions);
        free(vars->local_offsets);
    }
}

void adios_chunk_close (struct adios_file_struct * fd
                     ,struct adios_method_struct * method
                     )
{
    struct adios_MPI_data_struct * md = (struct adios_MPI_data_struct *)
                                                    method->method_data;
#ifdef DEBUG
    printf("%d: adios_close = %d\n", rank, varcnt);
#endif
    //write out the varaibles
    output_vars(header, varcnt, md, fd);
    //clean the counters
    varcnt=0;

    return;
}

void adios_chunk_get_write_buffer (struct adios_file_struct * fd
                                ,struct adios_var_struct * v
                                ,uint64_t * size
                                ,void ** buffer
                                ,struct adios_method_struct * method
                                )
{
}

void adios_chunk_finalize (int mype, struct adios_method_struct * method)
{
}

void adios_chunk_end_iteration (struct adios_method_struct * method)
{
}

void adios_chunk_start_calculation (struct adios_method_struct * method)
{
}

void adios_chunk_stop_calculation (struct adios_method_struct * method)
{
}

void copy_aggr_data (void *dst, void *src,
        int idim,
        int ndim,
        uint64_t* size_in_dset,
        uint64_t* ldims,
        const uint64_t * readsize,
        uint64_t dst_stride, //8 
        uint64_t src_stride, //9
        uint64_t dst_offset,
        uint64_t src_offset,
        uint64_t ele_num,
        int      size_of_type,
        int rank
        )
{
    unsigned int i, j;
    uint64_t dst_offset_new=0;
    uint64_t src_offset_new=0;
    uint64_t src_step, dst_step;


    if (ndim-1==idim) {
        for (i=0;i<size_in_dset[idim];i++) {
#ifdef DEBUG 
            printf("from: %d  copy to: %d ele_num=%d\n",i*src_stride+src_offset, i*dst_stride+dst_offset, ele_num);
#endif
            memcpy ((char *)dst + (i*dst_stride+dst_offset)*size_of_type,
                    (char *)src + (i*src_stride+src_offset)*size_of_type,
                    ele_num*size_of_type);
        }
        return;
    }

    for (i = 0; i<size_in_dset[idim];i++) {
        // get the different step granularity 
        // for each different reading pattern broke
        src_step = 1;
        dst_step = 1;
        for (j = idim+1; j <= ndim-1;j++) {
            src_step *= ldims[j];
            dst_step *= readsize[j];
        }
        src_offset_new  =src_offset + i * src_stride * src_step;
        dst_offset_new  = dst_offset + i * dst_stride * dst_step;

#ifdef DEBUG 
printf("src_offset =%d, dst-offset= %d\n", src_offset_new, dst_offset_new);
        printf("src_step=%d  dst_step=%d src_off=%d   dst_off=%d\n",\
                src_step, dst_step, src_offset_new, dst_offset_new);
#endif
        copy_aggr_data (dst, src, idim+1, ndim, size_in_dset, ldims,readsize, dst_stride, src_stride, dst_offset_new, src_offset_new, ele_num, size_of_type, rank);
    }
}

static void cal_gdims(int ndims, uint64_t *p_offsets, uint64_t *offsets, uint64_t *p_dims, uint64_t *ldims, uint64_t *gdims)
{
    if(ndims==1) {
        gdims[0]=ldims[0]+p_dims[0];
    }
    else if(ndims==2) {
        if(p_offsets[0]!=offsets[0] && p_offsets[1]==offsets[1]) 
               gdims[0]=ldims[0]+p_dims[0];
        else if(p_offsets[1]!=offsets[1] && p_offsets[0]==offsets[0]) 
            gdims[1]=ldims[1]+p_dims[1];
    }
    else if(ndims==3) {
        if(p_offsets[0]!=offsets[0] && p_offsets[1]==offsets[1] && p_offsets[2]==offsets[2]) 
               gdims[0]=ldims[0]+p_dims[0];
        else if(p_offsets[1]!=offsets[1] && p_offsets[0]==offsets[0] && p_offsets[2]==offsets[2])
            gdims[1]=ldims[1]+p_dims[1];
        else if(p_offsets[2]!=offsets[2] && p_offsets[0]==offsets[0] && p_offsets[1]==offsets[1])
            gdims[2]=ldims[2]+p_dims[2];
    }
}

static uint64_t do_spatial_aggr(int level, int *procs, int ndims, uint64_t *ldims, uint64_t *offsets, char *new_ldims, int rank,  void *data, uint64_t varsize, void *output, int type_size, MPI_Comm comm) 
{
    //struct adios_var_struct * v = g->vars;
    int i, j, k, client_cnt, lev;
    uint64_t aggrsize, buff_offset, tmpsize, alloc_size;
    uint64_t *tmp_dims, *tmp_offsets, *gdims, *ldims_list, *size_list;
    char   *tmpbuf, *recvbuf, *sendbuf;
    MPI_Status status;

    //store the local dimensions after aggregation
    gdims=malloc(ndims*sizeof(uint64_t));

    for(lev=0;lev<aggr_level;lev++) {
        /* aggregator */
        if(my_aggregator[ndims-1][lev] == rank) {
            if(lev==0) {
                alloc_size=varsize*(aggr_cnt[ndims-1][lev]+1);
#ifdef DEBUG 
                printf("%d: i am allocating %llu varsize=%llu\n", rank, alloc_size, varsize); 
#endif
                //XXX: hard coded memory allocation since realloc has seg fault
                tmpbuf=(char *)malloc(alloc_size*sizeof(char));
                //aggregator copies its own data first
                memcpy(tmpbuf, data, varsize);
                buff_offset=varsize;

                //allocate receive buffer for offsets and local_dims
                recvbuf=(char *)malloc(2*ndims*sizeof(uint64_t));

                tmp_dims=(uint64_t *)malloc(ndims*sizeof(uint64_t));
                tmp_offsets=(uint64_t *)malloc(ndims*sizeof(uint64_t));
            }
            
            else{
#ifdef DEBUG 
                printf("%d: reallocating %llu varsize=%llu\n", rank, (aggr_cnt[ndims-1][lev]+1)*varsize, varsize); 
#endif
                tmpbuf=(char *)realloc(tmpbuf, (aggr_cnt[ndims-1][lev]+1)*varsize);
                buff_offset=varsize;
            }

#ifdef DEBUG 
    printf("%d: ldims = %d %d %d  varsize=%llu nclient=%d \n", rank, ldims[0], ldims[1], ldims[2], varsize, aggr_cnt[ndims-1][lev]); 
#endif

            //store the local dimensions of all the chunks
            if(lev==0) {
#ifdef DEBUG
                printf("%d: ldims_list size=%d  size+list = %d\n", rank, (ndims*(aggr_cnt[ndims-1][lev]+1)), aggr_cnt[ndims-1][lev]+1);
#endif
                ldims_list=(uint64_t *)malloc(ndims*(aggr_cnt[ndims-1][lev]+1)*sizeof(uint64_t));
                size_list=(uint64_t *)malloc((aggr_cnt[ndims-1][lev]+1)*sizeof(uint64_t));
            }

            //prepare the local dimension list
            memcpy(ldims_list, ldims, ndims*sizeof(uint64_t));
            //initialize the gdims
            memcpy(gdims, ldims, ndims*sizeof(uint64_t));
            size_list[0]=varsize;

            k=1;
            //gather data from clients
            for(i=0;i<aggr_cnt[ndims-1][lev];i++) {
//                printf("%d: inside of loop  client=%d\n", rank, aggr_cnt[ndims-1][lev]);
#ifdef DEBUG 
                if(ndims==2) 
                    printf("rank %d: level=%d  waiting on................ %d \n", rank, lev, aggr2d_clients[lev][i].rank);
                else if (ndims==3)
                    printf("rank %d: level=%d  waiting on................ %d \n", rank, lev, aggr3d_clients[lev][i].rank);
#endif
                //memset(tmp_dims, 0x00, ndims*sizeof(uint64_t));
                //memset(tmp_offsets, 0x00, ndims*sizeof(uint64_t));
                //receive the ldims of the client process
                
                if(ndims==1) { 
                    MPI_Recv (recvbuf, 2*ndims*sizeof(uint64_t), MPI_BYTE, aggr1d_clients[lev][i].rank, 
                        aggr1d_clients[lev][i].rank, comm, &status);
                }
                else if(ndims==2) { 
                    MPI_Recv (recvbuf, 2*ndims*sizeof(uint64_t), MPI_BYTE, aggr2d_clients[lev][i].rank, 
                        aggr2d_clients[lev][i].rank, comm, &status);
                }
                else if(ndims==3) { 
                    MPI_Recv (recvbuf, 2*ndims*sizeof(uint64_t), MPI_BYTE, aggr3d_clients[lev][i].rank, 
                        aggr3d_clients[lev][i].rank, comm, &status);
                }

                //keep the ldims to the list
                memcpy(ldims_list+(i+1)*ndims, recvbuf, ndims*sizeof(uint64_t));
                memcpy(tmp_dims, recvbuf, ndims*sizeof(uint64_t));

#ifdef DEBUG 
                if(ndims==1) 
                    printf("%d: received from -- %d dims=%llu\n", rank, aggr1d_clients[lev][i].rank, tmp_dims[0]);
                else if(ndims==2) 
                    printf("%d: lev=%d received from -- %d dims=%llu   %llu\n", rank, lev, aggr2d_clients[lev][i].rank, tmp_dims[0], tmp_dims[1]); 
                else if(ndims==3) 
                    printf("%d: lev=%d received from -- %d dims=%llu   %llu\n", rank, lev, aggr3d_clients[lev][i].rank, tmp_dims[0], tmp_dims[1]); 
#endif

                //calculate the chunk size
                tmpsize=type_size;
                for(j=0;j<ndims;j++)
                    tmpsize*=tmp_dims[j];

                size_list[k]=tmpsize;
                k++;

                //calculate the aggregated chunk size
                //the chunk sizes may be different from the clients
                //XXX: FIXME harded for now
                //FIXME: incorrent algorithm, it only works for 3D domain
                //decomposition, probably we need a map
                
                //get the offsets
                memcpy(tmp_offsets, recvbuf+ndims*sizeof(uint64_t), ndims*sizeof(uint64_t));

                //calculate the aggregated chunk dimension
                //XXX: maybe better way to code this?
                cal_gdims(ndims, tmp_offsets, offsets, tmp_dims, ldims, gdims);

#ifdef DEBUG 
            if(ndims==2) 
                printf("rank %d: lev=%d waiting on................ %d size=%llu\n", rank, lev, aggr2d_clients[lev][i].rank, tmpsize); 
            if(ndims==3) 
                printf("rank %d: lev=%d waiting on................ %d size=%llu\n", rank, lev, aggr3d_clients[lev][i].rank, tmpsize); 
#endif

                //receive the data from the client process
                if(ndims==1) { 
                    MPI_Recv (tmpbuf+buff_offset, tmpsize, MPI_BYTE, 
                            aggr1d_clients[lev][i].rank, aggr1d_clients[lev][i].rank, comm, &status);
#ifdef DEBUG 
            printf("%d: received from -- %d, lev=%d size=%llu now offset=%llu\n", rank, aggr1d_clients[lev][i].rank, lev, tmpsize, buff_offset);
#endif
                }
                else if(ndims==2) { 
                    MPI_Recv (tmpbuf+buff_offset, tmpsize, MPI_BYTE, 
                            aggr2d_clients[lev][i].rank, aggr2d_clients[lev][i].rank, comm, &status);
#ifdef DEBUG 
            printf("%d: received from -- %d, lev=%d size=%llu now offset=%llu\n", rank, aggr2d_clients[lev][i].rank, lev, tmpsize, buff_offset);
#endif
                }
                if(ndims==3) { 
                    MPI_Recv (tmpbuf+buff_offset, tmpsize, MPI_BYTE, 
                            aggr3d_clients[lev][i].rank, aggr3d_clients[lev][i].rank, comm, &status);
#ifdef DEBUG 
            printf("%d: received from -- %d, lev=%d size=%llu now offset=%llu\n", rank, aggr3d_clients[lev][i].rank, lev, tmpsize, buff_offset);
#endif
                }

                //move the pointer
                buff_offset=buff_offset+tmpsize;
            }//end of for(i=0)

            for(i=0;i<ndims;i++) {
                if(i==0)
                    sprintf(new_ldims, "%llu", gdims[i]);
                else
                    sprintf(new_ldims, "%s,%llu", new_ldims, gdims[i]);
            }
                
            //aggregate the chunks
            //1D variable doesn't need to do this since the data is already
            //aligned in the buffer
            if(ndims>1) 
                aggr_chunks(&tmpbuf, procs, ndims, ldims_list, gdims, size_list, buff_offset, aggr_cnt[ndims-1][lev]+1, rank, lev, type_size);

            //update the local dimensions at the aggregators 
            memcpy(ldims, gdims, ndims*sizeof(uint64_t));
            //the buffer offset marks the current variable size
            varsize=buff_offset;
        }
        else{

            //the previous level clients don't need to do anything 
            if(lev>0 && my_aggregator[ndims-1][lev-1]!=rank)
                continue;
            else 
                sendbuf=(char *)malloc(2*ndims*sizeof(uint64_t));
#ifdef DEBUG 
            printf("rank %d: send...to %d, dims=%llu %llu %llu\n", rank, my_aggregator[ndims-1][lev], ldims[0], ldims[1], ldims[2]);
#endif
                
            //put in local dimensions
            memcpy(sendbuf, ldims, ndims*sizeof(uint64_t));
            //put in offsets
            memcpy(sendbuf+ndims*sizeof(uint64_t), offsets, ndims*sizeof(uint64_t));
            
            //clients send out the local dimension of data chunks
            MPI_Send(sendbuf, 2*ndims*sizeof(uint64_t), MPI_BYTE, my_aggregator[ndims-1][lev], rank, comm);
#ifdef DEBUG 
            printf("rank %d: sent...to %d, dims=%llu %llu %llu\n", rank, my_aggregator[ndims-1][lev], ldims[0], ldims[1], ldims[2]);
            if(lev==1)
                printf("%d:  sending to aggregator = %d size=%llu\n", rank, my_aggregator[ndims-1][lev], varsize);
#endif
            /*clients send out the data*/
            if(lev==0)
                MPI_Send(data, varsize, MPI_BYTE, my_aggregator[ndims-1][lev], rank, comm);
            else
                MPI_Send(tmpbuf, varsize, MPI_BYTE, my_aggregator[ndims-1][lev], rank, comm);

        } //end of if(rank==aggregator)
    }//end of for(lev=0)

    //the final aggregators need to return data
    if(my_aggregator[ndims-1][aggr_level-1]==rank) {
        memcpy(output, tmpbuf, buff_offset);
        //release the resources
        free(tmp_dims);
        free(gdims);
        free(ldims_list);
        free(size_list);
    }

    return buff_offset;
}


static void aggr_chunks(void **output, int *procs, int ndims, uint64_t *ldims_list,
    uint64_t *gdims, uint64_t *size_list, uint64_t totalsize, int nchunks, int rank, int level, int type_size)
{
    uint64_t count[3];
    uint64_t size;
    int i,j,k, m, cnt;
    uint64_t var_offset, dset_offset, buff_offset, size_in_dset[2];
    uint64_t datasize, dst_stride, src_stride;
    //cycles_t c1, c2;
    double tmem=0;
    char *input;
    int chunk_cnt;
    uint64_t prev_x, prev_y, prev_z;
    uint64_t m_offx, m_offy, m_offz;
    uint64_t offx, offy, offz;
    int ni, nj, nk;

#ifdef DEBUG
    printf("%d: chunk size = %d   %d  %d \n", rank, gdims[0],gdims[1], gdims[2]);
    printf("%d: totalsize = %llu nchunks=%d\n", rank, totalsize, nchunks);
#endif

    chunk_cnt=(int)pow(2, ndims);
    input=(char *)malloc(totalsize);
    memcpy(input, *output, totalsize); 

    dst_stride=1;
    src_stride=1;

    //for 3D variable, the number of chunks can only be 2, 4 or 8
    nk=1;
    nj=1;
    ni=1;

    m_offx=rank%procs[0];
    m_offy=rank/procs[0]%procs[1];
    m_offz=rank/(procs[0]*procs[1]);

    //determine the number of chunks on each dimension
    for(i=0;i<aggr_cnt[ndims-1][level];i++) {
        offx=offy=offz=0;
        if(ndims==1) { 
            offx=aggr1d_clients[level][i].rank%procs[0];
        }
        if(ndims==2) {
            offx=aggr2d_clients[level][i].rank%procs[0];
            offy=aggr2d_clients[level][i].rank/procs[0]%procs[1];
        }
        else if(ndims==3) {
            offx=aggr3d_clients[level][i].rank%procs[0];
            offy=aggr3d_clients[level][i].rank/procs[0]%procs[1];
            offz=aggr3d_clients[level][i].rank/(procs[0]*procs[1]);
        }

        if(offx!=m_offx && offy==m_offy) {
           if(ndims<3 || (ndims==3 && offz==m_offz))
               nk++;
        }
        else if(offy!=m_offy && offx==m_offx) {
           if(ndims==2 ||(ndims==3 && offz==m_offz)) 
               nj++; 
        }
        
        //only applies for 3D
        if(ndims==3 && offz!=m_offz && offx==m_offx && offy==m_offy)
            ni++;
    }
    
#ifdef DEBUG 
    printf("ni=%d  nj=%d  nk=%d\n", ni, nj, nk);
#endif
    cnt=0;
    prev_x=prev_y=prev_z=0;
    for(i=0;i<ni;i++) {
        for(j=0;j<nj;j++) {
            for(k=0;k<nk;k++) {
#ifdef DEBUG 
                if(ndims==2) 
                    printf("%d: chunk size = %llu   %llu  datasize=%llu\n", rank, ldims_list[ndims*cnt+0], ldims_list[ndims*cnt+1], size_list[cnt]); 
                if(ndims==3) 
                    printf("chunk size = %llu   %llu  %llu datasize=%llu\n", ldims_list[ndims*cnt+0], ldims_list[ndims*cnt+1], ldims_list[ndims*cnt+2], size_list[cnt]); 
#endif
                if(ndims==1) {
                    size_in_dset[0]=1;
                    size_in_dset[1]=ldims_list[ndims*cnt+0];
                }
                if(ndims==2) { 
                    size_in_dset[0]=ldims_list[ndims*cnt+1];
                    size_in_dset[1]=ldims_list[ndims*cnt+0];
                }
                else if(ndims==3) { 
                    size_in_dset[0]=ldims_list[ndims*cnt+2];
                    size_in_dset[1]=ldims_list[ndims*cnt+1];
                }

                datasize=ldims_list[ndims*cnt+0];
                src_stride=ldims_list[ndims*cnt+0];
                dst_stride=gdims[0];
                var_offset=0;

                if(cnt!=0) { 
                    //start offset in the aggregated chunk 
                    dset_offset=i*prev_z*gdims[0]*gdims[1]+j*prev_y*gdims[0]+k*prev_x;
                    //start offset in the input buffer 
                    buff_offset+=size_list[cnt-1];
                }
                else {
                    dset_offset=0;
                    buff_offset=0;
                }

#ifdef DEBUG 
        printf("rank=%d: ******input a new chunk, dst_offset=%llu  input_offset=%llu\n", rank, dset_offset, buff_offset);
#endif
#ifdef DEBUG 
        printf("size in dest= %llu %llu\n", size_in_dset[0], size_in_dset[1]);
        printf("gdims = %llu %llu %llu\n", gdims[0], gdims[1], gdims[2]);
        printf("dst_stride=%llu     src_stride=%llu     datasize=%llu\n", dst_stride, src_stride, datasize);
#endif
                copy_aggr_data(*output
                   ,input+buff_offset
                   ,0
                   ,ndims-1
                   ,size_in_dset
                   ,ldims_list+ndims*cnt
                   ,gdims
                   ,dst_stride
                   ,src_stride
                   ,dset_offset
                   ,var_offset
                   ,datasize
                   ,type_size
                   ,rank
                   );
                cnt++;
                prev_x=ldims_list[ndims*(cnt-1)+0];
            }
            prev_y=ldims_list[ndims*(cnt-1)+1];
        }//end of j
        if(ndims==3) { 
            prev_z=ldims_list[ndims*(cnt-1)+2];
        }
    } //end of i

    free(input);
}
