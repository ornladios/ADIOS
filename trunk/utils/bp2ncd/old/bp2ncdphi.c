#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "binpack-general.h"
#include "br-utils.h"
#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}
int makencd(char *ifname, char *ofname,char *destname);
int ncd_sizedset(int ncid,char *path,char *name,size_t *totalsize,\
                        int ranks, struct adios_bp_dimension_struct *dims);
int ncd_wdset(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type, int rank,              \
              struct adios_bp_dimension_struct *dims);
int ncd_wscalar(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type);

size_t start[2],count[2],totalsize[2];
size_t start_scalar[1],count_scalar[1];
int main(int argc, char** argv)
{
  char *ifname;
  char *ofname;
  char *dsetname;
  ofname=NULL;
  dsetname=NULL; 
  
  if(argc<2)
  {
     printf ("\nUSAGE: bp2cdf [-d dataset_name] XXXX.bp [XXXX.ncd]\n\n");
     return -1;
  }
  else if (argc<3)
  {
     ifname=argv[1];
  }
  else if (argc<4)
  {
     if(argv [1][0]& argv[1][0] == '-')
     {
        if(!strcmp(argv[1],"-d"))
        {
           dsetname=argv[2];
           ifname=argv[3];       
        }
     }
     else
     {
        ifname=argv[1];       
        ofname=argv[2];       
     }
  }  
  else if (argc<5)
  {
     if(argv [1][0]& argv[1][0] == '-')
        if(!strcmp(argv[1],"-d"))
        {
           dsetname=argv[2];
           ifname=argv[3];       
           ofname=argv[4];       
        }
  }
  makencd(ifname,ofname,dsetname);
  return 0;  
}

int makencd(char *ifname, char *ofname,char *destname)
{
  long long handle=0;
  int ncid,retval;
  char *val;
  char *tmpstr=NULL;
  unsigned long long DATALEN=100*1024*1024;
  struct adios_bp_element_struct *element=NULL;
  long long element_size=0; 
  if(ofname==NULL)
  {
     int size=strlen(ifname);
     tmpstr=malloc(strlen(ifname));
     strcpy(tmpstr,ifname);
     tmpstr[size-2]='n';
     tmpstr[size-1]='c';
  }
  else
     tmpstr=strdup(ofname);
  handle = br_fopen(ifname);
  if(!handle)
  {
     fprintf(stderr, "Failed to open %s file!\n",ifname);
     return -1;
  } 
  val=(char *)malloc(DATALEN);
  if(!val)
  {
     fprintf(stderr,"malloc failed!\n");
     return -1;
  } 

  char *dataname;
  nc_create(tmpstr,NC_CLOBBER | NC_64BIT_OFFSET,&ncid);
  int i,j;
  handle = br_fopen(ifname);
  while (element_size = br_get_next_element_general (handle ,val
                                                      ,DATALEN
                                                      ,&element
                                                      )
          )
  {
        switch (element->tag)
        {
            case DSTATRS_TAG:
//                hw_attr_str_ds (root_id, element->path, element->name, val);
                break;
              case DSTATRN_TAG:
//                hw_attr_num_ds (root_id, element->path, element->name, val \
//                               ,element->type				     \
//                               );         
                break;

            case GRPATRS_TAG:
//                hw_attr_str_gp (root_id, element->path, element->name, val);
                break;

            case GRPATRN_TAG:
//                hw_attr_num_gp (root_id, element->path, element->name, val \
                               ,element->type 				     \
                               );
                break;

            case SCR_TAG:
//                hw_scalar (root_id, element->path, element->name, val        \
                          ,element->type, 0                                  \
                          );
                ncd_wscalar(ncid,element->path,element->name,val,element->type);
                break;
            case DST_TAG:
                   for(i=0;i<strlen(element->path);i++)
                   {
                      j=*(element->path+i)-'0';
                      if(j<10 && j>-1)
                        break;
                   }
                   dataname=element->path+i;
                   start[0]=atoi(dataname);
                   start_scalar[0]=atoi(dataname);
                   //printf("%s,%d\n",element->name,start[0]);
                   ncd_wdset(ncid,element->path,element->name,val,element->type,
                        element->ranks, element->dims);
                break;
            default:
                break;
        }
        br_free_element (element);
     } // end of while
  br_fclose (handle);
  nc_close(ncid); 
#if DEBUG
  printf("out file: %s\n",tmpstr);
#endif
  if(!tmpstr)free(tmpstr); 
  return 0; 
}
int ncd_wdset(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type, int rank,              \
              struct adios_bp_dimension_struct *dims)
{
  int i,valid,retval;
  char fullname[100];
  char dimname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int dimids[2];
  dimids[1]=0;

  strcpy(fullname,name);
  if(start[0]==0)
  {
     retval=nc_redef(ncid);
     sprintf(dimname,"%s_%d",name,0);
     retval=nc_inq_unlimdim(ncid,&dimids[1]);
     if(dimids[1]==-1)retval=nc_def_dim(ncid,"rec_idx",NC_UNLIMITED,&dimids[0]);
     retval=nc_def_dim(ncid,dimname,dims[0].local_bound,&dimids[1]);
     if(type==bp_float)
        retval=nc_def_var(ncid,fullname,NC_FLOAT,2,dimids,&valid);
     if(type==bp_int)
        retval=nc_def_var(ncid,fullname,NC_INT,2,dimids,&valid);
     ERR(retval);
     printf("\t RANK=%d, DIMS:%s[0]=%d dimids[1]=%d\n",rank,dimname,dims[0].local_bound,dimids[1]);
     retval=nc_enddef(ncid);
  }
  else
     retval=nc_inq_varid (ncid, fullname, &valid);
  start[1]=0;//dims[0].local_bound;
 // count[0]=dims[1].local_bound;
  count[0]=1;
  count[1]=dims[0].local_bound;

  if(type==bp_float)
     retval=nc_put_vara_float(ncid,valid,start,count,val);
  if(type==bp_int)
     retval=nc_put_vara_int(ncid,valid,start,count,val);
  //retval=nc_put_vara_float(ncid,valid,start,count,val);
//#if DEBUG
  printf("create dataset:%s\n",fullname);
  printf("start:%dx%d. count:%dx%d\n",start[0],start[1],count[0],count[1]);
  printf("-------------------\n");
// iret = nf90_put_att(diag_flow_ncfileid, diag_flow_ncvarid(j,sp_type), 'units', var_units(j))
//#endif
  return;
}

int ncd_wscalar(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type)
{
  int i,valid,retval;
  char fullname[100];
  char dimname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int dimids[1];
  dimids[0]=0;

  strcpy(fullname,name);
  if(start[0]==0)
  {
     retval=nc_redef(ncid);
     sprintf(dimname,"%s_%d",name,0);
     retval=nc_inq_unlimdim(ncid,&dimids[0]);

     if(dimids[0]==-1)retval=nc_def_dim(ncid,"rec_idx",NC_UNLIMITED,&dimids[0]);
    
     if(type==bp_float)
        retval=nc_def_var(ncid,fullname,NC_FLOAT,1,dimids,&valid);
     if(type==bp_int)
        retval=nc_def_var(ncid,fullname,NC_INT,1,dimids,&valid);

     retval=nc_enddef(ncid);
  }

  else
     retval=nc_inq_varid (ncid, fullname, &valid);
  count_scalar[0]=1;
  if(type==bp_float)
     retval=nc_put_vara_float(ncid,valid,start_scalar,count_scalar,val);
  if(type==bp_int)
     retval=nc_put_vara_int(ncid,valid,start_scalar,count_scalar,val);
/*
//#if DEBUG
  printf("create scalar:%s\n",fullname);
  printf("start:%d. count:%d\n",start_scalar[0],count_scalar[0]);
  printf("-------------------\n");
// iret = nf90_put_att(diag_flow_ncfileid, diag_flow_ncvarid(j,sp_type), 'units', var_units(j))
//#endif
*/  return;
}


