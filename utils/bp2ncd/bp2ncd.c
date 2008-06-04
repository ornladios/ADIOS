#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "binpack-general.h"
#include "br-utils.h"
#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}
int makencd(char *ifname, char *ofname,char *destname);
int ncd_wdset(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type, int rank,              \
              struct adios_bp_dimension_struct *dims);
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

  nc_create(tmpstr,NC_CLOBBER | NC_64BIT_OFFSET,&ncid);
  while (element_size = br_get_next_element_general (handle
                                                      ,val
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
                break;
            case DST_TAG:
                //if(strcmp(element->name,"rdtemi")==0 && strcmp(element->path,"/node00000/param"))
                ncd_wdset (ncid,element->path, element->name, val
                        ,element->type, element->ranks, element->dims
                        );
                break;
            default:
                break;
        }
        br_free_element (element);
  }
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
  char dimname[100];
  char fullname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int *dimids=NULL;
  result = strtok(path,"/");
  nc_redef(ncid);
  while(result!=NULL)
  {
     sprintf(fullname,"%s_%s",fullname,result);
     result=strtok(NULL,"/");
  }
  dimids=(int *)malloc(sizeof(int)*rank);
  float *data=NULL;
//  int totalsize=1;
  sprintf(fullname,"%s_%s",fullname,name);
  for(i=0;i<rank;i++)
  {
     sprintf(dimname,"%s_%d",fullname,i);
     retval=nc_def_dim(ncid,dimname,dims[i].local_bound,&dimids[i]);
#if DEBUG
     printf("\tncid=%d, DIMS:%s=%d, dimids[%d]=%d\n",ncid,dimname,dims[i].local_bound,i,dimids[i]);
     ERR(retval);
#endif
//     totalsize=totalsize*dims[i].local_bound;
  } 
  //data=(float*)malloc(sizeof(float)*totalsize); 
  //for(i=0;i<totalsize;i++)
  //   data[i]=1.0*i; 
  retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
#if DEBUG
  ERR(retval);
#endif
  retval=nc_enddef(ncid);
#if DEBUG
  printf("create dataset:%s\n",fullname);
#endif
  //retval=nc_put_var_float(ncid,valid,&data[0]);
  retval=nc_put_var_float(ncid,valid,val);
  //retval=nc_put_var_float(ncid,valid,&((*((float*)val))[0]));
#if DEBUG
  printf("-------------------\n");
#endif
  if(dimids)free(dimids);
  //if(data)free(data);
  return;
}
