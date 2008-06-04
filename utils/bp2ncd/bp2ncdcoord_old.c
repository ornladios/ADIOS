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
int ncd_wgdset(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type, int rank,              \
              struct adios_bp_dimension_struct *dims);

size_t start[2],count[2],totalsize[2];
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

  char dataname[2];
  dataname[1]='\0'; 
  nc_create(tmpstr,NC_CLOBBER | NC_64BIT_OFFSET,&ncid);
  int i,j;
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
                if(strcmp(element->name,"X")==0 || 
                   strcmp(element->name,"Y")==0 ||
                   strcmp(element->name,"Z")==0)
                {
                   if(strcmp(element->name,"X")==0) 
                   { 
                      ncd_sizedset(ncid,element->path,element->name,totalsize,
                        element->ranks, element->dims);
#if DEBUG
                      printf("\t DIMS[0]=%d\n",totalsize[0]);
                      printf("\t DIMS[1]=%d\n",totalsize[1]);
#endif
                   }
                }
                else
                {
                   //ncd_wdset(ncid,element->path,element->name,val,element->type,
                   //     element->ranks, element->dims);
                }
                break;
            default:
                break;
        }
        br_free_element (element);
  }
  br_fclose (handle);
  for(i=0;i<3;i++)
  {
      
     for(j=0;j<2;j++)
        start[j]=0;
     //dataname[0]='X'+i; 
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
                break;
            case DST_TAG:
                //if(strcmp(element->name,dataname)==0) 
                   ncd_wgdset(ncid,element->path,element->name,val,element->type,
                        element->ranks, element->dims);
                break;
            default:
                break;
        }
        br_free_element (element);
  } // end of while
  br_fclose (handle);
  } // end of for
  nc_close(ncid); 
#if DEBUG
  printf("out file: %s\n",tmpstr);
#endif
  if(!tmpstr)free(tmpstr); 
  return 0; 
}
int ncd_sizedset(int ncid,char *path,char *name,size_t *totalsize,
                        int ranks, struct adios_bp_dimension_struct *dims)
{
  totalsize[0]+= dims[1].local_bound;
  totalsize[1]= dims[0].local_bound;
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
  int *dimids=NULL;
  dimids=(int*)malloc(sizeof(int)*rank);
/*  result = strtok(path,"/");
  while(result!=NULL)
  {
     sprintf(fullname,"%s_%s",fullname,result);
     result=strtok(NULL,"/");
  }
  sprintf(fullname,"%s_%s",fullname,name);
*/
  strcpy(fullname,name); 
  count[0]=dims[0].local_bound;
  {
     retval=nc_redef(ncid);
     for(i=0;i<rank;i++)
     {
        sprintf(dimname,"%s_%d",fullname,i);
#if DEBUG
        printf("\t RANK=%d, DIMS:%s[%d]=%d\n",rank,dimname,i,count[i]);
#endif
        retval=nc_def_dim(ncid,dimname,count[i],&dimids[i]);
#if DEBUG
        ERR(retval);
#endif
     }
#if DEBUG
     if(strcmp(name,"mtheta")==0)
        printf("mtheta: element->type=%d\n",type); 
#endif
     if(type==bp_float)
        retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
     if(type==bp_int)
        retval=nc_def_var(ncid,fullname,NC_INT,rank,dimids,&valid);
     //retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
#if DEBUG
     ERR(retval);
#endif
     retval=nc_enddef(ncid);
#if DEBUG
     ERR(retval);
#endif
  }
  
  if(type==bp_float)
     retval=nc_put_vara_float(ncid,valid,start,count,val);
  if(type==bp_int)
     retval=nc_put_vara_int(ncid,valid,start,count,val);
  start[0]+=dims[1].local_bound;

#if DEBUG
  printf("create dataset:%s\n",fullname);
  printf("-------------------\n");
#endif
  free(dimids);
  return;
}

int ncd_wgdset(int ncid,char *path,char *name, void *val,  \
              enum vartype_t type, int rank,              \
              struct adios_bp_dimension_struct *dims)
{
  int i,valid,retval;
  char fullname[100];
  char dimname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int *dimids=NULL;
  dimids=(int*)malloc(sizeof(int)*rank);
/*  result = strtok(path,"/");
  while(result!=NULL)
  {
     sprintf(fullname,"%s_%s",fullname,result);
     result=strtok(NULL,"/");
  }
  sprintf(fullname,"%s_%s",fullname,name);
*/
  strcpy(fullname,name);
  if(start[0]==0)
  {
     retval=nc_redef(ncid);
     for(i=0;i<rank;i++)
     {
#if DEBUG
     printf("\t RANK=%d, DIMS:%s[%d]=%d\n",rank,dimname,i,totalsize[i]);
#endif
     sprintf(dimname,"%s_%d",name,i);
     retval=nc_def_dim(ncid,dimname,totalsize[i],&dimids[i]);
#if DEBUG
     ERR(retval);
#endif
     }
     retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
#if DEBUG
     ERR(retval);
#endif
     retval=nc_enddef(ncid);
#if DEBUG
     ERR(retval);
#endif
  }
  else
     retval=nc_inq_varid (ncid, fullname, &valid);
  
  start[1]=0;//dims[0].local_bound;
  count[0]=dims[1].local_bound;
  count[1]=dims[0].local_bound;

  retval=nc_put_vara_float(ncid,valid,start,count,val);
  start[0]+=dims[1].local_bound;
//#if DEBUG
  printf("create dataset:%s\n",fullname);
  printf("start:%dx%d. count:%dx%d\n",start[0],start[1],count[0],count[1]);
  printf("-------------------\n");
//#endif
  free(dimids);
  return;
}

