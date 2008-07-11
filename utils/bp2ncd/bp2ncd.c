#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#include "binpack-general.h"
#include "br-utils.h"
#include "adios_types.h"

#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}

int makencd(char *ifname, char *ofname,char *destname);
int ncd_wscalar(int ncid,char *path,char *name, void *val,  \
              enum ADIOS_TYPES type);
int ncd_wdset(int ncid,char *path,char *name, void *val,  \
              enum ADIOS_TYPES type, int rank,              \
              struct adios_bp_dimension_struct *dims);
int ncd_attr_str_ds (int ncid, char *path, char *name, enum ADIOS_TYPES type,void *val);
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
                printf("attribute\n");
                ncd_attr_str_ds (ncid, element->path, element->name, element->type,val); 
                break;
              case DSTATRN_TAG:
                printf("attribute num\n");
                ncd_attr_str_ds (ncid, element->path, element->name, element->type,val); 
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
                 ncd_wscalar(ncid,element->path,element->name,val,element->type);
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
              enum ADIOS_TYPES type, int rank,              \
              struct adios_bp_dimension_struct *dims)
{
  int i,valid,retval,j;
  char tmpname[100];
  char dimname[100];
  char fullname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int *dimids=NULL;
  nc_redef(ncid);
  strcpy(tmpname,path);
  for(j=0;j<strlen(tmpname);j++)
      if(tmpname[j]=='[' || tmpname[j]==']')
           tmpname[j]='_';
  result = strtok(tmpname,"/");
  while(result!=NULL)
  {
     sprintf(fullname,"%s_%s",fullname,result);
     result=strtok(NULL,"/");
  }
  strcpy(tmpname,name);
  for(j=0;j<strlen(tmpname);j++)
      if(tmpname[j]=='[' || tmpname[j]==']')
           tmpname[j]='_';
  sprintf(dimname,"%s_%s",fullname,tmpname);
  strcpy(fullname,dimname+1);
  dimids=(int *)malloc(sizeof(int)*rank);
  float *data=NULL;
  for(i=0;i<rank;i++)
  {
     sprintf(dimname,"%s_%d",fullname,i);
     retval=nc_def_dim(ncid,dimname,dims[i].local_bound,&dimids[i]);
//     printf("\tfullname=%s, DIMS:%s=%d, dimids[%d]=%d\n",fullname,dimname,dims[i].local_bound,i,dimids[i]);
     ERR(retval);
  }
  switch(type)
  {
     case adios_real:
          retval=nc_def_var(ncid,fullname,NC_FLOAT,rank,dimids,&valid);
          retval=nc_enddef(ncid);
          retval=nc_put_var_float(ncid,valid,val);
          break;
     case adios_double:
          retval=nc_def_var(ncid,fullname,NC_DOUBLE,rank,dimids,&valid);
          retval=nc_enddef(ncid);
          retval=nc_put_var_double(ncid,valid,val);
          break;
     case adios_long:
          retval=nc_def_var(ncid,fullname,NC_LONG,rank,dimids,&valid);
          retval=nc_enddef(ncid);
          retval=nc_put_var_float(ncid,valid,val); 
          break;
     case adios_integer:
          retval=nc_def_var(ncid,fullname,NC_INT,rank,dimids,&valid);
          retval=nc_enddef(ncid);
          retval=nc_put_var_int(ncid,valid,val);
          break;
     default:
          retval=nc_enddef(ncid);
          break;
  }
#if DEBUG
  ERR(retval);
#endif
#if DEBUG
  printf("create dataset:%s\n",fullname);
#endif
#if DEBUG
  printf("-------------------\n");
#endif
  if(dimids)free(dimids);
  return;
}

int ncd_wscalar(int ncid,char *path,char *name, void *val,  \
              enum ADIOS_TYPES type)
{
  int i,valid,retval;
  char fullname[100];
  char dimname[100];
  fullname[0]='\0';
  char *result=NULL; 
  int dimids[1];
  dimids[0]=1;
  result = strtok(path,"/");
  while(result!=NULL)
  {
     sprintf(fullname,"%s_%s",fullname,result);
     result=strtok(NULL,"/");
  }
  sprintf(dimname,"%s_%s",fullname,name);
  strcpy(fullname,dimname+1); 
  retval=nc_inq_varid (ncid, fullname, &valid);
//  printf("fullname:%s,valid=%d\n",fullname,retval);
  if(retval<0)
  {
     retval=nc_redef(ncid);
     sprintf(dimname,"%s_%d",name,0);
     nc_def_dim(ncid,dimname,1,&dimids[0]);
     switch (type)
     {
        case adios_real:
             retval=nc_def_var(ncid,fullname,NC_FLOAT,1,dimids,&valid);
             break;
        case adios_integer:
             retval=nc_def_var(ncid,fullname,NC_INT,1,dimids,&valid);
             ERR(retval);
             break;
        case adios_long:
             retval=nc_def_var(ncid,fullname,NC_LONG,1,dimids,&valid);
             break;
        case adios_double:
             retval=nc_def_var(ncid,fullname,NC_DOUBLE,1,dimids,&valid);
             break;
        default:
             break;
     }
     retval=nc_enddef(ncid);
  //printf("fullname:%s,valid=%d\n",fullname,valid);
//  else
  count_scalar[0]=1;
  switch(type)
  {
    case adios_real:
       retval=nc_put_vara_float(ncid,valid,start_scalar,count_scalar,val);
       //ERR(retval);
       break;
    case adios_double:
       retval=nc_put_vara_double(ncid,valid,start_scalar,count_scalar,val);
       ERR(retval);
       break;
    case adios_integer:
       retval=nc_put_vara_long(ncid,valid,start_scalar,count_scalar,val);
       ERR(retval);
       break;
    case adios_long:
       retval=nc_put_vara_int(ncid,valid,start_scalar,count_scalar,val);
       ERR(retval);
       break;
    default:
       break;
  }//end of switch(type)
  }
  return 0;
}
int ncd_attr_str_ds (int ncid, char *path, char *name, enum ADIOS_TYPES type,void *val)
{
    int j;
    char tmpname[100],tmppath[100];
    int valid,retval;
    char *text=(char *)val;
    retval=nc_redef(ncid);
    retval=nc_inq_varid(ncid,path,&valid);
    if(retval>0)
    {
       retval=nc_put_att_text(ncid,valid,name,strlen(text),text);
       ERR(retval);
    }
    else 
    {
       if(path[strlen(path)-1]=='/')
       {
         if(path[0]=='/')
            strcpy(tmppath,path+1);
         else
            strcpy(tmppath,path);
         for(j=0;j<strlen(tmppath);j++)
         {
            if(tmppath[j]=='[' || tmppath[j]==']' || tmppath[j]=='/') 
               tmppath[j]='_';
         } 
         //printf("%s\n",tmppath);
         strcpy(tmpname,name);
         for(j=0;j<strlen(tmpname);j++)
         if(tmpname[j]=='[' || tmpname[j]==']')
           tmpname[j]='_';
         strcat(tmppath,tmpname); 
         if(type==adios_string)
         { 
            retval=nc_put_att_text(ncid,NC_GLOBAL,tmppath,strlen(text),text);
            ERR(retval);
         }
         else if(type==adios_integer)
         { 
            retval=nc_put_att_int(ncid,NC_GLOBAL,tmppath,NC_INT,1,text);
            ERR(retval);
         }
  
       }
    }
//    printf("varname=%s, valid=%d type=%d len=%d\n",path,valid,type,strlen(text));
       
    retval=nc_enddef(ncid);
    return 0;
}
