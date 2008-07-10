#include "hw-utils.h"
#include "hdf5.h"
void dump2ascii(char *argv[],int argc,int i);
double *hr_dataset(hid_t root_id, char* name,hsize_t *rank,hsize_t *dims,hid_t *type_id);
int main (int argc, char ** argv)
{
    if (argc < 4)
    {
        printf ("\nUSAGE: bp2ascii XXXX.bp  [XXXX.txt] -C[-c] data_name1  [data_name2...]\n\n");
        return -1;
    }
    int i;
    hw_makeh5(argv[1]);
    for(i=0;i<argc;i++)
      if((strcmp(argv[i],"-C")==0) || (strcmp(argv[i],"-c")==0))
        break;
   
    dump2ascii(argv,argc,i+1); 
    return 0;
}

void dump2ascii(char* argv[],int argc,int idx)
{
     char filename[100];
     int len,i,j,id;
     void *data;
     hid_t h5file_id,root_id,type_id;
     hsize_t rank,dims[3];
     strcpy(filename,argv[1]);
     len =strlen(filename);
     filename[len-2]='h';
     filename[len-1]='5';
     h5file_id=H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
     if(h5file_id<0)
       return;
     root_id=H5Gopen(h5file_id,"/");
     if(idx==3)
     { 
     strcpy(filename,argv[1]);
     len=strlen(argv[1]);
     filename[len-2]='t';
     filename[len-1]='x';
     filename[len]='t';
     filename[len+1]='\0';
     }
     if(idx==4)
     strcpy(filename,argv[2]);

     printf("----------------------------\n");
     printf("generating file: %s\n",filename);
     FILE* fid=fopen(filename,"w");
     for(id=idx;id<argc;id++)
     {
     printf("----------------------------\n");
     printf("write dataset: %s\n",argv[id]);
     data=hr_dataset(root_id,argv[id],&rank,dims,&type_id);
     if(rank==1)
     {
        if(H5Tequal(type_id,H5Tcopy(H5T_IEEE_F64LE)))
        {
           for(i=0;i<dims[0];i++)
           {
              fprintf(fid,"%e  ",*((double *)(data+8*i)));
              printf("\t%e  ",*((double *)(data+8*i)));
           }
           printf("\n");
        }
        else if(H5Tequal(type_id,H5Tcopy(H5T_IEEE_F32LE)))
        {
           {
              for(i=0;i<dims[0];i++)
              {
                  fprintf(fid,"%e  ",*((float *)(data+4*i)));
                  printf("\t%e  ",*((float *)(data+4*i)));
              }
              printf("\n");
           }
        }
        else if(H5Tequal(type_id,H5Tcopy(H5T_STD_I32LE)))
        {
           for(i=0;i<dims[0];i++)
           {
              fprintf(fid,"%d  ",*((int*)(data+4*i)));
              printf("\t%d  ",*((int *)(data+4*i)));
           }
           fprintf(fid,"%d  ",*((int*)(data+4*i)));
           printf("\n");
        }
            
     }
     if(rank==2)
     {
        if(H5Tequal(type_id,H5Tcopy(H5T_IEEE_F32LE)))
              for(i=0;i<dims[0];i++)
              {
                  for(j=0;j<dims[1];j++)
                      fprintf(fid,"%e  ",*((float*)(data+i*dims[1]+4*j)));
                  fprintf(fid,"\n"); 
              }
        else if(H5Tequal(type_id,H5Tcopy(H5T_IEEE_F64LE)))
              for(i=0;i<dims[0];i++)
              {
                  for(j=0;j<dims[1];j++)
                      fprintf(fid,"%e  ",*((double*)(data+i*dims[1]+8*j)));
                  fprintf(fid,"\n"); 
              }
        else if(H5Tequal(type_id,H5Tcopy(H5T_STD_I32LE)))
              for(i=0;i<dims[0];i++)
              {
                  for(j=0;j<dims[1];j++)
                      fprintf(fid,"%e  ",*((int*)(data+i*dims[1]+4*j)));
                  fprintf(fid,"\n"); 
              }
        else if(H5Tequal(type_id,H5Tcopy(H5T_STD_I64LE)))
              for(i=0;i<dims[0];i++)
              {
                  for(j=0;j<dims[1];j++)
                      fprintf(fid,"%e  ",*((long *)(data+i*dims[1]+8*j)));
                  fprintf(fid,"\n"); 
              }
     }
     }
     fclose(fid);
}
void *hr_dataset(hid_t root_id, char* name,hsize_t *rank,hsize_t *dims,hid_t *type_id)
{
   int i,j,size;
   void *temp;
   hid_t dsethist_id,filespace_id;
   //type_id = H5Tcopy(H5T_IEEE_F64LE);
   dsethist_id=H5Dopen(root_id,name);
   *type_id=H5Dget_type(dsethist_id);
   size=H5Tget_size(*type_id);  
   filespace_id=H5Dget_space(dsethist_id);
   *rank=H5Sget_simple_extent_ndims(filespace_id);
   H5Sget_simple_extent_dims(filespace_id,dims,NULL);
   if (*rank==1)
   {
      printf("\tDimension:%d \n",dims[0]);
      temp= malloc(size*dims[0]);
      H5Dread(dsethist_id,*type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
   }
   else if (*rank==2)
   {
      printf("\tDimension:%d x %d \n",dims[0],dims[1]);
      temp=malloc(size*dims[0]*dims[1]);
      H5Dread(dsethist_id,*type_id,H5S_ALL,H5S_ALL,H5P_DEFAULT,temp);
   }
  return  temp;
}
