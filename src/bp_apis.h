int bp_fread ( int64_t * fh,
	        char * fname,
		MPI_Comm comm
	      );
void bp_fclose ( int64_t fh);
void bp_inq_file ( int64_t fh_p, int *ngroup, 
		  int *nvar, int *nattr, int *nt, char **gnamelist); 
void bp_gopen ( int64_t fh,
		int64_t * gh, 
		char * grpname);
void bp_gclose ( int64_t fh);
void bp_inq_group (int64_t gh_p, int *nvar, char ** vnamelist);
void bp_get_var (int64_t fh,
		 char * varname,
		 void * var, 
		 int  *,
		 int  *, 
		 int);
void bp_inq_var (int64_t gh_p, char * varname,
		 int * type,
		 int * ndim,
		 int * dims);
