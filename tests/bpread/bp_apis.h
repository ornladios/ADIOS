int bp_fread ( int64_t * fh,
	        char * fname,
		MPI_Comm comm
	      );
void bp_fclose ( int64_t fh);
void bp_gopen ( int64_t fh,
		int64_t * gh, 
		char * grpname);
void bp_gclose ( int64_t fh);
void bp_get_var (int64_t fh,
		 char * varname,
		 void * var, 
		 int  *,
		 int  *, 
		 int);
