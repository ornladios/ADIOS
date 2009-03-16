
/*
header file for the subsetting read routines
March 2009, ORNL
*/

// C interface
/*
	IN:  fname
	     comm
	OUT: fh
*/

int bp_fread ( int64_t * fh,
	        char * fname,
		MPI_Comm comm
	      );

/*
	IN: fh 
*/
int bp_fclose ( int64_t fh);

/*
	IN:  fh_p 
	OUT: ngroup
	     nvar
	     nattr
	     nt
	     gnamelist
 */ 
int bp_inq_file ( int64_t fh_p, int *ngroup, 
		  int *nvar, int *nattr, int *nt, char **gnamelist);
 
/*
	IN:  fh
	     grpname 
	OUT: gh_p 
*/
int bp_gopen ( int64_t * gh_p,
		int64_t fh, 
		char * grpname);

/*
	IN:  fh
*/
int bp_gclose ( int64_t gh);

/*
	IN:  gh_p
*/
int bp_inq_group (int64_t gh, int *nvar, char ** vnamelist);

/*
	IN:  gh
	     varname 
	     start 
	     readsize
	     timestep 
	OUT: var 
*/
int bp_get_var (int64_t gh,
		 char * varname,
		 void * var, 
		 int  * start,
		 int  * readsize, 
		 int    timestep);

/*
	IN:  gh
	     varname 
	OUT: type 
	     ndim 
	     is_timebased
	     dims  
*/
int bp_inq_var (int64_t gh, char * varname,
		 int * type,
		 int * ndim,
		 int * is_timebased,
		 int * dims);
const char * bp_type_to_string (int type);

// Fortran interface

void bp_fread_ ( int64_t * fh,
	        char * fname,
		MPI_Comm comm,
		int * err,
		int fname_len
	      );

void bp_fclose_ ( int64_t * fh, int * err);

void bp_inq_file_ ( int64_t * fh_p, int *ngroup, 
		  int *nvar, int *nattr, int *nt, char **gnamelist, int * err); 

void bp_gopen_ (int64_t * gh_p,
		int64_t fh, 
		char * grpname);

void bp_gclose_ ( int64_t * gh);

void bp_inq_group_ (int64_t * gh, int *nvar, char ** vnamelist);

void bp_get_var_ (int64_t * gh,
		 char * varname,
		 void * var, 
		 int  * start,
		 int  * readsize, 
		 int timestep);

void bp_inq_var_ (int64_t * gh_p, char * varname,
		 int * type,
		 int * ndim,
		 int * is_timebased,
		 int * dims);
