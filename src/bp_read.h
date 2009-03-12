
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
void bp_fclose ( int64_t fh);

/*
	IN:  fh_p 
	OUT: ngroup
	     nvar
	     nattr
	     nt
	     gnamelist
 */ 
void bp_inq_file ( int64_t fh_p, int *ngroup, 
		  int *nvar, int *nattr, int *nt, char **gnamelist);
 
/*
	IN:  fh
	     grpname 
	OUT: gh_p 
*/
void bp_gopen ( int64_t * gh_p,
		int64_t fh, 
		char * grpname);

/*
	IN:  fh
*/
void bp_gclose ( int64_t gh);

/*
	IN:  gh_p
*/
void bp_inq_group (int64_t gh, int *nvar, char ** vnamelist);

/*
	IN:  gh
	     varname 
	     start 
	     readsize
	     timestep 
	OUT: var 
*/
void bp_get_var (int64_t gh,
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
void bp_inq_var (int64_t gh, char * varname,
		 int * type,
		 int * ndim,
		 int * is_timebased,
		 int * dims);

// Fortran interface

int bp_fread_ ( int64_t * fh,
	        char * fname,
		MPI_Comm comm
	      );

void bp_fclose_ ( int64_t fh);

void bp_inq_file_ ( int64_t fh_p, int *ngroup, 
		  int *nvar, int *nattr, int *nt, char **gnamelist); 

void bp_gopen_ (int64_t * gh_p,
		int64_t fh, 
		char * grpname);

void bp_gclose_ ( int64_t gh);

void bp_inq_group_ (int64_t gh, int *nvar, char ** vnamelist);

void bp_get_var_ (int64_t gh,
		 char * varname,
		 void * var, 
		 int  * start,
		 int  * readsize, 
		 int timestep);

void bp_inq_var_ (int64_t gh_p, char * varname,
		 int * type,
		 int * ndim,
		 int * is_timebased,
		 int * dims);
