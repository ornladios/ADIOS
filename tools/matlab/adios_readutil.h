/* Type: linked list of strings for dynamic building of a string list.*/
typedef struct adios_readutil_namelist adios_readutil_namelist;
struct adios_readutil_namelist {
    char * name;
    adios_readutil_namelist *next;
};

/********************************/
/*   Functions for adiosopenc   */
/********************************/
/* Open an adios file. 
   Input:  file name
   Output: fh file   handler
           ngroups   number of adios groups
           groupname linked list of strings
           tstart    first timestep
           tend      last timestep
           msg       error message
   Return: 0 ok, non-0 error
*/
int adios_readutil_open( const char *file, int64_t *fh, 
                         int *ngroups, adios_readutil_namelist **groupnames,
                         int *tstart, int *tend, char *msg);

/* Get names of variables and attributes of a group
   Input:  file handler
           group name
   Output: 
           nvars     number of variables
           varnames  linked list of strings
           nattrs    number of variables
           attrnames linked list of strings
           msg error message
   Return: 0 ok, non-0 error
*/
int adios_readutil_groupinfo (int64_t fh, const char *groupname, int64_t *gh, 
                              int *nvars, adios_readutil_namelist **varnames, 
                              int *nattrs, adios_readutil_namelist **attrnames, 
                              char *msg);

/* Get information about a variable
   Input: 
          gh        group handler
          path      var name
   Output:
          ndims     number of dimensions
          dims      dimension array (space should be provided by the caller!)
          type      adios type of variable
          timed     1 if variable has timesteps, 0 otherwise
          msg       error message
   Return: 
          0 OK, < 0 on error
*/
int adios_readutil_getvarinfo( int64_t gh, const char *path,
                               int *ndims, int *dims, int *type, int *timed, char *msg);

/* Get information about an attribute
   Input: 
          gh        group handler
          path      attr name
   Output:
          type      adios type of attribute
          size      string type: length of string, other types: 1
          msg       error message
   Return: 
          0 OK, < 0 on error
*/
int adios_readutil_getattrinfo( int64_t gh, const char *path,
                                int *type, int *size, char *msg);

/********************************/
/*   Functions for adiosclosec  */
/********************************/

void adios_readutil_fclose(int64_t fh);
void adios_readutil_gclose(int64_t gh);

/********************************/
/*   Functions for adiosreadc   */
/********************************/

/* Open an adios file. 
   Input:  file name
   Output: fh file handler
           msg error message
   Return: 0 ok, non-0 error
*/
int adios_readutil_fopen( const char *file, int64_t *fh, char *msg);

/* Open a group.
   Input:  group name
           fh file handler
   Output: gh group handler
           msg error message
   Return: 0 ok, non-0 error
*/
int adios_readutil_gopen_byname( const char *groupname, int64_t fh, int64_t *gh, char *msg);
int adios_readutil_gopen_byindex( int32_t groupidx, int64_t fh, int64_t *gh, char *msg);

/* name of an adios type (to print in logs) */
const char *adios_readutil_type_to_string(int type); 

/* interpret negative timesteps as -1 last step, -2 previous step, etc.
   positive timesteps are returned without modification
*/
int32_t adios_readutil_calctimestep(int64_t fh, int32_t timestep);

/* Get information about a variable/attribute.
   Input: 
          gh        group handler
          path      var/attr name
   Output:
          ndims     number of dimensions
          dims      dimension array (space should be provided by the caller!)
          type      adios type of variable
          msg       error message
          isvar     0 attribute, 1 variable was found for 'path'
   Return: 
          0 OK, < 0 on error
*/
int adios_readutil_getdatainfo( int64_t gh, const char *path, 
                                int *ndims, int *dims, int *type, char *msg, int *isvar);

/* Read in a variable/attribute.
   Input: 
          gh        group handler
          path      var/attr name
          isvar     1 read variable, 0 read attribute
          offsets   offset array for slicing
          counts    count array for counting
          timestep  (numbering starts from 1)
   Output:
          data (space should be provided by the caller!)
          msg error message
   Return: 
          the number of bytes read in, < 0 on error
*/

int adios_readutil_readdata(int64_t gh, const char *path, int isvar,
                            const int *offset, const int *count, int timestep,
                            void *data, char *msg);


/* Get the list of attributes that are directly below path, i.e.  attrpath = path/name.
   Input: 
         gh group handler
         path
   Output:
         attribute list, space will be allocated inside
   Return:
         number of attributes found
*/
int adios_readutil_getdirectattributes(int64_t gh, const char* path, adios_readutil_namelist **attrs);

/* Free list of strings returned by adios_readutil_getdirectattributes */
void adios_readutil_freenamelist(adios_readutil_namelist *namelist);

/* Set verbosity level for forthcoming calls this module */
void adios_readutil_setverbosity(int level);
