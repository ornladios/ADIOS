#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include "netcdf.h"
#include "adios_types.h"
#include "adios_bp_v1.h"
#include "adios_transport_hooks.h"
#include "adios_internals.h"
#include "bp2ncd.h"
#define ERR(e){if(e){printf("Error:%s\n",nc_strerror(e));return 2;}}
#define STRLEN 256
int main(int argc, char **argv)
{
	char 	*bp_fname = NULL;
	char 	*ncd_fname = NULL;
	uint32_t  from, to;
	int 	rc;
	enum verbose_level  verb;
	struct time_slice_struct time_slice;
	if (parse_cmdline(argc, argv, &bp_fname, &ncd_fname, &from, &to,  &verb)) {
		return -1;
	}
	if (ncd_fname==NULL) {
		int size;
		size = strlen(bp_fname);
		if (!strcmp(bp_fname+size-3,".bp")) {
			ncd_fname = strdup(bp_fname);
			ncd_fname[size - 2] = 'n';	
			ncd_fname[size - 1] = 'c';
		}
		else {
			ncd_fname = (char *)malloc(size+3);
		     	sprintf (ncd_fname,"%s%s",bp_fname,".nc");	
		}
	}
	fprintf(stderr, "bp2ncd\n\tconvert %s to %s\n", bp_fname, ncd_fname);
	if (to==0)
		fprintf(stderr, "\tfrom time step %d to the end\n", from);
	else
	fprintf(stderr, "\tfrom time step %d to %d\n", from, to);
	time_slice.from=from;
	time_slice.to=to;
	rc=makencd(bp_fname, ncd_fname, from, to);
	if (ncd_fname)
		free(ncd_fname);
	return rc;
}

int makencd(char *bp_fname, char *ncd_fname,int from, int to) 
{
	int size, i, j, rc;
	uint32_t version = 0;
	uint64_t element_size = 0;
	struct adios_bp_element_struct *element=NULL;
	struct adios_bp_buffer_struct_v1 *b, *b0;
	
	int ncid, retval;
	nc_create (ncd_fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
	
	b = (struct adios_bp_buffer_struct_v1 *) 
		malloc(sizeof(struct adios_bp_buffer_struct_v1));
	b0 = (struct adios_bp_buffer_struct_v1 *) 
		malloc(sizeof(struct adios_bp_buffer_struct_v1));
	adios_buffer_struct_init(b);
  	rc = adios_posix_open_read_internal(bp_fname,"",b);	
	adios_posix_read_version(b);
	adios_parse_version(b,&version);

	struct adios_index_process_group_struct_v1 *pg_root=0;
	struct adios_index_process_group_struct_v1 *pg=0;
	struct adios_index_var_struct_v1 *vars_root=0;
	struct adios_index_attribute_struct_v1 *attrs_root=0;
		
	adios_posix_read_index_offsets(b);
	adios_parse_index_offsets_v1(b);
	
	adios_posix_read_process_group_index(b);
	adios_parse_process_group_index_v1(b, &pg_root);
	
	adios_posix_read_vars_index(b);
	adios_parse_vars_index_v1(b, &vars_root);
	
	adios_posix_read_attributes_index(b);
	adios_parse_attributes_index_v1(b, &attrs_root);

	struct 	adios_process_group_header_struct_v1 pg_header;
	struct 	adios_vars_header_struct_v1 vars_header;
	struct 	adios_attributes_header_struct_v1 attrs_header;

	struct 	adios_var_header_struct_v1 var_header;
	struct 	adios_var_payload_struct_v1 var_payload;
	struct 	adios_attribute_struct_v1 attribute;
	struct 	var_dims_struct *var_dims = NULL; 
	int    	time_index_flag = 0;
	uint64_t element_num=0, var_dims_count=0, length_of_var;
	pg = pg_root;
	int vars_count=0, attrs_count=0;
	while(vars_root) {
		vars_count++;
		vars_root = vars_root->next;
	}
	while(attrs_root) {
		attrs_count++;
		attrs_root = attrs_root->next;
	}
	printf("TABLE size: %d, %d\n",vars_count,attrs_count);
	var_dims = (struct var_dims_struct *) 
		    malloc ((vars_count+attrs_count+1)*
				sizeof(struct var_dims_struct));
	if (pg-> time_index_name) {
		time_index_flag = 1;
		nc_redef(ncid);
		nc_def_dim(ncid,pg->time_index_name,
				NC_UNLIMITED, &var_dims[0].nc_dimid);
		nc_enddef(ncid);
		var_dims[0].id = 0;
		var_dims[0].rank = pg->time_index;
		strcpy(var_dims[0].varname,
				pg->time_index_name);
	}
	else{
		var_dims[0].id = 0;
		var_dims[0].rank = 0; 
		var_dims[0].nc_dimid = -1;
		strcpy(var_dims[0].varname,"");
	}
	var_dims_count = 1;
	while (pg) {
		b->read_pg_offset = pg->offset_in_file;
		if (pg->next) {
			b->read_pg_size = pg->next->offset_in_file
				- pg->offset_in_file;
		}
		else {
			b->read_pg_size = b->pg_index_offset
				- pg->offset_in_file;
		}	

		adios_posix_read_process_group(b);
		adios_parse_process_group_header_v1(b,&pg_header);
		adios_parse_vars_header_v1(b,&vars_header);
/*
		fprintf(stderr,"%s","XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
		fprintf(stderr,"%s","-----------------------\n");
		fprintf(stderr,"number of vars: %d\n",
				vars_header.count);
		fprintf(stderr,"%s","-----------------------\n");
*/
		if ((pg->time_index_name) &&
	  	    (pg->time_index<from || 
		     pg->time_index>to && to!=0)) {
			pg=pg->next;
		}
		else {	
		var_dims[0].rank = pg->time_index-from+1;
		if (element_num%2 == 0) {
			for (i=0;i<vars_header.count;i++) {
				var_payload.payload = 0;
				copy_buffer(b0, b);
				length_of_var = *(uint64_t*)(b->buff+b->offset);
				adios_parse_var_data_header_v1(b,&var_header);
/*
				fprintf(stderr, "var %d: %s %s\n",i,
						var_header.path,
						var_header.name);
				fprintf(stderr, "\tcount:%d\n",var_dims_count);
*/
				add_vardims_var(ncid,
						var_dims,
						&var_dims_count,
						&var_header, 
						&var_payload, 
						b);	

				copy_buffer(b,b0);
				b->offset = b->offset+length_of_var;
			}

			adios_parse_attributes_header_v1(b,&attrs_header);
/*
			printf("-----------------------\n");
			fprintf(stderr,"number of attrs: %d\n", 
					attrs_header.count);
			printf("-----------------------\n");
 */
			for (i=0;i<attrs_header.count;i++) {

				adios_parse_attribute_v1(b,&attribute);
				add_vardims_attribute(ncid,
						var_dims,
						&var_dims_count,
						&attribute);
			}
		}
		else {
			for (i=0;i<vars_header.count;i++) {
				var_payload.payload = 0;
				adios_parse_var_data_header_v1(b,&var_header);
				if (!var_payload.payload)
					var_payload.payload = 
						malloc(var_header.payload_size);
				adios_parse_var_data_payload_v1 (b, &var_header,
						&var_payload,
						var_header.payload_size
						);
				// scalar and is not used for variable dimension 
				if (var_header.dims) {
					//fprintf(stderr,"write var:%s\n",var_header.name);
					ncd_dataset(ncid, 
							&var_header,
							var_payload.payload, 
							var_dims, 
							var_dims_count);
				}
				else if (var_header.is_dim == adios_flag_no) {
					ncd_scalar(ncid,
							&var_header,
							var_payload.payload);
				}

				if (var_payload.payload)
					free(var_payload.payload);
			}
			adios_parse_attributes_header_v1(b,&attrs_header);
			for (i=0;i<attrs_header.count;i++) {
				adios_parse_attribute_v1(b,&attribute);
				ncd_attr(ncid, 
						&attribute,
						var_dims,
						var_dims_count);

			}
			pg = pg->next;
		}//else
		element_num++;
		element_num = element_num%2;
	} //from to
	} //while
	if (var_dims) 
		free(var_dims);
	nc_redef(ncid);
	retval=nc_put_att_int(ncid,NC_GLOBAL,"time_index_from",NC_INT,1,&from);
	ERR(retval);
	retval=nc_close(ncid);
	nc_enddef(ncid);
	ERR(retval);
	free(b0);
	free(b);
	return 0;	
}

int add_vardims_var(int ncid,
		struct var_dims_struct* var_dims,
		int	*var_dims_count, 
		struct adios_var_header_struct_v1 *var_header, 
		struct adios_var_payload_struct_v1 *var_payload, 
		struct adios_bp_buffer_struct_v1 *b) 
{
	int i,idx=*var_dims_count;
	char fullname[STRLEN];
	if (var_header->is_dim == adios_flag_yes && var_header->dims == 0) {
		for(i=0;i<idx;i++) {
			if(var_header->id == var_dims[i].id)
				return 0;
		}	
		ncd_gen_name(fullname, var_header->path, 
				var_header->name);
		if (!var_payload->payload) {
			var_payload->payload = malloc(var_header->payload_size);
			adios_parse_var_data_payload_v1(b, 
					var_header,
					var_payload,
					var_header->payload_size
					);
			var_dims[idx].id = var_header->id;
			var_dims[idx].rank = *(unsigned int *)
				var_payload->payload;
			if (var_payload->payload) {
				free(var_payload->payload);
				var_payload->payload = NULL;
			}
		}
		var_dims[idx].nc_dimid = -1;
		strcpy(var_dims[idx].varname,fullname);
		*var_dims_count=idx+1;	
		}
		return 0;
	}
	int add_vardims_attribute(int ncid,
			struct var_dims_struct* var_dims,
			int   *var_dims_count,
			struct adios_attribute_struct_v1 *attribute)
	{
		int i;
		int idx=*var_dims_count;
		char fullname[STRLEN];
		for(i=0;i<idx;i++) {
			if(attribute->id == var_dims[i].id)
				return 0;
		}
		ncd_gen_name(fullname, attribute->path, 
				attribute->name);
		//fprintf(stderr, "\tcount:%d attribute:%s\n",idx,fullname);
		if (attribute->is_var == adios_flag_yes) {
			var_dims[idx].id = attribute->var_id;
			var_dims[idx].rank = 0;
			for (i=0; i<idx;i++) {
				if (attribute->var_id == var_dims[i].id) {
					var_dims[idx].id = attribute->id;
					var_dims[idx].rank = var_dims[i].rank;
					i = idx;
				}
			}
		}
		else {
			var_dims[idx].id = attribute->id;
			var_dims[idx].rank = 0;
			assign_value_ulonglong(&(var_dims[idx].rank),
					attribute->type, 
					attribute->value);
		}
		var_dims[idx].nc_dimid = -1;
		strcpy(var_dims[idx].varname,fullname);
		*var_dims_count=idx+1;	
		return 0;
	}

	int assign_value_ulonglong(uint64_t *rank, 
			enum ADIOS_DATATYPES type, 
		void *value)
{
	if (!value) {
		fprintf(stderr, "data value is NULL in "
				"assign_value_ulonglong()\n");
		return -1;
	}	
	switch(type) {
		case adios_unsigned_short:
			*rank = (uint64_t) *((unsigned short *) value);
			break;
		case adios_unsigned_integer:
			*rank = (uint64_t) *((unsigned int *) value);
			break;
		case adios_unsigned_long:
			*rank = (uint64_t) *((unsigned long*) value);
			break;
		case adios_short:
			*rank = (uint64_t) *((short*) value);
			break;
		case adios_integer:
			*rank = (uint64_t) *((int *) value);
			break;
		case adios_long:
			*rank = (uint64_t) *((long *) value);
			break;
		case adios_byte:
		case adios_unsigned_byte:
		case adios_string:
			break;
		case adios_complex:
		case adios_double_complex:
			fprintf(stderr, "Error in mapping ADIOS Data "
					"Types to HDF5: complex not supported yet.\n");
			break;
		case adios_unknown:
		default:
			fprintf(stderr, "Error in mapping ADIOS Data "
					"Types to HDF5: unknown data type.\n");
	}
	return 0;
}

void print_usage() 
{
	printf("bp2ncd:\n");
	printf("NAME\n");
	printf("\tbp2ncd - Convert a bp file to a NetCDF file\n");
	printf("SYNOPSIS\n");
	printf("\tbp2ncd [OPTION] ... bp_file [ncd_file]\n");
	printf("DESCRIPTION\n");
	printf("\t-f, --from\n");
	printf("\t\twhich timestep the dataset is started from\n");
	printf("\t-t, --to\n");
	printf("\t\twhich timestep the dataset is ended to\n");
	printf("\t--help\n");
	printf("\t\tdisplay this help and exit\n");
	printf("\t--version\n");
	printf("\t\toutput version information and exit\n");
	printf("AUTHOR\n");
	printf("\tWritten by Chen Jin\n");
	printf("REPORTING BUGS\n");
	printf("\tReport bugs to <jinc@ornl.gov>\n");
	printf("COPYRIGHT\n");
}
int parse_cmdline(int argc,
		char **argv,
		char **bp_fname,
		  char **ncd_fname,
		  uint32_t *from_timeindex,
		  uint32_t *to_timeindex,
		  enum verbose_level *verb) 
{
	int i=1;
	*verb = NO_INFO;
 	*from_timeindex = 1;
	*to_timeindex= 0;
        int found_bp_file = 0;
	int found_ncd_file = 0;
	while (i < argc) {
 		if (!strcmp(argv[i], "--verbose") || !strcmp(argv[i], "-V")) {
			*verb = atoi(argv[i+1]);
			i++;
		}
		else if (!strcmp(argv[i], "--from") || !strcmp(argv[i], "-f")) {
			if (i+1<argc) {
				*from_timeindex = atoi(argv[i+1]);
				i++;	
			}
		}
		else if (!strcmp(argv[i], "--to") || !strcmp(argv[i], "-t")) {
			if (i+1<argc) {
				*to_timeindex = atoi(argv[i+1]);
				i++;	
			}
		}
		else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {
			print_usage();
			exit(0);
		}
		else if(!found_bp_file) {
			*bp_fname = argv[i];
			found_bp_file = 1;	
		}
		else if(found_bp_file && !found_ncd_file) {
			*ncd_fname = argv[i];
			found_ncd_file = 1;	
		}
		else {
			fprintf(stderr, 
				"Error in parsing command line: unknown argument %s\n\n", argv[i]);
            		print_usage();
			return -1;
		}
		i++;
	}
	if (!found_bp_file) {
		fprintf(stderr, "Error in parsing commnad line: bp_file is not provided\n");
		print_usage();
		return -1;	
	}	
	return 0;	
	
}

int ncd_gen_name (char *fullname, char *path, char *name) {
    int i;
    char *new_path = strdup (path);
    if ( path[0] == '/')
         new_path=new_path+1;

    for ( i = 0; i < strlen (new_path); i++) {
        if ( new_path[i] == '[' || new_path[i] == ']' || new_path[i] == '/' || new_path[i] == '\\')
            new_path[i] = '_';
    }
    if (*new_path != '\0') {
        if (new_path[i-1]!='_') {
            if (strcmp(name,"") )
                sprintf (fullname, "%s_%s", new_path, name);
            else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0';
            }
        }
        else
            if (strcmp(name,"") )
                sprintf (fullname, "%s%s", new_path, name);
            else {
                strcpy (fullname,new_path);
                fullname [strlen(fullname)] = '\0';
            }
    }
    else
        strcpy (fullname, name);
    return 0;
}

int ncd_scalar(int ncid, 
	       struct adios_var_header_struct_v1 *var_header,
	       void *var_payload) 
{
	int rank = 1, varid, retval, type=var_header->type;
   	void *val = var_payload;
	static int dimid = -1;
	char fullname[STRLEN];

	ncd_gen_name(fullname, var_header->path, 
			var_header->name);

	if (dimid==-1) {
		retval = nc_def_dim (ncid, "one", 1, &dimid);
		ERR(retval);
	}
	else {
		nc_redef(ncid);
		nc_inq_varid(ncid, fullname, &varid);
	}
	rank = 1;

	switch (type) {
		case adios_real:
			if (varid < 0 ) {
				fprintf(stderr,"\tncd-scalar-real: %d %d %s\n",
						&dimid,varid, fullname);
				retval=nc_def_var(ncid,fullname,NC_FLOAT,
						rank,&dimid,&varid);
				ERR(retval);
			}
			retval=nc_enddef(ncid);
			ERR(retval);
			retval=nc_put_var_float(ncid,varid,val);
			break;
		case adios_double:
			if (varid < 0 ) {
				fprintf(stderr,"\tncd-scalar: %d %d %s\n",
						&dimid,varid, fullname);
				retval=nc_def_var(ncid,fullname,NC_DOUBLE,
						rank,&dimid,&varid);
				ERR(retval);
				retval=nc_enddef(ncid);
				ERR(retval);
			}
			retval=nc_enddef(ncid);
			retval=nc_put_var_double(ncid,varid,val);
			break;
		case adios_long:
			if (varid < 0 ) {
				fprintf(stderr,"\tncd-scalar: %d %d %s\n",
						&dimid,varid, fullname);
				retval=nc_def_var(ncid,fullname,NC_LONG,
						rank,&dimid,&varid);
				ERR(retval);
			}
			retval=nc_enddef(ncid);
			ERR(retval);
			retval=nc_def_var(ncid,fullname,NC_LONG,
					rank,&dimid,&varid);
			retval=nc_enddef(ncid);
			retval=nc_put_var_long(ncid,varid,val);
			break;
		case adios_integer:
			if (varid < 0 ) {
				retval=nc_def_var(ncid,fullname,NC_INT,
						rank,&dimid,&varid);
				ERR(retval);
			}
			retval=nc_enddef(ncid);
			ERR(retval);
			retval=nc_put_var_int(ncid,varid,val);
			ERR(retval);
			break;
		default:
			retval=nc_enddef(ncid);
	}
	return 0;
}

int ncd_dataset(int 	ncid, 
		struct 	adios_var_header_struct_v1 *var_header, 
		void 	*val, 
		struct 	var_dims_struct *var_dims, 
		int 	var_dims_count)
{
	int  	maxrank=0, rank, time_dimrank, i,j, retval;
	struct 	adios_dimension_struct_v1 *dims = var_header->dims;
	int   	*dimids;
	size_t 	*start, *count;
	char 	fullname[STRLEN];
	int 	varid;

	ncd_gen_name(fullname,var_header->path,var_header->name);
	time_dimrank=-1;
	while (dims) {
		if (dims->dimension.time_index == adios_flag_yes) {
			time_dimrank = maxrank;
		}
		++maxrank;
		dims = dims->next;
	}

	dimids = (int *) malloc(sizeof(int)*maxrank);
	start = (size_t *) malloc(sizeof(size_t)*maxrank);
	count = (size_t *) malloc(sizeof(size_t)*maxrank);
	int idx=-1;
	dims = var_header->dims;

	if (dims->global_dimension.var_id != 0 
			|| dims->global_dimension.rank!=0) {
		printf("write global dataset!\n");
		return 0;
	}

	nc_redef(ncid);

	for(j=0;j<maxrank;j++) {
		if (time_dimrank==maxrank-1)
			rank = maxrank-1-j;
		else
			rank = j;
		start[rank]=0;
		if (dims->dimension.var_id != 0) {
			for (i=0;i<var_dims_count;i++) {
				if (dims->dimension.var_id == var_dims[i].id) {
					idx = i;
					i=var_dims_count;
				}
			}
			if (var_dims[idx].nc_dimid < 0 && idx<var_dims_count) {
				retval = nc_def_dim(ncid,
						    var_dims[idx].varname,
						    var_dims[idx].rank, 
						    &var_dims[idx].nc_dimid);
				ERR(retval);
			}
			dimids[rank]=var_dims[idx].nc_dimid;
			count[rank]=var_dims[idx].rank;
/*
			fprintf(stderr,"\tdims->varid:%d %d %s:%d\n",
				dims->dimension.var_id, idx,
					var_dims[idx].varname,
					var_dims[idx].nc_dimid);
			fprintf(stderr,"\tdim[%d]: c(%d):s(%d): dimid=%d\n"
				,rank
				,count[rank]
				,start[rank]
				,dimids[rank]
		       );
*/
		}
		else {
			if (dims->dimension.time_index == adios_flag_yes) {
				count[rank]=1;
				dimids[rank]=var_dims[0].nc_dimid;
				start[rank]=var_dims[0].rank-1;
/*
				fprintf(stderr,
						"\tset time-index dimension: %s %d\n",
						var_dims[0].varname,
						start[rank]);
*/
			}
			else {
/*
				fprintf(stderr,
						"\tset value dimension: %s\n",
						var_dims[i].varname);
*/
				count[rank]=dims->dimension.rank;
				dimids[rank]=var_dims[i].nc_dimid;
			}
		}
		dims=dims->next;
	}

/*
	for (rank=0;rank<maxrank;rank++)
		fprintf(stderr,"\t%s: dim[%d]: c(%d):s(%d): dimid=%d\n",
				fullname,
				rank,
				count[rank],
				start[rank],
				dimids[rank]
		       );
*/

	varid = -1;
	retval=nc_inq_varid(ncid, fullname, &varid);
	switch(var_header->type) {
		case adios_real:
			if ( varid<0) {                
				retval=nc_def_var(ncid,fullname,NC_FLOAT,
						maxrank,dimids,&varid);
				ERR(retval);
			}
			retval=nc_enddef(ncid);
			retval=nc_put_vara_float(ncid,varid,start,count,val);
			ERR(retval);
			break;
		case adios_double:
			if (varid<0) 
				retval=nc_def_var(ncid,fullname,NC_DOUBLE,
						maxrank,dimids,&varid);
			retval=nc_enddef(ncid); 
			ERR(retval);
			retval=nc_put_vara_double(ncid,varid,start,count,val);
			ERR(retval);
			break;
		case adios_long:
			if ( varid<0) 
				retval=nc_def_var(ncid,fullname,NC_LONG,
						maxrank,dimids,&varid);
			retval=nc_enddef(ncid);
			retval=nc_put_vara_long(ncid,varid,start,count,val);
			break; 
		case adios_unsigned_byte:
			if ( varid<0)
				retval=nc_def_var(ncid,fullname,NC_BYTE,
						maxrank,dimids,&varid);
			retval=nc_enddef(ncid);
			retval=nc_put_vara_uchar(ncid,varid,start,count,val);
			break;
		case adios_byte:
			if ( varid<0)
				retval=nc_def_var(ncid,fullname,NC_BYTE,
						maxrank,dimids,&varid);
			ERR (retval);
			retval=nc_enddef(ncid);
			retval=nc_put_vara_schar(ncid,varid,start,count,val);
			ERR (retval);
			break;
		case adios_integer:
			if (varid < 0) {
				retval = nc_def_var (ncid,fullname,NC_INT,
						maxrank,dimids,&varid);
			}
			retval = nc_enddef (ncid);
			ERR (retval);
			retval=nc_put_vara_int (ncid,varid,start,count,val);
			ERR (retval);
			break;
		default:
			retval=nc_enddef(ncid);
			break;
	}
	if (dimids)
		free(dimids);
	if (count)
		free(count);
	if (start)
		free(start);
	return 0;
}
int ncd_attr(int ncid, 
             struct adios_attribute_struct_v1 *attribute,
             struct var_dims_struct * var_dims,
	     int var_dims_count) {
	int i;
	char fullname[256];
	char *path = attribute->path;
	char *name = attribute->name;
	char *new_path;
	int  valid,retval,attid;

	ncd_gen_name (fullname, path, name);
	valid = -1;
	if (strcmp(path,"/")==0) {
		valid = NC_GLOBAL;
		strcpy(fullname, name);
	}
	else {
		ncd_gen_name (fullname, path, "");
		retval=nc_inq_varid(ncid,fullname,&valid);
		if(retval < 0)
			return; 
		else
			strcpy(fullname, name);
	}

	retval=nc_inq_attid(ncid,valid,fullname,&attid);

	if (retval == NC_NOERR ) 
		return;

	nc_redef(ncid);

	void *value = attribute->value;
	size_t len = 1; 
	enum ADIOS_DATATYPES type =  attribute->type;
	struct adios_var_header_struct_v1 var_header;
	struct adios_var_payload_struct_v1 var_payload;  
	uint64_t offset; 
	struct adios_index_var_struct_v1 * vars_root = 0;

	var_payload.payload = 0;

	if ( attribute->is_var == adios_flag_yes) {
		printf("%s(att) not implemented yet!\n", fullname);
		return;
	}
	//else
		//printf("\t      XML: ");  
 
	switch (type) {
		case adios_unsigned_byte:
			retval=nc_put_att_uchar(ncid,valid,fullname,NC_BYTE,len,value);
			break;
		case adios_byte:
			retval=nc_put_att_schar(ncid,valid,fullname,NC_BYTE,len,value);
			break;
		case adios_string:
			//printf("%s\n", (char *) value);    
			retval=nc_put_att_text(ncid,valid,fullname, strlen(value),value);
			break;
		case adios_short:
			retval=nc_put_att_short(ncid,valid,fullname,NC_SHORT,len,value);
			ERR(retval); 
			break;
		case adios_integer:
			//printf("%d\n", *((int *) value));    
			retval=nc_put_att_int(ncid,valid,fullname,NC_INT,len,value);
			break;
		case adios_long:
			retval=nc_put_att_long(ncid,valid,fullname,NC_LONG,len,value);
			break;
		case adios_real:
			retval=nc_put_att_float(ncid,valid,fullname,NC_FLOAT,len,value);
			break;
		case adios_double:
			retval=nc_put_att_double(ncid,valid,fullname,NC_DOUBLE,len,value);
			break;
		default:
			break;
	}
	ERR(retval);

	nc_enddef(ncid); 
}
