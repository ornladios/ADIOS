#ifndef __READ_XPMEM_H__
#define __READ_XPMEM_H__
// conditional libraries
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define BYTE_ALIGN 8
#define MINIFOOTER_SIZE 28


#define BUFREAD8(b,var)  do{var = (uint8_t) *(b->buff + b->offset);     \
		b->offset += 1;}while(0);

#define BUFREAD16(b,var) do{var = *(uint16_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_16(var); \
                         b->offset += 2;}while(0);

#define BUFREAD32(b,var) do{var = *(uint32_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_32(var); \
                         b->offset += 4;}while(0);

#define BUFREAD64(b,var) do{var = *(uint64_t *) (b->buff + b->offset);  \
                         if (b->change_endianness == adios_flag_yes) \
                             swap_64(var); \
                         b->offset += 8;}while(0);


typedef struct _xpmem_read_data
{
	xpmem_segid_t data_segid;
	xpmem_apid_t data_apid;
	shared_data *pg;
	char *data;
	uint64_t dsize;
	shared_data debug;
}xpmem_read_data;

typedef struct _xpmem_read_file
{
	uint32_t timestep;
	struct adios_bp_buffer_struct_v1 *b;
	struct bp_minifooter mfooter; 
	xpmem_read_data *fp;
	BP_PROC *bp;
	BP_FILE *fh;
}xpmem_read_file;

xpmem_read_data* fp = NULL;
static int show_hidden_attrs = 0; // don't show hidden attr by default

/* variety of functions taken from bp_utils to make them work with
 * shared memory style method
 */

static int print_shared(shared_data *pg)
{
		fprintf(stderr, "version %u size = %ll readcount = %u offset = %u reader = %u buffer=%p\n", pg->version,
	        pg->size,
	        pg->readcount,
	        pg->offset,
	        pg->reader,
	        pg->buffer);
}

static int xp_seek_to_step (ADIOS_FILE * fp, int tostep, int show_hidden_attrs)
{
	xpmem_read_file *xf = (xpmem_read_file*)fp->fh;
	BP_FILE *fh = xf->fh;

	int j, k, t, allstep;
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    struct adios_index_attribute_struct_v1 * attr_root;
    uint64_t i;
    int lenpath, lenname;

    /* Streaming starts with step 0. However, time index in BP file
     * starts with 1. If 'tostep' is -1, that means we want to get all steps.
     * If not, we seek to the specified step.
     */

    //TIME FIX
    // if (tostep == -1)
    // {
    //     allstep = 1;
    // }
    // else
    // {
    //     allstep = 0;
    //     t = get_time (var_root, tostep);
    // }

    allstep = 1;
    
    /* Prepare vars */
    fp->nvars = 0;
//    var_root = fh->vars_root;

    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (allstep || (!allstep && var_root->characteristics[i].time_index == t))
            {
                fp->nvars++;
                break;
            }
        }

        var_root = var_root->next;
    }

    fp->var_namelist = (char **) malloc (sizeof (char *) * fp->nvars);
    p->varid_mapping = (int *) malloc (fp->nvars * 4);

    var_root = fh->vars_root;
    j = 0;
    k = 0;
    while (var_root)
    {
        for (i = 0; i < var_root->characteristics_count; i++)
        {
            if (allstep || (!allstep && var_root->characteristics[i].time_index == t))
            {
                /* From 1.6, relative and full path (starts with /) are handled separately in search */
                // Full name of variable: concatenate var_path and var_name
                lenpath = strlen(var_root->var_path);
                lenname = strlen(var_root->var_name);
                if (lenpath > 0) {
                    fp->var_namelist [j] = (char *) malloc (lenname + lenpath + 1 + 1);
                                                                    // extra / and ending \0
                    strcpy(fp->var_namelist[j], var_root->var_path);
                    if (var_root->var_path[lenpath-1] != '/') {
                        fp->var_namelist[j][lenpath] = '/';
                        lenpath++;
                    }
                    strcpy(&(fp->var_namelist[j][lenpath]), var_root->var_name);
                }
                else {
                    fp->var_namelist[j] = (char *) malloc (lenname+1); 
                    strcpy(fp->var_namelist[j], var_root->var_name);
                }
                //printf ("Seek to step: Variable %d full path is [%s]\n", j, fp->var_namelist[j]);
                p->varid_mapping[j] = k;
                j++;

                break;
            }
        }

        k++;
        var_root = var_root->next;
    }

    /* Prepare attrs */
    fp->nattrs = 0;
    attr_root = fh->attrs_root;

    while (attr_root)
    {
        if (!show_hidden_attrs && strstr (attr_root->attr_path, "__adios__"))
        {
        }
        else
        {
            for (i = 0; i < attr_root->characteristics_count; i++)
            {
                if (allstep || (!allstep && attr_root->characteristics[i].time_index == t))
                {
                    fp->nattrs++;
                    break;
                }
            }
        }

        attr_root = attr_root->next;
    }

    fp->attr_namelist = (char **) malloc (sizeof (char *) * fp->nattrs);

    attr_root = fh->attrs_root;
    j = 0;
    while (attr_root)
    {
        if (!show_hidden_attrs && strstr (attr_root->attr_path, "__adios__"))
        {
        }
        else
        {
            for (i = 0; i < attr_root->characteristics_count; i++)
            {
                if (allstep || (!allstep && attr_root->characteristics[i].time_index == t))
                {
                    // Full name of attribute: concatenate attr_path and attr_name
                    lenpath = strlen(attr_root->attr_path);
                    lenname = strlen(attr_root->attr_name);
                    if (lenpath > 0) {
                        fp->attr_namelist [j] = (char *) malloc (lenname + lenpath + 1 + 1);
                                                                    // extra / and ending \0
                        strcpy(fp->attr_namelist[j], attr_root->attr_path);
                        if (attr_root->attr_path[lenpath-1] != '/') {
                            fp->attr_namelist[j][lenpath] = '/';
                            lenpath++;
                        }
                        strcpy(&(fp->attr_namelist[j][lenpath]), attr_root->attr_name);
                    }
                    else {
                        fp->attr_namelist[j] = (char *) malloc (lenname+1); 
                        strcpy(fp->attr_namelist[j], attr_root->attr_name);
                    }
                    //printf ("Seek to step: Attribute %d full path is [%s], path=[%s], name=[%s]\n", 
                    //        j, fp->attr_namelist[j], attr_root->attr_path, attr_root->attr_name);
                    j++;

                    break;
                }
            }
        }

        attr_root = attr_root->next;
    }

    fp->current_step = tostep;

    return 0;
}


static int xp_read_open(BP_FILE *fh, xpmem_read_data *f)
{
	//we don't need MPI_Comm or file name here
	struct adios_bp_buffer_struct_v1 *b;
	struct bp_minifooter * mh = &fh->mfooter;
	uint64_t attrs_end;
	uint64_t footer_size;	

	fh->sfh = NULL;
	fh->b = malloc(sizeof(struct adios_bp_buffer_struct_v1));
	fh->mfooter.file_size = f->dsize;
	fh->mpi_fh = 0;
	fh->fname = strdup("xpmem");
	fh->comm = 0;
	b = fh->b;
	
	log_info("file size = %llu\n", fh->b->file_size);

	adios_buffer_struct_init(fh->b);
	fh->b->file_size = f->dsize;

	//now we read the minifooter

	//allocate memory if we don't have it already
	if(!fh->b->buff)
	{
		bp_alloc_aligned(fh->b, MINIFOOTER_SIZE);
		memset(fh->b->buff, 0, MINIFOOTER_SIZE);
		fh->b->offset = 0;
	}

	memcpy(b->buff, &f->data[f->dsize - MINIFOOTER_SIZE],
	       MINIFOOTER_SIZE);

	b->offset = MINIFOOTER_SIZE - 4;

	adios_parse_version(b, &mh->version);
	mh->change_endianness = b->change_endianness;

	log_info("mh->version = %d\n", mh->version);

	//reset the offset for b
	b->offset = 0;

	//now we identify where each of the offsets are

	//pg index offset
	BUFREAD64(b, b->pg_index_offset);
	mh->pgs_index_offset = b->pg_index_offset;
    if (b->pg_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. PG index offset (%lld) > file size (%lld)\n",
                b->pg_index_offset, b->file_size);
        return -1;
    }

    //vars_index_offset
    BUFREAD64(b, b->vars_index_offset)
    mh->vars_index_offset = b->vars_index_offset;
    // validity check  
    if (b->vars_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) > file size (%lld)\n",
                b->vars_index_offset, b->file_size);
        return -1;
    }
    if (b->vars_index_offset < b->pg_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Variable index offset (%lld) < PG index offset (%lld)\n",
                b->vars_index_offset, b->pg_index_offset);
        return -1;
    }


    //attrs_index_offset
    BUFREAD64(b, b->attrs_index_offset)
    mh->attrs_index_offset = b->attrs_index_offset;
    // validity check  
    if (b->attrs_index_offset > b->file_size) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) > file size (%lld)\n",
                b->attrs_index_offset, b->file_size);
        return -1;
    }
    
    if (b->attrs_index_offset < b->vars_index_offset) {
        adios_error (err_file_open_error,
                "Invalid BP file detected. Attribute index offset (%lld) < Variable index offset (%lld)\n",
                b->attrs_index_offset, b->vars_index_offset);
        return -1;
    }

    //set the values in the buffer
    b->end_of_pgs = b->pg_index_offset;
    b->pg_size = b->vars_index_offset - b->pg_index_offset;
    b->vars_size = b->attrs_index_offset - b->vars_index_offset;
    
    attrs_end = b->file_size - MINIFOOTER_SIZE;
    b->attrs_size = attrs_end - b->attrs_index_offset;

    log_debug("offsets are %llu, %llu %llu\n",
              b->pg_index_offset, b->vars_index_offset,
              b->attrs_index_offset);

    //now figure out the rest of the footer

    //footer starts at pgs_index_offset and goes to end of file
    footer_size =  mh->file_size - mh->pgs_index_offset;
    bp_realloc_aligned (b, footer_size);

    //copy the footer into b starting from pgs_index_offset
    //to end of buffer
    
    memcpy(b->buff, &f->data[mh->pgs_index_offset], footer_size);

    //reset offset so we can read from buffer
    b->offset = 0;	

    bp_parse_pgs(fh);
    bp_parse_vars(fh);
    bp_parse_attrs(fh);
    
	
	return 0;
}

static ADIOS_VARINFO * xp_inq_var_byid (const ADIOS_FILE * fp, int varid)
{
	xpmem_read_file *xf = (xpmem_read_file *)fp->fh;
	BP_PROC * p = xf->bp;
	BP_FILE * fh = xf->fh;
	
    ADIOS_VARINFO * varinfo;
    int file_is_fortran, size, i;
    struct adios_index_var_struct_v1 * v;

    adios_errno = 0;

    v = bp_find_var_byid (fh, varid);

    varinfo = (ADIOS_VARINFO *) malloc (sizeof (ADIOS_VARINFO));
    assert (varinfo);

    /* Note: set varid as the real varid.
       Common read layer should convert it to the perceived id after the read me
thod returns.
    */
    varinfo->varid = varid;
    varinfo->type = v->type;
    file_is_fortran = is_fortran_file (fh);

    assert (v->characteristics_count);

    // Bugfix for block test (block.c). Actually changes
    // are actually inside bp_get_and_swap_dimensions.
    // For streaming mode, varinfo is built per steps.
    // Q. Liu 08/27/2014
    bp_get_and_swap_dimensions (fp, v, file_is_fortran,
                                &varinfo->ndim, &varinfo->dims,
                                &varinfo->nsteps,
                                file_is_fortran != futils_is_called_from_fortran()
                               );
    if (p->streaming)
    {
        varinfo->nsteps = 1;
    }

    // set value for scalar
    if (v->characteristics [0].value)
    {
        i = 0;

        if (p->streaming)
        {
            int time = fp->current_step + 1;
            i = 0;
            while (i < v->characteristics_count && v->characteristics[i].time_index != time)
            {
                i++;
            }

            if (i >= v->characteristics_count)
            {
                // shouldn't be here
            }
        }
        else
        {
            // keep i as 0
        }

        size = bp_get_type_size (v->type, v->characteristics [i].value);
        varinfo->value = (void *) malloc (size);
        assert (varinfo->value);

        memcpy (varinfo->value, v->characteristics [i].value, size);
    }
    else
    {
        varinfo->value = NULL;
    }

    varinfo->global = is_global_array (&(v->characteristics[0]));

    varinfo->nblocks = get_var_nblocks (v, varinfo->nsteps);

    assert (varinfo->nblocks);

    varinfo->sum_nblocks = (!p->streaming ? v->characteristics_count : varinfo->nblocks[0]) ;
    varinfo->statistics = 0;
    varinfo->blockinfo = 0;
    varinfo->meshinfo = 0;

    return varinfo;
}

#endif
