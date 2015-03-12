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

static int adios_wbidx_to_pgidx (const ADIOS_FILE * fp, read_request * r, int step_offset);
static ADIOS_VARCHUNK * read_var_bb (const ADIOS_FILE *fp, read_request * r);
static ADIOS_VARCHUNK * read_var (const ADIOS_FILE * fp, read_request * r);
static ADIOS_VARCHUNK * read_var_wb (const ADIOS_FILE * fp, read_request * r);

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
	BP_PROC *p = xf->bp;

	int j, k, t, allstep;
    struct adios_index_var_struct_v1 * var_root = fh->vars_root;
    struct adios_index_attribute_struct_v1 * attr_root;
    uint64_t i;
    int lenpath, lenname;

    /* Streaming starts with step 0. However, time index in BP file
     * starts with 1. If 'tostep' is -1, that means we want to get all steps.
     * If not, we seek to the specified step.
     */

    if (tostep == -1)
    {
        allstep = 1;
    }
    else
    {
        allstep = 0;
        t = get_time (var_root, tostep);
    }

    
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

static uint64_t get_req_datasize (const ADIOS_FILE * fp, read_request * r, struct adios_index_var_struct_v1 * v)
{
    ADIOS_SELECTION * sel = r->sel;
    uint64_t datasize = bp_get_type_size (v->type, "");
    int i, pgidx, ndims;

    if (sel->type == ADIOS_SELECTION_BOUNDINGBOX)
    {
        for (i = 0; i < sel->u.bb.ndim; i++)
        {
            datasize *=  sel->u.bb.count[i];
        }
    }
    else if (sel->type == ADIOS_SELECTION_POINTS)
    {
        datasize *= sel->u.points.npoints;
    }
    else if (sel->type == ADIOS_SELECTION_WRITEBLOCK)
    {
        //pgidx = adios_wbidx_to_pgidx (fp, r);
        // NCSU ALACRITY-ADIOS: Adding absoluet PG indexing
	    pgidx = sel->u.block.index;
        // pgidx = sel->u.block.is_absolute_index ?
        //             sel->u.block.index :
        //         adios_wbidx_to_pgidx (fp, r, 0);
        // NCSU ALACRITY-ADIOS: Adding sub-PG writeblock read support
        if (sel->u.block.is_sub_pg_selection) {
            datasize = sel->u.block.nelements;
        } else {
            // NCSU ALACRITY-ADIOS: This used to not be in this else block
            ndims = v->characteristics[pgidx].dims.count;
            for (i = 0; i < ndims; i++)
            {
                datasize *= v->characteristics[pgidx].dims.dims[i * 3];
            }
        }
    }

    return datasize;
}

/* This routine converts the write block index, which is of a particular step,
 * to the adios internal PG index.
 */
static int adios_wbidx_to_pgidx (const ADIOS_FILE * fp, read_request * r, int step_offset)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    int time, start_idx, stop_idx, c, idx;
    int mapped_varid, ridx;
    struct adios_index_var_struct_v1* v;

    if (r->sel->type != ADIOS_SELECTION_WRITEBLOCK)
    {
        return -1;
    }

    time = adios_step_to_time (fp, r->varid, r->from_steps + step_offset);
    mapped_varid = r->varid; //map_req_varid (fp, r->varid); // NCSU ALACRITY-ADIOS: Bugfix: r->varid has already been mapped
    v = bp_find_var_byid (fh, mapped_varid);

    start_idx = get_var_start_index (v, time);
    stop_idx = get_var_stop_index (v, time);
    if (start_idx < 0 || stop_idx < 0)
    {
        adios_error (err_no_data_at_timestep,
                     "No data at step %d\n",
                     r->from_steps);
    }

    ridx =  r->sel->u.block.index;
    c = -1;
    idx = start_idx;
    while (idx <= stop_idx)
    {
        if (v->characteristics[idx].time_index == time)
        {
            c++;
        }

        if (c < ridx)
        {
            idx++;
        }
        else
        {
            break;
        }
    }

    if (c != ridx)
    {
        log_debug ("Error in adios_wbidx_to_pgidx().\n");
    }

    return idx;
}

/* This routine processes a read request and returns data in ADIOS_VARCHUNK.
   If the selection type is not bounding box, convert it. The basic file reading
   functionality is implemented in read_var_bb() routine.
*/
static ADIOS_VARCHUNK * read_var (const ADIOS_FILE * fp, read_request * r)
{
    BP_PROC * p = GET_BP_PROC (fp);
    BP_FILE * fh = GET_BP_FILE (fp);

    int size_of_type;
    struct adios_index_var_struct_v1 * v;
    uint64_t i;
    read_request * nr;
    ADIOS_SELECTION * sel, * nsel;
    ADIOS_VARCHUNK * chunk;

    log_debug ("read_var()\n");
    sel = r->sel;
    v = bp_find_var_byid (fh, r->varid);

    switch (sel->type)
    {
        case ADIOS_SELECTION_BOUNDINGBOX:
            chunk = read_var_bb (fp, r);
            break;
        case ADIOS_SELECTION_POINTS:	        
	        break;
        case ADIOS_SELECTION_WRITEBLOCK:
            chunk = read_var_wb (fp, r);
            break;
        case ADIOS_SELECTION_AUTO:
            break;
        default:
            log_debug ("ADIOS selection type is wrong\n");
            break;
    }

    return chunk;
}

/* This routine reads in data for bounding box selection.
   If the selection is not bounding box, it should be converted to it.
   The data returned is saved in ADIOS_VARCHUNK.
 */
static ADIOS_VARCHUNK * read_var_bb (const ADIOS_FILE *fp, read_request * r)
{
	xpmem_read_file *xf = (xpmem_read_file *)fp->fh;	
	BP_PROC * p = GET_BP_PROC (fp);
	BP_FILE * fh = GET_BP_FILE (fp);
	xpmem_read_data *xd = xf->fp;

	ADIOS_SELECTION * sel;
	struct adios_index_var_struct_v1 * v;
	int i, j, t, time, nsteps;
	int64_t start_idx, stop_idx, idx;
	int ndim, file_is_fortran;
	uint64_t * dims, tmpcount;
	uint64_t ldims[32], gdims[32], offsets[32];
	uint64_t datasize, dset_stride,var_stride, total_size=0, items_read;
	uint64_t * count, * start;
	void * data;
	int dummy = -1, is_global = 0, size_of_type;
	uint64_t slice_offset, slice_size;
	MPI_Status status;
	ADIOS_VARCHUNK * chunk;
	struct adios_var_header_struct_v1 var_header;
	//struct adios_var_payload_struct_v1 var_payload;

//    log_debug ("read_var_bb()\n");
	file_is_fortran = is_fortran_file (fh);

	sel = r->sel;
	start = sel->u.bb.start;
	count = sel->u.bb.count;
	data = r->data;

	v = bp_find_var_byid (fh, r->varid);

	/* Get dimensions and flip if caller != writer language */
	/* Note: ndim below doesn't include time if there is any */
	// NCSU ALACRITY-ADIOS - Note: this function has been modified to return
	//   the "raw" dimensions (i.e., 1D byte array)
	bp_get_and_swap_dimensions (fp, v, file_is_fortran, &ndim, &dims, &nsteps, file_is_fortran);

	assert (ndim == sel->u.bb.ndim);
	ndim = sel->u.bb.ndim;

	/* Fortran reader was reported of Fortran dimension order so it gives counts and starts in that order.
	   We need to swap them here to read correctly in C order */
	if (futils_is_called_from_fortran ())
	{
		swap_order (ndim, start, &dummy);
		swap_order (ndim, count, &dummy);
	}

	/* items_read = how many data elements are we going to read in total (per timestep) */
	items_read = 1;
	for (i = 0; i < ndim; i++)
	{
		items_read *= count[i];
	}

	size_of_type = bp_get_type_size (v->type, v->characteristics [0].value);

//    log_debug ("read_var_bb: from_steps = %d, nsteps = %d\n", r->from_steps, r->nsteps);

	/* Note fp->current_step is always 0 for file mode. */
	for (t = fp->current_step + r->from_steps; t < fp->current_step + r->from_steps + r->nsteps; t++)
	{

		if (!p->streaming)
		{
			time = get_time (v, t);
		}
		else
		{
			// Fix: the assumption that for streaming mode, the time in file
			// always starts from 1 is not correct. So here we add fh->tidx_start to adjust
			// Q. Liu, 06/2013
			time = fh->tidx_start + t;
		}

		start_idx = get_var_start_index (v, time);
		stop_idx = get_var_stop_index (v, time);

		if (start_idx < 0 || stop_idx < 0)
		{
			adios_error (err_no_data_at_timestep,"Variable %s has no data at %d time step\n",
			             v->var_name, t);
			continue;
		}

		if (ndim == 0)
		{
			/* READ A SCALAR VARIABLE */
			/* Prepare slice_offset, slice_size and idx for the later macro:
			   MPI_FILE_READ_OPS1 and MPI_FILE_READ_OPS2
			*/
			idx = 0;
			slice_offset = v->characteristics[start_idx + idx].payload_offset;
			slice_size = size_of_type;

			if (v->type == adios_string)
			{
				// Note: strings are stored without \0 in file
				// size_of_type above includes \0 so decrease by one
				size_of_type--;
			}

			//copy the data
            
			memcpy ((char *)data, &xd->data[slice_offset], slice_size);

			if (fh->mfooter.change_endianness == adios_flag_yes)
			{
				log_error("endianness should not be changed\n");
			}

			if (v->type == adios_string)
			{
				// add \0 to the end of string
				// size_of_type here is the length of string
				// FIXME: how would this work for strings written over time?
				((char*)data)[size_of_type] = '\0';
			}

			//FIXME: Not clear why we read multiple timesteps like this. Ask norbert.
			//for scalars does it even make sense to do this?
			//and chunks will then point to the wrong pointer
			data = (char *)data + (v->type == adios_string?  size_of_type + 1 : size_of_type);
		}
		else
		{
			/* READ AN ARRAY VARIABLE */
			int * idx_table = (int *) malloc (sizeof (int) * (stop_idx - start_idx + 1));
			uint64_t write_offset = 0;

			// loop over the list of pgs to read from one-by-one
			for (idx = 0; idx < stop_idx - start_idx + 1; idx++)
			{
				int flag;
				datasize = 1;
				var_stride = 1;
				dset_stride = 1;
				idx_table[idx] = 1;
				uint64_t payload_size = size_of_type;

				is_global = bp_get_dimension_characteristics_notime (&(v->characteristics[start_idx + idx]),
				                                                     ldims, gdims, offsets, file_is_fortran);
				if (!is_global)
				{
					// we use gdims below, which is 0 for a local array; set to ldims here
					for (j = 0; j < ndim; j++)
					{
						gdims[j] = ldims[j];
					}
					// we need to read only the first PG, not all, so let's prevent a second loop
					stop_idx = start_idx;
				}
/*
  printf ("ldims   = "); for (j = 0; j<ndim; j++) printf ("%d ",ldims[j]); printf ("\n");
  printf ("gdims   = "); for (j = 0; j<ndim; j++) printf ("%d ",gdims[j]); printf ("\n");
  printf ("offsets = "); for (j = 0; j<ndim; j++) printf ("%d ",offsets[j]); printf ("\n");
*/
				for (j = 0; j < ndim; j++)
				{
					payload_size *= ldims [j];

					if ( (count[j] > gdims[j])
					     || (start[j] > gdims[j])
					     || (start[j] + count[j] > gdims[j]))
					{
						adios_error ( err_out_of_bound, "Error: Variable (id=%d) out of bound 1("
						              "the data in dimension %d to read is %llu elements from index %llu"
						              " but the actual data is [0,%llu])\n",
						              r->varid, j + 1, count[j], start[j], gdims[j] - 1);
						return 0;
					}

					/* check if there is any data in this pg and this dimension to read in */
					flag = (offsets[j] >= start[j]
					        && offsets[j] < start[j] + count[j])
						|| (offsets[j] < start[j]
						    && offsets[j] + ldims[j] > start[j] + count[j])
						|| (offsets[j] + ldims[j] > start[j]
						    && offsets[j] + ldims[j] <= start[j] + count[j]);
					idx_table [idx] = idx_table[idx] && flag;
				}

				if (!idx_table[idx])
				{
					continue;
				}

				/* determined how many (fastest changing) dimensions can we read in in one read */
				int hole_break;
				for (i = ndim - 1; i > -1; i--)
				{
					if (offsets[i] == start[i] && ldims[i] == count[i])
					{
						datasize *= ldims[i];
					}
					else
						break;
				}

				hole_break = i;
				slice_offset = 0;
				slice_size = 0;

				if (hole_break == -1)
				{
					/* The complete read happens to be exactly one pg, and the entire pg */
					/* This means we enter this only once, and npg=1 at the end */
					/* This is a rare case. FIXME: cannot eliminate this? */
					slice_size = payload_size;

					slice_offset = v->characteristics[start_idx + idx].payload_offset;

					memcpy ((char *)data, &xd->data[slice_offset], slice_size);
					if (fh->mfooter.change_endianness == adios_flag_yes)
					{
						log_error("endianness change not needed\n");
					}
				}
				else if (hole_break == 0)
				{
					/* The slowest changing dimensions should not be read completely but
					   we still need to read only one block */
					uint64_t isize;
					uint64_t size_in_dset = 0;
					uint64_t offset_in_dset = 0;
					uint64_t offset_in_var = 0;

					isize = offsets[0] + ldims[0];
					if (start[0] >= offsets[0])
					{
						// head is in
						if (start[0]<isize)
						{
							if (start[0] + count[0] > isize)
								size_in_dset = isize - start[0];
							else
								size_in_dset = count[0];
							offset_in_dset = start[0] - offsets[0];
							offset_in_var = 0;
						}
					}
					else
					{
						// middle is in
						if (isize < start[0] + count[0])
							size_in_dset = ldims[0];
						else
							// tail is in
							size_in_dset = count[0] + start[0] - offsets[0];
						offset_in_dset = 0;
						offset_in_var = offsets[0] - start[0];
					}

					slice_size = size_in_dset * datasize * size_of_type;
					write_offset = offset_in_var * datasize * size_of_type;

					slice_offset = v->characteristics[start_idx + idx].payload_offset
						+ offset_in_dset * datasize * size_of_type;
                    
					memcpy ((char *)data + write_offset, &xd->data[slice_offset], slice_size);

					if (fh->mfooter.change_endianness == adios_flag_yes)
					{
						log_error("endianness should not be changed\n");
					}

					//write_offset +=  slice_size;
				}
				else
				{
					uint64_t isize;
					uint64_t size_in_dset[10];
					uint64_t offset_in_dset[10];
					uint64_t offset_in_var[10];

					memset(size_in_dset, 0 , 10 * 8);
					memset(offset_in_dset, 0 , 10 * 8);
					memset(offset_in_var, 0 , 10 * 8);

					for (i = 0; i < ndim; i++)
					{
						isize = offsets[i] + ldims[i];
						if (start[i] >= offsets[i])
						{
							// head is in
							if (start[i]<isize)
							{
								if (start[i] + count[i] > isize)
									size_in_dset[i] = isize - start[i];
								else
									size_in_dset[i] = count[i];
								offset_in_dset[i] = start[i] - offsets[i];
								offset_in_var[i] = 0;
							}
							else
							{
							}
						}
						else
						{
							// middle is in
							if (isize < start[i] + count[i])
							{
								size_in_dset[i] = ldims[i];
							}
							else
							{
								// tail is in
								size_in_dset[i] = count[i] + start[i] - offsets[i];
							}
							offset_in_dset[i] = 0;
							offset_in_var[i] = offsets[i] - start[i];
						}
					}
					datasize = 1;
					var_stride = 1;

					for (i = ndim - 1; i >= hole_break; i--)
					{
						datasize *= size_in_dset[i];
						dset_stride *= ldims[i];
						var_stride *= count[i];
					}

					uint64_t start_in_payload = 0, end_in_payload = 0, s = 1;
					for (i = ndim - 1; i > -1; i--)
					{
						start_in_payload += s * offset_in_dset[i] * size_of_type;
						end_in_payload += s * (offset_in_dset[i] + size_in_dset[i] - 1) * size_of_type;
						s *= ldims[i];
					}

					slice_size = end_in_payload - start_in_payload + 1 * size_of_type;
					slice_offset =  v->characteristics[start_idx + idx].payload_offset
						+ start_in_payload;

					memcpy ((char *)data, &xd->data[slice_offset], slice_size);
						
					for (i = 0; i < ndim; i++)
					{
						offset_in_dset[i] = 0;
					}

					uint64_t var_offset = 0;
					uint64_t dset_offset = 0;

					for (i = 0; i < ndim; i++)
					{
						var_offset = offset_in_var[i] + var_offset * count[i];
						dset_offset = offset_in_dset[i] + dset_offset * ldims[i];
					}

					//FIXME find out what the heck is going on here
					copy_data (data
					           ,fh->b->buff + fh->b->offset
					           ,0
					           ,hole_break
					           ,size_in_dset
					           ,ldims
					           ,count
					           ,var_stride
					           ,dset_stride
					           ,var_offset
					           ,dset_offset
					           ,datasize
					           ,size_of_type
					           ,fh->mfooter.change_endianness
					           ,v->type
						);
				}
			}  // end for (idx ... loop over pgs

			free (idx_table);

			total_size += items_read * size_of_type;
			// shift target pointer for next read in
			data = (char *)data + (items_read * size_of_type);

		}
	} // end for t

	free (dims);

	chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
	assert (chunk);

	chunk->varid = r->varid;
	chunk->type = v->type;
	// NCSU ALACRITY-ADIOS - Added timestep information into varchunks
	chunk->from_steps = r->from_steps;
	chunk->nsteps = r->nsteps;
	chunk->sel = copy_selection (r->sel);
	chunk->data = r->data;
	return chunk;
}

/* This routine reads a write block. The 'index' value in the selection is
 * the block index within the context of the current step. Therefore, we
 * need to translate it to an absolute index.
 */
static ADIOS_VARCHUNK * read_var_wb (const ADIOS_FILE * fp, read_request * r)
{
	xpmem_read_file *xf = (xpmem_read_file *)fp->fh;	
	BP_PROC * p = GET_BP_PROC (fp);
	BP_FILE * fh = GET_BP_FILE (fp);
	xpmem_read_data *xd = xf->fp;

	struct adios_index_var_struct_v1 * v;
	int i = 0, j = 0, varid = 0, start_idx = 0, idx = 0;
	int ndim, has_subfile;
	uint64_t ldims[32], gdims[32], offsets[32];
	int size_of_type;
	uint64_t slice_offset, slice_size;
	void * data;
	ADIOS_VARCHUNK * chunk;
	MPI_Status status;
	const ADIOS_SELECTION_WRITEBLOCK_STRUCT *wb;// NCSU ALACRITY-ADIOS

	adios_errno = 0;

	has_subfile = has_subfiles (fh);
	data = r->data;
	varid = r->varid; //varid = map_req_varid (fp, r->varid); // NCSU ALACRITY-ADIOS: Bugfix: r->varid has already been mapped
	v = bp_find_var_byid (fh, varid);

	// NCSU ALACRITY-ADIOS: Add support for absolute PG index for efficiency
	//time = adios_step_to_time (fp, r->varid, r->from_steps);
	//idx = adios_wbidx_to_pgidx (fp, r);
	assert(r->sel->type == ADIOS_SELECTION_WRITEBLOCK);
	wb = &r->sel->u.block;

	idx = wb->is_absolute_index ? wb->index : adios_wbidx_to_pgidx (fp, r, i);
	assert (idx >= 0);

	ndim = v->characteristics [idx].dims.count;
	size_of_type = bp_get_type_size (v->type, v->characteristics [idx].value);

	if (ndim == 0)
	{
		r->datasize = size_of_type;
		slice_size = size_of_type;
		start_idx = 0; // OPS macros below need it

		if (v->type == adios_string)
		{
			size_of_type--;
		}

		slice_offset = v->characteristics[idx].payload_offset;

		//now copy the data from the buffer to
		//the user supplied data buffer
		//we copy from slice_offset and of size slice_size
		memcpy((char*)data, &xd->data[slice_offset], slice_size);

		if (fh->mfooter.change_endianness == adios_flag_yes)
		{
			log_error("should not have endianness issue in shared memory\n");
			// change_endianness ((char *)data,
			//                    size_of_type,
			//                    v->type
			// 	);
		}

		if (v->type == adios_string)
		{
			((char*)data)[size_of_type] = '\0';
		}

		//no need to advance since we are only reading 1 timestep
		// data = (char *) data + size_of_type;
	}
	else
	{
		// NCSU ALACRITY-ADIOS: Added sub-PG writeblock selection support
		// If this is a sub-PG selection, use nelements to compute slice_size
		// instead
		if (wb->is_sub_pg_selection) {
			// The start and end of the sub-PG selection must fall within the PG
			slice_size = wb->nelements * size_of_type;
		} else {
			// NCSU ALACRITY-ADIOS: This used to not be inside an else block
			// Else, do the old method of computing PG size from bounds
			slice_size = size_of_type;

			/* To get ldims for the chunk and then calculate payload size */
			bp_get_dimension_characteristics(&(v->characteristics[idx]),
			                                 ldims, gdims, offsets);

			for (j = 0; j < ndim; j++)
			{
				slice_size *= ldims [j];
			}
		}

		r->datasize = slice_size;
		/* Note: MPI_FILE_READ_OPS1 - for reading single BP file.
		 *       MPI_FILE_READ_OPS2 - for reading those with subfiles.
		 * Whenever to use OPS macro, start_idx and idx variable needs to be
		 * properly set.
		 */
		start_idx = 0;
		slice_offset = v->characteristics[idx].payload_offset;

		// NCSU ALACRITY-ADIOS: Added sub-PG writeblock selection support
		// If this is a sub-PG read, add the element_offset within the PG to the base offset in the file
		if (wb->is_sub_pg_selection) {
			slice_offset += wb->element_offset * size_of_type;
		}

		
		//now copy the data from the buffer to
		//the user supplied data buffer
		//we copy from slice_offset and of size slice_size
		memcpy((char*)data, &xd->data[slice_offset], slice_size);


		// NCSU ALACRITY-ADIOS: Reading directly to user buffer eliminates the need for this memcpy (profiling revealed it was hurting performance for transformed data)
		//memcpy ((char *)data, fh->b->buff + fh->b->offset, slice_size);
		if (fh->mfooter.change_endianness == adios_flag_yes)
		{
			log_error("should not have endianness issue in shared memory\n");
			// change_endianness ((char *)data, slice_size, v->type);
		}
	}

	chunk = (ADIOS_VARCHUNK *) malloc (sizeof (ADIOS_VARCHUNK));
	assert (chunk);

	chunk->varid = r->varid;
	chunk->type = v->type;
	// NCSU ALACRITY-ADIOS - Added timestep information into varchunks
	chunk->from_steps = r->from_steps;
	chunk->nsteps = r->nsteps;
	chunk->sel = copy_selection (r->sel);
	chunk->data = data;

	return chunk;
}


#endif
