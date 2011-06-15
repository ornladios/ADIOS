#ifndef __PORTAL_COMMON_H_
#define __PORTAL_COMMON_H_


//connection match

static int CONN_MATCH=0x0101010101;

static int LISTSIZE=1024*1024;


typedef struct portal_info
{
	ptl_process_id_t id;
	ptl_pt_index_t pt;

	//identify the destination array for requests
	ptl_match_bits_t match;
	ptl_hdr_data_t hdr;
	
	//ac index that we don't even use :P
	ptl_ac_index_t ac;	
} ptlinfo;

typedef struct info_
{
	ptl_interface_t iface;
	ptl_ni_limits_t nilimits;
	ptl_handle_ni_t nihandle;
	ptl_uid_t uid;
	ptl_process_id_t pid;
	ptl_handle_eq_t eqh; //eq for connection request
	ptl_handle_me_t meh; //match entry for connection request
	ptl_pt_index_t index;
	ptl_ac_index_t ac;	

}info;

typedef struct data_req
{
	uint32_t rank;

	ptl_size_t pg_size;
	ptl_match_bits_t pg_match;
	ptl_size_t pg_offset;
	
	ptl_size_t idx_size;
	ptl_match_bits_t idx_match;
	ptl_size_t idx_offset;
	
	char patch[MAXPATHLEN];
	uint32_t timestep;
}portal_req;

typedef struct conn_req_
{
	ptlinfo info;
	int rank;
	int size;
	int id;
}conn, conn_req;

	

typedef struct msg_
{
	ptl_handle_md_t h;
	ptl_md_t md;
	void *buffer;
	int size;
	ptl_match_bits_t match;
}msg;


static void portal_init_common(info *local_info)
{
	
	//initialize portals device

    int retval = 0;
    int max_interfaces;

#if  defined(__LIBCATAMOUNT__)
    ptl_interface_t ptl_iface = CRAY_QK_NAL;
#elif defined(__CRAYXT_COMPUTE_LINUX_TARGET) || defined(XT3)
    ptl_interface_t ptl_iface = CRAY_UK_SSNAL;
#else
    ptl_interface_t ptl_iface = PTL_NALTYPE_UTCP;
#endif
    ptl_ni_limits_t nilimits;
    ptl_handle_ni_t nihandle;
    ptl_uid_t uid;
    ptl_process_id_t pid;
    ptl_handle_eq_t eq_handle;
    ptl_handle_me_t me_handle;
    ptl_pt_index_t  ptl_index;
    ptl_handle_eq_t eq_handle_data;

    retval = PtlInit(&max_interface);
	if(retval != PTL_OK)
	{
		log_error("PORTALS: PtlInit \t %s\n", PtlErrorStr(retval));
		return retval;
	}

	log_debug("max interface = %d\n", max_interface);
	
	
	retval = PtlNIInit(ptl_iface, PTL_PID_ANY, NULL, &nilimits, &nihandle);
	if(retval != PTL_OK)
	{
		log_error("PORTALS PtlNIInit \t%s\n", PtlErrorStr(retval));
		return retval;
		
	}
	
	log_debug("%d %d %d %d %d %d %d %d\n",
         nilimits.max_mes, nilimits.max_mds,
         nilimits.max_eqs, nilimits.max_ac_index,
         nilimits.max_pt_index, nilimits.max_md_iovecs,
         nilimits.max_me_list, nilimits.max_getput_md);
    

    retval = PtlGetUid(nihandle, &uid);
    if (retval != PTL_OK)
    {
        log_error("PORTALS: PtlGetUid \t%s\n", PtlErrorStr(retval));
        return retval;
    }

    retval = PtlGetId(nihandle, &pid);
    if (retval != PTL_OK)
    {
        log_error("PORTALS: PtlGetId \t%s\n", PtlErrorStr(retval));
        return retval;
    }

	retval = PtlEQAlloc(nihandle, QSIZE, PTL_EQ_HANDLER_NONE, &local_info->eqh);
	if(retval != PTL_OK)
	{
		log_error("PORTALS: PtlEQAlloc\t\%s\n", PtlErrorStr(retval));
		return retval;
	}
	
	retval = PtlMEAttachAny(nihandle, &ptl_index, pidany, CONN_MATCH, 0, PTL_RETAIN, &local_info->meh);
	if(retval != PTL_OK)
	{
		log_error("PORTALS: PtlMEAttachAny\t\%s\n", PtlErrorStr(retval));
	 	return retval;
	}
	

	log_info("Allocating the memory and MD for connection requests\n");
	
	

	local_info->iface = ptl_iface;
	local_info->nilimits = nilimits;
	local_info->nihandle = nihandle;
	local_info->uid = uid;
	local_info->pid = pid;
	local_info->index = ptl_index;
	
	
	log_debug("iface = %d nihandle = %p uid = %d pid.pid = %d pid.nid = %d index = %d\n",
			  local_info->iface, local_info->nihandle, local_info->uid,
			  local_info->pid.pid, local_info->pid.nid, local_info->index);

}


#endif
