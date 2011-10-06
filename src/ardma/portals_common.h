#ifndef __PORTAL_COMMON_H_
#define __PORTAL_COMMON_H_

#if defined(__LIBCATAMOUNT__)
#include <portals/portals3.h>
#include <stdint.h>
typedef uint32_t ptl_time_t;
typedef void (*ptl_eq_handler_t)(ptl_event_t *event);
#define PtlEventKindStr(x) __FUNCTION__
#define PtlErrorStr(x) ptl_err_str[x]
#elif defined(__CRAYXT_COMPUTE_LINUX_TARGET) || defined(XT3)
#include <portals/portals3.h>
#include <stdint.h>
#define PtlEventKindStr(x) __FUNCTION__
#define PtlErrorStr(x) ptl_err_str[x]
#else
#include <portals3.h>
#include <p3nal_utcp.h>
#include <p3api/debug.h>
#endif


#ifndef PTL_EQ_HANDLER_NONE
#define PTL_EQ_HANDLER_NONE NULL
#endif

#ifndef PTL_NO_ACK_REQ
#define PTL_NO_ACK_REQ PTL_NOACK_REQ
#endif


#ifndef is_end_event
#define is_end_event(event) \
   ((event.type == PTL_EVENT_GET_END)           || \
    (event.type == PTL_EVENT_PUT_END)           || \
    (event.type == PTL_EVENT_REPLY_END)         || \
    (event.type == PTL_EVENT_SEND_END)          || \
    (event.type == PTL_EVENT_GETPUT_END))
#endif

//connection match

static int CONN_MATCH=0x0000001;
static int DATA_MATCH=0x0000010;
#define MAXPATHLEN 255

static int LISTSIZE=1024*1024;




typedef struct endpoint_info
{
    ptl_process_id_t id;
    ptl_pt_index_t pt;

    //identify the destination array for requests
    ptl_match_bits_t match;
    ptl_hdr_data_t hdr;
	
    //ac index that we don't even use :P
    ptl_ac_index_t ac;	
} einfo;

typedef struct nic_info_
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
    ptl_match_bits_t match;

}ninfo;

typedef struct data_req
{
    uint32_t rank;

    ptl_size_t pg_size;
    ptl_match_bits_t pg_match;
    ptl_size_t pg_offset;
	
    ptl_size_t idx_size;
    ptl_match_bits_t idx_match;
    ptl_size_t idx_offset;	
    char path[MAXPATHLEN];
    uint32_t timestep;
}req_data;

typedef struct conn_req_
{
    einfo info;
    int rank;
    int size;
    int id;
}req_conn;

	

typedef struct msg_
{
	ptl_handle_md_t h;
	ptl_md_t md;
	void *buffer;
	int size;
	ptl_match_bits_t match;
}msg;

extern int  portal_init_common(ninfo *local_info);


static int QSIZE = 8192;
ptl_process_id_t pidany = {PTL_PID_ANY,PTL_NID_ANY};

#endif
