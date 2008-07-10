#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* #define DEBUG */

#include "mpi-dummy.h"

#ifdef DEBUG
#define LOG(x) fprintf(stderr,"MPIDUMMY: %s\n",x);
#else
#define LOG(x) {};
#endif

/*
 * more advanced later versions -- actually have a message queue? 
 */

static short int sizes[MPI_MAXTYPES];
static short int freetype[MPI_MAXTYPES];
static int ntypes = 0;

int firstfree() {
	int i=0;
	while ((i<MPI_MAXTYPES) && (!freetype[i])) i++;

	if (i >= MPI_MAXTYPES) return -1;
	return i;
}

void printdata(char *s, void *data, MPI_DATATYPE dt, int count) {

	int i;


	fprintf(stderr,"%s =",s);

	switch (dt) {
		case MPI_CHARACTER:
			fprintf(stderr,"char: '");
			for (i=0; i<count; i++) 
				fprintf(stderr,"%c",((char *)data)[i]);
			fprintf(stderr,"'\n");
			break;
		case MPI_REAL:
			fprintf(stderr,"float: ");
			for (i=0; i<count; i++) 
				fprintf(stderr,"%f ",((float *)data)[i]);
			fprintf(stderr,"\n");
			break;
		case MPI_DOUBLE_PRECISION:
			fprintf(stderr,"dbl  : ");
			for (i=0; i<count; i++) 
				fprintf(stderr,"%lg ",((double *)data)[i]);
			fprintf(stderr,"\n");
			break;
		case MPI_2DOUBLE_PRECISION:
			fprintf(stderr,"2xdbl: ");
			for (i=0; i<count; i++) 
				fprintf(stderr,"(%lg,%lg) ",((double *)data)[2*i],
				                    ((double *)data)[2*i+1]);
			fprintf(stderr,"\n");
			break;
		case MPI_LOGICAL:
			fprintf(stderr,"bool : ");
			for (i=0; i<count; i++) 
				fprintf(stderr,"%s ",
                                 (((int *)data)[i] == 0 ? "FALSE" : "TRUE" ));
			fprintf(stderr,"\n");
			break;
		case MPI_INTEGER:
			fprintf(stderr,"int : ");
			for (i=0; i<count; i++) 
				fprintf(stderr,"%d ",((int *)data)[i]);
			fprintf(stderr,"\n");
			break;
                default:
			fprintf(stderr,"Can't print data of type %d\n",dt);
	}

	return;
}

void mpi_finalize_(int *ierr) {
        LOG("mpi_finalize_")
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_finalize(void) {
        LOG("mpi_finalize")
	return MPI_SUCCESS;
}

int mpi_abort(MPI_COMM c, int error ) {
        LOG("mpi_abort")
	fprintf(stderr, "MPI_Abort(C): error %d\n", error);
	mpi_finalize();
	exit(-1);
	return MPI_SUCCESS;
}

void mpi_abort_(MPI_COMM *c, int *error, int *ierr ) {
        LOG("mpi_abort_")
	fprintf(stderr, "MPI_Abort(F): error %d\n", *error);
        *ierr = MPI_SUCCESS;
	mpi_finalize();
	exit(-1);
        return;
}

void check_comm(MPI_COMM c) {
	if (c != MPI_COMM_WORLD) {
		fprintf(stderr,"Only supporting COMM_WORLD now!\n");
		fprintf(stderr,"  You tried to use communicator %d\n",c);
		mpi_abort(MPI_COMM_WORLD,1);
	}
}

void mpi_init_(int *ierr) {
	int i;
	*ierr = MPI_SUCCESS;

        LOG("mpi_init_")
        sizes[MPI_INTEGER] = sizeof(int); 
        sizes[MPI_REAL] = sizeof(float);
        sizes[MPI_LOGICAL] = 4*sizeof(char);
        sizes[MPI_DOUBLE_PRECISION] = 1*sizeof(double);
        sizes[MPI_2DOUBLE_PRECISION] = 2*sizeof(double);
        sizes[MPI_CHARACTER] = sizeof(char);
        ntypes = 6;

	for (i=0;i<ntypes;i++) freetype[i] = 0;
	for (i=ntypes;i<MPI_MAXTYPES;i++) freetype[i] = 1;

	return;
}

int mpi_init(int *argc, char **argv) {
	int foo;
        LOG("mpi_init")
	mpi_init_(&foo);
	return foo;
}

double mpi_wtime() {
        LOG("mpi_wtime")
	return (double)( ((double)clock())/CLOCKS_PER_SEC );
}

double mpi_wtime_() {
        LOG("mpi_wtime_")
	return (double)( ((double)clock())/CLOCKS_PER_SEC );
}

int mpi_barrier(MPI_COMM c) {
        LOG("mpi_barrier")
        check_comm(c);
	;
	return MPI_SUCCESS;
}

void mpi_barrier_(MPI_COMM *c, int *ierr) {
        LOG("mpi_barrier_")
        check_comm(*c);
	;
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_bcast(void *b, int count, MPI_DATATYPE dt, int root, MPI_COMM c) {

        LOG("mpi_bcast")
        check_comm(c);
        if (root != 0) {
		fprintf(stderr,"mpi_bcast(C): Invalid root %d\n",root);
                return MPI_FAILURE;
        }
	return MPI_SUCCESS;
}

void mpi_bcast_(void *b, int *count, MPI_DATATYPE *dt, int *root, MPI_COMM *c,
                int *ierr) {

        LOG("mpi_bcast_")
        check_comm(*c);
        if (*root != 0) {
		fprintf(stderr,"mpi_bcast(C): Invalid root %d\n",*root);
                *ierr = MPI_FAILURE;
	 	return;
        }
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_comm_rank(MPI_COMM c, int *rank) {
        LOG("mpi_comm_rank")
        check_comm(c);
	*rank = 0;
	return MPI_SUCCESS;
}

void mpi_comm_rank_(MPI_COMM *c, int *rank, int *ierr) {
        LOG("mpi_comm_rank_")
        check_comm(*c);
	*rank = 0;
	*ierr = MPI_SUCCESS;
	return;
}


int mpi_comm_size(MPI_COMM c, int *size) {
        LOG("mpi_comm_size")
        check_comm(c);
	*size = 1;
	return MPI_SUCCESS;
}

void mpi_comm_size_(MPI_COMM *c, int *size, int *ierr) {
        LOG("mpi_comm_size_")
        check_comm(*c);
	*size = 1;
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_allgather(void *sb, int sc, MPI_DATATYPE st,
                  void *rb, int rc, MPI_DATATYPE rt, MPI_COMM c) {

        LOG("mpi_allgather")
        check_comm(c);
	if (st > ntypes) {
		fprintf(stderr,"Error: allgatherc invalid stype %d!\n",st);
		return MPI_FAILURE;
	}
	if (rt > ntypes) {
		fprintf(stderr,"Error: allgatherc invalid rtype %d!\n",rt);
		return MPI_FAILURE;
	}

        fprintf(stderr,"mpi_allgather(C) called: %d x %d\n",sc,st);
        memcpy(rb, sb, sc*sizes[st] );
	
	return MPI_SUCCESS;
}

void mpi_allgather_(void *sb, int *sc, MPI_DATATYPE *st,
                    void *rb, int *rc, MPI_DATATYPE *rt, MPI_COMM *c, int *ierr)
{
        LOG("mpi_allgather_")
        check_comm(*c);
	if (*st > ntypes) {
		fprintf(stderr,"Error: allgatherf invalid stype %d!\n",*st);
		*ierr = MPI_FAILURE; return;
	}
	if (*rt > ntypes) {
		fprintf(stderr,"Error: allgatherf invalid rtype %d!\n",*rt);
		*ierr = MPI_FAILURE; return;
	}
	if (*st != *rt) {
		fprintf(stderr,"Error: allgatherf type Mistmatch %d != %d!\n",*st,*rt);
		*ierr =  MPI_FAILURE; return;
	}

#ifdef DEBUG
        fprintf(stderr,"mpi_allgather(F) called: %d x %d\n",*sc,*st);
        printdata("send: ", sb, *st, *sc); 
#endif
        memcpy(rb, sb, (*sc)*sizes[*st] );
#ifdef DEBUG
        printdata("recv: ", rb, *rt, *rc); 
#endif	
	*ierr = MPI_SUCCESS;
	return;
}


int mpi_gather(void *sb, int sc, MPI_DATATYPE st,
                  void *rb, int rc, MPI_DATATYPE rt, int root, MPI_COMM c) {

        LOG("mpi_gather")
        check_comm(c);
	if (st > ntypes) {
		fprintf(stderr,"Error: gatherc invalid stype %d!\n",st);
		return MPI_FAILURE;
	}
	if (rt > ntypes) {
		fprintf(stderr,"Error: gatherc invalid rtype %d!\n",rt);
		return MPI_FAILURE;
	}

        fprintf(stderr,"mpi_gather(C) called: %d x %d\n",sc,st);
        memcpy(rb, sb, sc*sizes[st] );
	
	return MPI_SUCCESS;
}

void mpi_gather_(void *sb, int *sc, MPI_DATATYPE *st,
                    void *rb, int *rc, MPI_DATATYPE *rt, int root, MPI_COMM *c, 
                    int *ierr)
{
        LOG("mpi_gather_")
        check_comm(*c);
	if (*st > ntypes) {
		fprintf(stderr,"Error: gatherf invalid stype %d!\n",*st);
		*ierr = MPI_FAILURE; return;
	}
	if (*rt > ntypes) {
		fprintf(stderr,"Error: gatherf invalid rtype %d!\n",*rt);
		*ierr = MPI_FAILURE; return;
	}
	if (*st != *rt) {
		fprintf(stderr,"Error: gatherf type Mistmatch %d != %d!\n",*st,*rt);
		*ierr =  MPI_FAILURE; return;
	}

#ifdef DEBUG
        fprintf(stderr,"mpi_gather(F) called: %d x %d\n",*sc,*st);
        printdata("send: ", sb, *st, *sc); 
#endif
        memcpy(rb, sb, (*sc)*sizes[*st] );
#ifdef DEBUG
        printdata("recv: ", rb, *rt, *rc); 
#endif	
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_waitall(int count, MPI_REQUEST *r, MPI_STATUS *s) {
        LOG("mpi_waitall")
	fprintf(stderr, "If you're waiting for messages in a 1-proc system,\n");
	fprintf(stderr, "you may be waiting a while.\n");
	fprintf(stderr, "waitallC\n");

	exit(-1);

	return MPI_FAILURE;
}

void mpi_waitall_(int *count, MPI_REQUEST *r, MPI_STATUS *s, int *ierr) {
        LOG("mpi_waitall_")
	fprintf(stderr, "If you're waiting for messages in a 1-proc system,\n");
	fprintf(stderr, "you may be waiting a while.\n");
	fprintf(stderr, "waitallF\n");

	exit(-1);

	*ierr = MPI_FAILURE;
	return;
}

int mpi_waitany(int count, MPI_REQUEST *r, int *index, MPI_STATUS *s) {
        LOG("mpi_waitany")
	fprintf(stderr, "If you're waiting for messages in a 1-proc system,\n");
	fprintf(stderr, "you may be waiting a while.\n");
	fprintf(stderr, "waitanyC\n");

	exit(-1);

	return MPI_FAILURE;
}

void mpi_waitany_(int *count, MPI_REQUEST *r, int *index, MPI_STATUS *s, 
        int *ierr) {
        LOG("mpi_waitany_")
	fprintf(stderr, "If you're waiting for messages in a 1-proc system,\n");
	fprintf(stderr, "you may be waiting a while.\n");
	fprintf(stderr, "waitanyF\n");

	exit(-1);

	*ierr = MPI_FAILURE;
	return;
}

int mpi_allreduce(void *sb, void *rb, int cout, MPI_DATATYPE dt,
                  MPI_OP op, MPI_COMM c) {
	int i;

        LOG("mpi_allreduce")
        check_comm(c);
	if (dt > ntypes) {
		fprintf(stderr,"Error: allreduceC invalid type %d!\n",dt);
		return MPI_FAILURE;
	}

#ifdef DEBUG
        fprintf(stderr,"mpi_allreduce(C) called: %d x %d (%d)\n",cout,dt,op);
#endif
        memcpy(rb, sb, cout*sizes[dt]  );

        if (dt == MPI_2DOUBLE_PRECISION) {
		for (i=0; i<cout;i++) {
			((double *)rb)[2*i+1] = (double)0;
		}
	}
	
	return MPI_SUCCESS;
}

void mpi_allreduce_(void *sb, void *rb, int *cout, MPI_DATATYPE *dt,
                  MPI_OP *op, MPI_COMM *c, int *ierr) {
	int i;

        LOG("mpi_allreduce_")
        check_comm(*c);
	if (*dt > ntypes) {
		fprintf(stderr,"Error: allreduceF invalid type %d!\n",*dt);
		*ierr = MPI_FAILURE;
		return;
	}

#ifdef DEBUG
        fprintf(stderr,"mpi_allreduce(F) called: %d x %d (%d)\n",*cout,*dt,*op);
        printdata("send: ", sb, *dt, *cout); 
#endif
        memcpy(rb, sb, (*cout)*(sizes[(*dt)])  );
/*        if (*dt == MPI_2DOUBLE_PRECISION) {
		for (i=0; i<(*cout);i++) {
			((double *)rb)[2*i+1] = (double)0;
		}
	} */
#ifdef DEBUG
        fprintf(stderr,"Copied %d bytes\n",(*cout)*sizes[(*dt)]);
        printdata("send: ", sb, *dt, *cout); 
        printdata("recv: ", rb, *dt, *cout); 
#endif	
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_iprobe(int source, int tag, MPI_COMM c, int *flag, MPI_STATUS *s) {
	*flag = 0;

        LOG("mpi_iprobe")
        check_comm(c);
        fprintf(stderr,"mpi_iprobe(C) called\n");
	return MPI_SUCCESS;
}

void mpi_iprobe_(int *source, int *tag, MPI_COMM *c, int *flag, MPI_STATUS *s, 
                 int *ierr) {
  
        LOG("mpi_iprobe_")
        check_comm(*c);
        fprintf(stderr,"mpi_iprobe(F) called\n");
        *flag = 0;
	*ierr = MPI_SUCCESS;
	return;
}

int mpi_reduce(void *sb, void *rb, int cout, MPI_DATATYPE dt,
                  MPI_OP op, int root, MPI_COMM c) {
        LOG("mpi_reduce")
        check_comm(c);
        fprintf(stderr,"mpi_reduce(C) called (op = %d) -- ",op);
        if (root != 0) {
           fprintf(stderr,"mpi_reduce(C): invalid root %d\n",root);
           return MPI_FAILURE;
        }
	return mpi_allreduce(sb,rb,cout,dt,op,c);
}

void mpi_reduce_(void *sb, void *rb, int *cout, MPI_DATATYPE *dt,
                  MPI_OP *op, int *root, MPI_COMM *c, int *ierr) {
        LOG("mpi_reduce_")
        check_comm(*c);
#ifdef DEBUG
        fprintf(stderr,"mpi_reduce(F) called (op = %d) -- ",*op);
#endif
        if (*root != 0) {
           fprintf(stderr,"mpi_reduce(C): invalid root %d\n",*root);
           *ierr = MPI_FAILURE;
           return;
        }
	mpi_allreduce_(sb,rb,cout,dt,op,c,ierr);
	return;
}

int mpi_scan(void *sb, void *rb, int cout, MPI_DATATYPE dt,
                  MPI_OP op, MPI_COMM c) {
        LOG("mpi_scan")
        check_comm(c);
#ifdef DEBUG
        fprintf(stderr,"mpi_scan(C) called (op = %d) -- ",op);
#endif
	return mpi_allreduce(sb,rb,cout,dt,op,c);
}

void mpi_scan_(void *sb, void *rb, int *cout, MPI_DATATYPE *dt,
                  MPI_OP *op, MPI_COMM *c, int *ierr) {
        LOG("mpi_scan_")
        check_comm(*c);
#ifdef DEBUG
        fprintf(stderr,"mpi_scan(F) called (op = %d) -- ",*op);
#endif
	mpi_allreduce_(sb,rb,cout,dt,op,c,ierr);
	return;
}

int mpi_ssend(void *b, int count, MPI_DATATYPE dt, int source, int tag, 
              MPI_COMM c){
        LOG("mpi_ssend")
        check_comm(c);
	fprintf(stderr,"Other processor %d not listening!\n",source);
        return MPI_FAILURE;     
}

void mpi_ssend_(void *b, int *count, MPI_DATATYPE *dt, int *source, int *tag, 
                MPI_COMM *c, int *ierr){
        LOG("mpi_ssend_")
        check_comm(*c);
	fprintf(stderr,"Other processor %d not listening!\n",*source);
        *ierr = MPI_FAILURE;     
	return;
}

int mpi_send(void *b, int count, MPI_DATATYPE dt, int source, int tag, 
             MPI_COMM c){
        LOG("mpi_send")
        check_comm(c);
	fprintf(stderr,"Other processor %d not listening!\n",source);
        return MPI_FAILURE;     
}

void mpi_send_(void *b, int *count, MPI_DATATYPE *dt, int *source, int *tag, 
                MPI_COMM *c, int *ierr){
        LOG("mpi_send_")
        check_comm(*c);
	fprintf(stderr,"Other processor %d not listening!\n",*source);
        *ierr = MPI_FAILURE;     
	return;
}

int mpi_recv(void *b, int count, MPI_DATATYPE dt, int source, int tag, 
             MPI_COMM c){
        LOG("mpi_recv")
        check_comm(c);
	fprintf(stderr,"Other processor %d not sending!\n",source);
        return MPI_FAILURE;     
}

void mpi_recv_(void *b, int *count, MPI_DATATYPE *dt, int *source, int *tag, 
                MPI_COMM *c, int *ierr){
        LOG("mpi_recv_")
        check_comm(*c);
	fprintf(stderr,"Other processor %d not sending!\n",*source);
        *ierr = MPI_FAILURE;     
	return;
}

int mpi_isend(void *b, int count, MPI_DATATYPE dt, int source, int tag, 
              MPI_COMM c){
        LOG("mpi_isend");
        check_comm(c);
	fprintf(stderr,"Other processor %d not listening!\n",source);
        return MPI_FAILURE;     
}

void mpi_isend_(void *b, int *count, MPI_DATATYPE *dt, int *source, int *tag, 
                MPI_COMM *c, int *ierr){
        LOG("mpi_isend_");
        check_comm(*c);
	fprintf(stderr,"Other processor %d not listening!\n",*source);
        *ierr = MPI_FAILURE;     
	return;
}

int mpi_irecv(void *b, int count, MPI_DATATYPE dt, int source, int tag, 
             MPI_COMM c){
        check_comm(c);

        LOG("mpi_irecv");
	fprintf(stderr,"Other processor %d not sending!\n",source);
        return MPI_FAILURE;     
}

void mpi_irecv_(void *b, int *count, MPI_DATATYPE *dt, int *source, int *tag, 
                MPI_COMM *c, int *ierr){
        check_comm(*c);

        LOG("mpi_irecv_");
	fprintf(stderr,"Other processor %d not sending!\n",*source);
        *ierr = MPI_FAILURE;     
	return;
}

int mpi_sendrecv(void *sb, int sc, MPI_DATATYPE st, int dest, int stag,
                 void *rb, int rc, MPI_DATATYPE rt, int src,  int rtag,
                 MPI_COMM c, MPI_STATUS *s) {

        check_comm(c);

        LOG("mpi_sendrecv");
        if ((dest == MPI_PROC_NULL) && (src == MPI_PROC_NULL)) {
            return  MPI_SUCCESS;
        }

	if (dest != 0 || src != 0) {
		fprintf(stderr, "sendrecvC: %d->%d?\n",src,dest);
		return MPI_FAILURE;
	}

        if (st > MPI_CHARACTER) {
		fprintf(stderr,"sendrecvC: sending derived type.\n");
        }
        if (rt > MPI_CHARACTER) {
		fprintf(stderr,"sendrecvC: recieving derived type.\n");
        }

#ifdef DEBUG
        fprintf(stderr,"mpi_sendrecv(C) called: %d x %d\n",sc,st);
#endif
        memcpy(rb, sb, sc*sizes[st]  );
	
	return MPI_SUCCESS;
}

void mpi_sendrecv_(void *sb, int *sc, MPI_DATATYPE *st, int *dest, int *stag,
                   void *rb, int *rc, MPI_DATATYPE *rt, int *src,  int *rtag,
                   MPI_COMM *c, MPI_STATUS *s, int *ierr) {

        LOG("mpi_sendrecv_");
        check_comm(*c);

        if ((*dest == MPI_PROC_NULL) && (*src == MPI_PROC_NULL)) {
            *ierr = MPI_SUCCESS;
            return; 
        }

        if ((*dest != 0) && (*src != 0)) {
            	fprintf(stderr,"sendrecvF: %d -> %d?\n",*src,*dest);
                *ierr = MPI_FAILURE;
                return;
        }

        if (*st > MPI_CHARACTER) {
		fprintf(stderr,"sendrecvF: sending derived type.\n");
        }
        if (*rt > MPI_CHARACTER) {
		fprintf(stderr,"sendrecvF: recieving derived type.\n");
        }

#ifdef DEBUG
        fprintf(stderr,"mpi_sendrecv(F) called: %d x %d\n",*sc,*st);
#endif
	*ierr = mpi_sendrecv(sb,*sc,*st,*dest,*stag,rb,*rc,*rt,*src,*rtag,*c,s);

	return;
}

void mpi_type_commit_(MPI_DATATYPE *dt, int *ierr) {
        LOG("mpi_type_commit_");
	*ierr = MPI_SUCCESS;	
	return;	
}

void mpi_type_free_(MPI_DATATYPE *dt, int *ierr) {
        int i;
        LOG("mpi_type_free_");
        if ((*dt < 0) || (*dt > MPI_MAXTYPES)) {
		fprintf(stderr,"MPI_Type_FreeF Invalid type %d!\n",*dt);
		*ierr = MPI_FAILURE;
		return;
	}
        if (freetype[*dt]) {
		fprintf(stderr,"MPI_Type_FreeF Invalid type %d!\n",*dt);
		*ierr = MPI_FAILURE;
		return;
	}
        if (*dt < 6) {
		fprintf(stderr,"MPI_Type_FreeF Can't Free Static Type %d!\n",*dt);
		*ierr = MPI_FAILURE;
		return;
	}

	freetype[*dt] = 1;
	ntypes--;
	*ierr = MPI_SUCCESS;	
	return;	
}

void mpi_type_contiguous_(int *count, MPI_DATATYPE *old, MPI_DATATYPE *new, 
                          int *ierr) {
	int i;

        LOG("mpi_type_contiguous_");
	if (*old >= ntypes) {
		fprintf(stderr,"MPI_Type_ContiguousF Invalid type %d!\n",*old);
		*ierr = MPI_FAILURE;
		return;
	}

        if ((*count < 0) || (*count > 999999)) {
		fprintf(stderr,"MPI_Type_ContiguousF Invalid count %d!\n", 
                        *count);
		*ierr = MPI_FAILURE;
		return;
	}

        i = 0;
        i = firstfree();

        if (i >= MPI_MAXTYPES || i < 0) {
		fprintf(stderr,"MPI_Type_ContiguousF FAILED: too many types\n");
                *ierr = MPI_FAILURE;
        } else  {
		*new = i;
                freetype[i] = 0;
		sizes[i] = (*count)*sizes[*old];
		ntypes++;
		*ierr = MPI_SUCCESS;	
	}
	return;	
}

void mpi_type_hvector_(int *count, int *blocklen, int *stride,
                       MPI_DATATYPE *old, MPI_DATATYPE *new, 
                       int *ierr) {
 	int i = 0;

        LOG("mpi_type_hvector_");
	if (*old >= ntypes) {
		fprintf(stderr,"MPI_Type_HVectorF Invalid type %d!\n",*old);
		*ierr = MPI_FAILURE;
		return;
	}

        i = firstfree();

        if (i >= MPI_MAXTYPES || i < 0) {
		fprintf(stderr,"MPI_Type_HVectorF FAILED: too many types\n");
                *ierr = MPI_FAILURE;
        } else  {
		*new = i;
                freetype[i] = 0;
		sizes[i] = (*count)*sizes[*old];
		ntypes++;
		*ierr = MPI_SUCCESS;	
	}
	return;	
}

void mpi_type_vector_(int *count, int *blocklen, int *stride,
                       MPI_DATATYPE *old, MPI_DATATYPE *new, 
                       int *ierr) {
	int i;

        LOG("mpi_type_vector_");
	if (*old >= ntypes) {
		fprintf(stderr,"MPI_Type_VectorF Invalid type %d!\n",*old);
		*ierr = MPI_FAILURE;
		return;
	}


        i = firstfree();

        if (i >= MPI_MAXTYPES || i < 0) {
		fprintf(stderr,"MPI_Type_VectorF FAILED: too many types\n");
                *ierr = MPI_FAILURE;
        } else  {
		*new = i;
                freetype[i] = 0;
		sizes[i] = (*count)*sizes[*old];
		ntypes++;
		*ierr = MPI_SUCCESS;	
	}
	return;	
}
