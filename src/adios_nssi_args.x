/* -------------------------------------------------------------------------- */
/**
 *   @file adios_nessie_args.x
 *
 *   @brief XDR argument structures for the adios proxy.
 *
 *   @author Todd Kordenbrock (thkorde\@sandia.gov).
 *
 */

/* Extra stuff to put at the beginning of the header file */
#ifdef RPC_HDR
/*%#include <sys/syslimits.h>  /* ADIOS_PATH_MAX, ... */
#endif

/* Extra stuff to put at the beginning of the C file */
#ifdef RPC_XDR
%#include <adios.h>
%#include "adios_nssi_args.h"
#endif

/**
 * Operation codes for the adios proxy.
 */
enum adios_opcode {
    ADIOS_NULL_OP = 10000,
    ADIOS_OPEN_OP,
    ADIOS_GROUP_SIZE_OP,
    ADIOS_READ_OP,
    ADIOS_WRITE_OP,
    ADIOS_END_ITER_OP,
    ADIOS_START_CALC_OP,
    ADIOS_STOP_CALC_OP,
    ADIOS_CLOSE_OP,
    ADIOS_FINALIZE_OP,

    ADIOS_READ_FOPEN_OP,
    ADIOS_READ_FCLOSE_OP,
    ADIOS_READ_GET_VARTYPE_SIZE_OP,
    ADIOS_READ_READ_VAR_OP
};


enum extra_errcodes {
    ADIOS_ENOTSUPP = -99
};

enum write_type {
    WRITE_DIRECT,
    WRITE_AGGREGATE_INDEPENDENT,
    WRITE_AGGREGATE_COLLECTIVE,
    WRITE_CACHING_INDEPENDENT,
    WRITE_CACHING_COLLECTIVE
};

enum adios_open_mode {
    ADIOS_MODE_READ,
    ADIOS_MODE_WRITE,
    ADIOS_MODE_APPEND,
    ADIOS_MODE_UPDATE

};

const ADIOS_PATH_MAX = 256;
const ADIOS_DIM_MAX = 16;

/* ********* ARGUMENTS FOR STUB FUNCTIONS ************* */

/**
 * Argument structure for adios_open
 */
struct adios_open_args {
    string          client_id<ADIOS_PATH_MAX>;
    string          gname<ADIOS_PATH_MAX>;
    string          fname<ADIOS_PATH_MAX>;
    adios_open_mode mode;
    int8_t          use_single_server;
};

/**
 * Structure for adios_open result
 */
struct adios_open_res {
    int64_t fd;
};

/**
 * Argument structure for adios_group_size
 */
struct adios_group_size_args {
    int64_t  fd;
    uint64_t data_size;
};

struct adios_var {
    uint8_t  is_anonymous;
    uint64_t vdata;
    uint64_t vdatasize;
    string   vpath<ADIOS_PATH_MAX>;
    string   vname<ADIOS_PATH_MAX>;
};

/**
 * Marshaled arguments for adios_read
 */
struct adios_read_args {
    int64_t  fd;
    string   vpath<ADIOS_PATH_MAX>;
    string   vname<ADIOS_PATH_MAX>;
    uint64_t max_read;
    uint16_t is_scalar;
    struct adios_var offsets<ADIOS_DIM_MAX>;
    struct adios_var ldims<ADIOS_DIM_MAX>;
    struct adios_var gdims<ADIOS_DIM_MAX>;
};

/**
 * Marshaled arguments for adios_read result
 */
struct adios_read_res {
    uint64_t bytes_read;
};

/**
 * Marshalled argument structure for adios_write
 */
struct adios_write_args {
    int64_t  fd;
    string   vpath<ADIOS_PATH_MAX>;
    string   vname<ADIOS_PATH_MAX>;
    uint64_t vsize;
    uint64_t atype;
    uint16_t is_scalar;
    int64_t  writer_rank;
    struct adios_var offsets<ADIOS_DIM_MAX>;
    struct adios_var ldims<ADIOS_DIM_MAX>;
    struct adios_var gdims<ADIOS_DIM_MAX>;
};

/**
 * Marshaled arguments for adios_write result
 */
struct adios_write_res {
    uint64_t bytes_written;
};

/**
 * Marshalled arguments for adios_end_iteration
 */
struct adios_end_iter_args {
    int64_t fd;
};

/**
 * Marshalled arguments for adios_start_calculation
 */
struct adios_start_calc_args {
    int64_t fd;
};

/**
 * Marshalled arguments for adios_stop_calculation
 */
struct adios_stop_calc_args {
    int64_t fd;
};

/**
 * Marshalled arguments for adios_close
 */
struct adios_close_args {
    string  fname<ADIOS_PATH_MAX>;
    int64_t fd;
};

/**
 * Marshalled arguments for adios_finalize
 */
struct adios_finalize_args {
    string client_id<ADIOS_PATH_MAX>;
};


/**
 * Argument structure for adios_read_fopen
 */
struct adios_read_fopen_args {
    string client_id<ADIOS_PATH_MAX>;
    string gname<ADIOS_PATH_MAX>;
    string fname<ADIOS_PATH_MAX>;
    uint64_t requested_timestep;
};

/**
 * Structure for adios_read_fopen result
 */
struct adios_read_fopen_res {
    int64_t fd;
};

/**
 * Marshaled arguments for adios_read
 */
struct adios_read_read_var_args {
    int64_t  fd;
    string   client_id<ADIOS_PATH_MAX>;
    uint64_t open_timestep;
    string   vpath<ADIOS_PATH_MAX>;
    string   vname<ADIOS_PATH_MAX>;
    uint64_t max_read;
    uint16_t is_scalar;

    uint64_t offsets[3];
    uint64_t counts[3];
};

/**
 * Marshaled arguments for adios_read result
 */
struct adios_read_read_var_res {
    uint64_t bytes_read;
};

/**
 * Marshaled arguments for adios_read
 */
struct adios_read_get_vartype_size_args {
    int64_t  fd;
    string   client_id<ADIOS_PATH_MAX>;
    uint64_t open_timestep;
    string   vpath<ADIOS_PATH_MAX>;
    string   vname<ADIOS_PATH_MAX>;
};

/**
 * Marshaled arguments for adios_read result
 */
struct adios_read_get_vartype_size_res {
    uint64_t vartype_size;
};

/**
 * Marshalled arguments for adios_fclose
 */
struct adios_read_fclose_args {
    int64_t  fd;
    uint64_t open_timestep;
    string   client_id<ADIOS_PATH_MAX>;
    string   fname<ADIOS_PATH_MAX>;
};
