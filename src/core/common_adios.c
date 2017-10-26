/*
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>  // gettimeofday
#include <unistd.h>

// xml parser
#include <mxml.h>

#include "core/adios_bp_v1.h"
#include "core/adios_internals.h"
#include "core/adios_internals_mxml.h"
#include "core/adios_logger.h"
#include "core/adios_timing.h"
#include "core/adios_transport_hooks.h"
#include "core/buffer.h"
#include "core/common_adios.h"
#include "core/qhashtbl.h"
#include "core/types.h"
#include "public/adios_error.h"

// NCSU ALACRITY-ADIOS
#include "core/transforms/adios_transforms_common.h"
#include "core/transforms/adios_transforms_read.h"
#include "core/transforms/adios_transforms_write.h"

#ifdef WITH_NCSU_TIMER
#include "timer.h"
#endif

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#include "adiost_callback_internal.h"

extern struct adios_transport_struct *adios_transports;
extern int adios_errno;

///////////////////////////////////////////////////////////////////////////////
int common_adios_init(const char *config, MPI_Comm comm) {
  // parse the config file
  if (comm == MPI_COMM_NULL) comm = MPI_COMM_SELF;
  adios_errno = err_no_error;
  adiost_pre_init();
  adios_parse_config(config, comm);
  adiost_post_init();
  ADIOST_CALLBACK(adiost_event_init, config, comm);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// all XML file pieces will be provided by another series of calls
int common_adios_init_noxml(MPI_Comm comm) {
  if (comm == MPI_COMM_NULL) comm = MPI_COMM_SELF;
  adios_errno = err_no_error;
  adiost_pre_init();
  adios_local_config(comm);
  adiost_post_init();
  ADIOST_CALLBACK(adiost_event_init_noxml, comm);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_finalize(int mype) {
  ADIOST_CALLBACK_ENTER(adiost_event_finalize, mype);
  struct adios_method_list_struct *m;
  struct adios_group_list_struct *g;

  // Time Aggregation: there may be time steps left in the time aggregated
  // buffer, which
  // needs to be dumped out before finalize
  for (g = adios_get_groups(); g; g = g->next) {
    if (TimeAggregationInProgress(g->group)) {
      // log_debug("time buffering enabled and data left in finalize... close
      // file\n");
      /* At this point the buffer of the last step has been closed and
       * the index has been built but it was not yet passed on to the
       * method's close function.
       * do_ts_finalize signals this special case to common_adios_close()
       */
      SetTimeAggregationFlush(g->group, 1);
      // => TimeAggregationIsFlushing => TimeAggregationLastStep
      common_adios_close(g->group->ts_fd);  // close file
      SetTimeAggregation(g->group, 0);      // turn off ts buffering
    }
  }

  adios_errno = err_no_error;
  for (m = adios_get_methods(); m; m = m->next) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_finalize_fn) {
      adios_transports[m->method->m].adios_finalize_fn(mype, m->method);
    }
  }

  adios_cleanup();

#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_finalize();
#endif

  ADIOST_CALLBACK_EXIT(adiost_event_finalize, mype);
  adiost_finalize();
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_allocate_buffer(
    enum ADIOS_BUFFER_ALLOC_WHEN adios_buffer_alloc_when,
    uint64_t buffer_size) {
  adios_errno = err_no_error;
  log_warn(
      "adios_allocate_buffer is not supported anymore. "
      "Use adios_set_max_buffer_size(size_in_MB) to set the maximum allowed "
      "buffer size for each adios_open()...adios_close() operation.\n");
  // adios_buffer_size_requested_set (buffer_size * 1024 * 1024);
  // adios_buffer_alloc_when_set (adios_buffer_alloc_when);
  // adios_set_buffer_size ();
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// Drew: used for experiments
static uint32_t pinned_timestep = 0;
void adios_pin_timestep(uint32_t ts) { pinned_timestep = ts; }

///////////////////////////////////////////////////////////////////////////////
static const char ADIOS_ATTR_PATH[] = "/__adios__";

int common_adios_open(int64_t *fd_p, const char *group_name, const char *name,
                      const char *file_mode, MPI_Comm comm) {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_start("adios_open_to_close");
  timer_start("adios_open");
#endif
  ADIOST_CALLBACK_ENTER(adiost_event_open, *fd_p, group_name, name, file_mode, comm);

  struct adios_group_struct *g = NULL;
  struct adios_method_list_struct *methods = NULL;
  struct adios_file_struct *fd = NULL;
  enum ADIOS_METHOD_MODE mode;

  // printf("beginning of file open... bytes_written=%llu\n",
  // fd->bytes_written);

  adios_errno = err_no_error;
  g = adios_common_get_group(group_name);
  if (!g) {
    adios_error(err_invalid_group,
                "adios_open: try to open file %s with undefined group: %s\n",
                name, group_name);
    *fd_p = 0;
    ADIOST_CALLBACK_EXIT(adiost_event_open, *fd_p, group_name, name, file_mode, comm);
    return adios_errno;
  }

  if (!strcasecmp(file_mode, "r"))
    mode = adios_mode_read;
  else if (!strcasecmp(file_mode, "w"))
    mode = adios_mode_write;
  else if (!strcasecmp(file_mode, "a"))
    mode = adios_mode_append;
  else if (!strcasecmp(file_mode, "u"))
    mode = adios_mode_update;
  else {
    adios_error(err_invalid_file_mode,
                "adios_open: unknown file mode: %s, supported r,w,a,u\n",
                file_mode);

    *fd_p = 0;

    return adios_errno;
  }

  // Time Aggregation: if file name suddenly changes for the group
  // we need to flush the buffer to the old file and start clean
  if (TimeAggregationInProgress(g) && strcmp(name, g->ts_fd->name)) {
    log_debug(
        "TimeAggr: new filename during aggregation. Flush and start buffering "
        "again\n");
    SetTimeAggregationFlush(
        g, 1);  // => TimeAggregationIsFlushing => TimeAggregationLastStep
    common_adios_close(g->ts_fd);   // close old file
    SetTimeAggregationFlush(g, 0);  // => NOT TimeAggregationIsFlushing
    g->ts_fd =
        NULL;  // => TimeAggregationJustBegan AND NOT TimeAggregationInProgress
  }

  // printf("ts_to buffer=%d max_tx=%d\n", g->ts_to_buffer, g->max_ts);
  // Time Aggregation: buffering time steps doesn't need to init file open every
  // time
  if (TimeAggregationInProgress(g)) {
    fd = g->ts_fd;  // continue writing to the previous file
    //   printf("skip file init... bytes_written=%llu file_handle_addr=%llu\n",
    //   fd->bytes_written, fd);
    //  printf("open: fd->pgs_written->pg_start_in_file=%lld,
    //  fd->current_pg->pg_start_in_file=%lld\n",
    //  fd->pgs_written->pg_start_in_file, fd->current_pg->pg_start_in_file);
    log_debug("TimeAggr: skip file name and group assignment\n");
  } else {
    log_debug("TimeAggr: new open... file struct init\n");
    fd = (struct adios_file_struct *)malloc(sizeof(struct adios_file_struct));
    adios_file_struct_init(fd);
    fd->name = strdup(name);
    fd->subfile_index = -1;  // subfile index is by default -1
    fd->group = g;
    fd->mode = mode;
    if (comm == MPI_COMM_NULL)
      fd->comm = MPI_COMM_NULL;
    else if (comm == MPI_COMM_SELF)
      fd->comm = MPI_COMM_SELF;
    else
      MPI_Comm_dup(comm, &fd->comm);
  }

  methods = g->methods;

  // Time Aggregation: if time aggregation is turned on, only open() at the
  // first time step
  if (TimeAggregationInProgress(g)) {
    *fd_p = (int64_t)fd;
    // printf("time buffering.... skip open\n");
    // printf("open: fd->offset=%llu, group_offset=%llu\n", fd->offset,
    // g->group_offset);
  } else {
    while (methods) {
      if (methods->method->m != ADIOS_METHOD_UNKNOWN &&
          methods->method->m != ADIOS_METHOD_NULL &&
          adios_transports[methods->method->m].adios_open_fn) {
        adios_transports[methods->method->m].adios_open_fn(fd, methods->method,
                                                           fd->comm);
      }

      methods = methods->next;
    }

    if (!adios_errno) {
      *fd_p = (int64_t)fd;
    } else {
      // free (fd_p);
      fd_p = 0L;
    }
  }

  /* Time index
   * It increases in write/append mode, it does not change in update/read mode.
   * This piece should be after calling the method's open function, which sets
   * the time index from the existing file in case of append/update. */
  if (mode == adios_mode_write || mode == adios_mode_append) {
    /* Traditionally, time=1 at the first step, and for subsequent file
       creations, time increases. Although, each file contains one step,
       the time index indicates that they are in a series.
     */
    g->time_index++;
  }
  /* time starts from 1 not from 0 in the BP design although it does not have
   * any meaning */
  if (g->time_index == 0) g->time_index = 1;

  // Drew: for experiments
  if (pinned_timestep > 0) g->time_index = pinned_timestep;

  if (!adios_errno && fd->mode != adios_mode_read) {
    /* Add ADIOS internal attributes now */
    if ((fd->group->process_id == 0 || fd->subfile_index != -1)) {
      struct timeval tp;
      char epoch[16];
      gettimeofday(&tp, NULL);
      sprintf(epoch, "%d", (int)tp.tv_sec);

      /* FIXME: this code works fine, it does not duplicate the attribute,
         but the index will still contain all copies and the read will see
         only the first one. Thus updating an attribute does not work
         in practice.
       */

      if (fd->group->time_index == 1) {
        log_debug(
            "Define ADIOS extra attributes, "
            "time = %d, rank = %d, epoch = %s subfile=%d\n",
            fd->group->time_index, fd->group->process_id, epoch,
            fd->subfile_index);

        adios_common_define_attribute((int64_t)fd->group, "version",
                                      ADIOS_ATTR_PATH, adios_string, VERSION,
                                      NULL);

        adios_common_define_attribute((int64_t)fd->group, "create_time_epoch",
                                      ADIOS_ATTR_PATH, adios_integer, epoch,
                                      NULL);
        adios_common_define_attribute((int64_t)fd->group, "update_time_epoch",
                                      ADIOS_ATTR_PATH, adios_integer, epoch,
                                      NULL);
        // id of last attribute is fd->group->member_count
        fd->group->attrid_update_epoch = fd->group->member_count;

      } else {
        // update attribute of update time (define would duplicate it)
        struct adios_attribute_struct *attr = adios_find_attribute_by_id(
            fd->group->attributes, fd->group->attrid_update_epoch);
        if (attr) {
          log_debug(
              "Update ADIOS extra attribute name=%s, "
              "time = %d, rank = %d, epoch = %s, subfile=%d\n",
              attr->name, fd->group->time_index, fd->group->process_id, epoch,
              fd->subfile_index);

          free(attr->value);
          adios_parse_scalar_string(adios_integer, (void *)epoch, &attr->value);
        }
      }
    }

    /* Add first PG to the group */
    if (NotTimeAggregated(g) || TimeAggregationJustBegan(g)) {
      assert(!fd->pgs_written);  // the start of the PG list
      assert(!fd->current_pg);   // the last of the PG list
    }
    // printf("rank %d: fd->bytes-written=%llu\n",fd->group->process_id,
    // fd->bytes_written);
    add_new_pg_written(fd);
    // keep a record of the PG offsets.
    // FIXME: Time Aggregation: Assume every process writes the same amount of
    // data
    if (TimeAggregated(g)) fd->current_pg->pg_start_in_file = fd->bytes_written;

#ifdef ADIOS_TIMERS
    /* Add timer variable definitions to this output */
    adios_add_timing_variables(fd);
#endif

    /* Now ask the methods if anyone wants common-layer BP formatted buffering
     */
    methods = g->methods;

    if (NotTimeAggregated(g) || TimeAggregationJustBegan(g)) {
      while (methods) {
        enum BUFFERING_STRATEGY bufstrat = no_buffering;
        if (methods->method->m != ADIOS_METHOD_UNKNOWN &&
            methods->method->m != ADIOS_METHOD_NULL &&
            adios_transports[methods->method->m].adios_should_buffer_fn) {
          bufstrat =
              adios_transports[methods->method->m].adios_should_buffer_fn(
                  fd, methods->method);
        }

        if (bufstrat != no_buffering) {
          fd->shared_buffer = adios_flag_yes;
          fd->bufstrat = bufstrat;
          /* FIXME: last method determines the value of buffering strategy here.
             This whole
             buffer overflow thing does not work if there are multiple methods
             called
             and they want something else (stop vs continue vs continue with new
             PG
           */
        }

        methods = methods->next;
      }
    }

    if (fd->bufstrat != no_buffering) {
      /* Allocate BP buffer with remembered size or max size or default size */
      uint64_t expected_bufsize;

      if (NotTimeAggregated(g)) {
        if (g->max_pg_size > 0)
          expected_bufsize = g->max_pg_size;
        else
          expected_bufsize = adios_databuffer_get_extension_size(fd);
      } else {
        /* Time Aggregation: set the user given buffer size at the first time
         * step to be buffered.
         * Later, check if the calculated number of timesteps knowing the size
         * encountered so far
         * would still fit. If not, try to increase the buffer to accommodate
         * it. The reason for this is
         * that we must try to avoid buffer overflow in time aggregation, which
         * situation is not handled.
         */
        if (TimeAggregationJustBegan(g)) {
          adios_databuffer_set_max_size(g->ts_buffsize);
          expected_bufsize = g->ts_buffsize;
        } else {
          // fd->bytes_written is the size of (g->max_ts - g->ts_to_buffer)
          // steps already buffered
          expected_bufsize =
              g->max_ts * fd->bytes_written / (g->max_ts - g->ts_to_buffer);
          if (expected_bufsize > fd->buffer_size) {
            adios_databuffer_set_max_size(
                expected_bufsize);  // allow more than what user specified
                                    // (g->ts_buffsize)
          }
        }
      }

      // printf("==open expected_bufsize=%llu\n", expected_bufsize);

      // Time Aggregation: when ts buffering is on, only the first step needs to
      // allocate the buffer here
      if (NotTimeAggregated(g) || TimeAggregationJustBegan(g))
        if (expected_bufsize > fd->buffer_size) {
          if (adios_databuffer_resize(fd, expected_bufsize)) {
            fd->bufstate = buffering_stopped;
            adios_error(err_no_memory,
                        "Cannot allocate %" PRIu64
                        " bytes for buffered output "
                        "of group %s in adios_open(). Output will fail.\n",
                        fd->buffer_size, g->name);
            return adios_errno;
          }
        }
      fd->bufstate = buffering_ongoing;

      // write the process group header
      adios_write_open_process_group_header_v1(fd);

      // setup for writing vars
      adios_write_open_vars_v1(fd);
    }

    // printf("end of adios_open fd->offset=%llu bytes_written=%llu\n",
    // fd->offset, fd->bytes_written);
  }

#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_stop("adios_open");
#endif
  ADIOST_CALLBACK_EXIT(adiost_event_open, *fd_p, group_name, name, file_mode, comm);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////

int common_adios_group_size(int64_t fd_p, uint64_t data_size,
                            uint64_t *total_size) {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_start("adios_group_size");
#endif
  // printf("enter adios_group_size datasize= %llu\n", data_size);
  ADIOST_CALLBACK_ENTER(adiost_event_group_size, fd_p, data_size, total_size);

  adios_errno = err_no_error;
  struct adios_file_struct *fd;
  fd = (struct adios_file_struct *)fd_p;

  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_group_size\n");
    ADIOST_CALLBACK_EXIT(adiost_event_group_size, fd_p, data_size, total_size);
    return adios_errno;
  }

  struct adios_method_list_struct *m = fd->group->methods;
  if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL) {
    *total_size = 0;
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_stop("adios_group_size");
#endif
    ADIOST_CALLBACK_EXIT(adiost_event_group_size, fd_p, data_size, total_size);
    return err_no_error;
  }

  if (fd->buffer_size == 0) {
    *total_size = 0;
    ADIOST_CALLBACK_EXIT(adiost_event_group_size, fd_p, data_size, total_size);
    return err_no_error;
  }

#ifdef ADIOS_TIMERS
  data_size += fd->group->tv_size;
#endif

  uint64_t overhead = adios_calc_overhead_v1(fd);
  *total_size = data_size + overhead;

  // NCSU ALACRITY-ADIOS - Current solution to group_size problem: find
  // the most "expansive" transform method used in the file, and assume
  // all of the data uses that method, for a very rough but safe upper bound.
  //
  // (see the comment in the top of adios_transforms.c, under section
  // 'The "group size" problem,' for more details)
  uint64_t wc_transformed_size =
      adios_transform_worst_case_transformed_group_size(data_size, fd);
  if (wc_transformed_size > data_size) {
    log_debug(
        "Computed worst-case bound on transformed data for a group size of "
        "%" PRIu64 " is %" PRIu64 "; increasing group size to match.\n",
        data_size, wc_transformed_size);

    *total_size += (wc_transformed_size - data_size);
  }

  if (*total_size > fd->buffer_size && fd->bufstate == buffering_ongoing) {
    if (adios_databuffer_resize(fd, *total_size)) {
      log_warn("Cannot reallocate data buffer to %" PRIu64
               " bytes "
               "for group %s in adios_group_size(). Continue buffering "
               "with buffer size %" PRIu64 " MB\n",
               *total_size, fd->group->name, fd->buffer_size / 1048576L);
    }
  }

#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_stop("adios_group_size");
#endif
  // each var will be added to the buffer by the adios_write calls
  // attributes will be added by adios_close

  ADIOST_CALLBACK_EXIT(adiost_event_group_size, fd_p, data_size, total_size);
  return adios_errno;
}

static int common_adios_write_transform_helper(struct adios_file_struct *fd,
                                               struct adios_var_struct *v) {
  int use_shared_buffer = (fd->bufstrat != no_buffering);
  int wrote_to_shared_buffer = 0;

  if (fd->bufstrat == no_buffering) {
    int ret = adios_transform_variable_data(fd, v, use_shared_buffer,
                                            &wrote_to_shared_buffer);
    assert(!wrote_to_shared_buffer);
    assert(v->data);
    return ret;
  } else if (fd->bufstate == buffering_ongoing) {
    // If we are using the shared buffer, transform the data directly into it
    uint16_t header_size = adios_calc_var_overhead_v1(v);
    uint64_t header_offset;
    uint64_t payload_offset;
    uint64_t end_offset;

    // Reserve space for the variable header (it will need to be written after
    // the transform to capture updated metadata)
    header_offset = fd->offset;
    fd->offset += header_size;
    payload_offset = fd->offset;

    // This function will either:
    // a) write to the shared buffer, leave v->data, v->data_size and
    //    v->free_data untouched, and return 1, OR
    // b) write to v->data, set v->data_size and v->free_data, and return 0
    //
    int success = adios_transform_variable_data(fd, v, use_shared_buffer,
                                                &wrote_to_shared_buffer);
    if (!success) {
      fd->offset = header_offset;
      return 0;
    }

    // Assumption: we don't change the header size, just contents, in
    // adios_transform_variable_data
    assert(adios_calc_var_overhead_v1(v) == header_size);

    // Store the ending offset of the payload write (if any)
    end_offset = fd->offset;

    // Rewind and write the header back where it should be
    fd->offset = header_offset;
    // var payload sent for sizing information
    adios_write_var_header_v1(fd, v);

    assert(fd->offset == payload_offset);

    // If the data was stored to the shared buffer, update v->data,
    // v->data_size and v->free_data. Else, write the payload to the shared
    // buffer (the other v->* fields have already been updated)
    if (wrote_to_shared_buffer) {
      v->adata = fd->buffer + payload_offset;
      v->data_size = end_offset - payload_offset;
      v->free_data = adios_flag_no;
      v->data = v->adata;

      // Update the buffer back to the end of the header+payload
      fd->offset = end_offset;
    } else {
      /* FIXME: This branch should not happen for memory reasons.
         The buffer either had enough space to store the transformation
         in place, or buffering has stopped.
         However, identity transform does no buffering, so we do it here.
       */
      // either no transformation happened (original data is in v->data), or
      // data was transformed into v->adata
      if (v->adata) {
        v->data = v->adata;
      }
      // write payload
      adios_write_var_payload_v1(fd, v);

      // fd->offset now points to the end of the header+payload
    }

    // Success!
    return 1;
  }
  /* else: shared buffer but buffering has stopped, nothing needs to be done */
  return 1;
}

///////////////////////////////////////////////////////////////////////////////
/* common_adios_write is just a partial implementation. It expects filled out
 * structures. This is because C and Fortran implementations of adios_write are
 * different for some part and this is the common part.
 */
int common_adios_write(struct adios_file_struct *fd, struct adios_var_struct *v,
                       const void *var) {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_start("adios_write");
#endif
  ADIOST_CALLBACK_WRITE_ENTER(adiost_event_write, fd, v);
  adios_errno = err_no_error;
  struct adios_method_list_struct *m = fd->group->methods;

  // First, before doing any transformation, compute variable statistics,
  // as we can't do this after the data is transformed
  adios_generate_var_characteristics_v1(fd, v);

  uint64_t vsize = 0;
  if (fd->bufstate == buffering_ongoing) {
    // Second, estimate the size needed for buffering and extend buffer when
    // needed
    // and handle the error if buffer cannot be extended
    vsize = adios_transform_worst_case_transformed_var_size(v);

    if (fd->buffer_size < fd->offset + vsize) {
      //    printf("adios_write fd->offset=%llu for variable= %s
      //    buffer_size=%llu\n", fd->offset, v->name, fd->buffer_size);
      //            printf("max_ts=%d ts_to_buffer=%d\n", fd->group->max_ts,
      //            fd->group->ts_to_buffer);
      /* Trouble: this variable does not fit into the current buffer */
      // First, try to realloc the buffer
      uint64_t extrasize = adios_databuffer_get_extension_size(fd);
      if (extrasize < vsize) extrasize = vsize;
      if (adios_databuffer_resize(fd, fd->buffer_size + extrasize)) {
        /* Second, let the method deal with it */
        log_debug(
            "adios_write(): buffer needs to be dumped before buffering "
            "variable %s/%s\n",
            v->path, v->name);
        // these calls don't extend the buffer but we will get a completed PG
        // here
        adios_write_close_vars_v1(fd);
        adios_write_close_process_group_header_v1(fd);

        /* Ask the method to do something with the current buffer then we either
           1. continue buffering from start with a new PG or
           2. skip buffering variables from now on, then method gets the same
           buffer again in close()
         */

        m = fd->group->methods;
        while (m) {
          if (m->method->m != ADIOS_METHOD_UNKNOWN &&
              m->method->m != ADIOS_METHOD_NULL &&
              adios_transports[m->method->m].adios_buffer_overflow_fn) {
            adios_transports[m->method->m].adios_buffer_overflow_fn(fd,
                                                                    m->method);
          }
          m = m->next;
        }

        if (fd->bufstrat == continue_with_new_pg) {
          // special case: fd->buffer_size is smaller than this single variable,
          // and the extension failed:
          // try to extend it to contain this single variable (plus headers) in
          // the next PG
          if (fd->buffer_size < vsize + 1024) {
            if (adios_databuffer_resize(fd, vsize + 1024)) {
              adios_error(
                  err_no_memory,
                  "adios_write(): buffer cannot accommodate variable %s/%s "
                  "with its storage size of %" PRIu64
                  " bytes at all. "
                  "No more variables will be written.\n",
                  v->path, v->name, vsize);
              //"This variable won't be written.\n", v->path, v->name, vsize);
              fd->bufstate = buffering_stopped;
              /* FIXME: This stops all writing, not just this variable! */
              /* FIXME: so maybe we should give the method a chance to write
               * this variable directly? */
            }
          }
          /* Start buffering from scratch (a new PG) */
          fd->offset = 0;
          adios_write_open_process_group_header_v1(fd);
          adios_write_open_vars_v1(fd);
          add_new_pg_written(fd);
        } else if (fd->bufstrat == stop_on_overflow) {
          fd->bufstate = buffering_stopped;
          if (!adios_errno) {
            // method is expected to throw an error in this case but we can
            // ensure it here
            // to signal error upward to not count this var in index
            adios_errno = err_buffer_overflow;
          }
        }
      }
    }
  }

  // If no transform is specified, do the normal thing (write to shared
  // buffer immediately, if one exists)
  if (v->transform_type == adios_transform_none) {
    if (fd->bufstate == buffering_ongoing) {
      /* Now buffer only if we have the buffer for it */
      if (fd->buffer_size > fd->offset + vsize) {
        // var payload sent for sizing information
        adios_write_var_header_v1(fd, v);

        // write payload
        adios_write_var_payload_v1(fd, v);
      }
    }
  } else  // Else, do a transform
  {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_start("adios_transform");
#endif
    ADIOST_CALLBACK_ENTER(adiost_event_transform, (int64_t)fd);
    int success = common_adios_write_transform_helper(fd, v);
    if (success) {
      // Make it appear as if the user had supplied the transformed data
      var = v->data;
    } else {
      log_error(
          "Error: unable to apply transform %s to variable %s; likely ran out "
          "of memory, check previous error messages\n",
          adios_transform_plugin_primary_xml_alias(v->transform_type), v->name);
      // FIXME: Reverse the transform metadata and write raw data as usual
    }
    ADIOST_CALLBACK_EXIT(adiost_event_transform, (int64_t)fd);
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_stop("adios_transform");
#endif
  }

  // now tell each transport attached that it is being written unless buffering
  // is stopped
  if (fd->bufstate == buffering_ongoing || fd->bufstrat == no_buffering) {
    m = fd->group->methods;
    while (m) {
      if (m->method->m != ADIOS_METHOD_UNKNOWN &&
          m->method->m != ADIOS_METHOD_NULL &&
          adios_transports[m->method->m].adios_write_fn) {
        adios_transports[m->method->m].adios_write_fn(fd, v, var, m->method);
      }

      m = m->next;
    }
  } else {
    adios_errno =
        err_buffer_overflow;  // signal error but don't print anything anymore
  }

  if (v->dimensions) {
    // NCSU ALACRITY-ADIOS - Free transform-method-allocated data buffers
    // TODO: Is this correct? Normally v->data is a user buffer, so we
    //   can't free it. However, after a transform, we probably do need
    //   to free it. We mark this by setting the free_data flag. However,
    //   as this flag is hardly ever used, I don't know whether this is
    //   using the flag correctly or not. Need verification with
    //   Gary/Norbert/someone knowledgable about ADIOS internals.
    if (v->transform_type != adios_transform_none &&
        v->free_data == adios_flag_yes && v->adata) {
      free(v->adata);
    }
    v->data = v->adata = 0;  // throw away pointers to user data in case of
                             // arrays to avoid trying
    // to free them in possible forthcoming adios_write() of the same variable
  }

  if (!adios_errno) {
    v->write_count++;
  }
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_stop("adios_write");
#endif
  ADIOST_CALLBACK_WRITE_EXIT(adiost_event_write, fd, v);
  // printf ("var: %s written %d\n", v->name, v->write_count);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_write_byid(struct adios_file_struct *fd,
                            struct adios_var_struct *v, const void *var) {
  struct adios_method_list_struct *m = fd->group->methods;

  ADIOST_CALLBACK_WRITE_ENTER(adiost_event_write_byid, fd, v);
  adios_errno = err_no_error;
  if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL) {
    ADIOST_CALLBACK_WRITE_EXIT(adiost_event_write_byid, fd, v);
    return adios_errno;
  }

  if (v->adata) {
    assert(
        v->dimensions ==
        NULL);  // this must happen only for scalars where we copied the value
    free(v->adata);
    v->adata = 0;
  }

  if (v->dimensions) {
    v->data = var;
  } else {
    uint64_t element_size = adios_get_type_size(v->type, var);

    switch (v->type) {
      case adios_byte:
      case adios_short:
      case adios_integer:
      case adios_long:
      case adios_unsigned_byte:
      case adios_unsigned_short:
      case adios_unsigned_integer:
      case adios_unsigned_long:
      case adios_real:
      case adios_double:
      case adios_long_double:
      case adios_complex:
      case adios_double_complex:
        v->adata = malloc(element_size);
        if (!v->adata) {
          adios_error(
              err_no_memory,
              "In adios_write, cannot allocate %lld bytes to copy scalar %s\n",
              element_size, v->name);
          ADIOST_CALLBACK_WRITE_EXIT(adiost_event_write_byid, fd, v);
          return adios_errno;
        }

        memcpy((char *)v->adata, var, element_size);
        v->data = v->adata;
        break;

      case adios_string:
        v->adata = malloc(element_size + 1);
        if (!v->adata) {
          adios_error(
              err_no_memory,
              "In adios_write, cannot allocate %lld bytes to copy string %s\n",
              element_size, v->name);
          ADIOST_CALLBACK_WRITE_EXIT(adiost_event_write_byid, fd, v);
          return adios_errno;
        }
        ((char *)v->adata)[element_size] = 0;
        memcpy((char *)v->adata, var, element_size);
        v->data = v->adata;
        break;

      default:
        v->data = 0;
        break;
    }
  }

  common_adios_write(fd, v, var);
  // v->data is set to NULL in the above call for arrays,
  // but it's pointing to v->adata, which is allocated in ADIOS, for scalars
  // to remember their value if used as dimension in forthcoming writes of
  // arrays.

  if (!adios_errno && fd->mode != adios_mode_read) {
    adios_copy_var_written(fd, v);
  }

  ADIOST_CALLBACK_WRITE_EXIT(adiost_event_write_byid, fd, v);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_get_write_buffer(int64_t fd_p, const char *name,
                                  uint64_t *size, void **buffer) {
  struct adios_file_struct *fd = (struct adios_file_struct *)fd_p;
  adios_errno = err_no_error;
  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_group_size\n");
    return adios_errno;
  }
  struct adios_var_struct *v = fd->group->vars;
  struct adios_method_list_struct *m = fd->group->methods;

  v = adios_find_var_by_name(fd->group, name);

  if (!v) {
    adios_error(err_invalid_varname, "Bad var name (ignored): '%s' (%c%c%c)\n",
                name, name[0], name[1], name[2]);
    return adios_errno;
  }

  if (fd->mode == adios_mode_read) {
    adios_error(err_invalid_file_mode,
                "write attempted on %s in %s. This was opened for read\n", name,
                fd->name);
    return adios_errno;
  }

  // since we are only getting one buffer, get it from the first
  // transport method that can provide it.
  while (m) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_get_write_buffer_fn) {
      adios_transports[m->method->m].adios_get_write_buffer_fn(
          fd, v, size, buffer, m->method);
      m = 0;
    } else
      m = m->next;
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////

int common_adios_read(int64_t fd_p, const char *name, void *buffer,
                      uint64_t buffer_size) {
  ADIOST_CALLBACK_ENTER(adiost_event_read, fd_p, name, buffer, buffer_size);
  struct adios_file_struct *fd = (struct adios_file_struct *)fd_p;
  adios_errno = err_no_error;
  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_group_size\n");
    ADIOST_CALLBACK_EXIT(adiost_event_read, fd_p, name, buffer, buffer_size);

    return adios_errno;
  }
  struct adios_var_struct *v;
  struct adios_method_list_struct *m = fd->group->methods;

  if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL) {
    // nothing to do so just return
    ADIOST_CALLBACK_EXIT(adiost_event_read, fd_p, name, buffer, buffer_size);
    return err_no_error;
  }

  if (!(fd->mode == adios_mode_read)) {
    adios_error(err_invalid_file_mode,
                "read attempted on %s which was opened for write\n", fd->name);

    ADIOST_CALLBACK_EXIT(adiost_event_read, fd_p, name, buffer, buffer_size);
    return adios_errno;
  }

  v = adios_find_var_by_name(fd->group, name);
  if (v) {
    // since can only read from one place into the buffer,
    // read from the first transport method that can
    while (m) {
      if (m->method->m != ADIOS_METHOD_UNKNOWN &&
          m->method->m != ADIOS_METHOD_NULL &&
          adios_transports[m->method->m].adios_read_fn) {
        adios_transports[m->method->m].adios_read_fn(fd, v, buffer, buffer_size,
                                                     m->method);
        m = 0;
      } else
        m = m->next;
    }
  } else {
    adios_error(err_invalid_varname, "var %s in file %s not found on read\n",
                name, fd->name);

    ADIOST_CALLBACK_EXIT(adiost_event_read, fd_p, name, buffer, buffer_size);
    return adios_errno;
  }

  ADIOST_CALLBACK_EXIT(adiost_event_read, fd_p, name, buffer, buffer_size);
  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// OBSOLETE, kept only for backward compatibility
int common_adios_set_path(int64_t fd_p, const char *path) {
  struct adios_file_struct *fd = (struct adios_file_struct *)fd_p;
  adios_errno = err_no_error;
  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_set_path\n");

    return adios_errno;
  }
  struct adios_group_struct *t = fd->group;
  struct adios_var_struct *v = t->vars;
  struct adios_attribute_struct *a = t->attributes;

  while (v) {
    if (v->path) {
      free(v->path);
    }

    v->path = strdup(path);

    v = v->next;
  }

  while (a) {
    // skip internal attributes
    if (a->path && strstr(a->path, "__adios__")) {
      a = a->next;
      continue;
    }

    if (a->path) {
      free(a->path);
    }

    a->path = strdup(path);

    a = a->next;
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// OBSOLETE, kept only for backward compatibility
// Inconsistent behavior with new ADIOS variable naming
// The variable is not replaced with the new path in the hash table here, so
// it is still found using the old path!
int common_adios_set_path_var(int64_t fd_p, const char *path,
                              const char *name) {
  struct adios_file_struct *fd = (struct adios_file_struct *)fd_p;
  adios_errno = err_no_error;
  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_set_path_var\n");

    return adios_errno;
  }
  struct adios_group_struct *t = fd->group;
  struct adios_var_struct *v = t->vars;

  // check for vars and then attributes
  v = adios_find_var_by_name(t, name);

  if (v) {
    if (v->path) {
      free(v->path);
    }

    v->path = strdup(path);

    /* Possible new behavior: replace the old path with the new path
     * in the hash table so that the variable is found by
     * the new path. Inconsistent with old codes that only used
     * the variable name without the path in adios_write()
     */
    // remove var from hash table by old fullpath...
    // t->hashtbl_vars->remove (t->hashtbl_vars, name);
    // and add it back with new fullpath
    // t->hashtbl_vars->put2 (t->hashtbl_vars, v->path, v->name, v);

  } else {
    adios_error(err_invalid_varname,
                "adios_set_path_var (path=%s, var=%s): var not found\n", path,
                name);

    return adios_errno;
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// hint that we reached the end of an iteration (for asynchronous pacing)
int common_adios_end_iteration() {
  struct adios_method_list_struct *m;

  adios_errno = err_no_error;
  for (m = adios_get_methods(); m; m = m->next) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_end_iteration_fn) {
      adios_transports[m->method->m].adios_end_iteration_fn(m->method);
    }
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// hint to start communicating
int common_adios_start_calculation() {
  struct adios_method_list_struct *m;

  adios_errno = err_no_error;
  for (m = adios_get_methods(); m; m = m->next) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_start_calculation_fn) {
      adios_transports[m->method->m].adios_start_calculation_fn(m->method);
    }
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
// hint to stop communicating
int common_adios_stop_calculation() {
  struct adios_method_list_struct *m;

  adios_errno = err_no_error;
  for (m = adios_get_methods(); m; m = m->next) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_stop_calculation_fn) {
      adios_transports[m->method->m].adios_stop_calculation_fn(m->method);
    }
  }

  return adios_errno;
}

///////////////////////////////////////////////////////////////////////////////
int common_adios_close(struct adios_file_struct *fd) {
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_start("adios_close");
#endif
  adios_errno = err_no_error;

  ADIOST_CALLBACK_ENTER(adiost_event_close, (int64_t)fd);
  if (!fd) {
    adios_error(err_invalid_file_pointer,
                "Invalid handle passed to adios_close\n");

    ADIOST_CALLBACK_EXIT(adiost_event_close, (int64_t)fd);
    return adios_errno;
  }

  struct adios_method_list_struct *m = fd->group->methods;
  if (m && m->next == NULL && m->method->m == ADIOS_METHOD_NULL) {
// nothing to do so just return
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
    timer_stop("adios_close");
    timer_stop("adios_open_to_close");
#endif
    ADIOST_CALLBACK_EXIT(adiost_event_close, (int64_t)fd);
    return 0;
  }

#ifdef ADIOS_TIMERS
  if (fd->mode != adios_mode_read) {
    adios_write_timing_variables(fd);
  }
#endif

  struct adios_attribute_struct *a = fd->group->attributes;
  struct adios_var_struct *v = fd->group->vars;

  if (fd->mode != adios_mode_read && !TimeAggregationIsFlushing(fd->group)) {
    // printf("=====common_adios_close offset=%llu %llu\n max_ts=%d
    // ts_to_buffer=%d \n", fd->pgs_written->pg_start_in_file,
    // fd->current_pg->pg_start_in_file, fd->group->max_ts,
    // fd->group->ts_to_buffer);
    //        fd->current_pg=fd->pgs_written;
    if (fd->bufstate == buffering_ongoing) {
      adios_write_close_vars_v1(fd);
    }

    /* Even if buffering has stopped on overflow, let's try to squeeze in the
       attributes
       to the end of the buffer */
    if (fd->bufstrat != no_buffering) {
      /* FIXME: this strategy writes all attributes defined in time step 0
         and thus duplicates them in the PGs and in the attribute index.
         For write mode, where files are new, this is good.
         For append/update it is unnecessary and replicates the attributes
         in the index.
         One should write the newly created attributes only in append mode.
       */
      uint64_t asize = 12;  // 12 bytes to store empty attribute list in a PG
                            // (the opening bytes)
      if (!fd->group->process_id || fd->subfile_index != -1) {
        // from ADIOS 1.4, only rank 0 writes attributes (or each process to its
        // subfile)
        asize +=
            adios_calc_attrs_overhead_v1(fd);  // size of full attribute list
      }
      if (fd->buffer_size < fd->offset + asize) {
        /* Trouble, attributes just don't fit into the end of buffer */
        // try to extend the buffer first
        log_debug(
            "Need more space for attributes in close(). "
            "Current buffer usage=%" PRIu64 " Attributes need %" PRIu64
            " bytes "
            "var_start offset=%" PRIu64 " and bytes_written=%" PRIu64 "\n",
            fd->offset, asize, fd->vars_start, fd->bytes_written);

        if (adios_databuffer_resize(fd, fd->offset + asize)) {
          /* FIXME:  Well, drop them in this case */
          log_error(
              "adios_close(): There is not enough buffer to write the "
              "attributes. "
              "They will be missing from the output\n");
        }
      }

      if (fd->buffer_size >= fd->offset + asize) {
        adios_write_open_attributes_v1(
            fd);  // this writes 12 bytes into the buffer
        if (!fd->group->process_id || fd->subfile_index != -1) {
          while (a) {
            adios_write_attribute_v1(fd, a);
            a = a->next;
          }
        }
        adios_write_close_attributes_v1(fd);
      }

      adios_write_close_process_group_header_v1(fd);
    }
  }

  if (TimeAggregationJustBegan(fd->group)) {
    // first close in time aggregation: we need to figure out how many
    // time steps we can buffer given the buffer size
    if (fd->group->ts_to_buffer > 0 && fd->group->ts_buffsize > 0) {
      fd->group->max_ts = fd->group->ts_buffsize / fd->bytes_written;
      // max_ts may be different on different processes, we need to use the
      // minimum of them
      int min_max_ts;
      MPI_Allreduce(&fd->group->max_ts, &min_max_ts, 1, MPI_INT, MPI_MIN,
                    fd->comm);
      fd->group->max_ts = min_max_ts;
      fd->group->ts_to_buffer =
          fd->group->max_ts - 1;  // we already have one in our hand
      //    printf("ts_buffsize=%llu bytes_writen=%llu\n",
      //    fd->group->ts_buffsize, fd->bytes_written);
      //    printf("max_ts=%d ts_to_buffer=%d\n", fd->group->max_ts,
      //    fd->group->ts_to_buffer);
    }
    // if we can only have one step in the buffer, now TimeAggregationLastStep()
    // holds.
    // TimeAggregationJustBegan() still holds (fd->group->ts_fd == NULL).

    fd->group->ts_fd = fd;
    // TimeAggregationJustBegan() now moved into TimeAggregationInProgress()
  }

  // in order to get the index assembled, we need to do it in the
  // transport once we have collected all of the pieces

  // now tell all of the transports to write the buffer during close
  for (; m; m = m->next) {
    if (m->method->m != ADIOS_METHOD_UNKNOWN &&
        m->method->m != ADIOS_METHOD_NULL &&
        adios_transports[m->method->m].adios_close_fn) {
      if (NotTimeAggregated(fd->group)) {
        adios_transports[m->method->m].adios_close_fn(fd, m->method);
      } else {
        if (!TimeAggregationIsFlushing(fd->group)) {
          // Time Aggregation: build index for current step and merge it in to
          // the aggregated index
          struct adios_index_struct_v1 *current_index;
          current_index = adios_alloc_index_v1(1);
          adios_build_index_v1(fd, current_index);
          // printf("current->index time=%d offset=%llu\n",
          // current_index->vars_root->characteristics[0].time_index,
          // current_index->vars_root->characteristics[0].offset);
          if (fd->group->index == NULL) {
            // first time step, we just built the index for the group
            fd->group->index = current_index;
          } else {
            adios_merge_index_v1(fd->group->index, current_index->pg_root,
                                 current_index->vars_root,
                                 current_index->attrs_root, 1);
            adios_free_index_v1(current_index);
          }
          // printf("fd->group->index time=%d offset=%llu\n",
          // fd->group->index->vars_root->characteristics[1].time_index,
          // fd->group->index->vars_root->characteristics[1].offset);
        }

        if (TimeAggregationLastStep(fd->group)) {
          // in last step, move the PG pointer to the first pg
          // to build correct starting offset of this process for all aggregated
          // steps
          fd->current_pg = fd->pgs_written;
          fd->group->built_index = 1;
          // printf("bytes-writen=%llu time index=%lu\n", fd->bytes_written,
          // fd->group->index->pg_root->time_index);

          adios_transports[m->method->m].adios_close_fn(fd, m->method);

          /* the method clears the full index data including the one built here.
           * we need to just free the container structure.
           */
          adios_free_index_v1(fd->group->index);
          fd->group->index = NULL;
        }
      }
    }
  }

  /* Synced groups for time-aggregation.
   * Force close on those groups that are time-aggregated and are synced with
   * this group
   */
  if (NotTimeAggregated(fd->group) || TimeAggregationLastStep(fd->group)) {
    if (TimeAggregationIsaSyncGroup(fd->group)) {
      struct adios_group_struct **syncedgroups;
      int ngroups, i;
      TimeAggregationGetSyncedGroups(fd->group, &syncedgroups, &ngroups);
      for (i = 0; i < ngroups; i++) {
        struct adios_group_struct *sg = syncedgroups[i];
        if (sg->ts_fd != NULL) {
          if (fd->group->process_id == 0) {
            log_info(
                "Sync flush group '%s' because we just wrote group '%s'. "
                "Synced group size is currently %" PRIu64
                " bytes holding %d steps\n",
                sg->name, fd->group->name, sg->ts_fd->bytes_written,
                sg->max_ts - sg->ts_to_buffer - 1);
          }
          SetTimeAggregationFlush(sg, 1);
          // => TimeAggregationIsFlushing => TimeAggregationLastStep
          common_adios_close(sg->ts_fd);  // close file for sync'd group
          SetTimeAggregationFlush(sg, 0);
        }
      }
    }
  }

  // Time Aggregation: only free the buffer if 1) we are not buffering time
  // steps;
  // 2) we have enough time steps buffered
  // printf("close: do_ts_aggr=%d  max_ts=%d ts_to_buffer=%d\n",
  // fd->group->do_ts_aggr, fd->group->max_ts, fd->group->ts_to_buffer);
  if (NotTimeAggregated(fd->group) || TimeAggregationLastStep(fd->group)) {
    while (v) {
      v->write_offset = 0;
      if (v->adata) {
        free(v->adata);
        v->data = v->adata = 0;
      }

      v = v->next;
    }

    /* clean-up all copied variables with statistics and data in all PGs
     * attached to this file */
    adios_free_pglist(fd);

    if (fd->name) {
      free(fd->name);
      fd->name = 0;
    }

    if (fd->comm != MPI_COMM_NULL && fd->comm != MPI_COMM_SELF) {
      MPI_Comm_free(&fd->comm);
    }
  }

#ifdef ADIOS_TIMER_EVENTS
  /* This has to come here before freeing the buffer that zeros out
   * fd->bytes_written */
  char *extension = ".perf";
  int name_len = strlen(fd->name);
  int fn_len = name_len + strlen(extension) + 1;
  char *fn = (char *)malloc(sizeof(char) * fn_len);

  sprintf(fn, "%s%s", fd->name, extension);

  adios_timing_write_xml_common(fd_p, fn);
#endif

  // Time Aggregation: add time step aggregation logic
  // ts aggregation will keep the buffer for next time steps
  // only free the buffer when ts aggregation is not required, or it
  // has reached the MAX_TS.
  // printf("close: bytes_written=%llu buffer_siz=%llu
  // fd->group->max_pg_size=%llu\n", fd->bytes_written, fd->buffer_size,
  // fd->group->max_pg_size);
  if (fd->bufstrat != no_buffering) {
    // printf("rank %d group %s last buffer size = %" PRIu64 " buffer_size = %"
    // PRIu64 "\n",
    //        fd->group->process_id, fd->group->name, fd->group->max_pg_size,
    //        fd->buffer_size);

    // Remember how much buffer we used for this group for helping buffer
    // allocation in the next cycle.
    // This only works for non-time-aggregation. In time aggregation,
    // bytes_written is the current size of all buffered steps
    if (NotTimeAggregated(fd->group)) {
      if (fd->group->max_pg_size < fd->bytes_written) {
        fd->group->max_pg_size = fd->bytes_written;
      }
    }
    if (NotTimeAggregated(fd->group) || TimeAggregationLastStep(fd->group)) {
      adios_databuffer_free(fd);
    }
  }

  if (TimeAggregated(fd->group)) {
    if (TimeAggregationLastStep(fd->group)) {
      // Time Aggregation: reset the ts counter when the last ts is in
      // printf("group %s reset ts_to_buffer to %d\n", fd->group->name,
      // fd->group->max_ts);
      fd->group->ts_to_buffer = fd->group->max_ts;
      fd->group->ts_fd = NULL;
      free(fd);
    } else {
      // fd->group->new_write_offset=fd->current_pg->pg_start_in_file+
      fd->group->ts_to_buffer--;
      // TimeAggregationInProgress() may have moved into
      //   TimeAggregationInProgress && TimeAggregationLastStep (i.e.
      //   ts_to_buffer==0)
      // printf("didn't free fd, continue buffering addr=%lld\n", fd_p);
    }
  } else {
    free(fd);
  }

// printf("close file==============\n");
/*
printf("close: do_ts_aggr=%d  max_ts=%d ts_to_buffer=%d\n",
        fd->group->do_ts_aggr, fd->group->max_ts, fd->group->ts_to_buffer);

printf("close: fd->offset=%lld, group_offset=%lld\n", fd->offset,
        fd->group->group_offset);

printf("close: fd->pgs_written->pg_start_in_file=%lld,
fd->current_pg->pg_start_in_file=%lld\n",
        fd->pgs_written->pg_start_in_file, fd->current_pg->pg_start_in_file);
 */
#if defined(WITH_NCSU_TIMER) && defined(TIMER_LEVEL) && (TIMER_LEVEL <= 0)
  timer_stop("adios_close");
  timer_stop("adios_open_to_close");
  //    printf ("Timers, ");
  //    printf ("%d, ", fd->group->process_id);
  //    printf ("%d, ", fd->group->time_index);
  //    printf ("%lf, ", timer_get_total_interval ("adios_open" ));
  //    printf ("%lf, ", timer_get_total_interval ("adios_group_size"));
  //    printf ("%lf, ", timer_get_total_interval ("adios_transform" ));
  //    printf ("%lf, ", timer_get_total_interval ("adios_write" ));
  //    printf ("%lf\n", timer_get_total_interval ("adios_close"     ));
  //    timer_reset_timers ();

  printf("[TIMERS] Proc: %d Time: %d ", fd->group->process_id,
         fd->group->time_index);
  int i;
  timer_result_t *results = timer_get_results_sorted();
  for (i = 0; i < timer_get_num_timers(); i++) {
    printf("%s: %0.4lf ", results[i].name, results[i].time);
  }
  printf("\n");
  free(results);

// timer_reset_timers ();
#endif

  ADIOST_CALLBACK_EXIT(adiost_event_close, (int64_t)fd);
  return adios_errno;
}

//////////////////////////////////////////////////////////////////////////////
// Methods normally only called by the XML parser
//////////////////////////////////////////////////////////////////////////////

// adios_common_declare_group is in adios_internals.c
// adios_common_define_var is in adios_internals.c
// adios_common_define_attribute is in adios_internals.c
// adios_common_select_method is in adios_internals.c
