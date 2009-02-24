#ifndef ADIOS_H
#define ADIOS_H

#include "adios_types.h"
#include <stdint.h>

// ADIOS - Adaptable IO System

// Global setup using the XML file
int adios_init (const char * config);
void adios_init_ (const char * config, int* err, int config_size);

// setup, but all XML file pieces will be provided by another series of calls
// yet to be worked out
// TODO
int adios_init_local (void);
void adios_init_local_ (int * err);

int adios_finalize (int mype);
void adios_finalize_ (int * mype, int * err);

int adios_allocate_buffer (void);
void adios_allocate_buffer_ (int * err);

// end user calls for each I/O operation
// modes = "r" = "read", "w" = "write", "a" = "append", "u" = "update"
int adios_open (int64_t * fd, const char * group_name, const char * name
               ,const char * mode
               );
void adios_open_ (int64_t * fd, const char * group_name, const char * name
                 ,const char * mode, int * err
                 ,int group_name_size, int name_size, int mode_size
                 );

int adios_group_size (int64_t fd_p, uint64_t data_size
                     ,uint64_t * total_size, void * comm
                     );
void adios_group_size_ (int64_t * fd_p, int64_t * data_size
                       ,int64_t * total_size, void * comm, int * err
                       );

int adios_write (int64_t fd_p, const char * name, void * var);
void adios_write_ (int64_t * fd_p, const char * name, void * var, int * err
                  ,int name_size
                  );

int adios_get_write_buffer (int64_t fd_p, const char * name
                           ,uint64_t * size
                           ,void ** buffer
                           );
void adios_get_write_buffer_ (int64_t * fd_p, const char * name
                             ,int64_t * size
                             ,void ** buffer, int * err, int name_size
                             );

int adios_read (int64_t fd_p, const char * name, void * buffer
               ,uint64_t buffer_size
               );
void adios_read_ (int64_t * fd_p, const char * name, void * buffer
                 ,int64_t * buffer_size, int * err, int name_size
                 );

int adios_set_path (int64_t fd_p, const char * path);
void adios_set_path_ (int64_t * fd_p, const char * path, int * err
                     ,int path_size
                     );

int adios_set_path_var (int64_t fd_p, const char * path, const char * name);
void adios_set_path_var_ (int64_t * fd_p, const char * path
                         ,const char * name, int * err
                         ,int path_size, int name_size
                         );

int adios_end_iteration (void);
void adios_end_iteration_ (int * err);

int adios_start_calculation (void);
void adios_start_calculation_ (int * err);

int adios_stop_calculation (void);
void adios_stop_calculation_ (int * err);

int adios_close (int64_t fd_p);
void adios_close_ (int64_t * fd_p, int * err);

// Generally internal use called when parsing the XML file
int adios_declare_group (int64_t * id, const char * name
                        ,const char * coordination_comm
                        ,const char * coordination_var
                        ,const char * time_index
                        );
void adios_declare_group_ (int64_t * id, const char * name
                          ,const char * coordination_comm
                          ,const char * coordination_var
                          ,const char * time_index, int * err
                          ,int name_size, int coordination_comm_size
                          ,int coordination_var_size, int time_index_size
                          );

int adios_define_var (int64_t group_id, const char * name
                     ,const char * path, int type
                     ,const char * dimensions
                     ,const char * global_dimensions
                     ,const char * local_offsets
                     );
void adios_define_var_ (int64_t * group_id, const char * name
                       ,const char * path, int * type
                       ,const char * dimensions
                       ,const char * global_dimensions
                       ,const char * local_offsets, int * err
                       ,int name_size, int path_size, int dimensions_size
                       ,int global_dimensions_size, int local_offsets_size
                       );

int adios_define_attribute (int64_t group, const char * name
                           ,const char * path, enum ADIOS_DATATYPES type
                           ,const char * value, const char * var
                           );
void adios_define_attribute_ (int64_t * group, const char * name
                             ,const char * path, int type, const char * value
                             ,const char * var, int * err
                             ,int name_size, int path_size, int value_size
                             ,int var_size
                             );

int adios_select_method (int priority, const char * method
                        ,const char * parameters, const char * type
                        ,const char * base_path, int iters
                        );
void adios_select_method_ (int * priority, const char * method
                          ,const char * parameters, const char * type
                          ,const char * base_path, int * iters, int * err
                          ,int method_size
                          ,int parameters_size, int type_size
                          ,int base_path_size
                          );

int adios_select_method (int priority, const char * method
                        ,const char * parameters, const char * type
                        ,const char * base_path, int iters
                        );
void adios_select_method_ (int * priority, const char * method
                          ,const char * parameters, const char * type
                          ,const char * base_path, int * iters, int * err
                          ,int method_size
                          ,int parameters_size, int type_size
                          ,int base_path_size
                          );

#endif
