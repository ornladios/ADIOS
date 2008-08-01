#ifndef ADIOS_TYPES_H
#define ADIOS_TYPES_H

// global defines needed for the type creation/setup functions
enum ADIOS_DATATYPES {adios_unknown = -1             // bp type         (size)

                     ,adios_byte = 0                 // bp_char           (1)
                     ,adios_short = 1                // bp_short          (2)
                     ,adios_integer = 2              // bp_int            (4)
                     ,adios_long = 4                 // bp_longlong       (8)

                     ,adios_unsigned_byte = 50       // bp_uchar          (1)
                     ,adios_unsigned_short = 51      // bp_ushort         (2)
                     ,adios_unsigned_integer = 52    // bp_uint           (4)
                     ,adios_unsigned_long = 54       // bp_ulonglong      (8)

                     ,adios_real = 5                 // bp_float          (4)
                     ,adios_double = 6               // bp_double         (8)
                     ,adios_long_double = 7          // bp_longdouble     (16)

                     ,adios_string = 9               // bp_string         (?)
                     ,adios_complex = 10             // bp_complex        (8)
                     ,adios_double_complex = 11      // bp_double_complex (16)
                     };

enum ADIOS_FLAG {adios_flag_unknown = 0
                ,adios_flag_yes = 1
                ,adios_flag_no = 2
                };

#endif
