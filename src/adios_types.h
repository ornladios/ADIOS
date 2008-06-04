/* global defines needed for the type creation/setup functions (see bottom) */
enum ADIOS_DATATYPES {adios_unknown = -1
                     ,adios_integer = 2  // bp_int
                     ,adios_long = 4     // bp_longlong
                     ,adios_real = 5     // bp_float
                     ,adios_string = 9   // bp_string
                     ,adios_double = 6   // bp_double
                     ,adios_complex = 10 // bp_complex
                     ,adios_byte = 50    // bp_uchar
                     };

enum ADIOS_FLAG {adios_flag_unknown = 0
                ,adios_flag_yes = 1
                ,adios_flag_no = 2
                };
