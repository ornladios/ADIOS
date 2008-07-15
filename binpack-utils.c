#include "binpack-utils.h"
#include "adios_types.h"

#include <stdint.h>

extern int adios_host_language_fortran;

char** bp_dirparser(char *str, int *nLevel)
{
  char **grp_name;
  char *pch;
  int idx = 0, len=0;
  char *tmpstr;
  tmpstr= (char *)malloc(sizeof(char)*(strlen(str)+1));
  strcpy(tmpstr,str);
  pch = strtok(tmpstr,"/");
  grp_name = (char **)malloc(NUM_GP*sizeof(char));
  while(pch!=NULL && *pch!=' ')
  {
     
     len = strlen(pch);
     grp_name[idx]  = (char *)malloc((len+1)*sizeof(char));
     grp_name[idx][0]='\0';
     strcat(grp_name[idx],pch);
     pch=strtok(NULL,"/");
     idx=idx+1;
  }
  *nLevel = idx;
  free(tmpstr);
  return grp_name;
}

// Returns the size in bytes of each possible vartype_t
int bp_getsize(enum vartype_t type, void * val)
{
    if (!val)
        return 0;

    switch (type)
    {
        case bp_char:
        case bp_uchar:
            return 1;

        case bp_string:
            return strlen ((char *) val);

        case bp_short:
        case bp_ushort:
            return 2;

        case bp_int:
        case bp_uint:
        case bp_long:
        case bp_ulong:
            return 4;

        case bp_longlong:
        case bp_ulonglong:
            return 8;

        case bp_float:
            return 4;

        case bp_double:
            return 8;

        case bp_long_double:
            return 16;

        case bp_pointer:
            return sizeof(char*);

        case bp_complex:
            return 2 * 4;

        case bp_double_complex:
            return 2 * 8;

        default:
            return -1;
    }
}

// Given the cell type, return the number of nodes per cell
// and the name of the cell type
int bp_getCellInfo(enum cell_t type, char* name)
{
  switch (type)
  {
    case point_cell:
      strcpy(name, "Point");
      return 1;
      break;
    case line_cell:
      strcpy(name, "Line");
      return 2;
      break;
    case tri_cell:
      strcpy(name, "Tri");
      return 3;
      break;
    case quad_cell:
      strcpy(name, "Quad");
      return 4;
      break;
    case tet_cell:
      strcpy(name, "Tet");
      return 4;
      break;
    case hex_cell:
      strcpy(name, "Hex");
      return 8;
      break;
    case pyr_cell:
      strcpy(name, "Pyr");
      return 5;
      break;
    case prism_cell:
      strcpy(name, "Prism");
      return 6;
      break;
    default:
      name[0] = '\0';
      return 50;
  }
}

// Read the next int. from packfile without advancing the stream position
int bp_tagPeek(FILE* packfile)
{
  if (packfile)
  {
    fpos_t position;
    int nexttag;
    
    if (fgetpos(packfile, &position) == 0) {
      fread(&nexttag, sizeof(int), 1, packfile);
      if (fsetpos(packfile, &position) != 0)
        return -2;
      else
        return nexttag;
    }
    else
      return -2;
  }
  else
    return -10;
}

void adios_var_element_count (int rank
                             ,struct adios_bp_dimension_struct * dims
                             ,uint64_t * use_count
                             ,uint64_t * total_count
                             )
{
    int i;

    *use_count = 1;
    *total_count = 1;

    for (i = 0; i < rank; i++)
    {
        *use_count *= dims [i].local_bound;
        *total_count *= dims [i].local_bound;
#if 0
        int use_size = dims [i].use_upper_bound - dims [i].use_lower_bound + 1;
        int total_size = dims [i].upper_bound - dims [i].lower_bound + 1;

        // adjust for the stride
        if (dims [i].stride <= use_size / 2)
        {
            if (use_size % dims [i].stride == 1) // correct fencepost error
                use_size = use_size / dims [i].stride + 1;
            else
                use_size = use_size / dims [i].stride;
        }
        else
        {
            if (dims [i].stride >= use_size)
                use_size = 1;
            else
                use_size = use_size / dims [i].stride + 1;  // maybe always 2?
        }

        // need to correct for empty/bad arrays
        if (   dims [i].use_upper_bound < dims [i].use_lower_bound
            || dims [i].upper_bound < dims [i].lower_bound
           )
        {
            use_size = 0;
            total_size = 0;
        }

        // number of items in the array
        *use_count *= use_size;
        *total_count *= total_size;
#endif
    }
}

// when reading or writing an array where we only want a subset, we need
// to know, for each entry, do we use it or not.  This function takes
// an initial description of the dimensions, builds a counter in the
// position parameter to iterate through all of the array entries.  On
// each call, it iterates one item, tests to see if it is in bounds to be
// used, and returns 1 for yes and 0 for no.
int adios_should_use_data (int element, int rank
                          ,struct adios_bp_dimension_struct * dims
                          ,int * position
                          )
{
#if 0
    int i;
    
    if (element == 0)
    {
        for (i = 0; i < rank; i++)
        {
            position [i] = dims [i].lower_bound;
        }   
    }   
    else  // increment our position
    {
        if (adios_host_language_fortran == adios_flag_yes)
        {
            i = 0;
            int done = 0;
            while (!done && i < rank)
            {
                // if less than max, just increment this dim
                if (position [i] < dims [i].upper_bound)
                {
                    position [i]++;
                    done = 1;
                }
                else  // reset dim and move to next to increment
                {
                    position [i] = dims [i].lower_bound;
                    i++;
                }
            }
        }
        else
        {
            i = rank - 1;
            int done = 0;
            while (!done && i >= 0)
            {
                // if less than max, just increment this dim
                if (position [i] < dims [i].upper_bound)
                {
                    position [i]++;
                    done = 1;
                }
                else  // reset dim and move to next to increment
                {
                    position [i] = dims [i].lower_bound;
                    i--;
                }
            }
        }
    }

    // check against bounds
    for (i = 0; i < rank; i++)
    {
        if (   position [i] < dims [i].use_lower_bound
            || position [i] > dims [i].use_upper_bound
           )
        {
            return 0;
        }
        else
        {
            // (pos - use lower) mod stride == 0 == use this element
            if (((position [i] - dims [i].use_lower_bound) % dims [i].stride) != 0)
                return 0;
        }
    }
#endif

    return 1;  // we only get here if we are within all bounds
}

const char * adios_tag_to_string (int tag)
{
    switch (tag)
    {
        case SCR_TAG:        return "Scalar";
        case DST_TAG:        return "Dataset";
        case GRP_TAG:        return "Group";
        case DIR_TAG:        return "Path";
        case VAL_TAG:        return "Value";
        case DSTVAL_TAG:     return "Dataset value";
        case DSTATRN_TAG:    return "Data attribute numeric";
        case DSTATRS_TAG:    return "Data attribute string";
        case GRPATRN_TAG:    return "Group attribute numeric";
        case GRPATRS_TAG:    return "Group attribute string";
        case NULL_TAG:       return "(null)";
        default:             return "(error)";
    }
}

const char * adios_type_to_string (int type)
{
    switch (type)
    {
        case bp_uchar:          return "unsigned byte";
        case bp_ushort:         return "unsigned short";
        case bp_uint:           return "unsigned integer";
        case bp_ulonglong:      return "unsigned long long";

        case bp_char:           return "byte";
        case bp_short:          return "short";
        case bp_int:            return "integer";
        case bp_longlong:       return "long long";

        case bp_float:          return "real";
        case bp_double:         return "double";
        case bp_long_double:    return "long double";

        case bp_string:         return "string";
        case bp_complex:        return "complex";
        case bp_double_complex: return "double complex";

        default:
        {
            static char buf [50];
            sprintf (buf, "(unknown: %d)", type);
            return buf;
        }
    }
}

const char * adios_file_mode_to_string (int mode)
{
    static char buf [50];

    switch (mode)
    {
        case 1: return "write";
        case 2: return "read";
        case 3: return "update";
        case 4: return "append";

        default:
            sprintf (buf, "(unknown: %d)", mode);
    }

    return buf;
}
