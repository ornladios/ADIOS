/* Auto-generated on Tue Feb 25 12:27:40 EST 2014 */

typedef struct icee_varinfo_rec {
    char* varname;
    int varid;
    int type;
    int typesize;
    int ndims;
    uint64_t* gdims;
    uint64_t* ldims;
    uint64_t* offsets;
    uint64_t varlen;
    char* data;
    struct icee_varinfo_rec * next;
} icee_varinfo_rec_t, *icee_varinfo_rec_ptr_t;

static FMField icee_varinfo_field_list[] =
{
    {"varname", "string", sizeof(char*), FMOffset(icee_varinfo_rec_ptr_t, varname)},
    {"varid", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, varid)},
    {"type", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, type)},
    {"typesize", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, typesize)},
    {"ndims", "integer", sizeof(int), FMOffset(icee_varinfo_rec_ptr_t, ndims)},
    {"gdims", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, gdims)},
    {"ldims", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, ldims)},
    {"offsets", "integer[ndims]", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, offsets)},
    {"varlen", "integer", sizeof(uint64_t), FMOffset(icee_varinfo_rec_ptr_t, varlen)},
    {"data", "char[varlen]", sizeof(char), FMOffset(icee_varinfo_rec_ptr_t, data)},
    {"next", "*icee_varinfo", sizeof(struct icee_varinfo_rec ), FMOffset(icee_varinfo_rec_ptr_t, next)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec icee_varinfo_format_list[] =
{
    {"icee_varinfo", icee_varinfo_field_list, sizeof(icee_varinfo_rec_t), NULL},
    {NULL, NULL}
};

/* Auto-generated on Tue Feb 25 12:27:40 EST 2014 */

typedef struct icee_fileinfo_rec {
    char* fname;
    int nvars;
    int rank;
    struct icee_varinfo_rec * varinfo;
    struct icee_fileinfo_rec * next;
} icee_fileinfo_rec_t, *icee_fileinfo_rec_ptr_t;

static FMField icee_fileinfo_field_list[] =
{
    {"fname", "string", sizeof(char*), FMOffset(icee_fileinfo_rec_ptr_t, fname)},
    {"nvars", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, nvars)},
    {"rank", "integer", sizeof(int), FMOffset(icee_fileinfo_rec_ptr_t, rank)},
    {"varinfo", "*icee_varinfo", sizeof(struct icee_varinfo_rec ), FMOffset(icee_fileinfo_rec_ptr_t, varinfo)},
    {"next", "*icee_fileinfo", sizeof(struct icee_fileinfo_rec ), FMOffset(icee_fileinfo_rec_ptr_t, next)},
    {NULL, NULL, 0, 0}
};

static FMStructDescRec icee_fileinfo_format_list[] =
{
    {"icee_fileinfo", icee_fileinfo_field_list, sizeof(icee_fileinfo_rec_t), NULL},
    {"icee_varinfo", icee_varinfo_field_list, sizeof(icee_varinfo_rec_t), NULL},
    {NULL, NULL}
};

