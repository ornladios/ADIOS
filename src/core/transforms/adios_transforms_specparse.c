/*
 * adios_transforms_specparse.c
 *
 *  Created on: May 20, 2013
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "core/transforms/adios_transforms_specparse.h"

inline static char * strsplit(char *input, char split) {
    char *pos = strchr(input, split);
    if (!pos)
        return NULL;

    *pos = '\0';
    return pos + 1;
}

inline static int strcount(char *input, char chr) {
    int count = 0;
    while ((input = strchr(input, chr))) {
        count++;
        input++;
    }
    return count;
}

// var is a pointer type
#define MALLOC_ARRAY(var, count) ((var) = (typeof(var))malloc(sizeof(*var) * (count)))
#define MALLOC_VAR(var) MALLOC_ARRAY(var, 1)
#define CALLOC_ARRAY(var, count) ((var) = (typeof(var))calloc((count), sizeof(*var)))
#define CALLOC_VAR(var) CALLOC_ARRAY(var, 1)

struct adios_transform_spec * adios_transform_parse_spec(const char *spec_str) {
    //struct adios_transform_spec *spec = (struct adios_transform_spec *)malloc(sizeof(struct adios_transform_spec));
    struct adios_transform_spec *spec;
    MALLOC_VAR(spec);

    *spec = (struct adios_transform_spec){
        .transform_type = adios_transform_none,
        .transform_type_str = "",
        .param_count = 0,
        .params = NULL,
        .backing_str = NULL,
        .backing_str_len = 0,
    };

    // If the spec string is null/empty, stop now with a blank spec
    if (!spec_str || strcmp(spec_str, "") == 0)
        return spec;
    assert(spec_str && strcmp(spec_str, "") != 0);

    // Duplicate the spec string so we can chop it up
    char *new_spec_str = strdup(spec_str);
    spec->backing_str = new_spec_str;
    spec->backing_str_len = strlen(new_spec_str);

    // Mark the transform method string in the spec string (the beginning)
    spec->transform_type_str = new_spec_str;

    // Split off the parameters if present
    char *param_list = strsplit(new_spec_str, ':');

    // Parse the transform method string
    spec->transform_type = adios_transform_find_type_by_xml_alias(spec->transform_type_str);

    // If the transform type is unknown (error) or none, stop now and return
    if (spec->transform_type == adios_transform_unknown ||
        spec->transform_type == adios_transform_none)
        return spec;
    assert(spec->transform_type != adios_transform_unknown &&
           spec->transform_type != adios_transform_none);

    // If there is no parameter list, we are done
    if (!param_list)
        return spec;
    assert(param_list);

    spec->param_count = strcount(param_list, ',') + 1;
    MALLOC_ARRAY(spec->params, spec->param_count);
    //spec->params = (typeof(spec->params))malloc(sizeof(*spec->params));

    struct adios_transform_spec_kv_pair *cur_kv = spec->params;
    char *cur_param;
    while (param_list) {
        cur_param = param_list;
        param_list = strsplit(param_list, ',');

        cur_kv->key = cur_param;
        cur_kv->value = strsplit(cur_param, '='); // NULL if no =

        cur_kv++;
    }

    return spec;
}

struct adios_transform_spec * adios_transform_spec_copy(struct adios_transform_spec *src) {
    struct adios_transform_spec *dst;
    CALLOC_VAR(dst);

    dst->transform_type = src->transform_type;

    // Duplicate the backing string if needed, then rebase all strings pointing into it
    if (src->backing_str) {
        // Duplicate the backing string
        dst->backing_str_len = src->backing_str_len;
        dst->backing_str = (char *)malloc(dst->backing_str_len + 1);
        memcpy(dst->backing_str, src->backing_str, src->backing_str_len + 1); // memcpy because it may have several \0's in it

        // Rebase the transform type string
        if (src->transform_type_str)
            dst->transform_type_str = src->transform_type_str - src->backing_str + dst->backing_str;
        else
            dst->transform_type_str = NULL;

        // Rebase all the parameters, if present
        if (src->params) {
            int i;
            dst->param_count = src->param_count;
            MALLOC_ARRAY(dst->params, dst->param_count);

            for (i = 0; i < dst->param_count; i++) {
                const struct adios_transform_spec_kv_pair *src_kv = &src->params[i];
                struct adios_transform_spec_kv_pair *dst_kv = &dst->params[i];

                if (src_kv->key)
                    dst_kv->key = src_kv->key - src->backing_str + dst->backing_str;
                else
                    dst_kv->key = NULL;

                if (src_kv->value)
                    dst_kv->value = src_kv->value - src->backing_str + dst->backing_str;
                else
                    dst_kv->value = NULL;
            }
        } else {
            dst->params = NULL;
        }
    } else {
        dst->backing_str = NULL;
    }

    return dst;
}

#define FREE(x) {if(x)free(x);(x)=NULL;}
void adios_transform_free_spec(struct adios_transform_spec **spec_ptr) {
    struct adios_transform_spec *spec = *spec_ptr;
    FREE(spec->params);
    FREE(spec->backing_str);
    FREE(*spec_ptr)
}
#undef FREE
