/*
 * transforms_specparse.c
 *
 * Tests the "specparse" functionality, which parses the string passed as transform="..." into
 * the transform ID and a list of key-value pairs.
 *
 *  Created on: Jul 25, 2013
 *      Author: David A. Boyuka II
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "core/transforms/adios_transforms_specparse.h"

void test1() {
    struct adios_transform_spec *origspec = adios_transform_parse_spec("identity:a=123,b,c=321,,,f=12321");
    struct adios_transform_spec *origspecToClear = origspec;

    struct adios_transform_spec *spec = adios_transform_spec_copy(origspec);
    adios_transform_free_spec(&origspec);
    memset(origspecToClear, 0, sizeof(*origspecToClear));

    assert(spec->transform_type == adios_transform_identity);
    assert(strcmp(spec->transform_type_str, "identity") == 0);

    assert(spec->param_count == 6);
    assert(spec->params);

    assert(strcmp(spec->params[0].key, "a") == 0);
    assert(strcmp(spec->params[0].value, "123") == 0);

    assert(strcmp(spec->params[1].key, "b") == 0);
    assert(spec->params[1].value == NULL);

    assert(strcmp(spec->params[2].key, "c") == 0);
    assert(strcmp(spec->params[2].value, "321") == 0);

    assert(strcmp(spec->params[3].key, "") == 0);
    assert(spec->params[3].value == NULL);

    assert(strcmp(spec->params[4].key, "") == 0);
    assert(spec->params[4].value == NULL);

    assert(strcmp(spec->params[5].key, "f") == 0);
    assert(strcmp(spec->params[5].value, "12321") == 0);
}

int main(int argc, char **argv) {
    test1();
}


