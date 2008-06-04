#include "hw-utils.h"

int main (int argc, char ** argv)
{
    if (argc < 2)
    {
        printf ("\nUSAGE: bp2h5 XXXX.bp (XXXX.h5)\n\n");
        return -1;
    }
    else if (argc < 3)
    {
        hw_makeh5 (argv [1],NULL);
        return 0;
    }
    else if (argc < 4)
    {
        hw_makeh5 (argv [1], argv [2]);
        return 0;
    }

}
