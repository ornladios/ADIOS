#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

void makedirectory_ (char *dirname){
  mkdir(dirname,0755);
}

void unlink_file_ (char *filename){
  unlink(filename);
}

