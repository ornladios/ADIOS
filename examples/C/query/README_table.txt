Example write_table/query_table demonstrates how to create/evaluate and use
a query on different columns of a 2D table variable.

Query: Get the Kinetic Energy of Carbon atoms whose Potential <= 96
i.e.  "column 1 of A is 0" AND "column 2 < 96"

1. Generate the example data set
=================================

$ ./write_table
$ bpls -l table.bp -d A -n 7 -f "%4d "
  integer  A         {10, 7} = -1 / 201 / 18.9714 / 40.7168 
    (0,0)       0    0   97   15    8    0    7 
    (1,0)       1    0   96   16    5   -1    7 
    (2,0)       2    1   32    1    7    1    7 
    (3,0)       3    2  200    8    6    0    7 
    (4,0)       4    0   94   14    1    3    5 
    (5,0)       5    0   96   13    2    4    5 
    (6,0)       6    2  201    2    3    3    5 
    (7,0)       7    1   37    9    4    2    5 
    (8,0)       8    0   98   15    5    1    5 
    (9,0)       9    0   99   14    6    2    5 


just FYI:

$ bpls -l table.bp -d -S Elements Columns
  byte     Elements  {3, 9} = 0 / 121 / 81.3704 / 40.0131 
    (0,0)    "Carbon  "
    (1,0)    "Nitrogen"
    (2,0)    "Oxygen  "

  byte     Columns   {7, 11} = 0 / 116 / 77.2597 / 39.3811 
    (0, 0)    "ID        "
    (1, 0)    "Element   "
    (2, 0)    "Potential "
    (3, 0)    "Kinetic E "
    (4, 0)    "Position X"
    (5, 0)    "Position Y"
    (6, 0)    "Position Z"


2. If ADIOS was built with FastBit, create the index file (optional)
====================================================================

$ adios_index_fastbit table.bp

$ ls -l table.*
-rw-rw-r-- 1 adios adios    1879 Dec 14 12:08 table.bp
-rw-rw-r-- 1 adios adios 1050846 Dec 14 12:14 table.idx

$ du -sh  table.*
4.0K	table.bp
8.0K	table.idx


3. Evaluate the query
=====================

$ ./query_table 
File : table.bp
Variable A has 10 rows and 7 columns

====== querying over columns of a table  =======
Query: ((A = 0) and (A <= 96))
Number of hits returned in batch 1 = 3 

Hit           i       j    Kinetic E
----------------------------------------------
    0         1       3      16
    1         4       3      14
    2         5       3      13


Read all rows that match the query:

ID         Element    Potential  Kinetic E  Position X Position Y Position Z 
----------------------------------------------------------------------------
    1          0         96         16          5         -1          7      
    4          0         94         14          1          3          5      
    5          0         96         13          2          4          5 


