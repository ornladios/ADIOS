Example write_vars/query_vars demonstrate how to create, 
evaluate and use a query over multiple queries.

Query: 80 < P < 90 AND V <= 50 
Read T where the query is true.

1. Generate the example data file
==================================

$ ./write_vars
$ bpls -ld vars.bp

  double   T     {5, 6} = 1.1 / 5.6 / 3.35 / 1.42449 
    (0,0)    1.1 1.2 1.3 1.4 1.5 1.6 
    (1,0)    2.1 2.2 2.3 2.4 2.5 2.6 
    (2,0)    3.1 3.2 3.3 3.4 3.5 3.6 
    (3,0)    4.1 4.2 4.3 4.4 4.5 4.6 
    (4,0)    5.1 5.2 5.3 5.4 5.5 5.6 

  double   P     {5, 6} = 31.6 / 95.5 / 66.6833 / 22.1717 
    (0,0)    41.1 61.2 81.3 81.4 91.5 31.6 
    (1,0)    42.1 62.2 82.3 82.4 92.5 32.6 
    (2,0)    43.1 63.2 83.3 83.4 93.5 33.6 
    (3,0)    44.1 64.2 84.3 84.4 94.5 34.6 
    (4,0)    45.1 65.2 85.3 85.4 95.5 35.6 

  double   V     {5, 6} = 41.1 / 55.6 / 49.15 / 5.30935 
    (0,0)    41.1 41.2 41.3 41.4 41.5 41.6 
    (1,0)    45.1 45.2 45.3 45.4 45.5 45.6 
    (2,0)    49.1 49.2 49.3 49.4 49.5 49.6 
    (3,0)    54.1 54.2 54.3 54.4 54.5 54.6 
    (4,0)    55.1 55.2 55.3 55.4 55.5 55.6 

2. If ADIOS was built with FastBit, create the index file (optional)
====================================================================

$ adios_index_fastbit vars.bp
$ ls -l vars.*
-rw-rw-r-- 1 adios adios    2253 Dec 14 10:52 vars.bp
-rw-rw-r-- 1 adios adios 1051182 Dec 14 10:52 vars.idx

Do not worry, the actual size of vars.idx is small:

$ du -sh vars.*
4.0K	vars.bp
8.0K	vars.idx


3. Run the query
================

The query and the data is chosen so that:
 - the query part of P is true for the two middle columns of P, 
 - the query part of V is true for the first three rows of V, 
 - and the combined query is true for 6 points 
   (the intersection of the two parts).

$ ./query_vars 

====== querying with one bound box for all variables =======
File : vars.bp
Query: (((P > 80.0) and (P < 90.0)) and (V <= 50.0))
Number of hits returned in batch 1 = 6 

Hit           i       j        P      V     T
----------------------------------------------
    0         0       2      81.3   41.3   1.3
    1         0       3      81.4   41.4   1.4
    2         1       2      82.3   45.3   2.3
    3         1       3      82.4   45.4   2.4
    4         2       2      83.3   49.3   3.3
    5         2       3      83.4   49.4   3.4







