cccccc   Edit arrays

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      iflprt: a print parameter regulating the print file created when
!       the calculation is terminated (ncycle = ncymax)
!
!          iflprt = 0:  no printout at termination
!          iflprt = 1:  configuration printout at termination
!          iflprt = 2:  full printout except editn at termination
!          iflprt = 3:  full printout except editng at termination
!          iflprt = 4:  full printout at termination
!-----------------------------------------------------------------------

iflprt                       3                                          iflprt

!-----------------------------------------------------------------------
!      intprt: the number of cycles between the closing of an old print file and the opening of a 
!       new print file. The print files are named modelxxx.d, where xxx is the right justified 
!       consecutive number of the print file, e.g., the fifth print file is named 'model005.d'.
!
!      nprt: a print counter giving the number of cycles since the last closing of a print file.
!-----------------------------------------------------------------------

intprt                     200                                          intprt
intprt                       0                                          nprt

!-----------------------------------------------------------------------
!      iprtbgn: a print parameter regulating the print file created when the calculation is initiated
!       (nrst = 0) or restarted (nrst ne 0);
!
!          iprtbgn = 0, nrst = 0 : abbreviated edit
!          iprtbgn = 0, nrst > 0 : configuration edit
!          iprtbgn = 1, nrst = 0 : abbreviated edit
!          iprtbgn = 1, nrst > 0 : abbreviated edit
!          iprtbgn = 2, nrst = 0 : full edit
!          iprtbgn = 2, nrst > 0 : configuration edit
!          iprtbgn = 3, nrst = 0 : full edit
!          iprtbgn = 3, nrst > 0 : abreviated edit
!          iprtbgn = 4, nrst = 0 : full edit
!          iprtbgn = 4, nrst > 0 : full edit
!-----------------------------------------------------------------------

intprt                       1                                          iprtbgn

!-----------------------------------------------------------------------
!        Obsolete parameter.
!-----------------------------------------------------------------------

noutfl                      19                                          noutfl

!-----------------------------------------------------------------------
!          0: print to a print file is bypassed.
!          1: print to a print file is implemented.
!-----------------------------------------------------------------------

output                       1                                          iprint

!-----------------------------------------------------------------------
!      nrstd1: the unit number for file 'rstdmp1.d' to which temporary restart dumps are written
!       every intrst cycles. These restart dumps are temporary in the sense that the next restart
!       dump written to file 'rstdmp1.d' writes over the preceding dump.
!
!      nrstd2: the unit number for file 'rstdmp2.d' to which temporary restart dumps can be written
!       every other intrst cycles alternating with file rstdmp1.d.
!-----------------------------------------------------------------------

output                      19                                          nrstd1
output                      21                                          nrstd2

!-----------------------------------------------------------------------
!      noutpmt: the unit number of the permanent restart dump file (typicallyfile 'restartxxxxxxx.d',
!       where xxxxxxx is the cycle number right-justified).
!-----------------------------------------------------------------------

output                      18                                          noutpmt

!-----------------------------------------------------------------------
!      nouttmp: the unit number of the temporary restart dump file (typically file 'rstdmp1.d' or
!       'rstdmp2.d').
!-----------------------------------------------------------------------

output                      19                                          nouttmp

!-----------------------------------------------------------------------
!          0: print to a plot file is bypassed.
!          1: print to a plot file is implemented.
!-----------------------------------------------------------------------

plot                         0                                          iplot

!-----------------------------------------------------------------------
!      nplotc: the unit number for downloading model configuration plot data.
!      nplote: the unit number for downloading electron neutrino plot data.
!      nplota: the unit number for downloading electron antineutrino plot data.
!      nplott: the unit number for downloading muon and tau neutrino plot data.
!-----------------------------------------------------------------------

plot                        25                                          nplotc
plot                        26                                          nplote
plot                        27                                          nplota
plot                        28                                          nplott

!-----------------------------------------------------------------------
!      intplf: the number of cycles between the printing to a plot file. 
!-----------------------------------------------------------------------

plot                   9000000                                          intplf

!-----------------------------------------------------------------------
!      npltf: the number of cycles since the last printing to a plot file.
!-----------------------------------------------------------------------

plot                         0                                          npltf

!-----------------------------------------------------------------------
!      nmodel: the number of the current print file.
!-----------------------------------------------------------------------

nmodel                       0                                          nmodel

!-----------------------------------------------------------------------
!      rhoprt (read in as rhoprint(i)): the ith central density at which a print file is created. 
!       The print file is named modeldxx.d, where xx is the value of i right justified.
!-----------------------------------------------------------------------

rhoprt             1      1.00000000E+11 
rhoprt             2      1.00000000E+12
rhoprt             3      1.00000000E+13
rhoprt             4      1.00000000E+14
rhoprt             5      4.26500000E+18
rhoprt             6      2.83400000E+18
rhoprt             7      1.39300000E+18
rhoprt             8      3.12300000E+18
rhoprt             9      6.00000000e+19

!-----------------------------------------------------------------------
!         Parameters for configuration edits
!
!      intedc(i): the number of cycles between the implementation of subsection i of subroutine editc. 
!
!      nedc(i): the number of cycles since the last implementation of subsection i of subroutine editc.
!
!      idxedc(i): a parameter used in subroutine editc. Subsection i of subroutine editc will print
!       data for every idxedc(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxedc(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intedc             1        50                                          intedc
intedc             1         0                                          nedc
intedc             1         1                                          idxedc

!-----------------------------------------------------------------------
!         Parameters for kinetic, internal, and gravitational energy edits
!
!      intede(i): the number of cycles between the implementation of subsection i of subroutine edit_e. 
!
!      nede(i): the number of cycles since the last implementation of subsection i of subroutine edit_e.
!
!      idxede(i): a parameter used in subroutine edit_e. Subsection i of subroutine editc will print
!       data for every idxede(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxede(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intede             1       200                                          intede
intede             1         0                                          nede
intede             1         1                                          idxede

!-----------------------------------------------------------------------
!        Parameters for editmi edits
!
!      intdmi(i): the number of cycles between the implementation of subsection i of subroutine editmi.
!
!      nedmi(i): the number of cycles since the last implementation of subsection i of subroutine editmi.
!
!      idxemi(i): a parameter used in subroutine editmi. Subsection i of subroutine editmi will print
!       data for every idxemi(i)'th radial zone, starting with the outermost zone. The data 
!       corresponding to the innermost radial zone will also be printed. Setting idxemi(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intdmi             1      1000                                          intdmi
intdmi             1         0                                          nedmi
intdmi             1         1                                          idxemi

!-----------------------------------------------------------------------
!         Parameters for mass average edits
!
!      intdma(i): the number of cycles between the implementation of subsection i of subroutine editma.
!
!      nedma(i): the number of cycles since the last implementation of subsection i of subroutine editma.
!
!      idxema(i): a parameter used in subroutine editma. Subsection i of subroutine editma will print
!       data for every idxema(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxema(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intdma             1      1000                                          intdma
intdma             1         0                                          nedma
intdma             1         1                                          idxema
intdma             2      1000                                          intdma
intdma             2         0                                          nedma
intdma             2         1                                          idxema
intdma             3      1000                                          intdma
intdma             3         0                                          nedma
intdma             3         1                                          idxema

!-----------------------------------------------------------------------
!         Parameters for hydro edits
!
!      intedh(i): the number of cycles between the implementation of subsection i of subroutine edith.
!
!      nedh(i): the number of cycles since the last implementation of subsection i of subroutine edith.
!
!      idxedh(i): a parameter used in subroutine edith. Subsection i of subroutine editma will print
!       data for every idxedh(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxedh(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intedh             1      1000                                          intedh
intedh             1         0                                          nedh
intedh             1         1                                          idxedh
intedh             2      1000                                          intedh
intedh             2         0                                          nedh
intedh             2         1                                          idxedh
intedh             3      1000                                          intedh
intedh             3         0                                          nedh
intedh             3         1                                          idxedh
intedh             4      1000                                          intedh
intedh             4         0                                          nedh
intedh             4         1                                          idxedh
intedh             5      1000                                          intedh
intedh             5         0                                          nedh
intedh             5         1                                          idxedh
intedh             6      1000                                          intedh
intedh             6         0                                          nedh
intedh             6         1                                          idxedh
intedh             7      1000                                          intedh
intedh             7         0                                          nedh
intedh             7         1                                          idxedh
intedh             8      1000                                          intedh
intedh             8         0                                          nedh
intedh             8         1                                          idxedh
intedh             9      1000                                          intedh
intedh             9         0                                          nedh
intedh             9         1                                          idxedh

!-----------------------------------------------------------------------
!         Parameters for pressure - stress edits
!
!      intdps(i): the number of cycles between the implementation of subsection i of subroutine editps.
!
!      nedps(i): the number of cycles since the last implementation of subsection i of subroutine editps.
!
!      idxeps(i): a parameter used in subroutine editps. Subsection i of subroutine editps will print
!       data for every idxeps(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxeps(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intdps             1      1000                                          intdps
intdps             1         0                                          nedps
intdps             1         1                                          idxeps

!-----------------------------------------------------------------------
!         Parameters for energy edits
!
!      intedu(i): the number of cycles between the implementation of subsection i of subroutine editu.
!
!      nedu(i): the number of cycles since the last implementation of subsection i of subroutine editu.
!
!      idxedu(i): a parameter used in subroutine editu. Subsection i of subroutine editu will print
!       data for every idxedu(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxedu(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intedu             1      1000                                          intedu
intedu             1         0                                          nedu
intedu             1         1                                          idxedu
intedu             2      1000                                          intedu
intedu             2         0                                          nedu
intedu             2         1                                          idxedu
intedu             3      1000                                          intedu
intedu             3         0                                          nedu
intedu             3         1                                          idxedu
intedu             4      1000                                          intedu
intedu             4         0                                          nedu
intedu             4         1                                          idxedu
intedu             5      1000                                          intedu
intedu             5         0                                          nedu
intedu             5         1                                          idxedu
intedu             6      1000                                          intedu
intedu             6         0                                          nedu
intedu             6         1                                          idxedu
intedu             7      1000                                          intedu
intedu             7         0                                          nedu
intedu             7         1                                          idxedu
intedu             8      1000                                          intedu
intedu             8         0                                          nedu
intedu             8         1                                          idxedu
intedu             9      1000                                          intedu
intedu             9         0                                          nedu
intedu             9         1                                          idxedu
intedu            10      1000                                          intedu
intedu            10         0                                          nedu
intedu            10         1                                          idxedu
intedu            11      1000                                          intedu
intedu            11         0                                          nedu
intedu            11         1                                          idxedu
intedu            12      1000                                          intedu
intedu            12         0                                          nedu
intedu            12         1                                          idxedu
intedu            13      1000                                          intedu
intedu            13         0                                          nedu
intedu            13         1                                          idxedu
intedu            14      1000                                          intedu
intedu            14         0                                          nedu
intedu            14         1                                          idxedu

!-----------------------------------------------------------------------
!         Parameters for composition edits
!
!      intedy(i): the number of cycles between the implementation of subsection i of subroutine edity.
!
!      nedy(i): the number of cycles since the last implementation of subsection i of subroutine edity.
!
!      idxedy(i): a parameter used in subroutine edity. Subsection i of subroutine edity will print
!       data for every idxedy(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxedy(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intedy             1      1000                                          intedy
intedy             1         0                                          nedy
intedy             1         1                                          idxedy
intedy             2      1000                                          intedy
intedy             2         0                                          nedy
intedy             2         1                                          idxedy
intedy             3      1000                                          intedy
intedy             3         0                                          nedy
intedy             3         1                                          idxedy
intedy             4      1000                                          intedy
intedy             4         0                                          nedy
intedy             4         1                                          idxedy
intedy             5      1000                                          intedy
intedy             5         0                                          nedy
intedy             5         1                                          idxedy
intedy             6      1000                                          intedy
intedy             6         0                                          nedy
intedy             6         1                                          idxedy

!-----------------------------------------------------------------------
!         Parameters for entropy and chemical potential edits
!
!      intdsc(i): the number of cycles between the implementation of subsection i of subroutine editsc.
!
!      nedsc(i): the number of cycles since the last implementation of subsection i of subroutine editsc.
!
!      idxesc(i): a parameter used in subroutine editsc. Subsection i of subroutine editsc will print 
!       data for every idxesc(i)'th radial zone, starting with the outermost zone. The data
!       corresponding to the innermost radial zone will also be printed. Setting idxesc(i) = 1 will
!       cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

intdsc             1      1000                                          intdsc
intdsc             1         0                                          nedsc
intdsc             1         1                                          idxesc
intdsc             2      1000                                          intdsc
intdsc             2         0                                          nedsc
intdsc             2         1                                          idxesc

!-----------------------------------------------------------------------
!         Parameters for diferential neutrino edits
!
!      intedn(n): the number of cycles between the implementation of subroutine editn for neutrinos 
!       of type n.
!
!      nedn(n): the number of cycles since the last implementation of subroutine editn for neutrinos 
!       of type n.
!
!      idxedn(n): a parameter used in subroutine editn. Subroutine editn will print the data
!       corresponding neutrinos of type n for every idxedn(n)'th radial zone, starting with the
!       outermost zone. The data corresponding to the innermost radial zone will also be printed.
!       Setting idxedy(i) = 1 will cause data for all radial zones to be printed.
!
!      niedn(i): a parameter used in subroutine editn. Subroutine editn will print the i'th datum
!       corresponding to neutrinos of type n every niedn(i) implementation of editn.
!
!      neden(i): the number of implementations of editn since the last printing of the i'th datum 
!       corresponding to neutrinos of type n.
!-----------------------------------------------------------------------

intedn             1      1000                                          intedn
intedn             1         0                                          nedn
intedn             1        20                                          idxedn
intedn             2      1000                                          intedn
intedn             2         0                                          nedn
intedn             2        20                                          idxedn
intedn             3      1000                                          intedn
intedn             3         0                                          nedn
intedn             3        20                                          idxedn
intedn             4      1000                                          intedn
intedn             4         0                                          nedn
intedn             4        20                                          idxedn
intedn             1         1                                          niedn
intedn             2         1                                          niedn
intedn             3         1                                          niedn
intedn             4         1                                          niedn
intedn             5         1                                          niedn
intedn             6         1                                          niedn
intedn             7         1                                          niedn
intedn             8         1                                          niedn
intedn             9         1                                          niedn
intedn            10         1                                          niedn
intedn            11         1                                          niedn
intedn            12         1                                          niedn
intedn            13         1                                          niedn
intedn            14         1                                          niedn
intedn            15         1                                          niedn
intedn            16         1                                          niedn
intedn            17         1                                          niedn
intedn            18         1                                          niedn
intedn            19         1                                          niedn
intedn            20         1                                          niedn
intedn            21         1                                          niedn
intedn            22         1                                          niedn
intedn            23         1                                          niedn
intedn            24         1                                          niedn
intedn            25         1                                          niedn
intedn            26         1                                          niedn
intedn            27         1                                          niedn
intedn            28         1                                          niedn
intedn            29         1                                          niedn
intedn            30         1                                          niedn
intedn            31         1                                          niedn
intedn            32         1                                          niedn
intedn            33         1                                          niedn
intedn            34         1                                          niedn
intedn            35         1                                          niedn
intedn            36         1                                          niedn
intedn            37         1                                          niedn
intedn            38         1                                          niedn
intedn            39         1                                          niedn
intedn            40         1                                          niedn
intedn             1         0                                          neden
intedn             2         0                                          neden
intedn             3         0                                          neden
intedn             4         0                                          neden
intedn             5         0                                          neden
intedn             6         0                                          neden
intedn             7         0                                          neden
intedn             8         0                                          neden
intedn             9         0                                          neden
intedn            10         0                                          neden
intedn            11         0                                          neden
intedn            12         0                                          neden
intedn            13         0                                          neden
intedn            14         0                                          neden
intedn            15         0                                          neden
intedn            16         0                                          neden
intedn            17         0                                          neden
intedn            18         0                                          neden
intedn            19         0                                          neden
intedn            20         0                                          neden
intedn            21         0                                          neden
intedn            22         0                                          neden
intedn            23         0                                          neden
intedn            24         0                                          neden
intedn            25         0                                          neden
intedn            26         0                                          neden
intedn            27         0                                          neden
intedn            28         0                                          neden
intedn            29         0                                          neden
intedn            30         0                                          neden
intedn            31         0                                          neden
intedn            32         0                                          neden
intedn            33         0                                          neden
intedn            34         0                                          neden
intedn            35         0                                          neden
intedn            36         0                                          neden
intedn            37         0                                          neden
intedn            38         0                                          neden
intedn            39         0                                          neden
intedn            40         0                                          neden

!-----------------------------------------------------------------------
!         Parameters for integral neutrino edits
!
!      intdng(i,n): the number of cycles between the implementation of subsection i of subroutine 
!       editng for neutrinos of type n.
!
!      nedng(i,n): the number of cycles since the last implementation of subsection i of subroutine 
!       editng for neutrinos of type n.
!
!      idxeng(i,n): a parameter used in subroutine editng. Subsection i of subroutine editng will 
!       print the data corresponding to every idxesc(i)'th radial zone, starting with the outermost
!       zone. The data corresponding to the innermost radial zone will also be printed. Setting
!       idxeng(i,n) = 1 will cause data for all radial zones to be printed. 
!-----------------------------------------------------------------------

intdng             1      5000         1                                intdng
intdng             1         0         1                                nedng
intdng             1         1         1                                idxeng
intdng             1      5000         2                                intdng
intdng             1         0         2                                nedng
intdng             1         1         2                                idxeng
intdng             1      5000         3                                intdng
intdng             1         0         3                                nedng
intdng             1         1         3                                idxeng
intdng             1      5000         4                                intdng
intdng             1         0         4                                nedng
intdng             1         1         4                                idxeng
intdng             2      5000         1                                intdng
intdng             2         0         1                                nedng
intdng             2         1         1                                idxeng
intdng             2      5000         2                                intdng
intdng             2         0         2                                nedng
intdng             2         1         2                                idxeng
intdng             2      5000         3                                intdng
intdng             2         0         3                                nedng
intdng             2         1         3                                idxeng
intdng             2      5000         4                                intdng
intdng             2         0         4                                nedng
intdng             2         1         4                                idxeng
intdng             3      5000         1                                intdng
intdng             3         0         1                                nedng
intdng             3         1         1                                idxeng
intdng             3      5000         2                                intdng
intdng             3         0         2                                nedng
intdng             3         1         2                                idxeng
intdng             3      5000         3                                intdng
intdng             3         0         3                                nedng
intdng             3         1         3                                idxeng
intdng             3      5000         4                                intdng
intdng             3         0         4                                nedng
intdng             3         1         4                                idxeng
intdng             4      5000         1                                intdng
intdng             4         0         1                                nedng
intdng             4         1         1                                idxeng
intdng             4      5000         2                                intdng
intdng             4         0         2                                nedng
intdng             4         1         2                                idxeng
intdng             4      5000         3                                intdng
intdng             4         0         3                                nedng
intdng             4         1         3                                idxeng
intdng             4      5000         4                                intdng
intdng             4         0         4                                nedng
intdng             4         1         4                                idxeng
intdng             5      5000         1                                intdng
intdng             5         0         1                                nedng
intdng             5         1         1                                idxeng
intdng             5      5000         2                                intdng
intdng             5         0         2                                nedng
intdng             5         1         2                                idxeng
intdng             5      5000         3                                intdng
intdng             5         0         3                                nedng
intdng             5         1         3                                idxeng
intdng             5      5000         4                                intdng
intdng             5         0         4                                nedng
intdng             5         1         4                                idxeng
intdng             6      5000         1                                intdng
intdng             6         0         1                                nedng
intdng             6         1         1                                idxeng
intdng             6      5000         2                                intdng
intdng             6         0         2                                nedng
intdng             6         1         2                                idxeng
intdng             6      5000         3                                intdng
intdng             6         0         3                                nedng
intdng             6         1         3                                idxeng
intdng             6      5000         4                                intdng
intdng             6         0         4                                nedng
intdng             6         1         4                                idxeng
intdng             7      5000         1                                intdng
intdng             7         0         1                                nedng
intdng             7         1         1                                idxeng
intdng             7      5000         2                                intdng
intdng             7         0         2                                nedng
intdng             7         1         2                                idxeng
intdng             7      5000         3                                intdng
intdng             7         0         3                                nedng
intdng             7         1         3                                idxeng
intdng             7      5000         4                                intdng
intdng             7         0         4                                nedng
intdng             7         1         4                                idxeng
intdng             8      5000         1                                intdng
intdng             8         0         1                                nedng
intdng             8         1         1                                idxeng
intdng             8      5000         2                                intdng
intdng             8         0         2                                nedng
intdng             8         1         2                                idxeng
intdng             8      5000         3                                intdng
intdng             8         0         3                                nedng
intdng             8         1         3                                idxeng
intdng             8      5000         4                                intdng
intdng             8         0         4                                nedng
intdng             8         1         4                                idxeng
intdng             9      5000         1                                intdng
intdng             9         0         1                                nedng
intdng             9         1         1                                idxeng
intdng             9      5000         2                                intdng
intdng             9         0         2                                nedng
intdng             9         1         2                                idxeng
intdng             9      5000         3                                intdng
intdng             9         0         3                                nedng
intdng             9         1         3                                idxeng
intdng             9      5000         4                                intdng
intdng             9         0         4                                nedng
intdng             9         1         4                                idxeng
intdng            10      5000         1                                intdng
intdng            10         0         1                                nedng
intdng            10         1         1                                idxeng
intdng            10      5000         2                                intdng
intdng            10         0         2                                nedng
intdng            10         1         2                                idxeng
intdng            10      5000         3                                intdng
intdng            10         0         3                                nedng
intdng            10         1         3                                idxeng
intdng            10      5000         4                                intdng
intdng            10         0         4                                nedng
intdng            10         1         4                                idxeng
intdng            11      5000         1                                intdng
intdng            11         0         1                                nedng
intdng            11         1         1                                idxeng
intdng            11      5000         2                                intdng
intdng            11         0         2                                nedng
intdng            11         1         2                                idxeng
intdng            11      5000         3                                intdng
intdng            11         0         3                                nedng
intdng            11         1         3                                idxeng
intdng            11      5000         4                                intdng
intdng            11         0         4                                nedng
intdng            11         1         4                                idxeng
intdng            12      5000         1                                intdng
intdng            12         0         1                                nedng
intdng            12         1         1                                idxeng
intdng            12      5000         2                                intdng
intdng            12         0         2                                nedng
intdng            12         1         2                                idxeng
intdng            12      5000         3                                intdng
intdng            12         0         3                                nedng
intdng            12         1         3                                idxeng
intdng            12      5000         4                                intdng
intdng            12         0         4                                nedng
intdng            12         1         4                                idxeng
intdng            13      5000         1                                intdng
intdng            13         0         1                                nedng
intdng            13         1         1                                idxeng
intdng            13      5000         2                                intdng
intdng            13         0         2                                nedng
intdng            13         1         2                                idxeng
intdng            13      5000         3                                intdng
intdng            13         0         3                                nedng
intdng            13         1         3                                idxeng
intdng            13      5000         4                                intdng
intdng            13         0         4                                nedng
intdng            13         1         4                                idxeng
intdng            14      5000         1                                intdng
intdng            14         0         1                                nedng
intdng            14         1         1                                idxeng
intdng            14      5000         2                                intdng
intdng            14         0         2                                nedng
intdng            14         1         2                                idxeng
intdng            14      5000         3                                intdng
intdng            14         0         3                                nedng
intdng            14         1         3                                idxeng
intdng            14      5000         4                                intdng
intdng            14         0         4                                nedng
intdng            14         1         4                                idxeng
intdng            15      5000         1                                intdng
intdng            15         0         1                                nedng
intdng            15         1         1                                idxeng
intdng            15      5000         2                                intdng
intdng            15         0         2                                nedng
intdng            15         1         2                                idxeng
intdng            15      5000         3                                intdng
intdng            15         0         3                                nedng
intdng            15         1         3                                idxeng
intdng            15      5000         4                                intdng
intdng            15         0         4                                nedng
intdng            15         1         4                                idxeng
intdng            16      5000         1                                intdng
intdng            16         0         1                                nedng
intdng            16         1         1                                idxeng
intdng            16      5000         2                                intdng
intdng            16         0         2                                nedng
intdng            16         1         2                                idxeng
intdng            16      5000         3                                intdng
intdng            16         0         3                                nedng
intdng            16         1         3                                idxeng
intdng            16      5000         4                                intdng
intdng            16         0         4                                nedng
intdng            16         1         4                                idxeng
intdng            17      5000         1                                intdng
intdng            17         0         1                                nedng
intdng            17         1         1                                idxeng
intdng            17      5000         2                                intdng
intdng            17         0         2                                nedng
intdng            17         1         2                                idxeng
intdng            17      5000         3                                intdng
intdng            17         0         3                                nedng
intdng            17         1         3                                idxeng
intdng            17      5000         4                                intdng
intdng            17         0         4                                nedng
intdng            17         1         4                                idxeng
intdng            18      5000         1                                intdng
intdng            18         0         1                                nedng
intdng            18         1         1                                idxeng
intdng            18      5000         2                                intdng
intdng            18         0         2                                nedng
intdng            18         1         2                                idxeng
intdng            18      5000         3                                intdng
intdng            18         0         3                                nedng
intdng            18         1         3                                idxeng
intdng            18      5000         4                                intdng
intdng            18         0         4                                nedng
intdng            18         1         4                                idxeng

cccccc   Restart file arrays
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!         Parameters for temperary restart dumps
!
!      intrst: the number of cycles between temporary restart dumps.
!
!      nnrst: the number of cycles since the last temporary restart dump.
!
!      nrstd1: the unit number for file 'rstdmp1.d' to which temporary restart dumps are written
!       every intrst cycles. These restart dumps are temporary in the sense that the next restart
!       dump written to file 'rstdmp1.d' writes over the preceding dump.
!
!      nrstd2: the unit number for file 'rstdmp2.d' to which temporary restart dumps can be written
!       every other intrst cycles alternating with file rstdmp1.d.
!
!      nrstfl: a temporary restart dump parameter -
!
!          2: restart dumps are alternated between file 'rstdmp1.d' and file 'rstdmp2.d'.
!          ne 2: restart dumps are written only to file 'rstdmp1.d'.
!
!      nouttmp: the unit number of the temporary restart dump file (typically file 'rstdmp1.d' or
!       'rstdmp2.d').
!-----------------------------------------------------------------------

intrst                     100                                          intrst

intrst                       0                                          nnrst

intrst                       2                                          nrstfl

!-----------------------------------------------------------------------
!         Parameters for permanent restart dumps
!
!      noutpmt: the unit number of the permanent restart dump file (typically
!       file 'restartxxxxxxx.d', where xxxxxxx is the cycle number right-justified).
!
!      intprm: the number of cycles between permanent restart dumps.
!
!      nprm: the number of cycles since the last permanent restart dump.
!
!      ncychg: at termination, write to restart file with
!          ncymax = ncymax + ncychg.
!
!      irstbgn: restart dump parameter -
!
!          1: temporary and permanent restart dumps are implemented at the
!              beginning of the run (ncycle = 0)
!          ne 1: emporary and permanent restart dumps are bypassed at the
!              beginning of the run(ncycle = 0)
!
!      ncyrst(i): the cycle number at which the ith prescribed permanent
!              restart dump is implemented.
!-----------------------------------------------------------------------

intprm                 1000000                                          intprm

intprm                       0                                          nprm

intrst                     600                                          ncychg

intrst                       0                                          irstbgn

ncyrst             1      2000                                          ncyrst
ncyrst             2      3000                                          ncyrst
ncyrst             3      4000                                          ncyrst
ncyrst             4      5000                                          ncyrst
ncyrst             5      6000                                          ncyrst
ncyrst             6      7000                                          ncyrst
ncyrst             7      8000                                          ncyrst
ncyrst             8      9000                                          ncyrst
ncyrst             9     10000                                          ncyrst
ncyrst            10     20000                                          ncyrst
ncyrst            11     30000                                          ncyrst
ncyrst            12     40000                                          ncyrst
ncyrst            13     50000                                          ncyrst
ncyrst            14     60000                                          ncyrst
ncyrst            15     70000                                          ncyrst
ncyrst            16     80000                                          ncyrst
ncyrst            17     90000                                          ncyrst
ncyrst            18    100000                                          ncyrst
ncyrst            19    120000                                          ncyrst
ncyrst            20    140000                                          ncyrst
ncyrst            21    160000                                          ncyrst
ncyrst            22    180000                                          ncyrst
ncyrst            23    200000                                          ncyrst
ncyrst            24    220000                                          ncyrst
ncyrst            25    240000                                          ncyrst
ncyrst            26    260000                                          ncyrst
ncyrst            27    280000                                          ncyrst
ncyrst            28    300000                                          ncyrst
ncyrst            29    320000                                          ncyrst
ncyrst            30    340000                                          ncyrst
ncyrst            31    360000                                          ncyrst
ncyrst            32    380000                                          ncyrst
ncyrst            33    400000                                          ncyrst
ncyrst            34    420000                                          ncyrst
ncyrst            35    440000                                          ncyrst
ncyrst            36    460000                                          ncyrst
ncyrst            37    480000                                          ncyrst
ncyrst            38    500000                                          ncyrst
ncyrst            39    550000                                          ncyrst
ncyrst            40    600000                                          ncyrst
ncyrst            41    650000                                          ncyrst
ncyrst            42    700000                                          ncyrst
ncyrst            43    750000                                          ncyrst
ncyrst            44    800000                                          ncyrst
ncyrst            45    850000                                          ncyrst
ncyrst            46    900000                                          ncyrst
ncyrst            47    950000                                          ncyrst
ncyrst            48   1000000                                          ncyrst

cccccc   Plot arrays
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!         Parameters for neutrino luminosity plot edits
!
!      ilumplt: luminosity plot edit switch
!
!          0: bypass subroutine lumplot.
!          1: enter subroutine lumplot.
!
!      nlumplt1: the unit number for editing lumplot data at r_lum.
!
!      nlumplt2: the unit number for editing lumplot data at d_lum.
!
!      nlum: number of cycles since the last implemtation of subroutine lumplot.
!
!      r_lum: the radius at which to evaluate luminosity data.
!
!      d_lum: the density at which to evaluate luminosity data.
!
!      intlum, ncylum:
!
!          ncycle < ncylum(1)             : subroutine lumplot is called every intlum(1) cycles.
!          ncylum(1) < ncycle < ncylum(2) : subroutine lumplot is called every intlum(2) cycles.
!          ncylum(2) < ncycle             : subroutine lumplot is called every intlum(3) cycles.
!-----------------------------------------------------------------------

lumplt                       0                                          ilumplt
lumplt            31        32                                          nlumplt
lumplt                              1.00000000e+08                      r_lum
lumplt                              1.00000000e+12                      d_lum
lumplt                       0                                          nlum
lumplt             1    900000                                          intlum
lumplt             2    900000                                          intlum
lumplt             3    900000                                          intlum
lumplt             1      2000                                          ncylum
lumplt             2     10000                                          ncylum

!-----------------------------------------------------------------------
!         Parameters for rms neutrino energy plot edits
!
!      ienuplt: rms neutrino energy plot edit switch
!
!          0: bypass subroutine enuvplot.
!          1: enter subroutine enuvplot.
!
!      nenuplt: the unit number for editing enuvplot data.
!
!      nenu: the number of cycles since the last implemtation of subroutine enuvplot.
!
!      r_e_rms: the radius at which to store neutrino energy data.
!
!      d_e_rms: the density at which to store neutrino energy data.
!
!      intenu, ncyenu:
!
!          ncycle < ncyenu(1)             : subroutine enuvplot is called every intenu(1) cycles.
!          ncyenu(1) < ncycle < ncyenu(2) : subroutine enuvplot is called every intenu(2) cycles.
!          ncyenu(2) < ncycle             : subroutine enuvplot is called every intenu(3) cycles.
!-----------------------------------------------------------------------

enuplt                       0                                          ienuplt
enuplt            35        36                                          nenuplt
enuplt                              1.00000000e+08                      r_e_rms
enuplt                              1.00000000e+12                      d_e_rms
enuplt                       0                                          nenu
enuplt             1    900000                                          intenu
enuplt             2    900000                                          intenu
enuplt             3    900000                                          intenu
enuplt             1      2000                                          ncyenu
enuplt             2     10000                                          ncyenu 

!-----------------------------------------------------------------------
!         Parameters for selected radii plot edits
!
!      irnuplt: selected radii plot edit switch
!
!          0: bypass subroutine rnuplot.
!          1: enter subroutine rnuplot.
!
!      nrnuplt: the unit number for editing rnuplot data.
!
!      nrum: the number of cycles since the last implemtation of subroutine rnuplot.
!
!      intrnu, ncyrnu:
!
!          ncycle < ncyrnu(1)             : subroutine rnuplot is called every intrnu(1) cycles.
!          ncyrnu(1) < ncycle < ncyrnu(2) : subroutine rnuplot is called every intrnu(2) cycles.
!          ncyrnu(2) < ncycle             : subroutine rnuplot is called every intrnu(3) cycles.
!-----------------------------------------------------------------------

rnuplt                       0                                          irnuplt
rnuplt                      37                                          nrnuplt
rnuplt                       0                                          nrnu
rnuplt             1    900000                                          intrnu
rnuplt             2    900000                                          intrnu
rnuplt             3    900000                                          intrnu
rnuplt             1      2000                                          ncyrnu
rnuplt             2     10000                                          ncyrnu 

!-----------------------------------------------------------------------
!         Parameters for bounday plot edits
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      dtimeplot: write to subroutine bnuplot every dtimeplot ms.
!-----------------------------------------------------------------------

bnuplt                              1.00000000e+00                      dtimeplot

!-----------------------------------------------------------------------
!      rinnerb: quantities are interpolated to rinnerb for ibound data.
!
!      nplotinnerb: the unit number for ibound data.
!
!      iplotinnerb: inner boundary plot switch
!
!          0: bypass writing ibound data.
!          ne 0: write ibound data.
!-----------------------------------------------------------------------

bnuplt            12         0      2.00000000e+06                      innerb

!-----------------------------------------------------------------------
!      routerb: quantities are interpolated to routerb for obound data.
!
!      nplotouterb: the unit number for obound data.
!
!      iplotouterb: outer boundary plot switch
!
!          0: bypass writing obound data.
!          ne 0: write obound data.
!-----------------------------------------------------------------------

bnuplt            13         0      1.00000000e+08                      outerb

!-----------------------------------------------------------------------
!      r_lumerms: luminosities and rms energies are interpolated to r_lumerms for lbound data.
!
!      nplotlum: the unit number for lbound data.
!
!      iplotlum: luminosities and rms energies plot switch
!
!          0: bypass writing lbound data.
!          ne 0: write lbound data.
!-----------------------------------------------------------------------

bnuplt            14         0      1.00000000e+08                      lumermsplt

!-----------------------------------------------------------------------
!      nplotshk: the unit number for shock data.
!
!      iplotshk: shock data plot switch
!
!          0: bypass writing shock data.
!          ne 0: write shock data.
!-----------------------------------------------------------------------

bnuplt            15         0                                          shkplt

!-----------------------------------------------------------------------
!      nplotcnv: the unit number for convection data.
!
!      iplotcnv: convection plot switch
!
!          0: bypass writing convection data.
!          ne 0: write convection data.
!-----------------------------------------------------------------------

bnuplt            16         0                                          cnvplt

!-----------------------------------------------------------------------
!         Parameters for configuration plot edits
!
!      ivarplt: configuration plot edit switch
!
!          0: bypass subroutine varplot.
!          1: enter subroutine varplot.
!
!      nvar: the number of cycles since the last implemtation of subroutine varplot.
!
!      nvarint: subroutine varplot is called every nvarint cycles.
!
!      nvarplt: the unit number for editing varplot data.
!
!      nvardump - varplot file number (varplot files are numbered sequentially.
!
!      dtvarplot : write to subroutine varplot every dtvarplot ms.
!-----------------------------------------------------------------------

varplt            17         0                                          ivarplt
varplt             0      1000         0                                nvarint
varplt                              1.00000000E+06                      dtvarplot

!-----------------------------------------------------------------------
!         Parameters for comparison plot edits
!
!      icomplt: comparison plot edit switch
!
!          0: bypass subroutine complot.
!          1: enter subroutine complot.
!
!      ncomplt: the unit number for editing complot data.
!
!      ncomdump: the complot file number (complot files are numbered sequentially.
!
!      dtcomplot: write to subroutine complot every dtcomplot ms.
!-----------------------------------------------------------------------

complt            17         0                                          icomplt
complt                       0      1.00000000E+06                      dtcomplot

!-----------------------------------------------------------------------
!         Parameters for lagrangian plot edits
!
!      nlagplt: the unit number for lagplot data.
!
!      ilagplt: lagrangian plot switch
!
!          0: bypass writing lagplot data.
!          ne 0: write nuplot data.
!
!      nlagdump: lagplot file number (lagplot files are numbered sequentially.
!
!      msslag: lagrangian mass of fluid element for lagplot data.
!-----------------------------------------------------------------------

lagplt            18         0         0                                ilagplt
lagplt                              0.00000000E+00                      msslag

!-----------------------------------------------------------------------
!         Parameters for r-lagrangian plot edits
!
!      nrlagplt: the unit number for rlagplot data.
!
!      irlagplt: r-lagrangian plot switch
!
!          0: bypass writing rlagplot data.
!          ne 0: write rlagplot data.
!
!      dmlag: the lagrangian mass difference between adjacent lagrangian
!       points to be plotted.
!-----------------------------------------------------------------------

rlagpt            31         0                                          irlagplt
rlagpt                              0.02500000E+00                      dmlag

!-----------------------------------------------------------------------
!         Parameters for energy-check plot edits
!
!      n_eplt: the unit number for e_chk data.
!
!      i_eplt: energy-check plot switch
!
!          0: bypass writing e_chk data.
!          ne 0: write e_chk data.
!
!      dt_eplot - write to file e_chk.d every dt_eplot ms.
!-----------------------------------------------------------------------

e_plt             17         1                                          i_eplt
e_plt                               5.00000000e+00                      dt_eplot

!-----------------------------------------------------------------------
!         Parameters for nurad plot edits
!
!      nplotnurad; the unit number for nuradplot data.
!
!      iplotnurad: nurad plot switch
!
!          0: bypass writing nuradplot data.
!          ne 0: write nuradplot data.
!
!      dtnuradplot: write to subroutine nuradplot every dtnuradplot ms.
!
!      r_nurad: the radius at which to evaluate neutrino radiation.
!
!      rho_nurad: the density at which to evaluate neutrino radiation.
!
!      nu_r(k,n): the number of neutrinos of energy group k radiated across r_nurad in time 
!       dtnuradplot.
!
!      nu_rt(k,n): the cumulative number of neutrinos of energy group k radiated across r_nurad.
!
!      nu_rho(k,n): the number of neutrinos of energy group k radiated across rho_nurad in time 
!       dtnuradplot.
!
!      nu_rhot(k,n): the cumulative number of neutrinos of energy group k radiated across rho_nurad.
!
!      unu_r(k,n): the energy of neutrinos of energy group k radiated across r_nurad in time
!       dtnuradplot.
!
!      unu_rt(k,n): the- cumulative energy of neutrinos of energy group k radiated across r_nurad.
!
!      unu_rho(k,n): the energy of neutrinos of energy group k radiated across rho_nurad in time
!       dtnuradplot.
!
!      unu_rhot(k,n); the cumulative energy of neutrinos of energy group k radiated across rho_nurad.
!-----------------------------------------------------------------------

nurad             18         0                                          inuplt
nurad                               1.00000000E+02                      dtnuradplot
nurad                               1.00000000E+08                      r_nurad
nurad                               1.00000000E+12                      rho_nurad

!-----------------------------------------------------------------------
!         Parameters for neutrino distribution plot edits
!
!      nnudata : the unit number for nuplot data.
!
!      inudata : neutrino distribution plot switch
!
!          0: bypass writing nuplot data.
!          ne 0: write nuplot data.
!
!      dtnudata : write to subroutine nuplot every dtnudata ms.
!
!      r_nudata : the radius at which psi0 and psi1 are evaluated for nuplot.
!
!      t_nudata : the time when the last data dump to nuplot was made.
!
!      psi0dat(k,n) : the time integrated psi0, to be time averaged on dumping.
!
!      psi1dat(k,n) : the time integrated psi1, to be time averaged on dumping.
!-----------------------------------------------------------------------

nudata            19         0                                          inudata
nudata                              1.00000000E+00                      dtnudata
nudata                              1.00000000E+08                      r_nudata

!-----------------------------------------------------------------------
!        Parameters for pinch parameter plot edits
!
!      r_pinch : the radius at which to evaluate the pinch parameters.
!
!      d_pinch : the density at which to evaluate the pinch parameters.
!-----------------------------------------------------------------------

pinch                               1.00000000E+08                      r_pinch
pinch                               1.00000000E+12                      d_pinch

cccccc   MD Edits
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Parameters for multidimensional edits
!
!      i_MDedit   : switch for dumping files for MD edits.
!
!      nd_MDedit  : the unit number for dumping files for MD edits.
!
!      n_MDedit   : MD edit counter.
!
!      dt_MDedit1 : dump files for MD edits every dt_MDedit1 ms before bounce.
!
!      dt_MDedit2 : dump files for MD edits every dt_MDedit1 ms after bounce.
!
!-----------------------------------------------------------------------

MDedit                       0                                          i_MDedit
MDedit                      31                                          nd_MDedit
MDedit                       0                                          n_MDedit
MDedit                              5.00000000E+00                      dt_MDedit1
MDedit                              5.00000000E+00                      dt_MDedit2
