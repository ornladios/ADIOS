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

intprt                       0                                          iprtbgn

!-----------------------------------------------------------------------
!        Obsolete parameter.
!-----------------------------------------------------------------------

noutfl                      19                                          noutfl

!-----------------------------------------------------------------------
!          0: print to a print file is bypassed.
!          1: print to a print file is implemented.
!-----------------------------------------------------------------------

output                       0                                          iprint

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

output                       6                                          nouttmp

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

intrst                 1000000                                          intrst

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

intprm                      50                                          intprm

intprm                       0                                          nprm

intrst                     600                                          ncychg

intrst                       0                                          irstbgn

ncyrst             1    900000                                          ncyrst
ncyrst             2    900010                                          ncyrst
ncyrst             3    900020                                          ncyrst
ncyrst             4    900030                                          ncyrst
ncyrst             5    905040                                          ncyrst
ncyrst             6    906050                                          ncyrst
ncyrst             7    907060                                          ncyrst
ncyrst             8    908070                                          ncyrst
ncyrst             9    909080                                          ncyrst
ncyrst            10    910000                                          ncyrst
ncyrst            11    912000                                          ncyrst
ncyrst            12    914000                                          ncyrst
ncyrst            13    916000                                          ncyrst
ncyrst            14    918000                                          ncyrst
ncyrst            15    920000                                          ncyrst
ncyrst            16    922000                                          ncyrst
ncyrst            17    924000                                          ncyrst
ncyrst            18    926000                                          ncyrst
ncyrst            19    928000                                          ncyrst
ncyrst            20    930000                                          ncyrst
ncyrst            21    935000                                          ncyrst
ncyrst            22    940000                                          ncyrst
ncyrst            23    945000                                          ncyrst
ncyrst            24    950000                                          ncyrst
ncyrst            25    960000                                          ncyrst
ncyrst            26    970000                                          ncyrst
ncyrst            27    980000                                          ncyrst
ncyrst            28    990000                                          ncyrst
ncyrst            29    100000                                          ncyrst
ncyrst            30    120000                                          ncyrst
ncyrst            31    140000                                          ncyrst
ncyrst            32    160000                                          ncyrst
ncyrst            33    180000                                          ncyrst
ncyrst            34    200000                                          ncyrst
ncyrst            35    220000                                          ncyrst
ncyrst            36    240000                                          ncyrst
ncyrst            37    260000                                          ncyrst
ncyrst            38    280000                                          ncyrst
ncyrst            39    300000                                          ncyrst
ncyrst            40    320000                                          ncyrst
ncyrst            41    340000                                          ncyrst
ncyrst            42    360000                                          ncyrst
ncyrst            43    380000                                          ncyrst
ncyrst            44    400000                                          ncyrst
ncyrst            45    420000                                          ncyrst
ncyrst            46    440000                                          ncyrst
ncyrst            47    460000                                          ncyrst
ncyrst            48    480000                                          ncyrst
ncyrst            49    500000                                          ncyrst
ncyrst            50    520000                                          ncyrst
ncyrst            51    540000                                          ncyrst
ncyrst            52    560000                                          ncyrst
ncyrst            53    580000                                          ncyrst
ncyrst            54    600000                                          ncyrst
ncyrst            55    620000                                          ncyrst
ncyrst            56    640000                                          ncyrst
ncyrst            57    660000                                          ncyrst
ncyrst            58    680000                                          ncyrst
ncyrst            59    700000                                          ncyrst
ncyrst            60    750000                                          ncyrst
ncyrst            61    800000                                          ncyrst
ncyrst            62    850000                                          ncyrst
ncyrst            63    900000                                          ncyrst
ncyrst            64    950000                                          ncyrst
ncyrst            65   1000000                                          ncyrst

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

lumplt                       1                                          ilumplt
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

enuplt                       1                                          ienuplt
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

rnuplt                       1                                          irnuplt
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

bnuplt            15         1                                          shkplt

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

rlagpt            31         1                                          irlagplt
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

nurad             18         1                                          inuplt
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

nudata            19         1                                          inudata
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
!        Parameters for M-dimensional u (radial) edits
!
!      iedMDu    : switch for performing MD u (radial) edits
!
!      nedMDu    : edit counter
!
!      intedMDu  : edit every intedMDu cycles
!
!      n_editMDu : model number counter
!
!-----------------------------------------------------------------------

MDedu                        0                                          iedMDu
MDedu                        0                                          nedMDu
MDedu                      200                                          intedMDu
MDedu                        0                                          n_editMDu

!-----------------------------------------------------------------------
!        Parameters for N-dimensional v (angular) edits
!
!      iedMDv    : switch for performing MD v (angular) edits
!
!      nedMDv    : edit counter
!
!      intedMDv  : edit every intedMDv cycles
!
!      n_editMDv : model number counter
!
!-----------------------------------------------------------------------

MDedv                        0                                          iedMDv
MDedv                        0                                          nedMDv
MDedv                      200                                          intedMDv
MDedv                        0                                          n_editMDv

!-----------------------------------------------------------------------
!        Parameters for N-dimensional w (azimuthal) edits
!
!      iedMDw    : switch for performing MD w (azimuthal) edits
!
!      nedMDw    : edit counter
!
!      intedMDw  : edit every intedMDw cycles
!
!      n_editMDw : model number counter
!
!-----------------------------------------------------------------------

MDedw                        0                                          iedMDw
MDedw                        0                                          nedMDw
MDedw                      200                                          intedMDw
MDedw                        0                                          n_editMDw

!-----------------------------------------------------------------------
!        Parameters for M-dimensional s (entropy) edits
!
!      iedMDs    : switch for performing MD s-edits
!
!      nedMDs    : edit counter
!
!      intedMDs  : edit every intedMDs cycles
!
!      n_editMDs : model number counter
!
!-----------------------------------------------------------------------

MDeds                        0                                          iedMDs
MDeds                        0                                          nedMDs
MDeds                      200                                          intedMDs
MDeds                        0                                          n_editMDs

!-----------------------------------------------------------------------
!        Parameters for M-dimensional d-edits
!
!      iedMDd    : switch for performing MD d-edits
!
!      nedMDd    : edit counter
!
!      intedMDd  : edit every intedMDd cycles
!
!      n_editMDd : model number counter
!
!-----------------------------------------------------------------------

MDedd                        0                                          iedMDd
MDedd                        0                                          nedMDd
MDedd                      200                                          intedMDd
MDedd                        0                                          n_editMDd

!-----------------------------------------------------------------------
!        Parameters for M-dimensional e-edits
!
!      iedMDe    : switch for performing MD e-edits
!
!      nedMDe    : edit counter
!
!      intedMDe  : edit every intedMDe cycles
!
!      n_editMDe : model number counter
!
!-----------------------------------------------------------------------

MDede                        0                                          iedMDe
MDede                        0                                          nedMDe
MDede                      200                                          intedMDe
MDede                        0                                          n_editMDe

!-----------------------------------------------------------------------
!        Parameters for M-dimensional p-edits
!
!      iedMDp    : switch for performing MD p-edits
!
!      nedMDp    : edit counter
!
!      intedMDp  : edit every intedMDp cycles
!
!      n_editMDp : model number counter
!
!-----------------------------------------------------------------------

MDedp                        0                                          iedMDp
MDedp                        0                                          nedMDp
MDedp                      200                                          intedMDp
MDedp                        0                                          n_editMDp

!-----------------------------------------------------------------------
!        Parameters for M-dimensional enu-edits
!
!      iedMDenu    : switch for performing MD enu-edits
!
!      nedMDenu    : edit counter
!
!      intedMDenu  : edit every intedMDenu cycles
!
!      n_editMDenu : model number counter
!
!-----------------------------------------------------------------------

MDedenu                      0                                          iedMDenu
MDedenu                      0                                          nedMDenu
MDedenu                    200                                          intedMDenu
MDedenu                      0                                          n_editMDenu

!-----------------------------------------------------------------------
!        Parameters for M-dimensional fnu-edits
!
!      iedMDfnu    : switch for performing MD fnu-edits
!
!      nedMDfnu    : edit counter
!
!      intedMDfnu  : edit every intedMDfnu cycles
!
!      n_editMDfnu : model number counter
!
!-----------------------------------------------------------------------

MDedfnu                      0                                          iedMDfnu
MDedfnu                      0                                          nedMDfnu
MDedfnu                    200                                          intedMDfnu
MDedfnu                      0                                          n_editMDfnu

!-----------------------------------------------------------------------
!        Parameters for M-dimensional a (mean mass number) edits
!
!      iedMDa    : switch for performing MD a (mean mass number) edits
!
!      nedMDa    : edit counter
!
!      intedMDa  : edit every intedMDa cycles
!
!      n_editMDa : model number counter
!
!-----------------------------------------------------------------------

MDeda                        0                                          iedMDa
MDeda                        0                                          nedMDa
MDeda                      200                                          intedMDa
MDeda                        0                                          n_editMDa

!-----------------------------------------------------------------------
!        Parameters for M-dimensional x (parameter) edits
!
!      iedMDx    : switch for performing MD x (parameter) edits
!
!      nedMDx    : edit counter
!
!      intedMDx  : edit every intedMDx cycles
!
!      n_editMDx : model number counter
!
!-----------------------------------------------------------------------

MDedx                        0                                          iedMDx
MDedx                        0                                          nedMDx
MDedx                      200                                          intedMDx
MDedx                        0                                          n_editMDx

!-----------------------------------------------------------------------
!        Parameters for M-dimensional ye (electron fraction) edits
!
!      iedMDye    : switch for performing MD ye (electron fraction) edits
!
!      nedMDye    : edit counter
!
!      intedMDye  : edit every intedMDye cycles
!
!      n_editMDye : model number counter
!
!-----------------------------------------------------------------------

MDedye                       0                                          iedMDye
MDedye                       0                                          nedMDye
MDedye                     200                                          intedMDye
MDedye                       0                                          n_editMDye

!-----------------------------------------------------------------------
!        Parameters for M-dimensional cm-edits
!
!      iedMDcm    : switch for performing MD cm-edits
!
!      nedMDcm    : edit counter
!
!      intedMDcm  : edit every intedMDcm cycles
!
!      n_editMDcm : model number counter
!
!-----------------------------------------------------------------------

MDedcm                       0                                          iedMDcm
MDedcm                       0                                          nedMDcm
MDedcm                     200                                          intedMDcm
MDedcm                       0                                          n_editMDcm

!-----------------------------------------------------------------------
!        Parameters for M-dimensional nu-edits
!
!      iedMDnu    : switch for performing MD nu-edits
!
!      nedMDnu    : edit counter
!
!      intedMDnu  : edit every intedMDnu cycles
!
!      n_editMDnu : model number counter
!
!-----------------------------------------------------------------------

MDednu                       0                                          iedMDnu
MDednu                       0                                          nedMDnu
MDednu                     200                                          intedMDnu
MDednu                       0                                          n_editMDnu

!-----------------------------------------------------------------------
!        Parameters for M-dimensional nc-edits
!
!      iedMDnc    : switch for performing MD nc-edits
!
!      nedMDnc    : edit counter
!
!      intedMDnc  : edit every intedMDnc cycles
!
!      n_editMDnc : model number counter
!
!-----------------------------------------------------------------------

MDednc                       0                                          iedMDnc
MDednc                       0                                          nedMDnc
MDednc                     200                                          intedMDnc
MDednc                       0                                          n_editMDnc

!-----------------------------------------------------------------------
!        Parameters for M-dimensional nl (neutrino luminosity) edits
!
!      iedMDnl    : switch for performing MD ye-edits
!
!      nedMDnl    : edit counter
!
!      intedMDnl  : edit every intedMDnl cycles
!
!      n_editMDnl : model number counter
!
!-----------------------------------------------------------------------

MDednl                       0                                          iedMDnl
MDednl                       0                                          nedMDnl
MDednl                     200                                          intedMDnl
MDednl                       0                                          n_editMDnl

!-----------------------------------------------------------------------
!        Parameters for M-dimensional ne (neutrino rms energy) edits
!
!      iedMDne    : switch for performing MD ne-edits
!
!      nedMDne    : edit counter
!
!      intedMDne  : edit every intedMDne cycles
!
!      n_editMDne : model number counter
!
!-----------------------------------------------------------------------

MDedne                       0                                          iedMDne
MDedne                       0                                          nedMDne
MDedne                     200                                          intedMDne
MDedne                       0                                          n_editMDne

!-----------------------------------------------------------------------
!        Parameters for M-dimensional gx (gravitational x-force) edits
!
!      iedMDgx    : switch for performing MD x-gravitational acceleration edits
!
!      nedMDgx    : edit counter
1
!      intedMDgx  : edit every intedMDgx cycles
!
!      n_editMDgx : model number counter
!
!-----------------------------------------------------------------------

MDedgx                       0                                          iedMDgx
MDedgx                       0                                          nedMDgx
MDedgx                     200                                          intedMDgx
MDedgx                       0                                          n_editMDgx

!-----------------------------------------------------------------------
!        Parameters for M-dimensional gy (gravitational y-force) edits
!
!      iedMDgy    : switch for performing MD y-gravitational acceleration edits
!
!      nedMDgy    : edit counter
1
!      intedMDgy  : edit every intedMDgy cycles
!
!      n_editMDgy : model number counter
!
!-----------------------------------------------------------------------

MDedgy                       0                                          iedMDgy
MDedgy                       0                                          nedMDgy
MDedgy                     200                                          intedMDgy
MDedgy                       0                                          n_editMDgy

!-----------------------------------------------------------------------
!        Parameters for M-dimensional gz (gravitational z-force) edits
!
!      iedMDgz    : switch for performing MD y-gravitational acceleration edits
!
!      nedMDgz    : edit counter
1
!      intedMDgz  : edit every intedMDgz cycles
!
!      n_editMDgz : model number counter
!
!-----------------------------------------------------------------------

MDedgz                       0                                          iedMDgz
MDedgz                       0                                          nedMDgz
MDedgz                     200                                          intedMDgz
MDedgz                       0                                          n_editMDgz

!-----------------------------------------------------------------------
!        Parameters for M-dimensional BVw (Brunt-Vaisala) edits
!
!      iedMDBVw    : switch for performing MD Brunt-Vaisala edits
!
!      nedMDBVw    : edit counter
1
!      intedMDBVw  : edit every intedMDBVw cycles
!
!      n_editMDBVw : model number counter
!
!-----------------------------------------------------------------------

MDedBV                       0                                          iedMDBVw
MDedBV                       0                                          nedMDBVw
MDedBV                     200                                          intedMDBVw
MDedBV                       0                                          n_editMDBVw

!-----------------------------------------------------------------------
!        Parameters for M-dimensional yl (lepton fraction) edits
!
!      iedMDyl    : switch for performing MD yl (lepton fraction) edits
!
!      nedMDyl    : edit counter
!
!      intedMDyl  : edit every intedMDyl cycles
!
!      n_editMDyl : model number counter
!
!-----------------------------------------------------------------------

MDedyl                       0                                          iedMDyl
MDedyl                       0                                          nedMDyl
MDedyl                     200                                          intedMDyl
MDedyl                       0                                          n_editMDyl

!-----------------------------------------------------------------------
!        Time parameters for M-D edits
!
!      i_editMD   : 2-D edit time switch
!      n_editMD   : model number counter
!      dt_MDedit1 : dump files for MD edits every dt_MDedit1 ms before bounce.
!      dt_MDedit2 : dump files for MD edits every dt_MDedit2 ms after bounce.
!-----------------------------------------------------------------------

MDedit                       1                                          i_editMD
MDedit                       0                                          n_editMD
MDedit                              5.00000000E+00                      dt_MDedit1
MDedit                              2.00000000E-01                      dt_MDedit2

cccccc   Global Edits
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Global edit cycle parameters
!
!      ied_global_n  : switch for performing global using cycle criterion
!      nned_global   : edit counter.
!      inted_global  : edit every inted_global cycles.
!      n_edit_global : model number counter.
!-----------------------------------------------------------------------

ed_gln                       0                                          ied_global_n
ed_gln                       0                                          nned_global
ed_gln                     200                                          inted_global
ed_gln                       0                                          n_edit_global

!-----------------------------------------------------------------------
!        Global time cycle parameters
!
!      ied_global_t  : switch for performing 2D edits at selected times
!      nted_global   : model number counter
!      dt_global_ed1 : dump files for global edits every dt_global_ed1 ms
!                       before bounce.
!      dt_global_ed2 : dump files for global edits every dt_global_ed2 ms
!                       after bounce.
!-----------------------------------------------------------------------

ed_glt                       1                                          ied_global_t
ed_glt                       0                                          nted_global
ed_glt                              5.00000000E+00                      dt_global_ed1
ed_glt                              2.00000000E+00                      dt_global_ed2

cccccc   HDF Edits
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Parameters for HDF edits
!
!      i_HDFedit   : switch for dumping files for HDF edits.
!
!      nd_HDFedit  : the unit number for dumping files for HDF edits.
!
!      n_HDFedit   : HDF edit counter.
!
!      dt_HDFedit1 : dump files for HDF edits every dt_HDFedit1 ms before bounce.
!
!      dt_HDFedit2 : dump files for HDF edits every dt_HDFedit1 ms after bounce.
!
!-----------------------------------------------------------------------

HDFedt                       0                                          i_HDFedit
HDFedt                      31                                          nd_HDFedit
HDFedt                       0                                          n_HDFedit
HDFedt                              5.00000000E+00                      dt_HDFedit1
HDFedt                              5.00000000E+00                      dt_HDFedit2
