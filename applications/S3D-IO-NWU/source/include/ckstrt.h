C       CVS Revision:$Revision: 1.1.1.1 $  created $Date: 2005/01/11 00:27:29 $
C	this is CHEMKIN-III file ckstrt.h V.3.1 December 2000;
C	it contains pointers for the gas-phase kinetics
C       subroutines' data storage arrays
C
C
C     include file for CHEMKIN-III cklib.f, dated: December 8, 2000
C
      INTEGER
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD, 
     4   MXORD, KEL,  NKKI,
     5   IcMM, IcKK,
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNSU,IcNS, IcNR, IcLT, 
     7   IcRL, IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, 
     8   IcEI, IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, 
     9   IcTD, IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcAVO,NcPER,
     2   NcECH,NcBOL,NcPI, NcKF, NcKR, NcRNU,
     3   NcRSU,NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4


      COMMON /CKSTRT/ 
C
C     Integer constants
C
C        0     1     2     3     4     5     6     7     8     9
     1   NMM,  NKK,  NII,  MXSP, MXTB, MXTP, NCP,  NCP1, NCP2, NCP2T,
C        10    11    12    13    14    15    16    17    18    19   
     2   NPAR, NLAR, NFAR, NLAN, NFAL, NREV, NTHB, NRLT, NWL,  NEIM,
C        20    21    22    23    24    25    26    27    28    29
     3   NJAN, NJAR, NFT1, NF1R, NEXC, NMOM, NXSM, NTDE, NRNU, NORD, 
C        30     31    32
     4   MXORD, KEL,  NKKI,
C
C     Integer pointers to character arrays in CCKWRK
C
C        33    34
     5   IcMM, IcKK,
C
C     Integer pointers to integer arrays in ICKWRK
C
C        35    36    37    38    39    40    41    42    43    44     
     6   IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNSU,IcNS, IcNR, IcLT, 
C        45    46    47    48    49    50    51    52    53    54
     7   IcRL, IcRV, IcWL, IcFL, IcFO, IcFT, IcKF, IcTB, IcKN, IcKT, 
C        55    56    57    58    59    60    61    62    63    64    
     8   IcEI, IcET, IcJN, IcF1, IcEX, IcMO, IcMK, IcXS, IcXI, IcXK, 
C        65    66    67    68    69    70    71    72    73   
     9   IcTD, IcTK, IcRNU,IcORD,IcKOR,IcKI, IcKTF,IcK1, IcK2,
C
C     Integer pointers to real variables and arrays in RCKWRK
C
C        74    75    76    77    78    79    80    81    82    83    
     *   NcAW, NcWT, NcTT, NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT,
C        84    85    86    87    88    89    90    91     92
     1   NcWL, NcJN, NcF1, NcEX, NcRU, NcRC, NcPA, NcAVO, NcPER,
C        93     94     95    96    97    98 
     2   NcECH, NcBOL, NcPI, NcKF, NcKR, NcRNU,
C        99    100   101   102   103   104   105   106   107   108
     3   NcRSU,NcKOR,NcK1, NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4


C Logical global variables
C
C      NPERT = TRUE if any of the perturbation factors for the rate
C              constants have been changed from the value of one.
C
       LOGICAL         LPERT
       COMMON /CKGLBL/ LPERT
C
C     END include file for cklib.f
C
