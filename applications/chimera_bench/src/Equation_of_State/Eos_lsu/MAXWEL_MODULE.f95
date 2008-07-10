!-----------------------------------------------------------------------
!    MODULE:       maxwel_module
!    Author:       S. W. Bruenn
!    Date:         11/20/02
!-----------------------------------------------------------------------
MODULE maxwel_module

USE kind_module
SAVE

!.........Number of T points in each boundary...........................

INTEGER, PARAMETER                              :: NUMTMP = 201

!.........Number of Ye points in each boundary..........................

INTEGER, PARAMETER                              :: NUMYE = 49

!.........Number of points in alpha-nuclei boundary.....................

INTEGER, PARAMETER                              :: NBPNTS = 101
INTEGER, PARAMETER                              :: NUMLOW = 51
INTEGER, PARAMETER                              :: NUMHI = 51

!.........Baryon densities..............................................

REAL(KIND=double), DIMENSION(NUMTMP,NUMYE)      :: BRYLOW
REAL(KIND=double), DIMENSION(NUMTMP,NUMYE)      :: BRYHI
REAL(KIND=double)                               :: LOWDNS, HIDNS, DNS_1, DNS_2
REAL(KIND=double)                               :: DNL_DT, DNH_DT, DNL_DY, DNH_DY
REAL(KIND=double), PARAMETER                    :: LNLOW = -8.92D+00
REAL(KIND=double), PARAMETER                    :: LNHI = -0.92D+00
REAL(KIND=double), PARAMETER                    :: LNCUT = -2.92D+00
REAL(KIND=double)                               :: LNMINS, LNPLUS
REAL(KIND=double)                               :: LOGBRY, LOGBCH, DLTLN1, DLTLN2, LNFRAC
REAL(KIND=double), PARAMETER                    :: YLOW = 0.16D+00
REAL(KIND=double), PARAMETER                    :: YHI = 0.51D+00
REAL(KIND=double), PARAMETER                    :: Y_CUT = 0.155D+00
REAL(KIND=double)                               :: TCHK_B, TCHK_N, T_MXWL, D_MXWL

!.........Highest temperatures where Coexistence........................
!.........of bulk and nuclear phases occurs.............................

REAL(KIND=double), DIMENSION(NUMYE)             :: T_H
REAL(KIND=double), DIMENSION(NUMYE)             :: D_H

!.........Arrays containing number of boundary points...................

REAL(KIND=double), DIMENSION(NBPNTS,NUMYE)      :: LBOUND
REAL(KIND=double), DIMENSION(NBPNTS,NUMYE)      :: UBOUND

!.........Minimum density code will work at.............................

REAL(KIND=double), PARAMETER                    :: MINDNS = 1.0D-10

REAL(KIND=double)                               :: NEWDNS
REAL(KIND=double)                               :: T_LOW, T_HI, TFRAC
REAL(KIND=double)                               :: Y_LOW, Y_LOW2, Y_HI, Y_HI2, YFRAC, YCH, YCUT
REAL(KIND=double)                               :: LNL, LNH, LNC

REAL(KIND=double)                               :: YMINUS, YPLUS, YINTRP, TMINUS, TPLUS, TINTRP
REAL(KIND=double)                               :: DELT_Y, DELT_T
INTEGER                                         :: I_MXWL, J_MXWL
INTEGER                                         :: I_BD, I_BNDY, J_BD, J_BNDY

!.........Total chemical potential......................................

REAL(KIND=double)                               :: MUTILD, MUTLOW, MUTHI

!.........Total pressure................................................

REAL(KIND=double)                               :: PRTILD, PRLOW, PRHI

!.........Muhat.........................................................

REAL(KIND=double)                               :: MUHLOW, MUHHI

!.........Electron chemical potential...................................

REAL(KIND=double)                               :: MUELOW, MUEHI

!.........Total entropy per baryon......................................

REAL(KIND=double)                               :: S_LOW, S_HI

!.........Total free energy density.....................................

REAL(KIND=double)                               :: F_LOW, F_HI

!.........Phase fraction................................................

REAL(KIND=double)                               :: PHASEF

END MODULE maxwel_module

