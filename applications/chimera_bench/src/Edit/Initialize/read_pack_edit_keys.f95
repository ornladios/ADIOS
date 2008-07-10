SUBROUTINE read_pack_edit_keys( nread, nprint, iskip, nez, nnu, i_edit_data, &
& d_edit_data, nrst )
!-----------------------------------------------------------------------
!
!    File:         read_pack_edit_keys
!    Module:       read_pack_edit_keys
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         10/26/04
!
!    Purpose:
!      To read in the edit keys determining what and when edits are to be done,
!       and to pack them into an integer and a real*8 array.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  nread       : unit number from which to read
!  nprint      : unit number from which to print
!  iskip       : echo data read flag
!  nez         : energy array dimension
!  nnu         : neutrino flavor dimension
!  nrst        : cycle number at start or restart
!
!    Output arguments:
!  i_edit_data : integer array of edit keys
!  d_edit_data : real*8 array of edit keys
!
!    Include files:
!  kind_module, numerical_module
!  edit_module
!
!-----------------------------------------------------------------------

USE kind_module, ONLY: double
USE numerical_module, ONLY : zero

USE edit_module, ONLY : nlog

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(in)              :: nread         ! unit number to read from
INTEGER, INTENT(in)              :: nprint        ! unit number to print to
INTEGER, INTENT(in)              :: iskip         ! echo data read flag
INTEGER, INTENT(in)              :: nez           ! energy array dimension
INTEGER, INTENT(in)              :: nnu           ! neutrino flavor dimension
INTEGER, INTENT(in)              :: nrst          ! cycle number at start or restart

!-----------------------------------------------------------------------
!        Output variables.
!-----------------------------------------------------------------------

INTEGER, INTENT(out), DIMENSION(1200+3*40*nnu) :: i_edit_data  ! integer array of edit keys

REAL(KIND=double), INTENT(out), DIMENSION(50)  :: d_edit_data  ! 64 bit real array of edit keys

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

CHARACTER (len=6)                :: type          ! 5 character variable to read in data chatagory type
CHARACTER (len=16)               :: name          ! 16 character variable to read in data name
CHARACTER (len=128)              :: line          ! 120 character variable to read in data line
CHARACTER (len=10)               :: var_name

INTEGER                          :: i             ! do index
INTEGER                          :: k             ! neutrino eneregy index
INTEGER                          :: n             ! neutrino flavor index
INTEGER                          :: istat         ! allocation status

INTEGER                          :: int           ! integer data variable to read in an interger datum

REAL(KIND=double)                :: rl            ! 64 bit real variable read in a real datum

!-----------------------------------------------------------------------
!
!               \\\\\ GENERAL EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  iflprt : a print parameter regulating the print file created when the
!   calculation is terminated (ncycle = ncymax)
!
!          iflprt = 0:  no printout at termination
!          iflprt = 1:  configuration printout at termination
!          iflprt = 2:  full printout except editn at termination
!          iflprt = 3:  full printout except editng at termination
!          iflprt = 4:  full printout at termination
!
!  iprint : the parameter used in the call to an edit subroutine;
!
!          0: print to a print file is bypassed.
!          1: print to a print file is implemented.
!-----------------------------------------------------------------------

INTEGER                                             :: iprint

INTEGER                                             :: iflprt
INTEGER                                             :: noutfl

!-----------------------------------------------------------------------
!  nrstd1 : the unit number for file 'rstdmp1.d' to which temporary restart dumps are written
!   every intrst cycles. These restart dumps are temporary in the sense that the next restart
!   dump written to file 'rstdmp1.d' writes over the preceding dump.
!
!  nrstd2 : the unit number for file 'rstdmp2.d' to which temporary restart dumps can be written
!   every other intrst cycles alternating with file rstdmp1.d.
!-----------------------------------------------------------------------

INTEGER                                             :: nrstd1
INTEGER                                             :: nrstd2

!-----------------------------------------------------------------------
!  noutpmt : the unit number of the permanent restart dump file (typicallyfile
!  'restartxxxxxxx.d',  where xxxxxxx is the cycle number right-justified).
!-----------------------------------------------------------------------

INTEGER                                             :: noutpmt

!-----------------------------------------------------------------------
!  nouttmp : the unit number of the temporary restart dump file (typically
!   file 'rstdmp1.d' or 'rstdmp2.d').
!-----------------------------------------------------------------------

INTEGER                                             :: nouttmp

!-----------------------------------------------------------------------
!  iplot  : a parameter used in the call to an edit subroutine -
!
!          0: print to a plot file is bypassed.
!          1: print to a plot file is implemented.
!
!  nplotc : the unit number for downloading model configuration plot data.
!
!  nplote : the unit number for downloading electron neutrino plot data.
!
!  nplota : the unit number for downloading electron antineutrino plot data.
!
!  nplott : the unit number for downloading muon and tau neutrino plot data.
!
!  intplf : the number of cycles between the printing to a plot file. 
!
!  npltf  : the number of cycles since the last printing to a plot file.
!-----------------------------------------------------------------------

INTEGER                                             :: iplot
INTEGER                                             :: nplotc
INTEGER                                             :: nplote
INTEGER                                             :: nplota
INTEGER                                             :: nplott
INTEGER                                             :: intplf
INTEGER                                             :: npltf

!-----------------------------------------------------------------------
!  nmodel : the number of the current print file.
!-----------------------------------------------------------------------

INTEGER                                             :: nmodel

!-----------------------------------------------------------------------
!  rhoprint(i): the ith central density at which a print file is created.
!   The print file is named modeldxx.d, where xx is the value of i right
!   justified.
!-----------------------------------------------------------------------

REAL(KIND=double), DIMENSION(10)                    :: rhoprint

!-----------------------------------------------------------------------
!  Parameters for configuration edits
!-----------------------------------------------------------------------
!  intedc(i) : the number of cycles between the implementation of
!   subsection i of subroutine editc. 
!
!  nedc(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine editc.
!
!  idxedc(i) : a parameter used in subroutine editc. Subsection i of
!   subroutine editc will print data for every idxedc(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedc(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedc
INTEGER, DIMENSION(20)                              :: nedc
INTEGER, DIMENSION(20)                              :: idxedc

!-----------------------------------------------------------------------
!  Parameters for kinetic, internal, and gravitational energy edits
!-----------------------------------------------------------------------
!  intede(i) : the number of cycles between the implementation of
!   subsection i of subroutine edit_e. 
!
!  nede(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine edit_e.
!
!  idxede(i) : a parameter used in subroutine edit_e. Subsection i of
!   subroutine editc will print data for every idxede(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxede(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intede
INTEGER, DIMENSION(20)                              :: nede
INTEGER, DIMENSION(20)                              :: idxede

!-----------------------------------------------------------------------
!  Parameters for editmi edits
!-----------------------------------------------------------------------
!  intdmi(i) : the number of cycles between the implementation of
!   subsection i of subroutine editmi.
!
!  nedmi(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editmi.
!
!  idxemi(i) : a parameter used in subroutine editmi. Subsection i of
!   subroutine editmi will print
!   data for every idxemi(i)'th radial zone, starting with the outermost
!   zone. The data corresponding to the innermost radial zone will also
!   be printed. Setting idxemi(i) = 1 will cause data for all radial
!   zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdmi
INTEGER, DIMENSION(20)                              :: nedmi
INTEGER, DIMENSION(20)                              :: idxemi

!-----------------------------------------------------------------------
!  Parameters for mass average edits
!-----------------------------------------------------------------------
!  intdma(i) : the number of cycles between the implementation of
!   subsection i of subroutine editma.
!
!  nedma(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editma.
!
!  idxema(i) : a parameter used in subroutine editma. Subsection i of
!   subroutine editma will print data for every idxema(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxema(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdma
INTEGER, DIMENSION(20)                              :: nedma
INTEGER, DIMENSION(20)                              :: idxema

!-----------------------------------------------------------------------
!  Parameters for hydrodynamic edits
!-----------------------------------------------------------------------
!  intedh(i) : the number of cycles between the implementation of
!   subsection i of subroutine edith.
!
!  nedh(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine edith.
!
!  idxedh(i) : a parameter used in subroutine edith. Subsection i of
!   subroutine editma will print data for every idxedh(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedh(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedh
INTEGER, DIMENSION(20)                              :: nedh
INTEGER, DIMENSION(20)                              :: idxedh

!-----------------------------------------------------------------------
!  Parameters for pressure - stress edits
!-----------------------------------------------------------------------
!  intdps(i) : the number of cycles between the implementation of
!   subsection i of subroutine editps.
!
!  nedps(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editps.
!
!  idxeps(i) : a parameter used in subroutine editps. Subsection i of
!   subroutine editps will print data for every idxeps(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxeps(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdps
INTEGER, DIMENSION(20)                              :: nedps
INTEGER, DIMENSION(20)                              :: idxeps

!-----------------------------------------------------------------------
!  Parameters for energy edits
!-----------------------------------------------------------------------
!  intedu(i) : the number of cycles between the implementation of
!   subsection i of subroutine editu.
!
!  nedu(i)   : the number of cycles since the last implementation of
!   subsection i of subroutine editu.
!
!  idxedu(i) : a parameter used in subroutine editu. Subsection i of
!   subroutine editu will print data for every idxedu(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedu(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedu
INTEGER, DIMENSION(20)                              :: nedu
INTEGER, DIMENSION(20)                              :: idxedu

!-----------------------------------------------------------------------
!  Parameters for composition edits
!-----------------------------------------------------------------------
!  intedy(i): the number of cycles between the implementation of
!   subsection i of subroutine edity.
!
!  nedy(i): the number of cycles since the last implementation of
!   subsection i of subroutine edity.
!
!  idxedy(i): a parameter used in subroutine edity. Subsection i of
!   subroutine edity will print data for every idxedy(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxedy(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intedy
INTEGER, DIMENSION(20)                              :: nedy
INTEGER, DIMENSION(20)                              :: idxedy

!-----------------------------------------------------------------------
!  Parameters for entropy and chemical potential edits
!-----------------------------------------------------------------------
!  intdsc(i) : the number of cycles between the implementation of
!   subsection i of subroutine editsc.
!
!  nedsc(i)  : the number of cycles since the last implementation of
!   subsection i of subroutine editsc.
!
!  idxesc(i) : a parameter used in subroutine editsc. Subsection i of
!   subroutine editsc will print data for every idxesc(i)'th radial zone,
!   starting with the outermost zone. The data corresponding to the
!   innermost radial zone will also be printed. Setting idxesc(i) = 1 will
!   cause data for all radial zones to be printed.
!-----------------------------------------------------------------------

INTEGER, DIMENSION(20)                              :: intdsc
INTEGER, DIMENSION(20)                              :: nedsc
INTEGER, DIMENSION(20)                              :: idxesc

!-----------------------------------------------------------------------
!  Parameters for diferential neutrino edits
!-----------------------------------------------------------------------
!  intedn(n) : the number of cycles between the implementation of
!   subroutine editn for neutrinos of type n.
!
!  nedn(n)   : the number of cycles since the last implementation of
!   subroutine editn for neutrinos of type n.
!
!  idxedn(n) : a parameter used in subroutine editn. Subroutine editn
!   will print the data corresponding neutrinos of type n for every
!   idxedn(n)'th radial zone, starting with the outermost zone. The data
!   corresponding to the innermost radial zone will also be printed.
!   Setting idxedn(i) = 1 will cause data for all radial zones to be
!   printed.
!
!  niedn(i)  : a parameter used in subroutine editn. Subroutine editn
!   will print the i'th datum corresponding to neutrinos of type n every
!   niedn(i) implementation of editn.
!
!  neden(i)  : the number of implementations of editn since the last
!   printing of the i'th datum corresponding to neutrinos of type n.
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:)                  :: intedn
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: nedn
INTEGER, ALLOCATABLE, DIMENSION(:)                  :: idxedn
INTEGER, DIMENSION(60)                              :: niedn
INTEGER, DIMENSION(60)                              :: neden

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------
!  intdng(i,n) : the number of cycles between the implementation of
!   subsection i of subroutine editng for neutrinos of type n.
!
!  nedng(i,n)  : the number of cycles since the last implementation
!   of subsection i of subroutine editng for neutrinos of type n.
!
!  idxeng(i,n) : a parameter used in subroutine editng. Subsection i of
!   subroutine editng will print the data corresponding to every
!   idxesc(i)'th radial zone, starting with the outermost zone. The
!   data corresponding to the innermost radial zone will also be printed.
!   Setting idxeng(i,n) = 1 will cause data for all radial zones to be
!   printed. 
!-----------------------------------------------------------------------

INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: intdng
INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: nedng
INTEGER, ALLOCATABLE, DIMENSION(:,:)                :: idxeng

!-----------------------------------------------------------------------
!  Parameters for edits at selected postbounce time intervals
!-----------------------------------------------------------------------
!  it_edit : key for editing at selected postbounce time intervals.
!
!  dt_edit : postbounce time intervals for editing (ms).
!-----------------------------------------------------------------------

INTEGER                                             :: it_edit

REAL(KIND=double)                                   :: dt_edit

!-----------------------------------------------------------------------
!  Parameters for closing and opening edit files
!-----------------------------------------------------------------------
!  intprt  : the number of cycles between the closing of an old print
!   file and the opening of a new print file. The print files are named
!   modelxxx.d, where xxx is the right justified consecutive number of
!   the print file, e.g., the fifth print file is named 'model005.d'.
!
!  nprt    : a print counter giving the number of cycles since the last
!   closing of a print file.
!
!  iprtbgn:  a print parameter regulating the print file created when the
!   calculation is initiated (nrst = 0) or restarted (nrst ne 0);
!
!     iprtbgn = 0, nrst = 0 : abbreviated edit
!     iprtbgn = 0, nrst > 0 : configuration edit
!     iprtbgn = 1, nrst = 0 : abbreviated edit
!     iprtbgn = 1, nrst > 0 : abbreviated edit
!     iprtbgn = 2, nrst = 0 : full edit
!     iprtbgn = 2, nrst > 0 : configuration edit
!     iprtbgn = 3, nrst = 0 : full edit
!     iprtbgn = 3, nrst > 0 : abreviated edit
!     iprtbgn = 4, nrst = 0 : full edit
!     iprtbgn = 4, nrst > 0 : full edit
!-----------------------------------------------------------------------

INTEGER                                             :: intprt
INTEGER                                             :: nprt
INTEGER                                             :: iprtbgn

!-----------------------------------------------------------------------
!  Parameters for temperary restart dumps
!-----------------------------------------------------------------------
!  intrst : the number of cycles between temporary restart dumps.
!
!  nnrst  : the number of cycles since the last temporary restart dump.
!
!  nrstfl : a temporary restart dump parameter -
!-----------------------------------------------------------------------

INTEGER                                             :: intrst
INTEGER                                             :: nnrst
INTEGER                                             :: nrstfl

!-----------------------------------------------------------------------
!  irstbgn : restart dump parameter -
!
!     1   : temporary and permanent restart dumps are implemented at the
!            beginning of the run (ncycle = 0)
!     ne 1: temporary and permanent restart dumps are bypassed at the
!            beginning of the run (ncycle = 0)
!-----------------------------------------------------------------------

INTEGER                                             :: irstbgn

!-----------------------------------------------------------------------
!  ncychg : at termination, write to restart file with
!          ncymax = ncymax + ncychg.
!-----------------------------------------------------------------------

INTEGER                                             :: ncychg

!-----------------------------------------------------------------------
!  Parameters for permanent restart dumps
!-----------------------------------------------------------------------
!  intprm    : the number of cycles between permanent restart dumps.
!
!  nprm      : the number of cycles since the last permanent restart dump.
!
!  ncyrst(i) : the cycle number at which the ith prescribed permanent
!   restart dump is implemented.
!-----------------------------------------------------------------------

INTEGER                                             :: intprm
INTEGER                                             :: nprm
INTEGER, DIMENSION(100)                             :: ncyrst

!-----------------------------------------------------------------------
!  Parameters for neutrino luminosity plot edits
!-----------------------------------------------------------------------
!  ilumplt  : luminosity plot edit switch
!
!          0: bypass subroutine lumplot.
!          1: enter subroutine lumplot.
!
!  nlumplt1 : the unit number for editing lumplot data at r_lum.
!
!  nlumplt2 : the unit number for editing lumplot data at d_lum.
!
!  nlum     : number of cycles since the last implemtation of subroutine lumplot.
!
!  r_lum    : the radius at which to evaluate luminosity data.
!
!  d_lum    : the density at which to evaluate luminosity data.
!
!  intlum, ncylum:
!
!     ncycle < ncylum(1)             : subroutine lumplot is called every
!                                       intlum(1) cycles.
!     ncylum(1) < ncycle < ncylum(2) : subroutine lumplot is called every
!                                       intlum(2) cycles.
!     ncylum(2) < ncycle             : subroutine lumplot is called every
!                                       intlum(3) cycles.                                                     !-----------------------------------------------------------------------

INTEGER                                             :: ilumplt
INTEGER                                             :: nlumplt1
INTEGER                                             :: nlumplt2
INTEGER                                             :: nlum
INTEGER, DIMENSION(3)                               :: intlum
INTEGER, DIMENSION(2)                               :: ncylum

REAL(KIND=double)                                   :: r_lum
REAL(KIND=double)                                   :: d_lum

!-----------------------------------------------------------------------
!  Parameters for rms neutrino energy plot edits
!-----------------------------------------------------------------------
!  ienuplt : rms neutrino energy plot edit switch
!
!          0: bypass subroutine enuvplot.
!          1: enter subroutine enuvplot.
!
!  nenuplt : the unit number for editing enuvplot data.
!
!  nenu    : the number of cycles since the last implemtation of subroutine
!   enuvplot.
!
!  r_e_rms :  the radius at which to store neutrino energy data.
!
!  d_e_rms : the density at which to store neutrino energy data.
!
!  intenu, ncyenu:
!
!     ncycle < ncyenu(1)             : subroutine enuvplot is called every                                          !                                       intenu(1) cycles.
!     ncyenu(1) < ncycle < ncyenu(2) : subroutine enuvplot is called every
!                                       intenu(2) cycles.
!     ncyenu(2) < ncycle             : subroutine enuvplot is called every
!                                       intenu(3) cycles.
!-----------------------------------------------------------------------

INTEGER                                             :: ienuplt
INTEGER                                             :: nenuplt1
INTEGER                                             :: nenuplt2
INTEGER                                             :: nenu
INTEGER, DIMENSION(3)                               :: intenu
INTEGER, DIMENSION(2)                               :: ncyenu

REAL(KIND=double)                                   :: r_e_rms
REAL(KIND=double)                                   :: d_e_rms

!-----------------------------------------------------------------------
!  Parameters for selected radii plot edits
!-----------------------------------------------------------------------
!  irnuplt : selected radii plot edit switch
!
!          0: bypass subroutine rnuplot.
!          1: enter subroutine rnuplot.
!
!  nrnuplt : the unit number for editing rnuplot data.
!
!  nrum    : the number of cycles since the last implemtation of subroutine
!   rnuplot.
!
!  intrnu, ncyrnu:
!
!     ncycle < ncyrnu(1)             : subroutine rnuplot is called every
!                                       intrnu(1) cycles.
!     ncyrnu(1) < ncycle < ncyrnu(2) : subroutine rnuplot is called every
!                                       intrnu(2) cycles.      
!     ncyrnu(2) < ncycle             : subroutine rnuplot is called every
!                                       intrnu(3) cycles.
!-----------------------------------------------------------------------

INTEGER                                             :: irnuplt
INTEGER                                             :: nrnuplt
INTEGER                                             :: nrnu
INTEGER, DIMENSION(3)                               :: intrnu
INTEGER, DIMENSION(2)                               :: ncyrnu

!-----------------------------------------------------------------------
!  Parameters for configuration plot edits
!-----------------------------------------------------------------------
!  ivarplt   : configuration plot edit switch
!
!     0 : bypass subroutine varplot.
!     1 : enter subroutine varplot.
!
!  nvar      : the number of cycles since the last implemtation of subroutine varplot.                                         !
!
!  nvarint   : subroutine varplot is called every nvarint cycles.
!
!  nvarplt   : the unit number for editing varplot data.
!
!  nvardump  : varplot file number (varplot files are numbered sequentially.                                               !
!  dtvarplot : write to subroutine varplot every dtvarplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: ivarplt
INTEGER                                             :: nvar
INTEGER                                             :: nvarint
INTEGER                                             :: nvarplt
INTEGER                                             :: nvardump

REAL(KIND=double)                                   :: dtvarplot

!-----------------------------------------------------------------------
!  Parameters for comparison plot edits
!-----------------------------------------------------------------------
!  icomplt   : comparison plot edit switch
!
!     0 : bypass subroutine complot.
!     1 : enter subroutine complot.
!
!  ncomplt   : the unit number for editing complot data.
!
!  ncomdump  : the complot file number (complot files are numbered
!   sequentially.                                               !
!
!  dtcomplot : write to subroutine complot every dtcomplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: icomplt
INTEGER                                             :: ncomplt
INTEGER                                             :: ncomdump

REAL(KIND=double)                                   :: dtcomplot

!-----------------------------------------------------------------------
!  Parameters for bounday plot edits
!-----------------------------------------------------------------------
!  dtimeplot   : write to subroutine bnuplot every dtimeplot ms.
!
!  rinnerb     : quantities are interpolated to rinnerb for ibound data.
!
!  nplotinnerb : the unit number for ibound data.
!
!  iplotinnerb : inner boundary plot switch
!
!     0    : bypass writing ibound data.
!     ne 0 : write ibound data.
!
!  routerb     : quantities are interpolated to routerb for obound data.
!
!  nplotouterb : the unit number for obound data.
!
!  iplotouterb : outer boundary plot switch
!
!     0   : bypass writing obound data.
!     ne 0: write obound data.
!
!  r_lumerms   : luminosities and rms energies are interpolated to
!   r_lumerms for lbound data.                                            
!
!  nplotlum    : the unit number for lbound data.
!
!  iplotlum    : luminosities and rms energies plot switch
!
!     0   : bypass writing lbound data.
!     ne 0: write lbound data.
!
!  nplotshk    : the unit number for shock data.
!
!  iplotshk    : shock data plot switch
!
!     0    : bypass writing shock data.
!     ne 0 : write shock data.
!
!  nplotcnv    : the unit number for convection data.
!
!  iplotcnv    : convection plot switch
!
!     0    : bypass writing convection data.
!     ne 0 : write convection data.
!
!  nplotmss    : the unit number for mass data.
!
!  iplotmss    : mass data plot switch
!
!     0    : bypass writing mass data.
!     ne 0 : write mass data.
!-----------------------------------------------------------------------

INTEGER                                             :: nplotinnerb
INTEGER                                             :: iplotinnerb
INTEGER                                             :: nplotouterb
INTEGER                                             :: iplotouterb
INTEGER                                             :: nplotlum
INTEGER                                             :: iplotlum
INTEGER                                             :: nplotshk
INTEGER                                             :: iplotshk
INTEGER                                             :: nplotcnv
INTEGER                                             :: iplotcnv

REAL(KIND=double)                                   :: dtimeplot
REAL(KIND=double)                                   :: rinnerb
REAL(KIND=double)                                   :: routerb
REAL(KIND=double)                                   :: r_lumerms

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------
!  nplotnurad   : the unit number for nuradplot data.
!
!  iplotnurad   : nurad plot switch
!
!     0    : bypass writing nuradplot data.
!     ne 0 : write nuradplot data.
!
!  dtnuradplot: write to subroutine nuradplot every dtnuradplot ms.
!
!  r_nurad      : the radius at which to evaluate neutrino radiation.
!
!  rho_nurad    : the density at which to evaluate neutrino radiation.
!-----------------------------------------------------------------------

INTEGER                                             :: nplotnurad
INTEGER                                             :: iplotnurad

REAL(KIND=double)                                   :: dtnuradplot
REAL(KIND=double)                                   :: r_nurad
REAL(KIND=double)                                   :: rho_nurad

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------
!  nnudata      : the unit number for nuplot data.
!
!  inudata      : neutrino distribution plot switch
!
!     0    : bypass writing nuplot data.
!     ne 0 : write nuplot data.
!
!  dtnudata     : write to subroutine nuplot every dtnudata ms.
!
!  r_nudata     : the radius at which psi0 and psi1 are evaluated for nuplot.
!
!  t_nudata     : the time when the last data dump to nuplot was made.
!-----------------------------------------------------------------------

INTEGER                                             :: nnudata
INTEGER                                             :: inudata

REAL(KIND=double)                                   :: dtnudata
REAL(KIND=double)                                   :: r_nudata
REAL(KIND=double)                                   :: t_nudata

!-----------------------------------------------------------------------
!  Parameters for lagrangian plot edits
!-----------------------------------------------------------------------
!  nlagplt  : the unit number for lagplot data.
!
!  ilagplt  : lagrangian plot switch
!
!     0    : bypass writing lagplot data.
!     ne 0 : write nuplot data.
!
!  nlagdump : lagplot file number (lagplot files are numbered sequentially.
!
!  msslag   : lagrangian mass of fluid element for lagplot data.
!-----------------------------------------------------------------------

INTEGER                                             :: nlagplt
INTEGER                                             :: ilagplt
INTEGER                                             :: nlagdump

REAL(KIND=double)                                   :: msslag

!-----------------------------------------------------------------------
!  Parameters for r-lagrangian plot edits
!-----------------------------------------------------------------------
!  nrlagplt : the unit number for rlagplot data.
!
!  irlagplt : r-lagrangian plot switch
!
!     0    : bypass writing rlagplot data.
!     ne 0 : write rlagplot data.
!
!  dmlag    : the lagrangian mass difference between adjacent lagrangian
!   points to be plotted.
!-----------------------------------------------------------------------

INTEGER                                             :: nrlagplt
INTEGER                                             :: irlagplt

REAL(KIND=double)                                   :: dmlag

!-----------------------------------------------------------------------
!  Parameters for energy-check plot edits
!-----------------------------------------------------------------------
!  n_eplt   : the unit number for e_chk data.
!
!  i_eplt   : energy-check plot switch
!
!     0    : bypass writing e_chk data.
!     ne 0 : write e_chk data.
!
!  dt_eplot : write to file e_chk.d every dt_eplot ms.
!-----------------------------------------------------------------------

INTEGER                                             :: n_eplt
INTEGER                                             :: i_eplt

REAL(KIND=double)                                   :: dt_eplot

!-----------------------------------------------------------------------
!  Parameters for pinch parameter plot edits
!-----------------------------------------------------------------------
!  r_pinch : the radius at which to evaluate the pinch parameters.
!
!  d_pinch : the density at which to evaluate the pinch parameters.
!-----------------------------------------------------------------------

REAL(KIND=double)                                   :: r_pinch
REAL(KIND=double)                                   :: d_pinch

!-----------------------------------------------------------------------
!  Parameters for M-D u-edits
!-----------------------------------------------------------------------
!
!  iedMDu    : switch for performing MD u-edits
!
!     0   : bypass
!     ne 0: edit when criteria satisfied
!
!  nedMDu    : edit counter
!  intedMDu  : edit every intedMDu cycles
!  n_editMDu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDu
INTEGER                                             :: nedMDu
INTEGER                                             :: intedMDu
INTEGER                                             :: n_editMDu

!-----------------------------------------------------------------------
!  Parameters for M-D v-edits
!-----------------------------------------------------------------------
!
!  iedMDv    : switch for performing MD v-edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  nedMDv    : edit counter
!  intedMDv  : edit every intedMDv cycles
!  n_editMDv : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDv
INTEGER                                             :: nedMDv
INTEGER                                             :: intedMDv
INTEGER                                             :: n_editMDv

!-----------------------------------------------------------------------
!  Parameters for M-D w-edits (z-velocity)
!-----------------------------------------------------------------------
!
!  iedMDw    : switch for performing MD w-edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  nedMDw    : edit counter
!  intedMDw  : edit every intedMDw cycles
!  n_editMDw : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDw
INTEGER                                             :: nedMDw
INTEGER                                             :: intedMDw
INTEGER                                             :: n_editMDw

!-----------------------------------------------------------------------
!  Parameters for M-D s-edits
!-----------------------------------------------------------------------
!
!  iedMDs    : switch for performing MD s-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDs    : edit counter
!  intedMDs  : edit every intedMDs cycles
!  n_editMDs : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDs
INTEGER                                             :: nedMDs
INTEGER                                             :: intedMDs
INTEGER                                             :: n_editMDs

!-----------------------------------------------------------------------
!  Parameters for M-D d-edits (density)
!-----------------------------------------------------------------------
!
!  iedMDd    : switch for performing MD d-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDd    : edit counter
!  intedMDd  : edit every intedMDd cycles
!  n_editMDd : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDd
INTEGER                                             :: nedMDd
INTEGER                                             :: intedMDd
INTEGER                                             :: n_editMDd

!-----------------------------------------------------------------------
!  Parameters for M-D e-edits (internal energy)
!-----------------------------------------------------------------------
!
!  iedMDe    : switch for performing MD e-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDe    : edit counter
!  intedMDe  : edit every intedMDe cycles
!  n_editMDe : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDe
INTEGER                                             :: nedMDe
INTEGER                                             :: intedMDe
INTEGER                                             :: n_editMDe

!-----------------------------------------------------------------------
!  Parameters for M-D e-edits (pressure)
!-----------------------------------------------------------------------
!
!  iedMDp    : switch for performing MD p-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDp    : edit counter
!  intedMDp  : edit every intedMDp cycles
!  n_editMDp : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDp
INTEGER                                             :: nedMDp
INTEGER                                             :: intedMDp
INTEGER                                             :: n_editMDp

!-----------------------------------------------------------------------
!  Parameters for M-D enu-edits (neutrino energy density)
!-----------------------------------------------------------------------
!
!  iedMDenu    : switch for performing MD enu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDenu    : edit counter
!  intedMDenu  : edit every intedMDenu cycles
!  n_editMDenu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDenu
INTEGER                                             :: nedMDenu
INTEGER                                             :: intedMDenu
INTEGER                                             :: n_editMDenu

!-----------------------------------------------------------------------
!  Parameters for M-D fnu-edits (neutrino flux)
!-----------------------------------------------------------------------
!
!  iedMDfnu    : switch for performing MD fnu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDfnu    : edit counter
!  intedMDfnu  : edit every intedMDfnu cycles
!  n_editMDfnu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDfnu
INTEGER                                             :: nedMDfnu
INTEGER                                             :: intedMDfnu
INTEGER                                             :: n_editMDfnu

!-----------------------------------------------------------------------
!  Parameters for M-D a-edits (atomic mass number)
!-----------------------------------------------------------------------
!
!  iedMDa    : switch for performing MD a-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDa    : edit counter
!  intedMDa  : edit every intedMDa cycles
!  n_editMDa : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDa
INTEGER                                             :: nedMDa
INTEGER                                             :: intedMDa
INTEGER                                             :: n_editMDa

!-----------------------------------------------------------------------
!  Parameters for M-D x-edits (parameter)
!-----------------------------------------------------------------------
!
!  iedMDx    : switch for performing MD a-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDx    : edit counter
!  intedMDx  : edit every intedMDx cycles
!  n_editMDx : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDx
INTEGER                                             :: nedMDx
INTEGER                                             :: intedMDx
INTEGER                                             :: n_editMDx

!-----------------------------------------------------------------------
!  Parameters for M-D ye-edits (electron fraction)
!-----------------------------------------------------------------------
!
!  iedMDye   : switch for performing MD ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDye    : edit counter
!  intedMDye  : edit every intedMDye cycles
!  n_editMDye : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDye
INTEGER                                             :: nedMDye
INTEGER                                             :: intedMDye
INTEGER                                             :: n_editMDye

!-----------------------------------------------------------------------
!  Parameters for M-D comp-edits (composition)
!-----------------------------------------------------------------------
!
!  iedMDcm   : switch for performing MD comp-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDcm    : edit counter
!  intedMDcm  : edit every intedMDcm cycles
!  n_editMDcm : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDcm
INTEGER                                             :: nedMDcm
INTEGER                                             :: intedMDcm
INTEGER                                             :: n_editMDcm

!-----------------------------------------------------------------------
!  Parameters for M-D dudt_nu-edits (neutrino energy deposition)
!-----------------------------------------------------------------------
!
!  iedMDnu   : switch for performing MD dudt_nu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnu    : edit counter
!  intedMDnu  : edit every intedMDnu cycles
!  n_editMDnu : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDnu
INTEGER                                             :: nedMDnu
INTEGER                                             :: intedMDnu
INTEGER                                             :: n_editMDnu

!-----------------------------------------------------------------------
!  Parameters for M-D dudt_nuc-edits (nuclear energy deposition)
!-----------------------------------------------------------------------
!
!  iedMDnc   : switch for performing MD dudt_nuc-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnc    : edit counter
!  intedMDnc  : edit every intedMDnc cycles
!  n_editMDnc : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDnc
INTEGER                                             :: nedMDnc
INTEGER                                             :: intedMDnc
INTEGER                                             :: n_editMDnc

!-----------------------------------------------------------------------
!  Parameters for M-D nl-edits (neutrino luminosity)
!-----------------------------------------------------------------------
!
!  iedMDnl   : switch for performing MD neutrino luminosity edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnl    : edit counter
!  intedMDnl  : edit every intedMDnl cycles
!  n_editMDnl : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDnl
INTEGER                                             :: nedMDnl
INTEGER                                             :: intedMDnl
INTEGER                                             :: n_editMDnl

!-----------------------------------------------------------------------
!  Parameters for M-D ne-edits (neutrino mean energies)
!-----------------------------------------------------------------------
!
!  iedMDne   : switch for performing MD neutrino rms energy edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDne    : edit counter
!  intedMDne  : edit every intedMDne cycles
!  n_editMDne : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDne
INTEGER                                             :: nedMDne
INTEGER                                             :: intedMDne
INTEGER                                             :: n_editMDne

!-----------------------------------------------------------------------
!  Parameters for M-D gx-edits
!-----------------------------------------------------------------------
!
!  iedMDgx   : switch for performing MD x-gravitational acceleration edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDgx    : edit counter
!  intedMDgx  : edit every intedMDgx cycles
!  n_editMDgx : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDgx
INTEGER                                             :: nedMDgx
INTEGER                                             :: intedMDgx
INTEGER                                             :: n_editMDgx

!-----------------------------------------------------------------------
!  Parameters for M-D gy-edits
!-----------------------------------------------------------------------
!
!  iedMDgy   : switch for performing MD y-gravitational acceleration edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDgy    : edit counter
!  intedMDgy  : edit every intedMDgy cycles
!  n_editMDgy : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDgy
INTEGER                                             :: nedMDgy
INTEGER                                             :: intedMDgy
INTEGER                                             :: n_editMDgy

!-----------------------------------------------------------------------
!  Parameters for M-D BVw-edits
!-----------------------------------------------------------------------
!
!  iedMDBVw   : switch for performing MD Brunt-Vaisala frequency edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDBVw    : edit counter
!  intedMDBVw  : edit every intedMDBVw cycles
!  n_editMDBVw : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDBVw
INTEGER                                             :: nedMDBVw
INTEGER                                             :: intedMDBVw
INTEGER                                             :: n_editMDBVw

!-----------------------------------------------------------------------
!  Parameters for M-D yl-edits
!-----------------------------------------------------------------------
!
!  iedMDyl   : switch for performing MD lwepton fraction edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDyl    : edit counter
!  intedMDyl  : edit every intedMDyl cycles
!  n_editMDyl : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: iedMDyl
INTEGER                                             :: nedMDyl
INTEGER                                             :: intedMDyl
INTEGER                                             :: n_editMDyl

!-----------------------------------------------------------------------
!  Time parameters for M-D edits
!-----------------------------------------------------------------------
!
!  i_editMD   : switch for performing MD edits at selected times
!  n_editMD   : model number counter
!  dt_MDedit1 : dump files for MD edits every dt_MDedit1 ms before bounce.
!  dt_MDedit2 : dump files for MD edits every dt_MDedit2 ms after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: i_editMD
INTEGER                                             :: n_editMD

REAL(KIND=double)                                   :: dt_MDedit1
REAL(KIND=double)                                   :: dt_MDedit2

!-----------------------------------------------------------------------
!  Cycle parameters for global edits
!-----------------------------------------------------------------------
!
!  ied_global_n  : switch for performing global using cycle criterion
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nned_global   : edit counter
!  inted_global  : edit every inted_global cycles
!  n_edit_global : model number counter
!-----------------------------------------------------------------------

INTEGER                                             :: ied_global_n
INTEGER                                             :: nned_global
INTEGER                                             :: inted_global
INTEGER                                             :: n_edit_global

!-----------------------------------------------------------------------
!  Time parameters for global edits
!-----------------------------------------------------------------------
!
!  ied_global_t  : switch for performing MD edits at selected times
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nted_global   : model number counter
!  dt_global_ed1 : dump files for global edits every dt_global_ed1 ms
!                   before bounce.
!  dt_global_ed2 : dump files for global edits every dt_global_ed2 ms
!                   after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: ied_global_t
INTEGER                                             :: nted_global

REAL(KIND=double)                                   :: dt_global_ed1
REAL(KIND=double)                                   :: dt_global_ed2

!-----------------------------------------------------------------------
!  Parameters for HDF M-D edits
!-----------------------------------------------------------------------
!
!  nd_HDFedit : the unit number for dumping files for HDF edits.
!
!  i_HDFedit : switch for dumping files for HDF edits.
!
!  n_HDFedit : HDF edit counter.
!
!     0    : bypass dumping files for HDF edits.
!     ne 0 : dump files for HDF edits.
!
!  dt_HDFedit1 : dump files for HDF edits every dt_HDFedit1 ms before bounce.
!  dt_HDFedit2 : dump files for HDF edits every dt_HDFedit1 ms after bounce.
!-----------------------------------------------------------------------

INTEGER                                             :: i_HDFedit
INTEGER                                             :: nd_HDFedit
INTEGER                                             :: n_HDFedit

REAL(KIND=double)                                   :: dt_HDFedit1
REAL(KIND=double)                                   :: dt_HDFedit2

    1 FORMAT (a128)
    2 FORMAT (1x,a128)
    3 FORMAT (/)
    5 FORMAT (10x,i10)
    7 FORMAT (' nrst=',i10)
    9 FORMAT (' nouttmp=',i7)
  101 FORMAT (a128)
  103 FORMAT (a6)
  105 FORMAT (1x,a128)
  111 FORMAT (20x,i10)
  113 FORMAT (1x,a6,14x,i10,42x,a16)
  121 FORMAT (10x,2i10)
  123 FORMAT (1x,a6,4x,2i10,42x,a16)
  131 FORMAT (20x,i10,5x,e15.8)
  133 FORMAT (1x,a6,14x,i10,5x,es15.8,22x,a16)
  141 FORMAT (35x,e15.8)
  143 FORMAT (1x,a6,29x,es15.8,22x,a16)
  151 FORMAT (10x,3i10)
  153 FORMAT (1x,a6,4x,3i10,32x,a16)
  161 FORMAT (10x,2i10,5x,e15.8)
  163 FORMAT (1x,a6,4x,2i10,5x,es15.8,22x,a16)
  211 FORMAT (10x,i10,5x,e15.8)
  213 FORMAT (1x,a6,4x,i10,5x,es15.8)
  251 FORMAT (10x,2i2,2i3,4(1x,e12.5))
  253 FORMAT (1x,a6,4x,2i2,2i3,4(1x,es12.5),a16)
  261 FORMAT (10x,2i2,2i3,2(1x,e12.5))
  263 FORMAT (1x,a6,4x,2i2,2i3,2(1x,es12.5),26x,a16)
  401 FORMAT (' The following data card was not recognized')
  403 FORMAT (1x,132a)
 1001 FORMAT (' Allocation problem for array ',a10,' in read_pack_edit_keys')
 2001 FORMAT (' Deallocation problem for array ',a10,' in read_pack_edit_keys')

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        Allocate arrays
!-----------------------------------------------------------------------

ALLOCATE (intedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intedn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedn      '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (idxedn(nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxedn    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (intdng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intdng    '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (nedng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedng     '; WRITE (nlog,1001) var_name; END IF
ALLOCATE (idxeng(40,nnu), STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxeng    '; WRITE (nlog,1001) var_name; END IF

!-----------------------------------------------------------------------
!
!                 \\\\\ INITIALIZE EDIT KEYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Initialize
!-----------------------------------------------------------------------

i_edit_data               = 0
d_edit_data               = zero

!-----------------------------------------------------------------------
!  General edit parameters
!-----------------------------------------------------------------------

iprint                    = 0
intprt                    = 0
nprt                      = 0
iprtbgn                   = 0
iflprt                    = 0
nmodel                    = 0
iplot                     = 0
intplf                    = 0
npltf                     = 0
nplotc                    = 0
nplote                    = 0
nplota                    = 0
nplott                    = 0

rhoprint                  = -1.d+00

!-----------------------------------------------------------------------
!  Parameters for configuration edits
!-----------------------------------------------------------------------

intedc                    = 9000000
nedc                      = 0
idxedc                    = 1

!-----------------------------------------------------------------------
!  Parameters for kinetic, internal, and gravitational energy edits
!-----------------------------------------------------------------------

intede                    = 9000000
nede                      = 0
idxede                    = 1

!-----------------------------------------------------------------------
!  Parameters for editmi edits
!-----------------------------------------------------------------------

intdmi                    = 9000000
nedmi                     = 0
idxemi                    = 1

!-----------------------------------------------------------------------
!  Parameters for mass average edits
!-----------------------------------------------------------------------

intdma                    = 9000000
nedma                     = 0
idxema                    = 1

!-----------------------------------------------------------------------
!  Parameters for hydro edits
!-----------------------------------------------------------------------

intedh                    = 9000000
nedh                      = 0
idxedh                    = 1

!-----------------------------------------------------------------------
!  Parameters for pressure - stress edits
!-----------------------------------------------------------------------

intdps                    = 9000000
nedps                     = 0
idxeps                    = 1

!-----------------------------------------------------------------------
!  Parameters for energy edits
!-----------------------------------------------------------------------

intedu                    = 9000000
nedu                      = 0
idxedu                    = 1

!-----------------------------------------------------------------------
!  Parameters for composition edits
!-----------------------------------------------------------------------

intedy                    = 9000000
nedy                      = 0
idxedy                    = 1

!-----------------------------------------------------------------------
!  Parameters for entropy and chemical potential edits
!-----------------------------------------------------------------------

intdsc                    = 9000000
nedsc                     = 0
idxesc                    = 1

!-----------------------------------------------------------------------
!  Parameters for diferential neutrino edit
!-----------------------------------------------------------------------

intedn                    = 9000000
nedn                      = 0
idxedn                    = 1
niedn                     = 9000000
neden                     = 0

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------

intdng                    = 9000000
nedng                     = 0
idxeng                    = 1

!-----------------------------------------------------------------------
!  Parameters editing at selected postbounce time intervals
!-----------------------------------------------------------------------

it_edit                   = 0
dt_edit                   = 1.d20

!-----------------------------------------------------------------------
!  Parameters for integral neutrino edits
!-----------------------------------------------------------------------

ilumplt                   = 0
nlumplt1                  = 0
nlumplt2                  = 0
nlum                      = 0
intlum                    = 0
ncylum                    = 0

r_lum                     = 1.d+08
d_lum                     = 1.d+12

!-----------------------------------------------------------------------
!  Parameters for rms neutrino energy plot edits
!-----------------------------------------------------------------------

ienuplt                   = 0
nenuplt1                  = 0
nenuplt2                  = 0
nenu                      = 0
intenu                    = 0
ncyenu                    = 0

r_e_rms                   = 1.d+08
d_e_rms                   = 1.d+12

!-----------------------------------------------------------------------
!  Parameters for selected radii plot edits
!-----------------------------------------------------------------------

irnuplt                   = 0
nrnuplt                   = 0
nrnu                      = 0
intrnu                    = 0
ncyrnu                    = 0

!-----------------------------------------------------------------------
!  Parameters for configuration plot edits
!-----------------------------------------------------------------------

ivarplt                   = 0
nvar                      = 0
nvarint                   = 0
nvarplt                   = 0
nvardump                  = 0

dtvarplot                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for comparison plot edits
!-----------------------------------------------------------------------

icomplt                   = 0
ncomplt                   = 0
ncomdump                  = 0

dtcomplot                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for bounday plot edits
!-----------------------------------------------------------------------

nplotinnerb               = 0
iplotinnerb               = 0
nplotouterb               = 0
iplotouterb               = 0
nplotlum                  = 0
iplotlum                  = 0
nplotshk                  = 0
iplotshk                  = 0
nplotcnv                  = 0
iplotcnv                  = 0

dtimeplot                 = 1.d+100
rinnerb                   = 1.d+100
routerb                   = 1.d+100
r_lumerms                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for nurad plot edits
!-----------------------------------------------------------------------

nplotnurad                = 0
iplotnurad                = 0

dtnuradplot               = 1.d+100
r_nurad                   = 1.d+100
rho_nurad                 = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for neutrino distribution plot edits
!-----------------------------------------------------------------------

nnudata                   = 0
inudata                   = 0

dtnudata                  = 1.d+100
r_nudata                  = 1.d+100
t_nudata                  = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for lagrangian plot edits
!-----------------------------------------------------------------------

nlagplt                   = 0
ilagplt                   = 0
nlagdump                  = 0

msslag                    = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for r-lagrangian plot edits
!-----------------------------------------------------------------------

nrlagplt                  = 0
irlagplt                  = 0

dmlag                     = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for energy-check plot edits
!-----------------------------------------------------------------------

n_eplt                    = 0
i_eplt                    = 0

dt_eplot                  = 1.d+100

!-----------------------------------------------------------------------
!  Parameters for M-D u-edits
!-----------------------------------------------------------------------
!
!  iedMDu    : switch for performing MD u-edits
!
!     0   : bypass
!     ne 0: edit when criteria satisfied
!
!  nedMDu    : edit counter
!  intedMDu  : edit every intedMDu cycles
!  n_editMDu : model number counter
!-----------------------------------------------------------------------

iedMDu                    = 0
nedMDu                    = 0
intedMDu                  = 200
n_editMDu                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D v-edits
!-----------------------------------------------------------------------
!
!  iedMDv    : switch for performing MD v-edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  nedMDv    : edit counter
!  intedMDv  : edit every intedMDv cycles
!  n_editMDv : model number counter
!-----------------------------------------------------------------------

iedMDv                    = 0
nedMDv                    = 0
intedMDv                  = 200
n_editMDv                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D w-edits (z-velocity)
!-----------------------------------------------------------------------
!
!  iedMDw    : switch for performing MD w-edits
!
!  0    : bypass
!  ne 0 : edit when criteria satisfied
!
!  nedMDw    : edit counter
!  intedMDw  : edit every intedMDw cycles
!  n_editMDw : model number counter
!-----------------------------------------------------------------------

iedMDw                    = 0
nedMDw                    = 0
intedMDw                  = 200
n_editMDw                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D s-edits
!-----------------------------------------------------------------------
!
!  iedMDs    : switch for performing MD s-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDs    : edit counter
!  intedMDs  : edit every intedMDs cycles
!  n_editMDs : model number counter
!-----------------------------------------------------------------------

iedMDs                    = 0
nedMDs                    = 0
intedMDs                  = 200
n_editMDs                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D d-edits (density)
!-----------------------------------------------------------------------
!
!  iedMDd    : switch for performing MD d-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDd    : edit counter
!  intedMDd  : edit every intedMDd cycles
!  n_editMDd : model number counter
!-----------------------------------------------------------------------

iedMDd                    = 0
nedMDd                    = 0
intedMDd                  = 200
n_editMDd                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D e-edits (internal energy)
!-----------------------------------------------------------------------
!
!  iedMDe    : switch for performing MD e-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDe    : edit counter
!  intedMDe  : edit every intedMDe cycles
!  n_editMDe : model number counter
!-----------------------------------------------------------------------

iedMDe                    = 0
nedMDe                    = 0
intedMDe                  = 200
n_editMDe                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D e-edits (pressure)
!-----------------------------------------------------------------------
!
!  iedMDp    : switch for performing MD p-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDp    : edit counter
!  intedMDp  : edit every intedMDp cycles
!  n_editMDp : model number counter
!-----------------------------------------------------------------------

iedMDp                    = 0
nedMDp                    = 0
intedMDp                  = 200
n_editMDp                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D enu-edits (neutrino energy density)
!-----------------------------------------------------------------------
!
!  iedMDenu    : switch for performing MD enu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDenu    : edit counter
!  intedMDenu  : edit every intedMDenu cycles
!  n_editMDenu : model number counter
!-----------------------------------------------------------------------

iedMDenu                  = 0
nedMDenu                  = 0
intedMDenu                = 200
n_editMDenu               = 0

!-----------------------------------------------------------------------
!  Parameters for M-D fnu-edits (neutrino flux)
!-----------------------------------------------------------------------
!
!  iedMDfnu    : switch for performing MD fnu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDfnu    : edit counter
!  intedMDfnu  : edit every intedMDfnu cycles
!  n_editMDfnu : model number counter
!-----------------------------------------------------------------------

iedMDfnu                  = 0
nedMDfnu                  = 0
intedMDfnu                = 200
n_editMDfnu               = 0

!-----------------------------------------------------------------------
!  Parameters for M-D a-edits (atomic mass number)
!-----------------------------------------------------------------------
!
!  iedMDa    : switch for performing MD a-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDa    : edit counter
!  intedMDa  : edit every intedMDa cycles
!  n_editMDa : model number counter
!-----------------------------------------------------------------------

iedMDa                    = 0
nedMDa                    = 0
intedMDa                  = 200
n_editMDa                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D x-edits (parameter)
!-----------------------------------------------------------------------
!
!  iedMDx    : switch for performing MD a-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDx    : edit counter
!  intedMDx  : edit every intedMDx cycles
!  n_editMDx : model number counter
!-----------------------------------------------------------------------

iedMDx                    = 0
nedMDx                    = 0
intedMDx                  = 200
n_editMDx                 = 0

!-----------------------------------------------------------------------
!  Parameters for M-D ye-edits
!-----------------------------------------------------------------------
!
!  iedMDye   : switch for performing MD ye-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDye    : edit counter
!  intedMDye  : edit every intedMDye cycles
!  n_editMDye : model number counter
!-----------------------------------------------------------------------

iedMDye                   = 0
nedMDye                   = 0
intedMDye                 = 200
n_editMDye                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D comp-edits (composition)
!-----------------------------------------------------------------------
!
!  iedMDcm   : switch for performing MD comp-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDcm    : edit counter
!  intedMDcm  : edit every intedMDcm cycles
!  n_editMDcm : model number counter
!-----------------------------------------------------------------------

iedMDcm                   = 0
nedMDcm                   = 0
intedMDcm                 = 200
n_editMDcm                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D dudt_nu-edits (neutrino energy deposition)
!-----------------------------------------------------------------------
!
!  iedMDnu   : switch for performing MD dudt_nu-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnu    : edit counter
!  intedMDnu  : edit every intedMDnu cycles
!  n_editMDnu : model number counter
!-----------------------------------------------------------------------

iedMDnu                   = 0
nedMDnu                   = 0
intedMDnu                 = 200
n_editMDnu                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D dudt_nuc-edits (nuclear energy deposition)
!-----------------------------------------------------------------------
!
!  iedMDnc   : switch for performing MD dudt_nuc-edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnc    : edit counter
!  intedMDnc  : edit every intedMDnc cycles
!  n_editMDnc : model number counter
!-----------------------------------------------------------------------

iedMDnc                   = 0
nedMDnc                   = 0
intedMDnc                 = 200
n_editMDnc                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D nl-edits
!-----------------------------------------------------------------------
!
!  iedMDnl   : switch for performing MD neutrino luminosity edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDnl    : edit counter
!  intedMDnl  : edit every intedMDnl cycles
!  n_editMDnl : model number counter
!-----------------------------------------------------------------------

iedMDnl                   = 0
nedMDnl                   = 0
intedMDnl                 = 200
n_editMDnl                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D ne-edits
!-----------------------------------------------------------------------
!
!  iedMDne   : switch for performing MD neutrino rms energy edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDne    : edit counter
!  intedMDne  : edit every intedMDne cycles
!  n_editMDne : model number counter
!-----------------------------------------------------------------------

iedMDne                   = 0
nedMDne                   = 0
intedMDne                 = 200
n_editMDne                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D gx-edits
!-----------------------------------------------------------------------
!
!  iedMDgx   : switch for performing MD x-gravitational acceleration edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDgx    : edit counter
!  intedMDgx  : edit every intedMDgx cycles
!  n_editMDgx : model number counter
!-----------------------------------------------------------------------

iedMDgx                   = 0
nedMDgx                   = 0
intedMDgx                 = 200
n_editMDgx                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D gy-edits
!-----------------------------------------------------------------------
!
!  iedMDgy   : switch for performing MD y-gravitational acceleration edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDgy    : edit counter
!  intedMDgy  : edit every intedMDgy cycles
!  n_editMDgy : model number counter
!-----------------------------------------------------------------------

iedMDgy                   = 0
nedMDgy                   = 0
intedMDgy                 = 200
n_editMDgy                = 0

!-----------------------------------------------------------------------
!  Parameters for M-D BVw-edits
!-----------------------------------------------------------------------
!
!  iedMDBVw   : switch for performing MD Brunt-Vaisala frequency edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDBVw    : edit counter
!  intedMDBVw  : edit every intedMDBVw cycles
!  n_editMDBVw : model number counter
!-----------------------------------------------------------------------

iedMDBVw                  = 0
nedMDBVw                  = 0
intedMDBVw                = 200
n_editMDBVw               = 0

!-----------------------------------------------------------------------
!  Parameters for M-D yl-edits
!-----------------------------------------------------------------------
!
!  iedMDyl   : switch for performing MD lwepton fraction edits
!
!     0    : bypass
!     ne 0 : edit when criteria satisfied
!
!  nedMDyl    : edit counter
!  intedMDyl  : edit every intedMDyl cycles
!  n_editMDyl : model number counter
!-----------------------------------------------------------------------

iedMDyl                   = 0
nedMDyl                   = 0
intedMDyl                 = 0
n_editMDyl                = 0

!-----------------------------------------------------------------------
!  Time parameters for M-D edits
!-----------------------------------------------------------------------

i_editMD                  = 0
n_editMD                  = 0
dt_MDedit1                = 1.d+20
dt_MDedit2                = 1.d+20

!-----------------------------------------------------------------------
!  Cycle parameters for global edits
!-----------------------------------------------------------------------

ied_global_n              = 0
nned_global               = 0
inted_global              = 200
n_edit_global             = 0

!-----------------------------------------------------------------------
!  Time parameters for global edits
!-----------------------------------------------------------------------

ied_global_t              = 0
nted_global               = 0
dt_global_ed1             = 1.d+20
dt_global_ed2             = 1.d+20

!-----------------------------------------------------------------------
!  Parameters for multidimensional edits
!-----------------------------------------------------------------------

i_HDFedit                  = 0
nd_HDFedit                 = 31
n_HDFedit                  = 0
dt_HDFedit1                = 1.d+20
dt_HDFedit2                = 1.d+20

!-----------------------------------------------------------------------
!  Parameters for temperary restart dumps
!-----------------------------------------------------------------------

intrst                    = 9000000
nnrst                     = 0
nrstd1                    = 0
nrstd2                    = 0
nrstfl                    = 0
nouttmp                   = 0

!-----------------------------------------------------------------------
!  Parameters for permanent restart dumps
!-----------------------------------------------------------------------

noutpmt                   = 0
intprm                    = 9000000
nprm                      = 0
ncyrst                    = -1
ncychg                    = 0
irstbgn                   = 0
noutfl                    = 0

!-----------------------------------------------------------------------
!  Parameters for HDF edits
!-----------------------------------------------------------------------

i_HDFedit                  = 0
nd_HDFedit                 = 31
n_HDFedit                  = 0

dt_HDFedit1                = 1.d100
dt_HDFedit2                = 1.d100

!-----------------------------------------------------------------------
!
!                   \\\\\ READ EDIT KEYS /////
!
!-----------------------------------------------------------------------

REWIND (nread)

READ: DO

!-----------------------------------------------------------------------
!  Read line
!-----------------------------------------------------------------------

  READ (nread,101,END=5000) line

  type               = line(1:6)
  name               = line(73:88)

!-----------------------------------------------------------------------
!  Ignore comment cards in the data stream
!-----------------------------------------------------------------------

  IF ( type == 'cccccc'  .or.  type(1:1) == '!'  .or.  type == '      ') THEN
    IF ( iskip == 0 ) WRITE (nprint,105) line
    CYCLE
  END IF ! type = 'cccccc' or type = '!'

!-----------------------------------------------------------------------
!
!                     \\\\\ EDIT ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  iflprt
!-----------------------------------------------------------------------

  IF ( type == 'iflprt' ) THEN
    READ (line ,111) iflprt
    IF ( iskip == 0 ) WRITE (nprint,113) type,iflprt,name
    CYCLE
  END IF ! type = 'iflprt'

!-----------------------------------------------------------------------
!  noutfl
!-----------------------------------------------------------------------

  IF ( type == 'noutfl' ) THEN
    READ (line ,111) noutfl
    IF ( iskip == 0 ) WRITE (nprint,113) type,noutfl,name
    CYCLE
  END IF ! type = 'noutfl'

!-----------------------------------------------------------------------
!  output
!-----------------------------------------------------------------------

  IF ( type == 'output' ) THEN
    READ (line ,111) int

    IF ( name == 'iprint  ') THEN
      iprint         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,iprint,name
      CYCLE
    END IF ! name = 'iprint  '

    IF ( name == 'nrstd1  ') THEN
      nrstd1         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrstd1,name
      CYCLE
    END IF ! name = 'nrstd1  '

    IF ( name == 'nrstd2  ') THEN
      nrstd2         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrstd2,name
      CYCLE
    END IF ! name = 'nrstd2  '

    IF ( name == 'noutpmt ') THEN
      noutpmt        = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,noutpmt,name
      CYCLE
    END IF ! name = 'noutpmt '

    IF ( name == 'nouttmp ') THEN
      nouttmp        = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nouttmp,name
      CYCLE
    END IF ! name = 'nouttmp '

  END IF ! type = 'output'

!-----------------------------------------------------------------------
!  plot
!-----------------------------------------------------------------------

  IF ( type == 'plot  ' ) THEN
    READ (line ,111) int

    IF ( name == 'iplot   ') THEN
      iplot          = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,iplot,name
      CYCLE
    END IF ! name = 'iplot   '

    IF ( name == 'nplotc  ') THEN
      nplotc         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nplotc,name
      CYCLE
    END IF ! name = 'nplotc  '

    IF ( name == 'nplote  ') THEN
      nplote         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nplote,name
      CYCLE
    END IF ! name = 'nplote  '

    IF ( name == 'nplota  ') THEN
      nplota         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nplota,name
      CYCLE
    END IF ! name = 'nplota  '

    IF ( name == 'nplott  ') THEN
      nplott         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nplott,name
      CYCLE
    END IF ! name = 'nplott  '

    IF ( name == 'intplf  ') THEN
      intplf         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intplf,name
      CYCLE
    END IF ! name = 'intplf  '

    IF ( name == 'npltf   ') THEN
      npltf          = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,npltf,name
      CYCLE
    END IF ! name = 'npltf   '

  END IF ! type = 'plot  '

!-----------------------------------------------------------------------
!  model number
!-----------------------------------------------------------------------

  IF ( type == 'nmodel' ) THEN
    READ (line ,111) nmodel
    IF ( iskip == 0 ) WRITE (nprint,113) type,nmodel,name
    CYCLE
  END IF ! type = 'nmodel'

!-----------------------------------------------------------------------
!  rhoprt
!-----------------------------------------------------------------------
                                                                     !
  IF ( type == 'rhoprt' ) THEN
    READ (line ,211) i,rl
    rhoprint(i)      = rl
    IF ( iskip == 0 ) WRITE (nprint,213) type,i,rhoprint(i)
    CYCLE
  END IF ! type = 'rhoprt'

!-----------------------------------------------------------------------
!  intedc
!-----------------------------------------------------------------------
 
  IF ( type == 'intedc' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intedc  ') THEN
      intedc(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intedc(i),name
      CYCLE
    END IF ! name = 'intedc  '

    IF ( name == 'nedc    ') THEN
      nedc(i)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedc(i),name
      CYCLE
    END IF ! name = 'nedc    '

    IF ( name == 'idxedc  ') THEN
      idxedc(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxedc(i),name
      CYCLE
    END IF ! name = 'idxedc  '

  END IF

!-----------------------------------------------------------------------
!  intede
!-----------------------------------------------------------------------

  IF ( type == 'intede' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intede  ') THEN
      intede(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intede(i),name
      CYCLE
    END IF ! name = 'intede  '

    IF ( name == 'nede    ') THEN
      nede(i)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nede(i),name
      CYCLE
    END IF ! name = 'nede    '

    IF ( name == 'idxede  ') THEN
      idxede(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxede(i),name
      CYCLE
    END IF ! name = 'idxede  '

  END IF

!-----------------------------------------------------------------------
!  intdmi
!-----------------------------------------------------------------------

  IF ( type == 'intdmi' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intdmi  ') THEN
      intdmi(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intdmi(i),name
      CYCLE
    END IF ! name = 'intdmi  '

    IF ( name == 'nedmi   ') THEN
      nedmi(i)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedmi(i),name
      CYCLE
    END IF ! name = 'nedmi   '

    IF ( name == 'idxemi  ') THEN
      idxemi(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxemi(i),name
      CYCLE
    END IF ! name = 'idxemi  '

  END IF ! type = 'intdmi'

!-----------------------------------------------------------------------
!  intdma
!-----------------------------------------------------------------------

  IF ( type == 'intdma' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intdma  ') THEN
      intdma(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intdma(i),name
      CYCLE
    END IF ! name = 'intdma  '

    IF ( name == 'nedma   ') THEN
      nedma(i)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedma(i),name
      CYCLE
    END IF ! name = 'nedma   '

    IF ( name == 'idxema  ') THEN
      idxema(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxema(i),name
      CYCLE
    END IF ! name = 'idxema  '

  END IF ! type = 'intdma'

!-----------------------------------------------------------------------
!  intedh
!-----------------------------------------------------------------------

  IF ( type == 'intedh' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intedh  ') THEN
      intedh(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intedh(i),name
      CYCLE
    END IF ! name = 'intedh  '

    IF ( name == 'nedh    ') THEN
      nedh(i)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedh(i),name
      CYCLE
    END IF ! name = 'nedh    '

    IF ( name == 'idxedh  ') THEN
      idxedh(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxedh(i),name
      CYCLE
    END IF ! name = 'idxedh  '

  END IF

!-----------------------------------------------------------------------
!  intdps
!-----------------------------------------------------------------------

  IF ( type == 'intdps' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intdps  ') THEN
      intdps(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intdps(i),name
      CYCLE
    END IF ! name = 'intdps  '

    IF ( name == 'nedps   ') THEN
      nedps(i)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedps(i),name
      CYCLE
    END IF ! name = 'nedps   '

    IF ( name == 'idxeps  ') THEN
      idxeps(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxeps(i),name
      CYCLE
    END IF ! name = 'idxeps  '

  END IF ! type = 'intdps'

!-----------------------------------------------------------------------
!  intedu
!-----------------------------------------------------------------------

  IF ( type == 'intedu' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intedu  ') THEN
      intedu(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intedu(i),name
      CYCLE
    END IF ! name = 'intedu  '

    IF ( name == 'nedu    ') THEN
      nedu(i)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedu(i),name
      CYCLE
    END IF ! name = 'nedu    '

    IF ( name == 'idxedu  ') THEN
      idxedu(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxedu(i),name
      CYCLE
    END IF ! name = 'idxedu  '

  END IF ! type = 'intedu'

!-----------------------------------------------------------------------
!  intedy
!-----------------------------------------------------------------------

  IF ( type == 'intedy' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intedy  ') THEN
      intedy(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intedy(i),name
      CYCLE
    END IF ! name = 'intedy  '

    IF ( name == 'nedy    ') THEN
      nedy(i)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedy(i),name
      CYCLE
    END IF ! name = 'nedy    '

    IF ( name == 'idxedy  ') THEN
      idxedy(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxedy(i),name
      CYCLE
    END IF ! name = 'idxedy  '

  END IF ! type = 'intedy'

!-----------------------------------------------------------------------
!  intdsc
!-----------------------------------------------------------------------

  IF ( type == 'intdsc' ) THEN
    READ (line ,121) i,int

    IF ( name == 'intdsc  ') THEN
      intdsc(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,intdsc(i),name
      CYCLE
    END IF ! name = 'intdsc  '

    IF ( name == 'nedsc   ') THEN
      nedsc(i)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,nedsc(i),name
      CYCLE
    END IF ! name = 'nedsc   '

    IF ( name == 'idxesc  ') THEN
      idxesc(i)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,i,idxesc(i),name
      CYCLE
    END IF ! name = 'idxesc  '

  END IF ! type = 'intdsc'

!-----------------------------------------------------------------------
!  intedn
!-----------------------------------------------------------------------

  IF ( type == 'intedn' ) THEN

    READ (line ,121) n,int

    IF ( name == 'intedn  ') THEN
     intedn(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,intedn(n),name
      CYCLE
    END IF ! name = 'intedn  '

    IF ( name == 'nedn    ') THEN
      nedn(n)        = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,nedn(n),name
      CYCLE
    END IF ! name = 'nedn    '

    IF ( name == 'idxedn  ') THEN
      idxedn(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,idxedn(n),name
      CYCLE
    END IF ! name = 'idxedn  '

    IF ( name == 'niedn   ') THEN
      niedn(n)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,niedn(n),name
      CYCLE
    END IF ! name = 'niedn   '

    IF ( name == 'neden   ') THEN
      neden(n)       = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,neden(n),name
      CYCLE
    END IF ! name = 'neden   '

  END IF ! type = 'intedn'

!-----------------------------------------------------------------------
!  intdng
!-----------------------------------------------------------------------

  IF ( type == 'intdng' ) THEN
    READ (line ,151) i,int,n

    IF ( name == 'intdng  ') THEN
      intdng(i,n)    = int
      IF ( iskip == 0 ) WRITE (nprint,153) type,i,intdng(i,n),n,name
      CYCLE
    END IF ! name = 'intdng  '

    IF ( name == 'nedng   ') THEN
      nedng(i,n)     = int
      IF ( iskip == 0 ) WRITE (nprint,153) type,i,nedng(i,n),n,name
      CYCLE
    END IF ! name = 'nedng   '

    IF ( name == 'idxeng  ') THEN
      idxeng(i,n)    = int
      IF ( iskip == 0 ) WRITE (nprint,153) type,i,idxeng(i,n),n,name
      CYCLE
    END IF ! name = 'idxeng  '
 
  END IF ! type = 'intdng'

!-----------------------------------------------------------------------
!  t_edit
!-----------------------------------------------------------------------

  IF ( type == 't_edit' ) THEN
    READ (line ,131) it_edit,dt_edit
    IF ( iskip == 0 ) WRITE (nprint,133) type,it_edit,dt_edit,name
  END IF ! type = 't_edit'

!-----------------------------------------------------------------------
!  intprt
!-----------------------------------------------------------------------

  IF ( type == 'intprt' ) THEN
    READ (line ,111) int

    IF ( name == 'intprt  ') THEN
      intprt     = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intprt,name
      CYCLE
    END IF ! name = 'intprt  '

    IF ( name == 'nprt    ') THEN
      nprt       = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nprt,name
      CYCLE
    END IF ! name = 'nprt    '

    IF ( name == 'iprtbgn ') THEN
      iprtbgn    = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,iprtbgn,name
      CYCLE
    END IF ! name = 'iprtbgn '

  END IF ! type = 'intprt'

!-----------------------------------------------------------------------
!
!                    \\\\\ RESTART ARAAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  intrst
!-----------------------------------------------------------------------

  IF ( type == 'intrst' ) THEN
    READ (line ,111) int

    IF ( name == 'intrst  ') THEN
      intrst         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intrst,name
      CYCLE
    END IF ! name = 'intrst  '

    IF ( name == 'nnrst   ') THEN
      nnrst          = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nnrst,name
      CYCLE
    END IF ! name = 'nnrst   '

    IF ( name == 'irstbgn ') THEN
      irstbgn        = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,irstbgn,name
      CYCLE
    END IF ! name = 'irstbgn '

    IF ( name == 'nrstfl  ') THEN
      nrstfl         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrstfl,name
      CYCLE
    END IF ! name = 'nrstfl  '

    IF ( name == 'ncychg  ') THEN
      ncychg         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,ncychg,name
      CYCLE
    END IF ! name = 'ncychg  '

  END IF ! type = 'intrst'

!-----------------------------------------------------------------------
!  intprm
!-----------------------------------------------------------------------

  IF ( type == 'intprm' ) THEN
    READ (line ,111) int

    IF ( name == 'intprm  ') THEN
      intprm         = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,intprm,name
      CYCLE
    END IF ! name = 'intprm  '

    IF ( name == 'nprm    ') THEN
      nprm           = int
      IF ( iskip == 0 ) WRITE (nprint,113) type,nprm,name
      CYCLE
    END IF ! name = 'nprm    '
  END IF ! type = 'intprm'

!-----------------------------------------------------------------------
!  ncyrst
!-----------------------------------------------------------------------

  IF ( type == 'ncyrst' ) THEN
    READ (line ,121) i,int
    ncyrst(i)       = int
    IF ( iskip == 0 ) WRITE (nprint,123) type,i,ncyrst(i),name
    CYCLE
  END IF ! type = 'ncyrst'

!-----------------------------------------------------------------------
!
!                      \\\\\ PLOT ARRAYS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  lumplt
!-----------------------------------------------------------------------

  IF ( type == 'lumplt' ) THEN

    IF ( name == 'ilumplt ') THEN
      READ (line ,111) ilumplt
      IF ( iskip == 0 ) WRITE (nprint,113) type,ilumplt,name
      CYCLE
    END IF ! name = 'ilumplt '

    IF ( name == 'nlumplt ') THEN
      READ (line ,121) nlumplt1,nlumplt2
      IF ( iskip == 0 ) WRITE (nprint,123) type,nlumplt1,nlumplt2,name
      CYCLE
    END IF ! name = 'nlumplt '

    IF ( name == 'r_lum    ') THEN
      READ (line ,141) r_lum
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_lum,name
      CYCLE
    END IF ! name = 'r_lum    '

    IF ( name == 'd_lum  ') THEN
      READ (line ,141) d_lum
      IF ( iskip == 0 ) WRITE (nprint,143) type,d_lum,name
      CYCLE
    END IF ! name = 'd_lum  '

    IF ( name == 'nlum    ') THEN
      READ (line ,111) nlum
      IF ( iskip == 0 ) WRITE (nprint,113) type,nlum,name
      CYCLE
    END IF ! name = 'nlum    '

    IF ( name == 'intlum  ') THEN
      READ (line ,121) n,int
      intlum(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,intlum(n),name
      CYCLE
    END IF ! name = 'intlum  '

    IF ( name == 'ncylum  ') THEN
      READ (line ,121) n,int
      ncylum(n)  = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,ncylum(n),name
      CYCLE
    END IF ! name = 'ncylum  '

  END IF

!-----------------------------------------------------------------------
!  enuplt
!-----------------------------------------------------------------------

  IF ( type == 'enuplt' ) THEN

    IF ( name == 'ienuplt ') THEN
      READ (line ,111) ienuplt
      IF ( iskip == 0 ) WRITE (nprint,113) type,ienuplt,name
      CYCLE
    END IF ! name = 'ienuplt '

    IF ( name == 'nenuplt ') THEN
      READ (line ,121) nenuplt1,nenuplt2
      IF ( iskip == 0 ) WRITE (nprint,123) type,nenuplt1,nenuplt2, name
      CYCLE
    END IF ! name = 'nenuplt '

    IF ( name == 'r_e_rms    ') THEN
      READ (line ,141) r_e_rms
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_e_rms,name
      CYCLE
    END IF ! name = 'r_e_rms    '

    IF ( name == 'd_e_rms  ') THEN
      READ (line ,141) d_e_rms
      IF ( iskip == 0 ) WRITE (nprint,143) type,d_e_rms,name
      CYCLE
    END IF ! name = 'd_e_rms  '

    IF ( name == 'nenu    ') THEN
      READ (line ,111) nenu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nenu,name
      CYCLE
    END IF ! name = 'nenu    '

    IF ( name == 'intenu  ') THEN
      READ (line ,121) n,int
      intenu(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,intenu(n),name
      CYCLE
    END IF ! name = 'intenu  '

    IF ( name == 'ncyenu  ') THEN
      READ (line ,121) n,int
      ncyenu(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,ncyenu(n),name
      CYCLE
    END IF ! name = 'ncyenu  '

  END IF ! type = 'enuplt'

!-----------------------------------------------------------------------
!  rnuplt
!-----------------------------------------------------------------------

  IF ( type == 'rnuplt' ) THEN

    IF ( name == 'irnuplt ') THEN
      READ (line ,111) irnuplt
      IF ( iskip == 0 ) WRITE (nprint,113) type,irnuplt,name
      CYCLE
    END IF ! name = 'irnuplt '

    IF ( name == 'nrnuplt ') THEN
      READ (line ,111) nrnuplt
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrnuplt,name
      CYCLE
    END IF ! name = 'nrnuplt '

    IF ( name == 'nrnu    ') THEN
      READ (line ,111) nrnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nrnu,name
      CYCLE

    END IF ! name = 'nrnu    '
     IF ( name == 'intrnu  ') THEN
      READ (line ,121) n,int
      intrnu(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,intrnu(n),name
      CYCLE
    END IF ! name = 'intrnu  '

    IF ( name == 'ncyrnu  ') THEN
      READ (line ,121) n,int
      ncyrnu(n)      = int
      IF ( iskip == 0 ) WRITE (nprint,123) type,n,ncyrnu(n),name
      CYCLE
    END IF ! name = 'ncyrnu  '

  END IF ! type = 'rnuplt'

!-----------------------------------------------------------------------
!  bnuplt
!-----------------------------------------------------------------------

  IF ( type == 'bnuplt' ) THEN

    IF ( name == 'dtimeplot') THEN
      READ (line ,141) dtimeplot
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtimeplot,name
      CYCLE
    END IF ! name = 'dtimeplot'

    IF ( name == 'innerb') THEN
      READ (line ,161) nplotinnerb,iplotinnerb,rinnerb
      IF ( iskip == 0 ) WRITE (nprint,163) type,nplotinnerb,iplotinnerb,rinnerb,name
      CYCLE
    END IF ! name = 'innerb'

    IF ( name == 'outerb') THEN
      READ (line ,161) nplotouterb,iplotouterb,routerb
      IF ( iskip == 0 ) WRITE (nprint,163) type,nplotouterb,iplotouterb,routerb,name
      CYCLE
    END IF ! name = 'outerb'

    IF ( name == 'lumermsplt') THEN
      READ (line ,161) nplotlum,iplotlum,r_lumerms
      IF ( iskip == 0 ) WRITE (nprint,163) type,nplotlum,iplotlum,r_lumerms,name
      CYCLE
    END IF ! name = 'lumermsplt'

    IF ( name == 'shkplt') THEN
      READ (line ,121) nplotshk,iplotshk
      IF ( iskip == 0 ) WRITE (nprint,123) type,nplotshk,iplotshk,name
      CYCLE 
    END IF ! name = 'shkplt'

    IF ( name == 'cnvplt') THEN
      READ (line ,121) nplotcnv,iplotcnv
      IF ( iskip == 0 ) WRITE (nprint,123) type,nplotcnv,iplotcnv,name
      CYCLE
    END IF ! name = 'cnvplt'

  END IF ! type = 'bnuplt'

!-----------------------------------------------------------------------
!  varplt
!-----------------------------------------------------------------------

  IF ( type == 'varplt' ) THEN

    IF ( name == 'ivarplt') THEN
      READ (line ,121) nvarplt,ivarplt
      IF ( iskip == 0 ) WRITE (nprint,123) type,ivarplt,nvarplt,name
      CYCLE
    END IF !  name = 'ivarplt'

    IF ( name == 'nvarint') THEN
      READ (line ,151) nvar,nvarint,nvardump
      IF ( iskip == 0 ) WRITE (nprint,153) type,nvar,nvarint,nvardump,name
      CYCLE
    END IF ! name = 'nvarint'

    IF ( name == 'dtvarplot' ) THEN
      READ (line ,141) dtvarplot
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtvarplot,name
      CYCLE
    END IF ! name = 'dtvarplot'

  END IF ! type = 'varplt'

!-----------------------------------------------------------------------
!  complt
!-----------------------------------------------------------------------

  IF ( type == 'complt' ) THEN

    IF ( name == 'icomplt') THEN
      READ (line ,121) ncomplt,icomplt
      IF ( iskip == 0 ) WRITE (nprint,123) type,icomplt,ncomplt,name
      CYCLE
    END IF !  name = 'icomplt'

    IF ( name == 'dtcomplot' ) THEN
      READ (line ,131) ncomdump,dtcomplot
      IF ( iskip == 0 ) WRITE (nprint,133) type,ncomdump,dtcomplot,name
      CYCLE
    END IF ! name = 'dtcomplot'

  END IF ! type = 'varplt'

!-----------------------------------------------------------------------
!  nurad
!-----------------------------------------------------------------------

  IF ( type == 'nurad ' ) THEN

    IF ( name == 'inuplt') THEN
      READ (line ,121) nplotnurad,iplotnurad
      IF ( iskip == 0 ) WRITE (nprint,123) type,nplotnurad,iplotnurad,name
      CYCLE
    END IF ! name = 'inuplt'

    IF ( name == 'dtnuradplot') THEN
      READ (line ,141) dtnuradplot
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnuradplot,name
      CYCLE
    END IF ! name = 'dtnuradplot'

    IF ( name == 'r_nurad') THEN
      READ (line ,141) r_nurad
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_nurad,name
      CYCLE
    END IF ! name = 'r_nurad'

    IF ( name == 'rho_nurad') THEN
      READ (line ,141) rho_nurad
      IF ( iskip == 0 ) WRITE (nprint,143) type,rho_nurad,name
      CYCLE
    END IF ! name = 'rho_nurad'

  END IF ! type = 'nurad '

!-----------------------------------------------------------------------
!  nudata
!-----------------------------------------------------------------------

  IF ( type == 'nudata' ) THEN

    IF ( name == 'inudata') THEN
      READ (line ,121) nnudata,inudata
      IF ( iskip == 0 ) WRITE (nprint,123) type,nnudata,inudata,name
      CYCLE
    END IF ! name = 'inudata'

    IF ( name == 'dtnudata') THEN
      READ (line ,141) dtnudata
      IF ( iskip == 0 ) WRITE (nprint,143) type,dtnudata,name
      CYCLE
    END IF ! name = 'dtnudata'

    IF ( name == 'r_nudata') THEN
      READ (line ,141) r_nudata
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_nudata,name
      CYCLE
    END IF ! name = 'r_nudata'

    IF ( name == 't_nudata') THEN
      READ (line ,141) t_nudata
      IF ( iskip == 0 ) WRITE (nprint,143) type,t_nudata,name
      CYCLE
    END IF ! name = 't_nudata'

  END IF ! type = 'nudata'

!-----------------------------------------------------------------------
!  lagplt
!-----------------------------------------------------------------------

  IF ( type == 'lagplt' ) THEN

    IF ( name == 'ilagplt') THEN
      READ (line ,151) nlagplt,ilagplt,nlagdump
      IF ( iskip == 0 ) WRITE (nprint,153) type,nlagplt,ilagplt,nlagdump,name
      CYCLE
    END IF ! name = 'ilagplt'

    IF ( name == 'msslag') THEN
      READ (line ,141) msslag
      IF ( iskip == 0 ) WRITE (nprint,143) type,msslag,name
       CYCLE
    END IF ! name = 'msslag'

  END IF ! type = 'lagplt' 

!-----------------------------------------------------------------------
!  rlagpt
!-----------------------------------------------------------------------

  IF ( type == 'rlagpt' ) THEN

    IF ( name == 'irlagplt') THEN
      READ (line ,121) nrlagplt,irlagplt
      IF ( iskip == 0 ) WRITE (nprint,123) type,nrlagplt,irlagplt,name
      CYCLE
    END IF ! name = 'irlagplt'

    IF ( name == 'dmlag') THEN
      READ (line ,141) dmlag
      IF ( iskip == 0 ) WRITE (nprint,143) type,dmlag,name
       CYCLE
    END IF ! name = 'dmlag'

  END IF ! type = 'rlagpt' 

!-----------------------------------------------------------------------
!  e_plt
!-----------------------------------------------------------------------

  IF ( type == 'e_plt ' ) THEN

    IF ( name == 'i_eplt') THEN
      READ (line ,121) n_eplt,i_eplt
      IF ( iskip == 0 ) WRITE (nprint,123) type,n_eplt,i_eplt,name
      CYCLE
    END IF ! name = 'i_eplt'

    IF ( name == 'dt_eplot') THEN
      READ (line ,141) dt_eplot
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_eplot,name
      CYCLE
    END IF ! name = 'dt_eplot'

  END IF ! type = 'e_plt ' 

!-----------------------------------------------------------------------
!  pinch
!-----------------------------------------------------------------------

  IF ( type == 'pinch ' ) THEN

    IF ( name == 'r_pinch') THEN
      READ (line ,141) r_pinch
      IF ( iskip == 0 ) WRITE (nprint,143) type,r_pinch,name
      CYCLE
    END IF ! name = 'r_pinch'

    IF ( name == 'd_pinch') THEN
      READ (line ,141) d_pinch
      IF ( iskip == 0 ) WRITE (nprint,143) type,d_pinch,name
      CYCLE
    END IF ! name = 'd_pinch'

  END IF ! type = 'pinch ' 

!-----------------------------------------------------------------------
!
!                    \\\\\ MD EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  MDeditu
!-----------------------------------------------------------------------

  IF ( type == 'MDedu ' ) THEN

    IF ( name == 'iedMDu') THEN
      READ (line ,111) iedMDu
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDu,name
      CYCLE
    END IF ! name = 'iedMDu'

    IF ( name == 'nedMDu') THEN
      READ (line ,111) nedMDu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDu,name
      CYCLE
    END IF ! name = 'nedMDu'

    IF ( name == 'intedMDu') THEN
      READ (line ,111) intedMDu
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDu,name
      CYCLE
    END IF ! name = 'intedMDu'

    IF ( name == 'n_editMDu') THEN
      READ (line ,111) n_editMDu
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDu,name
      CYCLE
    END IF ! name = 'n_editMDu'

  END IF ! type = 'MDedu ' 

!-----------------------------------------------------------------------
!  MDeditv
!-----------------------------------------------------------------------

  IF ( type == 'MDedv ' ) THEN

    IF ( name == 'iedMDv') THEN
      READ (line ,111) iedMDv
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDv,name
      CYCLE
    END IF ! name = 'iedMDv'

    IF ( name == 'nedMDv') THEN
      READ (line ,111) nedMDv
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDv,name
      CYCLE
    END IF ! name = 'nedMDv'

    IF ( name == 'intedMDv') THEN
      READ (line ,111) intedMDv
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDv,name
      CYCLE
    END IF ! name = 'intedMDv'

    IF ( name == 'n_editMDv') THEN
      READ (line ,111) n_editMDv
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDv,name
      CYCLE
    END IF ! name = 'n_editMDv'

  END IF ! type = 'MDedv ' 

!-----------------------------------------------------------------------
!  MDeditw
!-----------------------------------------------------------------------

  IF ( type == 'MDedw ' ) THEN

    IF ( name == 'iedMDw') THEN
      READ (line ,111) iedMDw
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDw,name
      CYCLE
    END IF ! name = 'iedMDw'

    IF ( name == 'nedMDw') THEN
      READ (line ,111) nedMDw
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDw,name
      CYCLE
    END IF ! name = 'nedMDw'

    IF ( name == 'intedMDw') THEN
      READ (line ,111) intedMDw
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDw,name
      CYCLE
    END IF ! name = 'intedMDw'

    IF ( name == 'n_editMDw') THEN
      READ (line ,111) n_editMDw
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDw,name
      CYCLE
    END IF ! name = 'n_editMDw'

  END IF ! type = 'MDedw ' 

!-----------------------------------------------------------------------
!  MDedits
!-----------------------------------------------------------------------

  IF ( type == 'MDeds ' ) THEN

    IF ( name == 'iedMDs') THEN
      READ (line ,111) iedMDs
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDs,name
      CYCLE
    END IF ! name = 'iedMDs'

    IF ( name == 'nedMDs') THEN
      READ (line ,111) nedMDs
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDs,name
      CYCLE
    END IF ! name = 'nedMDs'

    IF ( name == 'intedMDs') THEN
      READ (line ,111) intedMDs
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDs,name
      CYCLE
    END IF ! name = 'intedMDs'

    IF ( name == 'n_editMDs') THEN
      READ (line ,111) n_editMDs
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDs,name
      CYCLE
    END IF ! name = 'n_editMDs'

  END IF ! type = 'MDeds ' 

!-----------------------------------------------------------------------
!  MDeditd
!-----------------------------------------------------------------------

  IF ( type == 'MDedd ' ) THEN

    IF ( name == 'iedMDd') THEN
      READ (line ,111) iedMDd
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDd,name
      CYCLE
    END IF ! name = 'iedMDd'

    IF ( name == 'nedMDd') THEN
      READ (line ,111) nedMDd
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDd,name
      CYCLE
    END IF ! name = 'nedMDd'

    IF ( name == 'intedMDd') THEN
      READ (line ,111) intedMDd
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDd,name
      CYCLE
    END IF ! name = 'intedMDd'

    IF ( name == 'n_editMDd') THEN
      READ (line ,111) n_editMDd
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDd,name
      CYCLE
    END IF ! name = 'n_editMDd'

  END IF ! type = 'MDedd ' 

!-----------------------------------------------------------------------
!  MDedite
!-----------------------------------------------------------------------

  IF ( type == 'MDede ' ) THEN

    IF ( name == 'iedMDe') THEN
      READ (line ,111) iedMDe
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDe,name
      CYCLE
    END IF ! name = 'iedMDe'

    IF ( name == 'nedMDe') THEN
      READ (line ,111) nedMDe
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDe,name
      CYCLE
    END IF ! name = 'nedMDe'

    IF ( name == 'intedMDe') THEN
      READ (line ,111) intedMDe
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDe,name
      CYCLE
    END IF ! name = 'intedMDe'

    IF ( name == 'n_editMDe') THEN
      READ (line ,111) n_editMDe
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDe,name
      CYCLE
    END IF ! name = 'n_editMDe'

  END IF ! type = 'MDede ' 

!-----------------------------------------------------------------------
!  MDeditp
!-----------------------------------------------------------------------

  IF ( type == 'MDedp ' ) THEN

    IF ( name == 'iedMDp') THEN
      READ (line ,111) iedMDp
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDp,name
      CYCLE
    END IF ! name = 'iedMDp'

    IF ( name == 'nedMDp') THEN
      READ (line ,111) nedMDp
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDp,name
      CYCLE
    END IF ! name = 'nedMDp'

    IF ( name == 'intedMDp') THEN
      READ (line ,111) intedMDp
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDp,name
      CYCLE
    END IF ! name = 'intedMDp'

    IF ( name == 'n_editMDp') THEN
      READ (line ,111) n_editMDp
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDp,name
      CYCLE
    END IF ! name = 'n_editMDp'

  END IF ! type = 'MDedp ' 

!-----------------------------------------------------------------------
!  MDeditenu
!-----------------------------------------------------------------------

  IF ( type == 'MDeden' ) THEN

    IF ( name == 'iedMDenu') THEN
      READ (line ,111) iedMDenu
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDenu,name
      CYCLE
    END IF ! name = 'iedMDenu'

    IF ( name == 'nedMDenu') THEN
      READ (line ,111) nedMDenu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDenu,name
      CYCLE
    END IF ! name = 'nedMDenu'

    IF ( name == 'intedMDenu') THEN
      READ (line ,111) intedMDenu
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDenu,name
      CYCLE
    END IF ! name = 'intedMDenu'

    IF ( name == 'n_editMDenu') THEN
      READ (line ,111) n_editMDenu
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDenu,name
      CYCLE
    END IF ! name = 'n_editMDenu'

  END IF ! type = 'MDedenu ' 

!-----------------------------------------------------------------------
!  MDeditfnu
!-----------------------------------------------------------------------

  IF ( type == 'MDedfn' ) THEN

    IF ( name == 'iedMDfnu') THEN
      READ (line ,111) iedMDfnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDfnu,name
      CYCLE
    END IF ! name = 'iedMDfnu'

    IF ( name == 'nedMDfnu') THEN
      READ (line ,111) nedMDfnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDfnu,name
      CYCLE
    END IF ! name = 'nedMDfnu'

    IF ( name == 'intedMDfnu') THEN
      READ (line ,111) intedMDfnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDfnu,name
      CYCLE
    END IF ! name = 'intedMDfnu'

    IF ( name == 'n_editMDfnu') THEN
      READ (line ,111) n_editMDfnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDfnu,name
      CYCLE
    END IF ! name = 'n_editMDfnu'

  END IF ! type = 'MDedfnu ' 

!-----------------------------------------------------------------------
!  MDedita
!-----------------------------------------------------------------------

  IF ( type == 'MDeda ' ) THEN

    IF ( name == 'iedMDa') THEN
      READ (line ,111) iedMDa
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDa,name
      CYCLE
    END IF ! name = 'iedMDa'

    IF ( name == 'nedMDa') THEN
      READ (line ,111) nedMDa
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDa,name
      CYCLE
    END IF ! name = 'nedMDa'

    IF ( name == 'intedMDa') THEN
      READ (line ,111) intedMDa
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDa,name
      CYCLE
    END IF ! name = 'intedMDa'

    IF ( name == 'n_editMDa') THEN
      READ (line ,111) n_editMDa
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDa,name
      CYCLE
    END IF ! name = 'n_editMDa'

  END IF ! type = 'MDeda ' 

!-----------------------------------------------------------------------
!  MDeditx
!-----------------------------------------------------------------------

  IF ( type == 'MDedx ' ) THEN

    IF ( name == 'iedMDx') THEN
      READ (line ,111) iedMDx
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDx,name
      CYCLE
    END IF ! name = 'iedMDx'

    IF ( name == 'nedMDx') THEN
      READ (line ,111) nedMDx
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDx,name
      CYCLE
    END IF ! name = 'nedMDx'

    IF ( name == 'intedMDx') THEN
      READ (line ,111) intedMDx
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDx,name
      CYCLE
    END IF ! name = 'intedMDx'

    IF ( name == 'n_editMDx') THEN
      READ (line ,111) n_editMDx
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDx,name
      CYCLE
    END IF ! name = 'n_editMDx'

  END IF ! type = 'MDedx ' 

!-----------------------------------------------------------------------
!  MDeditye
!-----------------------------------------------------------------------

  IF ( type == 'MDedye' ) THEN

    IF ( name == 'iedMDye') THEN
      READ (line ,111) iedMDye
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDye,name
      CYCLE
    END IF ! name = 'iedMDye'

    IF ( name == 'nedMDye') THEN
      READ (line ,111) nedMDye
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDye,name
      CYCLE
    END IF ! name = 'nedMDye'

    IF ( name == 'intedMDye') THEN
      READ (line ,111) intedMDye
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDye,name
      CYCLE
    END IF ! name = 'intedMDye'

    IF ( name == 'n_editMDye') THEN
      READ (line ,111) n_editMDye
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDye,name
      CYCLE
    END IF ! name = 'n_editMDye'

  END IF ! type = 'MDedye' 

!-----------------------------------------------------------------------
!  MDeditcm
!-----------------------------------------------------------------------

  IF ( type == 'MDedcm' ) THEN

    IF ( name == 'iedMDcm') THEN
      READ (line ,111) iedMDcm
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDcm,name
      CYCLE
    END IF ! name = 'iedMDcm'

    IF ( name == 'nedMDcm') THEN
      READ (line ,111) nedMDcm
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDcm,name
      CYCLE
    END IF ! name = 'nedMDcm'

    IF ( name == 'intedMDcm') THEN
      READ (line ,111) intedMDcm
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDcm,name
      CYCLE
    END IF ! name = 'intedMDcm'

    IF ( name == 'n_editMDcm') THEN
      READ (line ,111) n_editMDcm
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDcm,name
      CYCLE
    END IF ! name = 'n_editMDcm'

  END IF ! type = 'MDedcm ' 

!-----------------------------------------------------------------------
!  MDeditnu
!-----------------------------------------------------------------------

  IF ( type == 'MDednu' ) THEN

    IF ( name == 'iedMDnu') THEN
      READ (line ,111) iedMDnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDnu,name
      CYCLE
    END IF ! name = 'iedMDnu'

    IF ( name == 'nedMDnu') THEN
      READ (line ,111) nedMDnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDnu,name
      CYCLE
    END IF ! name = 'nedMDnu'

    IF ( name == 'intedMDnu') THEN
      READ (line ,111) intedMDnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDnu,name
      CYCLE
    END IF ! name = 'intedMDnu'

    IF ( name == 'n_editMDnu') THEN
      READ (line ,111) n_editMDnu
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDnu,name
      CYCLE
    END IF ! name = 'n_editMDnu'

  END IF ! type = 'MDednu ' 

!-----------------------------------------------------------------------
!  MDeditnc
!-----------------------------------------------------------------------

  IF ( type == 'MDednc' ) THEN

    IF ( name == 'iedMDnc') THEN
      READ (line ,111) iedMDnc
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDnc,name
      CYCLE
    END IF ! name = 'iedMDnc'

    IF ( name == 'nedMDnc') THEN
      READ (line ,111) nedMDnc
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDnc,name
      CYCLE
    END IF ! name = 'nedMDnc'

    IF ( name == 'intedMDnc') THEN
      READ (line ,111) intedMDnc
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDnc,name
      CYCLE
    END IF ! name = 'intedMDnc'

    IF ( name == 'n_editMDnc') THEN
      READ (line ,111) n_editMDnc
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDnc,name
      CYCLE
    END IF ! name = 'n_editMDnc'

  END IF ! type = 'MDednc ' 

!-----------------------------------------------------------------------
!  MDeditnl
!-----------------------------------------------------------------------

  IF ( type == 'MDednl' ) THEN

    IF ( name == 'iedMDnl') THEN
      READ (line ,111) iedMDnl
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDnl,name
      CYCLE
    END IF ! name = 'iedMDnl'

    IF ( name == 'nedMDnl') THEN
      READ (line ,111) nedMDnl
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDnl,name
      CYCLE
    END IF ! name = 'nedMDnl'

    IF ( name == 'intedMDnl') THEN
      READ (line ,111) intedMDnl
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDnl,name
      CYCLE
    END IF ! name = 'intedMDnl'

    IF ( name == 'n_editMDnl') THEN
      READ (line ,111) n_editMDnl
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDnl,name
      CYCLE
    END IF ! name = 'n_editMDnl'

  END IF ! type = 'MDednl' 

!-----------------------------------------------------------------------
!  MDeditne
!-----------------------------------------------------------------------

  IF ( type == 'MDedne' ) THEN

    IF ( name == 'iedMDne') THEN
      READ (line ,111) iedMDne
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDne,name
      CYCLE
    END IF ! name = 'iedMDne'

    IF ( name == 'nedMDne') THEN
      READ (line ,111) nedMDne
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDne,name
      CYCLE
    END IF ! name = 'nedMDne'

    IF ( name == 'intedMDne') THEN
      READ (line ,111) intedMDne
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDne,name
      CYCLE
    END IF ! name = 'intedMDne'

    IF ( name == 'n_editMDne') THEN
      READ (line ,111) n_editMDne
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDne,name
      CYCLE
    END IF ! name = 'n_editMDne'

  END IF ! type = 'MDedne' 

!-----------------------------------------------------------------------
!  MDeditgx
!-----------------------------------------------------------------------

  IF ( type == 'MDedgx' ) THEN

    IF ( name == 'iedMDgx') THEN
      READ (line ,111) iedMDgx
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDgx,name
      CYCLE
    END IF ! name = 'iedMDgx'

    IF ( name == 'nedMDgx') THEN
      READ (line ,111) nedMDgx
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDgx,name
      CYCLE
    END IF ! name = 'nedMDgx'

    IF ( name == 'intedMDgx') THEN
      READ (line ,111) intedMDgx
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDgx,name
      CYCLE
    END IF ! name = 'intedMDgx'

    IF ( name == 'n_editMDgx') THEN
      READ (line ,111) n_editMDgx
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDgx,name
      CYCLE
    END IF ! name = 'n_editMDgx'

  END IF ! type = 'MDedgx' 

!-----------------------------------------------------------------------
!  MDeditgy
!-----------------------------------------------------------------------

  IF ( type == 'MDedgy' ) THEN

    IF ( name == 'iedMDgy') THEN
      READ (line ,111) iedMDgy
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDgy,name
      CYCLE
    END IF ! name = 'iedMDgy'

    IF ( name == 'nedMDgy') THEN
      READ (line ,111) nedMDgy
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDgy,name
      CYCLE
    END IF ! name = 'nedMDgy'

    IF ( name == 'intedMDgy') THEN
      READ (line ,111) intedMDgy
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDgy,name
      CYCLE
    END IF ! name = 'intedMDgy'

    IF ( name == 'n_editMDgy') THEN
      READ (line ,111) n_editMDgy
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDgy,name
      CYCLE
    END IF ! name = 'n_editMDgy'

  END IF ! type = 'MDedgy' 

!-----------------------------------------------------------------------
!  MDeditBVw
!-----------------------------------------------------------------------

  IF ( type == 'MDedBV' ) THEN

    IF ( name == 'iedMDBVw') THEN
      READ (line ,111) iedMDBVw
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDBVw,name
      CYCLE
    END IF ! name = 'iedMDBVw'

    IF ( name == 'nedMDBVw') THEN
      READ (line ,111) nedMDBVw
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDBVw,name
      CYCLE
    END IF ! name = 'nedMDBVw'

    IF ( name == 'intedMDBVw') THEN
      READ (line ,111) intedMDBVw
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDBVw,name
      CYCLE
    END IF ! name = 'intedMDBVw'

    IF ( name == 'n_editMDBVw') THEN
      READ (line ,111) n_editMDBVw
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDBVw,name
      CYCLE
    END IF ! name = 'n_editMDBVw'

  END IF ! type = 'MDedBV' 

!-----------------------------------------------------------------------
!  MDeditgy
!-----------------------------------------------------------------------

  IF ( type == 'MDedyl' ) THEN

    IF ( name == 'iedMDyl') THEN
      READ (line ,111) iedMDyl
      IF ( iskip == 0 ) WRITE (nprint,113) type,iedMDyl,name
      CYCLE
    END IF ! name = 'iedMDyl'

    IF ( name == 'nedMDyl') THEN
      READ (line ,111) nedMDyl
      IF ( iskip == 0 ) WRITE (nprint,113) type,nedMDyl,name
      CYCLE
    END IF ! name = 'nedMDyl'

    IF ( name == 'intedMDyl') THEN
      READ (line ,111) intedMDyl
      IF ( iskip == 0 ) WRITE (nprint,113) type,intedMDyl,name
      CYCLE
    END IF ! name = 'intedMDyl'

    IF ( name == 'n_editMDyl') THEN
      READ (line ,111) n_editMDyl
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMDyl,name
      CYCLE
    END IF ! name = 'n_editMDyl'

  END IF ! type = 'MDedyl' 

!-----------------------------------------------------------------------
!  Time parameters for MD edits
!-----------------------------------------------------------------------

  IF ( type == 'MDedit' ) THEN

    IF ( name == 'i_editMD') THEN
      READ (line ,111) i_editMD
      IF ( iskip == 0 ) WRITE (nprint,113) type,i_editMD,name
      CYCLE
    END IF ! name = 'i_editMD'

    IF ( name == 'n_editMD') THEN
      READ (line ,111) n_editMD
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_editMD,name
      CYCLE
    END IF ! name = 'n_editMD'

    IF ( name == 'dt_MDedit1') THEN
      READ (line ,141) dt_MDedit1
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_MDedit1,name
      CYCLE
    END IF ! name = 'dt_MDedit1'

    IF ( name == 'dt_MDedit2') THEN
      READ (line ,141) dt_MDedit2
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_MDedit2,name
      CYCLE
    END IF ! name = 'dt_MDedit2'

  END IF ! type = 'MDedit ' 

!-----------------------------------------------------------------------
!
!                 \\\\\ GLOBAL EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Cycle parameters for global edits
!-----------------------------------------------------------------------

  IF ( type == 'ed_gln' ) THEN

    IF ( name == 'ied_global_n') THEN
      READ (line ,111) ied_global_n
      IF ( iskip == 0 ) WRITE (nprint,113) type,ied_global_n,name
      CYCLE
    END IF ! name = 'ied_global_n'

    IF ( name == 'nned_global') THEN
      READ (line ,111) nned_global
      IF ( iskip == 0 ) WRITE (nprint,113) type,nned_global,name
      CYCLE
    END IF ! name = 'nned_global'

    IF ( name == 'inted_global') THEN
      READ (line ,111) inted_global
      IF ( iskip == 0 ) WRITE (nprint,113) type,inted_global,name
      CYCLE
    END IF ! name = 'inted_global'

    IF ( name == 'n_edit_global') THEN
      READ (line ,111) n_edit_global
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_edit_global,name
      CYCLE
    END IF ! name = 'n_edit_global'

  END IF ! type = 'ed_gln' 

!-----------------------------------------------------------------------
!  Time parameters for global edits
!-----------------------------------------------------------------------

  IF ( type == 'ed_glt' ) THEN

    IF ( name == 'ied_global_t') THEN
      READ (line ,111) ied_global_t
      IF ( iskip == 0 ) WRITE (nprint,113) type,ied_global_t,name
      CYCLE
    END IF ! name = 'ied_global_t'

    IF ( name == 'nted_global') THEN
      READ (line ,111) nted_global
      IF ( iskip == 0 ) WRITE (nprint,113) type,nted_global,name
      CYCLE
    END IF ! name = 'nted_global'

    IF ( name == 'dt_global_ed1') THEN
      READ (line ,141) dt_global_ed1
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_global_ed1,name
      CYCLE
    END IF ! name = 'dt_global_ed1'

    IF ( name == 'dt_global_ed2') THEN
      READ (line ,141) dt_global_ed2
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_global_ed2,name
      CYCLE
    END IF ! name = 'dt_global_ed2'

  END IF ! type = 'ed_glt ' 

!-----------------------------------------------------------------------
!
!                    \\\\\ HDF EDIT PARAMETERS /////
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  HDFedit
!-----------------------------------------------------------------------

  IF ( type == 'HDFedt' ) THEN

    IF ( name == 'i_HDFedit') THEN
      READ (line ,111) i_HDFedit
      IF ( iskip == 0 ) WRITE (nprint,113) type,i_HDFedit,name
      CYCLE
    END IF ! name = 'i_HDFedit'

    IF ( name == 'nd_HDFedit') THEN
      READ (line ,111) nd_HDFedit
      IF ( iskip == 0 ) WRITE (nprint,113) type,nd_HDFedit,name
      CYCLE
    END IF ! name = 'nd_HDFedit'

    IF ( name == 'n_HDFedit') THEN
      READ (line ,111) n_HDFedit
      IF ( iskip == 0 ) WRITE (nprint,113) type,n_HDFedit,name
      CYCLE
    END IF ! name = 'n_HDFedit'

    IF ( name == 'dt_HDFedit1') THEN
      READ (line ,141) dt_HDFedit1
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_HDFedit1,name
      CYCLE
    END IF ! name = 'dt_HDFedit1'

    IF ( name == 'dt_HDFedit2') THEN
      READ (line ,141) dt_HDFedit2
      IF ( iskip == 0 ) WRITE (nprint,143) type,dt_HDFedit2,name
      CYCLE
    END IF ! name = 'dt_HDFedit2'

  END IF ! type = 'e_plt ' 

!-----------------------------------------------------------------------
!  Unrecognized card
!-----------------------------------------------------------------------

  IF ( nrst == 0 ) THEN
    WRITE (nprint,401)
    WRITE (nprint,403) line
  END IF

END DO READ

!-----------------------------------------------------------------------
!
!                   \\\\\ PACK EDIT KEYS /////
!
!-----------------------------------------------------------------------

 5000 CONTINUE

i_edit_data(1)                    = iflprt
i_edit_data(2)                    = noutfl
i_edit_data(3)                    = iprint
i_edit_data(4)                    = nrstd1
i_edit_data(5)                    = nrstd2
i_edit_data(6)                    = noutpmt
i_edit_data(7)                    = nouttmp
i_edit_data(8)                    = iplot
i_edit_data(9)                    = nplotc
i_edit_data(10)                   = nplote
i_edit_data(11)                   = nplota
i_edit_data(12)                   = nplott
i_edit_data(13)                   = intplf
i_edit_data(14)                   = npltf
i_edit_data(15)                   = nmodel
i_edit_data(16)                   = intprt
i_edit_data(17)                   = nprt
i_edit_data(18)                   = iprtbgn
i_edit_data(19)                   = intrst
i_edit_data(20)                   = nnrst
i_edit_data(21)                   = nrstfl
i_edit_data(22)                   = irstbgn
i_edit_data(23)                   = ncychg
i_edit_data(24)                   = intprm
i_edit_data(25)                   = nprm
i_edit_data(26)                   = ilumplt
i_edit_data(27)                   = nlumplt1
i_edit_data(28)                   = nlumplt2
i_edit_data(29)                   = nlum
i_edit_data(30)                   = intlum(1)
i_edit_data(31)                   = intlum(2)
i_edit_data(32)                   = intlum(3)
i_edit_data(33)                   = ncylum(1)
i_edit_data(34)                   = ncylum(2)
i_edit_data(35)                   = ienuplt
i_edit_data(36)                   = nenuplt1
i_edit_data(37)                   = nenuplt2
i_edit_data(38)                   = nenu
i_edit_data(39)                   = intenu(1)
i_edit_data(40)                   = intenu(2)
i_edit_data(41)                   = intenu(3)
i_edit_data(42)                   = ncyenu(1)
i_edit_data(43)                   = ncyenu(2)
i_edit_data(44)                   = irnuplt
i_edit_data(45)                   = nrnuplt
i_edit_data(46)                   = nrnu
i_edit_data(47)                   = intrnu(1)
i_edit_data(48)                   = intrnu(2)
i_edit_data(49)                   = intrnu(3)
i_edit_data(50)                   = ncyrnu(1)
i_edit_data(51)                   = ncyrnu(1)
i_edit_data(52)                   = ivarplt
i_edit_data(53)                   = nvar
i_edit_data(54)                   = nvarint
i_edit_data(55)                   = nvarplt
i_edit_data(56)                   = nvardump
i_edit_data(57)                   = icomplt
i_edit_data(58)                   = ncomplt
i_edit_data(59)                   = ncomdump
i_edit_data(60)                   = nplotinnerb
i_edit_data(61)                   = iplotinnerb
i_edit_data(62)                   = nplotouterb
i_edit_data(63)                   = iplotouterb
i_edit_data(64)                   = nplotlum
i_edit_data(65)                   = iplotlum
i_edit_data(66)                   = nplotshk
i_edit_data(67)                   = iplotshk
i_edit_data(68)                   = nplotcnv
i_edit_data(69)                   = iplotcnv
i_edit_data(70)                   = nplotnurad
i_edit_data(71)                   = iplotnurad
i_edit_data(72)                   = nnudata
i_edit_data(73)                   = inudata
i_edit_data(74)                   = nlagplt
i_edit_data(75)                   = ilagplt
i_edit_data(76)                   = nlagdump
i_edit_data(77)                   = nrlagplt
i_edit_data(78)                   = irlagplt
i_edit_data(79)                   = n_eplt
i_edit_data(80)                   = i_eplt
i_edit_data(81)                   = iedMDu
i_edit_data(82)                   = nedMDu
i_edit_data(83)                   = intedMDu
i_edit_data(84)                   = n_editMDu
i_edit_data(85)                   = iedMDv
i_edit_data(86)                   = nedMDv
i_edit_data(87)                   = intedMDv
i_edit_data(88)                   = n_editMDv
i_edit_data(89)                   = iedMDw
i_edit_data(90)                   = nedMDw
i_edit_data(91)                   = intedMDw
i_edit_data(92)                   = n_editMDw
i_edit_data(93)                   = iedMDs
i_edit_data(94)                   = nedMDs
i_edit_data(95)                   = intedMDs
i_edit_data(96)                   = n_editMDs
i_edit_data(97)                   = iedMDd
i_edit_data(98)                   = nedMDd
i_edit_data(99)                   = intedMDd
i_edit_data(100)                  = n_editMDd
i_edit_data(101)                  = iedMDe
i_edit_data(102)                  = nedMDe
i_edit_data(103)                  = intedMDe
i_edit_data(104)                  = n_editMDe
i_edit_data(105)                  = iedMDp
i_edit_data(106)                  = nedMDp
i_edit_data(107)                  = intedMDp
i_edit_data(108)                  = n_editMDp
i_edit_data(109)                  = iedMDenu
i_edit_data(110)                  = nedMDenu
i_edit_data(111)                  = intedMDenu
i_edit_data(112)                  = n_editMDenu
i_edit_data(113)                  = iedMDfnu
i_edit_data(114)                  = nedMDfnu
i_edit_data(115)                  = intedMDfnu
i_edit_data(116)                  = n_editMDfnu
i_edit_data(117)                  = iedMDa
i_edit_data(118)                  = nedMDa
i_edit_data(119)                  = intedMDa
i_edit_data(120)                  = n_editMDa
i_edit_data(121)                  = iedMDx
i_edit_data(122)                  = nedMDx
i_edit_data(123)                  = intedMDx
i_edit_data(124)                  = n_editMDx
i_edit_data(125)                  = iedMDye
i_edit_data(126)                  = nedMDye
i_edit_data(127)                  = intedMDye
i_edit_data(128)                  = n_editMDye
i_edit_data(129)                  = iedMDcm
i_edit_data(130)                  = nedMDcm
i_edit_data(131)                  = intedMDcm
i_edit_data(132)                  = n_editMDcm
i_edit_data(133)                  = iedMDnu
i_edit_data(134)                  = nedMDnu
i_edit_data(135)                  = intedMDnu
i_edit_data(136)                  = n_editMDnu
i_edit_data(137)                  = iedMDnc
i_edit_data(138)                  = nedMDnc
i_edit_data(139)                  = intedMDnc
i_edit_data(140)                  = n_editMDnc
i_edit_data(141)                  = iedMDnl
i_edit_data(142)                  = nedMDnl
i_edit_data(143)                  = intedMDnl
i_edit_data(144)                  = n_editMDnl
i_edit_data(145)                  = iedMDne
i_edit_data(146)                  = nedMDne
i_edit_data(147)                  = intedMDne
i_edit_data(148)                  = n_editMDne
i_edit_data(149)                  = iedMDgx
i_edit_data(150)                  = nedMDgx
i_edit_data(151)                  = intedMDgx
i_edit_data(152)                  = n_editMDgx
i_edit_data(153)                  = iedMDgy
i_edit_data(154)                  = nedMDgy
i_edit_data(155)                  = intedMDgy
i_edit_data(156)                  = n_editMDgy
i_edit_data(157)                  = iedMDBVw
i_edit_data(158)                  = nedMDBVw
i_edit_data(159)                  = intedMDBVw
i_edit_data(160)                  = n_editMDBVw
i_edit_data(161)                  = iedMDyl
i_edit_data(162)                  = nedMDyl
i_edit_data(163)                  = intedMDyl
i_edit_data(164)                  = n_editMDyl
i_edit_data(165)                  = i_editMD
i_edit_data(166)                  = n_editMD
i_edit_data(167)                  = i_HDFedit
i_edit_data(168)                  = nd_HDFedit
i_edit_data(169)                  = n_HDFedit
i_edit_data(170)                  = it_edit
i_edit_data(171)                  = ied_global_n
i_edit_data(172)                  = nned_global
i_edit_data(173)                  = inted_global
i_edit_data(174)                  = n_edit_global
i_edit_data(175)                  = ied_global_t
i_edit_data(176)                  = nted_global

DO i = 1,20
  i_edit_data(200+i)              = intedc(i)
  i_edit_data(220+i)              = nedc(i)
  i_edit_data(240+i)              = idxedc(i)
  i_edit_data(260+i)              = intede(i)
  i_edit_data(280+i)              = nede(i)
  i_edit_data(300+i)              = idxede(i)
  i_edit_data(320+i)              = intdmi(i)
  i_edit_data(340+i)              = nedmi(i)
  i_edit_data(360+i)              = idxemi(i)
  i_edit_data(380+i)              = intdma(i)
  i_edit_data(400+i)              = nedma(i)
  i_edit_data(420+i)              = idxema(i)
  i_edit_data(440+i)              = intedh(i)
  i_edit_data(460+i)              = nedh(i)
  i_edit_data(480+i)              = idxedh(i)
  i_edit_data(500+i)              = intdps(i)
  i_edit_data(520+i)              = nedps(i)
  i_edit_data(540+i)              = idxeps(i)
  i_edit_data(560+i)              = intedu(i)
  i_edit_data(580+i)              = nedu(i)
  i_edit_data(600+i)              = idxedu(i)
  i_edit_data(620+i)              = intedy(i)
  i_edit_data(640+i)              = nedy(i)
  i_edit_data(660+i)              = idxedy(i)
  i_edit_data(680+i)              = intdsc(i)
  i_edit_data(700+i)              = nedsc(i)
  i_edit_data(720+i)              = idxesc(i)
END DO

DO i = 1,60
  i_edit_data(800+i)              = niedn(i)
  i_edit_data(900+i)              = neden(i)
END DO

DO i = 1,100
  i_edit_data(1000+i)             = ncyrst(i)
END DO

DO n = 1,nnu
  i_edit_data(1100+0*nnu+n)       = intedn(n)
  i_edit_data(1100+1*nnu+n)       = nedn(n)
  i_edit_data(1100+2*nnu+n)       = idxedn(n)
END DO

DO n = 1,nnu
  DO k = 1,40
    i_edit_data(1200+0*40+(n-1)*3*40+k) = intdng(k,n)
    i_edit_data(1200+1*40+(n-1)*3*40+k) = nedng(k,n)
    i_edit_data(1200+2*40+(n-1)*3*40+k) = idxeng(k,n)
  END DO
END DO

d_edit_data(1)                    = r_lum
d_edit_data(2)                    = d_lum
d_edit_data(3)                    = r_e_rms
d_edit_data(4)                    = d_e_rms
d_edit_data(5)                    = dtvarplot
d_edit_data(6)                    = dtcomplot
d_edit_data(7)                    = dtimeplot
d_edit_data(8)                    = rinnerb
d_edit_data(9)                    = routerb
d_edit_data(10)                   = r_lumerms
d_edit_data(11)                   = dtnuradplot
d_edit_data(12)                   = r_nurad
d_edit_data(13)                   = rho_nurad
d_edit_data(14)                   = dtnudata
d_edit_data(15)                   = r_nudata
d_edit_data(16)                   = t_nudata
d_edit_data(17)                   = msslag
d_edit_data(18)                   = dmlag
d_edit_data(19)                   = r_pinch
d_edit_data(20)                   = d_pinch
d_edit_data(21)                   = dt_eplot
d_edit_data(22)                   = dt_MDedit1
d_edit_data(23)                   = dt_MDedit2
d_edit_data(24)                   = dt_HDFedit1
d_edit_data(25)                   = dt_HDFedit2
d_edit_data(26)                   = dt_edit
d_edit_data(27)                   = dt_global_ed1
d_edit_data(28)                   = dt_global_ed2

DO i = 1,10
  d_edit_data(40+i)               = rhoprint(i)
END DO

!-----------------------------------------------------------------------
!  Deallocate arrays
!-----------------------------------------------------------------------

DEALLOCATE (intedn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intedn    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nedn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedn      '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (idxedn, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxedn    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (intdng, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'intdng    '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (nedng, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'nedng     '; WRITE (nlog,2001) var_name; END IF
DEALLOCATE (idxeng, STAT = istat)
  IF ( istat /= 0 ) THEN; var_name = 'idxeng    '; WRITE (nlog,2001) var_name; END IF

RETURN
END SUBROUTINE read_pack_edit_keys
