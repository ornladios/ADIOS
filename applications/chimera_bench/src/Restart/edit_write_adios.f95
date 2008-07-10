SUBROUTINE edit_write_adios( ndump, nez, nnu )
!-----------------------------------------------------------------------
!
!    File:         edit_write
!    Module:       edit_write
!    Type:         Subprogram
!    Author:       S. W. Bruenn, Dept of Physics, FAU,
!                  Boca Raton, FL 33431-0991
!
!    Date:         1/22/05
!
!    Purpose:
!        To dump the edit keys and current values of the edit counters.
!
!    Subprograms called:
!        none
!
!    Input arguments:
!  ndump      : unit number to print edit parameter dump
!
!    Output arguments:
!        none
!
!    Include files:
!      cycle_module, edit_module
!
!-----------------------------------------------------------------------

USE cycle_module, ONLY : ncycle
USE edit_module

IMPLICIT none
SAVE

!-----------------------------------------------------------------------
!        Input variables.
!-----------------------------------------------------------------------

INTEGER*8, INTENT(in)              :: ndump           ! unit number to print edit parameter dump
INTEGER, INTENT(in)              :: nez             ! neutrino energy array extent
INTEGER, INTENT(in)              :: nnu             ! neutrino flavor array extent

!-----------------------------------------------------------------------
!        Local variables
!-----------------------------------------------------------------------

INTEGER                          :: i               ! do index
INTEGER                          :: n               ! neutrino flavor index
INTEGER, PARAMETER               :: ij_ray = 1
INTEGER, PARAMETER               :: ik_ray = 1
  
    1 FORMAT ('!                  \\\\\                ///// ')
    2 FORMAT ('!                          EDIT KEYS')
    3 FORMAT ('!                  /////                \\\\\ ')

   13 FORMAT (/'!-----------------------------------------------------------------------')
   15 FORMAT ('!-----------------------------------------------------------------------'/)

!........data_path

   21 FORMAT ('!      data_path: the path directing output to the data directories.')
   22 FORMAT ('data_p',4x,a118)

!........iflprt

   26 FORMAT ('!      iflprt: a print parameter regulating the print file created when')
   27 FORMAT ('!       the calculation is terminated (ncycle = ncymax)')
   28 FORMAT ('iflprt',14x,i10,42x,'iflprt')

!........intprt

   31 FORMAT ('!      intprt: the number of cycles between the closing of')
   32 FORMAT ('!       an old opening of a new print file.')
   33 FORMAT ('!      nprt: number of cycles since the last closing of a print file.')
   34 FORMAT ('!      iprtbgn: initial print flag.')
   35 FORMAT ('intprt',14x,i10,42x,'intprt')
   36 FORMAT ('intprt',14x,i10,42x,'nprt')
   37 FORMAT ('intprt',14x,i10,42x,'iprtbgn')

!........output

   41 FORMAT ('!      iprint : print flag (1: print; 0: no print')
   42 FORMAT ('!      nrstd1 : unit number to dump temporary restart file')
   43 FORMAT ('!      nrstd2 : unit number to dump temporary restart file')
   44 FORMAT ('!      noutpmt : unit number to dump permanent restart files')
   45 FORMAT ('!      nouttmp : unit number to dump next temporary restart files')
   46 FORMAT ('output',14x,i10,42x,'iprint')
   47 FORMAT ('output',14x,i10,42x,'nrstd1')
   48 FORMAT ('output',14x,i10,42x,'nrstd2')
   49 FORMAT ('output',14x,i10,42x,'noutpmt')
   50 FORMAT ('output',14x,i10,42x,'nouttmp')

!........t_edit

   51 FORMAT ('!      it_edit : print flag for editing at specified time intervals after')
   52 FORMAT ('!      dt_edit : time interval (in ms) between time interval edits.')
   53 FORMAT ('t_edit',14x,i10,42x,'it_edit')
   54 FORMAT ('t_edit',29x,1pe15.8,22x,'dt_edit')

!........plot

   61 FORMAT ('!      iplot: plot file flag')
   62 FORMAT ('!      nplotc: unit number for downloading model configuration plot data.')
   63 FORMAT ('!      nplote: unit number for downloading electron neutrino plot data.')
   64 FORMAT ('!      nplota: unit number for downloading electron antineutrino plot data.')
   65 FORMAT ('!      nplott: unit number for downloading muon and tau neutrino plot data.')
   66 FORMAT ('!      intplf: the number of cycles between the printing to a plot file.')
   67 FORMAT ('!      npltf: the number of cycles since the last printing to a plot file.')
   68 FORMAT ('plot  ',14x,i10,42x,'iplot')
   69 FORMAT ('plot  ',14x,i10,42x,'nplotc')
   70 FORMAT ('plot  ',14x,i10,42x,'nplote')
   71 FORMAT ('plot  ',14x,i10,42x,'nplota')
   72 FORMAT ('plot  ',14x,i10,42x,'nplott')
   73 FORMAT ('plot  ',14x,i10,42x,'intplf')
   74 FORMAT ('plot  ',14x,i10,42x,'npltf')

!........model number

   81 FORMAT ('!      nmodel: the number of the current print file.')
   82 FORMAT ('nmodel',14x,i10,42x,'nmodel')

!........rhoprt

   91 FORMAT ('!      rhoprt(i) : the ith central density at which a print file is created.')
   92 FORMAT ('rhoprt',4x,i10,5x,1pe15.8)

!........intedc

  101 FORMAT ('!      intedc(i): number of cycles between subsection i of subroutine editc.')
  102 FORMAT ('!      nedc(i): number of cycles since the last subsection i of subroutine editc.')
  103 FORMAT ('!      idxedc(i): print every idxedc(i) lines of subsection i of subroutine editc.')
  104 FORMAT ('intedc',4x,2i10,42x,'intedc')
  105 FORMAT ('intedc',4x,2i10,42x,'nedc')
  106 FORMAT ('intedc',4x,2i10,42x,'idxedc')

!........intede

  111 FORMAT ('!      intede(i): number of cycles between subsection i of subroutine edit_e.')
  112 FORMAT ('!      nede(i): number of cycles since the last subsection i of subroutine edit_e.')
  113 FORMAT ('!      idxede(i): print every idxedc(i) lines of subsection i of subroutine edit_e.')
  114 FORMAT ('intede',4x,2i10,42x,'intede')
  115 FORMAT ('intede',4x,2i10,42x,'nede')
  116 FORMAT ('intede',4x,2i10,42x,'idxede')

!........intdmi

  121 FORMAT ('!      intdmi(i): number of cycles between subsection i of subroutine editmi.')
  122 FORMAT ('!      nedmi(i): number of cycles since the last subsection i of subroutine editmi.')
  123 FORMAT ('!      idxemi(i): print every idxedc(i) lines of subsection i of subroutine editmi.')
  124 FORMAT ('intdmi',4x,2i10,42x,'intdmi')
  125 FORMAT ('intdmi',4x,2i10,42x,'nedmi')
  126 FORMAT ('intdmi',4x,2i10,42x,'idxemi')

!........intdma

  131 FORMAT ('!      intdma(i): number of cycles between subsection i of subroutine editma.')
  132 FORMAT ('!      nedma(i): number of cycles since the last subsection i of subroutine editma.')
  133 FORMAT ('!      idxema(i): print every idxedc(i) lines of subsection i of subroutine editma.')
  134 FORMAT ('intdma',4x,2i10,42x,'intdma')
  135 FORMAT ('intdma',4x,2i10,42x,'nedma')
  136 FORMAT ('intdma',4x,2i10,42x,'idxema')

!........intedh

  141 FORMAT ('!      intedh(i): number of cycles between subsection i of subroutine edith.')
  142 FORMAT ('!      nedh(i): number of cycles since the last subsection i of subroutine edith.')
  143 FORMAT ('!      idxedh(i): print every idxedc(i) lines of subsection i of subroutine edith.')
  144 FORMAT ('intedh',4x,2i10,42x,'intedh')
  145 FORMAT ('intedh',4x,2i10,42x,'nedh')
  146 FORMAT ('intedh',4x,2i10,42x,'idxedh')

!........intdps

  151 FORMAT ('!      intdps(i): number of cycles between subsection i of subroutine editps.')
  152 FORMAT ('!      nedps(i): number of cycles since the last subsection i of subroutine editps.')
  153 FORMAT ('!      idxeps(i): print every idxedc(i) lines of subsection i of subroutine editps.')
  154 FORMAT ('intdps',4x,2i10,42x,'intdps')
  155 FORMAT ('intdps',4x,2i10,42x,'nedps')
  156 FORMAT ('intdps',4x,2i10,42x,'idxeps')

!........intedu

  161 FORMAT ('!      intedu(i): number of cycles between subsection i of subroutine editu.')
  162 FORMAT ('!      nedu(i): number of cycles since the last subsection i of subroutine editu.')
  163 FORMAT ('!      idxedu(i): print every idxedc(i) lines of subsection i of subroutine editu.')
  164 FORMAT ('intedu',4x,2i10,42x,'intedu')
  165 FORMAT ('intedu',4x,2i10,42x,'nedu')
  166 FORMAT ('intedu',4x,2i10,42x,'idxedu')

!........intedy

  171 FORMAT ('!      intedy(i): number of cycles between subsection i of subroutine edity.')
  172 FORMAT ('!      nedy(i): number of cycles since the last subsection i of subroutine edity.')
  173 FORMAT ('!      idxedy(i): print every idxedc(i) lines of subsection i of subroutine edity.')
  174 FORMAT ('intedy',4x,2i10,42x,'intedy')
  175 FORMAT ('intedy',4x,2i10,42x,'nedy')
  176 FORMAT ('intedy',4x,2i10,42x,'idxedy')

!........intdsc

  181 FORMAT ('!      intdsc(i): number of cycles between subsection i of subroutine editsc.')
  182 FORMAT ('!      nedsc(i): number of cycles since the last subsection i of subroutine editsc.')
  183 FORMAT ('!      idxesc(i): print every idxedc(i) lines of subsection i of subroutine editsc.')
  184 FORMAT ('intdsc',4x,2i10,42x,'intdsc')
  185 FORMAT ('intdsc',4x,2i10,42x,'nedsc')
  186 FORMAT ('intdsc',4x,2i10,42x,'idxesc')

!........intedn

  191 FORMAT ('!      intedn(i): number of cycles between subsection i of subroutine editn.')
  192 FORMAT ('!      nedn(i): number of cycles since the last subsection i of subroutine editn.')
  193 FORMAT ('!      idxedn(i): print every idxedc(i) lines of subsection i of subroutine editn.')
  194 FORMAT ('!      niedn(i): the ith datum is printed every niedn(i) implementation of editn.')
  195 FORMAT ('!      neden(i): number of implementations of editn since the last printing of the ith datum.')
  196 FORMAT ('intedn',4x,2i10,42x,'intedn')
  197 FORMAT ('intedn',4x,2i10,42x,'nedn')
  198 FORMAT ('intedn',4x,2i10,42x,'idxedn')
  199 FORMAT ('intedn',4x,2i10,42x,'niedn')
  200 FORMAT ('intedn',4x,2i10,42x,'neden')

!........intdng

  211 FORMAT ('!      intdng(i): number of cycles between subsection i of subroutine editng.')
  212 FORMAT ('!      nedng(i): number of cycles since the last subsection i of subroutine editng.')
  213 FORMAT ('!      idxeng(i): print every idxedc(i) lines of subsection i of subroutine editng.')
  214 FORMAT ('intdng',4x,3i10,32x,'intdng')
  215 FORMAT ('intdng',4x,3i10,32x,'nedng')
  216 FORMAT ('intdng',4x,3i10,32x,'idxeng')

!........t_edit

  221 FORMAT ('!      Parameters for edits at selected postbounce time intervals.')
  222 FORMAT ('!      it_edit : key for editing at selected postbounce time intervals.')
  223 FORMAT ('!      dt_edit : postbounce time intervals for editing (ms).')
  224 FORMAT ('t_edit',14x,i10,5x,1pe15.8,22x,'t_edit')

!cccccc   Plot arrays

  231 FORMAT (/'cccccc   Plot arrays'/)

!........lumplt

  232 FORMAT ('!      ilumplt : luminosity plot edit switch.')
  233 FORMAT ('!      nlumplt1 : the unit number for editing lumplot data at r_lum.')
  234 FORMAT ('!      nlumplt2 : the unit number for editing lumplot data at d_lum.')
  235 FORMAT ('!      nlum : number of cycles since the last implemtation of subroutine lumplot.')
  236 FORMAT ('!      r_lum : the radius at which to evaluate luminosity data.')
  237 FORMAT ('!      d_lum : the density at which to evaluate luminosity data.')
  238 FORMAT ('!      intlum : cycles between luminosity edits.')
  239 FORMAT ('!      ncylum : cycle number at which to change intlum.')
  240 FORMAT ('lumplt',14x,i10,42x,'ilumplt')
  241 FORMAT ('lumplt',4x,2i10,42x,'nlumplt')
  242 FORMAT ('lumplt',29x,1pe15.8,22x,'r_lum')
  243 FORMAT ('lumplt',29x,1pe15.8,22x,'d_lum')
  244 FORMAT ('lumplt',14x,i10,42x,'nlum')
  245 FORMAT ('lumplt',4x,2i10,42x,'intlum')
  246 FORMAT ('lumplt',4x,2i10,42x,'ncylum')

!........enuplt

  251 FORMAT ('!      ienuplt : rms neutrino energy plot edit switch.')
  252 FORMAT ('!      nenuplt1 : the unit number for editing enuvplot data at r_e_rms.')
  253 FORMAT ('!      nenuplt2 : the unit number for editing enuvplot data at d_e_rms.')
  254 FORMAT ('!      nenu: the number of cycles since the last implemtation of subroutine enuvplot.')
  255 FORMAT ('!      r_e_rms: the radius at which to store neutrino energy data.')
  256 FORMAT ('!      d_e_rms: the density at which to store neutrino energy data.')
  257 FORMAT ('!      intenu : cycles between rms energy edits.')
  258 FORMAT ('!      ncyenu : cycle number at which to change intenu.')
  259 FORMAT ('enuplt',14x,i10,42x,'ienuplt')
  260 FORMAT ('enuplt',4x,2i10,42x,'nenuplt')
  261 FORMAT ('enuplt',29x,1pe15.8,22x,'r_e_rms')
  262 FORMAT ('enuplt',29x,1pe15.8,22x,'d_e_rms')
  263 FORMAT ('enuplt',14x,i10,42x,'nenu')
  264 FORMAT ('enuplt',4x,2i10,42x,'intenu')
  265 FORMAT ('enuplt',4x,2i10,42x,'ncyenu')

!........rnuplt

  271 FORMAT ('!      irnuplt : selected radii plot edit switch.')
  272 FORMAT ('!      nrnuplt : tthe unit number for editing rnuplot data.')
  273 FORMAT ('!      nrnu: the number of cycles since the last implemtation of subroutine rnuplot.')
  274 FORMAT ('!      intrnu : cycles between radius edits.')
  275 FORMAT ('!      ncyrnu : cycle number at which to change intrnu.')
  276 FORMAT ('rnuplt',14x,i10,42x,'irnuplt')
  277 FORMAT ('rnuplt',14x,i10,42x,'nrnuplt')
  278 FORMAT ('rnuplt',14x,i10,42x,'nrnu')
  279 FORMAT ('rnuplt',4x,2i10,42x,'intrnu')
  280 FORMAT ('rnuplt',4x,2i10,42x,'ncyrnu')

!........bnuplt

  281 FORMAT ('!      dtimeplot : write to subroutine bnuplot every dtimeplot ms.')
  282 FORMAT ('!      rinnerb : quantities are interpolated to rinnerb for ibound data.')
  283 FORMAT ('!      nplotinnerb : the unit number for ibound data.')
  284 FORMAT ('!      iplotinnerb : inner boundary plot switch.')
  285 FORMAT ('!      routerb : quantities are interpolated to routerb for obound data.')
  286 FORMAT ('!      nplotouterb : the unit number for obound data.')
  287 FORMAT ('!      iplotouterb : outer boundary plot switch.')
  288 FORMAT ('!      r_lumerms : luminosities and rms energies are interpolated to r_lumerms for lbound data.')
  289 FORMAT ('!      nplotlum : the unit number for lbound data.')
  290 FORMAT ('!      iplotlum : luminosities and rms energies plot switch.')
  291 FORMAT ('!      nplotshk : the unit number for shock data.')
  292 FORMAT ('!      iplotshk : shock data plot switch.')
  293 FORMAT ('!      nplotcnv : the unit number for convection data.')
  294 FORMAT ('!      iplotcnv : convection plot switch.')
  295 FORMAT ('!      nplotmss : the unit number for mass data.')
  296 FORMAT ('!      iplotmss : mass data plot switch.')
  297 FORMAT ('bnuplt',29x,1pe15.8,22x,'dtimeplot')
  298 FORMAT ('bnuplt',4x,2i10,5x,1pe15.8,22x,'innerb')
  299 FORMAT ('bnuplt',4x,2i10,5x,1pe15.8,22x,'outerb')
  300 FORMAT ('bnuplt',4x,2i10,5x,1pe15.8,22x,'lumermsplt')
  301 FORMAT ('bnuplt',4x,2i10,42x,'shkplt')
  302 FORMAT ('bnuplt',4x,2i10,42x,'cnvplt')

!........varplt

  311 FORMAT ('!      ivarplt : configuration plot edit switch.')
  312 FORMAT ('!      nvar : the number of cycles since the last implemtation of subroutine varplot.')
  313 FORMAT ('!      nvarint : subroutine varplot is called every nvarint cycles.')
  314 FORMAT ('!      nvarplt : the unit number for editing varplot data.')
  315 FORMAT ('!      nvardump : varplot file number (varplot files are numbered sequentially.')
  316 FORMAT ('!      dtvarplot : write to subroutine varplot every dtvarplot ms.')
  317 FORMAT ('varplt',4x,2i10,42x,'ivarplt')
  318 FORMAT ('varplt',4x,3i10,32x,'nvarint')
  319 FORMAT ('varplt',29x,1pe15.8,22x,'dtvarplot')

!........complt

  321 FORMAT ('!      icomplt: comparison plot edit switch.')
  322 FORMAT ('!      ncomplt: the unit number for editing complot data.')
  323 FORMAT ('!      ncomdump: the complot file number (complot files are numbered sequentially).')
  324 FORMAT ('!      dtcomplot: write to subroutine complot every dtcomplot ms.')
  325 FORMAT ('complt',4x,2i10,42x,'icomplt')
  326 FORMAT ('complt',14x,i10,5x,1pe15.8,22x,'dtcomplot')

!........nurad

  331 FORMAT ('!      nplotnurad ; the unit number for nuradplot data.')
  332 FORMAT ('!      iplotnurad : nurad plot switch.')
  333 FORMAT ('!      dtnuradplot : write to subroutine nuradplot every dtnuradplot ms.')
  334 FORMAT ('!      r_nurad : the radius at which to evaluate neutrino radiation.')
  335 FORMAT ('!      rho_nurad : the density at which to evaluate neutrino radiation.')
  336 FORMAT ('!      nu_r(k,n) : number of neutrinos of energy group k radiated across r_nurad in time dtnuradplot.')
  337 FORMAT ('!      nu_rt(k,n): the cumulative number of neutrinos of energy group k radiated across r_nurad.')
  338 FORMAT ('!      nu_rho(k,n) : number of neutrinos of energy group k radiated across rho_nurad in time dtnuradplot.')
  339 FORMAT ('!      nu_rhot(k,n) : the cumulative number of neutrinos of energy group k radiated across rho_nurad.')
  340 FORMAT ('!      unu_r(k,n) : energy of neutrinos of energy group k radiated across r_nurad in time dtnuradplot.')
  341 FORMAT ('!      unu_rt(k,n) : the- cumulative energy of neutrinos of energy group k radiated across r_nurad.')
  342 FORMAT ('!      unu_rho(k,n) : energy of neutrinos of energy group k radiated across rho_nurad in time dtnuradplot.')
  343 FORMAT ('!      unu_rhot(k,n) : the cumulative energy of neutrinos of energy group k radiated across rho_nurad.')
  344 FORMAT ('nurad ',4x,2i10,42x,'inuplt')
  345 FORMAT ('nurad ',29x,1pe15.8,22x,'dtnuradplot')
  346 FORMAT ('nurad ',29x,1pe15.8,22x,'r_nurad')
  347 FORMAT ('nurad ',29x,1pe15.8,22x,'rho_nurad')
  348 FORMAT ('nurad ',4x,2i5,4(1x,1pe12.5),'nu_r')

!........nudata

  351 FORMAT ('!      nnudata : the unit number for nuplot data.')
  352 FORMAT ('!      inudata : neutrino distribution plot switch.')
  353 FORMAT ('!      dtnudata : write to subroutine nuplot every dtnudata ms.')
  354 FORMAT ('!      r_nudata : the radius at which psi0 and psi1 are evaluated for nuplot.')
  355 FORMAT ('!      t_nudata : the time when the last data dump to nuplot was made.')
  356 FORMAT ('!      psi0dat(k,n) : the time integrated psi0, to be time averaged on dumping.')
  357 FORMAT ('!      psi1dat(k,n) : the time integrated psi1, to be time averaged on dumping.')
  358 FORMAT ('nudata',4x,2i10,42x,'inudata')
  359 FORMAT ('nudata',29x,1pe15.8,22x,'dtnudata')
  360 FORMAT ('nudata',29x,1pe15.8,22x,'r_nudata')
  361 FORMAT ('nudata',29x,1pe15.8,22x,'t_nudata')
  362 FORMAT ('nudata',4x,2i5,2(1x,1pe12.5),26x,'nudata')

!........lagplt

  371 FORMAT ('!      nlagplt : the unit number for lagplot data.')
  372 FORMAT ('!      ilagplt : lagrangian plot switch.')
  373 FORMAT ('!      nlagdump : lagplot file number (lagplot files are numbered sequentially.')
  374 FORMAT ('!      msslag : lagrangian mass of fluid element for lagplot data.')
  375 FORMAT ('lagplt',4x,3i10,32x,'ilagplt')
  376 FORMAT ('lagplt',29x,1pe15.8,22x,'msslag')

!........rlagpt

  381 FORMAT ('!      nrlagplt : the unit number for rlagplot data.')
  382 FORMAT ('!      irlagplt : r-lagrangian plot switch.')
  383 FORMAT ('!      dmlag : the lagrangian mass difference between adjacent lagrangian points to be plotted.')
  384 FORMAT ('rlagpt',4x,2i10,42x,'irlagplt')
  385 FORMAT ('rlagpt',29x,1pe15.8,22x,'dmlag')

!........e_plt

  391 FORMAT ('!      n_eplt : the unit number for e_chk data.')
  392 FORMAT ('!      i_eplt : energy-check plot switch.')
  393 FORMAT ('!      dt_eplot : write to file e_chk.d every dt_eplot ms.')
  394 FORMAT ('e_plt ',4x,2i10,42x,'i_eplt')
  395 FORMAT ('e_plt ',29x,1pe15.8,22x,'dt_eplot')

!........pinch

  401 FORMAT ('!      r_pinch : the radius at which to evaluate the pinch parameters.')
  403 FORMAT ('!      d_pinch : the density at which to evaluate the pinch parameters.')
  405 FORMAT ('pinch ',29x,1pe15.8,22x,'r_pinch')
  407 FORMAT ('pinch ',29x,1pe15.8,22x,'d_pinch')

!cccccc   MD edit parameters

!........MDedu

  501 FORMAT (/'cccccc   MD u-edit parameters'/)
  502 FORMAT ('!      iedMDu : switch for performing MD u-edits.')
  503 FORMAT ('!      nedMDu : u-edit counter.')
  504 FORMAT ('!      intedMDu : u-edit cycles interval.')
  505 FORMAT ('!      n_editMDu : model number counter.')
  506 FORMAT ('MDedu ',14x,i10,42x,'iedMDu')
  507 FORMAT ('MDedu ',14x,i10,42x,'nedMDu')
  508 FORMAT ('MDedu ',14x,i10,42x,'intedMDu')
  509 FORMAT ('MDedu ',14x,i10,42x,'n_editMDu')

!........MDedv

  511 FORMAT (/'cccccc   MD v-edit parameters'/)
  512 FORMAT ('!      iedMDv : switch for performing MD v-edits.')
  513 FORMAT ('!      nedMDv : v-edit counter.')
  514 FORMAT ('!      intedMDv : v-edit cycles interval.')
  515 FORMAT ('!      n_editMDv : model number counter.')
  516 FORMAT ('MDedv ',14x,i10,42x,'iedMDv')
  517 FORMAT ('MDedv ',14x,i10,42x,'nedMDv')
  518 FORMAT ('MDedv ',14x,i10,42x,'intedMDv')
  519 FORMAT ('MDedv ',14x,i10,42x,'n_editMDv')

!........MDedw

  521 FORMAT (/'cccccc   MD w-edit parameters'/)
  522 FORMAT ('!      iedMDw : switch for performing MD w-edits.')
  523 FORMAT ('!      nedMDw : w-edit counter.')
  524 FORMAT ('!      intedMDw : w-edit cycles interval.')
  525 FORMAT ('!      n_editMDw : model number counter.')
  526 FORMAT ('MDedw ',14x,i10,42x,'iedMDw')
  527 FORMAT ('MDedw ',14x,i10,42x,'nedMDw')
  528 FORMAT ('MDedw ',14x,i10,42x,'intedMDw')
  529 FORMAT ('MDedw ',14x,i10,42x,'n_editMDw')

!........MDeds

  531 FORMAT (/'cccccc   MD s-edit parameters'/)
  532 FORMAT ('!      iedMDs : switch for performing MD s-edits.')
  533 FORMAT ('!      nedMDvs : s-edit counter.')
  534 FORMAT ('!      intedMDs : s-edit cycles interval.')
  535 FORMAT ('!      n_editMDs : model number counter.')
  536 FORMAT ('MDeds ',14x,i10,42x,'iedMDs')
  537 FORMAT ('MDeds ',14x,i10,42x,'nedMDs')
  538 FORMAT ('MDeds ',14x,i10,42x,'intedMDs')
  539 FORMAT ('MDeds ',14x,i10,42x,'n_editMDs')

!........MDedd

  541 FORMAT (/'cccccc   MD d-edit parameters'/)
  542 FORMAT ('!      iedMDd : switch for performing MD d-edits.')
  543 FORMAT ('!      nedMDd : d-edit counter.')
  544 FORMAT ('!      intedMDd : d-edit cycles interval.')
  545 FORMAT ('!      n_editMDd : model number counter.')
  546 FORMAT ('MDedd ',14x,i10,42x,'iedMDd')
  547 FORMAT ('MDedd ',14x,i10,42x,'nedMDd')
  548 FORMAT ('MDedd ',14x,i10,42x,'intedMDd')
  549 FORMAT ('MDedd ',14x,i10,42x,'n_editMDd')

!........MDede

  551 FORMAT (/'cccccc   MD e-edit parameters'/)
  552 FORMAT ('!      iedMDe : switch for performing MD e-edits.')
  553 FORMAT ('!      nedMDe : e-edit counter.')
  554 FORMAT ('!      intedMDe : e-edit cycles interval.')
  555 FORMAT ('!      n_editMDe : model number counter.')
  556 FORMAT ('MDede ',14x,i10,42x,'iedMDe')
  557 FORMAT ('MDede ',14x,i10,42x,'nedMDe')
  558 FORMAT ('MDede ',14x,i10,42x,'intedMDe')
  559 FORMAT ('MDede ',14x,i10,42x,'n_editMDe')

!........MDedp

  561 FORMAT (/'cccccc   MD p-edit parameters'/)
  562 FORMAT ('!      iedMDp : switch for performing MD p-edits.')
  563 FORMAT ('!      nedMDp : p-edit counter.')
  564 FORMAT ('!      intedMDp : p-edit cycles interval.')
  565 FORMAT ('!      n_editMDp : model number counter.')
  566 FORMAT ('MDedp ',14x,i10,42x,'iedMDp')
  567 FORMAT ('MDedp ',14x,i10,42x,'nedMDp')
  568 FORMAT ('MDedp ',14x,i10,42x,'intedMDp')
  569 FORMAT ('MDedp ',14x,i10,42x,'n_editMDp')

!........MDedenu

  571 FORMAT (/'cccccc   MD enu-edit parameters'/)
  572 FORMAT ('!      iedMDenu : switch for performing MD enu-edits.')
  573 FORMAT ('!      nedMDenu : enu-edit counter.')
  574 FORMAT ('!      intedMDenu : enu-edit cycles interval.')
  575 FORMAT ('!      n_editMDenu : model number counter.')
  576 FORMAT ('MDeden',14x,i10,42x,'iedMDenu')
  577 FORMAT ('MDeden',14x,i10,42x,'nedMDenu')
  578 FORMAT ('MDeden',14x,i10,42x,'intedMDenu')
  579 FORMAT ('MDeden',14x,i10,42x,'n_editMDenu')

!........MDedfnu

  581 FORMAT (/'cccccc   MD fnu-edit parameters'/)
  582 FORMAT ('!      iedMDfnu : switch for performing MD fnu-edits.')
  583 FORMAT ('!      nedMDfnu : fnu-edit counter.')
  584 FORMAT ('!      intedMDfnu : fnu-edit cycles interval.')
  585 FORMAT ('!      n_editMDfnu : model number counter.')
  586 FORMAT ('MDedfn',14x,i10,42x,'iedMDfnu')
  587 FORMAT ('MDedfn',14x,i10,42x,'nedMDfnu')
  588 FORMAT ('MDedfn',14x,i10,42x,'intedMDfnu')
  589 FORMAT ('MDedfn',14x,i10,42x,'n_editMDfnu')

!........MDeda

  591 FORMAT (/'cccccc   MD a-edit parameters'/)
  592 FORMAT ('!      iedMDa : switch for performing MD a-edits.')
  593 FORMAT ('!      nedMDa : a-edit counter.')
  594 FORMAT ('!      intedMDa : a-edit cycles interval.')
  595 FORMAT ('!      n_editMDa : model number counter.')
  596 FORMAT ('MDeda ',14x,i10,42x,'iedMDa')
  597 FORMAT ('MDeda ',14x,i10,42x,'nedMDa')
  598 FORMAT ('MDeda ',14x,i10,42x,'intedMDa')
  599 FORMAT ('MDeda ',14x,i10,42x,'n_editMDa')

!........MDedx

  601 FORMAT (/'cccccc   MD x-edit parameters'/)
  602 FORMAT ('!      iedMDx : switch for performing MD x-edits.')
  603 FORMAT ('!      nedMDx : x-edit counter.')
  604 FORMAT ('!      intedMDx : x-edit cycles interval.')
  605 FORMAT ('!      n_editMDx : model number counter.')
  606 FORMAT ('MDedx ',14x,i10,42x,'iedMDx')
  607 FORMAT ('MDedx ',14x,i10,42x,'nedMDx')
  608 FORMAT ('MDedx ',14x,i10,42x,'intedMDx')
  609 FORMAT ('MDedx ',14x,i10,42x,'n_editMDx')

!........MDedye

  611 FORMAT (/'cccccc   MD ye-edit parameters'/)
  612 FORMAT ('!      iedMDye : switch for performing MD ye-edits.')
  613 FORMAT ('!      nedMDvye : ye-edit counter.')
  614 FORMAT ('!      intedMDye : ye-edit cycles interval.')
  615 FORMAT ('!      n_editMDye : model number counter.')
  616 FORMAT ('MDedye',14x,i10,42x,'iedMDye')
  617 FORMAT ('MDedye',14x,i10,42x,'nedMDye')
  618 FORMAT ('MDedye',14x,i10,42x,'intedMDye')
  619 FORMAT ('MDedye',14x,i10,42x,'n_editMDye')

!........MDedcm

  621 FORMAT (/'cccccc   MD cm-edit parameters'/)
  622 FORMAT ('!      iedMDcm : switch for performing MD cm-edits.')
  623 FORMAT ('!      nedMDvcm : cm-edit counter.')
  624 FORMAT ('!      intedMDcm : cm-edit cycles interval.')
  625 FORMAT ('!      n_editMDcm : model number counter.')
  626 FORMAT ('MDedcm',14x,i10,42x,'iedMDcm')
  627 FORMAT ('MDedcm',14x,i10,42x,'nedMDcm')
  628 FORMAT ('MDedcm',14x,i10,42x,'intedMDcm')
  629 FORMAT ('MDedcm',14x,i10,42x,'n_editMDcm')

!........MDednu

  631 FORMAT (/'cccccc   MD nu-edit parameters'/)
  632 FORMAT ('!      iedMDnu : switch for performing MD nu-edits.')
  633 FORMAT ('!      nedMDvnu : nu-edit counter.')
  634 FORMAT ('!      intedMDnu : nu-edit cycles interval.')
  635 FORMAT ('!      n_editMDnu : model number counter.')
  636 FORMAT ('MDednu',14x,i10,42x,'iedMDnu')
  637 FORMAT ('MDednu',14x,i10,42x,'nedMDnu')
  638 FORMAT ('MDednu',14x,i10,42x,'intedMDnu')
  639 FORMAT ('MDednu',14x,i10,42x,'n_editMDnu')

!........MDednc

  641 FORMAT (/'cccccc   MD nc-edit parameters'/)
  642 FORMAT ('!      iedMDnc : switch for performing MD nc-edits.')
  643 FORMAT ('!      nedMDvnc : nc-edit counter.')
  644 FORMAT ('!      intedMDnc : nc-edit cycles interval.')
  645 FORMAT ('!      n_editMDnc : model number counter.')
  646 FORMAT ('MDednc',14x,i10,42x,'iedMDnc')
  647 FORMAT ('MDednc',14x,i10,42x,'nedMDnc')
  648 FORMAT ('MDednc',14x,i10,42x,'intedMDnc')
  649 FORMAT ('MDednc',14x,i10,42x,'n_editMDnc')

!........MDednl

  651 FORMAT (/'cccccc   MD nl-edit parameters'/)
  652 FORMAT ('!      iedMDnl : switch for performing MD nl-edits.')
  653 FORMAT ('!      nedMDvnl : nl-edit counter.')
  654 FORMAT ('!      intedMDnl : nl-edit cycles interval.')
  655 FORMAT ('!      n_editMDnl : model number counter.')
  656 FORMAT ('MDednl',14x,i10,42x,'iedMDnl')
  657 FORMAT ('MDednl',14x,i10,42x,'nedMDnl')
  658 FORMAT ('MDednl',14x,i10,42x,'intedMDnl')
  659 FORMAT ('MDednl',14x,i10,42x,'n_editMDnl')

!........MDedne

  661 FORMAT (/'cccccc   MD ne-edit parameters'/)
  662 FORMAT ('!      iedMDne : switch for performing MD ne-edits.')
  663 FORMAT ('!      nedMDvne : ne-edit counter.')
  664 FORMAT ('!      intedMDne : ne-edit cycles interval.')
  665 FORMAT ('!      n_editMDne : model number counter.')
  666 FORMAT ('MDedne',14x,i10,42x,'iedMDne')
  667 FORMAT ('MDedne',14x,i10,42x,'nedMDne')
  668 FORMAT ('MDedne',14x,i10,42x,'intedMDne')
  669 FORMAT ('MDedne',14x,i10,42x,'n_editMDne')

!........MDedgx

  671 FORMAT (/'cccccc   MD gx-edit parameters'/)
  672 FORMAT ('!      iedMDgx : switch for performing MD gx-edits.')
  673 FORMAT ('!      nedMDvgx : gx-edit counter.')
  674 FORMAT ('!      intedMDgx : gx-edit cycles interval.')
  675 FORMAT ('!      n_editMDgx : model number counter.')
  676 FORMAT ('MDedgx',14x,i10,42x,'iedMDgx')
  677 FORMAT ('MDedgx',14x,i10,42x,'gxdMDgx')
  678 FORMAT ('MDedgx',14x,i10,42x,'intedMDgx')
  679 FORMAT ('MDedgx',14x,i10,42x,'n_editMDgx')

!........MDedgy

  681 FORMAT (/'cccccc   MD gy-edit parameters'/)
  682 FORMAT ('!      iedMDgy : switch for performing MD gy-edits.')
  683 FORMAT ('!      nedMDvgy : gy-edit counter.')
  684 FORMAT ('!      intedMDgy : gy-edit cycles interval.')
  685 FORMAT ('!      n_editMDgy : model number counter.')
  686 FORMAT ('MDedgy',14x,i10,42x,'iedMDgy')
  687 FORMAT ('MDedgy',14x,i10,42x,'gydMDgy')
  688 FORMAT ('MDedgy',14x,i10,42x,'intedMDgy')
  689 FORMAT ('MDedgy',14x,i10,42x,'n_editMDgy')

!........MDedBV

  691 FORMAT (/'cccccc   MD BVw-edit parameters'/)
  692 FORMAT ('!      iedMDBVw : switch for performing MD BVw-edits.')
  693 FORMAT ('!      nedMDvBVw : BVw-edit counter.')
  694 FORMAT ('!      intedMDBVw : BVw-edit cycles interval.')
  695 FORMAT ('!      n_editMDBVw : model number counter.')
  696 FORMAT ('MDedBV',14x,i10,42x,'iedMDBVw')
  697 FORMAT ('MDedBV',14x,i10,42x,'gydMDBVw')
  698 FORMAT ('MDedBV',14x,i10,42x,'intedMDBVw')
  699 FORMAT ('MDedBV',14x,i10,42x,'n_editMDBVw')

!........MDedyl

  701 FORMAT (/'cccccc   MD yl-edit parameters'/)
  702 FORMAT ('!      iedMDyl : switch for performing MD yl-edits.')
  703 FORMAT ('!      nedMDvyl : yl-edit counter.')
  704 FORMAT ('!      intedMDyl : yl-edit cycles interval.')
  705 FORMAT ('!      n_editMDyl : model number counter.')
  706 FORMAT ('MDedyl',14x,i10,42x,'iedMDyl')
  707 FORMAT ('MDedyl',14x,i10,42x,'yldMDyl')
  708 FORMAT ('MDedyl',14x,i10,42x,'intedMDyl')
  709 FORMAT ('MDedyl',14x,i10,42x,'n_editMDyl')

!........MD edit time parameters

  761 FORMAT (/'cccccc   MD edit time parameters'/)
  762 FORMAT ('!      i_editMD   : switch for performing MD edits at selected times.')
  763 FORMAT ('!      n_editMD   : model number counter.')
  764 FORMAT ('!      dt_MDedit1 : dump files for MD edits every dt_MDedit1 ms before bounce.')
  765 FORMAT ('!      dt_MDedit2 : dump files for MD edits every dt_MDedit1 ms after bounce.')
  766 FORMAT ('MDedit',14x,i10,42x,'i_editMD')
  767 FORMAT ('MDedit',14x,i10,42x,'n_editMD')
  768 FORMAT ('MDedit',29x,1pe15.8,22x,'dt_MDedit1')
  769 FORMAT ('MDedit',29x,1pe15.8,22x,'dt_MDedit2')

!cccccc   Global edit parameters

!........Global edit cycle parameters

  771 FORMAT (/'cccccc   Global edit cycle parameters'/)
  772 FORMAT ('!      ied_global_n  : switch for performing global edits using cycle criterion.')
  773 FORMAT ('!      nned_global   : global edit counter.')
  774 FORMAT ('!      inted_global  : edit every inted_global cycles.')
  775 FORMAT ('!      n_edit_global : model number counter.')
  776 FORMAT ('ed_gln',14x,i10,42x,'ied_global_n')
  777 FORMAT ('ed_gln',14x,i10,42x,'nned_global')
  778 FORMAT ('ed_gln',14x,i10,42x,'inted_global')
  779 FORMAT ('ed_gln',14x,i10,42x,'n_edit_global')

!........Global time cycle parameters

  781 FORMAT (/'cccccc   Global edit time parameters'/)
  782 FORMAT ('!      ied_global_t  : switch for performing global edits at selected times.')
  783 FORMAT ('!      nted_global   : model number counter.')
  784 FORMAT ('!      dt_global_ed1 : dump files for global edits every dt_global_ed1 ms before bounce.')
  785 FORMAT ('!      dt_global_ed2 : dump files for global edits every dt_global_ed2 ms before bounce.')
  786 FORMAT ('ed_glt',14x,i10,42x,'ied_global_t')
  787 FORMAT ('ed_glt',14x,i10,42x,'nted_global')
  788 FORMAT ('ed_glt',29x,1pe15.8,22x,'dt_global_ed1')
  789 FORMAT ('ed_glt',29x,1pe15.8,22x,'dt_global_ed2')


!cccccc   HDF edit parameters

  790 FORMAT (/'cccccc   HDF edit parameters'/)

!........HDF edit parameters

  791 FORMAT ('!      i_HDFedit : switch for dumping files for HDF edits.')
  792 FORMAT ('!      nd_HDFedit : the unit number for dumping files for HDF edits.')
  793 FORMAT ('!      n_HDFedit : HDF edit counter.')
  794 FORMAT ('!      dt_HDFedit1 : dump files for HDF edits every dt_HDFedit1 ms before bounce.')
  795 FORMAT ('!      dt_HDFedit2 : dump files for HDF edits every dt_HDFedit1 ms after bounce.')
  796 FORMAT ('HDFedt',14x,i10,42x,'i_HDFedit')
  797 FORMAT ('HDFedt',14x,i10,42x,'nd_HDFedit')
  798 FORMAT ('HDFedt',14x,i10,42x,'n_HDFedit')
  799 FORMAT ('HDFedt',29x,1pe15.8,22x,'dt_HDFedit1')
  800 FORMAT ('HDFedt',29x,1pe15.8,22x,'dt_HDFedit2')

!cccccc   Restart file arrays

  801 FORMAT (/'cccccc   Restart file arrays'/)

!........intrst

  802 FORMAT ('!      intrst: the number of cycles between temporary restart dumps.')
  803 FORMAT ('!      nnrst: the number of cycles since the last temporary restart dump.')
  804 FORMAT ('!      irstbgn: restart dump parameter.')
  805 FORMAT ('!      nrstfl: a temporary restart dump parameter.')
  806 FORMAT ('!      ncychg: the number of cycles between temporary restart dumps.')
  807 FORMAT ('intrst',14x,i10,42x,'intrst')
  808 FORMAT ('intrst',14x,i10,42x,'nnrst')
  809 FORMAT ('intrst',14x,i10,42x,'irstbgn')
  810 FORMAT ('intrst',14x,i10,42x,'nrstfl')
  811 FORMAT ('intrst',14x,i10,42x,'ncychg')

!........intprm

  821 FORMAT ('!      intprm: the number of cycles between permanent restart dumps.')
  822 FORMAT ('!      nprm: the number of cycles since the last permanent restart dump.')
  823 FORMAT ('intprm',14x,i10,42x,'intprm')
  824 FORMAT ('intprm',14x,i10,42x,'nprm')

!........ncyrst

  831 FORMAT ('!      ncyrst(i): cycle number for the ith prescribed permanent restart dump.')
  832 FORMAT ('ncyrst',4x,2i10,42x,'ncyrst')

!........Document the dump

 1001 FORMAT (' ***Edit keys dump written at cycle       ',i10,' on unit',i5,'***')

!-----------------------------------------------------------------------
!  Edit arrays
!-----------------------------------------------------------------------

!........Header

WRITE (ndump,13)
WRITE (ndump,1)
WRITE (ndump,2)
WRITE (ndump,3)
WRITE (ndump,15)

!........data_path

WRITE (ndump,13)
WRITE (ndump,21)
WRITE (ndump,15)
WRITE (ndump,22) data_path

!........iflprt

WRITE (ndump,13)
WRITE (ndump,26)
WRITE (ndump,27)
WRITE (ndump,15)
WRITE (ndump,28) iflprt

!........intprt

WRITE (ndump,13)
WRITE (ndump,31)
WRITE (ndump,32)
WRITE (ndump,33)
WRITE (ndump,34)
WRITE (ndump,15)
WRITE (ndump,35) intprt
WRITE (ndump,36) nprt(ij_ray,ik_ray)
WRITE (ndump,37) iprtbgn

!........output

WRITE (ndump,13)
WRITE (ndump,41)
WRITE (ndump,42)
WRITE (ndump,43)
WRITE (ndump,44)
WRITE (ndump,45)
WRITE (ndump,15)
WRITE (ndump,46) iprint
WRITE (ndump,47) nrstd1
WRITE (ndump,48) nrstd2
WRITE (ndump,49) noutpmt
WRITE (ndump,50) nouttmp

!........t_edit

WRITE (ndump,13)
WRITE (ndump,51)
WRITE (ndump,52)
WRITE (ndump,15)
WRITE (ndump,53) it_edit
WRITE (ndump,54) dt_edit

!........plot

WRITE (ndump,13)
WRITE (ndump,61)
WRITE (ndump,62)
WRITE (ndump,63)
WRITE (ndump,64)
WRITE (ndump,65)
WRITE (ndump,66)
WRITE (ndump,67)
WRITE (ndump,15)
WRITE (ndump,68) iplot
WRITE (ndump,69) nplotc
WRITE (ndump,70) nplote
WRITE (ndump,71) nplota
WRITE (ndump,72) nplott
WRITE (ndump,73) intplf
WRITE (ndump,74) npltf

!........model number

WRITE (ndump,13)
WRITE (ndump,81)
WRITE (ndump,15)
WRITE (ndump,82) nmodel

!........rhoprt

WRITE (ndump,13)
WRITE (ndump,91)
WRITE (ndump,15)
WRITE (ndump,92) (i,rhoprint(i),i = 1,10)

!........intedc

WRITE (ndump,13)
WRITE (ndump,101)
WRITE (ndump,102)
WRITE (ndump,103)
WRITE (ndump,15)
WRITE (ndump,104) (i,intedc(i),i = 1,1)
WRITE (ndump,105) (i,nedc(i),i = 1,1)
WRITE (ndump,106) (i,idxedc(i),i = 1,1)

!........intede

WRITE (ndump,13)
WRITE (ndump,111)
WRITE (ndump,112)
WRITE (ndump,113)
WRITE (ndump,15)
WRITE (ndump,114) (i,intede(i),i = 1,1)
WRITE (ndump,115) (i,nede(i),i = 1,1)
WRITE (ndump,116) (i,idxede(i),i = 1,1)

!........intdmi

WRITE (ndump,13)
WRITE (ndump,121)
WRITE (ndump,122)
WRITE (ndump,123)
WRITE (ndump,15)
WRITE (ndump,124) (i,intdmi(i),i = 1,1)
WRITE (ndump,125) (i,nedmi(i),i = 1,1)
WRITE (ndump,126) (i,idxemi(i),i = 1,1)

!........intdma

WRITE (ndump,13)
WRITE (ndump,131)
WRITE (ndump,132)
WRITE (ndump,133)
WRITE (ndump,15)
WRITE (ndump,134) (i,intdma(i),i = 1,3)
WRITE (ndump,135) (i,nedma(i),i = 1,3)
WRITE (ndump,136) (i,idxema(i),i = 1,3)

!........intedh

WRITE (ndump,13)
WRITE (ndump,141)
WRITE (ndump,142)
WRITE (ndump,143)
WRITE (ndump,15)
WRITE (ndump,144) (i,intedh(i),i = 1,9)
WRITE (ndump,145) (i,nedh(i),i = 1,9)
WRITE (ndump,146) (i,idxedh(i),i = 1,9)

!........intdps

WRITE (ndump,13)
WRITE (ndump,151)
WRITE (ndump,152)
WRITE (ndump,153)
WRITE (ndump,15)
WRITE (ndump,154) (i,intdps(i),i = 1,1)
WRITE (ndump,155) (i,nedps(i),i = 1,1)
WRITE (ndump,156) (i,idxeps(i),i = 1,1)

!........intedu

WRITE (ndump,13)
WRITE (ndump,161)
WRITE (ndump,162)
WRITE (ndump,163)
WRITE (ndump,15)
WRITE (ndump,164) (i,intedu(i),i = 1,12)
WRITE (ndump,165) (i,nedu(i),i = 1,12)
WRITE (ndump,166) (i,idxedu(i),i = 1,12)

!........intedy

WRITE (ndump,13)
WRITE (ndump,171)
WRITE (ndump,172)
WRITE (ndump,173)
WRITE (ndump,15)
WRITE (ndump,174) (i,intedy(i),i = 1,6)
WRITE (ndump,175) (i,nedy(i),i = 1,6)
WRITE (ndump,176) (i,idxedy(i),i = 1,6)

!........intdsc

WRITE (ndump,13)
WRITE (ndump,181)
WRITE (ndump,182)
WRITE (ndump,183)
WRITE (ndump,15)
WRITE (ndump,184) (i,intdsc(i),i = 1,2)
WRITE (ndump,185) (i,nedsc(i),i = 1,2)
WRITE (ndump,186) (i,idxesc(i),i = 1,2)

!........intedn

WRITE (ndump,13)
WRITE (ndump,191)
WRITE (ndump,192)
WRITE (ndump,193)
WRITE (ndump,194)
WRITE (ndump,195)
WRITE (ndump,15)
WRITE (ndump,196) (n,intedn(n),n = 1,nnu)
WRITE (ndump,197) (n,nedn(n),n = 1,nnu)
WRITE (ndump,198) (n,idxedn(n),n = 1,nnu)
WRITE (ndump,199) (i,niedn(i),i = 1,40)
WRITE (ndump,200) (i,neden(i),i = 1,40)

!........intdng

WRITE (ndump,13)
WRITE (ndump,211)
WRITE (ndump,212)
WRITE (ndump,213)
WRITE (ndump,15)
DO i = 1,18
  WRITE (ndump,214) (i,intdng(i,n),n,n = 1,nnu)
  WRITE (ndump,215) (i,nedng(i,n),n,n = 1,nnu)
  WRITE (ndump,216) (i,idxeng(i,n),n,n = 1,nnu)
END DO

!........t_edit

WRITE (ndump,13)
WRITE (ndump,221)
WRITE (ndump,222)
WRITE (ndump,223)
WRITE (ndump,15)
WRITE (ndump,224) it_edit,dt_edit

!-----------------------------------------------------------------------
!  Plot arrays
!-----------------------------------------------------------------------

WRITE (ndump,231)

!........lumplt

WRITE (ndump,13)
WRITE (ndump,232)
WRITE (ndump,233)
WRITE (ndump,234)
WRITE (ndump,235)
WRITE (ndump,236)
WRITE (ndump,237)
WRITE (ndump,238)
WRITE (ndump,239)
WRITE (ndump,15)
WRITE (ndump,240) ilumplt
WRITE (ndump,241) nlumplt1,nlumplt2
WRITE (ndump,242) r_lum
WRITE (ndump,243) d_lum
WRITE (ndump,244) nlum
WRITE (ndump,245) (n,intlum(n),n = 1,3)
WRITE (ndump,246) (n,ncylum(n),n = 1,2)

!........enuplt

WRITE (ndump,13)
WRITE (ndump,251)
WRITE (ndump,252)
WRITE (ndump,253)
WRITE (ndump,254)
WRITE (ndump,255)
WRITE (ndump,256)
WRITE (ndump,257)
WRITE (ndump,258)
WRITE (ndump,15)
WRITE (ndump,259) ienuplt
WRITE (ndump,260) nenuplt1,nenuplt2
WRITE (ndump,261) r_e_rms
WRITE (ndump,262) d_e_rms
WRITE (ndump,263) nenu
WRITE (ndump,264) (n,intenu(n),n = 1,3)
WRITE (ndump,265) (n,ncyenu(n),n = 1,2)

!........rnuplt

WRITE (ndump,13)
WRITE (ndump,271)
WRITE (ndump,272)
WRITE (ndump,273)
WRITE (ndump,274)
WRITE (ndump,275)
WRITE (ndump,15)
WRITE (ndump,276) irnuplt
WRITE (ndump,277) nrnuplt
WRITE (ndump,278) nrnu
WRITE (ndump,279) (n,intrnu(n),n = 1,3)
WRITE (ndump,280) (n,ncyrnu(n),n = 1,2)

!........bnuplt

WRITE (ndump,13)
WRITE (ndump,281)
WRITE (ndump,282)
WRITE (ndump,283)
WRITE (ndump,284)
WRITE (ndump,285)
WRITE (ndump,286)
WRITE (ndump,287)
WRITE (ndump,288)
WRITE (ndump,289)
WRITE (ndump,290)
WRITE (ndump,291)
WRITE (ndump,292)
WRITE (ndump,293)
WRITE (ndump,294)
WRITE (ndump,295)
WRITE (ndump,296)
WRITE (ndump,15)
WRITE (ndump,297) dtimeplot
WRITE (ndump,298) nplotinnerb,iplotinnerb,rinnerb
WRITE (ndump,299) nplotouterb,iplotouterb,routerb
WRITE (ndump,300) nplotlum,iplotlum,r_lumerms
WRITE (ndump,301) nplotshk,iplotshk
WRITE (ndump,302) nplotcnv,iplotcnv

!........varplt

WRITE (ndump,13)
WRITE (ndump,311)
WRITE (ndump,312)
WRITE (ndump,313)
WRITE (ndump,314)
WRITE (ndump,315)
WRITE (ndump,316)
WRITE (ndump,15)
WRITE (ndump,317) nvarplt,ivarplt
WRITE (ndump,318) nvar,nvarint,nvardump
WRITE (ndump,319) dtvarplot

!........complt

WRITE (ndump,13)
WRITE (ndump,321)
WRITE (ndump,322)
WRITE (ndump,323)
WRITE (ndump,324)
WRITE (ndump,15)
WRITE (ndump,325) ncomplt,icomplt
WRITE (ndump,326) ncomdump,dtcomplot

!........nurad

WRITE (ndump,13)
WRITE (ndump,331)
WRITE (ndump,332)
WRITE (ndump,333)
WRITE (ndump,334)
WRITE (ndump,335)
WRITE (ndump,336)
WRITE (ndump,337)
WRITE (ndump,338)
WRITE (ndump,339)
WRITE (ndump,340)
WRITE (ndump,341)
WRITE (ndump,342)
WRITE (ndump,343)
WRITE (ndump,15)
WRITE (ndump,344) nplotnurad,iplotnurad
WRITE (ndump,345) dtnuradplot
WRITE (ndump,346) r_nurad
WRITE (ndump,347) rho_nurad

!........nudata

WRITE (ndump,13)
WRITE (ndump,351)
WRITE (ndump,352)
WRITE (ndump,353)
WRITE (ndump,354)
WRITE (ndump,355)
WRITE (ndump,356)
WRITE (ndump,357)
WRITE (ndump,15)
WRITE (ndump,358) nnudata,inudata
WRITE (ndump,359) dtnudata
WRITE (ndump,360) r_nudata
WRITE (ndump,361) t_nudata

!........lagplt

WRITE (ndump,13)
WRITE (ndump,371)
WRITE (ndump,372)
WRITE (ndump,373)
WRITE (ndump,374)
WRITE (ndump,15)
WRITE (ndump,375) nlagplt,ilagplt,nlagdump
WRITE (ndump,376) msslag

!........rlagpt

WRITE (ndump,13)
WRITE (ndump,381)
WRITE (ndump,382)
WRITE (ndump,383)
WRITE (ndump,15)
WRITE (ndump,384) nrlagplt,irlagplt
WRITE (ndump,385) dmlag

!........e_plt

WRITE (ndump,13)
WRITE (ndump,391)
WRITE (ndump,393)
WRITE (ndump,394)
WRITE (ndump,15)
WRITE (ndump,394) n_eplt,i_eplt
WRITE (ndump,395) dt_eplot

!........pinch

WRITE (ndump,13)
WRITE (ndump,401)
WRITE (ndump,403)
WRITE (ndump,15)
WRITE (ndump,405) r_pinch
WRITE (ndump,407) d_pinch

!-----------------------------------------------------------------------
!  MD edit parameters
!-----------------------------------------------------------------------

!........MD u-edit

WRITE (ndump,501)
WRITE (ndump,13)
WRITE (ndump,502)
WRITE (ndump,503)
WRITE (ndump,504)
WRITE (ndump,505)
WRITE (ndump,15)
WRITE (ndump,506) iedMDu
WRITE (ndump,507) nedMDu
WRITE (ndump,508) intedMDu
WRITE (ndump,509) n_editMDu

!........MD v-edit

WRITE (ndump,511)
WRITE (ndump,13)
WRITE (ndump,512)
WRITE (ndump,513)
WRITE (ndump,514)
WRITE (ndump,515)
WRITE (ndump,15)
WRITE (ndump,516) iedMDv
WRITE (ndump,517) nedMDv
WRITE (ndump,518) intedMDv
WRITE (ndump,519) n_editMDv

!........MD w-edit

WRITE (ndump,521)
WRITE (ndump,13)
WRITE (ndump,522)
WRITE (ndump,523)
WRITE (ndump,524)
WRITE (ndump,525)
WRITE (ndump,15)
WRITE (ndump,526) iedMDw
WRITE (ndump,527) nedMDw
WRITE (ndump,528) intedMDw
WRITE (ndump,529) n_editMDw

!........MD s-edit

WRITE (ndump,531)
WRITE (ndump,13)
WRITE (ndump,532)
WRITE (ndump,533)
WRITE (ndump,534)
WRITE (ndump,535)
WRITE (ndump,15)
WRITE (ndump,536) iedMDs
WRITE (ndump,537) nedMDs
WRITE (ndump,538) intedMDs
WRITE (ndump,539) n_editMDs

!........MD d-edit

WRITE (ndump,541)
WRITE (ndump,13)
WRITE (ndump,542)
WRITE (ndump,543)
WRITE (ndump,544)
WRITE (ndump,545)
WRITE (ndump,15)
WRITE (ndump,546) iedMDd
WRITE (ndump,547) nedMDd
WRITE (ndump,548) intedMDd
WRITE (ndump,549) n_editMDd

!........MD e-edit

WRITE (ndump,551)
WRITE (ndump,13)
WRITE (ndump,552)
WRITE (ndump,553)
WRITE (ndump,554)
WRITE (ndump,555)
WRITE (ndump,15)
WRITE (ndump,556) iedMDe
WRITE (ndump,557) nedMDe
WRITE (ndump,558) intedMDe
WRITE (ndump,559) n_editMDe

!........MD p-edit

WRITE (ndump,561)
WRITE (ndump,13)
WRITE (ndump,562)
WRITE (ndump,563)
WRITE (ndump,564)
WRITE (ndump,565)
WRITE (ndump,15)
WRITE (ndump,566) iedMDp
WRITE (ndump,567) nedMDp
WRITE (ndump,568) intedMDp
WRITE (ndump,569) n_editMDp

!........MD enu-edit

WRITE (ndump,571)
WRITE (ndump,13)
WRITE (ndump,572)
WRITE (ndump,573)
WRITE (ndump,574)
WRITE (ndump,575)
WRITE (ndump,15)
WRITE (ndump,576) iedMDenu
WRITE (ndump,577) nedMDenu
WRITE (ndump,578) intedMDenu
WRITE (ndump,579) n_editMDenu

!........MD fnu-edit

WRITE (ndump,581)
WRITE (ndump,13)
WRITE (ndump,582)
WRITE (ndump,583)
WRITE (ndump,584)
WRITE (ndump,585)
WRITE (ndump,15)
WRITE (ndump,586) iedMDfnu
WRITE (ndump,587) nedMDfnu
WRITE (ndump,588) intedMDfnu
WRITE (ndump,589) n_editMDfnu

!........MD a-edit

WRITE (ndump,591)
WRITE (ndump,13)
WRITE (ndump,592)
WRITE (ndump,593)
WRITE (ndump,594)
WRITE (ndump,595)
WRITE (ndump,15)
WRITE (ndump,596) iedMDa
WRITE (ndump,597) nedMDa
WRITE (ndump,598) intedMDa
WRITE (ndump,599) n_editMDa

!........MD x-edit

WRITE (ndump,601)
WRITE (ndump,13)
WRITE (ndump,602)
WRITE (ndump,603)
WRITE (ndump,604)
WRITE (ndump,605)
WRITE (ndump,15)
WRITE (ndump,606) iedMDx
WRITE (ndump,607) nedMDx
WRITE (ndump,608) intedMDx
WRITE (ndump,609) n_editMDx

!........MD ye-edit

WRITE (ndump,611)
WRITE (ndump,13)
WRITE (ndump,612)
WRITE (ndump,613)
WRITE (ndump,614)
WRITE (ndump,615)
WRITE (ndump,15)
WRITE (ndump,616) iedMDye
WRITE (ndump,617) nedMDye
WRITE (ndump,618) intedMDye
WRITE (ndump,619) n_editMDye

!........MD cm-edit

WRITE (ndump,621)
WRITE (ndump,13)
WRITE (ndump,622)
WRITE (ndump,623)
WRITE (ndump,624)
WRITE (ndump,625)
WRITE (ndump,15)
WRITE (ndump,626) iedMDcm
WRITE (ndump,627) nedMDcm
WRITE (ndump,628) intedMDcm
WRITE (ndump,629) n_editMDcm

!........MD nu-edit

WRITE (ndump,631)
WRITE (ndump,13)
WRITE (ndump,632)
WRITE (ndump,633)
WRITE (ndump,634)
WRITE (ndump,635)
WRITE (ndump,15)
WRITE (ndump,636) iedMDnu
WRITE (ndump,637) nedMDnu
WRITE (ndump,638) intedMDnu
WRITE (ndump,639) n_editMDnu

!........MD nc-edit

WRITE (ndump,641)
WRITE (ndump,13)
WRITE (ndump,642)
WRITE (ndump,643)
WRITE (ndump,644)
WRITE (ndump,645)
WRITE (ndump,15)
WRITE (ndump,646) iedMDnc
WRITE (ndump,647) nedMDnc
WRITE (ndump,648) intedMDnc
WRITE (ndump,649) n_editMDnc

!........MD nl-edit

WRITE (ndump,641)
WRITE (ndump,13)
WRITE (ndump,652)
WRITE (ndump,653)
WRITE (ndump,654)
WRITE (ndump,655)
WRITE (ndump,15)
WRITE (ndump,656) iedMDnl
WRITE (ndump,657) nedMDnl
WRITE (ndump,658) intedMDnl
WRITE (ndump,659) n_editMDnl

!........MD ne-edit

WRITE (ndump,661)
WRITE (ndump,13)
WRITE (ndump,662)
WRITE (ndump,663)
WRITE (ndump,664)
WRITE (ndump,665)
WRITE (ndump,15)
WRITE (ndump,666) iedMDne
WRITE (ndump,667) nedMDne
WRITE (ndump,668) intedMDne
WRITE (ndump,669) n_editMDne

!........MD gx-edit

WRITE (ndump,671)
WRITE (ndump,13)
WRITE (ndump,672)
WRITE (ndump,673)
WRITE (ndump,674)
WRITE (ndump,675)
WRITE (ndump,15)
WRITE (ndump,676) iedMDgx
WRITE (ndump,677) nedMDgx
WRITE (ndump,678) intedMDgx
WRITE (ndump,679) n_editMDgx

!........MD gy-edit

WRITE (ndump,681)
WRITE (ndump,13)
WRITE (ndump,682)
WRITE (ndump,683)
WRITE (ndump,684)
WRITE (ndump,685)
WRITE (ndump,15)
WRITE (ndump,686) iedMDgy
WRITE (ndump,687) nedMDgy
WRITE (ndump,688) intedMDgy
WRITE (ndump,689) n_editMDgy
!........2D BVw-edit

WRITE (ndump,691)
WRITE (ndump,13)
WRITE (ndump,692)
WRITE (ndump,693)
WRITE (ndump,694)
WRITE (ndump,695)
WRITE (ndump,15)
WRITE (ndump,696) iedMDBVw
WRITE (ndump,697) nedMDBVw
WRITE (ndump,698) intedMDBVw
WRITE (ndump,699) n_editMDBVw

!........MD gy-edit

WRITE (ndump,701)
WRITE (ndump,13)
WRITE (ndump,702)
WRITE (ndump,703)
WRITE (ndump,704)
WRITE (ndump,705)
WRITE (ndump,15)
WRITE (ndump,706) iedMDyl
WRITE (ndump,707) nedMDyl
WRITE (ndump,708) intedMDyl
WRITE (ndump,709) n_editMDyl


!........MD edit time parameters

WRITE (ndump,761)
WRITE (ndump,13)
WRITE (ndump,762)
WRITE (ndump,763)
WRITE (ndump,764)
WRITE (ndump,765)
WRITE (ndump,15)
WRITE (ndump,766) i_editMD
WRITE (ndump,767) n_editMD
WRITE (ndump,768) dt_MDedit1
WRITE (ndump,769) dt_MDedit2

!-----------------------------------------------------------------------
!  Global edit parameters
!-----------------------------------------------------------------------

!........Global edit cycle parameters

WRITE (ndump,771)
WRITE (ndump,13)
WRITE (ndump,772)
WRITE (ndump,773)
WRITE (ndump,774)
WRITE (ndump,775)
WRITE (ndump,15)
WRITE (ndump,776) ied_global_n
WRITE (ndump,777) nned_global
WRITE (ndump,778) inted_global
WRITE (ndump,779) n_edit_global

!........Global edit time parameters

WRITE (ndump,781)
WRITE (ndump,13)
WRITE (ndump,782)
WRITE (ndump,783)
WRITE (ndump,784)
WRITE (ndump,785)
WRITE (ndump,15)
WRITE (ndump,786) ied_global_t
WRITE (ndump,787) nted_global
WRITE (ndump,788) dt_global_ed1
WRITE (ndump,789) dt_global_ed2

!-----------------------------------------------------------------------
!  HDF edit parameters
!-----------------------------------------------------------------------

WRITE (ndump,790)

!........HDFedit

WRITE (ndump,13)
WRITE (ndump,791)
WRITE (ndump,792)
WRITE (ndump,793)
WRITE (ndump,794)
WRITE (ndump,795)
WRITE (ndump,15)
WRITE (ndump,796) i_HDFedit
WRITE (ndump,797) nd_HDFedit
WRITE (ndump,798) n_HDFedit
WRITE (ndump,799) dt_HDFedit1
WRITE (ndump,800) dt_HDFedit2

!-----------------------------------------------------------------------
!  Restart arrays
!-----------------------------------------------------------------------

WRITE (ndump,801)

!........intrst

WRITE (ndump,13)
WRITE (ndump,802)
WRITE (ndump,803)
WRITE (ndump,804)
WRITE (ndump,805)
WRITE (ndump,806)
WRITE (ndump,15)
WRITE (ndump,807) intrst
WRITE (ndump,808) nnrst
WRITE (ndump,809) irstbgn
WRITE (ndump,810) nrstfl
WRITE (ndump,811) ncychg

!........intprm

WRITE (ndump,13)
WRITE (ndump,821)
WRITE (ndump,822)
WRITE (ndump,15)
WRITE (ndump,823) intprm
WRITE (ndump,824) nprm

!........ncyrst

WRITE (ndump,13)
WRITE (ndump,831)
WRITE (ndump,15)
WRITE (ndump,832) (i,ncyrst(i),i=1,100)

!-----------------------------------------------------------------------
!  Record the dump
!-----------------------------------------------------------------------

WRITE (nprint,1001) ncycle,ndump

RETURN
END SUBROUTINE edit_write_adios
