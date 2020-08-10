   ! AutoHetoDet 0D Opt Model
   !  This is mostly being develop for the Entropy Manuscript for Jan 2020.  See notes,
   !  but basically this model is for a chemostat (0D) with phytoplankton (Auto),
   !  consumers (Het) and detritivours (Det, i.e., bacteria).  This is based
   !  on the Darwin-MEP-1D in the PIE-MEP-1D solution, but BiM will replace BACOLI
   !  since only an ODE solution is needed.
   !
   !  Also see MEPdarwin1 in SimpleBoxModels solution, as this is the one used in proposals
   !
   ! V1.0, 11-Dec-2019
   !  -  This is branched from AutoHetDet_V1.2 in the AutoHetDet_0D project
   !     This version is looking at changing the rediculous bland random search
   !     with an optimization based search.  This version is using the BOBYQA,
   !     and omg_cc is no longer treated at binary matrix, but values can
   !     vary between 0 and 1. There is no constraint on nLinks now, so
   !     some cc can be fully connected, while others can have none.
   !     MPI is not really being used here, but the MPI write are being
   !     left in for now, and there may be some use later (perhaps to investigate
   !     multiple starting points.
   ! V2.0, 15-Dec-2019
   !  -  Version 1 worked pretty well, so changing this to use hyperBOB which will
   !     search with MEP processes starting from different locations in a random
   !     latin hypercube.
   ! V2.1, 23-Dec-2019
   !  -  Removed "+ +" from wDen in Gz_bioS
   !  -  BUG: wDen used cc(1:ncc) instead of bac(1:nbac).  Not sure why this didn't generate error
   ! V3.0, 2-Jan-2020
   !  -  Same as V2.1, but adding two more "traits" to represent omg_pp and its time variation (freq and phase)
   !     Consequently, omg_pp(i) is replaced by f_pp(i) and phi_pp(i) (so just one increase in DOF compared to V2.1)
   !     That is, omg_pp(t) = (sin(2 pi f t + phi) + 1)/2
   !     where 0 <= f <= 2  (`1/d)
   !     and   0 <= phi < 2 pi
   ! V3.1, 6-Jan-2020
   !  -  Added time and date for error output associated with BiM
   !  -  Added warning on ODE integration failure while storing each solution.
   ! V3.2, 13-Jan-2020
   !  -  Added f_pp_max/min and omg_pp_max/min variables to parameters
   !  -  Using new v1.4 of hyperBOB (fixes problem with rhoBMin determination)
   !  -  Also, added scaling of phi_pp and f_pp to make it vary between 0 and 1 (not 0 and 2 pi, etc)
   ! V3.3, 14-Jan-2020
   !  -  Decided to noramlize all the control variables (i.e., traits)
   !  -  Note, while v3.2 forced f_pp to zero in code, this can be effectively accomplished by setting
   !     f_pp_min = 0 and f_pp_max = 1.e-10_mp (or something similarly small).
   ! V3.4, 16-Jan-2020
   !  -  Changed output of epp, ecc and ebac to allow tecplot plotting with n = 1
   ! V3.5, 22-Jan-2020
   !  -  Introduced two more control variables to alter omg_pp.  m_pp, which is the mean of omg_1,pp
   !     and amp_pp which specifies the amplitude of the time varying component of omg_1,pp.
   ! V3.6, 28-Jan-2019
   !  -  Returning to V3.4, as the addition of m_pp and amp_pp did not improve solutions
   !  -  Here, the objective is to keep the ration of C_pp to pp equal to a constant, k_pp,
   !     which is a new parameter.  Note, the feed, c_ppI, is adjusted to match k_pp, so is ignored.
   ! V3.7, 26-Feb-2020
   !  -  This is basically v3.6, but by setting k_pp to 0, the code runs v3.4 with the variable
   !     stoichiometry.
   !  -  Added ability to ramp up or down the dilution rate and the NH3 input.
   !  -  ***** Bug found.  For some reason dil was declared as an integer!!!! All previous versions
   !     have this error.  Since dil was set to 0.5, all dil have been run at dil = 0 then.  
   !  -  Added info on PDE fail and TO to output
   ! V3.8, 14-Apr-2020
   !  -  Changed the definition for sumSigmaRxn to NOT include the light dissipation by the photosynthetic part of 
   !     phytoplankton when growth is constrained.  Instead, add that to sumSigPart.  
   !  -  Also added a new vector, sumSigWeights, that allows the entropy production optimization to include any of the three terms
   !     rxn, part, h2o. 
   ! V3.9, 30-Apr-2020
   !  -  Removed C_D and C_P from particulate matter that intercepts light.
   ! V4.0, 9-May-2020
   !  -  Replacing the sin function on omg_pp with a square wave one.  The traits f_pp and phi_pp are replaced with
   !     three traits: tOn_pp, tOff_pp, omg1Amp_pp, and one parameter sigOmg_pp. 
   !     tOn_pp_min, tOn_pp_max, tOff_pp_min, tOff_pp_max and sigOmg_pp must be set in the *.ini file. 
   !  -  Aslo added storage of omp_pp in the misc file
   ! V4.1, 11-May-2020
   !  -  Change the control on omg1_pp so that it turns on when light level is above iOn_pp, where iOn_pp is in units of 
   !     uE/m2/s (note, the model uses I in uE/m2/d, so that is converted to 1/s.  sigOmg_pp is till used to set the step up or down., 
   ! V4.2, 14-May-2020
   !  -  This continues from V4.0 (not 4.1 branch). Adding additional parameter, k_chla, that use used in light extinction
   !     to capture the greater light capturing aspects of of Chl a.  See Wozniak2007 Light Absorption in Sea Water
   ! V4.3, 19-May-2020
   !  -  Added C_P back to light capture, but with k_p, not k_chla
   ! V4.4, 4-Jun-2020
   !  -  This version adds ability to maximize entropy production over a window within the overall simulation period. For instance
   !     the simulation may run from day 0 to day 730, but entropy could be maximized over just one days, such as 230 to 231.
   !  -  Corrected O2 balance for bacteria contribution (was being added instead of subtracted).
   !  -  Added summary file
   ! V4.5, 16-Jun-2020
   !  -  AutoHetDet_Surv_V1.2 was merged into this program so it can be used to generate 2D optimum surfaces too now.
   ! V4.6, 29-Jun-2020
   !  -  Now if optimize is set to FALSE, then random values are chosen for the traits and the ODE intergration run once based on the random value
   !  -  Added sigmaSI in stepInput routine to parameter file so it can be changed.  (8-Jul-2020)
   ! V4.7, 13-Jul-2020
   !  -  Introduce binaryOMG that when true turns omg_cc into a pure connection matrix, with only 0's or 1's.  To do this simply
   !     omg_cc is generated as before, but then anint is applied so that all values < 0.5 are driven to 0, and all >= 0.5 are set to 1.0 if binaryOMG is true.
   !  -  Fixed timeout in ODE integration.  It turns out there is a bug in BiM on line 1212, in that it returns if ierr is set in the ode, but it does not set
   !     idid to -6.  Also, in newton and jacobian calculations, BiM just assumes it's a one time error, and keeps going, when it should really just return.
   !     Consequently, instead of checking IDID, just calculated ODE time and issue a error based on that.  This way, all TO errors will be identified.

  module realPrec
      ! this sets the precision of reals
      implicit none
      integer, parameter:: sp = kind(1.0E0)
      integer, parameter:: dp = kind(1.0D0)
      integer, parameter:: mp = dp  ! Set to either sp or dp
      save      
   end module realPrec

   module interfaces
      interface
         subroutine omega_free (n, w, omg)
            use realPrec
            implicit none
            integer, intent(in)  :: n        ! Length of omg
            real(dp), intent(in) :: w(:)     ! lower case omega, control variables
            real(dp), intent(out):: omg(:)   ! upper case omega, these sum to 1, which is where the nth omega comes from
         end subroutine omega_free
      end interface
      interface
         pure subroutine insertion_sort(n, a)
            use realPrec
            implicit none
            integer, intent(in)                    :: n
            real(dp), intent(in out), dimension(:) :: a
         end subroutine insertion_sort
      end interface   
   end module interfaces
   
   module globalVars
      !  This contains a module for global variables
      use realPrec
      use mpi
      implicit none
      
      ! misc parameters
      character(len=2), parameter:: CRLF = char(13)//char(10) ! used by MPI writes for error output
      
      integer, parameter:: nconc = 6 ! number of concentration variables (dic, O2, NH3, C_L, C_D, N_D)
      
      integer mpiFHerr     ! Unit for MPI writing errors to for BiM, and is set by MPI_FILE_OPEN in the initialization routine.      
      integer iunit_state  ! Tecplot output file for state variables.  
      integer iunit_traits ! unit to be assigned the file containing the epsilon and omgega inputs
      integer iunit_rxn    ! units to be assigned for storing reaction rates
      integer iunit_sigDot ! unit for entropy production sumed for all pp and cc, but rxn, h2o and particles kept seprate
      integer iunit_sigTot ! unit for total entropy production from the three processes over time.
      integer iunit_epp    ! Saving iterations of epp
      integer iunit_omgpp  ! Saving iterations of omg_pp, but actuall save f_pp and phi_pp, which omg_pp is calculated from.
      integer iunit_ecc    ! saving ecc
      integer iunit_omgcc  ! omg_cc
      integer iunit_ebac   ! saving ebac
      integer iunit_omgbac ! omg_bac      
      integer iunit_misc   ! various output variables
      integer iunit_input  ! Saves how inputs change over time.  This is only saved once because it's the same for all solutions.
      integer iunit_summary! Save summmary info, such as wall time, function calls, sumSigTot
      integer iunit_EPsurf ! The 2D sigma_rxn surface      
   
      ! working variables
      integer nstates ! total number of state variables, which includes one for total integrated entropy production. 
      integer myRank ! for MPI
      type solInfo
         ! used to pass information about ODE integration
         integer fail ! = 1 if ODE integration failed; = 2 if time exceeded
         integer idid ! idid value from BACOLI
         real(mp) t0  ! Integration time where failure occured
      end type solInfo
   
      ! Control variables
      real(mp), allocatable:: epp(:), ecc(:), ebac(:) ! epsilon control variables for pp, cc and bac
      real(mp), allocatable:: omg_pp(:,:)  ! matrix of pp partitioning between co2 fixation and growth, which is now calculated!
      !real(mp), allocatable:: f_pp(:), phi_pp(:) ! New traits to calculated omg_pp over time.
      real(mp), allocatable:: tOn_pp(:), tOff_pp(:), omg1Amp_pp(:)
      real(mp), allocatable:: omg_cc(:,:)  ! matrix controlling predation. This is binary, either 0 or 1.  Could use sparce matrix, not yet. 
      real(mp), allocatable:: omg_bac(:,:) ! matrix controlling biomass synthesis, c_d or n_d decomposition (this must sum to unity)     
      ! state variables
      real(mp) dic, dic_lg
      real(mp) o2, o2_lg
      real(mp) nh3, nh3_lg
      real(mp) c_l, c_l_lg
      real(mp) c_d, c_d_lg
      real(mp) n_d, n_d_lg
      real(mp) P_D, P_D_lg ! not currently used (detrital P), but needs to be set in parameters
      !$OMP threadprivate (dic, dic_lg, o2, o2_lg, nh3, nh3_lg,c_l, c_l_lg, c_d, c_d_lg, n_d, n_d_lg, P_D, P_D_lg)  
      real(mp) sumSigTot(3) ! Time integrate entropy for Reactions, water and particles, integrated in BiM
      !$OMP threadprivate (sumSigTot)
      real(mp) hco3 ! obtained from carbonate routine
      real(mp) co2  ! Free CO2 (h2co3) obtained from carbonate routine
      real(mp) alk  ! currently appoximated from salinity
      real(mp) h3po4, h3po4_lg ! this is held fixed
      !$OMP threadprivate (alk, hco3, co2, h3po4, h3po4_lg)
      real(mp), allocatable:: pp(:), pp_lg(:) ! concentration of primary produces, and their max values
      real(mp), allocatable:: C_pp(:), C_pp_lg(:) ! internal C stores for each pp
      real(mp), allocatable:: cc(:), cc_lg(:) ! concentration of consumers, and their max values
      real(mp), allocatable:: bac(:), bac_lg(:)
      !$OMP threadprivate (pp, pp_lg, C_pp, C_pp_lg, cc, cc_lg)
      ! Reaction rate associated variables
      real(mp), allocatable:: rpp(:,:)! rpp reactions by produces
      real(mp), allocatable:: rcc(:,:) ! rcc reactions by consumers
      real(mp), allocatable:: rbac(:,:) ! rbac reactions by bacteria
      real(mp), allocatable:: n2_pp(:) ! stoichiometric coef in r2 equation of pp synthesis
      real(mp), allocatable:: aCi_cc(:) ! stoichiometric coef in cc synthesis
      real(mp) aA1_bac ! stoichiometric coef in bac synthesis
      !$OMP threadprivate (rpp, rcc, rbac, n2_pp, aCi_cc, aA1_bac)
      ! Variables for entropy production and its integrating over t and x
      real(mp) :: sumSigdot(3) ! entropy production from rxns and light absorption from water and particles, calculated in rxnRates
      !$OMP threadprivate (sumSigdot) ! Don't need to make these thread private
      real(mp) :: intSigDotXPrev(3) ! entropy production from rxns and light absorption from water and particles from previous time point
      real(mp) :: tPrev
      real(mp) :: bestEPever = 0._mp ! best solution found from all iterets
      ! Parameters
      real(mp) absZero ! value considered zero
      logical readeps
      real(mp) xc ! characteristic lenght scale (m)
      integer:: ompThreads = 1 ! in the inp file, this can be set to other values, but will default to 1 if not used.
   
      ! Initial and BCs
      ! Note, initial conditions are specified in the values given to the named state variables in the /params/ namelist.
      ! Boundary conditions
      real(mp) dicI ! Initial condition and feed concentration, dic (mmol/m3)
      real(mp) o2I  ! oxygen
      real(mp) nh3I ! ammonium
      real(mp) c_lI ! label organic carbon
      real(mp) c_dI ! C detritus
      real(mp) n_dI ! N detritus
      real(mp) alf, bet, gam, del ! elemental composition of phytoplankton, unit carbon based (H, O, N, P)
      real(mp) cell_F, delPsi
      real(mp) I0max, dLat
      real(mp) nuStar, kappa ! the values for maximum growht rate and 1/2 saturation for adaptive monod equation
      real(mp) nuDet ! maximum decomposition rate of detrital material
      real(mp) dGr_Ggamma ! Gibbs free energy for green photons J/mmol
      real(mp) k_w, k_p ! light absorption by water (1/m) and particles (1/m/uM)
      integer iseed ! used to set the seed for random_seed, and is also used by hyperBOB now
      integer, allocatable:: seed(:) ! used to set the random_seed.  
      
      namelist /params/ absZero
      namelist /params/ dicI, o2I, c_lI, c_dI, n_dI !  initial conditions nh3I, 
      namelist /params/ iseed
      namelist /params/ readeps ! if true, the the file basename.e is read to get epp and ecc.  
      namelist /params/ alf, bet, gam, del ! elemental composition of phytoplankton, unit carbon based (H, O, N, P).
      namelist /params/ cell_F ! concentration factor for C_P. Intracellular versus extracellular volume.  
                               ! In theory this changes as phyS conc does, but keep as constant for now.      
      namelist /params/ delPsi ! parameter in LaRowe2012 thermodynamic driver function (volts)
      namelist /params/ I0max, dLat ! Used for solar simulation.  I0Max is in mmol/m2/d
      namelist /params/ h3po4, P_d ! These are currently fixed parameters.
      namelist /params/ nuStar, kappa ! the values for maximum growht rate and 1/2 saturation for adaptive monod equation
      namelist /params/ nuDet ! maximum decomposition rate of detrital material
      namelist /params/ dGr_Ggamma ! Gibbs free energy for green photons J/mmol
      namelist /params/ k_w, k_p ! light absorption by water (1/m) and particles (1/m/uM)
      namelist /params/ ompThreads
      
      real(mp) dil ! dilution rate (1/d)
      !namelist /params/ dil
      
      ! v3.7, adding parameters and variables to allow stepping of dilution rate, nh3I and c_LI over time
      real(mp) dil_t0, dil_tf ! The initial and final values for dil over simulation
      integer dil_n ! number of steps between t0 and tf
      namelist /params/ dil_t0, dil_tf, dil_n
      real(mp) nh3I_t0, nh3I_tf
      integer nh3I_n 
      namelist /params/ nh3I_t0, nh3I_tf, nh3I_n
      
      real(mp) pCO2, pO2 ! partial pressures of CO2 and O2 in the air (atm)
      real(mp) pV_co2, pV_o2 ! piston velocities (m/d)
      namelist /params/ pCO2, pO2, pV_co2, pV_o2
      
      ! Variables associated with evolution of traits
      real(mp) maxTraitIters ! The maximum number of iterations to conduct
      real(mp) epp_min, epp_max ! min and max values for epp
      real(mp) ecc_min, ecc_max ! min and max values for ecc 
      real(mp) ebac_min, ebac_max ! min and max values for ebac 
      namelist /params/ maxTraitIters, epp_min, epp_max, ecc_min, ecc_max, ebac_min, ebac_max
      real(mp) sigEps_pp, sigEps_cc, sigEps_bac ! exponents for eps generation
      namelist /params/ sigEps_pp, sigEps_cc, sigEps_bac
      ! Associated with optimization-based solution
      !real(mp) epp_ini  ! initial value for all epp
      !real(mp) f_pp_ini ! initial value for all f_pp
      !real(mp) phi_pp_ini ! intial value for all phi_pp
      !real(mp) ecc_ini  ! initial value for ecc
      !real(mp) omg_cc_ini ! initial value for all omg_cc (there is no w version for omg_cc)
      !real(mp) ebac_ini ! initial value for all ebac
      !real(mp) w_bac1_ini ! initial value for all w_bac1 (connected to omg_bac)
      !real(mp) w_bac2_ini ! initial value for all w_bac2 (connected to omg_bac). NOTE w_bac2 >= w_bac1
      !namelist /params/ epp_ini, f_pp_ini, phi_pp_ini, ecc_ini, omg_cc_ini, ebac_ini, w_bac1_ini, w_bac2_ini
      
      ! v3.2 addition (removed in V4.0)
      !real(mp) f_pp_min, f_pp_max ! min and max bounds on f_pp
      !real(mp) phi_pp_min, phi_pp_max ! min and max bounds on phi_pp
      !namelist /params/ f_pp_min, f_pp_max, phi_pp_min, phi_pp_max
      
      ! v3.6
      real(mp):: k_pp = 0.0 ! ratio of c_pp to pp.  V3.7: Now, if k_pp is set to 0, then pp and C_pp are decoupled.
      namelist /params/ k_pp
      
      ! v3.7
      real(mp), dimension(3):: sumSigWeights = [1._mp, 0._mp, 0._mp] ! For EP optimization, these are the weights on Rxn, H20 and particles
      namelist /params/ sumSigWeights
      
      ! v4.0
      real(mp) tOn_pp_min, tOn_pp_max, tOff_pp_min, tOff_pp_max, sigOmg_pp
      namelist /params/ tOn_pp_min, tOn_pp_max, tOff_pp_min, tOff_pp_max, sigOmg_pp
      
      ! V4.2
      real(mp) k_chla ! light extinction coefficient for Chl a (1/m/mM C) Note, this is based on pp conc, not Chl a.
      namelist /params/ k_chla
      
      ! V4.4
      real(mp) t0_ep, tf_ep ! entropy will be maximized over t0_ep to tf_ep.  These times must be within or equal to t0 and tf
      namelist /params/ t0_ep, tf_ep
      
      ! V4.5 Parameters ported over from AutoHetDet_Surf_V1.2
      real(mp) :: bestEPfound = 0._mp ! best solution found from all iterets
      integer  :: bestI, bestJ      
      logical:: genSurf = .false. ! if set to true, a 2D surface based on supplied parameters will be produced instead of running the optimization
      integer nSurfPts    ! Number of points in the x and y dimension of the 2D surface to produce
      real(mp) reportTime ! How often to update during problem in min.
      integer whichPP     ! which of the possible pp's to generate the 2D surface over (if using more than one pp).
      namelist /params/ genSurf, nSurfPts, reportTime, whichPP
      
      ! V4.6 Add sigmaSI, which controls step up or downs in D and NH3I to parameter file input.  Default is old value
      real(mp):: sigmaSI = 5._mp ! Step will take about 1 day at sigmaSI = 10. units of (1/d) 
      namelist /params/ sigmaSI
      
      ! V4.7 Add binaryOMG to set omg_cc to binary matrix (only 0's or 1's)
      logical:: binaryOMG = .false. ! default
      namelist /params/ binaryOMG
            
      ! ODE constrain
      real(mp) maxODEtime ! maximum time allowed to solve the ODE = (tf-t0)/minComFac in days
      real(mp) ODE_t0 ! time at the start the ODE solution
      real(mp) minCompFac ! if a process takes longer than (tf-t0)/minCompFac, then it is terminated 
      namelist /params/ minCompFac
            
      ! BiM parameters   
      real(mp) t0, tDays ! the start time and number of days to run for ODE integration (d)
      integer maxstep_BiM ! maximum number of BiM iterations (set to 0 to use default of 100000)
      integer maxattempts ! number of attempts to solve ODEs before declaring failure
      namelist /params/ t0, tDays, maxstep_BiM, maxattempts
      integer useOmpJac ! set to 0 to have BiM calculate numerical gradient or 1 to use openMP-based numerical gradient
      namelist /params/ useOmpJac 
      real(mp) hmax_BiM ! largest integration step size (use to handle high frequency drivers) DEFAULT = (TEND-T0)/8
      namelist /params/ hmax_BiM
   
      integer npp ! number of primary produces in model
      integer ncc ! number of consumers (predators) in model 
      namelist /params/ npp, ncc ! this is placed in a separate namelist so that they can be used to allocate arrays.
      real(mp) ppI, c_ppI, ccI, bacI ! initial values and feed input for all pp and all cc and all bac
      real(mp) atol1, rtol1 ! scalar absolute and relative tolerance for BiM ODE integration
      namelist /params/ ppI, c_ppI, ccI, bacI
      namelist /params/ atol1, rtol1
      integer nbac ! number of bacterial specie
      namelist /params/ nbac
      
      ! Parameter needed for reaction rate calculations, etc
      real(mp) surA ! surface area (m^2).  If set to 1 m^2, then results can be though of as per m2
      real(mp) T_K ! get temperature in K
      real(mp) pH  ! get pH
      real(mp) depth ! average depth (m)
      real(mp) is ! use the ionic strength
      namelist /params/  surA, T_K, pH, depth, is
      
      ! parameters used by BOBYQA
      integer npt_bobyqa ! !Number of interpolation points between [N+2,(N+1)(N+2)/2]. best if < 2*N+1 In this version set to 2n+1 for now
      real(mp) rhobeg, rhoend ! initial and final values of a trust region radius
      integer iprint ! controls amount of printing (0, 1, 2 or 3)
      integer maxfun ! maximum number of calls to CALFUN
      logical optimize ! If true, MEP optimization occurs, otherwise just solve ODEs
      integer fcnUpdate ! After every fcnUpdate ODE integrations, infor is printed out.
      namelist /params/ rhobeg, rhoend, iprint, maxfun, optimize, fcnUpdate  
      integer fcnCalls, fcnCallsTot ! Not a parameter, but a counter
      integer ODEfailed, ODEtimedOut ! number of times the PDE integration failed or timedout
      real(mp):: fcnMax = 0._mp
      
      ! Global variables
      real(mp)  tf ! calculated end time.
      character basename*80 !filename w/o extenstion for input and output data.   
      ! calculate necessary thermodynamic variables, such as standard free energy's of formation (since many routines will need)  See thermoData module
      real(mp) dGf0_nh3
      real(mp) dGf0_h3po4
      real(mp) dGf0_h2o
      real(mp) dGf0_ch2o ! glucose, single carbon
      real(mp) dGf0_bioS ! free energies of formation of species at current temp (K), is (M), and pH. (kJ/mol)
      real(mp) dGf0_C_D ! Detritus C (use ch2o), N (use nh3) and P (use h3po4).
      real(mp) dGf0_N_D
      real(mp) dGf0_P_D
      real(mp) dGf0_o2aq
      real(mp) dGf0_dic ! or DIC
      !$OMP threadprivate (dGf0_nh3, dGf0_h3po4, dGf0_h2o, dGf0_ch2o, dGf0_bioS, dGf0_C_D, dGf0_N_D, dGf0_P_D, dGf0_o2aq, dGf0_dic)
      real(mp) kwp ! light absorption by all factors
      real(mp) Iht ! light a bottom of water column. This is a global variable for saving
      real(mp) aveI ! average light intensity (mmol photons /m^2 /d (Not micromoles))
      real(mp) aA2_pp, bA2_pp ! stoichiometric coef for reaction 2 of pp
      !$OMP threadprivate (kwp, Iht, aveI, aA2_pp, bA2_pp)  
      integer lenTS ! length of traitSpace      
      real(mp) WTime_it0 ! used as the global start time
      real(mp), parameter:: pi = 3.1415926535897932_mp
      real(mp):: aveODEtime = 0._mp ! Average PDE integration time
      
   contains
      subroutine initialize ()
         integer ioerr, ioE, mpierr, i
         integer npp_t, ncc_t, nbac_t         
         character (len=2000) string
         integer IUNITA(5) ! needed by SLATEC error handling routine.  See call to XSETUA below   
         ! This routine intializes and allocates variables and open/reads/writes files
         nstates = nconc + 2*npp + ncc + nbac + 3 ! total number of state variables (last three is integrated entropy production for reactions, water and particles)
         allocate( epp(npp), ecc(ncc), ebac(nbac) ) ! epsilon for pp and cc and bac, which does not depend on location.
         allocate( omg_pp(npp,2), omg_cc(ncc,npp+ncc+nbac) ) ! omg_cc specifies which predators each which prey, and includes cannabolism 
         !allocate( f_pp(npp), phi_pp(npp) ) ! New traits to calculate omg_pp over time.
         allocate( tOn_pp(npp), tOff_pp(npp), omg1Amp_pp(npp) )
         allocate( omg_bac(nbac,3) ) 
         !$OMP parallel
         ! Since these are threadprivate, they need to be allocated in those processes
         allocate( pp(npp), C_pp(npp) ) ! concentration of primary produces and associated internal C
         allocate( pp_lg(npp), C_pp_lg(npp) )
         allocate( cc(ncc) ) ! concentration of consumers
         allocate( cc_lg(ncc) )
         allocate( bac(nbac), bac_lg(nbac) )
         allocate( rpp(npp,2) ) ! reactions for PP, two for each
         allocate( rcc(ncc,npp+ncc+nbac) ) ! each consumer can potentially eat all pp, and all cc and all bac (stored in the matrix in that order)
         allocate( rbac(nbac,3) ) ! reactions for bac, three for each
         allocate( n2_pp(npp), aCi_cc(ncc) ) ! stoichiometric coef needed for ODE balance equations
         !$OMP end parallel
      
         ! BiM space
         ! open up file using MPI for error output from each process.  This is where mpiFHerr gets set.
         ! use MPI_MODE_WRONLY for write only, creat it if necessary, and for sequential access
         ! However, there is no easy way to clear out a file if it already exists, so open with close_on_delete, to get rid of trash, then reopen
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL+MPI_MODE_DELETE_ON_CLOSE, &
                  MPI_INFO_NULL, mpiFHerr, mpiErr)
         call MPI_FILE_CLOSE(mpiFHerr, mpiErr) ! This will delete the file.    
         call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(basename)//'_err.dat', MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_SEQUENTIAL, &
                  MPI_INFO_NULL, mpiFHerr, mpiErr)
         ! write a header to the file
         if (myRank == 0) then
            string = 'Errors associated with BiM integration:' 
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)
         end if   
      
         ! open file to store (or read if option set) epp, f_pp, phi_pp, ecc, omg_cc, ebac, omg_bac that are randomly generated
         ! Setup the random number generator so that it is reproducible. This is not needed in this version, but leave it here.
         call RANDOM_SEED (size=i) ! get number of integers random_seed uses
         allocate ( seed(i) )
         ! set a seed value based in iseed in the parameters file and myRank so that each process follows a different random 
         ! number sequence, but this allows a run to be repeated and all processes to generate different values.
         seed = iseed + myRank 
         call random_seed (put=seed(1:i))
         deallocate ( seed )
         
         if (readeps) then
            ! Open read only so all processes can read.
            open(newunit=iunit_traits, file=trim(basename)//'_traits.dat',action='read',status='old',iostat=ioerr)   
            if (ioerr == 0) then                
               ! read existing epp and ess from *.e file.  First line is just a header, so skipe that
               read(iunit_traits,'(a)') string
               read(iunit_traits,*) npp_t, ncc_t, nbac_t ! code should check that these match those in the inp file.
               if (npp_t /= npp .or. ncc_t /= ncc .or. nbac_t /= nbac) then
                  if (myRank == 0) write(*,'(a)') 'npp, ncc or nbac in traits file does not match that in params data file!'
                  ioerr = 1
               end if
               read(iunit_traits,*) epp
               !read(iunit_traits,*) f_pp
               !read(iunit_traits,*) phi_pp
               read(iunit_traits,*) tOn_pp
               read(iunit_traits,*) tOff_pp
               read(iunit_traits,*) omg1Amp_pp
               read(iunit_traits,*) ecc
               read(iunit_traits,*) omg_cc
               read(iunit_traits,*) ebac
               read(iunit_traits,*) omg_bac
            end if
            close(unit=iunit_traits)    
            ! see if any process messed up
            ioerr = abs(ioerr) ! remove negative sign, if present
            ! sum up ioerr across processes and place in ioE in ALL processes.  ioE should be 0 if no errors occured.
            call MPI_ALLREDUCE( ioerr, ioE, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)   
            ! check for error in opening the traits file by any of the processes.
            if (ioE /= 0) then
               ! an error occured reading the traits file by one or more processes
               if (myRank == 0) write(6,'(/a/)') 'Error opening or reading '//trim(basename)//'_traits.dat'//' file. Aborting.'
               call cleanUp()
               stop
            end if                  
         end if
            
         ! open output files. Note, only the root process saves output
         if (myRank /= 0) return
         open(newunit=iunit_state,   file=trim(basename)//'_state.dat'  ,status='unknown')
         open(newunit=iunit_rxn,     file=trim(basename)//'_rxn.dat'    ,status='unknown')
         open(newunit=iunit_sigDot,  file=trim(basename)//'_sigDot.dat' ,status='unknown')
         open(newunit=iunit_sigTot,  file=trim(basename)//'_sigTot.dat' ,status='unknown')
         open(newunit=iunit_epp,     file=trim(basename)//'_epp.dat'    ,status='unknown')
         open(newunit=iunit_omgpp,   file=trim(basename)//'_omgpp.dat'  ,status='unknown')
         open(newunit=iunit_ecc,     file=trim(basename)//'_ecc.dat'    ,status='unknown')
         open(newunit=iunit_omgcc,   file=trim(basename)//'_omgcc.dat'  ,status='unknown')
         open(newunit=iunit_ebac,    file=trim(basename)//'_ebac.dat'   ,status='unknown')
         open(newunit=iunit_omgbac,  file=trim(basename)//'_omgbac.dat' ,status='unknown')
         open(newunit=iunit_misc,    file=trim(basename)//'_misc.dat'   ,status='unknown') 
         open(newunit=iunit_traits,  file=trim(basename)//'_traits.dat' , status='unknown')   
         open(newunit=iunit_input,   file=trim(basename)//'_input.dat'  , status='unknown')   
         open(newunit=iunit_summary, file=trim(basename)//'_summary.dat', status='unknown')
         if (genSurf) open(newunit=iunit_EPsurf,  file=trim(basename)//'_EPsurf.dat' , status='unknown')   
         
         ! Control/trait variable techplot headers.  These are written with a call to storeTraits before ODE integration
         write(iunit_epp,  '(a)') 'Variables = "n" "epp"'
         !write(iunit_omgpp,'(a)') 'Variables = "n" "f_pp" "phi_pp"'
         write(iunit_omgpp,'(a)') 'Variables = "n" "tOn_pp" "tOff_pp" "omg1Amp_pp"'
         write(iunit_ecc,  '(a)') 'Variables = "n" "ecc"'
         write(iunit_omgcc,'(a)') 'Variables = "n"'
         do i=1, npp+ncc+nbac
            write(string,'(i)') i
            write(iunit_omgcc,'(a)') '"omg('//trim(adjustl(string))//',cc)"'
         end do         
         write(iunit_ebac,  '(a)') 'Variables = "n" "ebac"'
         write(iunit_omgbac,'(a)') 'Variables = "n" "omg(1,bac)" "omg(2,bac)" "omg(3,bac)"'  
         write(iunit_summary,'(a)') 'Variables = "Proc" "fcn calls" "tcpu (min)" "EP obj (J/K)" "sigTot_Rxn (J/K)" "sigTot_H2O (J/K)" & 
            &"sigTot_Part (J/K)" "sigTot (J/K)" "sigTot_Rxn (%)" "sigTot_H2O (%)" "sigTot_Part (%)"'
         if (genSurf) write(iunit_EPsurf,'(a)') 'Variables = "tOn_pp (d)" "tOff_pp (d)" "sigma_rxn (J/K)" "sigma_H2O (J/K)" "sigma_part (J/K)" "sigma_opt (J/K)"'
         return
      end subroutine initialize
   
      subroutine writeHeaders()
         ! This routine writes the headers to model output files
         ! local declarations
         integer i, j
         character (len=80) str80
          ! Write tecplot headers. Keep the p and c names on separate lines to avoid too many characters on a record.
         write(iunit_state,'(a)') 'Variables = "Date" "dic (<greek>m</greek>M)" "O2 (<greek>m</greek>M)" &
                                  &"NH3 (<greek>m</greek>M)"  "C_L (<greek>m</greek>M)"  "C_D (<greek>m</greek>M)" "N_D (<greek>m</greek>M)"'
         write(iunit_rxn,'(a)') 'Variables = "Date" '
         write(iunit_sigDot,'(a)') 'Variables = "Date" "sumSigDot-Rxn (J/d/K)" "sumSigDot-H2O (J/d/K)" "sumSigDot-Part (J/d/K)"&
                                 &"sumSigDot-Total (J/d/K)" "sumSigDot-Rxn (%)" "sumSigDot-H2O (%)" "sumSigDot-Part (%)"'
         write(iunit_sigTot,'(a)') 'Variables = "Date" "sumSigTot-Rxn (J/K)" "sumSigTot-H2O (J/K)" "sumSigTot-Part (J/K)"&
                                 &"sumSigTot-Total (J/K)" "sumSigTot-Rxn (%)" "sumSigTot-H2O (%)" "sumSigTot-Part (%)"'
         do i=1,npp
            write(str80,'(i)') i
            write(iunit_state,'(a)') '"p('//trim(adjustl(str80))//')"'
         end do
         do i=1,npp
            write(str80,'(i)') i
            write(iunit_state,'(a)') '"C_p('//trim(adjustl(str80))//')"'
         end do
         do i=1,ncc
            write(str80,'(i)') i
            write(iunit_state,'(a)') '"c('//trim(adjustl(str80))//')"'
         end do
         do i=1,nbac
            write(str80,'(i)') i
            write(iunit_state,'(a)') '"b('//trim(adjustl(str80))//')"'
         end do
         do i=1,npp
            do j=1,2
               write(str80,'(i0,a,i0)') i,',',j
               write(iunit_rxn,'(a)') '"rpp('//trim(adjustl(str80))//')"'
            end do
         end do
         do i=1,ncc
            do j=1,npp+ncc+nbac
               write(str80,'(i0,a,i0)') i,',',j
               write(iunit_rxn,'(a)') '"rcc('//trim(adjustl(str80))//')"'
            end do
         end do
         do i=1,nbac
            do j=1,3
               write(str80,'(i0,a,i0)') i,',',j
               write(iunit_rxn,'(a)') '"rbac('//trim(adjustl(str80))//')"'
            end do
         end do
         write(iunit_misc,'(a)') 'Variables = "Date" "I(h) (<greek>m</greek>moles/m2/s)"'
         do i=1,npp
            write(str80,'(i)') i
            write(iunit_misc,'(a)') '"omg_pp('//trim(adjustl(str80))//',1)" "omg_pp('//trim(adjustl(str80))//',2)"' 
         end do         
         return
      end subroutine writeHeaders
      
      subroutine writeZoneInfo(iter, valRank)
         ! This routine writes the tecplot zone info for a solution time or iteration
         integer , intent(in):: iter ! process that found the solution, and used for solutiontime
         real(mp), intent(in):: valRank(2) ! fnc value and process rank
         ! local declarations
         
         write(iunit_state ,'(a,i0,a,f0.2,a)') 'Zone T="State, proc: ',int(valRank(2)),',  fval = ',valRank(1),'"'
         write(iunit_state ,'(a,i0)') 'solutiontime = ', iter       
         write(iunit_rxn   ,'(a,i0,a,f0.2,a)') 'Zone T="Rates, proc: ',int(valRank(2)),',  fval = ',valRank(1),'"'
         write(iunit_rxn ,'(a,i0)') 'solutiontime = ', iter       
         write(iunit_sigDot,'(a,i0,a,f0.2,a)') 'Zone T="sigma dot, proc: ',int(valRank(2)),',  fval = ',valRank(1),'"'
         write(iunit_sigDot ,'(a,i0)') 'solutiontime = ', iter       
         write(iunit_sigTot,'(a,i0,a,f0.2,a)') 'Zone T="sigma total, proc: ',int(valRank(2)),',  fval = ',valRank(1),'"'  
         write(iunit_sigTot ,'(a,i0)') 'solutiontime = ', iter       
         write(iunit_misc,  '(a,i0,a,f0.2,a)') 'Zone T="misc variables, proc: ',int(valRank(2)),',  fval = ',valRank(1),'"'
         write(iunit_misc ,'(a,i0)') 'solutiontime = ', iter  
         ! The feed inputs are the same for all solutions, so rewind and just rewrite for all solutoins
         rewind(iunit_input)
         write(iunit_input, '(a)') 'Variables = "Date" "D (1/d)" "nh3I (<greek>m</greek>M)"'
         return
      end subroutine writeZoneInfo      
      
      subroutine cleanup ()
         ! This routine intializes and allocates variables
         integer mpierr
         deallocate( epp, ecc, ebac ) 
         deallocate( omg_pp, omg_cc, omg_bac )
         !deallocate( f_pp, phi_pp )
         deallocate( tOn_pp, tOff_pp, omg1Amp_pp )
         !$OMP parallel
         deallocate( pp, c_pp ) ! concentration of primary produces
         deallocate( pp_lg, c_pp_lg ) ! concentration of primary produces
         deallocate( cc ) ! concentration of consumers
         deallocate( cc_lg ) ! concentration of consumers
         deallocate( bac, bac_lg )
         deallocate( rpp, rcc, rbac )
         deallocate( n2_pp, aCi_cc ) 
         !$OMP end parallel
         call MPI_FILE_CLOSE(mpiFHerr, mpiErr)
         call MPI_FINALIZE(mpierr)    
         if (myRank /= 0) return
         close(unit=iunit_state)
         close(unit=iunit_rxn)
         close(unit=iunit_sigDot)
         close(unit=iunit_sigTot)
         close(unit=iunit_epp)
         close(unit=iunit_omgpp)
         close(unit=iunit_ecc)
         close(unit=iunit_omgcc)
         close(unit=iunit_ebac)
         close(unit=iunit_omgbac)
         close(unit=iunit_traits) 
         close(unit=iunit_input)
         close(unit=iunit_misc)
         close(unit=iunit_summary)
         if (genSurf) close(unit=iunit_EPsurf)
         return      
      end subroutine cleanup            
   end module globalVars

   module functions
      use realPrec
      implicit none
   contains
      real(mp) function fZero(xMin, x)
         ! This routine is used to prevent state variables from going negative or to zero, but in a continuous manner
         ! That is, while x >= 2 xMin, the function returns x, but as x decreases below 2 xMin, fZero also returns an number: xMin < fZero <= 2 xMin
         implicit none
         real(mp) xMin  ! Value that x must stay above or at
         real(mp) x     ! value of x that could be zero or less than zero that needs to be prevented from that

         if (x > 2._mp*xMin) then
            fZero = x
            return
         end if
         if (x < 0._mp) then
            fZero = xMin
            return
         end if
         ! x is between 0 and 2 xMin, so use quadratic for transition
         fZero = x**2/(4._mp*xMin) + xMin
         return
      end function fZero      
         
      real(mp) function F_Thermo (delG, ne, Tk)
         ! This function calculates the thermodynamic driver using LaRowe2012
         ! Function returns a unitless value between 0 and 1
         use globalVars, only: delPsi  ! The electric potential accross the cell wall (V).  Note, this is sigmodial 
                                       ! function, so F_Thermo is not zero if delG is less than it.  See LaRowe2012
         real(mp) delG  ! free energy of reaction (kJ/mol reaction)
         real(mp) ne    ! mole of electrons exchanged per mole of reaction extent
         real(mp) Tk    ! Temperature (K)
         ! Local declarations
         real(mp), parameter:: RkJ  = 8.3144598d-3   ! gas constant (kJ/(g-mol K) or J/(g-mmol K))
         real(mp), parameter:: F = 96485.3329 ! Faraday constant (C/mol-e)
         
         F_Thermo = 1.0_mp/( 1.0_mp + Exp( (delG/ne + F*delPsi/1000._mp)/(RkJ*Tk) ) )  ! Divide F*delPsi to get kJ/mol instead of J/mol
      end function  F_Thermo   
      
      real(mp) function I0(t)
         ! This fuction returns the light at the surface (mmol photons /m^2 /d (Not micromoles)) at time t (d)
         use globalVars, only: I0max, dLat
         real(mp) t
         real(mp) tJulian! since t is just linear time, it needs to be converted to a julian day, but don't bother with leap years.
         tJulian = mod(t,365._mp)
         call parSur (tJulian, dLat, I0max, I0)
! Note, this needs to be fixed for Tecplot which uses excel reference date  Need something that will convert excel date number to Ordinal Day
         return
      end function I0      
      
      real(mp) function stepUp(x, xm, sig)
         ! This routine generates a continuous sigmodial step up function from 0 to 1 around xm, with a gradient described by sig
         implicit none
         real(mp) x   ! value of x where function is evaluated
         real(mp) xm  ! value of x where setUp is 0.5
         real(mp) sig ! how steap the exponential function is around xm
         
         stepUp = 1._mp/( 1._mp + exp(-sig*(x - xm)) )
         return
      end function stepUp      
      
   end module functions   
     
   program AutoHetDet_Opt_Surf
      use mpi
      use omp_lib
      use realPrec
      use globalVars ! sets up space and stores parameters.  
      use hyperBOB_module
      use functions, only: fZero
      implicit none
      ! local variables
      integer i, ODEtimeOuts, ODEfailures
      type(solInfo) ODEinfo
      character(len=8) timeStr
      character(len=9) dateStr
      character(len=80) fmt80

      integer nU
      real(mp), allocatable:: traitSpace(:), tsL(:), tsU(:)
      real(mp) cpuTime
      ! MPI variables
      integer nameLen, noProc, mpierr, errBOB
      character (len=MPI_MAX_PROCESSOR_NAME) nodeName
      integer status(MPI_STATUS_SIZE)
      integer isent, sender,tag, irecv
      real(mp) WTime_0, WTime_f, WTime_t0
      ! openMP variables
      integer maxThreads, noThreads
      ! BOBYQA variables
      real(mp) maxTime ! This is not implemented in hyperBOB yet
      real(mp) EPbest ! entropy produced by best solution found
      external entropy
      ! hyperBOB related variables
      real(mp), allocatable:: EPsolns(:,:) ! this is the same as fVal_hyperBOB, but rank as second column, V4.4: +fncCalls
      integer:: mEP = 3 ! this is the column dimention of EPsolns

      ! Initialize MPI
      call MPI_INIT( mpierr )
      call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpierr) ! get rank of this process in world      
      if (myRank == 0) write(*,'(/a/)') 'AutHetDet_Opt Ver 4.7, 16-Jul-2020'
      ! initialize space via sysSetup
      call sysSetup ()
      ! Get info on threads
      maxThreads = OMP_GET_MAX_THREADS()
      call omp_set_num_threads(min(maxThreads, ompThreads)) ! Don't use more threads than available
      !$OMP parallel
      noThreads = OMP_GET_NUM_THREADS() ! Outside of omp parallel regions, this returns 1
      !$OMP end parallel
      call initialize () ! allocate variables. 
      tf = t0 + tDays ! set the end of the simulation 
      ! Get the number of processes running and the maxium number of threads
      ! Currently not using OpenMP
      call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)
      if (noProc == 1 .and. genSurf) then
         write(*,'(a)') 'Must run EP surface routine with more than 1 MPI process! Stopping.'
         stop
      end if
      if (myRank == 0) write(*,'(/,a,i4,a)') 'Running program with ',noProc, ' processes'     
      call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      call MPI_GET_PROCESSOR_NAME(nodeName, nameLen, mpiErr)
      do i=0,noProc-1
         if (i == myRank) then
            write(*,'(a,i4,2a,2(a,i4))') ' Process ', myRank, ' on node "', trim(nodeName), '" has ', noThreads, ' threads out of: ', maxThreads
            !write(*,'(a,i12,/)')  ' The thread stack size is set to: ', KMP_GET_STACKSIZE_S()   
         end if
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      end do  
      call MPI_Barrier(MPI_COMM_WORLD, mpiErr) 
      maxODEtime = (tf-t0)/minCompFac ! maximum allowed CPU wall time allowed to solve the ODE over t0 to tf.
      if (myRank == 0) then 
         write(*,'(//a)') 'Parameters being used for simulation:'
         write(*,nml=params)
         write(*,'(/,a,f0.2,a)') 'Maximum ODE solution time set to:  ', maxODEtime*24.0*60., ' (min)'
         write(*,'(a,f0.2,a,/)')   'Maximum possible optmization time: ', maxODEtime*real(maxFun), ' (days)'
      end if
      
      ! In V4.5, setting genSurf to true generates the 2D EP surfaces only
      if (genSurf) then
         if (readeps) then
            call EPsurfGen()
         else
            if (myRank == 0) write(*,'(a)') 'Parameter readeps must be set to .true. if 2D surface is to be created. Stopping.'
         end if         
         call cleanup ()      
         stop
      end if
      
      ! traitSpace stores epp, w_pp, ecc, omg_cc, ebac and w_bac.  Note, the rows of omg_pp and omg_bac must sum
      ! to one, so there is one DOF less than the column dimension of omg_pp and omg_bac (DOFs are 1 and 2, repsectively)
      ! NOTE, with V3.0, omg_pp is calculated from f_pp and phi_pp.  
      ! in version 4 omg_pp is calculated from tOn_pp, tOff_pp and omg1Amp_pp
      ! traitSpace storage: [epp(1:npp), w_pp(1:npp), ecc(1:ncc), omg_cc(1:ncc,1)...omg_cc(1:ncc,npp+ncc+nbac), ebac(1:nbac), w_bac(1:1:2) ... w_bac(nbac,1:2)]
      lenTS = 4*npp &                        ! This is for epp, tOn_pp, tOff_pp, omg1Amp_pp
            + ncc + ncc*(npp + ncc + nbac) & ! This is for ecc and omg_cc (note, no constraint on row sums due to biomass weighting (see Gz_bioS)
            + 3*nbac                         ! This is for ebac and w_bac(1:nbac,1:2), stored row-wise (rows must sum to 1) 
      allocate( traitSpace(lenTS), tsL(lenTS), tsU(lenTS), EPsolns(noProc, mEP) ) 
      call hyperBOB_initialize(lenTS) ! setup space for hyperBOB      
      ! ******************************************************************************************************      
      ! *** Begin optimization
      ! ******************************************************************************************************
      nU = lenTS
      ! setup npt_bobyqa which is Number of interpolation points between [N+2,(N+1)(N+2)/2]. best if < 2*N+1
      npt_bobyqa = 2*nU + 1  
      ! Initial conditions for bobyqa and set eps and omg vector/matrices; however, hyperBOB does not use the
      ! initial values for traitSpace assigned here, as it will generate it's own, but tsL and tsU are needed.
      ! Note, initializeTraits is not needed here, as all bounnds are not between 0 and 1.  
      ! call initializeTraits(traitSpace, tsL, tsU) 
      tsL = 0._mp
      tsU = 1._mp
      WTime_it0 = MPI_Wtime() ! Get the clock start time
      fcnCalls = 0
      ODEfailed = 0; ODEtimedOut = 0
      if (optimize) then
         ! run optimization
         call time(timeStr); call date(dateStr)         
         if (myRank == 0) write(*,'(a)') 'Begining hyperBOB optimization at '//timeStr//', '//dateStr//' ...'
         call hyperBOB(nU, npt_bobyqa, traitSpace, tsL, tsU, rhoBeg, rhoEnd, &
            iprint, maxFun, maxTime, entropy, iseed, EPbest, errBOB) 
         WTime_f = MPI_Wtime() ! get the clock time
         cpuTime = real(WTime_f-WTime_it0)/60.0/60. ! hrs
         ! copy fVal_hyperBOB to EPsolns, with rank
         EPsolns(1:noProc,1) = -fVal_hyperBOB(1:noProc) ! Take negative.  All ranks have the full fVal_hyperBOB vector 
         EPSolns(1:noProc,2) = [(real(i,mp),i=0,noProc-1)] ! this adds the rank corresponding to each fVal.
         EPSolns(myRank+1,3) = fcnCalls ! this column needs to be passed to all processes
         call MPI_allgather(MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, EPSolns(1:noProc,3), 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierr)
         ! sort EPsolutions (all processes do this)
         call dsortArray(noProc, mEP, EPsolns, 1) ! EP solutions now has the ascending order of best solution
         call MPI_ALLREDUCE(fcnCalls, fcnCallsTot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr) ! sum fcnCalls for all processes.  
         call time(timeStr); call date(dateStr)    
         if (myRank == 0) then
            write(*,'(/,a,f0.2,a)') 'hyperBOB optimization finished at at '//timeStr//', '//dateStr//' taking ', cpuTime, ' hrs'
            write(*,'(a,i0)') '   Total number of ODE integrations: ', fcnCallsTot
            write(*,'(a,i0,a,f0.0)') '   Process ', int(EPsolns(noProc,2)), ' found best solution: ', EPsolns(noProc,1)
            write(*,'(a)') ' '
            ! Store all the traitSpace solutions found in ascending order
            do i=1,noProc 
               call setTraits(xMat_hyperBOB(1:lenTS,int(EPsolns(i,2))+1))
               call storeTraits(i, EPsolns(i,1))
            end do
         end if
         call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
         do i=0,noProc-1
            if (i == myRank) then
               write(*,'(a,i0,a,f0.0,a,i0,a)') '   Process ', myRank, ' bestEP = ', -fVal_hyperBOB(myRank+1), ' after ', fcnCalls, ' fcn calls'
            end if
            call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
         end do  
      else 
         ! In this case, just randomly generate noProc guesses of the traits
         if (myRank == 0) then
            write(*,'(/,a)') 'Random solutions being generated. NO optimization done!'
            call RANDOM_NUMBER(xMat_hyperBOB) ! randomly generate all solutions
            ! Since the entropy production has not been calculated, just use a value of 0 for ep
            ! Note, the code could be expanded so that each process calculates EP for its random IC, but not really worth it.
            do i=1,noProc 
               EPsolns(i,1) = 0._mp
               EPsolns(i,2) = real(i-1)
               EPsolns(i,3) = 1._mp               
               call setTraits(xMat_hyperBOB(1:lenTS,i))
               call storeTraits(i, 0._mp)
            end do
         end if
      end if
      
      call MPI_Barrier(MPI_COMM_WORLD, mpiErr)
      if (myRank == 0) then
         write(*,'(/,a)') 'Storing all ODE solutions...'
         call writeHeaders() ! write tecplot headers.  
         ! Write the time hyperBOB took at a zone title in the summary file
         write(iunit_summary,'(a,f0.3,a)') 'Zone T = "Optimization time = ', cpuTime, ' (hr)"'
         do i=1,noProc
            ! Set traitSpace to the best one returned by hyperBOB for ODE solution next
            call setTraits( xMat_hyperBOB(1:lenTS,int(EPsolns(i,2))+1) ) ! Copies traitSpace to eps and omg vector/matrices for best solution            
            call writeZoneInfo(i, EPsolns(i,1:2))
            WTime_t0 = MPI_Wtime() ! Get the clock start time   
            ! set the state variables to their intitial conditions
            dic = dicI; o2 = o2I; nh3 = nh3I_t0; c_l = c_lI; c_d = c_dI; n_d = n_dI; pp = ppI; c_pp = c_ppI; bac = bacI; cc = ccI; sumSigTot = 0._mp
            call integrateODE(t0, tf, .true., ODEinfo)    
            WTime_f = MPI_Wtime() ! get the clock time
            fmt80 = ''
            if (ODEinfo%fail /= 0) fmt80 = ' ** Warning ODE integration failure!'
            cpuTime = real(WTime_f-WTime_t0)/60.0 ! min
            write(*,'(a,i0,a,f0.2,a)') 'Process ', int(EPsolns(i,2)), ' ODE solution time: ', cpuTime, ' (min)'//trim(fmt80) 
            ! write the summary out
            write(iunit_summary,'(2(f0.0,1x),f0.2,1x,8(f0.0,1x))') EPSolns(i,2), EPSolns(i,3), cpuTime, EPSolns(i,1), sumSigTot(1:3), sum(sumSigTot(1:3)), &
               100_mp*sumSigTot(1:3)/fZero(absZero,sum(sumSigTot(1:3)))
         end do         
      end if      
      ! clean things up.
      deallocate( traitSpace, tsL, tsU, EPsolns ) 
      call hyperBOB_cleanUp()      
      call cleanup ()      
   end program AutoHetDet_Opt_Surf
   
   subroutine EPsurfGen()
      use mpi
      use realPrec
      use globalVars ! sets up space and stores parameters.  
      use functions, only: fZero
      implicit none
      ! local variables
      integer i
      type(solInfo) ODEinfo
      character(len=8) timeStr
      character(len=9) dateStr
      character(len=80) fmt80
      ! EP surface declarations
      integer totalPts ! total number of function evaulations to produce 2D surface
      real(mp), allocatable:: epSurf(:,:,:) ! Need to store the EP surfaces, but on the master needs this.  Four surfaces are saved.
      real(mp) xGrid, yGrid
      integer ithPt ! keep track of what point in the 2D surface is being calculated.
      integer epI, epJ
      integer, parameter:: iresults = 6 ! this is the size of results vector
      real(mp) results(iresults)
      real(mp) EPsolns(2)
      real(mp) eProd(3)
      ! mpi declarations
      real(mp) WTime_0, WTime_f, WTime_t0, cpuTime
      integer noProc, mpierr
      integer status(MPI_STATUS_SIZE)
      integer isent, sender,tag, irecv

      ! get the number of processes
      call MPI_COMM_SIZE(MPI_COMM_WORLD, noProc, mpiErr)
      
      ! ******************************************************************************************************      
      ! *** Begin surface generation using MPI master/slave
      ! ******************************************************************************************************
      totalPts = nSurfPts*nSurfPts ! total number of function evaluations to get 2D surface
      WTime_it0 = MPI_Wtime() ! Get the global clock start time      
      WTime_0   = MPI_Wtime() ! Get the start time for output dispaly     
      if (myRank == 0) then
         ! The master process doling out slaves
         allocate( epSurf(nSurfPts, nSurfPts, 4) ) ! only the master needs this
         ! First, send jobs to all the processes except for the master
         call time(timeStr); call date(dateStr)         
         write(*,'(a,i0,a,i0,a)') 'Starting 2D sigma_rxn surface generation using ', nSurfPts,' x ',nSurfPts, ' grid at '//timeStr//' on '//dateStr         
         isent = 0
         do i=1, min(noProc-1,totalPts) ! The master process does not do ODE integrations
            tag = i
            call MPI_SEND (tag, 1, MPI_INTEGER, i, tag, MPI_COMM_WORLD, mpierr)   
            isent = isent + 1
         end do
         ! Now, start receiving jobs back, but also send out new ones after they finished if totalPts has not been reached.
         do i=1, totalPts
            ! results contains: 1) ithPt, 2) sigma_rxn, 3) error ID
            call MPI_RECV (results, iresults, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, &
               mpierr)
            sender = status(MPI_SOURCE) ! process number that has returned a solution
            ! store the sigma_rxn solution in the epSurf matrix
            ithPt = int(results(1))
            call ijPoint(nSurfPts, ithPt, epI, epJ, xGrid, yGrid)
            ! if an error occured, display that
            if (int(results(6)) /= 0) then
               write(*,'(5(a,i0),a)') '***Process: ', sender, ' returned BiM error: ', int(results(6)), ' (i=', epI,', j=', epJ,', ithPt=', ithPt, ')'
               results(5) = -1.0_mp ! set sigma_rxn to -1 so it can be easily observed in the plotting
            end if
            if (results(5) > bestEPfound) then
               bestEPfound = results(5)
               bestI = epI
               bestJ = epJ
            end if    
            epSurf(epI, epJ,1:4) = results(2:5) ! Note, if there was an error, results(5) = -1.0
            ! determine how much time has passed in min since last reset
            WTime_f = MPI_Wtime() ! get the clock time
            cpuTime = real(WTime_f-WTime_0)/60.0_mp ! min
            if (cpuTime >= reportTime) then
               cpuTime = real(WTime_f-WTime_it0)/60.0_mp/60._mp ! hrs
               call time(timeStr); call date(dateStr)         
               write(*,'(a,f0.2,a,i0,a,i0,a,f0.1,a)') '   After ', cpuTime, ' hrs, ', i, ' surface points have been determined (', &
                  int(100.0*real(i)/real(totalPts)),'%) at '//timeStr//' on '//dateStr//'. Time left: ', cpuTime*(real(totalPts)/real(i) - 1.0_mp),' hrs.'
               ! reset the report clock
               WTime_0 = MPI_Wtime()
            end if                           

            ! Since the sender just finished, send next job to it, unless totalPts has been found.          
            if (isent < totalPts) then
               ! keep going and send another job
               tag = isent + 1
               call MPI_SEND (tag, 1, MPI_INTEGER, sender, tag, MPI_COMM_WORLD, mpierr)   
               isent = isent + 1           
            else
               ! tell sender it's done.
               call MPI_SEND (MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, 0, MPI_COMM_WORLD, mpierr)
            end if            
         end do
      else
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! slave code
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         do      
            if (myRank > totalPts) exit ! more processes than needed
            call MPI_RECV (ithPt, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpierr)
            if (status(MPI_TAG) == 0) exit ! free the slaves        
            ! set the values of f_pp and phi_pp based on ithPt
            call ijPoint(nSurfPts, ithPt, epI, epJ, xGrid, yGrid)
            tOn_pp(whichPP)  = tOn_pp_min  + xGrid*(tOn_pp_max  - tOn_pp_min )
            tOff_pp(whichPP) = tOff_pp_min + yGrid*(tOff_pp_max - tOff_pp_min)
            ! Integrate ODE.
            dic = dicI; o2 = o2I; nh3 = nh3I_t0; c_l = c_lI; c_d = c_dI; n_d = n_dI; pp = ppI; c_pp = c_ppI; bac = bacI; cc = ccI; sumSigTot = 0._mp
            ODE: block ! this is being used to exit on error in any inteval integration
               ! First interval integration from t0 to t0_ep (if t0_ep = t0, then integrateODE just returns)
               call integrateODE(t0, t0_ep, .false., ODEinfo)
               if (ODEinfo%fail /= 0) exit ODE ! exit the block, as there is no need to continue integration
               eProd(1:3) = sumSigTot(1:3) ! Store cumulative entropy production up to t0_ep
               ! integrate the next interval over which EP is determined.
               call integrateODE(t0_ep, tf_ep, .false., ODEinfo)
               if (ODEinfo%fail /= 0) exit ODE ! exit the block, as there is no need to continue integration
               eProd(1:3) = sumSigTot(1:3) - eProd(1:3) ! This is the total entropy produced over t0_ep to tf_ep
               ! There is no need to inegrate to tf because it does not change entropy prod over the specified EP interval.
            end block ODE
            ! send results back to master process
            ! results contains: 1) ithPt, 2) sigma_rxn, 3) error ID
            results(1) = real(ithPt)
            results(2) = eProd(1) ! This is entropy production from reactions only
            results(3) = eProd(2) ! This is entropy production from H2O only
            results(4) = eProd(3) ! This is entropy production from particles only
            results(5) = dot_product(sumSigWeights,eProd) ! sumSigWeights weights each of the three terms: Rxn, H2O, and Particle
            results(6) = real(ODEinfo%fail)
            tag = status(MPI_TAG) ! note, this is the same as i
            call MPI_SEND (results, iresults, MPI_DOUBLE_PRECISION, 0, tag, MPI_COMM_WORLD, mpierr) 
         end do
      end if         

      if (myRank == 0) then
         WTime_f = MPI_Wtime() ! get the clock time
         cpuTime = real(WTime_f - WTime_it0)/60.0_mp/60._mp ! hrs  
         ithPt = (bestI-1)*nSurfPts + bestJ
         call ijPoint(nSurfPts, ithPt, epI, epJ, xGrid, yGrid)
         tOn_pp(whichPP)  = tOn_pp_min  + xGrid*(tOn_pp_max  - tOn_pp_min )
         tOff_pp(whichPP) = tOff_pp_min + yGrid*(tOff_pp_max - tOff_pp_min)
         write(*,'(/,a)')             'Finished generating 2D sigma_rxn surface'
         write(*,'(a,f0.1,a)')        '   Total time: ', cpuTime, ' hrs'
         write(*,'(a,f0.0,2(a,i0))')  '   Best EP solution: ', bestEPfound, ' (J/K) at i = ', epI, ', j = ', epJ
         write(*,'(2(a,f0.6))')       '   tOn_pp = ', tOn_pp(whichPP),', tOff_pp = ', tOff_pp(whichPP) 
         write(*,'(/,a)')             'Saving best solution and EP surfaces...'
         call writeHeaders() ! write tecplot headers.  
         call storeTraits(1, bestEPfound)
         EPsolns = [bestEPfound, 0._mp] 
         call writeZoneInfo(1, EPsolns)
         WTime_t0 = MPI_Wtime() ! Get the clock start time   
         dic = dicI; o2 = o2I; nh3 = nh3I_t0; c_l = c_lI; c_d = c_dI; n_d = n_dI; pp = ppI; c_pp = c_ppI; bac = bacI; cc = ccI; sumSigTot = 0._mp
         call integrateODE(t0, tf, .true., ODEinfo)         
         WTime_f = MPI_Wtime() ! get the clock time
         fmt80 = ''
         if (ODEinfo%fail /= 0) fmt80 = ' ** Warning ODE integration failure!'
         cpuTime = real(WTime_f-WTime_t0)/60.0 ! min
         write(*,'(a,i0,a,f0.2,a)') 'Process ', 0, ' ODE solution time: ', cpuTime, ' (min)'//trim(fmt80) 
         ! write the summary out
         write(iunit_summary,'(2(f0.0,1x),f0.2,1x,8(f0.0,1x))') 0., 1., cpuTime, bestEPfound, sumSigTot(1:3), sum(sumSigTot(1:3)), &
            100_mp*sumSigTot(1:3)/fZero(absZero,sum(sumSigTot(1:3)))   
         call writeEPsurf(whichPP, epSurf) ! save the EP surfaces. Note, this changes the values for tOn and tOff, so doing last here.         
      end if      
      ! clean things up.
      if (myRank == 0) deallocate( epSurf ) 
      return
   end subroutine EPsurfGen
   
   subroutine entropy (n, x, f)
      ! this routine is called by BOBYQA to find the function value, in this case entropy production over simulation time t0_ep to tf_ep
      ! Note, t0 <= t0_ep < tf and t0 < tf_ep <= tf.  Previous version to V4.4 implicitly set t0_ep = t0 and tf_ep = tf.
      ! This version does potentially tree enegrations
      use realPrec
      use globalVars
      use mpi
      use kind_module
      implicit none
      ! Dummy variables
      ! ****** NOTE, the function must be defined EXACTLY as shown by the interface in bobyqa.  *******
      ! that is, don't instead declare as "real(wp) x(n)", etc as that will cause problems      
      integer,intent(in)               :: n ! number of control variables
      real(wp),dimension(:),intent(in) :: x ! control variables
      real(wp),intent(out)             :: f ! objective function, integrated entropy production
      ! Local declarations
      integer i, j    
      type(solInfo) ODEinfo ! info on ODE success, etc.
      real(mp) WTime_f, cpuTime, WTime_s, intTime, eProd(3)
      character(len=8) timeStr
      character(len=9) dateStr
      
      fcnCalls = fcnCalls + 1
      call setTraits(x)  ! set eps and omg, etc  
      WTime_s = MPI_Wtime() ! get the clock time
      ! set the state variables to their intitial conditions
      dic = dicI; o2 = o2I; nh3 = nh3I_t0; c_l = c_lI; c_d = c_dI; n_d = n_dI; pp = ppI; c_pp = c_ppI; bac = bacI; cc = ccI; sumSigTot = 0._mp
      ODE: block ! this is being used to exit on error in any inteval integration
         ! First interval integration from t0 to t0_ep (if t0_ep = t0, then integrateODE just returns)
         call integrateODE(t0, t0_ep, .false., ODEinfo)
         if (ODEinfo%fail /= 0) exit ODE ! exit the block, as there is no need to continue integration
         eProd(1:3) = sumSigTot(1:3) ! Store cumulative entropy production up to t0_ep
         ! integrate the next interval over which EP is determined.
         call integrateODE(t0_ep, tf_ep, .false., ODEinfo)
         if (ODEinfo%fail /= 0) exit ODE ! exit the block, as there is no need to continue integration
         eProd(1:3) = sumSigTot(1:3) - eProd(1:3) ! This is the total entropy produced over t0_ep to tf_ep
         ! There is no need to inegrate to tf because it does not change entropy prod over the specified EP interval.
         ! However, it is possible that the integration fails within tf_ep and tf, but for now, assume that is rare.
         !call integrateODE(tf_ep, tf, .false., ODEinfo)
      end block ODE
      WTime_f = MPI_Wtime() ! get the clock time
      f = -dot_product(sumSigWeights,eProd) ! sumSigWeights weights each of the three terms: Rxn, H2O, and Particle
      if (ODEinfo%fail /= 0) then
         f = 0._mp ! BOBYQA does not handle failures, so just set to 0.
         if (ODEinfo%fail == 1) ODEfailed   = ODEfailed   + 1
         if (ODEinfo%fail == 2) ODEtimedOut = ODEtimedOut + 1
      end if
      if (-f > fcnMax) fcnMax = -f
      cpuTime = real(WTime_f-WTime_it0)/60.0/60. ! hrs
      intTime = real(WTime_f-WTime_s)/60.0 ! min
      aveODEtime = aveODEtime + intTime ! determine the average PDE integration time over fcnUpdate steps
      if (mod(fcnCalls, fcnUpdate) == 0._mp) then
         call time(timeStr); call date(dateStr)         
         write(*,'(a,i0,a,f0.1,3(a,i0),a,f0.0,a,f0.2,a)')  'Process ', myRank, ' after ', cpuTime, ' hrs, ODE integtrated ', &
            fcnCalls, ' times (',ODEfailed, ' fail, ', ODEtimedOut, ' TO), fcnMax = ', fcnMax, &
            ' (ODE time: ', aveODEtime/real(fcnUpdate),' min; '//timeStr//', '//dateStr//')'
         aveODEtime = 0._mp         
      end if      
      return
   end subroutine entropy  
   
   subroutine integrateODE(t0_ode, tf_ode, strSoln, ODEinfo)
      ! This routine integrates the ADR equation over t0 to tf
      use mpi
      use realPrec
      use globalVars
      implicit none     
      real(mp),      intent(in) :: t0_ode ! start of integration time (d)
      real(mp),      intent(in) :: tf_ode ! end of integration time (d)
      logical,       intent(in) :: strSoln ! if set to true, the solution is stored.
      type(solInfo), intent(out):: ODEinfo 
      ! This derived type has three variables
      ! If the integration fails or times out, ifail is set to non zero value
      !  ifail = 0 ODE integrated without error and within the max time allowed
      !  ifail = 1 ODE generated an error
      !  ifail = 2 ODE took longer than it's maximum time allowed.
      !  idid: value of idid from BACOLI
      !  t0: value of t0 (where the solution got to)
      
      ! local declarations
      integer i, j, mpierr
      real(mp) ODE_ti ! current wall time at iteration
      logical repeat
      integer attempts, iout
      character(len=200) string
      character(len=8) timeStr
      character(len=9) dateStr

      ! Local declarations needed for BiM
      real(mp) x(nStates), tEnd          
      integer lenw, leniw
      integer  ierflg, ijac, mljac, mujac, idid, ierr, ipar(1)
      real(mp) rtolL, atolL, h_BiM, rpar(1), t0L
      real(mp), allocatable:: wk(:)
      integer, allocatable:: iwk(:)
      external ODEs, ompJAC
      
      ODEinfo%fail = 0 
      ODEinfo%t0 = t0_ode
      ODEinfo%idid = 0
      if(tf_ode <= t0_ode) return ! nothing to integrate.
      
      ! Space needed by BIM
      lenw = 14+10+8*nstates+4*10*nstates+2*nstates**2
      leniw = nstates+37
      allocate (  wk(lenw) )
      allocate ( iwk(leniw) )
      
            
      repeat = .True.
      attempts = 1
      rtolL = rtol1 ! relative and absolute tolerances in parameter module
      atolL = atol1
      
      !************* Begin ODE Solution from t0_ode to tf_ode ***********************
      ODE_t0 = MPI_Wtime() ! start clock for this solution
      t0L = t0_ode ! on output from BiM, t0_ode is where the solution got to, so make this local variable
      ! Initialize x at t0_ode. Note, these global variables must have been set before calling integrateODE   
      !x(1:nstates) = [dicI, o2I, nh3I_t0, c_lI, c_dI, n_dI, (ppI, i=1,npp), (c_ppI, i=1,npp), (ccI, i=1,ncc), (bacI, i=1,nbac), 0._mp,  0._mp,  0._mp]      
      x(1:nstates) = [dic, o2, nh3, c_l, c_d, n_d, pp(1:npp), c_pp(1:npp), cc(1:ncc), bac(1:nbac), sumSigTot(1:3)]      
      Do While (repeat)
         repeat = .False. 
         h_BiM = 0. ! use default initial BiM stepsize
         ijac = useOmpJac ! BiM: 0=calculate jacobian numerically, otherwise use analytical routine provided.
         iout = 0
         if (strSoln) iout = 1 ! BiM: set to 1 to call solout.  BiM determines how frequently to save solution points
         mljac = nstates; mujac = nstates !BiM: banding aspect of jacobian? just set to nstates.
         iwk(1:8) = 0
         wk(1:14) = 0.0d0
         iwk(1) = maxstep_BiM ! default iterations 100000 used if set to zero
         wk(1) = epsilon(wk) ! set true precision
         ! note, after email exchanges with Cecilia Magherini, hmax actual refers to internal mesh spacing, that depends on the order of the method
         ! At maximum orde (12), there can be 10 interval, so you need to divide hmax by 10 if an output point needs to be < hmax
         wk(2) = hmax_BiM/10._mp ! largest integration step size (use to handle high frequency drivers) DEFAULT = (TEND-T0)/8
         idid = 0
         CALL BiM (nstates, ODEs, t0L, tf_ode, x, h_BiM, rtolL, atolL, ompJAC, ijac, mljac, mujac, &
                     wk, lenw, iwk, leniw, rpar, ipar, iout, idid) 
         ODE_ti = MPI_Wtime() ! get the clock time
         ODEinfo%t0 = t0L
         ODEinfo%idid = idid
         call mapx2names(x,nstates) ! this is mostly used to set sumSigTot(1:3) at t0L (=tf_ode if no failures)
         ! First check if maximum ODE time has been exceeded.
         if ((ODE_ti-ODE_t0)/86400._mp > maxODEtime) then
            ! Maximum solution time exceeded. Note, BiM does not handle ierr in the ODE subrouinte correctly, so it can return idid = 0 when it should have been -6.
            ODEinfo%fail = 2
            call time(timeStr); call date(dateStr)         
            write(string,'(a,i0,a,f0.2,a)') 'Process: ', myRank, ' ODE TIME EXCEEDED (t0 = ',  t0L, ') ('//timeStr//', '//dateStr//')'    
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)            
            return ! exits this ODE solution loop
         end if
         ! Check for other integration errors
         If (idid < 0) Then             
            ! Save the error info
            call time(timeStr); call date(dateStr)         
            write(string, '(a,i4,a,f0.2,2(a,i0),a)') 'Error in BiM::idid = ', idid, ' at time: ', t0L, ' Attempts: ', attempts, &
               ' for Process: ', myRank,' ('//timeStr//', '//dateStr//')'
            call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)            
            rtolL = 10.*rtolL
            atolL = 10.*atolL
            repeat = .True.
            if (attempts > maxattempts) then
               write(string, '(a,i4,a,f0.2,a,i0,a)') '**MAX** Attempts in BiM::idid = ', idid, ' at time: ', t0L, ' for Process: ', &
                  myRank,' ('//timeStr//', '//dateStr//')'
               call MPI_FILE_WRITE_SHARED(mpiFHerr, trim(string)//CRLF, len_trim(string)+2, MPI_CHARACTER, MPI_STATUS_IGNORE, mpiErr)               
               ODEinfo%fail = 1
               return
            end if
            attempts = attempts + 1
         End If
      End Do            
      deallocate (  wk )
      deallocate ( iwk )
      return
   end subroutine integrateODE

   subroutine solout(m, t, x, f, k, ord, irtrn, rpar, ipar) 
      ! This routine is called by BiM.  It turns out that BiM runs much faster
      ! if you allow it to determine when to output points.  It will store points
      ! often when the state is changing fast, and slowley when not.  Makes for the
      ! smoothest output with the least amount of points (and it's faster)
      ! NOTE must use bim_SOLOUTmod.f with this SOLOUT routine.
      ! This routine is just a wrapper for saveSolution
      use realPrec
      use globalVars
      implicit none
      integer , intent(in):: m
      integer , intent(in):: k
      integer , intent(in):: ord
      integer , intent(in):: irtrn
      real(mp), intent(in):: t(k)
      real(mp), intent(in):: x(m,k)
      real(mp), intent(in):: f(m,k)
      integer , intent(in):: ipar(1)
      real(mp), intent(in):: rpar(1)
      ! Local declarations
      
      if (myRank == 0) call storeSoln (t(k), x(1:nconc,k))            
      return
   end subroutine solout  
   
   subroutine sysSetup ()
      use globalVars
      use mpi
      implicit none
      ! Local declarations
      integer iunit, ioerr, ioE, mpierr
      integer narg

      ! First try getting basename from commandline
      ! if that fails, then have processes zero get name and broadcast it to other processes 
      narg = command_argument_count () ! see if a filename has been included in command line
      if (narg == 1) then
         call get_command_argument(narg,value=basename)
      else 
         if (myRank == 0) then
            write(6,'(a,$)') 'Enter parameter filename base, no extension: '
            read(5,*) basename
         end if
         ! send basename to all processes
         call MPI_BCAST(basename,80,MPI_CHARACTER,0,MPI_COMM_WORLD,mpierr)
      end if
      ! open a readonly unit to the inp file. This has to be readonly so that each process can open
      ! it without conflict. This file must be available for all processes.
      open(newunit=iUnit,file=trim(basename)//'.inp',action='read',status='old',iostat=ioerr)   
      if (ioerr /= 0) then 
         write(6,'(a,i5,a)') 'Process ',myRank,' could not open INP file' ! this may or maynot print.
      else   
         read(iunit,nml=params) ! read parameters. see module globalVars for items in list inputData
         close(unit=iunit)   
      end if
      ! with V3.7, set c_ppI to be consistent with k_pp if k_pp is not equal to 0.  
      ! If k_pp equals 0 (or is very small), then use C_ppI from parameter file.
      if (k_pp >= epsilon(k_pp)) c_ppI = k_pp*ppI
      ! see if any process messed up
      ioerr = abs(ioerr) ! remove negative sign, if present
      ! sum up ioerr across processes and place in ioE in ALL processes.  ioE should be 0 if no errors occured.
      call MPI_ALLREDUCE( ioerr, ioE, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)   
      ! check for error in opening the PRM file by any of the processes.
      if (ioE /= 0) then
         ! an error occured reading the prm file by one or more processes
         if (myRank == 0) then
            write(6,'(/a/)') 'Error opening or reading '//trim(basename)//'.inp file. Aborting.'
         end if
         call cleanUp()
         stop
      end if         
      return
   end subroutine sysSetup

   subroutine mapx2names(x,n)
      ! This routine copies the x vector into meaningful names 
      ! pointer might be better here for pp, c_pp and cc, or simply elsewhere
      use realPrec
      use globalVars
      use functions
      implicit none
      integer , intent(in):: n
      real(mp), intent(in):: x(n)
      ! local declarations
      integer i
      real(mp), parameter:: ln106 = -6._mp*log(10._mp) ! this is the natural logorithm of 10^-6 to convert uM to M

      dic   = fZero(absZero, x(1)); dic_lg   = log(dic)  + ln106
      o2    = fZero(absZero, x(2)); o2_lg    = log(o2)   + ln106
      nh3   = fZero(absZero, x(3)); nh3_lg   = log(nh3)  + ln106
      c_l   = fZero(absZero, x(4)); c_l_lg   = log(c_l)  + ln106
      c_d   = fZero(absZero, x(5)); c_d_lg   = log(c_d)  + ln106
      n_d   = fZero(absZero, x(6)); n_d_lg   = log(n_d)  + ln106
      
      ! Since h3po4 and P_D are held constant, these only need to be calculated once, but put here anyway
      h3po4_lg = log(h3po4) + ln106
      p_d_lg   = log(P_D)   + ln106

      do i=1, npp
         pp(i)   = fZero(absZero, x(nconc+i));       pp_lg(i)   = log(pp(i)) + ln106
      end do
      do i=1, npp
         c_pp(i) = fZero(absZero, x(nconc+npp+i));   c_pp_lg(i) = log(c_pp(i)) + ln106
      end do      
      do i=1, ncc
         cc(i)   = fZero(absZero, x(nconc+2*npp+i)); cc_lg(i)   = log(cc(i)) + ln106
      end do
      do i=1, nbac
         bac(i)  = fZero(absZero, x(nconc+2*npp+ncc+i)); bac_lg(i)   = log(bac(i)) + ln106
      end do
      sumSigTot(1:3) = x(n-2:n) ! integrated entropy production for reactions, water and particles to current time (J/K)
      return
   end subroutine mapx2names

   subroutine storeSoln (t, x)
      ! This routine stores the output.  It is assumed that the file has already been opened and
      ! that the solution has been mapped to variables names.
      use realPrec
      use globalVars
      use functions
      implicit none
      real(mp)          , intent(in) :: t ! Time solution is valid for
      real(mp)          , intent(in) :: x(nstates) ! solution from BACOLI at time t
      ! local declarations
      integer i,j,k
      if (myRank /= 0) return ! only root process saves data
      ! make sure everything is up to date.  This routine only assumes that all state variables have been defined.
      ! everything else is calculated here.
      ! These are needed by by rxnRates
      call mapx2names(x, nstates)
      call rxnRates(t) ! Calculate reaction rates and sumSigDot for rxns, h2o and paricles
      ! Note, tecplot will read across line breaks, so there is no need to put everythink on one record (line)
      write(iunit_state ,'(20(g15.7,1x))') t, dic, o2, nh3, c_l, c_d, n_d, pp(1:npp), c_pp(1:npp), cc(1:ncc), bac(1:nbac)
      write(iunit_rxn   ,'(20(g15.7,1x))') t, ((rpp(k,j),j=1,2),k=1,npp), ((rcc(k,j),j=1,npp+ncc+nbac),k=1,ncc), ((rbac(k,j),j=1,3),k=1,nbac)  
      write(iunit_sigDot, '(9(g15.7,1x))') t, sumSigDot(1:3), &
      sum(sumSigDot(1:3)), 100_mp*sumSigDot(1:3)/fZero(absZero,sum(sumSigDot(1:3)))
      write(iunit_misc,'(12(g15.7,1x))') t, iht*1000._mp/86400._mp, (omg_pp(i,1), omg_pp(i,2),i=1,npp) ! saves light at bottome of column converted to umole/m2/s and omg_pp

      ! store sumSigTot solution
      write(iunit_sigTot, '(8(g15.7,1x))') t, sumSigTot(1:3), sum(sumSigTot(1:3)), 100_mp*sumSigTot(1:3)/fZero(absZero,sum(sumSigTot(1:3)))
      ! Store inputs.  Note, there are the same for all solutions, so this file gets reround each time.
      call stepInput(t, dil_t0, dil_tf, dil_n, dil)         
      call stepInput(t, nh3I_t0, nh3I_tf, nh3I_n, nh3I) 
      write(iunit_input,'(3(g15.7,1x))') t, dil, nh3I 
      return   
   end subroutine storeSoln   

   subroutine storeTraits(iter, fVal)
      ! This routine saves the trait variables 
      ! Eact trait matrix/vector is staved in its own file using Tecplot syntax.  
      use realPrec
      use globalVars
      use functions, only: fZero
      implicit none
      integer,  intent(in)::  iter  ! The current iteration, but not really used in this version
      real(mp), intent(in)::  fVal  ! Value of objective function
      ! local declarations
      character(len=4) iterStr, nppStr, nccStr, nbacStr
      integer i
      ! 
      if (myRank /= 0) return ! only root process saves data
      
      write(iterStr,'(i4)') iter
      write(nppStr ,'(i4)') npp
      write(nccStr ,'(i4)') ncc
      write(nbacStr,'(i4)') nbac
    
      ! Save epp
      write(iunit_epp,'(a,f0.2,a)') 'Zone T = "epp iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_epp,'(a)') 'i = '//trim(adjustL(nppStr))//', j = 1, k = 1'
      write(iunit_epp,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1, npp
         write(iunit_epp,*) i, epp(i)
      end do      
      ! Save f_pp and phi_pp used to calculate omg_pp
      write(iunit_omgpp,'(a,f0.2,a)') 'Zone T = "omg_pp iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_omgpp,'(a)') 'i = '//trim(adjustL(nppStr))//', j = 1, k = 1'
      write(iunit_omgpp,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1,npp
         !write(iunit_omgpp,'(i0,2(1x,f0.8))') i, f_pp(i), phi_pp(i)
         write(iunit_omgpp,'(i0,3(1x,f0.8))') i, tOn_pp(i), tOff_pp(i), omg1Amp_pp(i) 
      end do
      ! Save ecc
      write(iunit_ecc,'(a,f0.2,a)') 'Zone T = "ecc iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_ecc,'(a)') 'i = '//trim(adjustL(nccStr))//', j = 1, k = 1'
      write(iunit_ecc,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1,ncc
         write(iunit_ecc,*) i, ecc(i)
      end do      
      ! Save omg_cc
      write(nppStr,'(i4)') npp+ncc+nbac
      write(iunit_omgcc,'(a,f0.2,a)') 'Zone T = "omg_cc iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_omgcc,'(a)') 'i = '//trim(adjustL(nccStr))//', j = 1, k = 1'
      write(iunit_omgcc,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1,ncc
         write(iunit_omgcc,*) i, omg_cc(i,1:npp+ncc+nbac)
      end do      
      ! Save ebac
      write(iunit_ebac,'(a,f0.2,a)') 'Zone T = "ebac iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_ebac,'(a)') 'i = '//trim(adjustL(nbacStr))//', j = 1, k = 1'
      write(iunit_ebac,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1, nbac
         write(iunit_ebac,*) i, ebac(i)
      end do      
      ! Save omg_bac
      write(iunit_omgbac,'(a,f0.2,a)') 'Zone T = "omg_bac iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal,'"'
      write(iunit_omgbac,'(a)') 'i = '//trim(adjustL(nbacStr))//', j = 1, k = 1'
      write(iunit_omgbac,'(a)') 'solutiontime = '//trim(adjustL(iterStr))
      do i=1, nbac
         write(iunit_omgbac,*) i, omg_bac(i,1:3)
      end do   
      if (genSurf) return ! don't store these, as they were already read in if a surface was requested.
      ! Aslo, write this to traits file so that they can be read easily if needed to another run.
      write (iunit_traits,'(a,f0.2)') 'All traits for iteration: '//trim(adjustL(iterStr))//'; fVal = ', fVal
      write (iunit_traits,*) npp, ncc, nbac
      write (iunit_traits,*) epp
      write (iunit_traits,*) ' '
      write (iunit_traits,*) tOn_pp
      write (iunit_traits,*) ' '
      write (iunit_traits,*) tOff_pp
      write (iunit_traits,*) ' '
      write (iunit_traits,*) omg1Amp_pp
      write (iunit_traits,*) ' '
      write (iunit_traits,*) ecc
      write (iunit_traits,*) ' '
      write (iunit_traits,*) omg_cc
      write (iunit_traits,*) ' '
      write (iunit_traits,*) ebac
      write (iunit_traits,*) ' '
      write (iunit_traits,*) omg_bac
      write (iunit_traits,*) ' '      
      return
   end subroutine storeTraits      

   subroutine initializeTraits(traitSpace, tsL, tsU)
      ! This routine is used to set traitSpace and it's lower and upper bounds from initial conditions
      ! It can also set traitSpace if the full eps and omg vectors/matrices had been read in from
      ! a previous stored simulation.
      ! This routine also sets the eps and omg vectors/matrices, which is only needed if the optimization is not run
      use globalVars
      implicit none
      real(mp), intent(out):: traitSpace(lenTS) ! traits in vector
      real(mp), intent(out):: tsL(lenTS) ! lower bounds on traitSpace
      real(mp), intent(out):: tsU(lenTS) ! upper bounds on ts
      ! local declarations
      ! nLinks is number of prey a cc can have, randomly generated here, but is a global variable
      integer iPrey, i, j, links
      real(mp) rnd, omg23
         
      if (readeps) then
         ! traits where already read in for first iteration.  These are actually be omg and not w, so need
         ! to extract w from omg_pp and omg_bac (not needed for omg_cc)
         ! Note, if omg has a dimenstion of three then, 
         !  w1 = omg1
         !  w2 = omg1 + omg2 (omg3 is not needed)
         ! In V3.0, omg_pp is now calculated from two other traits, f_pp and phi_pp (freq and phase), so only those are stored.
         ! In V4.0, omg_pp is now calculated from three traits, tOn_pp, t_Off_pp and omg1Amp_pp
         traitSpace(1:lenTS) = [(epp(1:npp)-epp_min)/(epp_max-epp_min), (tOn_pp(1:npp)-tOn_pp_min)/(tOn_pp_max-tOn_pp_min), (tOff_pp(1:npp)-tOff_pp_min)/(tOff_pp_max-tOff_pp_min),  &
                    omg1Amp_pp(1:npp), &
                    (ecc(1:ncc)-ecc_min)/(ecc_max-ecc_min), omg_cc(1:ncc,1:npp+ncc+nbac), (ebac(1:nbac)-ebac_min)/(ebac_max-ebac_min), &
                    ([omg_bac(i,1),omg_bac(i,1)+omg_bac(i,2)],i=1,nbac)]  ! not sure if this will work      
      else
         ! Note, only traitSpace needs to be set here, not eps or omg, but this may be hepful to understand wtf is going on.
         ! this is now longer used with hyperBOB version
         !epp  = epp_ini
         !f_pp(1:npp) = f_pp_ini
         !phi_pp(1:npp) = phi_pp_ini
         !ecc  = ecc_ini
         !omg_cc = omg_cc_ini
         !ebac = ebac_ini
         !omg_bac(1:nbac,1) = w_bac1_ini
         !omg_bac(1:nbac,2) = w_bac2_ini - w_bac1_ini ! Note, w_bac2_ini >= w_bac1_ini
         !omg_bac(1:nbac,3) = 1._mp - w_bac2_ini
         !traitSpace(1:lenTS) = [(epp(1:npp)-epp_min)/(epp_max-epp_min), (tOn_pp(1:npp)-tOn_pp_min)/(tOn_pp_max-tOn_pp_min), (tOff_pp(1:npp)-tOff_pp_min)/(tOff_pp_max-tOff_pp_min),  &
         !           omg1Amp_pp(1:npp), &
         !           (ecc(1:ncc)-ecc_min)/(ecc_max-ecc_min), omg_cc(1:ncc,1:npp+ncc+nbac),  (ebac(1:nbac)-ebac_min)/(ebac_max-ebac_min), &
         !           ([w_bac1_ini, w_bac2_ini], i=1,nbac)]        
      end if 
      ! Set the lower and upper boundaries on traitSpace
      tsL = 0._mp
      tsU = 1._mp
      return
   end subroutine initializeTraits   
   
   subroutine setTraits(traitSpace)
      ! This routine copies the traits held in traitSpace to the global vectors/matrices used in ODE and rxnRates, etc
      use realPrec
      use globalVars
      use interfaces
      implicit none
      integer i
      real(mp), intent(in):: traitSpace(lenTS)
      epp    = (epp_max-epp_min)*traitSpace(1:npp) + epp_min
      !f_pp   = (f_pp_max-f_pp_min)*traitSpace(npp+1:2*npp) + f_pp_min
      !phi_pp = (phi_pp_max-phi_pp_min)*traitSpace(2*npp+1:3*npp) + phi_pp_min
      tOn_pp   = (tOn_pp_max-tOn_pp_min)*traitSpace(npp+1:2*npp) + tOn_pp_min
      tOff_pp = (tOff_pp_max-tOff_pp_min)*traitSpace(2*npp+1:3*npp) + tOff_pp_min
      omg1Amp_pp = traitSpace(3*npp+1:4*npp)
      ecc    = (ecc_max-ecc_min)*traitSpace(4*npp+1:4*npp+ncc) + ecc_min
      omg_cc = reshape(traitSpace(4*npp+ncc+1 : 4*npp+ncc + ncc*(npp + ncc + nbac)),[ncc,npp+ncc+nbac])   
      if (binaryOMG) omg_cc = anint(omg_cc) ! This was added in V4.7 to force binary matrix
      ebac   = (ebac_max-ebac_min)*traitSpace(4*npp+ncc+ncc*(npp+ncc+nbac) + 1 : 4*npp+ncc+ncc*(npp + ncc + nbac) + nbac) + ebac_min
      do i=1,nbac
         call omega_free(3, traitSpace(4*npp+ncc+ncc*(npp+ncc+nbac)+nbac + 2*i-1 : 4*npp+ncc+ncc*(npp+ncc+nbac)+nbac + 2*i), omg_bac(i,1:3))
      end do      
      return
   end subroutine setTraits   
   
   subroutine omega_free (n, w, omg)
      ! this routine calclates omega (capital) from the lower case omega based on sorting
      ! of w then calculating distances between w_i.  
      ! Note, if there is only 1 sub reaction, then omg = w = 1, so there is no need
      ! to call this routine, but it is robust to such a call.
      ! This routine requires that sum(omg(1:n)) = 1
      use realPrec
      use interfaces
      implicit none
      integer, intent(in)  :: n              ! size of omg, which is one more than w
      real(dp), intent(in) :: w(:)           ! lower case omega, control variables
      real(dp), intent(out):: omg(:)         ! upper case omega, these sum to 1, which is where the nth omega comes from
      ! Local declaration
      integer i, l
      ! w1(1:n+1) = [0._dp, w(1:n-1), 1._dp] ! add 0 and 1 to ends
      real(dp) w1(n+1) ! Code may run faster if a parameter is used to specify size here.

      if (n == 1) then
         omg = 1.0_mp ! no need to have called this routine
         return
      end if      
      w1(1) = 0_mp
      w1(2:n) = w(1:n-1)
      w1(n+1) = 1_mp
      call insertion_sort(n, w1) ! sort w1
      omg(1:n) = w1(2:n+1) - w1(1:n)
      return
   end subroutine omega_free
   
   pure subroutine insertion_sort(n, a)
      ! uses insertion sorting
      ! from here https://rosettacode.org/wiki/Category:Sorting_Algorithms      
      ! but also see https://en.wikipedia.org/wiki/Sorting_algorithm
      use realPrec
      implicit none
      integer, intent(in)                    :: n
      real(dp), intent(in out), dimension(:) :: a
      real(dp) :: temp
      integer :: i, j
 
      do i = 2, n
         j = i - 1
         temp = a(i)
         do while (j>=1 .and. a(j)>temp) ! Strange, if I remove j>=1, it runs slower
            a(j+1) = a(j)
            j = j - 1
            if (j<=0) exit ! this prevents do while test of a(0)
         end do
         a(j+1) = temp
      end do
   end subroutine insertion_sort  
   
   subroutine rxnRates(t)
      ! This calculates all reaction rates
      use globalVars
      use functions
      use thermoData
      implicit none
      real(mp), intent(in) :: t ! current time (d)
      ! Local declrations
      integer i
      real(mp) I0t, sigDotRxn, sigDotPart, sumSigDotRxn, sumSigDoth2o, sumSigDotPart

      ! NOTE, mapx2names(x,n) must be called before this routine to update values of state variables from BALCOLI95
      ! ALSO, must call, geom, tempK, and pHtx.
      co2         = freeCO2 (T_K,is,pH)*dic ! Free CO2 concentration in mmol/m3 
      hco3        = freeHCO3(T_K,is,pH)*dic ! Free HCO3 concentration in mmol/m3   
      ! calculate necessary thermodynamic variables, such as standard free energy's of formation (since many routines will need)  See thermoData module
      dGf0_nh3    = dGf0(ammoniasp, T_K,is,pH)
      dGf0_h3po4  = dGf0(pisp,      T_K,is,pH)
      dGf0_h2o    = dGf0(h2osp,T_K,is,pH)
      dGf0_ch2o   = dGf0(ch2osp,    T_K,is,pH)   ! glucose, single carbon
      dGf0_bioS   = dGf0(yeastsp,   T_K,is,pH)  ! free energies of formation of species at current temp (K), is (M), and pH. (kJ/mol)
      dGf0_C_D    = dGf0_ch2o  ! Detritus C (use ch2o), N (use nh3) and P (use h3po4).
      dGf0_N_D    = dGf0_nh3
      dGf0_P_D    = dGf0_h3po4
      dGf0_o2aq   = dGf0(o2aqsp,    T_K,is,pH)
      dGf0_dic    = dGf0(co2totsp,  T_K,is,pH) ! or DIC
      !dGf0_hno3   = dGf0(nitratesp, T_K,is,pH) ! Not using yet
      
      ! V3.0: calculated omg_pp at current time
      !omg_pp(1:npp,1) = (sin(2_mp*pi*f_pp(1:npp)*t + phi_pp(1:npp)) + 1_mp)/2._mp
      ! V4.0: calculated omg_pp at current time using square wave
      do i=1,npp      
         omg_pp(i,1) = max( omg1Amp_pp(i)*( stepUp(mod(t,1._mp), tOn_pp(i), sigOmg_pp) + stepUp(mod(t,1._mp), tOff_pp(i), -sigOmg_pp) - 1._mp ), 0._mp )
      end do         
      omg_pp(1:npp,2) = 1._mp - omg_pp(1:npp,1)
      
      ! Calculate light at bottom, and average delI_phy.  Note, assume N_d is part of C_d, so don't double count it.
      !kwp = k_w + k_p*( sum(pp(1:npp)) + sum(c_pp(1:npp)) + sum(cc(1:ncc)) + sum(bac(1:nbac)) + c_d ) ! light absorption by all factors
      !kwp = k_w + k_p*( dot_product(omg_pp(1:npp,2),pp(1:npp)) + sum(cc(1:ncc)) + sum(bac(1:nbac)) ) + k_chla*dot_product(omg_pp(1:npp,1),pp(1:npp))! light absorption by all factors
      kwp = k_w + k_p*( dot_product(omg_pp(1:npp,2),pp(1:npp)) + sum(c_pp(1:npp)) + sum(cc(1:ncc)) + sum(bac(1:nbac)) ) + k_chla*dot_product(omg_pp(1:npp,1),pp(1:npp))! light absorption by all factors
      I0t = I0(t)
      Iht = I0t*exp(-kwp*depth) ! light a bottom of water column. This is a global variable for saving
      aveI = (I0t - Iht)/(kwp*depth) ! average light intensity for rectangular channel (mmol photons /m^2 /d (Not micromoles))
      !aveI = 2._mp*(Iht - I0t + I0t*depth*kwp)/(depth*kwp)**2 ! Average light intensity for a triangular channel
            
      ! Get the reaction rates associated with both pp and cc, as well as stoichiometric coef.
      sumSigDotRxn  = 0._mp  
      sumSigDotPart = 0._mp
      do i=1,npp
         call Phy_bioS(t, i, sigDotRxn, sigDotPart)
         sumSigDotRxn  = sumSigDotRxn  + sigDotRxn ! from all reactions, including the light harvesting part.
         sumSigDotPart = sumSigDotPart + sigDotPart ! entropy production from light absorption
      end do
      do i=1,ncc
         call Gz_bioS(t, i, sigDotRxn, sigDotPart)
         sumSigDotRxn  = sumSigDotRxn  + sigDotRxn ! from all cc reactions
         sumSigDotPart = sumSigDotPart + sigDotPart ! entropy production from light absorption
      end do      
      do i=1,nbac
         call Bac_BioS(t, i, sigDotRxn, sigDotPart)
         sumSigDotRxn  = sumSigDotRxn  + sigDotRxn ! from all bac reactions
         sumSigDotPart = sumSigDotPart + sigDotPart ! entropy production from light absorption
      end do            
      ! calculate entropy contributions (J/d/K) Note, if surA is set to 1, then units could be J/m2/d/K
      sumSigDotH2O  = - aveI*dGr_Ggamma*k_w*surA*depth/T_K ! absorption by water
      ! Particles are everything except the photoactive component of pp (i.e., Chl a, etc) which is already accounted for in sumSigmaDotRxn
      !sumSigDotPart = - aveI*dGr_Ggamma*k_p*( dot_product(pp(1:npp),omg_pp(1:npp,2)) + sum(c_pp(1:npp)) + sum(cc(1:ncc)) + sum(bac(1:nbac)) + c_d )*surA*depth/T_K
      ! Account for other matter that absorbs light, in this case only C_D remains.  This was removed in V3.9
      !sumSigDotPart = sumSigDotPart - aveI*dGr_Ggamma*k_p*c_d*surA*depth/T_K
      ! Saves this to global variable, which is saved for all t.
      sumSigDot(1) = sumSigDotRxn
      sumSigDot(2) = sumSigDotH2O
      sumSigDot(3) = sumSigDotPart
      return
   end subroutine rxnRates
   
   subroutine Phy_bioS(t, i_pp, sigDotRxn, sigDotPart)
      ! *** This is for V3.7, where the ratio of c_pp to pp is kept constant at k_pp, so both reactions 
      ! *** must be adjusted to insure a fixed ratio, but only if k_pp = 0; otherwise, the two reactions are not coupled.
      !!DIR$ OMP DECLARE SIMD
      ! This routine is use to model aerobic phytoplankton (1: Phy)
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only phosphate limitation as found in Siders Pond.
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      ! An alternative would be in add a derived type that could package all the parameters, but that would be for later development
      ! This version does not include a mortality reaction, so only co2 fixation to ch2o (cL) and converstion of that to phyS
      use realPrec
      use globalVars
      use functions , only: F_Thermo, fZero
      use thermoData, only: RkJ
      implicit none
      real(mp), intent(in)  :: t     ! time (d)
      integer , intent(in)  :: i_pp  ! The primary producter that needs reaction and entropy prod calculated.     
      real(mp), intent(out) :: sigDotRxn ! entropy production by imperfect light use and reactions.  
      real(mp), intent(out) :: sigDotPart ! entropy production associated with light absorption by particles
      ! local declarations
      real(mp) delI_pp ! Light captured for photosynthesis by the component of phytoplankton allocated to photosynthesis, given by omg(1)
      real(mp) dG0_r1A, dG_r1A ! free energy for reaction 1, CO2 fixation (kJ/mol or J/mmol)
      real(mp) dG0_r2A, dG_r2A, dG0_r2C, dG_r2C
      real(mp) concFac ! testing a concentration factor
      real(mp) ne, F_T, F_K, dG_r1, dG_r2, n1_pp
      real(mp) rxnRatio, sigDotBoth
      
      ! Begin calculations
      ! Reaction 1, CO2 fixation
      ! eps dic + n(1) hv -> eps(C_P + o2)  NOTE, now considering dic to include CO2 and biocarbonate.
      ! Calculate light captured by Phy photosynthetic apparatus only.
      delI_pp = k_chla*omg_pp(i_pp,1)*pp(i_pp)*aveI  ! note, only omg(1) fraction of phyS is allocated to light harvesting for chemosynthesis.
      ! calculate free energy associated with CH2O3 -> CH2O + O2
      dG0_r1A = (dGf0_ch2o + dGf0_o2aq) - dGf0_dic ! CO2 + H2O fixation into glucose and oxygen (unit carbon)   
      dG_r1A  = dG0_r1A + RkJ*T_K*( c_pp_lg(i_pp) + o2_lg - dic_lg  )
      ! calculate the mmol of photons needed to fix 1 mmol co2
      n1_pp = -dG_r1A/dGr_Ggamma ! this is the n1 coef.
      dG_r1 = -(1_mp - epp(i_pp))*dG_r1A ! free energy of reaction for reaction 1, CO2 fixation
      ne = 4 ! electrons transfered in catabolic reaction, this has been returned to 4
      F_T = F_Thermo(dG_r1, ne, T_K)
      F_K =  (co2+hco3)/(co2+hco3 + kappa*epp(i_pp)**4) ! added bicarbonate uptake to kinetics
      !F_K =  (co2)/(co2 + kappa*epp(i_pp)**4) ! only free CO2 uptake kinetics.  Huge effect on simulation
      rpp(i_pp,1) = delI_pp/n1_pp*F_T*F_K
         
      ! Calculate for phototroph biosynthesis reaction, r2
      ! (1+eps n2)C_Phy + eps(gam NH3 + del H2PO4) + [1 + eps(aA2 + n2 -1)]O2 -> eps PhyS + eps bA2 H2O + (1 + eps(n2-1))dic
      ! calculate reaction coefficients, but only anabolic reactions have them
      aA2_pp = (      - alf + 2._mp*bet + 3._mp*gam - 5._mp*del)/4._mp
      bA2_pp = (2._mp - alf                + 3._mp*gam + 3._mp*del)/2._mp
      ! anabolic reaction thermodynamics
      dG0_r2A = (dGf0_bioS +  bA2_pp*dGf0_h2o) - (dGf0_ch2o + gam*dGf0_nh3 + del*dGf0_h3po4 + aA2_pp*dGf0_o2aq)
      dG_r2A  = dG0_r2A + RkJ*T_K*( pp_lg(i_pp) - c_pp_lg(i_pp) - gam*nh3_lg - del*h3po4_lg - aA2_pp*o2_lg  )
      ! catabolic reaction thermodynamics
      dG0_r2C = dGf0_dic - ( dGf0_ch2o + dGf0_o2aq )
      dG_r2C = dG0_r2C + RkJ*T_K*( dic_lg - c_pp_lg(i_pp) - o2_lg )
      ! coupling of anabolic to catabolic reaction
      n2_pp(i_pp) = - dG_r2A/dG_r2C
      ! Free energy of combined whole reaction
      dG_r2 = (1._mp - epp(i_pp))*dG_r2C
      ! Reaction (2) rate
      F_T = F_Thermo(dG_r2, ne, T_K)
      F_K = ( c_pp(i_pp)*cell_F/(c_pp(i_pp)*cell_F + kappa*epp(i_pp)**4) )*( nh3/gam/(nh3/gam + kappa*epp(i_pp)**4) )*( o2/(o2 + kappa*epp(i_pp)**4)  )
      rpp(i_pp,2) = nuStar*omg_pp(i_pp,2)*epp(i_pp)**2*pp(i_pp)*F_T*F_K
      
      ! for v3.7, insure the ratio of r1 to r2 equals rxnRatio if k_pp /= 0
      if (k_pp >= epsilon(k_pp)) then
         rxnRatio = (1._mp/epp(i_pp) - n2_pp(i_pp) + k_pp) 
         if (rpp(i_pp,1) <= rxnRatio*rpp(i_pp,2)) then
            ! make r1 the driver
            rpp(i_pp,2) = rpp(i_pp,1)/rxnRatio
         else
            ! make r2 the driver
            rpp(i_pp,1) = rxnRatio*rpp(i_pp,2)
         end if
      end if
      ! Entropy production from reaction rates and delG
      ! for r1
      sigDotBoth = dG_r1A*( delI_pp/n1_pp - epp(i_pp)*rpp(i_pp,1) )*surA*depth/T_K  ! Entropy production by both rxn and particles
      sigDotRxn = rpp(i_pp,1)/fZero(absZero, delI_pp/n1_pp)*sigDotBoth ! if constraints are placed on rpp(i_pp,1), then treat as particle absorption
      sigDotPart = sigDotBoth - sigDotRxn ! Most of the time, entropy production is associated with reaction except when limited by kinetics.
      ! for r2
      sigDotRxn = sigDotRxn - rpp(i_pp,2)*dG_r2*surA*depth/T_K
      !sigDotPart = sigDotPart - aveI*dGr_Ggamma*k_p*(pp(i_pp)*omg_pp(i_pp,2) + c_pp(i_pp))*surA*depth/T_K ! For the non photosynthetic machinery and c_pp, it behaves as a particle      
      !sigDotPart = sigDotPart - aveI*dGr_Ggamma*k_p*(pp(i_pp)*omg_pp(i_pp,2))*surA*depth/T_K ! For the non photosynthetic machinery behaves as a particle	(V3.9 removed c_pp)
      sigDotPart = sigDotPart - aveI*dGr_Ggamma*k_p*(pp(i_pp)*omg_pp(i_pp,2) + c_pp(i_pp))*surA*depth/T_K ! For the non photosynthetic machinery behaves as a particle	(V4.3 added back c_pp)
      return
   end subroutine Phy_bioS

   subroutine Gz_bioS(t, i_cc, sigDotRxn, sigDotPart)
      !!DIR$ OMP DECLARE SIMD
      ! This routine is use to model aerobic grazers (3: Gz)
      ! See Notes for Frontiers MS
      ! Input and output variables are largely handled by declared types, so see module typeVars for details.
      use realPrec
      use functions, only: F_Thermo, fZero
      use globalVars
      use thermoData, only: RkJ
      implicit none
      real(mp), intent(in)  :: t     ! time (d)
      integer , intent(in)  :: i_cc  ! The consumer that needs reaction and entropy prod calculated.      
      real(mp), intent(out) :: sigDotRxn ! entropy production summed for all reactions.
      real(mp), intent(out) :: sigDotPart ! entropy production due to light absorption
      
      ! local declarations
      integer j
      real(mp) dG0_ri, dG_ri ! Free energy for reaction i, which is the consumption rate of BioS(i) (kJ/mol or J/mmol)
      real(mp) F_K, F_T
      real(mp) bCi, nei, c2ppj
      real(mp) wt, wDen ! weighting on prey uptake based on omg_cc and prey concentrations.
      
      ! zero out sigDotRxn
      sigDotRxn = 0._mp
      ! Catabolic stiochiometric coef. Note, aCi_cc depends on ecc, so must save for reaction ODE/PDE
      aCi_cc(i_cc) = (  4._mp + alf - 2._mp*bet - 3._mp*gam + 5._mp*del - 4._mp*ecc(i_cc))/4._mp ! global variable
      bCi          = (- 2._mp + alf             - 3._mp*gam - 3._mp*del                  )/2._mp
      
      ! loop over all prey and calculate rxn rates and entropy production
      ! First get denominator of weighting (prey preference) for i_cc consumer
      wDen = dot_product(omg_cc(i_cc,1:npp),pp(1:npp)) + dot_product(omg_cc(i_cc,npp+1:npp+ncc),cc(1:ncc)) + dot_product(omg_cc(i_cc,npp+ncc+1:npp+ncc+nbac),bac(1:nbac))
      wDen = fZero(absZero, wDen) ! unlikely, but make sure this is not zero
      do j=1,npp ! first loop over pp
         c2ppj = c_pp(j)/pp(j) ! ratio of C_pp to pp.  This ratio can get ugly, may wish to constrain
         ! Free energy of reactions.  Here, C_pp to just compusted to dic instead of C_L like in the Siders model.
         dG0_ri = ecc(i_cc)*dGf0_bioS + (1._mp - ecc(i_cc))**2*dGf0_dic + ecc(i_cc)*(1._mp - ecc(i_cc))*dGf0_c_d  &
                + del*(1._mp - ecc(i_cc))*dGf0_h3po4 + gam*(1._mp - ecc(i_cc))*dGf0_nh3 + bCi*(1._mp - ecc(i_cc))*dGf0_h2o + c2ppj*dGf0_dic &
                - ( dGf0_bioS + c2ppj*dGf0_ch2o + (aCi_cc(i_cc)*(1._mp - ecc(i_cc)) + c2ppj)*dGf0_o2aq )
         ! Account for the concentrations.
         dG_ri  = dG0_ri + RkJ*T_K*( ecc(i_cc)*cc_lg(i_cc) + (1._mp - ecc(i_cc))**2*dic_lg + ecc(i_cc)*(1._mp - ecc(i_cc))*C_D_lg &
                + del*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*h3po4_lg + ecc(i_cc)*P_D_lg) &
                + gam*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*nh3_lg   + ecc(i_cc)*N_D_lg) &
                + c2ppj*dic_lg &
                - (pp_lg(j) + c2ppj*c_pp_lg(j) + (aCi_cc(i_cc)*(1._mp - ecc(i_cc)) + c2ppj)*o2_lg ) )
         nei = 4._mp ! Electron transfer same for all reactions. it was 4._mp*aCi_cc(i_cc), but this is fueled by CH2O oxidation, so set at 4.         
         F_T = F_Thermo(dG_ri, nei, T_K) ! Thermo driver         
         F_K = o2/(o2 + kappa*ecc(i_cc)**4)*pp(j)/(pp(j) + kappa*ecc(i_cc)**4) ! kinetic driver
         wt = omg_cc(i_cc,j)*pp(j)/wDen ! weighting for uptake of prey pp(j)
         rcc(i_cc,j) = nuStar*ecc(i_cc)**2*wt*cc(i_cc)*F_K*F_T ! rate of pp(j) uptake
         ! Entropy Production
         sigDotRxn = sigDotRxn - rcc(i_cc,j)*dG_ri*surA*depth/T_K
      end do      
      ! now loop over cc(i_cc) eating all other cc including itself.
      do j=1,ncc ! first loop over pp
         ! Free energy of reactions.  Here, C_pp to just compusted to dic instead of C_L like in the Siders model.
         dG0_ri = ecc(i_cc)*dGf0_bioS + (1._mp - ecc(i_cc))**2*dGf0_dic + ecc(i_cc)*(1._mp - ecc(i_cc))*dGf0_c_d  &
                + del*(1._mp - ecc(i_cc))*dGf0_h3po4 + gam*(1._mp - ecc(i_cc))*dGf0_nh3 + bCi*(1._mp - ecc(i_cc))*dGf0_h2o  &
                - ( dGf0_bioS + aCi_cc(i_cc)*(1._mp - ecc(i_cc))*dGf0_o2aq )
         ! Account for the concentrations.
         dG_ri  = dG0_ri + RkJ*T_K*( ecc(i_cc)*cc_lg(i_cc) + (1._mp - ecc(i_cc))**2*dic_lg + ecc(i_cc)*(1._mp - ecc(i_cc))*C_D_lg &
                + del*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*h3po4_lg + ecc(i_cc)*P_D_lg) &
                + gam*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*nh3_lg   + ecc(i_cc)*N_D_lg) &
                - (cc_lg(j) + aCi_cc(i_cc)*(1._mp - ecc(i_cc))*o2_lg ) )
         nei = 4._mp ! Electron transfer same for all reactions. it was 4._mp*aCi_cc(i_cc), but this is fueled by CH2O oxidation, so set at 4.         
         F_T = F_Thermo(dG_ri, nei, T_K) ! Thermo driver         
         F_K = o2/(o2 + kappa*ecc(i_cc)**4)*cc(j)/(cc(j) + kappa*ecc(i_cc)**4) ! kinetic driver
         wt = omg_cc(i_cc,npp+j)*cc(j)/wDen ! weighting for uptake of prey cc(j)
         rcc(i_cc,npp+j) = nuStar*ecc(i_cc)**2*wt*cc(i_cc)*F_K*F_T ! reaction rate of cc(j) uptake
         ! Entropy Production
         sigDotRxn = sigDotRxn - rcc(i_cc,npp+j)*dG_ri*surA*depth/T_K
      end do
      ! now loop over cc(i_cc) eating all bac.
      do j=1,nbac ! loop over bac
         ! Free energy of reactions.  Here, C_pp to just compusted to dic instead of C_L like in the Siders model.
         dG0_ri = ecc(i_cc)*dGf0_bioS + (1._mp - ecc(i_cc))**2*dGf0_dic + ecc(i_cc)*(1._mp - ecc(i_cc))*dGf0_c_d  &
                + del*(1._mp - ecc(i_cc))*dGf0_h3po4 + gam*(1._mp - ecc(i_cc))*dGf0_nh3 + bCi*(1._mp - ecc(i_cc))*dGf0_h2o  &
                - ( dGf0_bioS + aCi_cc(i_cc)*(1._mp - ecc(i_cc))*dGf0_o2aq )
         ! Account for the concentrations.
         dG_ri  = dG0_ri + RkJ*T_K*( ecc(i_cc)*cc_lg(i_cc) + (1._mp - ecc(i_cc))**2*dic_lg + ecc(i_cc)*(1._mp - ecc(i_cc))*C_D_lg &
                + del*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*h3po4_lg + ecc(i_cc)*P_D_lg) &
                + gam*(1._mp - ecc(i_cc))*((1._mp - ecc(i_cc))*nh3_lg   + ecc(i_cc)*N_D_lg) &
                - (bac_lg(j) + aCi_cc(i_cc)*(1._mp - ecc(i_cc))*o2_lg ) )
         nei = 4._mp ! Electron transfer same for all reactions. it was 4._mp*aCi_cc(i_cc), but this is fueled by CH2O oxidation, so set at 4.         
         F_T = F_Thermo(dG_ri, nei, T_K) ! Thermo driver         
         F_K = o2/(o2 + kappa*ecc(i_cc)**4)*bac(j)/(bac(j) + kappa*ecc(i_cc)**4) ! kinetic driver
         wt = omg_cc(i_cc,npp+ncc+j)*bac(j)/wDen ! weighting for uptake of prey bac(j)
         rcc(i_cc,npp+ncc+j) = nuStar*ecc(i_cc)**2*wt*cc(i_cc)*F_K*F_T ! reaction rate of bac(j) uptake
         ! Entropy Production
         sigDotRxn = sigDotRxn - rcc(i_cc,npp+ncc+j)*dG_ri*surA*depth/T_K
      end do
      ! Entropy production due to light absorption
      sigDotPart = - aveI*dGr_Ggamma*k_p*cc(i_cc)*surA*depth/T_K
      
      return
   end subroutine Gz_BioS   
   
   subroutine Bac_BioS(t, i_bac, sigDotRxn, sigDotPart)
      ! This routine is use to model heterotrophic bacteria
      ! See C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\MEPphotoautotrophGrowh_v2.docx
      ! Also see Notes for Frontiers MS
      ! This model will assume only NH3 limitation.
      ! Input and output variables are largely handled by global variables, see globalVars.
      ! Three reactions are considered here: 1) Bacterial growth, 2) C_D decomposition, 3) N_D decomposition
      use realPrec
      use globalVars
      use functions , only: F_Thermo, stepUp
      use thermoData, only: RkJ
      implicit none
      real(mp), intent(in)  :: t     ! time (d)
      integer , intent(in)  :: i_bac  ! The bacteria that needs reaction and entropy prod calculated.      
      real(mp), intent(out) :: sigDotRxn ! entropy production summed for all reactions. 
      real(mp), intent(out) :: sigDotPart ! entropy production due to light absorption. 
      
      ! local declarations
      real(mp) dG0_rA1, dG_rA1, dG0_rC1, dG_rC1 ! free energy for anabolic and catabolic reactions for r1
      real(mp), parameter:: ln106 = -6._mp*log(10._mp) ! this is the natural logorithm of 10^-6 to convert uM to M
      real(mp) F_T, F_K, ne, dG, bA1_bac

      ! Calculate heterogrophic growth, r1
      ! C_L + eps(gam NH3 + del H2PO4) + (1-eps)o2 -> eps aA1 Bac + eps bA1 H2O + [2 - eps(aA1-1)]dic
      ! calculate reaction coefficients, but only anabolic reactions have them
      
      aA1_bac = (4._mp + 3._mp*gam - 5._mp*del)/(4._mp + alf - 2._mp*bet) ! This one is needed by the ODE routine, so make it globalVar
      bA1_bac = (4._mp - 2._mp*alf + 9._mp*gam - 3._mp*bet*gam + del + 4._mp*alf*del - 3._mp*bet*del)/(4._mp + alf - 2._mp*bet)
      ! anabolic reaction thermodynamics
      dG0_rA1 = (aA1_bac*dGf0_bioS + (1._mp-aA1_bac)*dGf0_dic + bA1_bac*dGf0_h2o) - (dGf0_ch2o + gam*dGf0_nh3 + del*dGf0_h3po4)
      dG_rA1  = dG0_rA1 + RkJ*T_K*( aA1_bac*bac_lg(i_bac) + (1._mp-aA1_bac)*dic_lg - c_l_lg - gam*nh3_lg - del*h3po4_lg )

      ! catabolic reaction thermodynamics
      dG0_rC1 = dGf0_dic - ( dGf0_ch2o + dGf0_o2aq )
      dG_rC1  = dG0_rC1 + RkJ*T_K*( dic_lg - c_l_lg - o2_lg )
      ! Free energy of combined whole reaction
      dG = ebac(i_bac)*dG_rA1 + (1._mp - ebac(i_bac))*dG_rC1
      ! Reaction (1) rate
      ne = 4. ! electrons transfered in catabolic reaction
      F_T = F_Thermo(dG, ne, T_K)
      F_K = ( C_L/(C_L + kappa*ebac(i_bac)**4) )*( nh3/gam/(nh3/gam + kappa*ebac(i_bac)**4)  )*( o2/(o2 + kappa*ebac(i_bac)**4)  ) ! don't include NH3 limatations yet
      rbac(i_bac,1) = nuStar*ebac(i_bac)**2*omg_bac(i_bac,1)*bac(i_bac)*F_T*F_K
      ! Entropy production
      sigDotRxn = - rbac(i_bac,1)*dG*surA*depth/T_K
      
      ! Reaction 2, C_D decompositoin to C_L.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! C_D -> C_L
      dG = RkJ*T_K*log(C_L/C_D) ! This will be negligable, but easy to calculate so keep
      F_K = C_D/(C_D + kappa*ebac(i_bac)**4)
      F_T = stepUp(-dG, 0.3_mp, 20._mp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      rbac(i_bac,2) = nuDet*ebac(i_bac)**2*omg_bac(i_bac,2)*bac(i_bac)*F_T*F_K
      sigDotRxn = sigDotRxn - rbac(i_bac,2)*dG*surA*depth/T_K
      
      ! Reaction 3, N_D decompositoin to nh3.  Note, there is not thermodyamic driver for this, BUT detritus degridation uses different nuStar: nuDet
      ! N_D -> NH3
      dG = RkJ*T_K*log(nh3/N_D) ! This will be negligable, but easy to calculate so keep
      F_K = N_D/(N_D + kappa*ebac(i_bac)**4)
      F_T = stepUp(-dG, 0.3_mp, 20._mp) ! want F_T to go to zero when delG >= 0.  However, use a smooth step.
      rbac(i_bac,3) = nuDet*ebac(i_bac)**2*omg_bac(i_bac,3)*bac(i_bac)*F_T*F_K
      sigDotRxn = sigDotRxn - rbac(i_bac,3)*dG*surA*depth/T_K
      ! Entropy production due to light absorption
      sigDotPart = - aveI*dGr_Ggamma*k_p*bac(i_bac)*surA*depth/T_K      
      return
   end subroutine Bac_bioS      
   
   !-----------------------------------------------------------------------
   ! Routines below are provide the definition of the ODEs for BiM
   !-----------------------------------------------------------------------   
   subroutine ODEs(nx, t, x ,xdot, ierr, rpar, ipar)
      ! This routine contains the state equations
      ! for which a solution is being seeked.
      use mpi
      use realPrec
      use globalVars 
      use thermoData, only: co2g2l, solO2
      implicit none
      integer,  intent(in) :: nx       ! number of differential equations
      real(mp), intent(in) :: t        ! time (d)
      real(mp), intent(in) :: x(nx)    ! state vector
      real(mp), intent(out):: xdot(nx) ! derivative of x set by this routine
      integer,  intent(out):: ierr     ! set to /= 0 to indicate an error and have BiM terminate
      real(mp), intent(in) :: rpar(*)  ! user passed real vector
      integer,  intent(in) :: ipar(*)  ! user passed integer vector
         
      ! local declarations
      integer i, j
      real(mp) o2Eq, co2Eq
      real(mp) ODE_ti
      
      ierr = 0 ! set default to no errors
      ! First see if maximum allowed time has been exceeded
      ODE_ti = MPI_Wtime() ! note, ODE_t0 is a global variable set in integratePDE routine
      if ((ODE_ti-ODE_t0)/86400._mp > maxODEtime) then
         ierr = 1 ! set to 1 to exit BiM integration        
         return
      end if
      
      ! map the x vector into useful names, and others that are needed by rxnRates and here
      call mapx2names(x, nx)
      ! V3.7, get current concentrations of nh3 and c_L in the input stream, as well as the dilution rate
      ! Calculate reaction rates and entropy production for reactions and light absorption by water and particles
      call rxnRates(t)
      call stepInput(t, dil_t0, dil_tf, dil_n, dil) ! get current dilution rate
      call stepInput(t, nh3I_t0, nh3I_tf, nh3I_n, nh3I)
      ! Define ODE equations.  
      ! dic(1)
      co2Eq = pCO2*co2g2l(T_K, is, pH) ! this is the CO2 concentration (mmol/m3) in equilibrium with air (pCO2 in atm)
      xdot(1) = pV_co2*( co2Eq - co2 )/depth & ! CO2 gas exchange
              + dil*(dicI - dic) & ! chemostat/diffusion transport
              - dot_product(epp(1:npp), rpp(1:npp,1)) & ! co2 fixation by pp
              + dot_product(1._mp + epp(1:npp)*(n2_pp(1:npp) - 1._mp), rpp(1:npp,2)) ! pp respiration                 
      do j=1,npp
         ! cc respiration by consuming pp_j
         xdot(1) = xdot(1) + dot_product((1._mp - ecc(1:ncc))**2 + c_pp(j)/pp(j), rcc(1:ncc,j) )
      end do
      do j=1,ncc
         ! cc respiration by consuming cc_j
         xdot(1) = xdot(1) + dot_product((1._mp - ecc(1:ncc))**2, rcc(1:ncc,npp+j) )                  
      end do      
      do j=1,nbac
         ! cc respiration by consuming bac_j
         xdot(1) = xdot(1) + dot_product((1._mp - ecc(1:ncc))**2, rcc(1:ncc,npp+ncc+j) )                  
      end do     
      ! bac respiration by consuming cL
      xdot(1) = xdot(1) + dot_product((2._mp - ebac(1:nbac)*(aA1_bac + 1._mp)), rbac(1:nbac,1) )                   
      
      ! o2(2)
      o2Eq = pO2*solO2(T_K, is, pH) ! O2 in equilbrium with atmosphere O2 pressure set in pO2      
      xdot(2) = pV_o2*( o2Eq - o2 )/depth & ! O2 gas exchange
              + dil*(o2I - o2) & ! chemostat/diffusion transport
              + dot_product(epp(1:npp), rpp(1:npp,1)) & ! o2 production from autotrophs
              - dot_product(1._mp + epp(1:npp)*(aA2_pp + n2_pp(1:npp) - 1._mp), rpp(1:npp,2)) ! pp respiration
      do j=1,npp
         ! cc respiration of pp_j
         xdot(2) = xdot(2) - dot_product(aCi_cc(1:ncc)*(1._mp - ecc(1:ncc)) + c_pp(j)/pp(j), rcc(1:ncc,j) )
      end do
      do j=1,ncc
         ! cc respiration of cc_j
         xdot(2) = xdot(2) - dot_product(aCi_cc(1:ncc)*(1._mp - ecc(1:ncc)), rcc(1:ncc,npp+j) )                  
      end do   
      do j=1,nbac
         ! cc respiration of bac_j
         xdot(2) = xdot(2) - dot_product(aCi_cc(1:ncc)*(1._mp - ecc(1:ncc)), rcc(1:ncc,npp+ncc+j) )                  
      end do         
      ! bac respiration by consuming cL
      xdot(2) = xdot(2) - dot_product(1._mp - ebac(1:nbac), rbac(1:nbac,1) )                   
      
      ! nh3(3)
      xdot(3) = dil*(nh3I - nh3) & ! chemostat/diffusion transport
              - dot_product(epp(1:npp)*gam, rpp(1:npp,2)) ! uptake by pp
      do j=1,npp+ncc+nbac
         ! remin by cc
         xdot(3) = xdot(3) + dot_product(gam*(1._mp - ecc(1:ncc))**2, rcc(1:ncc,j) )
      end do
      ! uptake by bac, and remin of Nd by bac
      xdot(3) = xdot(3) - gam*dot_product( ebac(1:nbac),rbac(1:nbac,1) ) + sum(rbac(1:nbac,3))
      
      ! c_L(4)
      xdot(4) = dil*(c_LI - c_L) & ! chemostat/diffusion transport
              + sum(rbac(1:nbac,2) - rbac(1:nbac,1)) ! remin of c_D and uptake by bac
      
      ! c_d(5)
      xdot(5) = dil*(c_dI - c_d) ! chemostat/diffusion transport
      do j=1,npp+ncc+nbac
         ! cc consuming pp_j, cc_j and bac_j
         xdot(5) = xdot(5) + dot_product(ecc(1:ncc)*(1._mp - ecc(1:ncc)), rcc(1:ncc,j) )
      end do
      xdot(5) = xdot(5) - sum(rbac(1:nbac,2)) ! bac remin of c_d
      
      ! n_d(6)
      xdot(6) = dil*(n_dI - n_d)  ! chemostat/diffusion transport
      do j=1,npp+ncc+nbac
         xdot(6) = xdot(6) + dot_product(gam*ecc(1:ncc)*(1._mp - ecc(1:ncc)), rcc(1:ncc,j) )
      end do
      xdot(6) = xdot(6) - sum(rbac(1:nbac,3)) ! bac remin of n_d

      ! pp(i)
      do i=1,npp
         xdot(nconc+i) = dil*(ppI - pp(i)) + epp(i)*rpp(i,2) - sum(rcc(1:ncc,i))                        
      end do
   
      ! c_pp(i)
      do i=1,npp
         xdot(nconc+npp+i) = dil*(c_ppI - c_pp(i)) + epp(i)*rpp(i,1) &
                           - (1._mp + epp(i)*n2_pp(i))*rpp(i,2) - c_pp(i)/pp(i)*sum(rcc(1:ncc,i))
      end do
      
      ! cc(i)
      do i=1,ncc
         xdot(nconc+2*npp+i) = dil*(ccI - cc(i)) + ecc(i)*sum(rcc(i,1:npp+ncc+nbac)) &
                             - sum(rcc(1:ncc,npp+i))
      end do
      
      ! bac(i)
      do i=1,nbac
         xdot(nconc+2*npp+ncc+i) = dil*(bacI - bac(i)) + ebac(i)*aA1_bac*rbac(i,1) - sum(rcc(1:ncc,npp+ncc+i))
      end do
      
      ! Last three states are integrated intropy production for reactions and light absorption by water and particles
      xdot(nx-2:nx) = sumSigDot(1:3) 
      
      !if (any(isnan(xdot))) then
      !   write(*,*) 'NaN in xdot in PDEs routine'
      !   read(*,'(a)') junk
      !end if      
      return
   end subroutine ODEs

   subroutine ompJAC(n,t,x,jacmat,ldjac,ierr,rpar,ipar)
      ! This calculates a numberical jacobian using the same method in BiM code
      ! but uses OpenMP to speed things up.
      use realPrec
      implicit none
      integer,  intent(in)    :: n              ! Size of system
      real(mp), intent(in)    :: t              ! time
      real(mp), intent(inout) :: x(n)           ! state variables. Not changed, but do need to temp alter.
      real(mp), intent(out)   :: jacmat(ldjac,n)! Jacobian output
      integer,  intent(in)    :: ldjac          ! Size of jacmat
      integer,  intent(out)   :: ierr           ! if errors occur
      real(mp), intent(inout) :: rpar(*)        ! User defined vector
      integer,  intent(inout) :: ipar(*)        ! User defined vector
      ! Local declarations
      integer i, j
      real(mp) xsave, delt, uround 
      real(mp) xdot(n), xdot0(n)
  
      ! Values +1 is still 1. BiM uses 1.0d-16, but could use epsilon here.
      uround = 1.0d-16 
      
      ! get the value of xdot at current x
      call ODEs(n, t, x ,xdot0, ierr, rpar, ipar)
      continue
      ! Begin parallel calculation of jacobian matrix
      !$OMP PARALLEL DO firstprivate(x) private(i,j,xsave,delt,xdot) schedule(static)      
      do i=1,n
         xsave = x(i)
         delt = sqrt(uround*max(1.D-5,abs(xsave)))
         x(i) = xsave+delt
         call ODEs(n, t, x ,xdot, ierr, rpar, ipar)
         do j=1,n
            jacmat(j,i)=(xdot(j)-xdot0(j))/delt
         end do
         x(i)=xsave
      end do   
      !$OMP END PARALLEL DO      
      ierr = 0
      return
   end subroutine ompJAC

   subroutine stepInput(t, f_t0, f_tf, f_n, f_t)
      ! This routine steps up (or down) a feed input based on
      ! the starting and ending values and the number of steps
      ! that are to occur between t0 and tf of the simulation      
      use realPrec
      use globalVars, only: t0, tf, sigmaSI
      implicit none
      real(mp), intent(in) :: t     ! current time (d)
      real(mp), intent(in) :: f_t0  ! value of variable at t0 (units depend on variable passed)
      real(mp), intent(in) :: f_tf  ! value of variable at tf 
      integer , intent(in) :: f_n   ! number of steps between t0 and tf.  0 means no steps (f_t = f_t0)
      real(mp), intent(out):: f_t   ! value of f at time t
      ! local declartions
      integer i
      real(mp) tStep, tDel, fStep
      
      f_t = f_t0
      if (f_n == 0) return ! no steps in f
      
      tDel  = (tf - t0)/(real(f_n) + 1.0_mp) ! number of days between steps in f
      fStep = (f_tf - f_t0)/real(f_n) ! increase in f at each step
      tStep = 0.0_mp
      do i=1, f_n
         tStep = tStep + tDel
         f_t = f_t + fStep/(1._mp + exp(-(sigmaSI*(t - tStep))))
      end do      
      return
   end subroutine stepInput         
   
   Subroutine parSur (time,dlat,solarC,par)
      use realPrec
      Implicit None
      real(mp) time, dlat, par, solarC
      !
      !     This routine the light level at the surface of the earth given
      !     time and latitude (and some other stuff)
      !
      !     USAGE - Call parSur (time,dlat,solarC,par)
      !
      !     INPUT
      !        TIME     REAL*8. Time of year in days (Julian days)
      !        dlat     Real*8  Latitude in degrees
      !        solarC   Real*8  Solar constant, in whatever units you want par in.
      !                         e.g.  1353. W/m2
      !                               6216. microEinstein/s/m2 Note, 0.75 of this seems to give more reasonable values (2000 vs 2650 on day 177.5)
      !
      !     OUTPUT
      !        par      REAL*8. photosynthetic available radiation.
      !
      !     MODIFIED:   31Jan97, 7Mar2016

      !     Local Declarations
      real(mp) pi, rlat, hour, dec, ha, cosz, aa           
      !real(mp) sigma ! this is an exponential weighting coef. to allow a smooth transition to night.
      !      
      !     Declination and cos of Zenith angle.
      !     Ref: Brock (1981) Ecol.Model. 14:1-19, changed to radians
      pi = 3.141592654D0
      rlat = dlat*pi/180.D0
      hour = (time - aint(time))*24.D0
      dec = 0.4093D0*sin(2.*pi*(284.D0+time)/365.D0)
      ha = (hour - 12.D0)*pi/12.D0                                           
      cosz = Dsin(dec)*Dsin(rlat) + Dcos(dec)*Dcos(rlat)*Dcos(ha)           
      !     See if it is night time    
      !If (cosz <= 0._mp) cosz = 0.0D0

      ! see if exponential weighting that removes the discontinuity helps numerical solution
      ! see C:\Users\jvallino\Projects\NSF-GeoBio_Apr2015\Modeling\Mma\SmoothPARfcn.nb
      !sigma = 50._mp
      !If (cosz <= 0._mp) then
      !   cosz = 0.0D0
      !else
      !   cosz = cosz*( 1._mp - exp(-sigma*cosz) )/( 1._mp + exp(-sigma*cosz) )
      !end if            
      cosz = max(0._mp, cosz) ! the above exponential function did not seem to matter, so just use this
      
      ! PAR to Total radiation at surface: a = Epar(0+)/Etot(0+)
      ! Ref: Baker & Frouin (1987) L&O 32:1370-1377     
      aa = 0.45 
      par = solarC*aa*cosz      
      return
   end subroutine parSur

	subroutine carbonate(dic, alk, t, s, h2co3, co3, h, hco3)
	   implicit none
	   real(8) dic, alk, t, s, h2co3, co3, h, hco3

	   ! This routine solves the unknowns in the carbonate system given DIC and alkalinity.
	   ! The equations solved (via Mathematica) are:
	   !
	   ! Solve[{h2co3 k1 == hco3 h, k2 hco3 == co3 h, alk == hco3 + 2co3,
	   !         dic == hco3 + co3 + h2co3}, {h2co3, hco3, co3, h}]
	   !
	   ! Input
	   !	dic	DIC concentration (uM)
	   !	alk	Alkalinity (uM eq)
	   !	t		Temperature (K)
	   !	s		Salinity (S)
	   !
	   ! Output
	   !	h2co3	uM
	   !	co3	uM
	   !	h		uM Note, to get pH, first divide by 10^6 to get into molar conc.
	   !	hco3	uM
	   !
	   ! History
	   !	31 Jan 02 - First written.
	   ! Local declarations
	   real(8) k1, k2

	   ! Get the equilibrium constants at T and S
	   call  k1k2(t,s,k1,k2)

	   h2co3 =  (-alk + 2.*dic + (dic*k1)/(-k1 + 4.*k2) & 
                - Sqrt(4.*dic**2*k1**2 + 4.*alk*(-alk + 2.*dic)*k1*(-k1 + 4.*k2))/(2.*(-k1 + 4.*k2)))/2.

	   co3 = (alk + (dic*k1)/(-k1 + 4*k2) - Sqrt(4*dic**2*k1**2 + 4*alk*(-alk + 2*dic)*k1*(-k1 + 4*k2))/(2.*(-k1 + 4*k2)))/2.

	   h = (-(alk*k1) + 2.*dic*k1 + (dic*k1**2)/(-k1 + 4.*k2) - (4.*dic*k1*k2)/(-k1 + 4.*k2) & 
		   - (k1*Sqrt(4.*dic**2*k1**2 + 4.*alk*(-alk + 2.*dic)*k1*(-k1 + 4.*k2)))/(2.*(-k1 + 4.*k2)) & 
		   + (2.*k2*Sqrt(4.*dic**2*k1**2 + 4.*alk*(-alk + 2.*dic)*k1*(-k1 + 4.*k2)))/(-k1 + 4.*k2))/(2.*alk)

	   hco3 = (-2.*dic*k1 + Sqrt(4.*dic**2*k1**2 + 4.*alk*(-alk + 2.*dic)*k1*(-k1 + 4.*k2)))/(2.*(-k1 + 4.*k2))
	   return
   end subroutine carbonate
   
   Subroutine k1k2(T,S,K1,K2)
	   implicit none
	   real(8) T, S, K1, K2 
	   ! This routine returns the K1 and K2 equilibrium constants for the carbonate
	   ! system at a given termperature and salinity. F. J. Millero. Thermodynamics of the carbon dioxide system in 
	   ! the oceans. Geochim.Cosmochim.Acta 59 (4):661-677, 1995. Used Equations (41)-(44).
	   !
	   ! Input
	   !	T		temperature (K)
	   !	S		salinity (PSU)
	   !
	   ! Output
	   !	K1		Equilibrium constant for CO2 + H2O -> H(+) + HCO3(-) (uM)
	   !	K2		Equilibrium constant for HCO3(-) -> H+ + CO3(2-) (uM)

	   ! Local declarations

	   if (s >= 5) then
		   ! Use normal salinity equations (41) and (42)
		   K1 = 1.0d6*EXP( 2.83655D0 - 2307.1266D0/T - 1.5529413D0*log(T) &
				   + (-0.20760841D0 - 4.0484D0/T)*sqrt(S) + 0.08468345D0*S - 0.00654208D0*S*sqrt(S) )
		   K2 = 1.0d6*EXP( -9.226508d0 - 3351.6106d0/T - 0.2005743d0*log(T) &
				   + (-0.106901773d0 - 23.9722d0/T)*sqrt(S) + 0.1130822d0*S - 0.00846934*S*sqrt(S) )
	   else 
		   ! Use low salinity equations (43) and (44)
		   K1 = 1.0d6*EXP( 290.9097d0 - 14554.21d0/T - 45.0575d0*log(T) &
				   + (-228.39774d0 + 9714.36839d0/T + 34.485796d0*log(T))*sqrt(S) &
				   + (54.20871d0 - 2310.48919d0/T - 8.19515*log(T))*S &
				   + (-3.969101d0 + 170.22169d0/T + 0.603627d0*log(T))*S*sqrt(S) - 0.00258768*S**2 )
		   K2 = 1.0d6*EXP( 207.6548d0 - 11843.79d0/T - 33.6485d0*log(T) &
				   + (-167.69908d0 + 6551.35253d0/T + 25.928788d0*log(T))*sqrt(S) &
				   + (39.75854d0 - 1566.13883d0/T - 6.171951d0*log(T))*S &
				   + (-2.892532d0 + 116.270079d0/T + 0.45788501d0*log(T))*S*sqrt(S) - 0.00613142*S**2 )
	   end if
	   return
   end subroutine k1k2
   
   subroutine isort2D (n, y)
      ! This sorts an integer 2D array, y(n,2) in ascending order based on the
      ! first column. The second column corresponds to an index for the first
      implicit none
      integer, intent(in)   :: n ! Row dimension of y
      integer, intent(inout):: y(n,2) ! The first column is sorted and the second column moves accordingly.
      ! local declarations
      integer i,j
      integer ymin(2)
    
      j = 1
      do while (j < n)
         ymin = y(j,:)
         do i=j+1,n
            if (y(i,1) < ymin(1)) then
               y(j,:) = y(i,:)
               y(i,:) = ymin
               ymin = y(j,:)
            end if
         end do
         j = j+1
      end do
      return
   end subroutine isort2D
   
   subroutine dsort2D (n, y)
      ! This sorts an real(mp) 2D array, y(n,2) in ascending order based on the
      ! first column. The second column corresponds to an index for the first
      use realPrec
      implicit none
      integer , intent(in)   :: n ! Row dimension of y
      real(mp), intent(inout):: y(n,2) ! The first column is sorted and the second column moves accordingly.
      ! local declarations
      integer i,j
      real(mp) ymin(2)
    
      j = 1
      do while (j < n)
         ymin = y(j,:)
         do i=j+1,n
            if (y(i,1) < ymin(1)) then
               y(j,:) = y(i,:)
               y(i,:) = ymin
               ymin = y(j,:)
            end if
         end do
         j = j+1
      end do
      return
   end subroutine dsort2D  
   
   subroutine dsortArray (n, m, y, icol)
      ! This sorts an real(mp) nxm array, y(n,m) in ascending order based on
      ! icol column. The second column corresponds to an index for the first
      use realPrec
      implicit none
      integer , intent(in)   :: n ! Row dimension of y
      integer , intent(in)   :: m ! column dimention of y
      real(mp), intent(inout):: y(n,m) ! The icol column is sorted and the corresponding row moves accordingly.
      integer , intent(in)   :: icol ! which column is used to sort the rows.
      ! local declarations
      integer i,j
      real(mp) rowMin(m)
    
      j = 1
      do while (j < n)
         rowMin = y(j,1:m)
         do i=j+1,n
            if (y(i,icol) < rowMin(icol)) then
               y(j,1:m) = y(i,1:m)
               y(i,1:m) = rowMin
               rowMin = y(j,1:m)
            end if
         end do
         j = j+1
      end do
      return
   end subroutine dsortArray  
   
   subroutine writeEPsurf(whichPP, epSurf)
      ! This save the sigma_rxn surface as a function of f_pp (x) and phi_pp (y)
      use realPrec
      use globalVars, only: nSurfPts, iunit_EPsurf, tOn_pp_min, tOn_pp_max, tOff_pp_min, tOff_pp_max, tOn_pp, tOff_pp
      implicit none
      integer,  intent(in):: whichPP ! which of the pp is manipulated wrt f_pp and phi_pp
      real(mp), intent(in):: epSurf(nSurfPts, nSurfPts,4)
      ! local declarations
      integer i, epI, epJ
      real(mp) xGrid, yGrid
      character(80) string
      string = 'Zone T = "sigma surfaces for p('
      write(string(len_trim(string)+1:),'(i0)') whichPP
      string = trim(string)//') (J/K)"'
      write(iunit_EPsurf, '(a)') trim(string)
      write(iunit_EPsurf, '(2(a,i0),a)') 'I=', nSurfPts, ' J = ', nSurfPts, ' K = 1'
      
      do i=1, nSurfPts**2
         call ijPoint(nSurfPts, i, epI, epJ, xGrid, yGrid)
         tOn_pp(whichPP)  = tOn_pp_min  + xGrid*(tOn_pp_max  - tOn_pp_min )
         tOff_pp(whichPP) = tOff_pp_min + yGrid*(tOff_pp_max - tOff_pp_min)
         write(iunit_EPsurf, *) tOn_pp(whichPP), tOff_pp(whichPP), epSurf(epI, epJ, 1:4)
      end do      
      return
   end subroutine writeEPsurf   

   subroutine ijPoint(nGrid, ithPt, i, j, x, y)
      ! This routine translates the ithPt to the (i,j) grid, and returns x and y values of the two parameters
      ! Note, matrix numbering is row oriented, NOT like Fortran, which is column oriented.
      ! Example, nGrid = 5
      ! 1  2  3  4  5
      ! 6  7  8  9  10
      ! ...
      ! 21 22 23 24 25
      use realPrec
      implicit none
      integer , intent(in) :: nGrid ! Number of grid points in each dimension
      integer , intent(in) :: ithPt ! the serial value of the grid point (i.e., each point has a unique number)
      integer , intent(out):: i     ! The value of the ith row
      integer , intent(out):: j     ! the value of the jth column
      real(mp), intent(out):: x     ! the value of the x variable between 0.0 and 1.0 (i row)
      real(mp), intent(out):: y     ! the value of the y variable between 0.0 and 1.0 (j column)
      ! local declarations
      integer lastJ, modJ
      ! begin
      if (ithPt > nGrid*nGrid) then
         write(*,'(a)') 'Error::ijPoint: ithPt > nGrid^2'
         i = nGrid; j = nGrid; x = 1.0_mp; y = 1.0_mp
         return
      end if      
      modJ = mod(ithPt, nGrid) ! give the j value, unless it's the last column, then gives 0
      lastJ = 0
      if (modJ == 0) lastJ = 1
      
      ! determine the i value
      i = ithPt/nGrid + 1 - lastJ
      ! determine the j value
      j = modJ + lastJ*nGrid
      
      ! calculate the real spacing [0,1]
      x = (real(i) - 1.0_mp)/(real(nGrid) - 1.0_mp)
      y = (real(j) - 1.0_mp)/(real(nGrid) - 1.0_mp)      
      return
   end subroutine ijPoint      
   
