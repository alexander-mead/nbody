MODULE constants

  ! TODO: Can I do calculations in this file? For example, for non-fundamental constants like sigma_T
  IMPLICIT NONE

  !!

  ! Mathematical constants
  REAL, PARAMETER :: pi=3.14159265359 ! pi  
  REAL, PARAMETER :: em=0.5772156649  ! Euler–Mascheroni
  REAL, PARAMETER :: zero=0.          ! zero
  REAL, PARAMETER :: one=1.           ! one
  REAL, PARAMETER :: rad2deg=180./pi  ! radians-to-degrees conversion

  ! Mathematical derived
  REAL, PARAMETER :: twopi=2.*pi        ! 2pi or tau
  REAL, PARAMETER :: deg2rad=1./rad2deg ! degress-to-radians conversion

  !!

  !!

  ! Physical constants
  REAL, PARAMETER :: kB=1.38064852e-23        ! Boltzmann constant [m^2 kg s^-2 K^-1]  
  REAL, PARAMETER :: mp=1.6726219e-27         ! proton mass [kg]
  REAL, PARAMETER :: me=9.10938356e-31        ! electron mass [kg]
  REAL, PARAMETER :: bigG=6.67408e-11         ! Gravitational constant [kg^-1 m^3 s^-2]
  REAL, PARAMETER :: eV=1.60218e-19           ! electronvolt [kg m^2 s^-2]  
  REAL, PARAMETER :: SBconst=5.670367e-8      ! Steffan-Boltzmann constant [kg s^-3 K^-4]
  REAL, PARAMETER :: c_light=2.99792458e8     ! speed of light [m/s]
  REAL, PARAMETER :: sigma_T=6.6524587158e-29 ! Thompson-scatter cross section [m^2]
  REAL, PARAMETER :: h_Planck=6.62607004e-34  ! Planck constant [kg m^2/s]

  ! Derived physical constants
  REAL, PARAMETER :: me_eV=me*(c_light)**2/eV ! Electron rest mass in eV [~511 keV]
  REAL, PARAMETER :: mp_eV=mp*(c_light)**2/eV ! Proton rest mass in eV [~938 MeV]

  !!

  !!
  
  ! Astronomy constants
  REAL, PARAMETER :: Htime=9.7776                   ! Hubble time (1/H0) [Gyrs/h] ! TODO: relate this to fundamental constants
  REAL, PARAMETER :: H0=3.243e-18                   ! Hubble parameter, H0 [s]
  REAL, PARAMETER :: critical_density=2.7755e11     ! Universal critical density at (equal to 3*H0^2 / 8piG) [(M_sun/h)/(Mpc/h)^3] ! TODO: relate this to fundamental constants
  REAL, PARAMETER :: dc0=(3./20.)*(12.*pi)**(2./3.) ! Einstein-de Sitter linear collapse density ~1.686
  REAL, PARAMETER :: Dv0=18.*pi**2                  ! Einsten-de Sitter virialised collapse threshold ~178
  REAL, PARAMETER :: Msun=1.989e30                  ! Solar mass [kg]
  REAL, PARAMETER :: Mpc=3.086e22                   ! Mpc [m]
  REAL, PARAMETER :: SI_to_Jansky=1e26              ! [W m^-2 Hz^-1 / Jy]

  ! Derived astronomy constants
  REAL, PARAMETER :: Hdist=c_light/1e5              ! Hubble parameter distance (c/H0) [Mpc/h]
  REAL, PARAMETER :: yfac=sigma_T/(me*c_light**2)   ! sigma_T/m_e*c^2 [kg^-1 s^2], prefactor of Compton-y integral over *pressure*

  !!
  
END MODULE constants