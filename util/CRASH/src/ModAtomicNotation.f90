!^CFG COPYRIGHT UM
module CRASH_ModAtomicNotation
 !Chemistry and spectroscopy
  !\
  ! The Mendeleev table
  !/
  character(LEN=2),parameter:: NameElement_I(-1:54)=(/&
     'Ay' ,'Pl',                                                                               & !-1:0
     'H_' ,                                                                               'He',& ! 1:2
     'Li','Be',                                                  'B_','C_','N_','O_','F_','Ne',& ! 3:10
     'Na','Mg',                                                  'Al','Si','P_','S_','Cl','Ar',& !11:18
     'K_','Ca','Sc','Ti','V_','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',& !19:36
     'Rb','Sr','Y_','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I_','Xe'/) !37:54
  !Orbital quantum number: transform a numerical value to a symbol:
  character(LEN=1),parameter::TypeL_I(0:9) = (/'s','p','d','f','g','h','i','k','l','m'/)
contains
    integer function l_orbital(TypeL)
    character(LEN=1),intent(in)::TypeL
    !--------------------------------!
    l_orbital=-1
    select case(TypeL)
    case('s','S')
       l_orbital = 0
    case('p','P')
       l_orbital = 1
    case('d','D')
       l_orbital = 2
    case('f','F')
       l_orbital = 3
    case('g','G')
       l_orbital = 4
    case('h')
       l_orbital = 5
    case('i')
       l_orbital = 6
    case('k')
       l_orbital = 7
    case('l')
       l_orbital = 8
    case('m')
       l_orbital = 9
    case default
       call CON_stop('The spectroscopy symbol '//TypeL//' is not implemented')
    end select
  end function l_orbital
  !====================================================================

end module CRASH_ModAtomicNotation
