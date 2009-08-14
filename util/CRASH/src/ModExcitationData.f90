module CRASH_ModExcitationData
  use CRASH_ModIonization,ONLY : get_ioniz_potential
  implicit none
  PRIVATE !except

  !public methods
  public:: n_ground  !Principal quantum number as a function of the charge state and the element number
  public:: get_excitation_energy
  public:: get_virial_coeff4_energy
  public:: get_degeneracy
contains
  !=====================================================================================
  integer function n_ground(iZ,nZ)
    integer,intent(in)::iZ,nZ
    !The principal quantum number of the outermost electron in bounded with an ion
    !with I electrons
    integer,parameter :: nGround0_I(100) = (/ &
         1, 1, &                                                                    !  2
         2, 2, 2, 2, 2, 2, 2, 2, &                                                  !  8
         3, 3, 3, 3, 3, 3, 3, 3, &                                                  !  8
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &                    ! 18
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, &                    ! 18
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &                          !
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &                          ! 32
         7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 /)                                !
    integer,parameter :: nGround_I(100) = (/ &
         1, 1, &                                                                    !  2
         2, 2, 2, 2, 2, 2, 2, 2, &                                                  !  8
         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &                    ! 18
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &                    ! 
         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       &                    ! 32
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &                    !
         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       &                    !
         5, 5, 5, 5, 5, 5, 5, 5    /)                                               ! 40
    !-------------------------------------------------------------------------------!
    if(iZ==0)then 
       n_ground = nGround0_I(nZ)
    else
       n_ground = nGround_I(nZ - iZ)
    end if
  end function n_ground
  !======================================================================================
  !Fill in the array ExcitationEnergy_III with tabulated or calculted excitation energies
  !for quantum states defined by n, l of atoms or ions of a particular element defined
  !by their charge state
  subroutine get_excitation_energy(nExcitation, nZ, ExcitationEnergy_III)

    ! From SPECTR-w3
    ! Be, nZ = 4
    !             S  !    P       !    D      !       F       !     G      !  
    real, parameter, dimension(0:4, 1:5, 0:3) :: cExcitationBe_III = reshape((/&
         !Ionization potential = 9.32263, nGround = 2 
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    2.724652,         0.0,         0.0,         0.0,&  !n = 2
           6.455745,    7.288175,    7.693182,         0.0,         0.0,&  !n = 3
           7.949928,    8.220561,    8.366949,         0.0,         0.0,&  !n = 4
           8.484647,    8.611297,    8.680418,         0.0,         0.0,&  !n = 5
         !Ionization potential = 18.21116, nGround = 2 
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    3.956807,         0.0,         0.0,         0.0,&  !n = 2
          10.938703,   11.963904,   12.156649,         0.0,         0.0,&  !n = 3
          14.315214,   14.724362,   14.806191,         0.0,         0.0,&  !n = 4
          15.787526,   15.990240,   16.031031,         0.0,         0.0,&  !n = 5
         !Ionization potential = 153.897, nGround = 1
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         118.590874,  121.923569,         0.0,         0.0,         0.0,&  !n = 2
         139.008589,  139.891357,  140.273228,         0.0,         0.0,&  !n = 3
         145.718613,  146.074448,  146.233148,         0.0,         0.0,&  !n = 4
         148.731429,  148.912446,  148.991796,  148.996755,         0.0,&  !n = 5
         !Ionization potential = 217.713, nGround = 1
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         163.285338,  163.284594,         0.0,         0.0,         0.0,&  !n = 2
         193.527201,  193.526977,  193.530424,         0.0,         0.0,&  !n = 3
         204.111433,  204.111346,  204.112797,  204.113280,         0.0,&  !n = 4
         209.010234,  209.010185,  209.010928,  209.011176,  209.011300 &  !n = 5 
         /), (/5, 5, 4/))

    ! From SPECTR-w3 and
    ! http://physics.nist.gov/PhysRefData/Handbook/Tables/carbontable5.htm
    real, parameter, dimension(0:4, 1:5, 0:5) :: cExcitationC_III = reshape((/&
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    0.002480,         0.0,         0.0,         0.0,&  !n = 2
           7.479966,    8.537097,    9.631092,         0.0,         0.0,&  !n = 3
           9.683165,         0.0,   10.352680,         0.0,         0.0,&  !n = 4
          10.382436,         0.0,   10.686197,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    0.007439,         0.0,         0.0,         0.0,&  !n = 2
          14.448869,   16.331197,   18.045898,         0.0,         0.0,&  !n = 3
          19.494034,   20.149786,   20.844221,   20.950848,         0.0,&  !n = 4
          21.492411,   21.732940,   22.130557,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    6.488266,         0.0,         0.0,         0.0,&  !n = 2
          29.538612,   32.096530,   33.476226,         0.0,         0.0,&  !n = 3
          38.409805,   39.376633,   39.849757,   39.922908,         0.0,&  !n = 4
          41.969887,   42.494836,   42.734249,   43.042350,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    7.994252,         0.0,         0.0,         0.0,&  !n = 2
          37.548611,   39.679981,   40.189981,         0.0,         0.0,&  !n = 3
          49.736504,   50.623983,   50.875175,   50.886829,         0.0,&  !n = 4
          55.200363,   55.651294,   55.779245,   55.785445,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         298.958108,  304.399774,         0.0,         0.0,         0.0,&  !n = 2
         352.060534,  353.530987,  354.260014,         0.0,         0.0,&  !n = 3
         369.911777,  370.509381,  370.809423,  370.826781,         0.0,&  !n = 4
         378.019103,  378.316665,  378.472885,  378.480324,  378.481564,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         367.477261,  367.474025,         0.0,         0.0,         0.0,&  !n = 2
         435.547716,  435.546749,  435.564169,         0.0,         0.0,&  !n = 3
         459.370211,  459.369802,  459.377154,  459.379596,         0.0,&  !n = 4
         470.395666,  470.395455,  470.399212,  470.400464,  470.401096 &  !n = 5
         /), (/5, 5, 6/))

    ! From SPECTR-w3
    real, parameter, dimension(0:4, 1:5, 0:6) :: cExcitationN_III = reshape((/&
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    2.382976,         0.0,         0.0,         0.0,&  !n = 2
          10.325403,   11.844209,   12.971226,         0.0,         0.0,&  !n = 3
          12.847241,   19.148118,   13.661817,         0.0,         0.0,&  !n = 4
          13.614703,   19.629176,   13.980457,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    0.006199,         0.0,         0.0,         0.0,&  !n = 2
          18.462485,   20.409037,   23.196201,         0.0,         0.0,&  !n = 3
          24.367852,   25.065883,   26.028000,         0.0,         0.0,&  !n = 4
          26.558652,         0.0,   27.338513,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    0.021077,         0.0,         0.0,         0.0,&  !n = 2
          27.437700,   30.458699,   33.133410,         0.0,         0.0,&  !n = 3
          37.329159,   38.644755,   39.395975,   39.710895,         0.0,&  !n = 4
          41.374763,   41.898472,   42.395896,   42.495580,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    8.338842,         0.0,         0.0,         0.0,&  !n = 2
          46.776754,   50.142428,   52.079557,         0.0,         0.0,&  !n = 3
          60.455929,   62.448355,   63.343149,   64.045271,         0.0,&  !n = 4
          67.437230,   68.081080,   68.411498,   68.811223,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
                0.0,    9.975768,         0.0,         0.0,         0.0,&  !n = 2
          56.552287,   59.235304,   60.057940,         0.0,         0.0,&  !n = 3
          75.149667,   76.267632,   76.611440,   76.629542,         0.0,&  !n = 4
          83.530874,   84.098473,   84.273291,   84.283210,   84.284449,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         419.796815,  426.294826,         0.0,         0.0,         0.0,&  !n = 2
         494.927512,  496.684367,  497.600611,         0.0,         0.0,&  !n = 3
         520.336831,  521.052219,  521.432851,  521.453928,         0.0,&  !n = 4
         531.913234,  532.272788,  532.458764,  532.474882,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 1
         500.252260,  500.246681,         0.0,         0.0,         0.0,&  !n = 2
         592.926844,  592.925108,  592.957344,         0.0,         0.0,&  !n = 3
         625.358875,  625.358131,  625.371769,  625.376233,         0.0,&  !n = 4
         640.368400,  640.368028,  640.375096,  640.377327,  640.378567 &  !n = 5
         /), (/5, 5, 7/))

    ! Partly from SPECTR-w3, partly from
    ! http://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5.htm and
    ! http://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable6.htm
    real, parameter, dimension(0:4, 2:5, 0:7) :: cExcitationO_III = reshape((/&
                0.0,    0.019837,         0.0,         0.0,         0.0,&  !n = 2
           9.146313,   10.740224,   12.078539,         0.0,         0.0,&  !n = 3
          11.838010,         0.0,   12.754253,         0.0,         0.0,&  !n = 4
          12.661265,         0.0,   13.069173,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
          22.966245,   25.285619,   28.677062,         0.0,         0.0,&  !n = 3
          29.586031,         0.0,         0.0,         0.0,         0.0,&  !n = 4
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 5

                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 2
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 3
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 4
                0.0,         0.0,         0.0,         0.0,         0.0,&  !n = 5

                0.0,    0.047920,         0.0,         0.0,         0.0,&  !n = 2
          44.337985,   48.373670,   52.015085,         0.0,         0.0,&  !n = 3
          60.233997,   61.927001,   63.301366,   63.627444,         0.0,&  !n = 4
          66.873102,   67.775459,   68.442990,   68.501263,         0.0,&  !n = 5

                0.0,   10.174961,         0.0,         0.0,         0.0,&  !n = 2
          67.818358,   72.004808,   74.498006,         0.0,         0.0,&  !n = 3
          89.425082,   91.082006,   91.951755,         0.0,         0.0,&  !n = 4
          98.605987,   99.423290,   99.847936,         0.0,         0.0,&  !n = 5

                0.0,   11.948852,         0.0,         0.0,         0.0,&  !n = 2
          79.353970,   82.587230,   83.642955,         0.0,         0.0,&  !n = 3
         105.688583,  107.038647,  107.478667,  107.503092,         0.0,&  !n = 4
         117.601604,  118.287732,  118.518963,         0.0,         0.0,&  !n = 5

         560.983806,  568.544362,         0.0,         0.0,         0.0,&  !n = 2
         661.929250,  664.018384,  665.104485,         0.0,         0.0,&  !n = 3
         696.337342,  697.113483,  697.597021,  697.561065,         0.0,&  !n = 4
         711.991585,  712.343700,  712.613985,  712.561912,         0.0,&  !n = 5

         653.502789,  653.493738,         0.0,         0.0,         0.0,&  !n = 2
         774.581653,  774.578926,  774.633975,         0.0,         0.0,&  !n = 3
         816.952257,  816.951141,  816.974326,  816.982137,         0.0,&  !n = 4
         836.560728,  836.560232,  836.572134,  836.575978,  836.577962 &  !n = 5
         /), (/5, 4, 8/))

    ! E.B.Saloman, Energy Levels and Observed Spectral Lines of Xenon,
    ! Xe I through Xe LIV
    real, parameter, dimension(0:4, 4:10) :: cExcitationXe0_III = reshape((/&
               0.0,       0.0,       0.0,  11.26270,       0.0,&  !n = 4
               0.0,       0.0,   9.89038,  11.57549,  11.58275,&  !n = 5
           8.31532,   9.58015,  10.97149,  11.74563,       0.0,&  !n = 6
          10.56206,  10.90157,  11.43877,  11.84806,       0.0,&  !n = 7
          11.25833,  11.42555,  11.66993,  11.91442,       0.0,&  !n = 8
          11.57991,  11.66281,  11.80076,  11.95984,       0.0,&  !n = 9
          11.74873,  11.79764,  11.88912,  11.99231,       0.0 &  !n = 10
         /), (/5, 7/))

    real, parameter, dimension(0:4, 4:8) :: cExcitationXe1_III = reshape((/&
               0.0,       0.0,       0.0,  17.22982,       0.0,&  !n = 4
               0.0,  11.26692,  11.82769,  18.60964,  18.77645,&  !n = 5
          11.53901,  13.86046,  16.80076,  19.32030,  19.44976,&  !n = 6
          16.43024,  17.24977,       0.0,       0.0,  19.85531,&  !n = 7
          18.26418,       0.0,       0.0,       0.0,       0.0 &  !n = 8
         /), (/5, 5/))

    real, parameter, dimension(0:3, 4:7) :: cExcitationXe2_III = reshape((/&
               0.0,       0.0,       0.0,  20.62542,&  !n = 4
               0.0,  12.18299,  13.83731,  24.46339,&  !n = 5
          15.06110,  18.19858,  22.60701,       0.0,&  !n = 6
          22.62497,       0.0,       0.0,       0.0 &  !n = 7
         /), (/4, 4/))

    integer :: nExcitation
    integer :: nZ
    real,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: ExcitationEnergy_III
    real,dimension(1:nZ) :: IonizPotential_I

    integer :: iL, iN, iZ, nGround
    !--------------
    ExcitationEnergy_III = 0.0

    call get_ioniz_potential(nZ, IonizPotential_I(1:nZ))

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround+1, nExcitation
          ExcitationEnergy_III(0:iN-1,iN,iZ) = IonizPotential_I(iZ+1) * &
               (1.0 - (real(nGround) / iN)**2.0)

          do iL = 0, iN-1
             select case (nZ)
                case (4)
                   if (iL <= 4 .and. iN <= 5 .and. cExcitationBe_III(iL,iN,iZ) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationBe_III(iL,iN,iZ)

                case (6)
                   if (iL <= 4 .and. iN <= 5 .and. cExcitationC_III(iL,iN,iZ) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationC_III(iL,iN,iZ)

                case (7)
                   if (iL <= 4 .and. iN <= 5 .and. cExcitationN_III(iL,iN,iZ) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationN_III(iL,iN,iZ)

                case (8)
                   if (iL <= 4 .and. iN >= 2 .and. iN <= 5 .and. cExcitationO_III(iL,iN,iZ) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationO_III(iL,iN,iZ)

                case (54)
                   if (iZ == 0 .and. iL <= 4 .and. iN >= 4 .and. iN <= 10 .and. &
                        cExcitationXe0_III(iL,iN) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationXe0_III(iL,iN)
                   if (iZ == 1 .and. iL <= 4 .and. iN >= 4 .and. iN <= 8 .and. &
                        cExcitationXe1_III(iL,iN) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationXe1_III(iL,iN)
                   if (iZ == 2 .and. iL <= 3 .and. iN >= 4 .and. iN <= 7 .and. &
                        cExcitationXe2_III(iL,iN) /= 0.0)&
                        ExcitationEnergy_III(iL,iN,iZ) = cExcitationXe2_III(iL,iN)

             end select
          end do
       end do
    end do
  end subroutine get_excitation_energy
  !======================================================================================
  subroutine get_virial_coeff4_energy(nExcitation, nZ, VirialCoeff4Energy_III)
    integer :: nExcitation
    integer :: nZ
    real,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: VirialCoeff4Energy_III
    real,dimension(1:nZ) :: IonizPotential_I

    integer :: iZ, iN, iL, nGround
    !--------------
    VirialCoeff4Energy_III = 0.0

    call get_ioniz_potential(nZ, IonizPotential_I(1:nZ))

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround, nExcitation
          do iL = 0, iN-1
             VirialCoeff4Energy_III(iL,iN,iZ) = IonizPotential_I(iZ+1) * &
                  real(nGround)**2 * real(iN - iL)**2 / real(iZ+1)**2
          end do
       end do
    end do
  end subroutine get_virial_coeff4_energy
  !======================================================================================
  subroutine get_degeneracy(nExcitation, nZ, Degeneracy_III)
    integer :: nExcitation
    integer :: nZ
    integer,dimension(0:nExcitation-1,nExcitation,0:nZ-1),intent(out) :: Degeneracy_III
    integer :: iMix

    integer :: iZ, iN, iL, nGround


    !Degeneracies of ground states for positive ions of the first 10 elements.
    !The data for the first 3 columns are taken from the book 
    !"Allen' Astrophysical Quantities, Edn iV" p. 33
    !by Allen, Clabon W. Editor: Cox, Arthur N.
    !Publisher: Springer (2000) (g_0 for Y=1, Y=2, Y=3

    !First 10 elements - full ionizations
    integer, parameter, dimension(10,10) :: cDegeneracy10_II = reshape(  (/   &
         !   I  !  II  !  III !  IV  !   V  !  VI  !  VII ! VIII !  IX  !   X  !
         2 ,     1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 1 - H
         1 ,     2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 2 - He
         2 ,     1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 3 - Li
         1 ,     2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 4 - Be
         6 ,     1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1 ,    1     ,&  ! 5 - B
         9 ,     6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1 ,    1     ,&  ! 6 - C
         4 ,     9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1 ,    1     ,&  ! 7 - N
         9 ,     4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1 ,    1     ,&  ! 8 - O
         6 ,     9 ,    4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2 ,    1     ,&  ! 9 - F
         1 ,     6 ,    9 ,    4 ,    9 ,    6 ,    1 ,    2 ,    1 ,    2   /),&  !10 - Ne
                                                                                  (/10,10/))
    !--------------
    Degeneracy_III = 0

    do iZ = 0, nZ-1
       nGround = n_ground(iZ, nZ)

       do iN = nGround, nExcitation
          do iL = 0, iN-1
             if (iN == nGround .and. nZ <= 10) then
                if (iL == 0)&
                     Degeneracy_III(iZ,iN,iL) = cDegeneracy10_II(iZ+1,nZ)
             else
                Degeneracy_III(iZ,iN,iL) = 2 * (2*iL + 1)

             end if
          end do
       end do
    end do
  end subroutine get_degeneracy
  !======================================================================================
end module CRASH_ModExcitationData
