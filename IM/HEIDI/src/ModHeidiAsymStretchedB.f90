module get_gradB_components

  use ModHeidiMain,      ONLY: DipoleFactor
  use ModHeidiInput,     ONLY: StretchingFactorA, StretchingFactorB

contains
  !========================================================================
  subroutine get_dBdx(x, y, z, dBdx)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: dBdx
    ! Local variables used for optimization
    real :: t1,   t2,   t3,   t4,   t5,   t6,   t7,   t8,   t9,   t10
    real :: t13,  t14,  t15,  t18,  t19,  t23,  t24,  t25,  t28,  t29
    real :: t32,  t33,  t34,  t37,  t38,  t39,  t42,  t43,  t46,  t47  
    real :: t51,  t52,  t56,  t60,  t61,  t65,  t67,  t69,  t72,  t74
    real :: t77,  t80,  t96,  t99,  t104, t112, t115, t116, t117, t118 
    real :: t120, t126, t127, t128, t130, t132, t138, t139, t146, t149
    real :: t155, t156, t159, t161, t165, t187, t203, t211, t214, t249
    real :: t252, t255, t279
    real :: a, b

    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB
    
    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/

    
    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + t2
    t4 = t3 ** 2
    t5 = t4 * t3
    t6 = z ** 2
    t7 = t1 ** 2
    t8 = t7 * t6
    t9 = a ** 2
    t10 = t2 * t9
    t13 = t1 * t6
    t14 = t2 ** 2
    t15 = t14 * t9
    t18 = b ** 2
    t19 = t18 * t6
    t23 = t7 * t1
    t24 = t9 * t23
    t25 = t2 * t18
    t28 = t7 ** 2
    t29 = a * t28
    t32 = sqrt(0.1D1 / t1 * t3)
    t33 = t18 * b
    t34 = t33 * t32
    t37 = t9 * a
    t38 = t37 * t28
    t39 = b * t32
    t42 = t37 * t23
    t43 = t2 * t39
    t46 = t32 * t6
    t47 = a * b
    t51 = t32 * a
    t52 = t2 * b
    t56 = t14 * b
    t60 = a * t46
    t61 = t7 * b
    t65 = t18 ** 2
    t67 = t23 * t6
    t69 = t14 * t2
    t72 = t9 ** 2
    t74 = t6 ** 2
    
    t77 = -dble(8 * t10 * t8) - dble(4 * t15 * t13) - 0.4D1 * t2 * t7 * t19 + &
         dble(6 * t25 * t24) + 0.4D1 * t34 * t29 + 0.4D1 * t39 * t38 + 0.4D1 * &
         t43 * t42 - 0.8D1 * t23 * t47 * t46 + 0.4D1 * t23 * t52 * t51 + 0.4D1 * &
         t7 * t56 * t51 - 0.8D1 * t2 * t61 * t60 + t65 * t28 + 0.9D1 * t67 +  &
         0.5D1 * t69 * t6 + t72 * t28 + 0.4D1 * t7 * t74

    t80 = t7 * t14
    t96 = t9 * t7
    t99 = t9 * t1
    t104 = t14 * t18
    t112 = t2 * t74
    t115 = t14 ** 2
    
    t116 = 0.4D1 * t14 * t74 + t80 + 0.2D1 * t1 * t69 + 0.23D2 * t2 * dble(t8) + &
         0.19D2 * t14 * dble(t13) - 0.4D1 * t9 * t67 - 0.4D1 * t23 * t19 + 0.6D1  &
         * t18 * t9 * t28 + 0.2D1 * t2 * dble(t24) + 0.4D1 * t14 * t96 + 0.2D1 * &
         t69 * t99 + 0.2D1 * t23 * dble(t25) + 0.2D1 * t7 * t104 + 0.2D1 * t2 * &
         t72 * t23 + t14 * t72 * t7 + 0.8D1 * t1 * t112 + t115
    
    t117 = t77 + t116
    t118 = t117 * t5
    t120 = a * t7
    
    t126 = t96 + t2 * t99 + 0.2D1 * t39 * t120 + t18 * t7 + t1 * t2 + t14 + dble(t13) + t2 * t6
    t127 = t126 ** 2
    t128 = t127 ** 2
    t130 = 0.1D1 / t128 / t126
    t132 = sqrt(t130 * t118)
    t138 = t1 * x
    t139 = t7 * t138
    t146 = t138 * t6
    t149 = x * t6
    t155 = t7 * x
    t156 = t9 * t155
    t159 = 0.1D1 / t32
    t161 = a * t159 * t6
    t165 = 0.1D1 / x - 0.1D1 / t138 * t3
    t187 = 0.2D1 * t165 * b * t159
    t203 = b * t159 * a
    t211 = 0.32D2 * t34 * a * t139 + 0.32D2 * t39 * t37 * t139 - 0.32D2 * dble(t10) * &
         t146 - 0.8D1 * dble(t15) * t149 - 0.16D2 * t2 * t138 * t19 + 0.36D2 *        &
         dble(t25) * t156 - 0.8D1 * t165 * t2 * t61 * t161 + 0.24D2 * t43 * t37 *     &
         t155 - 0.48D2 * t155 * t47 * t46 + 0.24D2 * t155 * t52 * t51 + 0.16D2 *      &
         t138 * t56 * t51 + 0.4D1 * t165 * t33 * t159 * t29 + 0.2D1 * t187 * t38 -    &
         0.32D2 * t2 * t138 * b * t60 + 0.4D1 * t165 * t52 * t159 * t42 - 0.8D1 *     &
         t165 * t23 * b * t161 + 0.4D1 * t165 * t23 * t2 * t203 + 0.4D1 * t165 * t80 * t203
    
    t214 = t155 * t6
    t249 = t9 * t138
    t252 = t9 * x
    
    t255 = 0.8D1 * t65 * t139 + 0.54D2 * t214 + 0.8D1 * t72 * t139 + 0.16D2 * t138 * &
         t74 + 0.4D1 * t138 * t14 + 0.4D1 * x * t69 + 0.12D2 * t2 * t72 * t155 +     &
         0.48D2 * t18 * t9 * t139 + 0.16D2 * x * t112 + 0.12D2 * t155 * dble(t25) +  &
         0.92D2 * t2 * t146 + 0.12D2 * t2 * t156 + 0.4D1 * t14 * t72 * t138 -        &
         0.24D2 * t155 * t19 - 0.24D2 * t9 * t214 + 0.8D1 * t138 * t104 + 0.38D2 *   &
         t14 * t149 + 0.16D2 * t14 * t249 + 0.4D1 * t69 * t252
    
    t279 = (0.6D1 * x * t130 * t117 * t4 + t130 * (t211 + t255) * t5 - 0.5D1 *  &
         (0.4D1 * t249 + 0.2D1 * t2 * t252 + 0.8D1 * t39 * a * t138 + t187 *    &
         t120 + 0.4D1 * t18 * t138 + 0.2D1 * x * t2 + 0.2D1 * t149) / t128 /    &
         t127 * t118) / t132 / 0.2D1




    dBdx = abs(DipoleFactor) * t279

  end subroutine get_dBdx
  !========================================================================
  subroutine get_dBdy(x, y, z, dBdy)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: dBdy
    ! Local variables used for optimization
    real :: t1,   t2,    t3,   t4,   t5,   t6,   t7,  t8,  t9,   t13
    real :: t14,  t18,   t19, t23,  t24,  t25,  t28,  t32, t33,  t37
    real :: t39,  t42,   t43, t44,  t48,  t56,  t65,  t67, t69,  t72
    real :: t74,  t77,   t96, t99,  t107, t110, t115, t116, t117, t118
    real :: t126, t127, t128, t130, t132, t138, t140, t141, t145, t159
    real :: t168, t172, t175, t180, t191, t196, t228, t250
    real :: a, b

    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB

    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/

    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + dble(t2)
    t4 = t3 ** 2
    t5 = t4 * t3
    t6 = z ** 2
    t7 = t1 ** 2
    t8 = t7 * t6
    t9 = a ** 2
    t13 = t1 * t6
    t14 = t2 ** 2
    t18 = b ** 2
    t19 = t18 * t6
    t23 = t7 * t1
    t24 = dble(t9) * t23
    t25 = dble(t2) * t18
    t28 = t7 ** 2
    t32 = sqrt(0.1D1 / t1 * t3)
    t33 = t18 * b
    t37 = dble(t9) * a
    t39 = b * t32
    t42 = t32 * t6
    t43 = a * t42
    t44 = t7 * b
    t48 = t37 * t23
    t56 = t32 * a
    t65 = t18 ** 2
    t67 = t23 * t6
    t69 = dble(t14) * dble(t2)
    t72 = t9 ** 2
    t74 = t6 ** 2
    
    t77 = -dble(8 * t2 * t9 * t8) - dble(4 * t14 * t9 * t13) - 0.4D1 * dble(t2) * t7 * t19 + &
         dble(6 * t25 * t24) + 0.4D1 * t33 * t32 * a * t28 + 0.4D1 * t39 * t37 * t28 -       &
         0.8D1 * dble(t2) * t44 * t43 + 0.4D1 * dble(t2) * t39 * t48 - 0.8D1 * t23 * b * a * &
         t42 + 0.4D1 * t23 * dble(t2) * b * t56 + 0.4D1 * t7 * dble(t14) * b * t56 + t65 *  &
         t28 + 0.9D1 * t67 + 0.5D1 * t69 * t6 + t72 * t28 + 0.4D1 * t7 * t74

    t96 = dble(t9) * t7
    t99 = dble(t9) * t1
    t107 = t72 * t23
    t110 = t72 * t7
    t115 = t14 ** 2
    
    t116 = 0.4D1 * dble(t14) * t74 + t7 * dble(t14) + 0.2D1 * t1 * t69 + dble(23 * t2 * t8) + &
         dble(19 * t14 * t13) - 0.4D1 * dble(t9) * t67 - 0.4D1 * t23 * t19 + 0.6D1 * t18 *    &
         dble(t9) * t28 + dble(2 * t2 * t24) + 0.4D1 * dble(t14) * t96 + 0.2D1 * t69 * t99 +  &
         0.2D1 * t23 * dble(t25) + 0.2D1 * t7 * dble(t14) * t18 + 0.2D1 * dble(t2) * t107 +   &
         dble(t14) * t110 + 0.8D1 * t1 * dble(t2) * t74 + dble(t115)
    
    t117 = t77 + t116
    t118 = t117 * t5
    
    t126 = t96 + dble(t2) * t99 + 0.2D1 * t39 * a * t7 + t18 * t7 + t1 * &
         dble(t2) + dble(t14) + dble(t13) + dble(t2) * t6
    
    t127 = t126 ** 2
    t128 = t127 ** 2
    t130 = 0.1D1 / t128 / t126
    t132 = sqrt(t130 * t118)
    t138 = 0.1D1 / t32
    t140 = a * t138 * t6
    t141 = y * t44
    t145 = dble(t2) * y
    t159 = t7 * t145 * b
    t168 = b * t138
    t172 = t138 * a
    t175 = dble(t14) * y
    t180 = y * t168
    
    t191 = -0.8D1 * t141 * t140 - 0.8D1 * t145 * t1 * b * t140 - 0.16D2 * t141 * &
         t43 + 0.8D1 * y * t39 * t48 + 0.8D1 * t23 * y * b * t56 + 0.16D2 * t159 * &
         t56 + 0.4D1 * y * t33 * t138 * a * t23 + 0.4D1 * t145 * t168 * t37 * t7 + &
         0.4D1 * t159 * t172 + 0.4D1 * t1 * t175 * b * t172 + 0.4D1 * t180 * t48 + &
         0.30D2 * t175 * t6 + 0.16D2 * t145 * t74 + 0.4D1 * t7 * t145 + 0.12D2 * t1 * t175
    
    t196 = y * t18
    
    t228 = 0.4D1 * y * t107 + 0.76D2 * t145 * dble(t13) + 0.4D1 * t23 * t196 + 0.4D1 *&
         t145 * t110 + 0.46D2 * y * dble(t8) + 0.16D2 * t145 * t96 + 0.12D2 * t175 *  &
         t99 + 0.16D2 * t1 * y * t74 + 0.4D1 * y * dble(t24) + 0.8D1 * t7 * t145 *    &
         t18 + 0.8D1 * dble(t14) * t145 - 0.16D2 * y * dble(t9) * dble(t8) - 0.16D2 * &
         t145 * dble(t9) * dble(t13) - 0.8D1 * y * t7 * t19 + 0.12D2 * t196 * dble(t24)
    
    t250 = (0.6D1 * y * t130 * t117 * t4 + t130 * (t191 + t228) * t5 - 0.5D1 *&
         (0.2D1 * y * t99 + 0.2D1 * t180 * a * t1 + 0.2D1 * t1 * y + 0.4D1 * &
         t145 + 0.2D1 * y * t6) / t128 / t127 * t118) / t132 / 0.2D1

    dBdy = abs(DipoleFactor) * t250

  end subroutine get_dBdy
  !========================================================================
  subroutine get_dBdz(x, y, z, dBdz)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: dBdz
    ! Local variables used for optimization
    real :: t1,   t2,   t3,   t4,   t5,   t6,   t7,   t8,   t9,   t10
    real :: t13,  t14,  t15,  t18,  t19,  t20,  t23,  t24,  t25,  t28
    real :: t32,  t37,  t39,  t42,  t45,  t49,  t52,  t65,  t67,  t69
    real :: t72,  t74,  t77,  t83,  t102, t105, t116, t118, t126, t127
    real :: t128, t130, t132, t134, t137, t140, t143, t151, t156, t169
    real :: t181
   
    real :: a, b


    real :: t11, t22
    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB

    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/

    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + t2
    t4 = t3 ** 2
    t5 = t4 * t3
    t6 = z ** 2
    t7 = t1 ** 2
    t8 = t7 * t6
    t9 = a ** 2
    t10 = t2 * t9
    t13 = t1 * t6
    t14 = t2 ** 2
    t15 = t14 * t9
    t18 = b ** 2
    t19 = t18 * t6
    t20 = t2 * t7
    t23 = t7 * t1
    t24 = t9 * t23
    t25 = t2 * t18
    t28 = t7 ** 2
    t32 = sqrt(0.1D1 / t1 * t3)
    t37 = t9 * a
    t39 = b * t32
    t42 = t32 * t6
    t45 = t2 * t7 * b
    t49 = t23 * b * a
    t52 = t32 * a
    t65 = t18 ** 2
    t67 = t23 * t6
    t69 = t14 * t2
    t72 = t9 ** 2
    t74 = t6 ** 2
    
    t77 = -dble(8 * t10 * t8) - dble(4 * t15 * t13) - dble(4 * t20 * t19) + dble(6 * t25 * t24) + &
         0.4D1 * t18 * b * t32 * a * t28 + 0.4D1 * t39 * t37 * t28 - 0.8D1 * t45 * a * t42 -      &
         0.8D1 * t49 * t42 + 0.4D1 * t23 * t2 * b * t52 + 0.4D1 * t7 * t14 * b * t52 + 0.4D1 *    &
         t2 * t39 * t37 * t23 + t65 * t28 + 0.9D1 * t67 + 0.5D1 * t69 * t6 + t72 * t28 + 0.4D1 * t7 * t74

    t83 = t14 ** 2
    t102 = t9 * t7
    t105 = t9 * t1
    
    t116 = 0.4D1 * t14 * t74 + t7 * t14 + 0.2D1 * t1 * t69 + t83 + 0.2D1 * t23 * dble(t25) + &
         0.23D2 * t2 * dble(t8) + 0.19D2 * t14 * dble(t13) - 0.4D1 * t9 * t67 - 0.4D1 *      &
         t23 * dble(t19) + 0.6D1 * t18 * t9 * t28 + 0.8D1 * t1 * t2 * t74 + 0.2D1 * t2 *     &
         dble(t24) + 0.4D1 * t14 * t102 + 0.2D1 * t69 * t105 + 0.2D1 * t7 * t14 * t18 +      &
         0.2D1 * t2 * t72 * t23 + t14 * t72 * t7
    
    t118 = (t77 + t116) * t5
    t126 = t102 + t2 * t105 + 0.2D1 * t39 * a * t7 + t18 * t7 + t1 * t2 + t14 + dble(t13) + t2 * t6
    
    t127 = t126 ** 2
    t128 = t127 ** 2
    t130 = 0.1D1 / t128 / t126
    t132 = sqrt(t130 * t118)
    t134 = t7 * z
    t137 = t1 * z
    t140 = t18 * z
    t143 = t32 * z
    t151 = t23 * z
    t156 = t6 * z
    
    t169 = -0.16D2 * dble(t10) * t134 - 0.8D1 * dble(t15) * t137 - 0.8D1 * dble(t20) * t140 - &
         0.16D2 * t45 * a * t143 + 0.46D2 * t2 * t134 + 0.38D2 * t14 * t137 - 0.8D1 * t9 *    &
         t151 - 0.8D1 * t23 * t140 + 0.32D2 * t1 * t2 * t156 + 0.18D2 * t151 + 0.10D2 * t69 * &
         z + 0.16D2 * t7 * t156 + 0.16D2 * t14 * t156 - 0.16D2 * t49 * t143
    
    t181 = (t130 * t169 * t5 - 0.10D2 * (t137 + t2 * z) / t128 / t127 * t118) / t132 / 0.2D1

   dBdz = abs(DipoleFactor) * t181




  end subroutine get_dBdz

  !========================================================================
  subroutine get_GradB2CrossB_x(x, y, z, Vx)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: Vx
    ! Local variables used for optimization
    real :: t1,   t2,   t3,   t4,   t5,    t6,   t7,   t8,   t9,   t12
    real :: t13,  t14,  t17,  t18,  t19,   t22,  t23,  t24,  t27,  t29
    real :: t31,  t32,  t36,  t38,  t41,   t45,  t47,  t50,  t59,  t60
    real :: t61,  t64,  t66,  t68,  t71,   t73,  t76,  t86,  t89,  t103
    real :: t106, t114, t115, t116, t125,  t126, t127, t129, t133, t137
    real :: t144, t155, t159, t164, t165,  t172, t175, t180, t184, t191
    real :: t224, t229, t231, t250, t254,  t255, t258, t259, t260, t262
    real :: t264, t268, t271, t274, t278,  t281, t299, t313
    real :: a, b
    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB
    
    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/
   
    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + t2
    t4 = t3 ** 2
    t5 = z ** 2
    t6 = t1 ** 2
    t7 = t6 * t5
    t8 = a ** 2
    t9 = t2 * t8
    t12 = t1 * t5
    t13 = t2 ** 2
    t14 = t13 * t8
    t17 = b ** 2
    t18 = t17 * t5
    t19 = t2 * t6
    t22 = t6 * t1
    t23 = t8 * t22
    t24 = t2 * t17
    t27 = t6 ** 2
    t29 = 0.1D1 / t1
    t31 = sqrt(t29 * t3)
    t32 = t17 * b
    t36 = t8 * a
    t38 = b * t31
    t41 = t36 * t22
    t45 = t31 * t5
    t47 = t22 * b * a
    t50 = t31 * a
    t59 = a * t45
    t60 = t6 * b
    t61 = t2 * t60
    t64 = t17 ** 2
    t66 = t22 * t5
    t68 = t13 * t2
    t71 = t8 ** 2
    t73 = t5 ** 2

    t76 = -dble(8 * t9 * t7) - dble(4 * t14 * t12) - dble(4 * t19 * t18) + dble(6 * t24 * t23) + &
         0.4D1 * t32 * t31 * a * t27 + 0.4D1 * t38 * t36 * t27 + 0.4D1 * t2 * t38 * t41 - 0.8D1 *&
         t47 * t45 + 0.4D1 * t22 * t2 * b * t50 + 0.4D1 * t6 * t13 * b * t50 - 0.8D1 * t61 *     &
         t59 + t64 * t27 + 0.9D1 * t66 + 0.5D1 * t68 * t5 + t71 * t27 + 0.4D1 * t6 * t73
    
    t86 = t71 * t22
    t89 = t71 * t6
    t103 = t8 * t6
    t106 = t8 * t1
    t114 = t13 ** 2
    
    t115 = 0.4D1 * t13 * t73 + t6 * t13 + 0.2D1 * t1 * t68 + 0.23D2 * t2 * dble(t7) + &
         0.19D2 * t13 * dble(t12) + 0.2D1 * t2 * t86 + t13 * t89 + 0.8D1 * t1 * t2 *  &
         t73 - 0.4D1 * t8 * t66 - 0.4D1 * t22 * dble(t18) + 0.6D1 * t17 * t8 * t27 +  &
         0.2D1 * t2 * dble(t23) + 0.4D1 * t13 * t103 + 0.2D1 * t68 * t106 + 0.2D1 *   &
         t22 * dble(t24) + 0.2D1 * t6 * t13 * t17 + t114
    
    t116 = t76 + t115
    t125 = t103 + t2 * t106 + 0.2D1 * t38 * a * t6 + t17 * t6 + t1 * t2 + t13 + dble(t12) + t2 * t5
    
    t126 = t125 ** 2
    t127 = t126 ** 2
    t129 = 0.1D1 / t127 / t125
    t133 = t4 * t3
    t137 = t2 * y
    t144 = y * t17
    t155 = t6 * t137 * b
    t159 = 0.1D1 / t31
    t164 = b * t159
    t165 = y * t164
    t172 = t159 * a
    t175 = t13 * y
    t180 = y * t60
    t184 = a * t159 * t5

    t191 = -0.16D2 * y * t8 * dble(t7) - 0.16D2 * t137 * t8 * dble(t12) - 0.8D1 * y * t6 * &
         dble(t18) + 0.12D2 * t144 * dble(t23) + 0.8D1 * y * t38 * t41 + 0.8D1 * t22 * y * &
         b * t50 + 0.16D2 * t155 * t50 + 0.4D1 * y * t32 * t159 * a * t22 + 0.4D1 * t165 * &
         t41 + 0.4D1 * t137 * t164 * t36 * t6 + 0.4D1 * t155 * t172 + 0.4D1 * t1 * t175 *  &
         b * t172 - 0.16D2 * t180 * t59 - 0.8D1 * t137 * t1 * b * t184 - 0.8D1 * t180 * t184
    
    t224 = 0.30D2 * t175 * t5 + 0.16D2 * t137 * t73 + 0.4D1 * t6 * t137 + 0.12D2 * t1 *   &
         t175 + 0.16D2 * t1 * y * t73 + 0.12D2 * t175 * t106 + 0.4D1 * y * dble(t23) +    &
         0.4D1 * y * t86 + 0.4D1 * t22 * t144 + 0.4D1 * t137 * t89 + 0.8D1 * t6 * t137 *  &
         t17 + 0.16D2 * t137 * t103 + 0.46D2 * y * dble(t7) + 0.76D2 * t137 * dble(t12) + &
         0.8D1 * t13 * t137
    
    t229 = t116 * t133
    t231 = 0.1D1 / t127 / t126
    t250 = sqrt(0.1D1 + t29 * t2)
    t254 = (a + 0.1D1 / t250 * b) ** 2
    t255 = t254 * t1
    t258 = t255 + t2 + t5
    t259 = t258 ** 2
    t260 = sqrt(t258)
    t262 = 0.1D1 / t260 / t259
    t264 = t31 * z
    t268 = t6 * z
    t271 = t1 * z
    t274 = t5 * z
    t278 = t22 * z
    t281 = t17 * z

    t299 = -0.16D2 * t61 * a * t264 + 0.46D2 * t2 * t268 + 0.38D2 * t13 * t271 + 0.32D2 * t1 * &
         t2 * t274 - 0.8D1 * t8 * t278 - 0.8D1 * t22 * t281 + 0.18D2 * t278 + 0.10D2 * t68 *   &
         z + 0.16D2 * t6 * t274 + 0.16D2 * t13 * t274 - 0.16D2 * dble(t9) * t268 - 0.8D1 *     &
         dble(t14) * t271 - 0.8D1 * dble(t19) * t281 - 0.16D2 * t47 * t264
    
    t313 = -t262 * (-0.2D1 * t5 + t255 + t2) * (0.3D1 * y * t129 * t116 * t4 + t129 * &
         (t191 + t224) * t133 / 0.2D1 - 0.5D1 / 0.2D1 * (0.2D1 * y * t106 + 0.2D1 *   &
         t165 * a * t1 + 0.2D1 * t1 * y + 0.4D1 * t137 + 0.2D1 * y * t5) * t231 * t229) -&
         0.3D1 * t262 * y * z * (t129 * t299 * t133 / 0.2D1 - 0.5D1 * (t271 + t2 * z) * t231 * t229)

    Vx = t313


  end subroutine get_GradB2CrossB_x
  !========================================================================
  subroutine get_GradB2CrossB_y(x, y, z, Vy)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: Vy
    ! Local variables used for optimization
    real :: t1,   t2,   t3,   t4,   t5,   t6,   t8,   t9,   t11,  t12
    real :: t13,  t16,  t17,  t18,  t21,  t23,  t24,  t27,  t28,  t33
    real :: t36,  t39,  t42,  t43,  t49,  t52,  t55,  t58,  t60,  t61
    real :: t63,  t64,  t69,  t71,  t72,  t73,  t75,  t78,  t79,  t82
    real :: t87,  t90,  t91,  t94,  t95,  t96,  t97,  t100, t101, t104
    real :: t105, t110, t111, t115, t119, t121, t125, t127, t130, t133
    real :: t136, t158, t166, t167, t168, t169, t171, t181, t185, t186
    real :: t187, t188, t189, t191, t199, t201, t203, t206, t224, t232
    real :: t239, t242, t248, t249, t269, t272, t275, t290, t298, t316
    real :: t342
    real :: a, b
    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB
    
    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/
    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + t2
    t4 = t3 ** 2
    t5 = t4 * t3
    t6 = 0.1D1 / t1
    t8 = sqrt(t6 * t3)
    t9 = t8 * z
    t11 = t1 ** 2
    t12 = t11 * b
    t13 = t2 * t12
    t16 = a * b
    t17 = t11 * t1
    t18 = t17 * t16
    t21 = t17 * z
    t23 = t2 ** 2
    t24 = t23 * t2
    t27 = z ** 2
    t28 = t27 * z
    t33 = t11 * z
    t36 = t1 * z
    t39 = a ** 2
    t42 = b ** 2
    t43 = t42 * z
    t49 = t2 * t39
    t52 = t23 * t39
    t55 = t2 * t11
    
    t58 = -0.16D2 * t13 * a * t9 - 0.16D2 * t18 * t9 + 0.18D2 * t21 + 0.10D2 * t24 * z + &
         0.16D2 * t11 * t28 + 0.16D2 * t23 * t28 + 0.46D2 * t2 * t33 + 0.38D2 * t23 *    &
         t36 - 0.8D1 * t39 * t21 - 0.8D1 * t17 * t43 + 0.32D2 * t1 * t2 * t28 - 0.16D2 * &
         t49 * t33 - 0.8D1 * t52 * t36 - 0.8D1 * t55 * t43
    
    t60 = t39 * t11
    t61 = t39 * t1
    t63 = a * t11
    t64 = b * t8
    t69 = t1 * t27
    
    t71 = t60 + t2 * t61 + 0.2D1 * t64 * t63 + t42 * t11 + t1 * t2 + t23 + t69 + t2 * t27
    t72 = t71 ** 2
    t73 = t72 ** 2
    t75 = 0.1D1 / t73 / t71
    t78 = t8 * t27
    t79 = a * t78
    t82 = t11 * t27
    t87 = t42 * t27
    t90 = t39 * t17
    t91 = t2 * t42
    t94 = t11 ** 2
    t95 = a * t94
    t96 = t42 * b
    t97 = t96 * t8
    t100 = t39 * a
    t101 = t100 * t94
    t104 = t100 * t17
    t105 = t2 * t64
    t110 = t8 * a
    t111 = t2 * b
    t115 = t23 * b
    t119 = t42 ** 2
    t121 = t17 * t27
    t125 = t39 ** 2
    t127 = t27 ** 2

    t130 = -0.8D1 * t13 * t79 - 0.8D1 * t49 * t82 - 0.4D1 * t52 * t69 - 0.4D1 * t55 * t87 + &
         0.6D1 * t91 * t90 + 0.4D1 * t97 * t95 + 0.4D1 * t64 * t101 + 0.4D1 * t105 * t104 - &
         0.8D1 * t18 * t78 + 0.4D1 * t17 * t111 * t110 + 0.4D1 * t11 * t115 * t110 + t119 * &
         t94 + 0.9D1 * t121 + 0.5D1 * t24 * t27 + t125 * t94 + 0.4D1 * t11 * t127
    
    t133 = t11 * t23
    t136 = t2 * t127
    t158 = t23 * t42
    t166 = t23 ** 2
    
    t167 = 0.4D1 * t23 * t127 + t133 + 0.2D1 * t1 * t24 + 0.8D1 * t1 * t136 + 0.2D1 * t17 * &
         t91 + 0.19D2 * t23 * t69 + 0.23D2 * t2 * t82 - 0.4D1 * t39 * t121 - 0.4D1 * t17 *  &
         t87 + 0.6D1 * t42 * t39 * t94 + 0.2D1 * t2 * t90 + 0.4D1 * t23 * t60 + 0.2D1 *     &
         t24 * t61 + 0.2D1 * t11 * t158 + 0.2D1 * t2 * t125 * t17 + t23 * t125 * t11 + t166
    
    t168 = t130 + t167
    t169 = t168 * t5
    t171 = 0.1D1 / t73 / t72
    t181 = sqrt(0.1D1 + t6 * t2)
    t185 = (a + 0.1D1 / t181 * b) ** 2
    t186 = t185 * t1
    t187 = t186 + t2 + t27
    t188 = t187 ** 2
    t189 = sqrt(t187)
    t191 = 0.1D1 / t189 / t188
    t199 = 0.1D1 / t8
    t201 = a * t199 * t27
    t203 = t1 * x
    t206 = 0.1D1 / x - 0.1D1 / t203 * t3
    t224 = b * t199 * a
    t232 = t11 * t203
    t239 = t203 * t27
    t242 = x * t27
    t248 = t11 * x
    t249 = t39 * t248
    t269 = 0.2D1 * t206 * b * t199
    
    t272 = -0.8D1 * t206 * t2 * t12 * t201 - 0.32D2 * t2 * t203 * b * t79 + 0.4D1 * t206 * &
         t111 * t199 * t104 - 0.8D1 * t206 * t17 * b * t201 + 0.4D1 * t206 * t17 * t2 *    &
         t224 + 0.4D1 * t206 * t133 * t224 + 0.32D2 * t97 * a * t232 + 0.32D2 * t64 *      &
         t100 * t232 - 0.32D2 * t49 * t239 - 0.8D1 * t52 * t242 - 0.16D2 * t2 * t203 *     &
         t87 + 0.36D2 * t91 * t249 + 0.24D2 * t105 * t100 * t248 - 0.48D2 * t248 * t16 *   &
         t78 + 0.24D2 * t248 * t111 * t110 + 0.16D2 * t203 * t115 * t110 + 0.4D1 * t206 *  &
         t96 * t199 * t95 + 0.2D1 * t269 * t101
    
    t275 = t248 * t27
    t290 = t39 * t203
    t298 = t39 * x

    t316 = 0.8D1 * t119 * t232 + 0.54D2 * t275 + 0.8D1 * t125 * t232 + 0.16D2 * t203 * &
         t127 + 0.4D1 * t203 * t23 + 0.4D1 * x * t24 - 0.24D2 * t248 * t87 + 0.4D1 *   &
         t23 * t125 * t203 + 0.16D2 * t23 * t290 + 0.38D2 * t23 * t242 + 0.12D2 * t2 * &
         t125 * t248 + 0.4D1 * t24 * t298 + 0.8D1 * t203 * t158 + 0.16D2 * x * t136 +  &
         0.48D2 * t42 * t39 * t232 + 0.12D2 * t248 * t91 + 0.92D2 * t2 * t239 -        &
         0.24D2 * t39 * t275 + 0.12D2 * t2 * t249
    
    t342 = 0.3D1 * t191 * z * x * (t75 * t58 * t5 / 0.2D1 - 0.5D1 * (t36 + t2 * z) *  &
         t171 * t169) + t191 * (-0.2D1 * t27 + t186 + t2) * (0.3D1 * x * t75 * t168 * &
         t4 + t75 * (t272 + t316) * t5 / 0.2D1 - 0.5D1 / 0.2D1 * (0.4D1 * t290 +      &
         0.2D1 * t2 * t298 + 0.8D1 * t64 * a * t203 + t269 * t63 + 0.4D1 * t42 *      &
         t203 + 0.2D1 * x * t2 + 0.2D1 * t242) * t171 * t169)

    Vy = t342

  end subroutine get_GradB2CrossB_y
  !========================================================================
  subroutine get_GradB2CrossB_z(x, y, z, Vz)

    implicit none

    real, intent(IN) :: x, y, z
    real, intent(OUT):: Vz
    ! Local variables used for optimization
    real :: t1,   t2,   t3,   t4,   t5,   t6,   t7,   t8,   t9,   t12
    real :: t13,  t14,  t17,  t18,  t22,  t23,  t24,  t27,  t28,  t29
    real :: t31,  t32,  t33,  t36,  t37,  t38,  t41,  t42,  t45,  t46
    real :: t50,  t51,  t55,  t59,  t60,  t64,  t65,  t68,  t75,  t78
    real :: t86,  t89,  t90,  t95,  t96,  t99,  t101, t111, t114, t115
    real :: t116, t117, t119, t125, t126, t127, t129, t133, t134, t135
    real :: t142, t145, t151, t152, t155, t157, t161, t178, t182, t183
    real :: t198, t199, t207, t215, t225, t228, t251, t256, t258, t279
    real :: t283, t285, t286, t287, t290, t298, t305, t313, t316, t325
    real :: t337, t346, t379, t402
    real :: a, b
    !------------------------------------------------------------------------
    a = StretchingFactorA
    b = StretchingFactorB
    
    !\
    ! Calculations are done by Maple using CodeGeneration with optimize flag.
    !/

    t1 = x ** 2
    t2 = y ** 2
    t3 = t1 + t2
    t4 = t3 ** 2
    t5 = z ** 2
    t6 = t1 ** 2
    t7 = t6 * t5
    t8 = a ** 2
    t9 = t2 * t8
    t12 = t1 * t5
    t13 = t2 ** 2
    t14 = t13 * t8
    t17 = b ** 2
    t18 = t17 * t5
    t22 = t6 * t1
    t23 = t8 * t22
    t24 = t2 * t17
    t27 = t6 ** 2
    t28 = a * t27
    t29 = 0.1D1 / t1
    t31 = sqrt(t29 * t3)
    t32 = t17 * b
    t33 = t32 * t31
    t36 = t8 * a
    t37 = t36 * t27
    t38 = b * t31
    t41 = t36 * t22
    t42 = t2 * t38
    t45 = t31 * t5
    t46 = a * b
    t50 = t31 * a
    t51 = t2 * b
    t55 = t13 * b
    t59 = a * t45
    t60 = t6 * b
    t64 = t5 ** 2
    t65 = t2 * t64
    t68 = t13 * t17
    t75 = t22 * t5
    
    t78 = -dble(8 * t9 * t7) - dble(4 * t14 * t12) - 0.4D1 * t2 * t6 * t18 + &
         dble(6 * t24 * t23) + 0.4D1 * t33 * t28 + 0.4D1 * t38 * t37 +       &
         0.4D1 * t42 * t41 - 0.8D1 * t22 * t46 * t45 + 0.4D1 * t22 * t51 *   &
         t50 + 0.4D1 * t6 * t55 * t50 - 0.8D1 * t2 * t60 * t59 + 0.8D1 *     &
         t1 * t65 + 0.2D1 * t6 * t68 + 0.23D2 * t2 * dble(t7) + 0.19D2 *     &
         t13 * dble(t12) - 0.4D1 * t8 * t75
    
    t86 = t8 * t6
    t89 = t8 * t1
    t90 = t13 * t2
    t95 = t8 ** 2
    t96 = t95 * t22
    t99 = t95 * t6
    t101 = t17 ** 2
    t111 = t6 * t13
    t114 = t13 ** 2
    
    t115 = -0.4D1 * t22 * t18 + 0.6D1 * t17 * t8 * t27 + 0.2D1 * t2 * dble(t23) + &
         0.4D1 * t13 * t86 + 0.2D1 * t90 * t89 + 0.2D1 * t22 * dble(t24) + 0.2D1 *&
         t2 * t96 + t13 * t99 + t101 * t27 + 0.9D1 * t75 + 0.5D1 * t90 * t5 +     &
         t95 * t27 + 0.4D1 * t6 * t64 + 0.4D1 * t13 * t64 + t111 + 0.2D1 * t1 * t90 + t114
    
    t116 = t78 + t115
    t117 = t116 * t4
    t119 = a * t6
    
    t125 = t86 + t2 * t89 + 0.2D1 * t38 * t119 + t17 * t6 + t1 * t2 + t13 + dble(t12) + t2 * t5
    t126 = t125 ** 2
    t127 = t126 ** 2
    t129 = 0.1D1 / t127 / t125
    t133 = t4 * t3
    t134 = t1 * x
    t135 = t6 * t134
    t142 = t134 * t5
    t145 = x * t5
    t151 = t6 * x
    t152 = t8 * t151
    t155 = 0.1D1 / t31
    t157 = a * t155 * t5
    t161 = 0.1D1 / x - 0.1D1 / t134 * t3
    t178 = t32 * t155
    t182 = b * t155
    t183 = 0.2D1 * t161 * t182
    t198 = t155 * a
    t199 = b * t198
    
    t207 = 0.32D2 * t33 * a * t135 + 0.32D2 * t38 * t36 * t135 - 0.32D2 * dble(t9) * &
         t142 - 0.8D1 * dble(t14) * t145 - 0.16D2 * t2 * t134 * t18 + 0.36D2 *       &
         dble(t24) * t152 - 0.8D1 * t161 * t2 * t60 * t157 + 0.24D2 * t42 * t36 *    &
         t151 - 0.48D2 * t151 * t46 * t45 + 0.24D2 * t151 * t51 * t50 + 0.16D2 *     &
         t134 * t55 * t50 + 0.4D1 * t161 * t178 * t28 + 0.2D1 * t183 * t37 - 0.32D2 *&
         t2 * t134 * b * t59 + 0.4D1 * t161 * t51 * t155 * t41 - 0.8D1 * t161 * t22 *&
         b * t157 + 0.4D1 * t161 * t22 * t2 * t199 + 0.4D1 * t161 * t111 * t199
    
    t215 = t8 * t134
    t225 = t8 * x
    t228 = t151 * t5

    t251 = 0.12D2 * t2 * t95 * t151 + 0.38D2 * t13 * t145 + 0.12D2 * t151 * dble(t24) + &
         0.16D2 * t13 * t215 + 0.4D1 * t13 * t95 * t134 + 0.8D1 * t134 * t68 + 0.12D2 * &
         t2 * t152 + 0.4D1 * t90 * t225 - 0.24D2 * t8 * t228 + 0.16D2 * x * t65 +       &
         0.92D2 * t2 * t142 - 0.24D2 * t151 * t18 + 0.48D2 * t17 * t8 * t135 + 0.8D1 *  &
         t101 * t135 + 0.54D2 * t228 + 0.8D1 * t95 * t135 + 0.16D2 * t134 * t64 +       &
         0.4D1 * t134 * t13 + 0.4D1 * x * t90
    
    t256 = t116 * t133
    t258 = 0.1D1 / t127 / t126
    t279 = sqrt(0.1D1 + t29 * t2)
    t283 = (a + 0.1D1 / t279 * b) ** 2
    t285 = t283 * t1 + t2 + t5
    t286 = t285 ** 2
    t287 = sqrt(t285)
    t290 = 0.1D1 / t287 / t286 * z
    t298 = t2 * y
    t305 = y * t17
    t313 = t6 * t298 * b
    t316 = t13 * y
    t325 = y * t182
    t337 = y * t60
    
    t346 = -0.16D2 * y * t8 * dble(t7) - 0.16D2 * t298 * t8 * dble(t12) - 0.8D1 * y * &
         t6 * t18 + 0.12D2 * t305 * dble(t23) + 0.4D1 * t298 * t182 * t36 * t6 +      &
         0.4D1 * t313 * t198 + 0.4D1 * t1 * t316 * b * t198 + 0.4D1 * y * t178 * a *  &
         t22 + 0.4D1 * t325 * t41 + 0.8D1 * y * t38 * t41 + 0.8D1 * t22 * y * b *     &
         t50 + 0.16D2 * t313 * t50 - 0.8D1 * t337 * t157 - 0.8D1 * t298 * t1 * b *    &
         t157 - 0.16D2 * t337 * t59
    
    t379 = 0.76D2 * t298 * dble(t12) + 0.4D1 * t298 * t99 + 0.12D2 * t316 * t89 + &
         0.4D1 * y * t96 + 0.4D1 * t22 * t305 + 0.16D2 * t1 * y * t64 + 0.8D1 *   &
         t6 * t298 * t17 + 0.4D1 * y * dble(t23) + 0.16D2 * t298 * t86 + 0.46D2 * &
         y * dble(t7) + 0.30D2 * t316 * t5 + 0.16D2 * t298 * t64 + 0.4D1 * t6 *   &
         t298 + 0.12D2 * t1 * t316 + 0.8D1 * t13 * t298
    
    t402 = 0.3D1 * (t290 * y * (0.3D1 * x * t129 * t117 + t129 * (t207 + t251) * &
         t133 / 0.2D1 - 0.5D1 / 0.2D1 * (0.4D1 * t215 + 0.2D1 * t2 * t225 +      &
         0.8D1 * t38 * a * t134 + t183 * t119 + 0.4D1 * t17 * t134 + 0.2D1 * x * &
         t2 + 0.2D1 * t145) * t258 * t256) - t290 * x * (0.3D1 * y * t129 *      &
         t117 + t129 * (t346 + t379) * t133 / 0.2D1 - 0.5D1 / 0.2D1 * (0.2D1 *   &
         y * t89 + 0.2D1 * t325 * a * t1 + 0.2D1 * t1 * y + 0.4D1 * t298 +       &
         0.2D1 * y * t5) * t258 * t256)) 
    

    Vz = t402

  end subroutine get_GradB2CrossB_z
  !========================================================================
  end module get_gradB_components
