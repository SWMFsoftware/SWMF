

 Below is a template for a subroutine that comes after a CONTAINS statement, 
 i.e., routines internal to a module or another subroutine/function.
 The internal routine documentation will become the subsection of 
 the module's LaTeX section, thus saving the conceptual hierarchy.

 The template should be followed as closely as possible, but do not 
 leave in blank parts: 

 If no modules are used, omit the !USES: tag.
 If there are no output arguments, omit the !OUTPUT ARGUMENTS: tag.
 If local variables are not worth documenting, omit the !LOCAL VARIABLES: tag.
 If the source code is not short, do not use the !BOC ... !EOC tags.
 If the revision history is maintained for the module/subroutine containing, 
    this subroutine, omit the !REVISION HISTORY: tag.

The template (can be cut and pasted and EDITED):

  !BOP ======================================================================
  !IROUTINE: put_name_here - put short description here

  !INTERFACE:
  subroutine put_name_here(input_var, output_var, input_output_var)

    !USES:
    use ModSomething
    use CON_Something

    !INPUT ARGUMENTS: 
    real, intent(in) :: input_var           ! short description of input_var

    !OUTPUT ARGUMENTS:
    logical, intent(out) :: output_var      ! short description of output_var

    !INPUT/OUTPUT ARGUMENTS: 
    real, intent(inout) :: input_output_var ! short descr. of input_output_var

    !DESCRIPTION: 
    ! Long description of subroutine in Latex format
    ! Do not repeat information already provided by the above tags.
    !

    !LOCAL VARIABLES:
    real :: AnImportantLocalVariable
    !
    !REVISION HISTORY: 
    ! 04/01/04 My Name - description of change
    !EOP  

    ! local variables not worth of documenting come here

    character(len=*), parameter:: NameSub=NameMod//'put_name_here'

    !------------------------------------------------------------------------
    !BOC
    write(*,*)'Executable statement worth documenting come here'
    !EOC

    write(*,*)'This part should not appear in the documentation'

    end subroutine put_name_here
