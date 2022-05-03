!--------------------------------------------------------------------------
subroutine read_oddalpha_file (file_odd_alpha, number_alpha, alpha_orb, spread_orb, &
                               delta_spread_orb, center_orb, delta_center_orb, iflag_alpha)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode_id, ionode
  USE mp,         ONLY : mp_bcast
  USE kinds, only : DP
  
  !
  implicit none
  !
  character (len=*)  :: file_odd_alpha
  !
  integer :: number_alpha, iflag_alpha
  !
  real(DP):: alpha_orb(number_alpha), spread_orb(number_alpha), & 
             delta_spread_orb(number_alpha), center_orb (number_alpha, 3), &
             delta_center_orb(number_alpha)
  !
  integer :: iun_alphafile, ia, ipol, ios
  !
  if (ionode) then
    !
    if (file_odd_alpha == ' ') call errore ('file_odd_alpha', 'filename missing', 1)
    !
    iun_alphafile = 4
    !
    write( stdout, '(5x,"Reading orbital alpha from file  ",a)') TRIM(file_odd_alpha)
    !
    open (unit = iun_alphafile, file = file_odd_alpha, form = 'formatted', &
          status = 'old', iostat = ios)
    !
    rewind (iun_alphafile)
    !
    read (iun_alphafile, *) number_alpha 
    !
    do ia = 1, number_alpha
       read (iun_alphafile, * ) alpha_orb(ia), spread_orb(ia), delta_spread_orb(ia), &
                       ( center_orb(ia,ipol), ipol=1,3 ), delta_center_orb(ia), &
                       iflag_alpha
    enddo 
    !
    close (unit = iun_alphafile)
    !
  endif
  !
  call mp_bcast( number_alpha, ionode_id )
  call mp_bcast( iflag_alpha, ionode_id )
  call mp_bcast( alpha_orb, ionode_id )
  call mp_bcast( spread_orb, ionode_id )
  call mp_bcast( delta_spread_orb, ionode_id )
  call mp_bcast( center_orb, ionode_id )
  call mp_bcast( delta_center_orb, ionode_id )
  call mp_bcast( iflag_alpha, ionode_id )
  !
  return
  !
end subroutine read_oddalpha_file

