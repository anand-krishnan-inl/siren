module transient
use kind, only : rk, ik
implicit none

interface
  real(c_double) function transient_albedo(time, albedo_coeff, pow_lhs, pow_rhs) &
      bind(C, name='transient_albedo')
    use, intrinsic :: iso_c_binding, only : c_double
    real(c_double), value, intent(in) :: time, albedo_coeff, pow_lhs, pow_rhs
  endfunction transient_albedo
endinterface

private

public :: transient_read, transient_summary, transient_solve, transient_cleanup

contains

  subroutine transient_read(fname, dnd)
    use fileio, only : fileio_open_read
    use exception_handler, only : exception_fatal
    use delayed_neutron_data, only : DelayedNeutronData
    character(*), intent(in) :: fname
    type(DelayedNeutronData), intent(out) :: dnd

    integer, parameter :: iounit = 15
    character(*), parameter :: comment_char = '#'

    character(1024) :: line, card
    integer :: ios

    call fileio_open_read(trim(adjustl(fname)), iounit)

    do
      read(iounit, '(a)', iostat=ios) line
      if (ios /= 0) then
        exit
      endif

      if (line == '') then
        cycle
      endif

      line = trim(adjustl(line))
      if (line(1:1) == comment_char) then
        cycle
      endif

      read(line, *) card
      backspace(iounit)

      select case (card)
        case ('nd')
          read(iounit, *) card, dnd%nd
          allocate(dnd%beta(dnd%nd))
          allocate(dnd%lamd(dnd%nd))
          dnd%beta = 0.0_rk
          dnd%lamd = 0.0_rk
        case ('beta')
          read(iounit, *) card, dnd%beta
        case ('lambda')
          read(iounit, *) card, dnd%lamd
        case ('ng')
          read(iounit, *) card, dnd%ng
          allocate(dnd%vel(dnd%ng))
          dnd%vel = 0.0_rk
        case ('vel')
          read(iounit, *) card, dnd%vel
        case ('deltat')
          read(iounit, *) card, dnd%deltat
        case ('tend')
          read(iounit, *) card, dnd%tend
        case ('reference')
          read(iounit, *) card, dnd%reference
        case ('tedit')
          read(iounit, *) card, dnd%ntime
          backspace(iounit)
          allocate(dnd%tedit(dnd%ntime))
          dnd%tedit = 0.0_rk
          read(iounit, *) card, dnd%ntime, dnd%tedit
        case default
          call exception_fatal('unknown transient input card: ' // &
            trim(adjustl(card)))
      endselect
    enddo

    close(iounit)
  endsubroutine transient_read

  subroutine transient_summary(dnd)
    use output, only : output_write
    use delayed_neutron_data, only : DelayedNeutronData
    type(DelayedNeutronData), intent(in) :: dnd

    character(2**16) :: line, tmp
    integer(ik) :: d, g, i

    call output_write('=== TRANSIENT ===')

    write(line, '(a,i0)') 'Number of delayed neutron precursor families: ', dnd%nd
    call output_write(line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%beta(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('beta = ' // line)

    line = ''
    do d = 1,dnd%nd
      write(tmp, '(es9.2)') dnd%lamd(d)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! d = 1,dnd%nd
    call output_write('lambda [1/s] = ' // line)

    write(line, '(a,i0)') 'Number of energy groups in transient data: ', dnd%ng
    call output_write(line)
    
    line = ''
    do g = 1,dnd%ng
      write(tmp, '(es9.2)') dnd%vel(g)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! g = dnd%ng
    call output_write('velocity [cm/s] = ' // line)

    write(line, '(a,es9.2)') 'Time-step size [s] = ', dnd%deltat
    call output_write(line)
    write(line, '(a,es9.2)') 'Tend [s] = ', dnd%tend
    call output_write(line)
    
    call output_write('Transient reference case: ' &
      // trim(adjustl(dnd%reference)))

    write(line, '(a,i0)') 'Number of time edit points: ', dnd%ntime
    call output_write(line)

    line = ''
    do i = 1,dnd%ntime
      write(tmp, '(es9.2)') dnd%tedit(i)
      line = trim(adjustl(line)) // ' ' // trim(adjustl(tmp))
    enddo ! i = 1,dnd%ntime
    call output_write('Edit times [s] = ' // trim(adjustl(line)))

    call output_write('')
  endsubroutine transient_summary

  subroutine transient_init_precursors(nx, mat_map, xslib, dnd, kcrit, flux, prec)
    use xs, only : XSLibrary
    use delayed_neutron_data, only : DelayedNeutronData
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: kcrit
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: prec(:,:) ! (nx,nd)

    integer(ik) :: i, mthis
    real(rk) :: fsrc

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsrc = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))/kcrit
        prec(i,:) = dnd%beta/dnd%lamd * fsrc
      else
        prec(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx
  endsubroutine transient_init_precursors

  subroutine transient_update_precursors(nx, mat_map, xslib, dnd, kcrit, flux, &
    prec_old, prec)
    use xs, only : XSLibrary
    use delayed_neutron_data, only : DelayedNeutronData
    integer(ik), intent(in) :: nx
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: kcrit
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(in) :: prec_old(:,:) ! (nx,nd)
    real(rk), intent(out) :: prec(:,:) ! (nx,nd)

    integer(ik) :: i, mthis
    real(rk) :: fsrc

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsrc = sum(xslib%mat(mthis)%nusf(:) * flux(i,:))/kcrit
        prec(i,:) = dnd%deltat / (1.0_rk + dnd%deltat * dnd%lamd(:)) &
          * (dnd%beta(:)*fsrc + prec_old(i,:)/dnd%deltat)
      else
        prec(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx
  endsubroutine transient_update_precursors

  subroutine transient_build_kinsrc(nx, dx, mat_map, xslib, dnd, flux, prec, kinsrc)
    use xs, only : XSLibrary
    use delayed_neutron_data, only : DelayedNeutronData
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(in) :: prec(:,:) ! (nx,nd)
    real(rk), intent(out) :: kinsrc(:,:) ! (nx,ngroup)

    integer(ik) :: i, mthis
    real(rk) :: prec_src

    do i = 1,nx
      mthis = mat_map(i)
      kinsrc(i,:) = flux(i,:) / (dnd%vel(:) * dnd%deltat) * dx(i)
      if (xslib%mat(mthis)%is_fiss) then
        prec_src = sum((dnd%lamd * prec(i,:)) / (1.0_rk + dnd%lamd*dnd%deltat))
        kinsrc(i,:) = kinsrc(i,:) &
          + xslib%mat(mthis)%chi_delay * prec_src * dx(i)
      endif
    enddo ! i = 1,nx
  endsubroutine transient_build_kinsrc

  subroutine transient_build_fsource(nx, dx, mat_map, xslib, dnd, flux, fsrc)
    use xs, only : XSLibrary
    use delayed_neutron_data, only : DelayedNeutronData
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    real(rk), intent(in) :: flux(:,:) ! (nx,ngroup)
    real(rk), intent(out) :: fsrc(:,:) ! (nx,ngroup)

    integer(ik) :: i
    integer(ik) :: mthis
    real(rk) :: beff, fsum
    real(rk), allocatable :: chi_tilde(:)

    allocate(chi_tilde(xslib%ngroup))

    beff = sum(dnd%beta)

    do i = 1,nx
      mthis = mat_map(i)
      if (xslib%mat(mthis)%is_fiss) then
        fsum = sum(xslib%mat(mthis)%nusf(:) * flux(i,:)) * dx(i)
        chi_tilde = (1.0 - beff) * xslib%mat(mthis)%chi &
          + xslib%mat(mthis)%chi_delay * sum((dnd%lamd * dnd%beta * dnd%deltat)&
          /(1.0_rk + dnd%lamd * dnd%deltat))
        fsrc(i,:) = chi_tilde * fsum
      else
        fsrc(i,:) = 0.0_rk
      endif
    enddo ! i = 1,nx

    deallocate(chi_tilde)
  endsubroutine transient_build_fsource

  subroutine transient_build_diagonal(nx, dx, mat_map, xslib, dnd, &
      boundary_left, boundary_right, albedo_alpha, time, dia)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    use delayed_neutron_data, only : DelayedNeutronData
    use transient_mms, only : transient_mms_sigma_a
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xslib
    type(DelayedNeutronData), intent(in) :: dnd
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: albedo_alpha
    real(rk), intent(in) :: time
    real(rk), intent(out) :: dia(:,:) ! (nx,ngroup)

    integer(ik) :: i, g
    integer(ik) :: mprev, mthis, mnext
    real(rk) :: dprev, dnext
    real(rk) :: xt, xleft

    ! BC at x=0, i=1
    mthis = mat_map(1)
    mnext = mat_map(2)
    do g = 1,xslib%ngroup
      dnext = 2 &
        * (xslib%mat(mthis)%diffusion(g) / dx(1) * xslib%mat(mnext)%diffusion(g) / dx(2)) &
        / (xslib%mat(mthis)%diffusion(g) / dx(1) + xslib%mat(mnext)%diffusion(g) / dx(2))
      if (dnd%reference == 'mms') then
        xt = transient_mms_sigma_a(0.5_rk * dx(1), time) &
          + xslib%mat(mthis)%scatter(g,g,1)
      elseif (dnd%reference == 'bc-cont' .or. dnd%reference == 'bc-cont-capi') then
        xt = transient_boundary_control_continuous_sigma_t(xslib, 0.5_rk * dx(1), time)
      else
        xt = xslib%mat(mthis)%sigma_t(g)
      endif
      select case (boundary_left)
        case ('zero')
          dia(1,g) = dnext &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(1) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(1) &
            + 2 * xslib%mat(mthis)%diffusion(g) / dx(1)
        case ('mirror')
          dia(1,g) = dnext &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(1) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(1)
        case ('albedo')
          dia(1,g) = dnext &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(1) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(1) &
            + 1.0_rk / (1.0_rk / albedo_alpha + 0.5_rk * dx(1) / xslib%mat(mthis)%diffusion(g))
        case default
          call exception_fatal(&
            'Unknown boundary_left in transient_build_diagonal: ' &
            // trim(adjustl(boundary_left)))
      endselect
    enddo ! g = 1,xslib%ngroup

    do g = 1,xslib%ngroup
      xleft = dx(1)
      do i = 2,nx-1

        mprev = mat_map(i-1)
        mthis = mat_map(i)
        mnext = mat_map(i+1)

        dprev = 2 &
          * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mprev)%diffusion(g) / dx(i-1)) &
          / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mprev)%diffusion(g) / dx(i-1))
        dnext = 2 &
          * (xslib%mat(mthis)%diffusion(g) / dx(i) * xslib%mat(mnext)%diffusion(g) / dx(i+1)) &
          / (xslib%mat(mthis)%diffusion(g) / dx(i) + xslib%mat(mnext)%diffusion(g) / dx(i+1))

        if (dnd%reference == 'mms') then
          xt = transient_mms_sigma_a(xleft + 0.5_rk * dx(i), time) &
            + xslib%mat(mthis)%scatter(g,g,1)
          xleft = xleft + dx(i)
        elseif (dnd%reference == 'bc-cont' .or. dnd%reference == 'bc-cont-capi') then
          xt = transient_boundary_control_continuous_sigma_t(xslib, &
            xleft + 0.5_rk * dx(i), time)
          xleft = xleft + dx(i)
        else
          xt = xslib%mat(mthis)%sigma_t(g)
        endif

        dia(i,g) = dprev + dnext &
          + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(i) &
          + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(i)

      enddo ! i = 2,nx-1
    enddo ! g = 1,xslib%ngroup

    ! BC at x=L, i=n
    mprev = mat_map(nx-1)
    mthis = mat_map(nx)
    do g = 1,xslib%ngroup
      dprev = 2 &
        * (xslib%mat(mthis)%diffusion(g) / dx(nx) * xslib%mat(mprev)%diffusion(g) / dx(nx-1)) &
        / (xslib%mat(mthis)%diffusion(g) / dx(nx) + xslib%mat(mprev)%diffusion(g) / dx(nx-1))
      if (dnd%reference == 'mms') then
        xt = transient_mms_sigma_a(sum(dx(1:nx-1)) + 0.5_rk*dx(nx), time) &
          + xslib%mat(mthis)%scatter(g,g,1)
      elseif (dnd%reference == 'bc-cont' .or. dnd%reference == 'bc-cont-capi') then
        xt = transient_boundary_control_continuous_sigma_t(xslib, &
          sum(dx(1:nx-1)) + 0.5_rk*dx(nx), time)
      else
        xt = xslib%mat(mthis)%sigma_t(g)
      endif
      select case (boundary_right)
        case ('zero')
          dia(nx,g) = dprev &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(nx) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(nx) &
            + 2 * xslib%mat(mthis)%diffusion(g) / dx(nx)
        case ('mirror')
          dia(nx,g) = dprev &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(nx) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(nx)
        case ('albedo')
          dia(nx,g) = dprev &
            + (xt - xslib%mat(mthis)%scatter(g,g,1)) * dx(nx) &
            + (1.0_rk / (dnd%vel(g) * dnd%deltat)) * dx(nx) &
            + 1.0_rk / (1.0_rk / albedo_alpha + 0.5_rk * dx(nx) / xslib%mat(mthis)%diffusion(g))
        case default
          call exception_fatal(&
            'Unknown boundary_right in transient_build_diagonal: ' &
            // trim(adjustl(boundary_right)))
      endselect
    enddo ! g = 1,xslib%ngroup
  endsubroutine transient_build_diagonal

  subroutine transient_solve(fname_stub, nx, dx, mat_map, xs, dnd, &
      boundary_left, boundary_right, albedo_coeff_in, &
      phi_tol, max_iter, keff, flux)
    use xs, only : XSLibrary
    use output, only : output_write
    use power, only : power_calculate, power_total
    use diffusion, only : diffusion_build_matrix, &
      diffusion_build_upscatter, diffusion_build_downscatter
    use linalg, only : trid
    use fileio, only : fileio_open_write
    use output, only : output_power_csv, output_flux_csv
    use albedo, only : albedo_calculate_alpha
    use delayed_neutron_data, only : DelayedNeutronData
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: fname_stub
    integer(ik), intent(in) :: nx
    real(rk), intent(in) :: dx(:) ! (nx)
    integer(ik), intent(in) :: mat_map(:) ! (nx)
    type(XSLibrary), intent(in) :: xs
    type(DelayedNeutronData), intent(in) :: dnd
    character(*), intent(in) :: boundary_left, boundary_right
    real(rk), intent(in) :: albedo_coeff_in ! A \in [-1,1]
    real(rk), intent(in) :: phi_tol
    integer(ik), intent(in) :: max_iter
    real(rk), intent(in) :: keff
    real(rk), intent(inout) :: flux(:,:) ! (nx,ngroup)

    real(rk), allocatable :: sub(:,:), dia(:,:), sup(:,:) ! (nx,ngroup)
    real(rk), allocatable :: dia_copy(:) ! (nx)
    real(rk), allocatable :: fsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: upsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: downsource(:) ! (nx)
    real(rk), allocatable :: kinsource(:,:) ! (nx,ngroup)
    real(rk), allocatable :: q(:) ! (nx)

    real(rk), allocatable :: prec(:,:), prec_old(:,:) ! (nx,nd)

    real(rk), allocatable :: flux_old(:,:) ! (nx,ngroup)
    real(rk), allocatable :: power(:) ! (nx)

    integer(ik) :: step, iter, g, idx
    real(rk) :: tfinal
    real(rk) :: phi_conv

    real(rk) :: albedo_alpha ! \alpha \in [0,\infty]
    real(rk) :: albedo_coeff
    real(rk) :: pow_lhs, pow_rhs

    character(1024) :: line

    type(XSLibrary) :: xslib

    integer(ik), parameter :: iout = 19
    real(rk), parameter :: omega = 1.9_rk ! over-relaxation
    real(rk), parameter :: power_init = 1.0_rk ! solving for relative power

    character(1024) :: fname_kin, fname_tpower, fname_tflux

    ! store a mutable copy so I can modify it during the transient
    xslib = xs
    albedo_coeff = albedo_coeff_in

    allocate(sub(nx,xslib%ngroup))
    allocate(dia(nx,xslib%ngroup))
    allocate(sup(nx,xslib%ngroup))
    allocate(dia_copy(nx))
    sub = 0.0_rk
    dia = 0.0_rk
    sup = 0.0_rk
    dia_copy = 0.0_rk

    allocate(fsource(nx,xslib%ngroup))
    allocate(upsource(nx,xslib%ngroup))
    allocate(downsource(nx))
    allocate(kinsource(nx,xslib%ngroup))
    allocate(q(nx))
    fsource = 0.0_rk
    upsource = 0.0_rk
    downsource = 0.0_rk
    kinsource = 0.0_rk
    q = 0.0_rk

    allocate(flux_old(nx,xslib%ngroup))
    flux_old = 0.0_rk

    ! compute intital power and normalize
    allocate(power(nx))
    power = 0.0_rk
    call power_calculate(nx, mat_map, xslib, flux, power)
    flux = flux * (power_init / power_total(dx, power))
    ! now, recompute with normalized flux
    call power_calculate(nx, mat_map, xslib, flux, power)

    allocate(prec(nx,dnd%nd))
    allocate(prec_old(nx,dnd%nd))
    prec = 0.0_rk
    prec_old = 0.0_rk
    call transient_init_precursors(nx, mat_map, xslib, dnd, keff, flux, prec)

    fname_kin = trim(adjustl(fname_stub)) // trim(adjustl('_transient.csv'))
    call fileio_open_write(fname_kin, iout)

    call output_write('=== TRANSIENT CALCULATION ===')
    call output_write('elapt [s] , deltat [s] , rel. power')
    write(line, '(es13.6, " , ", es13.6, " , ", es13.6)') &
      0.0_rk, dnd%deltat, power_total(dx, power)
    call output_write(line)
    write(iout, '(a)') trim(adjustl(line))

    ! only the diagonal needs to be modified compared to steady-state solution
    ! the off-diagonal never change (if diffusion coeff never changes)
    albedo_alpha = albedo_calculate_alpha(albedo_coeff)
    call diffusion_build_matrix(&
      nx, dx, mat_map, xslib, boundary_left, boundary_right, albedo_alpha, &
      sub, dia, sup)

    ! edit requested before first step
    idx = 0
    tfinal = 0.0_rk
    if (allocated(dnd%tedit)) then
      if (any(dnd%tedit < 0.5_rk*dnd%deltat)) then
        idx = idx + 1
        write(fname_tpower, '(a,i0,a)') trim(adjustl(fname_stub)) // '_power_t', idx, '.csv'
        call output_power_csv(fname_tpower, nx, dx, power)
        write(line, '(a,es13.6,a)') ' --- writing power at time t=', tfinal, ' on ' // trim(adjustl(fname_tpower))
        call output_write(line)
        write(fname_tflux, '(a,i0,a)') trim(adjustl(fname_stub)) // '_flux_t', idx, '.csv'
        call output_flux_csv(fname_tflux, nx, xslib%ngroup, dx, flux)
        write(line, '(a,es13.6,a)') ' --- writing flux  at time t=', tfinal, ' on ' // trim(adjustl(fname_tflux))
        call output_write(line)
      endif
    endif

    ! TIME LOOP
    step = 0
    do
      ! This works for uniform time steps.
      ! It is more accurate to use a multiplication rather than floating-point
      ! addition which will accumulate round-off.
      step = step + 1
      tfinal = step * dnd%deltat
      
      prec_old = prec

      call transient_build_kinsrc(nx, dx, mat_map, xslib, dnd, flux, prec, kinsource)

      ! recompute an albedo coefficient
      ! this is a bit verbose to support boundary control
      ! the boundary controller may need the power at the boundaries
      pow_lhs = -(power(2)-power(1))/(dx(2)+dx(1))*dx(1) + power(1)
      pow_rhs = (power(nx)-power(nx-1))/(dx(nx-1)+dx(nx))*dx(nx) + power(nx)
      albedo_coeff = &
        transient_update_albedo(dnd%reference, tfinal, albedo_coeff, &
        pow_lhs, pow_rhs)
      albedo_alpha = albedo_calculate_alpha(albedo_coeff)

      ! update xs and re-build the diagonal
      ! the rest of the entries in the matrix do not change
      call transient_update_xs(dnd%reference, tfinal, xslib)
      call transient_build_diagonal(&
        nx, dx, mat_map, xslib, dnd, &
        boundary_left, boundary_right, albedo_alpha, tfinal, dia)

      ! iterative scheme is necesary for one-group at-a-time problem
      do iter = 1,max_iter
        flux_old = flux

        call transient_build_fsource(nx, dx, mat_map, xslib, dnd, flux, fsource)
        call diffusion_build_upscatter(nx, dx, mat_map, xslib, flux, upsource)

        do g = 1,xslib%ngroup
          call diffusion_build_downscatter(nx, dx, mat_map, xslib, flux, g, &
            downsource)
          q = fsource(:,g)/keff + upsource(:,g) + kinsource(:,g) + downsource

          ! Note that we only trample over the diagonal
          dia_copy = dia(:,g)
          call trid(nx, sub(:,g), dia_copy, sup(:,g), q, flux(:,g))

          ! over-relaxation
          flux(:,g) = flux_old(:,g) + omega * (flux(:,g) - flux_old(:,g))
        enddo ! g = 1,xslib%ngroup

        call transient_update_precursors(nx, mat_map, xslib, dnd, keff, flux, &
          prec_old, prec)

        phi_conv = maxval(abs(flux - flux_old))/maxval(flux)

        if ((phi_conv < phi_tol) .and. (iter > 1)) then
          exit
        endif
      enddo ! iter = 1,max_iter

      if (iter > max_iter) then
        call exception_fatal('Maximum iterations exceeded.')
      endif

      call power_calculate(nx, mat_map, xslib, flux, power)
      write(line, '(es13.6, " , ", es13.6, " , ", es13.6)') &
        tfinal, dnd%deltat, power_total(dx, power)
      call output_write(line)
      write(iout, '(a)') trim(adjustl(line))

      if (allocated(dnd%tedit)) then
        if (any(abs(tfinal - dnd%tedit) < 0.5_rk*dnd%deltat)) then
          idx = idx + 1
          write(fname_tpower, '(a,i0,a)') trim(adjustl(fname_stub)) // '_power_t', idx, '.csv'
          call output_power_csv(fname_tpower, nx, dx, power)
          write(line, '(a,es13.6,a)') ' --- writing power at time t=', tfinal, ' on ' // trim(adjustl(fname_tpower))
          call output_write(line)
          write(fname_tflux, '(a,i0,a)') trim(adjustl(fname_stub)) // '_flux_t', idx, '.csv'
          call output_flux_csv(fname_tflux, nx, xslib%ngroup, dx, flux)
          write(line, '(a,es13.6,a)') ' --- writing flux  at time t=', tfinal, ' on ' // trim(adjustl(fname_tflux))
          call output_write(line)
        endif
      endif

      if (tfinal > dnd%tend) then
        call output_write('Transient Completed!')
        exit
      endif
    enddo ! TIME LOOP

    call output_write('')
    flush(iout)
    close(iout)

    deallocate(sub, dia, sup)
    deallocate(dia_copy)
    deallocate(fsource, upsource, downsource, kinsource, q)
    deallocate(prec, prec_old)
    deallocate(power)
    deallocate(flux_old)
  endsubroutine transient_solve

  subroutine transient_cleanup(dnd)
    use delayed_neutron_data, only : DelayedNeutronData
    type(DelayedNeutronData), intent(out) :: dnd
    if (allocated(dnd%beta)) then
      deallocate(dnd%beta)
    endif
    if (allocated(dnd%lamd)) then
      deallocate(dnd%lamd)
    endif
    if (allocated(dnd%vel)) then
      deallocate(dnd%vel)
    endif
  endsubroutine transient_cleanup

  subroutine transient_update_xs(name, time, xs)
    use xs, only : XSLibrary
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: name
    real(rk), intent(in) :: time
    type(XSLibrary), intent(inout) :: xs

    logical, save :: first = .true.
    real(rk), save :: sigma0
    real(rk) :: fac
    character(len=64) :: buf
    integer :: stat

    select case (name)
      case ('anl-slab-6-a1')
        if (first) then
          sigma0 = xs%mat(1)%sigma_t(2)
          first = .false.
        endif
        if (time <= 1.0_rk) then
          xs%mat(1)%sigma_t(2) = sigma0 + time * 0.03_rk * sigma0
        endif
      case ('anl-slab-6-a2')
        if (first) then
          sigma0 = xs%mat(1)%sigma_t(2)
          first = .false.
        endif
        if (time <= 1.0_rk) then
          xs%mat(1)%sigma_t(2) = sigma0 - time * 0.01_rk * sigma0
        endif
      case ('anl-slab-6-a3', 'anl-slab-6-a4')
        if (first) then
          sigma0 = xs%mat(1)%sigma_t(2)
          first = .false.
        endif
        if (time <= 0.01_rk) then
          xs%mat(1)%sigma_t(2) = sigma0 - time/0.01_rk * 0.05_rk * sigma0
        endif
      case ('albedo-pid', 'albedo-pid-capi')
        if (first) then
          xs%mat(1)%sigma_t(1) = 0.99_rk * xs%mat(1)%sigma_t(1)
          first = .false.
        endif
      case ('bc-disc', 'bc-disc-capi')
        if (first) then
          ! SIREN_BC_DISC_FRAC sizes the step (fractional sigma_t drop in F1);
          ! unset keeps the original 1% step exactly
          fac = 0.99_rk
          call get_environment_variable('SIREN_BC_DISC_FRAC', buf, status=stat)
          if (stat == 0) then
            read(buf, *) fac
            fac = 1.0_rk - fac
          endif
          xs%mat(1)%sigma_t(1) = fac * xs%mat(1)%sigma_t(1)
          first = .false.
        endif
      case ('bc-cont', 'bc-cont-capi')
        ! do nothing here
        ! this is handeled in transient_build_diagonal
      case ('mms')
        ! do nothing here
        ! this is handled in transient_build_diagonal
      case ('null', 'albedo')
        ! do nothing
      case default
        call exception_fatal('Unknown transient reference name: ' &
          // trim(adjustl(name)))
    endselect
  endsubroutine transient_update_xs

  real(rk) function transient_update_albedo(name, time, albedo_coeff, &
      pow_lhs, pow_rhs)
    use exception_handler, only : exception_fatal
    character(*), intent(in) :: name
    real(rk), intent(in) :: time
    real(rk), intent(in) :: albedo_coeff
    real(rk), intent(in) :: pow_lhs, pow_rhs
    select case (name)
      case ('null', &
          'anl-slab-6-a1', 'anl-slab-6-a2', 'anl-slab-6-a3', 'anl-slab-6-a4', &
          'mms', 'bc-disc', 'bc-cont')
        ! do nothing
        ! transparent pass-through
        transient_update_albedo = albedo_coeff
      case ('albedo')
        transient_update_albedo = max(1.0_rk - (time/64.0_rk), 0.0_rk)
      case ('albedo-pid')
        transient_update_albedo = &
          transient_albedo_pid(time, albedo_coeff, pow_rhs)
      case ('albedo-pid-capi')
        transient_update_albedo = &
          transient_albedo_capi(time, albedo_coeff, pow_lhs, pow_rhs)
      case ('bc-disc-capi', 'bc-cont-capi')
        transient_update_albedo = &
          transient_albedo_capi(time, albedo_coeff, pow_lhs, pow_rhs)
      case default
        transient_update_albedo = 0.0_rk
        call exception_fatal('Unknown transient albedo reference name: ' &
          // trim(adjustl(name)))
    endselect
  endfunction transient_update_albedo

  real(rk) function transient_albedo_pid(time, albedo_coeff, pboundary)
    real(rk), intent(in) :: time ! [s]
    real(rk), intent(in) :: albedo_coeff
    real(rk), intent(in) :: pboundary

    logical, save :: first = .true.

    real(rk) :: dt
    real(rk) :: error, derror_dt
    real(rk), save :: p0
    real(rk), save :: prev_error, prev_time
    real(rk), save :: interror

    real(rk), parameter :: kp = 0.004_rk
    real(rk), parameter :: ki = 0.0_rk
    real(rk), parameter :: kd = 0.001_rk

    if (first) then
      p0 = pboundary
      transient_albedo_pid = albedo_coeff
      error = 0.0_rk
      interror = 0.0_rk
      first = .false.
    else
      ! maintain constant pboundary
      error = p0 - pboundary
      dt = time - prev_time
      derror_dt = (error - prev_error) / dt
      interror = interror + error * dt
      transient_albedo_pid = albedo_coeff &
        + kp * error + kd * derror_dt + ki * interror
      ! clamp
      transient_albedo_pid = min(max(transient_albedo_pid, -1.0_rk), 1.0_rk)
      ! write(*,'(a,1x,es9.2,1x,a,1x,es9.2,1x,a,1x,es9.2)') &
      !   'error', error, 'albedo', transient_albedo_pid, 'pboundary', pboundary
    endif
    prev_error = error
    prev_time = time
  endfunction transient_albedo_pid

  real(rk) function transient_albedo_capi(time, albedo_coeff, pow_lhs, pow_rhs)
    use, intrinsic :: iso_c_binding, only : c_double
    real(rk), intent(in) :: time
    real(rk), intent(in) :: albedo_coeff
    real(rk), intent(in) :: pow_lhs, pow_rhs

    real(c_double) :: c_time
    real(c_double) :: c_albedo_coeff
    real(c_double) :: c_pow_lhs ! power at x=0 (left-hand-side)
    real(c_double) :: c_pow_rhs ! power at x=L (right-hand-side)
    real(c_double) :: c_alb

    c_time = real(time, c_double)
    c_albedo_coeff = real(albedo_coeff, c_double)
    c_pow_lhs = real(pow_lhs, c_double)
    c_pow_rhs = real(pow_rhs, c_double)

    c_alb = transient_albedo(c_time, c_albedo_coeff, c_pow_lhs, c_pow_rhs)

    transient_albedo_capi = real(c_alb, rk)
  endfunction transient_albedo_capi

  real(rk) function transient_boundary_control_continuous_sigma_t(xslib, x, t)
    use xs, only : XSLibrary
    type(XSLibrary), intent(in) :: xslib
    real(rk), intent(in) :: x ! [cm]
    real(rk), intent(in) :: t ! [t]

    logical, save :: first = .true.
    real(rk), save :: xs0
    real(rk), save :: thold = huge(1.0_rk)
    real(rk), parameter :: alpha = 2e-8_rk
    character(len=64) :: buf
    integer :: stat

    if (first) then
      xs0 = xslib%mat(1)%sigma_t(1)
      ! optional ramp-and-hold: freeze the perturbation at t = SIREN_BC_CONT_THOLD [s]
      call get_environment_variable('SIREN_BC_CONT_THOLD', buf, status=stat)
      if (stat == 0) read(buf, *) thold
      first = .false.
    endif

    transient_boundary_control_continuous_sigma_t = &
      xs0 - alpha * (50._rk - x)**2 * min(t, thold)**2
  endfunction transient_boundary_control_continuous_sigma_t

endmodule transient
