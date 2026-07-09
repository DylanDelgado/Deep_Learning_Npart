      subroutine c_urqmd_main
      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
      include 'freezeout.f'
      include 'outcom.f'
      include 'curqmd_inc.f'

      character frame*8,proj*8,targ*8,cwin*15
c     save frame, proj, targ, cwin, mevent, nstep, bmin, bmax, dt
      integer i, idx_a
      double precision :: nonpart, cpx, cpy, cpz, cpt, cpt2,
     $ cpcm, cycm ,ceta, cecm, cmass, cpcm2, cchge,
     $ crx, cry, crz, crt, crft
      integer, intent(in) :: c_i, iev, inv
      integer, intent(out) :: iout
      double precision, intent(out) :: dout
      double precision, intent(in) :: din
      integer, intent(inout) :: cmevent, cnstep, cnpart
      double precision, intent(inout) :: cbmin, cbmax, cdt, cmentry,
     $ cmxcoll
      character, dimension(*),intent(inout)::ccwin, cframe, cproj, ctarg
      double precision, dimension(*), intent(inout) :: cd_ary
      real, dimension(*), intent(inout) :: px_ary, py_ary, pz_ary,
     $ e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary,
     $ x_ary, y_ary, z_ary, t_ary, ft_ary
      integer, dimension(*), intent(inout) :: ks_ary, pid_ary
      integer in_test, ou_test
      character in_char_t

      integer j,k,steps,ii,ocharge,ncharge, it1,it2

      integer pdgid

      c_use_external_seed = .false.
      c_quiet = 0
      c_ufile9  = .false.
      c_ufile10 = .false.
      c_ufile13 = .true.
      c_ufile14 = .true.
      c_ufile15 = .true.
      c_ufile16 = .true.
      c_ufile17 = .false.
      c_ufile18 = .true.
      c_ufile19 = .true.
      c_ufile20 = .true.

c From uinit
c Reset of all indexed variables
      call set0
c Setting of global paramters
      call params

c.. display logo
      call urqmdlogo

c init input file
      call c_input(0)
c     call input(0)
c     call uinit(0)
      return

      entry c_urqmd_fin
      write(6,301)'no. of collisions = ',mc/dble(nevents),' (per event)'
      write(6,301)'final particles   = ',mp/dble(nevents),' (per event)'
      write(6,302)'empty events      : ', noc,noc*1d2/dble(nevents)
 301  format(a19,f8.1,a12)
 302  format(a19,i8,' = ',f5.1,'%')
      return

      entry c_urqmd_init
      call c_init_getparams
      call c_set_params_finish
      call c_set_ftn
      if (.not. c_ufile10) then
        call c_post_inputs
      end if
      call c_urqmd_uinit
      write(*,*) Ap, Zp, At, Zt
      write(*,*) CTOption(5), bmin, bdist
      mc=0
      mp=0
      noc=0
      return

      entry c_urqmd_evt(iev)
      call curqmd_evt(iev)
      return

      entry c_on_external_seed()
      c_use_external_seed = .true.
      return

      entry c_off_external_seed()
      c_use_external_seed = .false.
      return

      entry c_set_evt_seed(inv)
      c_external_seed = inv
      return

      entry c_get_evt_seed(iout)
      iout = c_current_event_seed
      return

      entry c_s_cto(inv, din)
      call c_set_cto(inv, int(din))
      return

      entry c_g_cto(inv, dout)
      dout = dble(CTOption(inv))
      return

      entry c_s_ctp(inv, din)
      call c_set_ctp(inv, real(din, 8))
      return

      entry c_g_ctp(inv, dout)
      dout = dble(CTParam(inv))
      return

      entry c_s_nevents(din)
      nevents = int(din)
      return
      entry c_g_nevents(dout)
      dout = dble(nevents)
      return

      entry c_s_Ap(inv)
      Ap=inv
      return
      entry c_g_Ap(iout)
      iout=Ap
      return

      entry c_s_Zp(inv)
      Zp=inv
      return
      entry c_g_Zp(iout)
      iout=Zp
      return

      entry c_wdata(pid_ary, px_ary, py_ary,
     $ pz_ary)
c     do i = 1, nspec
c       pid_ary(i) = pdgid(sityp(i),siso3(i))
c       px_ary(i) = real(pxs(i))
c       py_ary(i) = real(pys(i))
c       pz_ary(i) = real(pzs(i))
c       e_ary(i) = real(p0s(i))
c       m_ary(i) = real(sfmass(i))
c       scharge is integer
c       chge_ary(i) = real(scharge(i))
c     end do
c     do i = 1, npart
c       idx_a = i + nspec
c       pid_ary(idx_a) = pdgid(ityp(i),iso3(i))
c       px_ary(idx_a) = real(px(i))
c       py_ary(idx_a) = real(py(i))
c       pz_ary(idx_a) = real(pz(i))
c       e_ary(idx_a) = real(p0(i))
c       m_ary(idx_a) = real(fmass(i))
c       charge is integer
c       chge_ary(idx_a) = real(charge(i))
c     end do
      do i = 1, npart
        pid_ary(i) = pdgid(ityp(i),iso3(i))
        px_ary(i)  = real(px(i)+ffermpx(i))
        py_ary(i)  = real(py(i)+ffermpy(i))
        pz_ary(i)  = real(pz(i)+ffermpz(i))
      end do
      return

      entry c_t1(c_i, iev, inv, iout, din, dout, cmevent, cnstep,
     $ cnpart, cbmin, cbmax, cdt, cmentry, cmxcoll, ccwin, cframe,
     $ cproj,  ctarg, cd_ary,
     $ px_ary, py_ary, pz_ary, e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary,
     $ x_ary, y_ary, z_ary, t_ary, ft_ary, pid_ary, ks_ary)
      in_test = din
      in_test = c_i
      in_test = iev
      in_test = inv
      iout = ou_test
      dout = ou_test
      in_test = cmevent
      in_test = cnstep
      in_test = cnpart
      in_test = cbmin
      in_test = cbmax
      in_test = cdt
      in_test = cmentry
      in_test = cmxcoll
      in_test = cnpart
      in_test = cbmin
      in_test = cbmax
      in_test = cdt
      in_char_t =  ccwin(1)
      in_char_t = cframe(1)
      in_char_t =  cproj(1)
      in_char_t =  ctarg(1)

      in_test = cd_ary(1)
      in_test = px_ary(1)
      in_test = py_ary(1)
      in_test = pz_ary(1)
      in_test = e_ary(1)
      in_test = m_ary(1)
      in_test = chge_ary(1)
      in_test = pt_ary(1)
      in_test = ycm_ary(1)
      in_test = eta_ary(1)
      in_test = x_ary(1)
      in_test = y_ary(1)
      in_test = z_ary(1)
      in_test = t_ary(1)
      in_test = ft_ary(1)
      in_test = ks_ary(1)
      in_test = pid_ary(1)
      return

      end ! subroutine curqmd
