
      subroutine curqmd_out(opt)
      implicit none
      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'

      include 'curqmd_inc.f'
c
      integer iunit,i,ttime,iu,app,att,zpp,ztt,l
      integer iiunit,isunit, id, pdgid
      integer timestep,itotcoll,iinelcoll
      real*8 sigmatot,ptsigtot,stot,otime
      common /outco2/sigmatot

      integer opt

      character*3 abox3, abox4, abox5
      character*5 abox6
      character*4 reffram
      character*20 aa,ah,ai,ak,abox1,abox2
      character*36 ae,abt
      character*31 aee
      character*15 ab,aj,al,am
      character*13 ac,ag,pds,tds
      character*12 ad
      character*7 af
      character*9 ag2
      character*1 add
      character*171 apa14,apa15,apav,line
      character*2 apa,aop
      character*35 aboxhead

c file15out
      integer ind,ind1,ind2,nin
      integer istr,ich,ii,iid

      real*8 sqrts, sigpart, colldens, cdens,cdens_
      logical bdum,paulibl

      include 'outcom.f'
     
      integer fchg,strit
      character*1 echar
c     temporary arrays for CTO and CTP when read in from old event
c     (before they are overwritten)
      integer CTOtmp(numcto)
      real*8 CTPtmp(numctp)
c     CTOtc and CTPtc are the temporary fields for CTOdc and CTPdc.
      character CTOtc(numcto)*2
      character CTPtc(numctp)*2
      integer ctoforty, ctofoone
      integer ctolines,ctplines
      
      integer iou(13:20)

chp variables necessary for hydro evolution output in f15
      real*8 thydro_start,thydro
chp hydro flag for visualization output
      logical hydro_flag
chp new variable for vis-out to count correct npart
      integer npart_form
      integer, intent(in) :: iin1
      double precision, intent(in) :: din1
      integer, intent(out) :: iou1
      double precision, intent(out) :: dou1


      save 

      data iou/13,14,15,16,17,18,19,20/

c     copy projectile/target info to local vars
      app=ap
      zpp=zp
      att=at
      ztt=zt

c     c_val_proxies
      entry c_s_npart(iin1)
      npart = iin1
      return
      entry c_g_npart(iou1)
      iou1 = npart
      return
      entry c_s_nbar(iin1)
      nbar = iin1
      return
      entry c_g_nbar(iou1)
      iou1 = nbar
      return
      entry c_s_nmes(iin1)
      nmes = iin1
      return
      entry c_g_nmes(iou1)
      iou1 = nmes
      return
      entry c_s_time(din1)
      time = din1
      return
      entry c_g_time(dou1)
      dou1 = time
      return
      entry c_s_acttime(din1)
      acttime = din1
      return
      entry c_g_acttime(dou1)
      dou1 = acttime
      return
      entry c_s_nsteps(iin1)
      nsteps = iin1
      return
      entry c_g_nsteps(iou1)
      iou1 = nsteps
      return
      entry c_s_bimp(din1)
      bimp = din1
      return
      entry c_g_bimp(dou1)
      dou1 = bimp
      return
      entry c_s_nspec(iin1)
      nspec = iin1
      return
      entry c_g_nspec(iou1)
      iou1 = nspec
      return
      entry c_s_app(iin1)
      App = iin1
      return
      entry c_g_app(iou1)
      iou1 = App
      return
      entry c_s_att(iin1)
      Att = iin1
      return
      entry c_g_att(iou1)
      iou1 = Att
      return
      entry c_s_sigmatot(din1) ! in mbarn
      sigmatot = din1
      return
      entry c_g_sigmatot(dou1)
      dou1 = sigmatot
      return

      entry c_s_cnpartit(iin1)
      cnpartit = iin1
      return
      entry c_g_cnpartit(iou1)
      iou1 = cnpartit
      return

      entry c_s_cmul(iin1)
      cmul = iin1
      return
      entry c_g_cmul(iou1)
      iou1 = cmul
      return
c     end of c_val_proxies

      entry out_test(opt, iin1, din1,
     $iou1, dou1)
      i=opt
      i=iin1
      i=int(din1)
      iou1=i
      dou1=dble(i)
      return

      end subroutine ! curqmd_out
