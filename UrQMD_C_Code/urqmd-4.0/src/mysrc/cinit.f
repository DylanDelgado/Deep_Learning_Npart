      subroutine c_urqmd_uinit
c
c Revision : 1.0
c
cinput io : flag for call to {\tt input(io)} 
c
c This subroutine calls initialization procedures for uqmd
c i.e. random generater, etc. 
c Routines called before the first (physical) event of uqmd should
c enter here. an exception is the subroutine {\tt init}. {\tt init} should NOT
c be included in uinit. 
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'boxinc.f'
      include 'inputs.f'
      integer io

      io = 0
c     {\tt strini} calculates mixing angles for the meson-multipletts
      call strini

c     print the header into output file
      call output(19)

      firstseed=.true.
      fixedseed=ranseed.gt.0
      if(fixedseed)write(6,*)'fixed random number:',ranseed
      call sseed(ranseed)
      call loginit
      if(CTOption(33).eq.0.or.CTOption(9).eq.0) call loadwtab(io)
c in case of CASCADE mode, the potentials need not be calculated
      if(EoS.ne.0) then
         if(logSky) call potdww
         if(logYuk) call potYuk
         if(logPau) call potPau
         if(logCb)  call potCb
      endif
c
c  calculate the normalization of resonances distribution...
c
      call norm_init

      if(io.ne.0) return

c     do not initialize projectile and target if old event is read in
      if(CTOption(40).ne.0) return

      if(boxflag.eq.0) then
c
c initialize nuclei (projectile and target) and store them
c
c initialize normal projectile
         if(prspflg.eq.0) then
c         if(eos.eq.0) then
            call cascinit(Zp,Ap,1)
c         else
c            write(6,*)'illegal EOS in init.'
c            stop 137
c         endif
         endif
c initialize normal target
         if(At.ne.0) then
            if(trspflg.eq.0) then
c            if(eos.eq.0) then
               call cascinit(Zt,At,2)
c            else
c               write(6,*)'illegal EOS in init.'
c               stop 137
c            endif
            endif
         endif
      endif

      return
      end ! c_urqmd_init
