      Program NewCT
c.......................................................................
c.... Program to study a multi-turn extraction based on island-capture
c.... It allows:
c.... - Simulation of evolution of beam distribution with graphical 
c....   output
c.... - Computation of initial and final beam distribution
c.... - Multi-extractions
c.... The detailed options are selected via the parameter Icomp, which 
c.... can assume the following values:
c.... 0 - Distribution evolution, plotting
c.... 1 - Particle identification
c.... 2 - Computation of distribution parameters (initial and final), 
c....     including optical mismatch
c.... 3 - Distribution evolution and plotting (histogram) with new option 
c....     of dumping final beam distribution
c.... 4 - Multiple extractions via rescaling of the nonlinearities
c.... 5 - Multiple extractions via emittance blow-up
c.... 6 - Reversibility study
c.... 7 - Multi-turn injection
c.... 8 - Multi-turn injection with varying octupole strength
c.... 9 - Distribution evolution with decapoles
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Dimension Istabp(Max_part),ici(5)
c.......................................................................
      Real hmal,x1coor,x2coor,xgsiz,ygsiz,scft,scfw
      Real Sexg,Octg,Tunx,Tuny,Laps1,Laps2,Nsurv,Nlost
      Real rxvec(Max_part),rxpvec(Max_part)
      Real ryvec(Max_part),rypvec(Max_part)
      Real x1tmp(1),x2tmp(2)
      Real hpos,vpos,vincr
c.......................................................................
      Character*3 ext
      Character*4 cmap
      Character*8 txtstring
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Logical Pcheck
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Momen/xmean(Max_boxes,2),xpmean(Max_boxes,2),
     .             x2(Max_boxes,2),xp2(Max_boxes,2),xxp(Max_boxes,2),
     .             x3(Max_boxes,2),xp3(Max_boxes,2),x2xp(Max_boxes,2),
     .             xxp2(Max_boxes,2),x4(Max_boxes,2),xp4(Max_boxes,2),
     .             x2xp2(Max_boxes,2),xxp3(Max_boxes,2),
     .             x3xp(Max_boxes,2),
     .             emit(Max_boxes,2),beta(Max_boxes,2),
     .             alpha(Max_boxes,2),ellipse(Max_boxes,Max_dim),
     .             surf(Max_boxes,Max_dim),halop1(Max_boxes,2),
     .             halop2(Max_boxes,2),Idist1D(Max_bin,Max_boxes,2)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Distribution/Radius,Sigma,Sigmac,Amean,Bmean,Radv,
     .                    Sigv,Sigvc,Ameav,Bmeav,Ameanc,Bmeanc,
     .                    Ameanvc,Bmeanvc,xsep,Sigman,Sigmah(Max_bord),
     .                    Akick(Max_bord),rxvec,rxpvec,
     .                    ryvec,rypvec,Namp,Nangles,Nampc,Nanglesc
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/Plot/Sexg(Max_dat),Octg(Max_dat),
     .            Tunx(Max_dat),Tuny(Max_dat),Laps1(Max_dat),
     .            Laps2(Max_dat),Nsurv(Max_dat),Nlost(Max_dat),
     .            Icount1,Icount2
      Common/Multext/Ibord(0:Max_bord,2),Ixbordmax,Ncorners,
     .               Iarray(Max_part)
      Common/pawc/hmal(nplot)
c.......................................................................
      Data xpmin,ypmin/2*1d38/
      Data xpmax,ypmax/2*-1d38/
c.......................................................................
      do ipart=1,Max_part
        Istabp(ipart)=1
      enddo
c.......................................................................
      pi2=8d0*datan(1d0)                  !..initialisation of variables
      Iout=0
      Iden=7
      Isym=0
      Icount1=0
      Icount2=0
      Ixbordmax=Max_bin+1
c..............................Setting-up of the parameters for graphics
      Call hlimit(nplot)                 !..declares size of common PAWC
c.......................................................................
      ici(1)=1                                        !..defines colours
      ici(2)=2
      ici(3)=3
      ici(4)=7
      ici(5)=6
c.......................................................................
      cmap(0)='COL '
      cmap(1)='COLZ'
      cmap(2)='    '
      cmap(3)='S   '
c.......................................................................
      Nxz(1)=1    !..Number of div. along x (used for multi-projections)
      Nxz(2)=2
      Nxz(3)=2
      Nxz(4)=2
      Nxz(5)=2
      Nxz(6)=2
      Nxz(7)=1
      Nyz(1)=1    !..Number of div. along y (used for multi-projections)
      Nyz(2)=1
      Nyz(3)=2
      Nyz(4)=2
      Nyz(5)=3
      Nyz(6)=3
      Nyz(7)=3
c.......................................................................
      xgsiz(1)=21.0 !..Size of plot along x (used for multi-projections)
      xgsiz(2)=22.4
      xgsiz(3)=22.4
      xgsiz(4)=22.4
      xgsiz(5)=22.7
      xgsiz(6)=22.7
      xgsiz(7)=21.0
      ygsiz(1)=19.0 !..Size of plot along y (used for multi-projections)
      ygsiz(2)=10.2
      ygsiz(3)=18.4
      ygsiz(4)=18.4
      ygsiz(5)=27.0
      ygsiz(6)=27.0
      ygsiz(7)=27.0
c.......................................................................
      scft(1)=0.90       !..Text sc. factor (used for multi-projections)
      scft(2)=0.48
      scft(3)=0.48 
      scft(4)=0.48 
      scft(5)=0.38 
      scft(6)=0.38 
      scft(7)=0.42 
c.......................................................................
      scfw(1)=1.00   !..H. window sc factor (used for multi-projections)
      scfw(2)=2.00
      scfw(3)=2.00 
      scfw(4)=2.00 
      scfw(5)=2.00 
      scfw(6)=2.00 
      scfw(7)=1.00
c.......................................................................
      chtit(1,1)='X '
      chtit(2,1)='X'''
      chtit(3,1)='Y '
      chtit(4,1)='Y'''
      chtit(1,2)='x   (mm^1/2?)'
      chtit(2,2)='x'' (mm^1/2?)'
      chtit(3,2)='y   (mm^1/2?)'
      chtit(4,2)='y'' (mm^1/2?)'
      chtit(1,3)='x   (mm)'
      chtit(2,3)='x'' (mrad)'
      chtit(3,3)='y   (mm)'
      chtit(4,3)='y'' (mrad)'
c.......................................................................
      write(*,*) 'Reads parameters'
      Call Read_par
c.......................................................................
      if ((Icomp.eq.4).or.(Icomp.eq.5)) then
        Sigman=Sigma        !..special definitions for multi-extractions
        Sigmah(1)=Sigma
        do ipart=1,Ntot
          Iarray(ipart)=ipart
        enddo
      endif
c.......................................................................
      Call Find_script
c.......................................................................
      write(*,*) 'Generates distribution'
      Call Gen_dist
c.......................................................................
      write(*,*) 'Starts computations'
c.......................................................................
      if (Icomp.eq.0) then                     !..Graphics and evolution
c.......................................................................
        iturn=0
        if (Ianim.ge.1) Iout=Iout+1
        ext='ps '
        Call Gen_name(ext)
        open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
        Call Init_graph
        Call Opt_graph(Ianim)
c.......................................................................
        Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')       !..sets zone def.
        do jproj=1,iproj
c.......................................................................
          Call Graph2_book(Ianim,iturn,jproj)
          Call Graph_plot(jproj)
c.......................................................................
          Iplane=1+2*(jproj-1)
          do iprint=1,Ntot                            !..graphics output
            do kdim=1,Idim           !..copies co-ordinates into xtransf
              xtransf(kdim)=x(iprint,kdim)
            enddo
            Call Transf                       !..transforms co-ordinates
            x1tmp(1)=real(xtransf(icpr(1)))
            x2tmp(1)=real(xtransf(icpr(2)))
            Call ipm(1,x1tmp,x2tmp)
          enddo
c.......................................................................
        enddo
        Call hplzon(1,1,1,' ')                       !..resets zone def.
c.......................................................................
        if (Ianim.ge.1) then              !..terminates graphics package
          Call End_graph
          close(Iden)
          Call system(script(Ianim)(iscriptl(Ianim):
     .               iscriptu(Ianim))//fnameo(Idim/2))
          Call system('rm -f '//fnameo(Idim/2))            !..deletes ps
        endif
c.......................................................................
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=0
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
c.......................................................................
              if (Ianim.ge.1) then !..initialisation of graphics package
                Iout=Iout+1
                ext='ps '
                Call Gen_name(ext)
                open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
                Call Init_graph
                Call Opt_graph(Ianim)
              endif
c.......................................................................
              Call Graph2_book(Ianim,iturn,iproj)
              Call Graph_plot(iproj)
c.......................................................................
              do iprint=1,Ntot                        !..graphics output
                if (Istab(iprint).eq.1) then 
                  do kdim=1,Idim     !..copies co-ordinates into xtransf
                    xtransf(kdim)=x(iprint,kdim)
                  enddo
                  Call Transf                 !..transforms co-ordinates
                  x1tmp(1)=real(xtransf(icpr(1)))
                  x2tmp(1)=real(xtransf(icpr(2)))
                  Call ipm(1,x1tmp,x2tmp)
                endif
              enddo 
c.......................................................................
              if (Ianim.ge.1) then        !..terminates graphics package
                Call End_graph
                close(Iden)
                Call system(script(Ianim)(iscriptl(Ianim)
     .                    :iscriptu(Ianim))//fnameo(Idim/2))
                Call system('rm -f '//fnameo(Idim/2))      !..deletes ps
              endif
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
c.......................................................................
          enddo
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext                          
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
c.......................................................................
              if (Ianim.ge.1) then !..initialisation of graphics package
                Iout=Iout+1
                ext='ps '
                Call Gen_name(ext)
                open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
                Call Init_graph
                Call Opt_graph(Ianim)
              endif
c.......................................................................
              Call hplzon(Nxz(iproj),Nyz(iproj),1,' ') !..sets zone def.
              do jproj=1,iproj
c.......................................................................
                Call Graph2_book(Ianim,iturn,jproj)
                Call Graph_plot(jproj)
c.......................................................................
                Iplane=1+2*(jproj-1)
                do iprint=1,Ntot                !..fills graphics output
                  do kdim=1,Idim     !..copies co-ordinates into xtransf
                    xtransf(kdim)=x(iprint,kdim)
                  enddo
                  Call Transf                 !..transforms co-ordinates
                  x1tmp(1)=real(xtransf(icpr(1)))
                  x2tmp(1)=real(xtransf(icpr(2)))
                  Call ipm(1,x1tmp,x2tmp)
                enddo 
c.......................................................................
              enddo
              Call hplzon(1,1,1,' ')                 !..resets zone def.
c.......................................................................
              if (Ianim.ge.1) then        !..terminates graphics package
                Call End_graph
                close(Iden)
                Call system(script(Ianim)(iscriptl(Ianim)
     .                    :iscriptu(Ianim))//fnameo(Idim/2))
                Call system('rm -f '//fnameo(Idim/2))      !..deletes ps
              endif
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  x(ipart,3)=0d0
                  x(ipart,4)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
        endif
c.......................................................................
      elseif (Icomp.eq.1) then               !..Particles identification
c.......................................................................
        ext='ps '
        Call Gen_name(ext)
        open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
c.......................................................................
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=0
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
          enddo
c.......................................................................
          tproj=iproj
          iproj=5                        !..resets number of projections
c.......................................................................
          do jbox=1,Nboxes(1)       !..defines boxes to look for islands
c.......................................................................
            bordc1=0.5d0*(bound(1,jbox+1)+bound(3,jbox+1))
            bordc2=0.5d0*(bound(2,jbox+1)+bound(4,jbox+1))
c.......................................................................
            bordw1=0.5d0*(bound(3,jbox+1)-bound(1,jbox+1))
            bordw2=0.5d0*(bound(4,jbox+1)-bound(2,jbox+1))
c.......................................................................
            bord(1,1,jbox+1)=bordc1-bordw1
            bord(2,1,jbox+1)=bordc2-bordw2
            bord(3,1,jbox+1)=bordc1+bordw1
            bord(4,1,jbox+1)=bordc2+bordw2
          enddo
c.......................................................................
          tmpbbox1=bbox(1)                      !..stores initial b. box
          tmpbbox2=bbox(2)
          tmpbbox3=bbox(3)
          tmpbbox4=bbox(4)
c.......................................................................
          bound(1,1)=-xsep                           !..redefines b. box
          bound(2,1)=-xsep
          bound(3,1)=xsep
          bound(4,1)=xsep
          bbox(1)=-xsep                         
          bbox(2)=-xsep
          bbox(3)=xsep
          bbox(4)=xsep
          Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin)   !..step s. C-S
c.......................................................................
          do ipoint=1,4      !..transforms b. box in proper co-ordinates
            anew=ipoint-1.9       !..fancy def. for proper coor. couples
            ianew=anew/2          !..fancy def. for proper coor. couples
            Isgn=(anew+Abs(anew))/anew-1                !..sign function
            Ivar1=ianew-(Isgn-1)/2        !..fancy def. for proper coor. 
            xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
            Ivar2=ipoint/3+1
            xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
            xtransf(3)=xtransf(1)                               
            xtransf(4)=xtransf(2)                               
            Call Transf             !..transforms corners of initial box
            xpmin=min(xpmin,xtransf(2))    !..finds min angle in H plane
            xpmax=max(xpmax,xtransf(2))    !..finds max angle in H plane
            ypmin=min(ypmin,xtransf(4))    !..finds min angle in V plane
            ypmax=max(ypmax,xtransf(4))    !..finds max angle in V plane
          enddo
c.......................................................................
          frame(1,1)=Tmat(1,1,1)*bound(1,1)  !..defines min frame limits
          frame(1,2)=xpmin
          frame(1,3)=Tmat(1,1,2)*bound(1,1)
          frame(1,4)=ypmin
c.......................................................................
          frame(2,1)=Tmat(1,1,1)*bound(3,1)  !..defines max frame limits
          frame(2,2)=xpmax
          frame(2,3)=Tmat(1,1,2)*bound(3,1)
          frame(2,4)=ypmax
c.......................................................................
          Call Init_graph
          Call Opt_graph(Ianim)
          Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')     !..sets zone def.
          iproj=1                        !..resets number of projections
c.......................................................................
          if (Idist.eq.1) then
c.......................................................................
            Rstep=Radius/Dfloat(Namp)
            Astep=pi2/Dfloat(Nangles)
c.......................................................................
            Call Graph2_book(Ianim,Nturn,iproj) !..identifies part. ori.
            Call Graph_plot(iproj)
            do iprint=1,Ntot                          !..graphics output
              if (Istab(iprint).eq.1) then
c.......................................................................
                do jbox=1,Nboxes(1)
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
                    Icli=jbox
                    goto 10
                  endif
c.......................................................................
                enddo
c.......................................................................
 10             Ind1=iprint/Nangles
                Ind2=iprint-Ind1*Nangles
                Amp=Rstep*Dfloat(Ind1)
                Ang=Astep*Dfloat(Ind2)
                xtransf(1)=Amp*Dcos(Ang)+Amean         !..restores init.
                xtransf(2)=Amp*Dsin(Ang)+Bmean
                Call Transf                   !..transforms co-ordinates
                x1tmp(1)=real(xtransf(icpr(1)))
                x2tmp(1)=real(xtransf(icpr(2)))
                Call ispmci(ici(Icli))
                Call ipm(1,x1tmp,x2tmp)
c.......................................................................
              endif
            enddo
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              hpos=0.45*real(bbox(3))
              vpos=real(bbox(4))
              vincr=jbox*0.3*real(bbox(4)-bbox(2))/Nboxes(1)
              Call istxci(ici(jbox))
              write(txtstring(1:7),'(a7)') 'Island '
              write(txtstring(8:8),'(i1)') Ibox(jbox)
              Call itx(hpos,vpos-vincr,txtstring)
c.......................................................................
            enddo
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              Call istxci(ici(1))                       !..resets colour
              Call Graph2_book(Ianim,Nturn,iproj)    !..identifies part.
              Call Graph_plot(iproj)
c.......................................................................
              do iprint=1,Ntot                        !..graphics output
                if (Istab(iprint).eq.1) then
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
c.......................................................................
                    Ind1=iprint/Nangles
                    Ind2=iprint-Ind1*Nangles
                    Amp=Rstep*Dfloat(Ind1)
                    Ang=Astep*Dfloat(Ind2)
c.......................................................................
                    xtransf(1)=Amp*Dcos(Ang)+Amean     !..restores init.
                    xtransf(2)=Amp*Dsin(Ang)+Bmean
                    Call Transf               !..transforms co-ordinates
                    x1tmp(1)=real(xtransf(icpr(1)))
                    x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                    Call ispmci(ici(jbox))
                    Call ipm(1,x1tmp,x2tmp)
c.......................................................................
                  endif
c.......................................................................
                endif
c.......................................................................
              enddo
c.......................................................................
            enddo
c.......................................................................
          elseif (Idist.eq.2) then
c.......................................................................
            Call Graph2_book(Ianim,Nturn,iproj) !..identifies part. ori.
            Call Graph_plot(iproj)
c.......................................................................
            do iprint=1,Ntot                          !..graphics output
              if (Istab(iprint).eq.1) then
c.......................................................................
                do jbox=1,Nboxes(1)
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
                    Icli=jbox
                    goto 20
                  endif
c.......................................................................
                enddo
c.......................................................................
 20             xtransf(1)=Sigma*rxvec(iprint)+Amean   !..restores init.
                xtransf(2)=Sigma*rxpvec(iprint)+Bmean
                Call Transf                   !..transforms co-ordinates
                x1tmp(1)=real(xtransf(icpr(1)))
                x2tmp(1)=real(xtransf(icpr(2)))
                Call ispmci(ici(Icli))
                Call ipm(1,x1tmp,x2tmp)
c.......................................................................
              endif
            enddo 
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              hpos=0.45*real(bbox(3))
              vpos=real(bbox(4))
              vincr=jbox*0.3*real(bbox(4)-bbox(2))/Nboxes(1)
              Call istxci(ici(jbox))
              write(txtstring(1:7),'(a7)') 'Island '
              write(txtstring(8:8),'(i1)') Ibox(jbox)
              Call itx(hpos,vpos-vincr,txtstring)
c.......................................................................
            enddo
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              Call istxci(ici(1))                       !..resets colour
              Call Graph2_book(Ianim,Nturn,iproj)    !..identifies part.
              Call Graph_plot(iproj)
c.......................................................................
              do iprint=1,Ntot                        !..graphics output
                if (Istab(iprint).eq.1) then
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
c.......................................................................
                    xtransf(1)=Sigma*rxvec(iprint)+Amean     !..restores
                    xtransf(2)=Sigma*rxpvec(iprint)+Bmean
                    Call Transf               !..transforms co-ordinates
                    x1tmp(1)=real(xtransf(icpr(1)))
                    x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                    Call ispmci(ici(jbox))
                    Call ipm(1,x1tmp,x2tmp)
c.......................................................................
                  endif
c.......................................................................
                endif
c.......................................................................
              enddo
c.......................................................................
            enddo
c.......................................................................
          elseif (Idist.ge.3) then
c.......................................................................
            Call Graph2_book(Ianim,Nturn,iproj) !..identifies part. ori.
            Call Graph_plot(iproj)
c.......................................................................
            do iprint=1,Ntot                          !..graphics output
              if (Istab(iprint).eq.1) then
c.......................................................................
                do jbox=1,Nboxes(1)
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
                    Icli=jbox
                    goto 30
                  endif
c.......................................................................
                enddo
c.........................................................restores init.
 30             xtransf(1)=(Sigma*rxvec(iprint)+Amean)*
     .                                           Dcos(pi2*rxpvec(ipart))   
                xtransf(2)=(Sigma*rxvec(iprint)+Bmean)*
     .                                           Dsin(pi2*rxpvec(ipart))
                Call Transf                   !..transforms co-ordinates
                x1tmp(1)=real(xtransf(icpr(1)))
                x2tmp(1)=real(xtransf(icpr(2)))
                Call ispmci(ici(Icli))
                Call ipm(1,x1tmp,x2tmp)
c.......................................................................
              endif
            enddo 
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              hpos=0.45*real(bbox(3))
              vpos=real(bbox(4))
              vincr=jbox*0.3*real(bbox(4)-bbox(2))/Nboxes(1)
              Call istxci(ici(jbox))
              write(txtstring(1:7),'(a7)') 'Island '
              write(txtstring(8:8),'(i1)') Ibox(jbox)
              Call itx(hpos,vpos-vincr,txtstring)
c.......................................................................
            enddo
c.......................................................................
            do jbox=1,Nboxes(1)
c.......................................................................
              Call istxci(ici(1))                       !..resets colour
              Call Graph2_book(Ianim,Nturn,iproj)    !..identifies part.
              Call Graph_plot(iproj)
c.......................................................................
              do iprint=1,Ntot                        !..graphics output
                if (Istab(iprint).eq.1) then
c.......................................................................
                  Pcheck=(bord(1,1,jbox+1).le.x(iprint,1)).and.
     .                   (bord(3,1,jbox+1).gt.x(iprint,1)).and.
     .                   (bord(2,1,jbox+1).le.x(iprint,2)).and.
     .                   (bord(4,1,jbox+1).gt.x(iprint,2))
                  if (Pcheck) then
c.........................................................restores init.
                    xtransf(1)=(Sigma*rxvec(iprint)+Amean)*
     .                                           Dcos(pi2*rxpvec(ipart))   
                    xtransf(2)=(Sigma*rxvec(iprint)+Bmean)*
     .                                           Dsin(pi2*rxpvec(ipart))
                    Call Transf               !..transforms co-ordinates
                    x1tmp(1)=real(xtransf(icpr(1)))
                    x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                    Call ispmci(ici(jbox))
                    Call ipm(1,x1tmp,x2tmp)
c.......................................................................
                  endif
c.......................................................................
                endif
c.......................................................................
              enddo
c.......................................................................
            enddo
c.......................................................................
          endif
c.......................................................................
          bound(1,1)=tmpbbox1                         !..restores b. box
          bound(2,1)=tmpbbox2
          bound(3,1)=tmpbbox3
          bound(4,1)=tmpbbox4
          bbox(1)=tmpbbox1                      
          bbox(2)=tmpbbox2
          bbox(3)=tmpbbox3
          bbox(4)=tmpbbox4
          Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin)   !..step s. C-S
c.......................................................................
          do ipoint=1,4      !..transforms b. box in proper co-ordinates
            anew=ipoint-1.9       !..fancy def. for proper coor. couples
            ianew=anew/2          !..fancy def. for proper coor. couples
            Isgn=(anew+Abs(anew))/anew-1                !..sign function
            Ivar1=ianew-(Isgn-1)/2        !..fancy def. for proper coor. 
            xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
            Ivar2=ipoint/3+1
            xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
            xtransf(3)=xtransf(1)                               
            xtransf(4)=xtransf(2)                               
            Call Transf             !..transforms corners of initial box
            xpmin=min(xpmin,xtransf(2))    !..finds min angle in H plane
            xpmax=max(xpmax,xtransf(2))    !..finds max angle in H plane
            ypmin=min(ypmin,xtransf(4))    !..finds min angle in V plane
            ypmax=max(ypmax,xtransf(4))    !..finds max angle in V plane
          enddo
c.......................................................................
          frame(1,1)=Tmat(1,1,1)*bound(1,1)    !..defines min frame lim.
          frame(1,2)=xpmin
          frame(1,3)=Tmat(1,1,2)*bound(1,1)
          frame(1,4)=ypmin
c.......................................................................
          frame(2,1)=Tmat(1,1,1)*bound(3,1)    !..defines max frame lim.
          frame(2,2)=xpmax
          frame(2,3)=Tmat(1,1,2)*bound(3,1)
          frame(2,4)=ypmax
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
        endif
c.......................................................................
        iproj=tproj                    !..restores number of projections
c.......................................................................
      elseif (Icomp.eq.2) then      !..computations of beam distribution
c.......................................................................
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=0
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn
            Isurv=0
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
          enddo
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn
            Isurv=0
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
               Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        iflag=1
        Call Comp_dist(iflag,iturn)    !..computes/prints i/e beam dist.
c.......................................................................
      elseif (Icomp.eq.3) then
c.......................................................................
        if (Ianim.eq.0) then       !..initialisation of graphics package
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        if (Idump.eq.1) then  !..initialisation of file for dumping data
          ext='dat'
          Call Gen_name(ext)
          open(Iden+2,file=fnameo(9),form='formatted',status='unknown')
        endif
c.......................................................................
        iturn=0
        iflag=2
        Call Comp_dist(iflag,iturn)     !..computes 2D beam distribution
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
              Call Comp_dist(iflag,iturn)      !..computes 2D beam dist.
c.......................................................................
              if ((Idump.eq.1).and.(iturn.eq.Nturn)) 
     .                         then  !..dumps final co-ordinates on file
c.......................................................................
                do ipart=1,Ntot
                  if (Istab(ipart).eq.1) write(Iden+2,'(4(1x,e12.6))') 
     .                (x(ipart,irow),irow=1,Idim) !..writes co-ordinates 
                enddo
c.......................................................................
                close(Iden+2)
c.......................................................................
              endif
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
          enddo
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
              Call Comp_dist(iflag,iturn)      !..computes 4D beam dist.
c.......................................................................
              if ((Idump.eq.1).and.(iturn.eq.Nturn)) 
     .                         then  !..dumps final co-ordinates on file
c.......................................................................
                do ipart=1,Ntot
                  if (Istab(ipart).eq.1) write(Iden+2,'(4(1x,e12.6))') 
     .                (x(ipart,irow),irow=1,Idim) !..writes co-ordinates 
                enddo
c.......................................................................
                close(Iden+2)
c.......................................................................
              endif
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  x(ipart,3)=0d0
                  x(ipart,4)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
        endif
c.......................................................................
        ext='ps '
        Call Gen_name(ext)
        Call system('mv -f fort.7 '//fnameo(Idim/2))
c.......................................................................
      elseif ((Icomp.eq.4).or.(Icomp.eq.5)) then
c.......................................................................
        Ikick=Icomp-4   !..switch to blow-up with kick or nonlinearities
c.......................................................................
        if (Ianim.eq.0) then       !..initialisation of graphics package
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        iturn=0
        iflag=3
        Call Comp_dist(iflag,iturn)     !..computes 2D beam distribution
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do imul=1,Mult                         !..multiple extractions
            do iturn=1,Nturn
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
                Itest=Iuser_map2(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=iturn-1+Iper*(imul-1)
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=3
            Call Comp_dist(iflag,iturn-1+Iper*(imul-1)) !..2D beam dist.
c.......................................................................
            Sigmah(imul+1)=Sigman               !..new sigma becomes old
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,iturn+3+Iper*(imul-1)) !..2D beam dist.
c.......................................................................
            do iturn=ipar(5),ipar(6)
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
                Itest=Iuser_map2(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount2=Icount2+1                !..writes tunes/nlin. grad.
            Laps2(Icount2)=iturn-1+Iper*(imul-1)
            Tunx(Icount2)=Real(par(24)/pi2)
            Tuny(Icount2)=Real(par(23)/pi2)
            Sexg(Icount2)=Real(par(22))
            Octg(Icount2)=Real(par(21))
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=iturn-1+Iper*(imul-1)
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,iturn-1+Iper*(imul-1))       !..2D beam
c.......................................................................
            if (Ikick.eq.0) then
              Slope=par(4)*(Sigmah(imul)/Sigmah(imul+1)-1d0)/
     .                                                 (ipar(7)-ipar(6))
              Akick(imul)=0d0
            else
              Slope=0d0
              Akick(imul)=0.5d0*Dsqrt(Sigmah(1)*Sigmah(1)-
     .                                Sigmah(imul+1)*Sigmah(imul+1)/
     .                                          Dfloat(ipar(7)-ipar(6)))
            endif
            Vinit=par(4)
            do iturn=ipar(6)+1,ipar(7)
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
c.......................................................................
                par(4)=Vinit+Slope*Dfloat(iturn-ipar(6))        !..ramps
c.......................................................................
                x(Iarray(ipart),2)=x(Iarray(ipart),2)+Akick(imul)*Ikick
c.......................................................................
                Itest=Iuser_map2(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount2=Icount2+1                !..writes tunes/nlin. grad.
            Laps2(Icount2)=Iper*imul
            Tunx(Icount2)=Real(par(24)/pi2)
            Tuny(Icount2)=Real(par(23)/pi2)
            Sexg(Icount2)=Real(par(22))
            Octg(Icount2)=Real(par(21))
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=Iper*imul
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,Iper*imul)             !..2D beam dist.
c.......................................................................
          enddo
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do imul=1,Mult                         !..multiple extractions
            do iturn=1,Nturn
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
                Itest=Iuser_map4(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=iturn-1+Iper*(imul-1)
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=3
            Call Comp_dist(iflag,iturn-1+Iper*(imul-1)) !..4D beam dist.
c.......................................................................
            Sigmah(imul+1)=Sigman               !..new sigma becomes old
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,iturn+3+Iper*(imul-1)) !..4D beam dist.
c.......................................................................
            do iturn=ipar(5),ipar(6)
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
                Itest=Iuser_map4(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount2=Icount2+1                !..writes tunes/nlin. grad.
            Laps2(Icount2)=iturn-1+Iper*(imul-1)
            Tunx(Icount2)=Real(par(24)/pi2)
            Tuny(Icount2)=Real(par(23)/pi2)
            Sexg(Icount2)=Real(par(22))
            Octg(Icount2)=Real(par(21))
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=iturn-1+Iper*(imul-1)
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,iturn-1+Iper*(imul-1))       !..4D beam
c.......................................................................
            if (Ikick.eq.0) then
              Slope=par(4)*(Sigmah(imul)/Sigmah(imul+1)-1d0)/
     .                                                 (ipar(7)-ipar(6))
              Akick(imul)=0d0
            else
              Slope=0d0
              Akick(imul)=0.5d0*Dsqrt(Sigmah(1)*Sigmah(1)-
     .                                Sigmah(imul+1)*Sigmah(imul+1)/
     .                                          Dfloat(ipar(7)-ipar(6)))
            endif
            Vinit=par(4)
            do iturn=ipar(6)+1,ipar(7)
c.......................................................................
              Isurv=0
c.......................................................................
              do ipart=1,Npmult
c.......................................................................
                par(4)=Vinit+Slope*Dfloat(iturn-ipar(6))     !..ramps up
c.......................................................................
                x(Iarray(ipart),2)=x(Iarray(ipart),2)+Akick(imul)*Ikick
c.......................................................................
                Itest=Iuser_map2(Iarray(ipart),iturn)       !..iteration
                Istab(Iarray(ipart))=Istabp(Iarray(ipart))-Itest !..lost
                Isurv=Isurv+Istab(Iarray(ipart))   !..counts surv. part.
                Istabp(Iarray(ipart))=Istab(Iarray(ipart))
              enddo
              if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
                Icount1=Icount1+1                    !..writes rem./lost
                Laps1(Icount1)=iturn+Iper*(imul-1)
                Nsurv(Icount1)=Isurv
                Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
                Icount2=Icount2+1            !..writes tunes/nlin. grad.
                Laps2(Icount2)=iturn+Iper*(imul-1)
                Tunx(Icount2)=Real(par(24)/pi2)
                Tuny(Icount2)=Real(par(23)/pi2)
                Sexg(Icount2)=Real(par(22))
                Octg(Icount2)=Real(par(21))
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            Icount2=Icount2+1                !..writes tunes/nlin. grad.
            Laps2(Icount2)=Iper*imul
            Tunx(Icount2)=Real(par(24)/pi2)
            Tuny(Icount2)=Real(par(23)/pi2)
            Sexg(Icount2)=Real(par(22))
            Octg(Icount2)=Real(par(21))
c.......................................................................
            Icount1=Icount1+1                        !..writes rem./lost
            Laps1(Icount1)=Iper*imul
            Nsurv(Icount1)=Isurv
            Nlost(Icount1)=(Npmult-Isurv)
c.......................................................................
            iflag=4
            Call Comp_dist(iflag,Iper*imul)             !..4D beam dist.
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
        endif
c.......................................................................
      elseif (Icomp.eq.6) then
c.......................................................................
        if (Ianim.eq.0) then       !..initialisation of graphics package
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        iturn=0
        iflag=2
        Call Comp_dist(iflag,iturn)     !..computes 2D beam distribution
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
            if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
          Icount1=Icount1+1                          !..writes rem./lost
          Laps1(Icount1)=iturn-1
          Nsurv(Icount1)=Isurv
          Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
          Call Comp_dist(iflag,iturn-1)          !..2D beam distribution
c.......................................................................
          do iturn=ipar(5),ipar(6)
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
            if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
          Icount2=Icount2+1                  !..writes tunes/nlin. grad.
          Laps2(Icount2)=iturn-1
          Tunx(Icount2)=Real(par(24)/pi2)
          Tuny(Icount2)=Real(par(23)/pi2)
          Sexg(Icount2)=Real(par(22))
          Octg(Icount2)=Real(par(21))
c.......................................................................
          Icount1=Icount1+1                          !..writes rem./lost
          Laps1(Icount1)=iturn-1
          Nsurv(Icount1)=Isurv
          Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
          Call Comp_dist(iflag,iturn-1)          !..2D beam distribution
          Call End_graph                      !..closes graphics package
          close(Iden)
c.......................................................................
          iflag=0
          Call Comp_dist(iflag,iturn-1) !..computes/prints i/e beam dist
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
            if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
          Icount1=Icount1+1                          !..writes rem./lost
          Laps1(Icount1)=iturn
          Nsurv(Icount1)=Isurv
          Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
          Call Comp_dist(iflag,iturn-1)          !..4D beam distribution
c.......................................................................
          do iturn=ipar(5),ipar(6)
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
            if (mod(iturn,2*ipar(9)).eq.0) then 
c.......................................................................
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
          Icount2=Icount2+1                  !..writes tunes/nlin. grad.
          Laps2(Icount2)=iturn-1
          Tunx(Icount2)=Real(par(24)/pi2)
          Tuny(Icount2)=Real(par(23)/pi2)
          Sexg(Icount2)=Real(par(22))
          Octg(Icount2)=Real(par(21))
c.......................................................................
          Icount1=Icount1+1                          !..writes rem./lost
          Laps1(Icount1)=iturn-1
          Nsurv(Icount1)=Isurv
          Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
          Call Comp_dist(iflag,iturn-1)          !..4D beam distribution
          Call End_graph                      !..closes graphics package
          close(Iden)
c.......................................................................
          iflag=0
          Call Comp_dist(iflag,iturn-1) !..computes/prints i/e beam dist
c.......................................................................
        endif
c.......................................................................
      elseif (Icomp.eq.7) then
c.......................................................................
        iturn=Ninj
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
        iflag=2
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn-(Ninj+Icentre)
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Ninj+Icentre+Noff,Nout).eq.0) then   !..checks 
              Call Comp_dist(iflag,iturn+Ninj+Icentre+Noff)   !..comp 2D
c.......................................................................
            endif
            if (mod(iturn+Ninj+Icentre,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn+Ninj+Icentre
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn+Ninj+Icentre
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
          enddo
c.......................................................................
          if (Ianim.eq.0) then
c.......................................................................
            tmpbbox1=bbox(1)                    !..stores initial b. box
            tmpbbox2=bbox(2)
            tmpbbox3=bbox(3)
            tmpbbox4=bbox(4)
c.......................................................................
            bound(1,1)=-xsep                         !..redefines b. box
            bound(2,1)=-xsep
            bound(3,1)=xsep
            bound(4,1)=xsep
            bbox(1)=-xsep                         
            bbox(2)=-xsep
            bbox(3)=xsep
            bbox(4)=xsep
            Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin) !..step s. C-S
c.......................................................................
            do ipoint=1,4    !..transforms b. box in proper co-ordinates
              anew=ipoint-1.9     !..fancy def. for proper coor. couples
              ianew=anew/2        !..fancy def. for proper coor. couples
              Isgn=(anew+Abs(anew))/anew-1              !..sign function
              Ivar1=ianew-(Isgn-1)/2      !..fancy def. for proper coor. 
              xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
              Ivar2=ipoint/3+1
              xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
              xtransf(3)=xtransf(1)                               
              xtransf(4)=xtransf(2)                               
              Call Transf           !..transforms corners of initial box
              xpmin=min(xpmin,xtransf(2))  !..finds min angle in H plane
              xpmax=max(xpmax,xtransf(2))  !..finds max angle in H plane
              ypmin=min(ypmin,xtransf(4))  !..finds min angle in V plane
              ypmax=max(ypmax,xtransf(4))  !..finds max angle in V plane
            enddo
c.......................................................................
            frame(1,1)=Tmat(1,1,1)*bound(1,1)!..defines min frame limits
            frame(1,2)=xpmin
            frame(1,3)=Tmat(1,1,2)*bound(1,1)
            frame(1,4)=ypmin
c.......................................................................
            frame(2,1)=Tmat(1,1,1)*bound(3,1)!..defines max frame limits
            frame(2,2)=xpmax
            frame(2,3)=Tmat(1,1,2)*bound(3,1)
            frame(2,4)=ypmax
c.......................................................................
            tproj=iproj
            iproj=Ninj+1+Icentre         !..resets number of projections
            Call Opt_graph(Ianim)
            Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')   !..sets zone def.
            iproj=1                      !..resets number of projections
c.......................................................................
            Call Graph2_book(Ianim,Nturn,iproj)      !..identifies part.
            Call Graph_plot(iproj)
            do iprint=1,Ntot                          !..graphics output
              if (Istab(iprint).eq.1) then
c.......................................................................
                xtransf(1)=x(iprint,1)            !..copies co-ordinates
                xtransf(2)=x(iprint,2)
                Call Transf                   !..transforms co-ordinates
                x1tmp(1)=real(xtransf(icpr(1)))
                x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                if ((Icentre.eq.0).or.(iprint.le.Namp*Nangles)) then
                  Icli=iprint/(Namp*Nangles)+1
                else
                  Icli=Ninj+Icentre
                endif
c.......................................................................
                Call ispmci(ici(Icli))
                Call ipm(1,x1tmp,x2tmp)
c.......................................................................
              endif
            enddo
c.......................................................................
            do jbox=1,Ninj+Icentre
c.......................................................................
              hpos=0.45*real(bbox(3))
              vpos=real(bbox(4))
              vincr=jbox*0.3*real(bbox(4)-bbox(2))/(Ninj+Icentre)
              Call istxci(ici(jbox))
              write(txtstring(1:7),'(a7)') '  Turn '
              write(txtstring(8:8),'(i1)') jbox
              Call itx(hpos,vpos-vincr,txtstring)
c.......................................................................
            enddo
c.......................................................................
            do jbox=1,Ninj+Icentre
c.......................................................................
              Call istxci(ici(1))                       !..resets colour
              Call Graph2_book(Ianim,Nturn,iproj)    !..identifies part.
              Call Graph_plot(iproj)
c.......................................................................
              iextremum=Namp*Nangles*(1-jbox/(Ninj+Icentre))+
     .                  Nampc*Nanglesc*jbox/(Ninj+Icentre)
              do iprn=1,iextremum                     !..graphics output
                iprint=(jbox-1)*Namp*Nangles+iprn
                if (Istab(iprint).eq.1) then
c.......................................................................
                  xtransf(1)=x(iprint,1)          !..copies co-ordinates
                  xtransf(2)=x(iprint,2)
                  Call Transf                 !..transforms co-ordinates
                  x1tmp(1)=real(xtransf(icpr(1)))
                  x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                  Call ispmci(ici(jbox))
                  Call ipm(1,x1tmp,x2tmp)
c.......................................................................
                endif
c.......................................................................
              enddo
c.......................................................................
            enddo
c.......................................................................
            iproj=tproj                !..restores number of projections
c.......................................................................
            xpmin=1d38
            xpmax=-1d38
            ypmin=1d38
            ypmax=-1d38
            bound(1,1)=tmpbbox1                       !..restores b. box
            bound(2,1)=tmpbbox2
            bound(3,1)=tmpbbox3
            bound(4,1)=tmpbbox4
            bbox(1)=tmpbbox1                      
            bbox(2)=tmpbbox2
            bbox(3)=tmpbbox3
            bbox(4)=tmpbbox4
            Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin) !..step s. C-S
c.......................................................................
            do ipoint=1,4    !..transforms b. box in proper co-ordinates
              anew=ipoint-1.9     !..fancy def. for proper coor. couples
              ianew=anew/2        !..fancy def. for proper coor. couples
              Isgn=(anew+Abs(anew))/anew-1              !..sign function
              Ivar1=ianew-(Isgn-1)/2      !..fancy def. for proper coor. 
              xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
              Ivar2=ipoint/3+1
              xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
              xtransf(3)=xtransf(1)                               
              xtransf(4)=xtransf(2)                               
              Call Transf           !..transforms corners of initial box
              xpmin=min(xpmin,xtransf(2))  !..finds min angle in H plane
              xpmax=max(xpmax,xtransf(2))  !..finds max angle in H plane
              ypmin=min(ypmin,xtransf(4))  !..finds min angle in V plane
              ypmax=max(ypmax,xtransf(4))  !..finds max angle in V plane
            enddo
c.......................................................................
            frame(1,1)=Tmat(1,1,1)*bound(1,1)  !..defines min frame lim.
            frame(1,2)=xpmin
            frame(1,3)=Tmat(1,1,2)*bound(1,1)
            frame(1,4)=ypmin
c.......................................................................
            frame(2,1)=Tmat(1,1,1)*bound(3,1)  !..defines max frame lim.
            frame(2,2)=xpmax
            frame(2,3)=Tmat(1,1,2)*bound(3,1)
            frame(2,4)=ypmax
c.......................................................................
          endif
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn-(Ninj+Icentre)
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Ninj+Icentre+Noff,Nout).eq.0) then   !..checks
              Call Comp_dist(iflag,iturn+Ninj+Icentre+Noff)  !..comp. 2D
c.......................................................................
            endif
            if (mod(iturn+Ninj+Icentre,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn+Ninj+Icentre
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn+Ninj+Icentre
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
c.......................................................................
          iflag=0
          Call Comp_dist(iflag,iturn)  !..computes/prints i/e beam dist.
        endif
c.......................................................................
      elseif (Icomp.eq.8) then
c.......................................................................
        iturn=Ninj
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
        iflag=2
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn-Ninj
c.......................................................................
            if ((iturn.le.ipar(4))) then                 !..Changes oct.
              Xctk_ft=par(37)
              Octk_ft=par(27)
            elseif ((iturn.gt.ipar(4)).and.(iturn.le.ipar(1))) then 
              Xctk_ft=par(37)+par(35)*
     .                       (Dexp(par(36)*Dlog(Dfloat(iturn)))-1d0)
              Octk_ft=par(27)+par(25)*
     .                       (Dexp(par(26)*Dlog(Dfloat(iturn)))-1d0)
            elseif ((iturn.gt.ipar(1))) then
              Xctk_ft=0d0
              Octk_ft=0d0
            endif
            par(34)=Xctk_ft
            par(4)=Octk_ft
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Ninj+Noff,Nout).eq.0) then !..checks for print
              Call Comp_dist(iflag,iturn+Ninj+Noff) !..computes 2D dist.
c.......................................................................
            endif
            if (mod(iturn+Ninj,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn+Ninj
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn+Ninj
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
          enddo
c.......................................................................
          if (Ianim.eq.0) then
c.......................................................................
            tmpbbox1=bbox(1)                    !..stores initial b. box
            tmpbbox2=bbox(2)
            tmpbbox3=bbox(3)
            tmpbbox4=bbox(4)
c.......................................................................
            bound(1,1)=-xsep                         !..redefines b. box
            bound(2,1)=-xsep
            bound(3,1)=xsep
            bound(4,1)=xsep
            bbox(1)=-xsep                         
            bbox(2)=-xsep
            bbox(3)=xsep
            bbox(4)=xsep
            Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin) !..step s. C-S
c.......................................................................
            do ipoint=1,4    !..transforms b. box in proper co-ordinates
              anew=ipoint-1.9     !..fancy def. for proper coor. couples
              ianew=anew/2        !..fancy def. for proper coor. couples
              Isgn=(anew+Abs(anew))/anew-1              !..sign function
              Ivar1=ianew-(Isgn-1)/2      !..fancy def. for proper coor. 
              xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
              Ivar2=ipoint/3+1
              xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
              xtransf(3)=xtransf(1)                               
              xtransf(4)=xtransf(2)                               
              Call Transf           !..transforms corners of initial box
              xpmin=min(xpmin,xtransf(2))  !..finds min angle in H plane
              xpmax=max(xpmax,xtransf(2))  !..finds max angle in H plane
              ypmin=min(ypmin,xtransf(4))  !..finds min angle in V plane
              ypmax=max(ypmax,xtransf(4))  !..finds max angle in V plane
            enddo
c.......................................................................
            frame(1,1)=Tmat(1,1,1)*bound(1,1)!..defines min frame limits
            frame(1,2)=xpmin
            frame(1,3)=Tmat(1,1,2)*bound(1,1)
            frame(1,4)=ypmin
c.......................................................................
            frame(2,1)=Tmat(1,1,1)*bound(3,1)!..defines max frame limits
            frame(2,2)=xpmax
            frame(2,3)=Tmat(1,1,2)*bound(3,1)
            frame(2,4)=ypmax
c.......................................................................
            tproj=iproj
            iproj=Ninj+1                 !..resets number of projections
            Call Opt_graph(Ianim)
            Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')   !..sets zone def.
            iproj=1                      !..resets number of projections
c.......................................................................
            Call Graph2_book(Ianim,Nturn,iproj)      !..identifies part.
            Call Graph_plot(iproj)
            do iprint=1,Ntot                          !..graphics output
              if (Istab(iprint).eq.1) then
c.......................................................................
                xtransf(1)=x(iprint,1)            !..copies co-ordinates
                xtransf(2)=x(iprint,2)
                Call Transf                   !..transforms co-ordinates
                x1tmp(1)=real(xtransf(icpr(1)))
                x2tmp(1)=real(xtransf(icpr(2)))
                Icli=iprint/(Namp*Nangles)+1
                Call ispmci(ici(Icli))
                Call ipm(1,x1tmp,x2tmp)
c.......................................................................
              endif
            enddo
c.......................................................................
            do jbox=1,Ninj
c.......................................................................
              hpos=0.45*real(bbox(3))
              vpos=real(bbox(4))
              vincr=jbox*0.3*real(bbox(4)-bbox(2))/Ninj
              Call istxci(ici(jbox))
              write(txtstring(1:7),'(a7)') '  Turn '
              write(txtstring(8:8),'(i1)') jbox
              Call itx(hpos,vpos-vincr,txtstring)
c.......................................................................
            enddo
c.......................................................................
            do jbox=1,Ninj
c.......................................................................
              Call istxci(ici(1))                       !..resets colour
              Call Graph2_book(Ianim,Nturn,iproj)    !..identifies part.
              Call Graph_plot(iproj)
c.......................................................................
              do iprn=1,Namp*Nangles                  !..graphics output
                iprint=(jbox-1)*Namp*Nangles+iprn
                if (Istab(iprint).eq.1) then
c.......................................................................
                  xtransf(1)=x(iprint,1)          !..copies co-ordinates
                  xtransf(2)=x(iprint,2)
                  Call Transf                 !..transforms co-ordinates
                  x1tmp(1)=real(xtransf(icpr(1)))
                  x2tmp(1)=real(xtransf(icpr(2)))
c.......................................................................
                  Call ispmci(ici(jbox))
                  Call ipm(1,x1tmp,x2tmp)
c.......................................................................
                endif
c.......................................................................
              enddo
c.......................................................................
            enddo
c.......................................................................
            iproj=tproj                !..restores number of projections
c.......................................................................
            xpmin=1d38
            xpmax=-1d38
            ypmin=1d38
            ypmax=-1d38
            bound(1,1)=tmpbbox1                       !..restores b. box
            bound(2,1)=tmpbbox2
            bound(3,1)=tmpbbox3
            bound(4,1)=tmpbbox4
            bbox(1)=tmpbbox1                      
            bbox(2)=tmpbbox2
            bbox(3)=tmpbbox3
            bbox(4)=tmpbbox4
            Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin) !..step s. C-S
c.......................................................................
            do ipoint=1,4    !..transforms b. box in proper co-ordinates
              anew=ipoint-1.9     !..fancy def. for proper coor. couples
              ianew=anew/2        !..fancy def. for proper coor. couples
              Isgn=(anew+Abs(anew))/anew-1              !..sign function
              Ivar1=ianew-(Isgn-1)/2      !..fancy def. for proper coor. 
              xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)
              Ivar2=ipoint/3+1
              xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)
              xtransf(3)=xtransf(1)                               
              xtransf(4)=xtransf(2)                               
              Call Transf           !..transforms corners of initial box
              xpmin=min(xpmin,xtransf(2))  !..finds min angle in H plane
              xpmax=max(xpmax,xtransf(2))  !..finds max angle in H plane
              ypmin=min(ypmin,xtransf(4))  !..finds min angle in V plane
              ypmax=max(ypmax,xtransf(4))  !..finds max angle in V plane
            enddo
c.......................................................................
            frame(1,1)=Tmat(1,1,1)*bound(1,1)  !..defines min frame lim.
            frame(1,2)=xpmin
            frame(1,3)=Tmat(1,1,2)*bound(1,1)
            frame(1,4)=ypmin
c.......................................................................
            frame(2,1)=Tmat(1,1,1)*bound(3,1)  !..defines max frame lim.
            frame(2,2)=xpmax
            frame(2,3)=Tmat(1,1,2)*bound(3,1)
            frame(2,4)=ypmax
c.......................................................................
          endif
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          do iturn=1,Nturn-Ninj
c.......................................................................
            if ((iturn.le.ipar(4))) then                 !..Changes oct.
              Octk_ft=par(27)
            elseif ((iturn.gt.ipar(4)).and.(iturn.le.ipar(1))) then 
              Octk_ft=par(27)+par(25)*Dexp(par(26)*
     .                                Dlog(Dfloat(ipar(1)-iturn)))
            elseif ((iturn.gt.ipar(1))) then
              Octk_ft=0d0
            endif
            par(4)=Octk_ft
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,iturn)                 !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Ninj+Noff,Nout).eq.0) then !..checks for print
              Call Comp_dist(iflag,iturn+Ninj+Nout) !..computes 4D dist.
c.......................................................................
            endif
            if (mod(iturn+Ninj,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn+Ninj
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn+Ninj
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
c.......................................................................
          iflag=0
          Call Comp_dist(iflag,iturn)  !..computes/prints i/e beam dist.
        endif
c.......................................................................
      elseif (Icomp.eq.9) then
c.......................................................................
        if (Ianim.eq.0) then       !..initialisation of graphics package
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        iturn=0
        iflag=2
        Call Comp_dist(iflag,iturn)     !..computes 2D beam distribution
        Icount1=Icount1+1                            !..writes rem./lost
        Laps1(Icount1)=iturn
        Nsurv(Icount1)=Ntot
        Nlost(Icount1)=0
        Icount2=Icount2+1                    !..writes tunes/nlin. grad.
        Laps2(Icount2)=0
        Tunx(Icount2)=Real(par(1)/pi2)
        Tuny(Icount2)=Real(par(14)/pi2)
        Sexg(Icount2)=Real(1d0)
        Octg(Icount2)=Real(par(4))
        if (ipar(4).ne.0) then
          Sexg(Icount2)=Real(0d0)
          Octg(Icount2)=Real(0d0)
        endif
c.......................................................................
        if (Idim.eq.2) then                           !..2D computations
c.......................................................................
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map2_s(ipart,iturn)               !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
              Call Comp_dist(iflag,iturn)      !..computes 2D beam dist.
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
          enddo
c.......................................................................
        elseif (Idim.eq.4) then                       !..4D computations
c.......................................................................
          write(*,*) '*** Option not active yet!'
          write(*,*) '*** Program stops'
          stop
          do iturn=1,Nturn+4*Iext
c.......................................................................
            Isurv=0
c.......................................................................
            do ipart=1,Ntot
              Itest=Iuser_map4_s(ipart,iturn)               !..iteration
              Istab(ipart)=Istabp(ipart)-Itest  !..labels lost particles
              Isurv=Isurv+Istab(ipart)     !..counts surviving particles
              Istabp(ipart)=Istab(ipart)
            enddo
c.......................................................................
            if (mod(iturn+Noff,Nout).eq.0) then   !..checks for printout
              Call Comp_dist(iflag,iturn)      !..computes 4D beam dist.
c.......................................................................
            endif
            if (mod(iturn,2*ipar(9)).eq.0) then 
              Icount1=Icount1+1                      !..writes rem./lost
              Laps1(Icount1)=iturn
              Nsurv(Icount1)=Isurv
              Nlost(Icount1)=(Ntot-Isurv)
c.......................................................................
              Icount2=Icount2+1              !..writes tunes/nlin. grad.
              Laps2(Icount2)=iturn
              Tunx(Icount2)=Real(par(24)/pi2)
              Tuny(Icount2)=Real(par(23)/pi2)
              Sexg(Icount2)=Real(par(22))
              Octg(Icount2)=Real(par(21))
            endif
c.......................................................................
            if (iturn.gt.Nturn) then                !..Starts extraction
              Nout=1          !..Redefines parameter to output each turn
              do ipart=1,Ntot     !..Kills particles beyond septum blade
                if (x(ipart,1).gt.xsep) then
                  x(ipart,1)=0d0
                  x(ipart,2)=0d0
                  x(ipart,3)=0d0
                  x(ipart,4)=0d0
                  Istabp(ipart)=0
                endif
              enddo
            endif
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        if (Ianim.eq.0) then                  !..closes graphics package
          Call End_graph
          close(Iden)
        endif
c.......................................................................
      endif
c.......................................................................
      Call Plot_par(Ierr)   !..generates plot with simulation parameters
c.......................................................................
      if (Ierr.eq.0) Call End_graph       !..terminates graphics package
      close(Iden+10)
c................................................generates ps file names
      ext='ps '
      Call Gen_name(ext)
      Call system('mv -f fort.17 '//fnameo(2+Idim/2))
c....................................generates text file with parameters
      ext='out'
      Call Gen_name(ext)
      open(Iden+11,file=fnameo(7),form='formatted',status='unknown')
      open(Iden+12,file=fnameo(8),form='formatted',status='unknown')
      Call Write_par        !..generates file with simulation parameters
      close(Iden+11)
      close(Iden+12)
c.......................................................................
      End
c.......................................................................

      Integer Function Iuser_map2(ipart,n)
c.......................................................................
c.... Function to compute the evolution of particles. It is a 2D
c.... polynomial map
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
c.....par(1) is the initial Hor. frequency*pi2
c.....par(2) is the intermediate Hor. frequency*pi2
c.....par(3) is the final Hor. frequency*pi2
c.....par(4) is the octupole coefficient
c.....par(5) is the first Hor. slope
c.....par(6) is the second Hor. slope
c.....par(7) is the exponent of the first slope
c.....par(8) is the slope to ramp the sextupole strength
c.....par(9) is the slope to ramp the octupole strength
c.....par(10) is the slope to ramp up to intial tune for multiple 
c.....extractions
c.....par(11) is the beta ratio at sextupole location
c.....par(12) is the 3*beta ratio at the octupole location
c.....par(13) is the square of the beta ratio at the octupole location
c.....par(14) is the initial Ver. frequency*pi2
c.....par(15) is the intermediate Ver. frequency*pi2
c.....par(16) is the final Ver. frequency*pi2
c.....par(17) is the first Ver. slope
c.....par(18) is the second Ver. slope
c.....ipar(1) is the time between par(1) and par(2)
c.....ipar(2) is the time at end of par(2)
c.....ipar(3) is the time at end of par(3)
c.....ipar(5) is the time at end of par(3)+4 turns for extraction
c.....ipar(6) is the time at end of par(3)+ time to go back for 
c.....multiple extractions
c.....ipar(7) is the time at end of par(6)+ time to blow-up the beam 
c.....emittance for multiple extractions
c.....ipar(4) is the time to end the ramp of the nonlinearities
c.......................................................................
      Character*29 script
c.......................................................................
      Logical Dcheck
c.......................................................................
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
c.......................................................................
      turns=Dfloat(n-1)
      ripq=0d0
      do irip=1,Nriplq                        !..computes ripple on tune
        Arg=pi2*Dmod(Freqrq(irip)*turns+Phaserq(irip),1d0)
        ripq=ripq+Amprq(irip)*Dcos(Arg)
      enddo
c.......................................................................
      ripo=0d0
      do irip=1,Nriplo  !..computes ripple on other elements (sex.+oct.)
        Arg=pi2*Dmod(Freqro(irip)*turns+Phasero(irip),1d0)
        ripo=ripo+Ampro(irip)*Dcos(Arg)
      enddo
c.......................................................................
      if (n.le.ipar(4)) then                      !..nonlinearities ramp
        sexk=turns*par(8)*(1d0+ripo)
        octk=turns*par(9)*(1d0+ripo)
      else
        sexk=par(34)*(1d0+ripo)
        octk=par(4)*(1d0+ripo)
      endif
c.......................................................................
      par(22)=sexk
      par(21)=octk
c.......................................................................
      sexk=sexk*par(28)        !..rescales kick with the number of kicks
      octk=octk*par(28)
c.......................................................................
      if (n.le.ipar(4)) then                                !..tune ramp
        freq=par(1)
      elseif ((n.gt.ipar(4)).and.(n.le.ipar(1))) then
        freq=par(2)+par(5)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
      elseif ((n.gt.ipar(1)).and.(n.le.ipar(2))) then
        freq=par(2)
      elseif ((n.gt.ipar(2)).and.(n.le.ipar(3))) then
        freq=par(2)+par(6)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
      elseif ((n.gt.ipar(3)).and.(n.le.ipar(5))) then
        freq=par(3)
      elseif ((n.gt.ipar(5)).and.(n.le.ipar(6))) then
        freq=par(3)+par(10)*(n-ipar(5))
      elseif ((n.gt.ipar(6)).and.(n.le.ipar(7))) then
        freq=par(1)
      endif
c.......................................................................
      freq=freq*(1d0+ripq)                                !..adds ripple
      par(24)=freq
c.......................................................from Sept -> Oct
      x1= csepoctH*x(ipart,1)+ssepoctH*x(ipart,2)
      y1=-ssepoctH*x(ipart,1)+csepoctH*x(ipart,2)
c.......................................................................
      do jkick=1,ipar(8)
        y1=y1+octk*x1*x1*x1                                  !..Oct kick
      enddo
      octnl=y1
c.......................................................................
      x2= coctsexH*x1+soctsexH*octnl                  !..from Oct -> Sex
      y2=-soctsexH*x1+coctsexH*octnl
c.......................................................................
      do jkick=1,ipar(8)
        y2=y2+sexk*x2*x2                                     !..Sex kick
      enddo
      sexnl=y2
c...............................defines transfer matrix Sex/Sep: H plane
      PhAdv=AoffH+freq-FtuneH         !..freq does not contain int. part
      csexsepH=Dcos(PhAdv)
      ssexsepH=Dsin(PhAdv)
c.......................................................................
      t1= csexsepH*x2+ssexsepH*sexnl                 !..from Sex -> Sept
      t2=-ssexsepH*x2+csexsepH*sexnl
c.......................................................................
      Dcheck=(Dabs(t1).le.bbox(3)).and.(Dabs(t2).le.bbox(4))
      if (Dcheck) then                            !..checks for overflow
c.......................................................................
        x(ipart,1)=t1                  !..copies back final co-ordinates
        x(ipart,2)=t2
        Iuser_map2=0
c.......................................................................
      else
c.......................................................................
        x(ipart,1)=0d0
        x(ipart,2)=0d0
        Iuser_map2=1
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Integer Function Iuser_map2_s(ipart,n)
c.......................................................................
c.... Function to compute the evolution of particles. It is a 2D
c.... polynomial map including decapoles and condition for setting 
c.... detuning with amplitude to zero with octupoles
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
c.....par(1) is the initial Hor. frequency*pi2
c.....par(2) is the intermediate Hor. frequency*pi2
c.....par(3) is the final Hor. frequency*pi2
c.....par(4) is the octupole coefficient
c.....par(5) is the first Hor. slope
c.....par(6) is the second Hor. slope
c.....par(7) is the exponent of the first slope
c.....par(8) is the slope to ramp the sextupole strength
c.....par(9) is the slope to ramp the octupole strength
c.....par(10) is the slope to ramp up to intial tune for multiple 
c.....extractions
c.....par(11) is the beta ratio at sextupole location
c.....par(12) is the 3*beta ratio at the octupole location
c.....par(13) is the square of the beta ratio at the octupole location
c.....par(14) is the initial Ver. frequency*pi2
c.....par(15) is the intermediate Ver. frequency*pi2
c.....par(16) is the final Ver. frequency*pi2
c.....par(17) is the first Ver. slope
c.....par(18) is the second Ver. slope
c.....ipar(1) is the time between par(1) and par(2)
c.....ipar(2) is the time at end of par(2)
c.....ipar(3) is the time at end of par(3)
c.....ipar(5) is the time at end of par(3)+4 turns for extraction
c.....ipar(6) is the time at end of par(3)+ time to go back for 
c.....multiple extractions
c.....ipar(7) is the time at end of par(6)+ time to blow-up the beam 
c.....emittance for multiple extractions
c.....ipar(4) is the time to end the ramp of the nonlinearities
c.......................................................................
      Character*29 script
c.......................................................................
      Logical Dcheck
c.......................................................................
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
c.......................................................................
      turns=Dfloat(n-1)
      ripq=0d0
      do irip=1,Nriplq                        !..computes ripple on tune
        Arg=pi2*Dmod(Freqrq(irip)*turns+Phaserq(irip),1d0)
        ripq=ripq+Amprq(irip)*Dcos(Arg)
      enddo
c.......................................................................
      ripo=0d0
      do irip=1,Nriplo  !..computes ripple on other elements (sex.+oct.)
        Arg=pi2*Dmod(Freqro(irip)*turns+Phasero(irip),1d0)
        ripo=ripo+Ampro(irip)*Dcos(Arg)
      enddo
c.......................................................................
      if (n.le.ipar(4)) then                      !..nonlinearities ramp
        sexk=turns*par(8)*(1d0+ripo)
        octk=turns*par(9)*(1d0+ripo)
      else
        sexk=par(34)*(1d0+ripo)
        octk=par(4)*(1d0+ripo)
      endif
c.......................................................................
      par(22)=sexk
c.......................................................................
      if (n.le.ipar(4)) then                                !..tune ramp
        freq=par(1)
      elseif ((n.gt.ipar(4)).and.(n.le.ipar(1))) then
        freq=par(2)+par(5)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
      elseif ((n.gt.ipar(1)).and.(n.le.ipar(2))) then
        freq=par(2)
      elseif ((n.gt.ipar(2)).and.(n.le.ipar(3))) then
        freq=par(2)+par(6)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
      elseif ((n.gt.ipar(3)).and.(n.le.ipar(5))) then
        freq=par(3)
      elseif ((n.gt.ipar(5)).and.(n.le.ipar(6))) then
        freq=par(3)+par(10)*(n-ipar(5))
      elseif ((n.gt.ipar(6)).and.(n.le.ipar(7))) then
        freq=par(1)
      endif
c.......................................................................
      deck=par(4)                           !..defines decapole gradient
c......................defines octupole gradient to set detuning to zero
      octk=-1d0/6d0*(3d0*Dcos(0.5d0*freq)/Dsin(0.5d0*freq)+
     .                   Dcos(1.5d0*freq)/Dsin(1.5d0*freq)) 
      par(21)=octk
c.......................................................................
      freq=freq*(1d0+ripq)                                !..adds ripple
      par(24)=freq
c.......................................................................
      cc=Dcos(freq)
      ss=Dsin(freq)
c.......................................................................
      x1=x(ipart,2)
      y1=x(ipart,2)
c.......................................................................
      polnl=y1+x1*x1*(sexk+x1*(octk+deck*x1))          !..nonlinear kick
c.......................................................................
      t1= cc*x1+ss*polnl                                    !..iteration
      t2=-ss*x1+cc*polnl
c.......................................................................
      Dcheck=(Dabs(t1).le.bbox(3)).and.(Dabs(t2).le.bbox(4))
      if (Dcheck) then                            !..checks for overflow
c.......................................................................
        x(ipart,1)=t1                  !..copies back final co-ordinates
        x(ipart,2)=t2
        Iuser_map2_s=0
c.......................................................................
      else
c.......................................................................
        x(ipart,1)=0d0
        x(ipart,2)=0d0
        Iuser_map2_s=1
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Integer Function Iuser_map4(ipart,n)
c.......................................................................
c.... Function to compute the evolution of particles. It is a 4D
c.... polynomial map
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
c.....par(1) is the initial frequency*pi2
c.....par(2) is the intermediate frequency*pi2
c.....par(3) is the final frequency*pi2
c.....par(4) is the octupole coefficient
c.....par(5) is the first slope
c.....par(6) is the second slope
c.....par(7) is the exponent of the first slope
c.....par(8) is the slope to ramp the sextupole strength
c.....par(9) is the slope to ramp the octupole strength
c.....par(10) is the slope to ramp up to intial tune for multiple 
c.....extractions
c.....par(11) is the beta ratio at sextupole location
c.....par(12) is the 3*beta ratio at the octupole location
c.....par(13) is the square of the beta ratio at the octupole location
c.....par(14) is the initial Ver. frequency*pi2
c.....par(15) is the intermediate Ver. frequency*pi2
c.....par(16) is the final Ver. frequency*pi2
c.....par(17) is the first Ver. slope
c.....par(18) is the second Ver. slope
c.....par(19) is the slope to ramp up to intial Ver. tune for multiple 
c.....extractions
c.....ipar(1) is the time between par(1) and par(2)
c.....ipar(2) is the time at end of par(2)
c.....ipar(3) is the time at end of par(3)
c.....ipar(5) is the time at end of par(3)+4 turns for extraction
c.....ipar(6) is the time at end of par(3)+ time to go back for 
c.....multiple extractions
c.....ipar(7) is the time at end of par(6)+ time to blow-up the beam 
c.....emittance for multiple extractions
c.....ipar(4) is the time to end the ramp of the nonlinearities
c.......................................................................
      Character*29 script
c.......................................................................
      Logical Dcheck
c.......................................................................
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
c.......................................................................
      turns=Dfloat(n-1)
      ripq=0d0
      do irip=1,Nriplq                                !..computes ripple
        Arg=pi2*Dmod(Freqrq(irip)*turns+Phaserq(irip),1d0)
        ripq=ripq+Amprq(irip)*Dcos(Arg)
      enddo
c.......................................................................
      ripo=0d0
      do irip=1,Nriplo  !..computes ripple on other elements (sex.+oct.)
        Arg=pi2*Dmod(Freqro(irip)*turns+Phasero(irip),1d0)
        ripo=ripo+Ampro(irip)*Dcos(Arg)
      enddo
c.......................................................................
      if (n.le.ipar(4)) then                      !..nonlinearities ramp
        sexk=turns*par(8)*(1d0+ripo)
        octk=turns*par(9)*(1d0+ripo)
      else
        sexk=par(34)*(1d0+ripo)
        octk=par(4)*(1d0+ripo)
      endif
c.......................................................................
      par(22)=sexk
      par(21)=octk
c.......................................................................
      sexk=sexk*par(28)        !..rescales kick with the number of kicks
      octk=octk*par(28)
c.......................................................................
      if (n.le.ipar(4)) then                  !..hor. and ver. tune ramp
        freqx=par(1)
        freqy=par(14)
      elseif ((n.gt.ipar(4)).and.(n.le.ipar(1))) then
        freqx=par(2)+par(5)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
        freqy=par(15)+par(17)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
      elseif ((n.gt.ipar(1)).and.(n.le.ipar(2))) then
        freqx=par(2)
        freqy=par(15)
      elseif ((n.gt.ipar(2)).and.(n.le.ipar(3))) then
        freqx=par(2)+par(6)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
        freqy=par(15)+par(18)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
      elseif ((n.gt.ipar(3)).and.(n.le.ipar(5))) then
        freqx=par(3)
        freqy=par(16)
      elseif ((n.gt.ipar(5)).and.(n.le.ipar(6))) then
        freqx=par(3)+par(10)*(n-ipar(5))
        freqy=par(16)+par(19)*(n-ipar(5))
      elseif ((n.gt.ipar(6)).and.(n.le.ipar(7))) then
        freqx=par(1)
        freqy=par(14)
      endif
c.......................................................................
      freqx=freqx*(1d0+ripq)                              !..adds ripple
      par(24)=freqx
c.......................................................................
      freqy=freqy*(1d0+ripq)                              !..adds ripple
      par(23)=freqy
c.......................................................from Sept -> Oct
      x1 = csepoctH*x(ipart,1)+ssepoctH*x(ipart,2)
      px1=-ssepoctH*x(ipart,1)+csepoctH*x(ipart,2)
      y1 = csepoctV*x(ipart,3)+ssepoctV*x(ipart,4)
      py1=-ssepoctV*x(ipart,3)+csepoctV*x(ipart,4)
c.......................................................................
      do jkick=1,ipar(8)
        px1=px1+octk*(        x1*x1*x1-par(12)*x1*y1*y1)     !..Oct kick
        py1=py1+octk*(par(12)*x1*x1*y1-par(13)*y1*y1*y1)
      enddo
      octnlx=px1
      octnly=py1
c.......................................................................
      x2 = coctsexH*x1+soctsexH*octnlx                !..from Oct -> Sex
      px2=-soctsexH*x1+coctsexH*octnlx
      y2 = coctsexV*y1+soctsexV*octnly
      py2=-soctsexV*y1+coctsexV*octnly
c.......................................................................
      do jkick=1,ipar(8)
        px2=px2+sexk*(x2*x2-par(11)*y2*y2)                   !..Sex kick
        py2=py2-2d0*par(11)*sexk*x2*y2
      enddo
      sexnlx=px2
      sexnly=py2
c...............................defines transfer matrix Sex/Sep: H plane
      PhAdvH=AoffH+freqx-FtuneH       !..freq does not contain int. part
      csexsepH=Dcos(PhAdvH)
      ssexsepH=Dsin(PhAdvH)
c.......................................................................
      PhAdvV=AoffV+freqy-FtuneV       !..freq does not contain int. part
      csexsepV=Dcos(PhAdvV)
      ssexsepV=Dsin(PhAdvV)
c.......................................................................
      t1= csexsepH*x2+ssexsepH*sexnlx                !..from Sex -> Sept
      t2=-ssexsepH*x2+csexsepH*sexnlx
      t3= csexsepV*y2+ssexsepV*sexnly
      t4=-ssexsepV*y2+csexsepV*sexnly
c.......................................................................
      Dcheck=(Dabs(t1).le.bbox(3)).and.(Dabs(t2).le.bbox(4)).and.
     .       (Dabs(t3).le.bbox(3)).and.(Dabs(t4).le.bbox(4))
      if (Dcheck) then                            !..checks for overflow
c.......................................................................
        x(ipart,1)=t1                  !..copies back final co-ordinates
        x(ipart,2)=t2
        x(ipart,3)=t3
        x(ipart,4)=t4
        Iuser_map4=0
c.......................................................................
      else
c.......................................................................
        x(ipart,1)=0d0
        x(ipart,2)=0d0
        x(ipart,3)=0d0
        x(ipart,4)=0d0
        Iuser_map4=1
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Integer Function Iuser_map4_s(ipart,n)
c.......................................................................
c.... Function to compute the evolution of particles. It is a 4D
c.... polynomial map
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
c.....par(1) is the initial frequency*pi2
c.....par(2) is the intermediate frequency*pi2
c.....par(3) is the final frequency*pi2
c.....par(4) is the octupole coefficient
c.....par(5) is the first slope
c.....par(6) is the second slope
c.....par(7) is the exponent of the first slope
c.....par(8) is the slope to ramp the sextupole strength
c.....par(9) is the slope to ramp the octupole strength
c.....par(10) is the slope to ramp up to intial tune for multiple 
c.....extractions
c.....par(11) is the beta ratio at sextupole location
c.....par(12) is the 3*beta ratio at the octupole location
c.....par(13) is the square of the beta ratio at the octupole location
c.....par(14) is the initial Ver. frequency*pi2
c.....par(15) is the intermediate Ver. frequency*pi2
c.....par(16) is the final Ver. frequency*pi2
c.....par(17) is the first Ver. slope
c.....par(18) is the second Ver. slope
c.....par(19) is the slope to ramp up to intial Ver. tune for multiple 
c.....extractions
c.....ipar(1) is the time between par(1) and par(2)
c.....ipar(2) is the time at end of par(2)
c.....ipar(3) is the time at end of par(3)
c.....ipar(5) is the time at end of par(3)+4 turns for extraction
c.....ipar(6) is the time at end of par(3)+ time to go back for 
c.....multiple extractions
c.....ipar(7) is the time at end of par(6)+ time to blow-up the beam 
c.....emittance for multiple extractions
c.....ipar(4) is the time to end the ramp of the nonlinearities
c.......................................................................
      Character*29 script
c.......................................................................
      Logical Dcheck
c.......................................................................
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
c.......................................................................
      turns=Dfloat(n-1)
      ripq=0d0
      do irip=1,Nriplq                                !..computes ripple
        Arg=pi2*Dmod(Freqrq(irip)*turns+Phaserq(irip),1d0)
        ripq=ripq+Amprq(irip)*Dcos(Arg)
      enddo
c.......................................................................
      ripo=0d0
      do irip=1,Nriplo  !..computes ripple on other elements (sex.+oct.)
        Arg=pi2*Dmod(Freqro(irip)*turns+Phasero(irip),1d0)
        ripo=ripo+Ampro(irip)*Dcos(Arg)
      enddo
c.......................................................................
      if (n.le.ipar(4)) then                      !..nonlinearities ramp
        sexk=turns*par(8)*(1d0+ripo)
        octk=turns*par(9)*(1d0+ripo)
      else
        sexk=par(34)*(1d0+ripo)
        octk=par(4)*(1d0+ripo)
      endif
c.......................................................................
      par(22)=sexk
      par(21)=octk
c.......................................................................
      if (n.le.ipar(4)) then                  !..hor. and ver. tune ramp
        freqx=par(1)
        freqy=par(14)
      elseif ((n.gt.ipar(4)).and.(n.le.ipar(1))) then
        freqx=par(2)+par(5)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
        freqy=par(15)+par(17)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-n)))
      elseif ((n.gt.ipar(1)).and.(n.le.ipar(2))) then
        freqx=par(2)
        freqy=par(15)
      elseif ((n.gt.ipar(2)).and.(n.le.ipar(3))) then
        freqx=par(2)+par(6)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
        freqy=par(15)+par(18)*Dexp(par(31)*Dlog(Dfloat(n-ipar(2))))
      elseif ((n.gt.ipar(3)).and.(n.le.ipar(5))) then
        freqx=par(3)
        freqy=par(16)
      elseif ((n.gt.ipar(5)).and.(n.le.ipar(6))) then
        freqx=par(3)+par(10)*(n-ipar(5))
        freqy=par(16)+par(19)*(n-ipar(5))
      elseif ((n.gt.ipar(6)).and.(n.le.ipar(7))) then
        freqx=par(1)
        freqy=par(14)
      endif
c.......................................................................
      freqx=freqx*(1d0+ripq)                              !..adds ripple
      par(24)=freqx
c.......................................................................
      freqy=freqy*(1d0+ripq)                              !..adds ripple
      par(23)=freqy
c.......................................................from Sept -> Oct
      x1 = csepoctH*x(ipart,1)+ssepoctH*x(ipart,2)
      px1=-ssepoctH*x(ipart,1)+csepoctH*x(ipart,2)
      y1 = csepoctV*x(ipart,3)+ssepoctV*x(ipart,4)
      py1=-ssepoctV*x(ipart,3)+csepoctV*x(ipart,4)
c.......................................................................
      octnlx=px1+octk*(        x1*x1*x1-par(12)*x1*y1*y1)    !..Oct kick
      octnly=py1+octk*(par(12)*x1*x1*y1-par(13)*y1*y1*y1)
c.......................................................................
      x2 = coctsexH*x1+soctsexH*octnlx                !..from Oct -> Sex
      px2=-soctsexH*x1+coctsexH*octnlx
      y2 = coctsexV*y1+soctsexV*octnly
      py2=-soctsexV*y1+coctsexV*octnly
c.......................................................................
      sexnlx=px2+sexk*(x2*x2-par(11)*y2*y2)                  !..Sex kick
      sexnly=py2-2d0*par(11)*sexk*x2*y2
c...............................defines transfer matrix Sex/Sep: H plane
      PhAdvH=AoffH+freqx-FtuneH       !..freq does not contain int. part
      csexsepH=Dcos(PhAdvH)
      ssexsepH=Dsin(PhAdvH)
c.......................................................................
      PhAdvV=AoffV+freqy-FtuneV       !..freq does not contain int. part
      csexsepV=Dcos(PhAdvV)
      ssexsepV=Dsin(PhAdvV)
c.......................................................................
      t1= csexsepH*x2+ssexsepH*sexnlx                !..from Sex -> Sept
      t2=-ssexsepH*x2+csexsepH*sexnlx
      t3= csexsepV*y2+ssexsepV*sexnly
      t4=-ssexsepV*y2+csexsepV*sexnly
c.......................................................................
      Dcheck=(Dabs(t1).le.bbox(3)).and.(Dabs(t2).le.bbox(4)).and.
     .       (Dabs(t3).le.bbox(3)).and.(Dabs(t4).le.bbox(4))
      if (Dcheck) then                            !..checks for overflow
c.......................................................................
        x(ipart,1)=t1                  !..copies back final co-ordinates
        x(ipart,2)=t2
        x(ipart,3)=t3
        x(ipart,4)=t4
        Iuser_map4_s=0
c.......................................................................
      else
c.......................................................................
        x(ipart,1)=0d0
        x(ipart,2)=0d0
        x(ipart,3)=0d0
        x(ipart,4)=0d0
        Iuser_map4_s=1
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Read_par
c.......................................................................
c.... Subroutine to read input parameters. It allows to put comments in 
c.... the input file. Comments should start with a '!' sign.
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
      Real rxvec(Max_part),rxpvec(Max_part)
      Real ryvec(Max_part),rypvec(Max_part)
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
      Character*80 aline
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Distribution/Radius,Sigma,Sigmac,Amean,Bmean,Radv,
     .                    Sigv,Sigvc,Ameav,Bmeav,Ameanc,Bmeanc,
     .                    Ameanvc,Bmeanvc,xsep,Sigman,Sigmah(Max_bord),
     .                    Akick(Max_bord),rxvec,rxpvec,
     .                    ryvec,rypvec,Namp,Nangles,Nampc,Nanglesc
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      Data AmatSepOct,AmatOctSex,Tmat,Tmati/48*0d0/
c.......................................................................
      do idiag=1,2
        do ii=1,2
          Tmat(idiag,idiag,ii)=1d0                           !..Identity

        enddo
      enddo
c.......................................................................
      read(1,'(a80)',end=10) aline
      ilimit=index(aline,'!')-1
      read(aline(1:ilimit),*) Ntur,Nout,Noff,Idim,Iext,Iorl,Ioru,Mult,
     .                        Idump
c.......................................................................
      read(1,'(a80)',end=10) aline
      ilimit=index(aline,'!')-1
      read(aline(1:ilimit),*) iproj
      read(aline(1:ilimit),*) idummy,(icpr(ii),ii=1,2*iproj)
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Om1,Om2,Om3,Oct,Nkick,Omv1,Omv2,Omv3,
     .                                                     Itwiss,Bratio
c.......................................................................
      do iplace=1,4 !..Reads Twiss parameters (used only if Itwiss.eq.1)
        read(1,'(a80)',end=10) aline 
        ilimit=index(aline,'!')-1 
        read(aline(1:ilimit),*) BetaH(iplace),AlphaH(iplace),  !oct. 1st
     .                          AmuH(iplace),                  !sex. 2nd
     .                          BetaV(iplace),AlphaV(iplace),  !sep. 3rd
     .                          AmuV(iplace)                   !end  4th
      enddo
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Ico,emx,emy
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Nriplq,Nriplo
      do irip=1,Nriplq                                 !..ripple on tune
        read(1,'(a80)',end=10) aline 
        ilimit=index(aline,'!')-1 
        read(aline(1:ilimit),*) Amprq(irip),Fr,Ph  !..Fr,Ph units of 2pi
        Freqrq(irip)=2.2*Fr*1d-6                 !..converts Hz -> turns
        Phaserq(irip)=Ph
      enddo
      do irip=1,Nriplo           !..ripple on other elements (sex.+oct.)
        read(1,'(a80)',end=10) aline 
        ilimit=index(aline,'!')-1 
        read(aline(1:ilimit),*) Ampro(irip),Fr,Ph  !..Fr,Ph units of 2pi
        Freqro(irip)=2.2*Fr*1d-6                 !..converts Hz -> turns
        Phasero(irip)=Ph
      enddo
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Nt1,Nt2,Nt3,Nt4,Nt5,Nt6
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) expon,expon1,expxct,expoct
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Icomp,Idist,Ninj
c.......................................................................
      read(1,'(a80)',end=10) aline        !..computation of distribution
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Nbin,Nbx,Nby 
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      read(aline(1:ilimit),*) Nboxes(1),(Ibox(jbox),jbox=1,Nboxes(1))
c.......................................................................
      read(1,'(a80)',end=10) aline 
      ilimit=index(aline,'!')-1 
      if (Idist.eq.1) then                       !..Uniform distribution
        read(aline(1:ilimit),*) Radius,Amean,Bmean,Radv,Ameav,Bmeav,
     .                          Namp,Nangles,
     .                          Icentre,Sigmac,Ameanc,Bmeanc,Sigvc,
     .                          Ameavc,Bmeavc,Nampc,Nanglesc
        SSigma=Radius/Dsqrt(6d0)
      elseif (Idist.ge.2) then           !..Gaussian/hollow distribution
        read(aline(1:ilimit),*) Sigma,Amean,Bmean,Sigv,Ameav,Bmeav,Namp,
     .                          Nangles,
     .                          Icentre,Sigmac,Ameanc,Bmeanc,Sigvc,
     .                          Ameavc,Bmeavc,Nampc,Nanglesc
        SSigma=Sigma
      endif
c.......................................................................
      read(1,'(a80)',end=10) aline
      ilimit=index(aline,'!')-1
      read(aline(1:ilimit),*) bbox(1),bbox(2),bbox(3),bbox(4),xsep,Ianim
      read(1,'(a80)',end=10) aline
      ilimit=index(aline,' ')-1
      fnamei=aline(1:ilimit)                   !..file name output files
c.......................................................................
 10   close(1)
c.......................................................................
      Nboxes(2)=1
c.......................................................................
      ItuneH=AmuH(4)             !..defines integer part of machine tune
      ItuneV=AmuV(4)
      FtuneH=pi2*(AmuH(4)-ItuneH)  !..defines frac. part of machine tune
      FtuneV=pi2*(AmuV(4)-ItuneV)
      if (Itwiss.eq.0) then        !..defines the appropriate parameters
c.......................................................................
        par(11)=Bratio                !..defines beta ratio for 4D kicks
        par(12)=3*par(11)
        par(13)=par(11)*par(11)
c.......................................................................
      elseif(Itwiss.eq.1) then
c.......................................................................
        par(11)=BetaV(2)/BetaH(2)     !..defines beta ratio for 4D kicks
        par(12)=3*BetaV(1)/BetaH(1)
        par(13)=(BetaV(1)/BetaH(1))*(BetaV(1)/BetaH(1))
c...............................defines transfer matrix Sep/Oct: H plane
        AmatSepOct(1,1)=Dsqrt(BetaH(1)/BetaH(3))*
     .  (Dcos(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))+AlphaH(3)*
     .   Dsin(pi2*(AmuH(1)+AmuH(4)-AmuH(3))))
        AmatSepOct(1,2)=Dsqrt(BetaH(1)*BetaH(3))*
     .   Dsin(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))
        AmatSepOct(2,1)=1d0/Dsqrt(BetaH(1)*BetaH(3))*
     .  ((AlphaH(1)-AlphaH(3))*Dcos(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))+
     .  (1d0+AlphaH(3)*AlphaH(1))*Dsin(pi2*(AmuH(1)+AmuH(4)-AmuH(3))))
        AmatSepOct(2,2)=Dsqrt(BetaH(3)/BetaH(1))*
     .  (Dcos(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))-AlphaH(1)*
     .   Dsin(pi2*(AmuH(1)+AmuH(4)-AmuH(3))))
c...............................defines transfer matrix Sep/Oct: V plane
        AmatSepOct(3,3)=Dsqrt(BetaV(1)/BetaV(3))*
     .  (Dcos(pi2*(AmuV(1)-AmuV(3)))+AlphaV(3)*
     .   Dsin(pi2*(AmuV(1)-AmuV(3))))
        AmatSepOct(3,4)=Dsqrt(BetaV(1)*BetaV(3))*
     .   Dsin(pi2*(AmuV(1)-AmuV(3)))
        AmatSepOct(4,3)=1d0/Dsqrt(BetaV(1)*BetaV(3))*
     .  ((AlphaV(1)-AlphaV(3))*Dcos(pi2*(AmuV(1)-AmuV(3)))+
     .  (1d0+AlphaV(3)*AlphaV(1))*Dsin(pi2*(AmuV(1)-AmuV(3))))
        AmatSepOct(4,4)=Dsqrt(BetaV(3)/BetaV(1))*
     .  (Dcos(pi2*(AmuV(1)-AmuV(3)))-AlphaV(1)*
     .   Dsin(pi2*(AmuV(1)-AmuV(3))))
c...............................defines transfer matrix Oct/Sex: H plane
        AmatOctSex(1,1)=Dsqrt(BetaH(2)/BetaH(1))*
     .  (Dcos(pi2*(AmuH(2)-AmuH(1)))+AlphaH(1)*
     .   Dsin(pi2*(AmuH(2)-AmuH(1))))
        AmatOctSex(1,2)=Dsqrt(BetaH(2)*BetaH(1))*
     .   Dsin(pi2*(AmuH(2)-AmuH(1)))
        AmatOctSex(2,1)=1d0/Dsqrt(BetaH(2)*BetaH(1))*
     .  ((AlphaH(2)-AlphaH(1))*Dcos(pi2*(AmuH(2)-AmuH(1)))+
     .  (1d0+AlphaH(1)*AlphaH(2))*Dsin(pi2*(AmuH(2)-AmuH(1))))
        AmatOctSex(2,2)=Dsqrt(BetaH(1)/BetaH(2))*
     .  (Dcos(pi2*(AmuH(2)-AmuH(1)))-AlphaH(2)*
     .   Dsin(pi2*(AmuH(2)-AmuH(1))))
c...............................defines transfer matrix Oct/Sex: V plane
        AmatOctSex(3,3)=Dsqrt(BetaV(2)/BetaV(1))*
     .  (Dcos(pi2*(AmuV(2)-AmuV(1)))+AlphaV(1)*
     .   Dsin(pi2*(AmuV(2)-AmuV(1))))
        AmatOctSex(3,4)=Dsqrt(BetaV(2)*BetaV(1))*
     .   Dsin(pi2*(AmuV(2)-AmuV(1)))
        AmatOctSex(4,3)=1d0/Dsqrt(BetaV(2)*BetaV(1))*
     .  ((AlphaV(2)-AlphaV(1))*Dcos(pi2*(AmuV(2)-AmuV(1)))+
     .  (1d0+AlphaV(1)*AlphaV(2))*Dsin(pi2*(AmuV(2)-AmuV(1))))
        AmatOctSex(4,4)=Dsqrt(BetaV(1)/BetaV(2))*
     .  (Dcos(pi2*(AmuV(2)-AmuV(1)))-AlphaV(2)*
     .   Dsin(pi2*(AmuV(2)-AmuV(1))))
c.......................................................................
      endif
c.......................................................................
      Alambda=1d0 !..defines normalisation form special C-S co-ordinates
      Alambdai=1d0
      if (Ico.ne.0) then
        Alambda=SSigma/Dsqrt(emx)
        if (Alambda.ne.0d0) Alambdai=1d0/Alambda
      endif
c.......................................................................
      if (Ico.eq.2) then
c.......................................................................
        Tmat(1,1,1)=Dsqrt(BetaH(3))  !..defines transformation C-S/phys.
        Tmat(2,1,1)=-AlphaH(3)/Dsqrt(BetaH(3))
        Tmat(2,2,1)=1d0/Dsqrt(BetaH(3))
        Tmat(1,1,2)=Dsqrt(BetaV(3))
        Tmat(2,1,2)=-AlphaV(3)/Dsqrt(BetaV(3))
        Tmat(2,2,2)=1d0/Dsqrt(BetaV(3))
        Tmati(1,1,1)=1d0/Dsqrt(BetaH(3))  !..defines transfor. phys./C-S
        Tmati(2,1,1)=AlphaH(3)/Dsqrt(BetaH(3))
        Tmati(2,2,1)=Dsqrt(BetaH(3))
        Tmati(1,1,2)=1d0/Dsqrt(BetaV(3))
        Tmati(2,1,2)=AlphaV(3)/Dsqrt(BetaV(3))
        Tmati(2,2,2)=Dsqrt(BetaV(3))
c.......................................................................
      endif
c.......................................................................
      if (Ico.gt.0) then        !..transforms from special to normal C-S
         do irow=1,2
           do icol=1,2
             do ii=1,2
               Tmat(irow,icol,ii)=Alambdai*Tmat(irow,icol,ii)
               Tmati(irow,icol,ii)=Alambda*Tmati(irow,icol,ii)
             enddo
           enddo
         enddo
      endif
c.......................................................................
      Dett=Tmat(1,1,1)*Tmat(2,2,1)-Tmat(1,2,1)*Tmat(2,1,1)
      Detti= Tmati(1,1,1)*Tmati(2,2,1)-Tmati(1,2,1)*Tmati(2,1,1)
c...........................................transforms 2nd order moments
      do ii=1,2
        Tmom(1,1,ii)=Tmat(1,1,ii)*Tmat(1,1,ii)
        Tmom(1,2,ii)=2d0*Tmat(1,1,ii)*Tmat(1,2,ii)
        Tmom(1,3,ii)=Tmat(1,2,ii)*Tmat(1,2,ii)
        Tmom(2,1,ii)=Tmat(1,1,ii)*Tmat(2,1,ii)
        Tmom(2,2,ii)=Tmat(1,1,ii)*Tmat(2,2,ii)+
     .                  Tmat(1,2,ii)*Tmat(2,1,ii)
        Tmom(2,3,ii)=Tmat(1,2,ii)*Tmat(2,2,ii)
        Tmom(3,1,ii)=Tmat(2,1,ii)*Tmat(2,1,ii)
        Tmom(3,2,ii)=2d0*Tmat(2,1,ii)*Tmat(2,2,ii)
        Tmom(3,3,ii)=Tmat(2,2,ii)*Tmat(2,2,ii)
c..............................................inverse transform moments
        Tmomi(1,1,ii)=Tmati(1,1,ii)*Tmati(1,1,ii)
        Tmomi(1,2,ii)=2d0*Tmati(1,1,ii)*Tmati(1,2,ii)
        Tmomi(1,3,ii)=Tmati(1,2,ii)*Tmati(1,2,ii)
        Tmomi(2,1,ii)=Tmati(1,1,ii)*Tmati(2,1,ii)
        Tmomi(2,2,ii)=Tmati(1,1,ii)*Tmati(2,2,ii)+
     .                   Tmati(1,2,ii)*Tmati(2,1,ii)
        Tmomi(2,3,ii)=Tmati(1,2,ii)*Tmati(2,2,ii)
        Tmomi(3,1,ii)=Tmati(2,1,ii)*Tmati(2,1,ii)
        Tmomi(3,2,ii)=2d0*Tmati(2,1,ii)*Tmati(2,2,ii)
        Tmomi(3,3,ii)=Tmati(2,2,ii)*Tmati(2,2,ii)
c.......................................................................
      enddo
c.......................................................................
      csepoctH=Dcos(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))   !..rotation sep/oct
      ssepoctH=Dsin(pi2*(AmuH(1)+AmuH(4)-AmuH(3)))
      csepoctV=Dcos(pi2*(AmuV(1)+AmuV(4)-AmuV(3)))
      ssepoctV=Dsin(pi2*(AmuV(1)+AmuV(4)-AmuV(3)))
c.......................................................................
      coctsexH=Dcos(pi2*(AmuH(2)-AmuH(1)))           !..rotation oct/sex
      soctsexH=Dsin(pi2*(AmuH(2)-AmuH(1)))
      coctsexV=Dcos(pi2*(AmuV(2)-AmuV(1)))
      soctsexV=Dsin(pi2*(AmuV(2)-AmuV(1)))
c.......................................................................
      AoffH=pi2*(AmuH(3)-AmuH(2))              !..rotation angle sex/sep
      AoffV=pi2*(AmuV(3)-AmuV(2))
c.......................................................................
      Ntot=Namp*Nangles
      Npmult=Ntot
c.......................................................................
      if (Ntot.gt.Max_part) then
        write(*,*) 
     .    '*** Error: max. number of particles exceeded. Program stops.'
        stop
      endif
c.......................................................................
      if ((Icomp.eq.7).and.(Ninj*Ntot.gt.Max_part)) then
        write(*,*) 
     .    '*** Error: max. number of particles exceeded. Program stops.'
        write(*,*) 
     .    '*** Ninj*Npart has to be smaller than ',Max_part
        stop
      endif
c.......................................................................
      if (Ntur.ne.(Nt1+Nt2+Nt3)) then
        write(*,*) '*** Error: inconsistency in the definition of tune',
     .             ' curve. Program stops.'
        stop
      endif
c.......................................................................
      Nturn=Ntur+Nt4
c.......................................................................
      ipar(1)=Nt1+Nt4                            !..fills in ipar vector
      ipar(2)=ipar(1)+Nt2
      ipar(3)=ipar(2)+Nt3
      ipar(4)=Nt4
      ipar(5)=ipar(3)+4*Iext
      ipar(6)=ipar(5)+Nt5
      ipar(7)=ipar(6)+Nt6
      ipar(8)=Nkick
      ipar(9)=Nturn/1000
      iramp=ipar(6)-ipar(5)
      iblowup=ipar(7)-ipar(6)
      Iper=Nturn+iramp+iblowup
c.......................................................................
      par(1)=pi2*Om1                             ! ..fills in par vector
      par(2)=pi2*Om2
      par(3)=pi2*Om3
      par(4)=Oct
      par(34)=1d0
      Slopeh1=(par(1)-par(2))/Dexp(expon*Dlog(Dfloat(ipar(1)-Nt4)))
      Slopeh2=(par(3)-par(2))/Dexp(expon1*Dlog(Dfloat(ipar(3)-ipar(2))))
      par(5)=Slopeh1
      par(6)=Slopeh2
      par(7)=expon
      par(31)=expon1
      par(8)=1d0/Dfloat(Nt4-1)
      par(9)=Oct/Dfloat(Nt4-1)
      par(10)=0d0
      par(19)=0d0
      par(14)=pi2*Omv1
      par(15)=pi2*Omv2
      par(16)=pi2*Omv3
      Slopev1=(par(14)-par(15))/Dexp(expon*Dlog(Dfloat(ipar(1)-Nt4)))
      Slopev2=(par(16)-par(15))/
     .                     Dexp(expon1*Dlog(Dfloat(ipar(3)-ipar(2))))
      par(17)=Slopev1
      par(18)=Slopev2
      if (Nt5.ne.0) then
        par(10)=(par(1)-par(3))/Dfloat(Nt5)
        par(19)=(par(14)-par(16))/Dfloat(Nt5)
      endif
      if (Icomp.eq.8) then
        if (expoct.eq.0d0) then
          Slopeoct=0d0
        else
          Slopeoct=par(4)/(1-Dexp(expoct*Dlog(Dfloat(ipar(1)-Nt4))))
        endif
        par(25)=Slopeoct
        par(26)=expoct
        par(27)=par(4)
        if (expxct.eq.0d0) then
          Slopexct=0d0
        else
          Slopexct=1d0/(1-Dexp(expxct*Dlog(Dfloat(ipar(1)-Nt4))))
        endif
        par(35)=Slopexct
        par(36)=expxct
        par(37)=par(34)
      endif
      par(28)=1d0/Nkick
c.......................................................................
      Call Gen_boxes                                 !..limits for boxes
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Find_script
c.......................................................................
c.... Subroutine to read the directory where the scripts to perform file
c.... conversion are installed.
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Character*29 script
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
c.......................................................................
      Call system('which pstops | cat > tmp.txt')
      open(1,file='tmp.txt',status='old')
      read(1,'(a29)') script(0)
      close(1)
      Call system('rm -f tmp.txt')
      iscriptl(0)=index(script(0),'/')
      iscriptu(0)=index(script(0),'pstops')+6
c.......................................................................
      script(1)='$HOME/private/script/pstogif '
      iscriptl(1)=1
      iscriptu(1)=29
c.......................................................................
      Call system('which ps2pdf | cat > tmp.txt')
      open(1,file='tmp.txt',status='old')
      read(1,'(a29)') script(2)
      close(1)
      Call system('rm -f tmp.txt')
      iscriptl(2)=index(script(2),'/')
      iscriptu(2)=index(script(2),'ps2pdf')+6
c.......................................................................
      script(3)='$HOME/private/script/pstogif '
      iscriptl(3)=1
      iscriptu(3)=29
c.......................................................................
      Call system('which ps2epsi | cat > tmp.txt')
      open(1,file='tmp.txt',status='old')
      read(1,'(a29)') script(4)
      close(1)
      Call system('rm -f tmp.txt')
      iscriptl(4)=index(script(4),'/')
      iscriptu(4)=index(script(4),'ps2epsi')+7
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Gen_dist
c.......................................................................
c.... Subroutine to generate the initial distribution of particles.
c.... Idist = 1: Uniform distribution
c.... Idist = 2: Gaussian distribution
c.... Idist = 3: Hollow distribution in both H and V phase space
c.... Idist = 4: Hollow distribution in H and Gaussian in V phase space
c....
c.... Author: M. Giovannozzi
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Dimension avg(2),smom(3)
      Dimension Ipoin(Max_part)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
      Real rxvec(Max_part),rxpvec(Max_part),radxv(Max_part)
      Real ryvec(Max_part),rypvec(Max_part),radyv(Max_part)
c.......................................................................
      Character*3 ext
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Momen/xmean(Max_boxes,2),xpmean(Max_boxes,2),
     .             x2(Max_boxes,2),xp2(Max_boxes,2),xxp(Max_boxes,2),
     .             x3(Max_boxes,2),xp3(Max_boxes,2),x2xp(Max_boxes,2),
     .             xxp2(Max_boxes,2),x4(Max_boxes,2),xp4(Max_boxes,2),
     .             x2xp2(Max_boxes,2),xxp3(Max_boxes,2),
     .             x3xp(Max_boxes,2),
     .             emit(Max_boxes,2),beta(Max_boxes,2),
     .             alpha(Max_boxes,2),ellipse(Max_boxes,Max_dim),
     .             surf(Max_boxes,Max_dim),halop1(Max_boxes,2),
     .             halop2(Max_boxes,2),Idist1D(Max_bin,Max_boxes,2)
      Common/Distribution/Radius,Sigma,Sigmac,Amean,Bmean,Radv,
     .                    Sigv,Sigvc,Ameav,Bmeav,Ameanc,Bmeanc,
     .                    Ameanvc,Bmeanvc,xsep,Sigman,Sigmah(Max_bord),
     .                    Akick(Max_bord),rxvec,rxpvec,
     .                    ryvec,rypvec,Namp,Nangles,Nampc,Nanglesc
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
C.......................................................................
      external ranlux
C.......................................................................
      Data Idist1D/Max_e1D*0/
      Data Id1Dmax/Max_dim*-1/,Id2Dmax/Max_pro2*-1/
c.......................................................................
      isw=1     !..used to switch on surface computation via convex hull
      isw1=0     !..used to switch on inj. mode in particles' generation
      if (Ninj.ne.0) isw1=1
c.......................................................................
      do jdim=1,Idim/2                                 !..initialisation
        xmean(1,jdim)=0d0
        xpmean(1,jdim)=0d0
        x2(1,jdim)=0d0
        xxp(1,jdim)=0d0
        xp2(1,jdim)=0d0
        x3(1,jdim)=0d0      !..higher-order moments used to compute halo
        x2xp(1,jdim)=0d0
        xxp2(1,jdim)=0d0 
        xp3(1,jdim)=0d0
        x4(1,jdim)=0d0
        x2xp2(1,jdim)=0d0
        x3xp(1,jdim)=0d0
        xxp3(1,jdim)=0d0
        xp4(1,jdim)=0d0
      enddo
c.......................................................................
      cminx=1d38                                       !..initialisation
      cminxp=1d38
      cmaxx=-1d38
      cmaxxp=-1d38
c.......................................................................
      if (Idim.eq.4) then                              !..initialisation
        cminy=1d38
        cminyp=1d38
        cmaxy=-1d38
        cmaxyp=-1d38
      endif
c.......................................................................
      if (Idist.eq.1) then                       !..uniform distribution
c.......................................................................
        lux=2                                         !..2D distribution
        iseed=7675039
        NRNGen=Ntot*(1+isw1*(Ninj-1))+Icentre*Namp*Nangles
        Call rluxgo(lux,iseed,0,0)        !..initializes random sequence
        Call ranlux(rxvec,NRNGen)             !..random number generator
        Call ranlux(rxpvec,NRNGen)
c.......................................................................
        do ipart=1,Ntot
          Ang=pi2*rxpvec(ipart)
          x(ipart,1)=Radius*Sqrt(rxvec(ipart))*Sin(Ang)+Amean
          x(ipart,2)=Radius*Sqrt(rxvec(ipart))*Cos(Ang)+Bmean
c.......................................................................
          x2t=x(ipart,1)*x(ipart,1)
          x3t=x2t*x(ipart,1)
          x4t=x2t*x2t
          xp2t=x(ipart,2)*x(ipart,2)
          xp3t=xp2t*x(ipart,2)
          xp4t=xp2t*xp2t
          x2xpt=x2t*x(ipart,2)
          xxp2t=x(ipart,1)*xp2t
          x3xpt=x3t*x(ipart,2)
          xxp3t=x(ipart,1)*xp3t
          x2xp2t=x2t*xp2t
c.......................................................................
          xmean(1,1)=xmean(1,1)+x(ipart,1)                 !..mean/sigma
          xpmean(1,1)=xpmean(1,1)+x(ipart,2)               !..mean/sigma
          x2(1,1)=x2(1,1)+x2t
          xxp(1,1)=xxp(1,1)+x(ipart,1)*x(ipart,2)
          xp2(1,1)=xp2(1,1)+xp2t
          x3(1,1)=x3(1,1)+x3t
          x4(1,1)=x4(1,1)+x4t
          xp3(1,1)=xp3(1,1)+xp3t
          xp4(1,1)=xp4(1,1)+xp4t
          x2xp(1,1)=x2xp(1,1)+x2xpt
          xxp2(1,1)=xxp2(1,1)+xxp2t
          xxp3(1,1)=xxp3(1,1)+xxp3t
          x3xp(1,1)=x3xp(1,1)+x3xpt
          x2xp2(1,1)=x2xp2(1,1)+x2xp2t
c.......................................................................
          ipointer=((x(ipart,1)-frame(1,1))*Stpi(1))+1           !..dist
          Idist1D(ipointer,1,1)=Idist1D(ipointer,1,1)+1
          Id1Dmax(1)=Max(Id1Dmax(1),Idist1D(ipointer,1,1))
c.......................................................................
          cminx=Dmin1(cminx,x(ipart,1))
          cminxp=Dmin1(cminxp,x(ipart,2))
          cmaxx=Dmax1(cmaxx,x(ipart,1))
          cmaxxp=Dmax1(cmaxxp,x(ipart,2))
c.......................................................................
          Ipoin(ipart)=ipart            !..defines position of particles
c.......................................................................
        enddo
c.......................................................................
        xaxe=0.5d0*(cmaxx-cminx)       !..Computes surface using ellipse
        xpaxe=0.5d0*(cmaxxp-cminxp)
        ellipse(1,1)=0.5d0*pi2*xaxe*xpaxe*Dett
        surf(1,1)=Convex_hull(isw,1,Ntot,Ipoin)*Dett
c.......................................................................
        if (Idim.eq.4) then                           !..4D distribution
c.......................................................................
          lux=2
          iseed=7678443
          Call rluxgo(lux,iseed,0,0)      !..initializes random sequence
          Call ranlux(ryvec,NRNGen)           !..random number generator
          Call ranlux(rypvec,NRNGen)
c.......................................................................
          do ipart=1,Ntot
            Ang=pi2*rypvec(ipart)
            x(ipart,3)=Radv*Sqrt(ryvec(ipart))*Sin(Ang)+Ameav
            x(ipart,4)=Radv*Sqrt(ryvec(ipart))*Cos(Ang)+Bmeav
c.......................................................................
            y2t=x(ipart,3)*x(ipart,3)
            y3t=y2t*x(ipart,3)
            y4t=y2t*y2t
            yp2t=x(ipart,4)*x(ipart,4)
            yp3t=yp2t*x(ipart,4)
            yp4t=yp2t*yp2t
            y2ypt=y2t*x(ipart,4)
            yyp2t=x(ipart,3)*yp2t
            y3ypt=y3t*x(ipart,4)
            yyp3t=x(ipart,3)*yp3t
            y2yp2t=y2t*yp2t
c.......................................................................
            xmean(1,2)=xmean(1,2)+x(ipart,3)               !..mean/sigma
            xpmean(1,2)=xpmean(1,2)+x(ipart,4)             !..mean/sigma
            x2(1,2)=x2(1,2)+y2t
            xxp(1,2)=xxp(1,2)+x(ipart,3)*x(ipart,4)
            xp2(1,2)=xp2(1,2)+yp2t
            x3(1,2)=x3(1,2)+y3t
            x4(1,2)=x4(1,2)+y4t
            xp3(1,2)=xp3(1,2)+yp3t
            xp4(1,2)=xp4(1,2)+yp4t
            x2xp(1,2)=x2xp(1,2)+y2ypt
            xxp2(1,2)=xxp2(1,2)+yyp2t
            xxp3(1,2)=xxp3(1,2)+yyp3t
            x3xp(1,2)=x3xp(1,2)+y3ypt
            x2xp2(1,2)=x2xp2(1,2)+y2yp2t
c.......................................................................
            ipointer=((x(ipart,3)-frame(1,3))*Stpi(3))+1
            Idist1D(ipointer,1,2)=Idist1D(ipointer,1,2)+1
            Id1Dmax(2)=Max(Id1Dmax(2),Idist1D(ipointer,1,2))
c.......................................................................
            cminy=Dmin1(cminy,x(ipart,3))
            cminyp=Dmin1(cminyp,x(ipart,4))
            cmaxy=Dmax1(cmaxy,x(ipart,3))
            cmaxyp=Dmax1(cmaxyp,x(ipart,4))
c.......................................................................
            Ipoin(ipart)=ipart          !..defines position of particles
c.......................................................................
          enddo
c.......................................................................
          yaxe=0.5d0*(cmaxy-cminy)     !..Computes surface using ellipse
          ypaxe=0.5d0*(cmaxyp-cminyp)
          ellipse(1,2)=0.5d0*pi2*yaxe*ypaxe*Dett
          surf(1,2)=Convex_hull(isw,2,Ntot,Ipoin)*Dett
c.......................................................................
        endif
c.......................................................................
      elseif (Idist.eq.2) then                  !..Gaussian distribution
c.......................................................................
        lux=2                                         !..2D distribution
        iseed=7675039
        NRNGen=Ntot*(1+isw1*(Ninj-1))+Icentre*Nampc*Nanglesc
        Call rluxgo(lux,iseed,0,0)        !..initializes random sequence
        Call rnormx(rxvec,NRNGen,ranlux)      !..random number generator
        Call rnormx(rxpvec,NRNGen,ranlux)
c.......................................................................
        do ipart=1,Ntot
          x(ipart,1)=Sigma*rxvec(ipart)+Amean
          x(ipart,2)=Sigma*rxpvec(ipart)+Bmean
c.......................................................................
          x2t=x(ipart,1)*x(ipart,1)
          x3t=x2t*x(ipart,1)
          x4t=x2t*x2t
          xp2t=x(ipart,2)*x(ipart,2)
          xp3t=xp2t*x(ipart,2)
          xp4t=xp2t*xp2t
          x2xpt=x2t*x(ipart,2)
          xxp2t=x(ipart,1)*xp2t
          x3xpt=x3t*x(ipart,2)
          xxp3t=x(ipart,1)*xp3t
          x2xp2t=x2t*xp2t
c.......................................................................
          xmean(1,1)=xmean(1,1)+x(ipart,1)                 !..mean/sigma
          xpmean(1,1)=xpmean(1,1)+x(ipart,2)               !..mean/sigma
          x2(1,1)=x2(1,1)+x2t
          xxp(1,1)=xxp(1,1)+x(ipart,1)*x(ipart,2)
          xp2(1,1)=xp2(1,1)+xp2t
          x3(1,1)=x3(1,1)+x3t
          x4(1,1)=x4(1,1)+x4t
          xp3(1,1)=xp3(1,1)+xp3t
          xp4(1,1)=xp4(1,1)+xp4t
          x2xp(1,1)=x2xp(1,1)+x2xpt
          xxp2(1,1)=xxp2(1,1)+xxp2t
          xxp3(1,1)=xxp3(1,1)+xxp3t
          x3xp(1,1)=x3xp(1,1)+x3xpt
          x2xp2(1,1)=x2xp2(1,1)+x2xp2t
c.......................................................................
          ipointer=((x(ipart,1)-frame(1,1))*Stpi(1))+1           !..dist
          Idist1D(ipointer,1,1)=Idist1D(ipointer,1,1)+1
          Id1Dmax(1)=Max(Id1Dmax(1),Idist1D(ipointer,1,1))
c.......................................................................
          cminx=Dmin1(cminx,x(ipart,1))
          cminxp=Dmin1(cminxp,x(ipart,2))
          cmaxx=Dmax1(cmaxx,x(ipart,1))
          cmaxxp=Dmax1(cmaxxp,x(ipart,2))
c.......................................................................
          Ipoin(ipart)=ipart            !..defines position of particles
c.......................................................................
        enddo
c.......................................................................
        xaxe=0.5d0*(cmaxx-cminx)       !..Computes surface using ellipse
        xpaxe=0.5d0*(cmaxxp-cminxp)
        ellipse(1,1)=0.5d0*pi2*xaxe*xpaxe*Dett
        surf(1,1)=Convex_hull(isw,1,Ntot,Ipoin)*Dett
c.......................................................................
        if (Idim.eq.4) then                           !..4D distribution
c.......................................................................
          lux=2
          iseed=7678443
          Call rluxgo(lux,iseed,0,0)      !..initializes random sequence
          Call rnormx(ryvec,NRNGen,ranlux)    !..random number generator
          Call rnormx(rypvec,NRNGen,ranlux)
c.......................................................................
          do ipart=1,Ntot
            x(ipart,3)=Sigv*ryvec(ipart)+Ameav
            x(ipart,4)=Sigv*rypvec(ipart)+Bmeav
c.......................................................................
            y2t=x(ipart,3)*x(ipart,3)
            y3t=y2t*x(ipart,3)
            y4t=y2t*y2t
            yp2t=x(ipart,4)*x(ipart,4)
            yp3t=yp2t*x(ipart,4)
            yp4t=yp2t*yp2t
            y2ypt=y2t*x(ipart,4)
            yyp2t=x(ipart,3)*yp2t
            y3ypt=y3t*x(ipart,4)
            yyp3t=x(ipart,3)*yp3t
            y2yp2t=y2t*yp2t
c.......................................................................
            xmean(1,2)=xmean(1,2)+x(ipart,3)               !..mean/sigma
            xpmean(1,2)=xpmean(1,2)+x(ipart,4)             !..mean/sigma
            x2(1,2)=x2(1,2)+y2t
            xxp(1,2)=xxp(1,2)+x(ipart,3)*x(ipart,4)
            xp2(1,2)=xp2(1,2)+yp2t
            x3(1,2)=x3(1,2)+y3t
            x4(1,2)=x4(1,2)+y4t
            xp3(1,2)=xp3(1,2)+yp3t
            xp4(1,2)=xp4(1,2)+yp4t
            x2xp(1,2)=x2xp(1,2)+y2ypt
            xxp2(1,2)=xxp2(1,2)+yyp2t
            xxp3(1,2)=xxp3(1,2)+yyp3t
            x3xp(1,2)=x3xp(1,2)+y3ypt
            x2xp2(1,2)=x2xp2(1,2)+y2yp2t
c.......................................................................
            ipointer=((x(ipart,3)-frame(1,3))*Stpi(3))+1         !..dist
            Idist1D(ipointer,1,2)=Idist1D(ipointer,1,2)+1
            Id1Dmax(2)=Max(Id1Dmax(2),Idist1D(ipointer,1,2))
c.......................................................................
            cminy=Dmin1(cminy,x(ipart,3))
            cminyp=Dmin1(cminyp,x(ipart,4))
            cmaxy=Dmax1(cmaxy,x(ipart,3))
            cmaxyp=Dmax1(cmaxyp,x(ipart,4))
          enddo
c.......................................................................
          yaxe=0.5d0*(cmaxy-cminy)     !..Computes surface using ellipse
          ypaxe=0.5d0*(cmaxyp-cminyp)
          ellipse(1,2)=0.5d0*pi2*yaxe*ypaxe*Dett
          surf(1,2)=Convex_hull(isw,2,Ntot,Ipoin)*Dett
c.......................................................................
        endif
c.......................................................................
      elseif (Idist.eq.3) then         !..Hollow distribution in H and V
c.......................................................................
        nsigma=4
        ainf=Dmax1(0d0,Amean-nsigma*Sigma)
        asup=Amean+nsigma*Sigma
        Sigma2=Sigma*Sigma
        Amean2=Amean*Amean
        sigmafact=1d0/(8d0*Sigma2)
        Sigmafact1=1d0/(2d0*Sigma2)
        fmax=0.5d0*(Amean+Dsqrt(Amean2+4d0*Sigma2))*
     .       Dexp(-(Amean-Dsqrt(Amean2+4d0*Sigma2))*
     .                       (Amean-Dsqrt(Amean2+4d0*Sigma2))*sigmafact)
        fmaxinv=1d0/fmax
        surfmax=fmax*(asup-ainf)
c.......................................................................
        kpart=0
c.......................................................................
 10     NRNGen=Ntot*(1+isw1*(Ninj-1))+Icentre*Nampc*Nanglesc
        lux=2                                         !..2D distribution
        iseed=7675039+kpart
        Call rluxgo(lux,iseed,0,0)        !..initializes random sequence
        Call ranlux(rxvec,NRNGen)             !..random number generator
        Call ranlux(rxpvec,NRNGen)
        Call ranlux(radxv,NRNGen)
c.......................................................................
        do ipart=1,NRNGen
          radius=rxvec(ipart)*surfmax*fmaxinv+ainf
          Deltar=radius-Amean
          fdistrib=radius*Dexp(-Deltar*Deltar*sigmafact1)
          if (rxpvec(ipart)*fmax.lt.fdistrib) then
c.......................................................................
            kpart=kpart+1
            x(kpart,1)=radius*Dcos(pi2*radxv(kpart))
            x(kpart,2)=radius*Dsin(pi2*radxv(kpart))
c.......................................................................
            x2t=x(kpart,1)*x(kpart,1)
            x3t=x2t*x(kpart,1)
            x4t=x2t*x2t
            xp2t=x(kpart,2)*x(kpart,2)
            xp3t=xp2t*x(kpart,2)
            xp4t=xp2t*xp2t
            x2xpt=x2t*x(kpart,2)
            xxp2t=x(kpart,1)*xp2t
            x3xpt=x3t*x(kpart,2)
            xxp3t=x(kpart,1)*xp3t
            x2xp2t=x2t*xp2t
c.......................................................................
            xmean(1,1)=xmean(1,1)+x(kpart,1)               !..mean/sigma
            xpmean(1,1)=xpmean(1,1)+x(kpart,2)             !..mean/sigma
            x2(1,1)=x2(1,1)+x2t
            xxp(1,1)=xxp(1,1)+x(kpart,1)*x(kpart,2)
            xp2(1,1)=xp2(1,1)+xp2t
            x3(1,1)=x3(1,1)+x3t
            x4(1,1)=x4(1,1)+x4t
            xp3(1,1)=xp3(1,1)+xp3t
            xp4(1,1)=xp4(1,1)+xp4t
            x2xp(1,1)=x2xp(1,1)+x2xpt
            xxp2(1,1)=xxp2(1,1)+xxp2t
            xxp3(1,1)=xxp3(1,1)+xxp3t
            x3xp(1,1)=x3xp(1,1)+x3xpt
            x2xp2(1,1)=x2xp2(1,1)+x2xp2t
c.......................................................................
            ipointer=((x(kpart,1)-frame(1,1))*Stpi(1))+1         !..dist
            Idist1D(ipointer,1,1)=Idist1D(ipointer,1,1)+1
            Id1Dmax(1)=Max(Id1Dmax(1),Idist1D(ipointer,1,1))
c.......................................................................
            cminx=Dmin1(cminx,x(kpart,1))
            cminxp=Dmin1(cminxp,x(kpart,2))
            cmaxx=Dmax1(cmaxx,x(kpart,1))
            cmaxxp=Dmax1(cmaxxp,x(kpart,2))
c.......................................................................
            Ipoin(kpart)=kpart          !..defines position of particles
c.......................................................................
          endif
c.......................................................................
          if (kpart.gt.Ntot) goto 20
c.......................................................................
        enddo
c.......................................................................
        if (kpart.lt.Ntot) goto 10
c.......................................................................
 20     xaxe=0.5d0*(cmaxx-cminx)       !..Computes surface using ellipse
        xpaxe=0.5d0*(cmaxxp-cminxp)
        ellipse(1,1)=0.5d0*pi2*xaxe*xpaxe*Dett
        surf(1,1)=Convex_hull(isw,1,Ntot,Ipoin)*Dett
c.......................................................................
        if (Idim.eq.4) then                           !..4D distribution
c.......................................................................
          ainfv=Dmax1(0d0,Ameav-nsigma*Sigv)
          asupv=Ameav+nsigma*Sigv 
          Sigv2=Sigv*Sigv
          Ameav2=Ameav*Ameav
          sigvfact=1d0/(8d0*Sigv2)
          sigvfact1=1d0/(2d0*Sigv2)
          fmaxv=0.5d0*(Ameav+Dsqrt(Ameav2+4d0*Sigv2))*
     .          Dexp(-(Ameav-Dsqrt(Ameav2+4d0*Sigv2))*
     .                         (Ameav-Dsqrt(Ameav2+4d0*Sigv2))*sigvfact)
          fmaxvinv=1d0/fmaxv
          surfmaxv=fmaxv*(asupv-ainfv)
c.......................................................................
          lpart=0
c.......................................................................
 30       lux=2
          iseed=7678443+lpart
          Call rluxgo(lux,iseed,0,0) !..initializes random sequence
          Call ranlux(ryvec,NRNGen)           !..random number generator
          Call ranlux(rypvec,NRNGen)
          Call ranlux(radyv,NRNGen)
c.......................................................................
          do ipart=1,NRNGen
            radius=ryvec(ipart)*surfmaxv*fmaxvinv+ainfv
            Deltar=radius-Ameav
            fdistrib=radius*Dexp(-Deltar*Deltar*sigvfact1)
            if (rypvec(ipart)*fmaxv.lt.fdistrib) then
c.......................................................................
              lpart=lpart+1
              x(lpart,3)=radius*Dcos(pi2*radyv(lpart))
              x(lpart,4)=radius*Dsin(pi2*radyv(lpart))
c.......................................................................
              y2t=x(lpart,3)*x(lpart,3)
              y3t=y2t*x(lpart,3)
              y4t=y2t*y2t
              yp2t=x(lpart,4)*x(lpart,4)
              yp3t=yp2t*x(lpart,4)
              yp4t=yp2t*yp2t
              y2ypt=y2t*x(lpart,4)
              yyp2t=x(lpart,3)*yp2t
              y3ypt=y3t*x(lpart,4)
              yyp3t=x(lpart,3)*yp3t
              y2yp2t=y2t*yp2t
c.......................................................................
              xmean(1,2)=xmean(1,2)+x(lpart,3)             !..mean/sigma
              xpmean(1,2)=xpmean(1,2)+x(lpart,4)           !..mean/sigma
              x2(1,2)=x2(1,2)+y2t
              xxp(1,2)=xxp(1,2)+x(lpart,3)*x(lpart,4)
              xp2(1,2)=xp2(1,2)+yp2t
              x3(1,2)=x3(1,2)+y3t
              x4(1,2)=x4(1,2)+y4t
              xp3(1,2)=xp3(1,2)+yp3t
              xp4(1,2)=xp4(1,2)+yp4t
              x2xp(1,2)=x2xp(1,2)+y2ypt
              xxp2(1,2)=xxp2(1,2)+yyp2t
              xxp3(1,2)=xxp3(1,2)+yyp3t
              x3xp(1,2)=x3xp(1,2)+y3ypt
              x2xp2(1,2)=x2xp2(1,2)+y2yp2t
c.......................................................................
              ipointer=((x(lpart,3)-frame(1,3))*Stpi(3))+1       !..dist
              Idist1D(ipointer,1,2)=Idist1D(ipointer,1,2)+1
              Id1Dmax(2)=Max(Id1Dmax(2),Idist1D(ipointer,1,2))
c.......................................................................
              cminy=Dmin1(cminy,x(lpart,3))
              cminyp=Dmin1(cminyp,x(lpart,4))
              cmaxy=Dmax1(cmaxy,x(lpart,3))
              cmaxyp=Dmax1(cmaxyp,x(lpart,4))
c.......................................................................
            endif
c.......................................................................
            if (lpart.gt.Ntot) goto 40
          enddo
c.......................................................................
          if (lpart.lt.Ntot) goto 30
c.......................................................................
 40       yaxe=0.5d0*(cmaxy-cminy)     !..Computes surface using ellipse
          ypaxe=0.5d0*(cmaxyp-cminyp)
          ellipse(1,2)=0.5d0*pi2*yaxe*ypaxe*Dett
          surf(1,2)=Convex_hull(isw,2,Ntot,Ipoin)*Dett
c.......................................................................
        endif
c.......................................................................
      elseif (Idist.eq.4) then    !..Hollow dist. in H and Gaussian in V
c.......................................................................
        nsigma=4
        ainf=Dmax1(0d0,Amean-nsigma*Sigma)
        asup=Amean+nsigma*Sigma
        Sigma2=Sigma*Sigma
        Amean2=Amean*Amean
        sigmafact=1d0/(8d0*Sigma2)
        Sigmafact1=1d0/(2d0*Sigma2)
        fmax=0.5d0*(Amean+Dsqrt(Amean2+4d0*Sigma2))*
     .       Dexp(-(Amean-Dsqrt(Amean2+4d0*Sigma2))*
     .                       (Amean-Dsqrt(Amean2+4d0*Sigma2))*sigmafact)
        fmaxinv=1d0/fmax
        surfmax=fmax*(asup-ainf)
c.......................................................................
        kpart=0
c.......................................................................
 50     NRNGen=Ntot*(1+isw1*(Ninj-1))+Icentre*Nampc*Nanglesc
        lux=2                                         !..2D distribution
        iseed=7675039+kpart
        Call rluxgo(lux,iseed,0,0)        !..initializes random sequence
        Call ranlux(rxvec,NRNGen)             !..random number generator
        Call ranlux(rxpvec,NRNGen)
        Call ranlux(radxv,NRNGen)
c.......................................................................
        do ipart=1,NRNGen
          radius=rxvec(ipart)*surfmax*fmaxinv+ainf
          Deltar=radius-Amean
          fdistrib=radius*Dexp(-Deltar*Deltar*sigmafact1)
          if (rxpvec(ipart)*fmax.lt.fdistrib) then
c.......................................................................
            kpart=kpart+1
            x(kpart,1)=radius*Dcos(pi2*radxv(kpart))
            x(kpart,2)=radius*Dsin(pi2*radxv(kpart))
c.......................................................................
            x2t=x(kpart,1)*x(kpart,1)
            x3t=x2t*x(kpart,1)
            x4t=x2t*x2t
            xp2t=x(kpart,2)*x(kpart,2)
            xp3t=xp2t*x(kpart,2)
            xp4t=xp2t*xp2t
            x2xpt=x2t*x(kpart,2)
            xxp2t=x(kpart,1)*xp2t
            x3xpt=x3t*x(kpart,2)
            xxp3t=x(kpart,1)*xp3t
            x2xp2t=x2t*xp2t
c.......................................................................
            xmean(1,1)=xmean(1,1)+x(kpart,1)               !..mean/sigma
            xpmean(1,1)=xpmean(1,1)+x(kpart,2)             !..mean/sigma
            x2(1,1)=x2(1,1)+x2t
            xxp(1,1)=xxp(1,1)+x(kpart,1)*x(kpart,2)
            xp2(1,1)=xp2(1,1)+xp2t
            x3(1,1)=x3(1,1)+x3t
            x4(1,1)=x4(1,1)+x4t
            xp3(1,1)=xp3(1,1)+xp3t
            xp4(1,1)=xp4(1,1)+xp4t
            x2xp(1,1)=x2xp(1,1)+x2xpt
            xxp2(1,1)=xxp2(1,1)+xxp2t
            xxp3(1,1)=xxp3(1,1)+xxp3t
            x3xp(1,1)=x3xp(1,1)+x3xpt
            x2xp2(1,1)=x2xp2(1,1)+x2xp2t
c.......................................................................
            ipointer=((x(kpart,1)-frame(1,1))*Stpi(1))+1         !..dist
            Idist1D(ipointer,1,1)=Idist1D(ipointer,1,1)+1
            Id1Dmax(1)=Max(Id1Dmax(1),Idist1D(ipointer,1,1))
c.......................................................................
            cminx=Dmin1(cminx,x(kpart,1))
            cminxp=Dmin1(cminxp,x(kpart,2))
            cmaxx=Dmax1(cmaxx,x(kpart,1))
            cmaxxp=Dmax1(cmaxxp,x(kpart,2))
c.......................................................................
            Ipoin(kpart)=kpart          !..defines position of particles
c.......................................................................
          endif
c.......................................................................
          if (kpart.gt.Ntot) goto 60
c.......................................................................
        enddo
c.......................................................................
        if (kpart.lt.Ntot) goto 50
c.......................................................................
 60     xaxe=0.5d0*(cmaxx-cminx)       !..Computes surface using ellipse
        xpaxe=0.5d0*(cmaxxp-cminxp)
        ellipse(1,1)=0.5d0*pi2*xaxe*xpaxe*Dett
        surf(1,1)=Convex_hull(isw,1,Ntot,Ipoin)*Dett
c.......................................................................
        if (Idim.eq.4) then                           !..4D distribution
c.......................................................................
          lux=2
          iseed=7678443
          Call rluxgo(lux,iseed,0,0)      !..initializes random sequence
          Call rnormx(ryvec,NRNGen,ranlux)    !..random number generator
          Call rnormx(rypvec,NRNGen,ranlux)
c.......................................................................
          do ipart=1,Ntot
            x(ipart,3)=Sigv*ryvec(ipart)+Ameav
            x(ipart,4)=Sigv*rypvec(ipart)+Bmeav
c.......................................................................
            y2t=x(ipart,3)*x(ipart,3)
            y3t=y2t*x(ipart,3)
            y4t=y2t*y2t
            yp2t=x(ipart,4)*x(ipart,4)
            yp3t=yp2t*x(ipart,4)
            yp4t=yp2t*yp2t
            y2ypt=y2t*x(ipart,4)
            yyp2t=x(ipart,3)*yp2t
            y3ypt=y3t*x(ipart,4)
            yyp3t=x(ipart,3)*yp3t
            y2yp2t=y2t*yp2t
c.......................................................................
            xmean(1,2)=xmean(1,2)+x(ipart,3)               !..mean/sigma
            xpmean(1,2)=xpmean(1,2)+x(ipart,4)             !..mean/sigma
            x2(1,2)=x2(1,2)+y2t
            xxp(1,2)=xxp(1,2)+x(ipart,3)*x(ipart,4)
            xp2(1,2)=xp2(1,2)+yp2t
            x3(1,2)=x3(1,2)+y3t
            x4(1,2)=x4(1,2)+y4t
            xp3(1,2)=xp3(1,2)+yp3t
            xp4(1,2)=xp4(1,2)+yp4t
            x2xp(1,2)=x2xp(1,2)+y2ypt
            xxp2(1,2)=xxp2(1,2)+yyp2t
            xxp3(1,2)=xxp3(1,2)+yyp3t
            x3xp(1,2)=x3xp(1,2)+y3ypt
            x2xp2(1,2)=x2xp2(1,2)+y2yp2t
c.......................................................................
            ipointer=((x(ipart,3)-frame(1,3))*Stpi(3))+1         !..dist
            Idist1D(ipointer,1,2)=Idist1D(ipointer,1,2)+1
            Id1Dmax(2)=Max(Id1Dmax(2),Idist1D(ipointer,1,2))
c.......................................................................
            cminy=Dmin1(cminy,x(ipart,3))
            cminyp=Dmin1(cminyp,x(ipart,4))
            cmaxy=Dmax1(cmaxy,x(ipart,3))
            cmaxyp=Dmax1(cmaxyp,x(ipart,4))
          enddo
c.......................................................................
          yaxe=0.5d0*(cmaxy-cminy)     !..Computes surface using ellipse
          ypaxe=0.5d0*(cmaxyp-cminyp)
          ellipse(1,2)=0.5d0*pi2*yaxe*ypaxe*Dett
          surf(1,2)=Convex_hull(isw,2,Ntot,Ipoin)*Dett
c.......................................................................
        endif
c.......................................................................
      endif
c.......................................................................
      do jdim=1,Idim/2                           !..loop over H/V planes
c.....................................................transforms moments
        avg(1)=Tmat(1,1,jdim)*xmean(1,jdim)+Tmat(1,2,jdim)*
     .         xpmean(1,jdim)
        avg(2)=Tmat(2,1,jdim)*xmean(1,jdim)+Tmat(2,2,jdim)*
     .         xpmean(1,jdim)
        xmean(1,jdim)=avg(1)/Dfloat(Ntot)                        !..mean
        xpmean(1,jdim)=avg(2)/Dfloat(Ntot)                       !..mean
c.......................................................................
        smom(1)=Tmom(1,1,jdim)*x2(1,jdim)+Tmom(1,2,jdim)*xxp(1,jdim)+
     .          Tmom(1,3,jdim)*xp2(1,jdim)
        smom(2)=Tmom(2,1,jdim)*x2(1,jdim)+Tmom(2,2,jdim)*xxp(1,jdim)+
     .          Tmom(2,3,jdim)*xp2(1,jdim)
        smom(3)=Tmom(3,1,jdim)*x2(1,jdim)+Tmom(3,2,jdim)*xxp(1,jdim)+
     .          Tmom(3,3,jdim)*xp2(1,jdim)
c.................higher-order moments are not transformed back!!!!!!!!!
        tmpm3=x3(1,jdim)/Dfloat(Ntot)
        tmpm4=x4(1,jdim)/Dfloat(Ntot)
        x3(1,jdim)=tmpm3-3d0*x2(1,jdim)/Dfloat(Ntot)*xmean(1,jdim)+
     .                     2d0*xmean(1,jdim)*xmean(1,jdim)*xmean(1,jdim)
        x4(1,jdim)=tmpm4+6d0*x2(1,jdim)/Dfloat(Ntot)*xmean(1,jdim)*
     .                       xmean(1,jdim)-
     .                   4d0*tmpm3*xmean(1,jdim)-3d0*xmean(1,jdim)*
     .                       xmean(1,jdim)*xmean(1,jdim)*xmean(1,jdim)
        dtest=x2(1,jdim)/Dfloat(Ntot)-xmean(1,jdim)*xmean(1,jdim)
        if (dtest.ne.0d0) halop1(1,jdim)=x4(1,jdim)/(dtest*dtest)-2d0
c.......................................................................
        tmpmp3=xp3(1,jdim)/Dfloat(Ntot)
        tmpmp4=xp4(1,jdim)/Dfloat(Ntot)
        tmpm2p1=x2xp(1,jdim)/Dfloat(Ntot)
        tmpm1p2=xxp2(1,jdim)/Dfloat(Ntot)
        tmpm1p3=xxp3(1,jdim)/Dfloat(Ntot)
        tmpm3p1=x3xp(1,jdim)/Dfloat(Ntot)
        tmpm2p2=x2xp2(1,jdim)/Dfloat(Ntot)
        x2xp(1,jdim)=tmpm2p1-x2(1,jdim)/Dfloat(Ntot)*xpmean(1,jdim)-
     .               2d0*xxp(1,jdim)/Dfloat(Ntot)*xmean(1,jdim)+
     .               2d0*xmean(1,jdim)*xmean(1,jdim)*xpmean(1,jdim)
        xxp2(1,jdim)=tmpm1p2-xp2(1,jdim)/Dfloat(Ntot)*xmean(1,jdim)-
     .               2d0*xxp(1,jdim)/Dfloat(Ntot)*xpmean(1,jdim)+
     .               2d0*xpmean(1,jdim)*xpmean(1,jdim)*xmean(1,jdim)
        xp3(1,jdim)=tmpmp3-3d0*xp2(1,jdim)/Dfloat(Ntot)*xpmean(1,jdim)+
     .                  2d0*xpmean(1,jdim)*xpmean(1,jdim)*xpmean(1,jdim)
        xp4(1,jdim)=tmpmp4+6d0*xp2(1,jdim)/Dfloat(Ntot)*xpmean(1,jdim)*
     .      xpmean(1,jdim)-4d0*tmpmp3*xpmean(1,jdim)-3d0*xpmean(1,jdim)*
     .                      xpmean(1,jdim)*xpmean(1,jdim)*xpmean(1,jdim)
        xxp3(1,jdim)=tmpm1p3-xmean(1,jdim)*tmpmp3-3d0*xpmean(1,jdim)*
     .     tmpm1p2+3d0*xmean(1,jdim)*xpmean(1,jdim)*xp2(1,jdim)/
     .     Dfloat(Ntot)+3d0*xpmean(1,jdim)*xpmean(1,jdim)*xxp(1,jdim)/
     .     Dfloat(Ntot)-3d0*xmean(1,jdim)*xpmean(1,jdim)*xpmean(1,jdim)*
     .     xpmean(1,jdim)
        x3xp(1,jdim)=tmpm3p1-xpmean(1,jdim)*tmpm3-3d0*xmean(1,jdim)*
     .     tmpm2p1+3d0*xpmean(1,jdim)*xmean(1,jdim)*x2(1,jdim)/
     .     Dfloat(Ntot)+3d0*xmean(1,jdim)*xmean(1,jdim)*xxp(1,jdim)/
     .     Dfloat(Ntot)-3d0*xpmean(1,jdim)*xmean(1,jdim)*xmean(1,jdim)*
     .     xmean(1,jdim)
        x2xp2(1,jdim)=tmpm2p2-2d0*tmpm2p1*xpmean(1,jdim)-2d0*tmpm1p2*
     .     xmean(1,jdim)+x2(1,jdim)/Dfloat(Ntot)*xpmean(1,jdim)*
     .     xpmean(1,jdim)+4d0*xmean(1,jdim)*xpmean(1,jdim)*xxp(1,jdim)/
     .     Dfloat(Ntot)+xp2(1,jdim)/Dfloat(Ntot)*xmean(1,jdim)*
     .     xmean(1,jdim)-3d0*xmean(1,jdim)*xmean(1,jdim)*xpmean(1,jdim)*
     .     xpmean(1,jdim)
c.......................................................................
        x2(1,jdim)=smom(1)/Dfloat(Ntot)-xmean(1,jdim)*xmean(1,jdim)
        xxp(1,jdim)=smom(2)/Dfloat(Ntot)-xmean(1,jdim)*xpmean(1,jdim)
        xp2(1,jdim)=smom(3)/Dfloat(Ntot)-xpmean(1,jdim)*xpmean(1,jdim)
c.......................................................................
        etest=x2(1,jdim)*xp2(1,jdim)-xxp(1,jdim)*xxp(1,jdim)
        if (etest.gt.0d0) emit(1,jdim)=Dsqrt(etest)
        beta(1,jdim)=x2(1,jdim)/emit(1,jdim)
        alpha(1,jdim)=xxp(1,jdim)/emit(1,jdim)
        ftest=3d0*(x4(1,jdim)*xp4(1,jdim)+3d0*x2xp2(1,jdim)*
     .             x2xp2(1,jdim)-4d0*xxp3(1,jdim)*x3xp(1,jdim))
        if (ftest.gt.0d0) halop2(1,jdim)=dsqrt(ftest)/(2d0*etest)-2d0
c.......................................................................
      enddo
c.......................................................................
      if ((Icomp.eq.7).or.(Icomp.eq.8)) then     !..multi-turn injection
c.......................................................................
        if (Ianim.eq.0) then
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        do kinj=1,Ninj-1
c.......................................................................
          Ntot=Namp*Nangles*kinj  !..redefines total number of particles
          iflag=2
          Call Comp_dist(iflag,kinj)          !..computes 2D beam distr.
c.......................................................................
          if (Idim.eq.2) then
            do ipart=1,Ntot
              Itest=Iuser_map2(ipart,1)
            enddo
          elseif (Idim.eq.4) then
            do ipart=1,Ntot
              Itest=Iuser_map4(ipart,1)
            enddo
          endif
c.......................................................................
          if (Idist.eq.1) then                   !..uniform distribution
c.......................................................................
            do ipart=1,Ntot
              Ang=pi2*rxpvec(ipart)
              x(ipart+Ntot,1)=Radius*Sqrt(rxvec(ipart))*Sin(Ang)+Amean
              x(ipart+Ntot,2)=Radius*Sqrt(rxvec(ipart))*Cos(Ang)+Bmean
c.......................................................................
              Ipoin(ipart+Ntot)=ipart+Ntot            !..particles' pos.
c.......................................................................
            enddo
c.......................................................................
            if (Idim.eq.4) then                       !..4D distribution
c.......................................................................
              do ipart=1,Ntot
                Ang=pi2*rypvec(ipart)
                x(ipart+Ntot,3)=Radv*Sqrt(ryvec(ipart))*Sin(Ang)+Ameav
                x(ipart+Ntot,4)=Radv*Sqrt(ryvec(ipart))*Cos(Ang)+Bmeav
c.......................................................................
              enddo
c.......................................................................
            endif
c.......................................................................
          elseif (Idist.eq.2) then              !..Gaussian distribution
c.......................................................................
            do ipart=1,Namp*Nangles
              x(ipart+Ntot,1)=Sigma*rxvec(ipart+Ntot)+Amean
              x(ipart+Ntot,2)=Sigma*rxpvec(ipart+Ntot)+Bmean
c.......................................................................
              Ipoin(ipart+Ntot)=ipart+Ntot         !..defines part. pos.
c.......................................................................
            enddo
c.......................................................................
            if (Idim.eq.4) then                       !..4D distribution
c.......................................................................
              do ipart=1,Namp*Nangles
                x(ipart+Ntot,3)=Sigv*ryvec(ipart+Ntot)+Ameav
                x(ipart+Ntot,4)=Sigv*rypvec(ipart+Ntot)+Bmeav
              enddo
c.......................................................................
            endif
c.......................................................................
          endif
c.......................................................................
        enddo
c.......................................................................
        Ntot=Namp*Nangles*Ninj    !..redefines total number of particles
c.......................................................................
        if (Icentre.eq.1) then !..injects beam in the centre
c.......................................................................
          if (Idist.eq.1) then                   !..uniform distribution
c.......................................................................
            do iamp=1,Nampc                           !..2D distribution
              Amp=Sigmac*Dfloat(iamp)
              do iang=1,Nanglesc
                Ang=Astep*Dfloat(iang)
                ipart=(iamp-1)*Nampc+iang
                x(ipart+Ntot,1)=Amp*Dcos(Ang)+Ameanc
                x(ipart+Ntot,2)=Amp*Dsin(Ang)+Bmeanc
c.......................................................................
                Ipoin(ipart+Ntot)=ipart+Ntot          !..particles' pos.
c.......................................................................
              enddo
            enddo
c.......................................................................
            if (Idim.eq.4) then                       !..4D distribution
c.......................................................................
              do iamp=1,Nampc
                Amp=Sigvc*Dfloat(iamp)
                do iang=1,Nanglesc
                  Ang=Astep*Dfloat(iang)
                  ipart=(iamp-1)*Nampc+iang
                  x(ipart+Ntot,3)=Amp*Dcos(Ang)+Ameanvc
                  x(ipart+Ntot,4)=Amp*Dsin(Ang)+Bmeanvc
c.......................................................................
                enddo
c.......................................................................
              enddo
c.......................................................................
            endif
c.......................................................................
          elseif (Idist.eq.2) then              !..Gaussian distribution
c.......................................................................
            do ipart=1,Nampc*Nanglesc
              x(ipart+Ntot,1)=Sigmac*rxvec(ipart+Ntot)+Ameanc
              x(ipart+Ntot,2)=Sigmac*rxpvec(ipart+Ntot)+Bmeanc
c.......................................................................
              Ipoin(ipart+Ntot)=ipart+Ntot         !..defines part. pos.
c.......................................................................
            enddo
c.......................................................................
            if (Idim.eq.4) then                       !..4D distribution
c.......................................................................
              do ipart=1,Nampc*Nanglesc
                x(ipart+Ntot,3)=Sigvc*ryvec(ipart+Ntot)+Ameanc
                x(ipart+Ntot,4)=Sigvc*rypvec(ipart+Ntot)+Bmeanc
              enddo
c.......................................................................
            endif
c.......................................................................
          endif
c.......................................................................
        endif
c.......................................................................
        Ntot=Ntot+Icentre*Nampc*Nanglesc        !..tot numb of particles
        iflag=2
        Call Comp_dist(iflag,Ninj+Icentre)    !..computes 2D beam distr.
c.......................................................................
        do ipart=1,Ntot-Namp*Nangles*(1-Icentre)    !..computes profiles
c.......................................................................
          ipointer=((x(ipart,1)-frame(1,1))*Stpi(1))+1
          Idist1D(ipointer,1,1)=Idist1D(ipointer,1,1)+1
          Id1Dmax(1)=Max(Id1Dmax(1),Idist1D(ipointer,1,1))
c.......................................................................
          if (Idim.eq.4) then
            ipointer=((x(ipart,3)-frame(1,3))*Stpi(3))+1
            Idist1D(ipointer,1,2)=Idist1D(ipointer,1,2)+1
            Id1Dmax(2)=Max(Id1Dmax(2),Idist1D(ipointer,1,2))
          endif
c.......................................................................
        enddo
      endif
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Gen_boxes
c.......................................................................
c.... Subroutine to generate the boxes used in the computation of the
c.... beam distribution.
c....
c.... Author: M. Giovannozzi
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      Data xpmin,ypmin/2*1d38/
      Data xpmax,ypmax/2*-1d38/
c.......................................................................
      bound(1,1)=bbox(1)                             !..generates b. box
      bound(2,1)=bbox(2)
      bound(3,1)=bbox(3)
      bound(4,1)=bbox(4)
      Stepx(1)=(bound(3,1)-bound(1,1))/Dfloat(Nbin) !..step in spec. C-S
c.......................................................................
      do ipoint=1,4    !..transforms bounding box in proper co-ordinates
        anew=ipoint-1.9           !..fancy def. for proper coor. couples
        ianew=anew/2              !..fancy def. for proper coor. couples
        Isgn=(anew+Abs(anew))/anew-1                    !..sign function
        Ivar1=ianew-(Isgn-1)/2    !..fancy def. for proper coor. couples
        xtransf(1)=bound(1,1)*Ivar1+bound(3,1)*(1-Ivar1)      !..H plane
        Ivar2=ipoint/3+1
        xtransf(2)=bound(2,1)*(2-Ivar2)+bound(4,1)*(Ivar2-1)  !..H plane
        xtransf(3)=xtransf(1)                                 !..V plane
        xtransf(4)=xtransf(2)                                 !..V plane
        Call Transf                 !..transforms corners of initial box
        xpmin=min(xpmin,xtransf(2))        !..finds min angle in H plane
        xpmax=max(xpmax,xtransf(2))        !..finds max angle in H plane
        ypmin=min(ypmin,xtransf(4))        !..finds min angle in V plane
        ypmax=max(ypmax,xtransf(4))        !..finds max angle in V plane
      enddo
c.......................................................................
      frame(1,1)=Tmat(1,1,1)*bound(1,1)!..defines min frame limits (H/V) 
      frame(1,2)=xpmin
      frame(1,3)=Tmat(1,1,2)*bound(1,1)
      frame(1,4)=ypmin
c.......................................................................
      frame(2,1)=Tmat(1,1,1)*bound(3,1)!..defines max frame limits (H/V) 
      frame(2,2)=xpmax
      frame(2,3)=Tmat(1,1,2)*bound(3,1)
      frame(2,4)=ypmax
c.......................................................................
      fmul=1d0
      if (Ico.eq.1) fmul=1d1
      if (Ico.gt.0) then
        do ii=1,2
          do jj=1,Max_dim
            frame(ii,jj)=Int(frame(ii,jj)*fmul)/fmul  !..truncates scale
          enddo
        enddo
      endif
c.......................................................................
      do ii=1,4
        Stp(ii)=(frame(2,ii)-frame(1,ii))/
     .          Dfloat(Nbin)          !..transformed steps for Comp_dist
        Stpi(ii)=1d0/Stp(ii)                             !..inverse step
      enddo
c.......................................................................
      Deltax=(bbox(3)-bbox(1))/Dfloat(Nbx)   !..steps x, x' in spec. C-S
      Deltay=(bbox(4)-bbox(2))/Dfloat(Nby)
c.......................................................................
      do jbox=1,Nboxes(1)
c.......................................................................
        irow=Ibox(jbox)/Nbx+1                       !..finds coordinates
        if (mod(Ibox(jbox),Nbx).eq.0) irow=irow-1
        icol=Ibox(jbox)-(irow-1)*Nbx
c.......................................................................
        bound(1,jbox+1)=bbox(1)+(icol-1)*Deltax      !..generates b. box
        bound(2,jbox+1)=bbox(4)-irow*Deltay
        bound(3,jbox+1)=bbox(1)+icol*Deltax
        bound(4,jbox+1)=bbox(4)-(irow-1)*Deltay
        Stepx(jbox+1)=(bound(3,jbox+1)-bound(1,jbox+1))/Dfloat(Nbin)
c.......................................................................
      enddo
c.......................................................................
      do jcomp=1,4                     !..defines border for initial box
        do jdim=1,2
          bord(jcomp,jdim,1)=bound(jcomp,1)
        enddo
      enddo
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Comp_dist(iflag,iturn)
c.......................................................................
c.... Subroutine to compute the beam distribution. It computes the mean
c.... and sigma of the different beam parts. 
c....
c.... iflag = 0: Final beam distribution. The beam core is used to 
c....            compute the beam distribution. The results are printed 
c....            out
c.... iflag = 1: Final beam distribution. The different beam pieces are
c....            used to compute the different beam distributions. The
c....            results are printed out
c.... iflag = 2: Histogram of the beam distribution is computed for 
c....            graphical output. The results are printed out
c.... iflag = 3: Histogram of the beam distribution is computed for 
c....            graphical output. The results are printed out. However, 
c....            islands are removed from the final distribution (for 
c....            multi-extraction purposes).
c.... iflag = 4: Same as 2 but to be used in case of multi-extractions.
c....            It considers that only a subset (specified by the array
c....            Iarray) of the initial particles is to be used in the 
c....            computations. The results are printed out
c....
c.... iturn: turn number
c....
c.... Author: M. Giovannozzi
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Dimension bordc(2,Max_dim,Max_boxes),bordw(2,Max_dim,Max_boxes)
      Dimension Idist2D(Max_bin,Max_bin,Max_pro),Npart(Max_boxes)
      Dimension avg(2),smom(3)
      Dimension Ipoin(Max_part),Itemp11(0:Max_bord),Itemp12(0:Max_bord),
     .  iptest(2),Itemp21(Max_bord),Itemp22(Max_bord),Iscratch(Max_part)
c.......................................................................
      Real rxvec(Max_part),rxpvec(Max_part)
      Real ryvec(Max_part),rypvec(Max_part)
      Real hmal,xc,yc,zc,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*3 ext
      Character*4 cmap
      Character*12 headmom(2,14)
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
      Character*100 command
c.......................................................................
      Logical Pcheck,inside
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Momen/xmean(Max_boxes,2),xpmean(Max_boxes,2),
     .             x2(Max_boxes,2),xp2(Max_boxes,2),xxp(Max_boxes,2),
     .             x3(Max_boxes,2),xp3(Max_boxes,2),x2xp(Max_boxes,2),
     .             xxp2(Max_boxes,2),x4(Max_boxes,2),xp4(Max_boxes,2),
     .             x2xp2(Max_boxes,2),xxp3(Max_boxes,2),
     .             x3xp(Max_boxes,2),
     .             emit(Max_boxes,2),beta(Max_boxes,2),
     .             alpha(Max_boxes,2),ellipse(Max_boxes,Max_dim),
     .             surf(Max_boxes,Max_dim),halop1(Max_boxes,2),
     .             halop2(Max_boxes,2),Idist1D(Max_bin,Max_boxes,2)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/Distribution/Radius,Sigma,Sigmac,Amean,Bmean,Radv,
     .                    Sigv,Sigvc,Ameav,Bmeav,Ameanc,Bmeanc,
     .                    Ameanvc,Bmeanvc,xsep,Sigman,Sigmah(Max_bord),
     .                    Akick(Max_bord),rxvec,rxpvec,
     .                    ryvec,rypvec,Namp,Nangles,Nampc,Nanglesc
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
      Common/Multext/Ibord(0:Max_bord,2),Ixbordmax,Ncorners,
     .               Iarray(Max_part)
      Common/pawc/hmal(nplot)
c.......................................................................
      headmom(1,1)='     <x>    '
      headmom(1,2)='    <x''>    '
      headmom(1,3)='    <x^2>   '
      headmom(1,4)='    <x x''>  '
      headmom(1,5)='   <x''^2>   '
      headmom(1,6)='    <x^3>   '
      headmom(1,7)='  <x^2 x''>  '
      headmom(1,8)='  <x x''^2>  '
      headmom(1,9)='   <x''^3>   '
      headmom(1,10)='    <x^4>   '
      headmom(1,11)='  <x^3 x''>  '
      headmom(1,12)=' <x^2 x''^2> '
      headmom(1,13)='  <x x''^3>  '
      headmom(1,14)='   <x''^4>   '
c.......................................................................
      headmom(2,1)='     <y>    '
      headmom(2,2)='    <y''>    '
      headmom(2,3)='    <y^2>   '
      headmom(2,4)='    <y y''>  '
      headmom(2,5)='   <y''^2>   '
      headmom(2,6)='    <y^3>   '
      headmom(2,7)='  <y^2 y''>  '
      headmom(2,8)='  <y y''^2>  '
      headmom(2,9)='   <y''^3>   '
      headmom(2,10)='    <y^4>   '
      headmom(2,11)='  <y^3 y''>  '
      headmom(2,12)=' <y^2 y''^2> '
      headmom(2,13)='  <y y''^3>  '
      headmom(2,14)='   <y''^4>   '
c.......................................................................
      do ii=1,Max_bin
        do jj=1,Max_bin
          do kk=1,Max_pro
            Idist2D(ii,jj,kk)=0                        !..initialisation
          enddo
        enddo
      enddo
c.......................................................................
      if (iflag.eq.0) then        !..final beam distribution (core only)
c.......................................................................
        isw=1   !..used to switch on surface computation via convex hull
c.......................................................................
        do jdim=1,Idim/2                          !..loop over H/V plane
          Iplane=1+2*(jdim-1)
          Ic=2-jdim
c          Id1Dmax(jdim)=-1
c.......................................................................
          limitl(1,jdim)=(frame(1,Iplane)-frame(1,Iplane))*Stpi(Iplane)
     .                   +1
          limitu(1,jdim)=(frame(2,Iplane)-frame(1,Iplane))*Stpi(Iplane)
     .                   +1
          Npart(1)=Ntot
c.......................................................................
          limitl(2,jdim)=limitl(1,jdim)  !..index 2 used for final dist.
          limitu(2,jdim)=limitu(1,jdim)
          Npart(2)=0
c.......................................................................
          ext='dis'
          Call Gen_name(ext)                  !..file with distributions
          open(6+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          ext='opt'
          Call Gen_name(ext)                 !..file with optical param.
          open(7+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          write(7+Iplane,'(1x,a3,1x,a7,7(a13))') 'Isl#','  Nb   ',
     .      '    Emit     ','    Beta     ','    Alpha    ',
     .      ' Ellipse app.',' Convex hull', ' h param.  ', ' H param.  '
          ext='mom'
          Call Gen_name(ext)                        !..file with moments
          open(10+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          write(10+Iplane,'(1x,a3,1x,a7,14(a12))') 'Isl#','  Nb   ',
     .                                  (headmom(jdim,imom),imom=1,14)
c.......................................................................
          write(7+Iplane,'(1x,i3,1x,i7,7(1x,e12.6))') 
     .           0,Npart(1),emit(1,jdim),beta(1,jdim),alpha(1,jdim),
     .           ellipse(1,jdim),surf(1,jdim),halop1(1,jdim),
     .           halop2(1,jdim)
c.......................................................................
          write(10+Iplane,'(1x,i3,1x,i7,14(1x,e11.5))') 
     .           0,Npart(1),xmean(1,jdim),xpmean(1,jdim),
     .           x2(1,jdim),xxp(1,jdim),xp2(1,jdim),x3(1,jdim),
     .           x2xp(1,jdim),xxp2(1,jdim),xp3(1,jdim),x4(1,jdim),
     .           x3xp(1,jdim),x2xp2(1,jdim),xxp3(1,jdim),xp4(1,jdim)
c.......................................................................
          xmean(2,jdim)=0d0 
          xpmean(2,jdim)=0d0
          x2(2,jdim)=0d0
          xxp(2,jdim)=0d0
          xp2(2,jdim)=0d0
          x3(2,jdim)=0d0    !..higher-order moments used to compute halo
          x2xp(2,jdim)=0d0
          xxp2(2,jdim)=0d0 
          xp3(2,jdim)=0d0
          x4(2,jdim)=0d0
          x2xp2(2,jdim)=0d0
          x3xp(2,jdim)=0d0
          xxp3(2,jdim)=0d0
          xp4(2,jdim)=0d0
          do ii=1,Max_bin
            Idist1D(ii,2,jdim)=0
          enddo
c.......................................................................
          cminx=1d38                                   !..initialisation
          cminxp=1d38
          cmaxx=-1d38
          cmaxxp=-1d38
c.......................................................................
          do ipart=1,Ntot                         !..loop over particles
c.......................................................................
            if (Istab(ipart).eq.1) then
c.......................................................................
              Npart(2)=Npart(2)+1
              Ipoin(Npart(2))=ipart
c.......................................................................
              x2t=x(ipart,Iplane)*x(ipart,Iplane)
              x3t=x2t*x(ipart,Iplane)
              x4t=x2t*x2t
              xp2t=x(ipart,Iplane+1)*x(ipart,Iplane+1)
              xp3t=xp2t*x(ipart,Iplane+1)
              xp4t=xp2t*xp2t
              x2xpt=x2t*x(ipart,Iplane+1)
              xxp2t=x(ipart,Iplane)*xp2t
              x3xpt=x3t*x(ipart,Iplane+1)
              xxp3t=x(ipart,Iplane)*xp3t
              x2xp2t=x2t*xp2t
c.......................................................................
              xmean(2,jdim)=xmean(2,jdim)+x(ipart,Iplane)        !..mean
              xpmean(2,jdim)=xpmean(2,jdim)+x(ipart,Iplane+1)    !..mean
              x2(2,jdim)=x2(2,jdim)+x2t
              xxp(2,jdim)=xxp(2,jdim)+x(ipart,Iplane)*x(ipart,Iplane+1)
              xp2(2,jdim)=xp2(2,jdim)+xp2t
              x3(2,jdim)=x3(2,jdim)+x3t
              x4(2,jdim)=x4(2,jdim)+x4t
              xp3(2,jdim)=xp3(2,jdim)+xp3t
              xp4(2,jdim)=xp4(2,jdim)+xp4t
              x2xp(2,jdim)=x2xp(2,jdim)+x2xpt
              xxp2(2,jdim)=xxp2(2,jdim)+xxp2t
              xxp3(2,jdim)=xxp3(2,jdim)+xxp3t
              x3xp(2,jdim)=x3xp(2,jdim)+x3xpt
              x2xp2(2,jdim)=x2xp2(2,jdim)+x2xp2t
c.......................................................................
              do kdim=1,Idim         !..copies co-ordinates into xtransf
                xtransf(kdim)=x(ipart,kdim)
              enddo
              Call Transf                     !..transforms co-ordinates
              x1coor=xtransf(Iplane)
c.......................................................................
              ipointer=((x1coor-frame(1,Iplane))/Stepx(1))+1
              Idist1D(ipointer,2,jdim)=Idist1D(ipointer,2,jdim)+1
              Id1Dmax(jdim)=Max(Id1Dmax(jdim),Idist1D(ipointer,1,jdim))
c.......................................................................
              cminx=Dmin1(cminx,x(ipart,Iplane))
              cminxp=Dmin1(cminxp,x(ipart,Iplane+1))
              cmaxx=Dmax1(cmaxx,x(ipart,Iplane))
              cmaxxp=Dmax1(cmaxxp,x(ipart,Iplane+1))
c.......................................................................
            endif
c.......................................................................
          enddo
c.....................................................transforms moments
          avg(1)=Tmat(1,1,jdim)*xmean(2,jdim)+Tmat(1,2,jdim)*
     .                                                    xpmean(2,jdim)
          avg(2)=Tmat(2,1,jdim)*xmean(2,jdim)+Tmat(2,2,jdim)*
     .                                                    xpmean(2,jdim)
          xmean(2,jdim)=avg(1)/Dfloat(Npart(2))                  !..mean
          xpmean(2,jdim)=avg(2)/Dfloat(Npart(2))                 !..mean
c.......................................................................
          smom(1)=Tmom(1,1,jdim)*x2(2,jdim)+
     .            Tmom(1,2,jdim)*xxp(2,jdim)+Tmom(1,3,jdim)*xp2(2,jdim)
          smom(2)=Tmom(2,1,jdim)*x2(2,jdim)+
     .            Tmom(2,2,jdim)*xxp(2,jdim)+Tmom(2,3,jdim)*xp2(2,jdim)
          smom(3)=Tmom(3,1,jdim)*x2(2,jdim)+
     .            Tmom(3,2,jdim)*xxp(2,jdim)+Tmom(3,3,jdim)*xp2(2,jdim)
c.................higher-order moments are not transformed back!!!!!!!!!
          tmpm3=x3(2,jdim)/Dfloat(Npart(2))
          tmpm4=x4(2,jdim)/Dfloat(Npart(2))
          x3(2,jdim)=tmpm3-3d0*x2(2,jdim)/Dfloat(Npart(2))*xmean(2,jdim)
     .                    +2d0*xmean(2,jdim)*xmean(2,jdim)*xmean(2,jdim)
          x4(2,jdim)=tmpm4+6d0*x2(2,jdim)/Dfloat(Npart(2))*xmean(2,jdim)
     .          *xmean(2,jdim)-4d0*tmpm3*xmean(2,jdim)-3d0*xmean(2,jdim)
     .                        *xmean(2,jdim)*xmean(2,jdim)*xmean(2,jdim)
          dtest=x2(2,jdim)/Dfloat(Npart(2))-xmean(2,jdim)*xmean(2,jdim)
          if (dtest.ne.0d0) halop1(2,jdim)=x4(2,jdim)/(dtest*dtest)-2d0
c.......................................................................
          tmpmp3=xp3(2,jdim)/Dfloat(Npart(2))
          tmpmp4=xp4(2,jdim)/Dfloat(Npart(2))
          tmpm2p1=x2xp(2,jdim)/Dfloat(Npart(2))
          tmpm1p2=xxp2(2,jdim)/Dfloat(Npart(2))
          tmpm1p3=xxp3(2,jdim)/Dfloat(Npart(2))
          tmpm3p1=x3xp(2,jdim)/Dfloat(Npart(2))
          tmpm2p2=x2xp2(2,jdim)/Dfloat(Npart(2))
          x2xp(2,jdim)=tmpm2p1-x2(2,jdim)/Dfloat(Ntot)*xpmean(2,jdim)-
     .               2d0*xxp(2,jdim)/Dfloat(Ntot)*xmean(2,jdim)+
     .               2d0*xmean(2,jdim)*xmean(2,jdim)*xpmean(2,jdim)
          xxp2(2,jdim)=tmpm1p2-xp2(2,jdim)/Dfloat(Ntot)*xmean(2,jdim)-
     .               2d0*xxp(2,jdim)/Dfloat(Ntot)*xpmean(2,jdim)+
     .               2d0*xpmean(2,jdim)*xpmean(2,jdim)*xmean(2,jdim)
          xp3(2,jdim)=tmpmp3-3d0*xp2(2,jdim)/Dfloat(Npart(2))*
     .                 xpmean(2,jdim)+2d0*xpmean(2,jdim)*xpmean(2,jdim)*
     .                 xpmean(2,jdim)
          xp4(2,jdim)=tmpmp4+6d0*xp2(2,jdim)/Dfloat(Npart(2))*
     .      xpmean(2,jdim)*xpmean(2,jdim)-4d0*tmpmp3*xpmean(2,jdim)-3d0*
     .      xpmean(2,jdim)*xpmean(2,jdim)*xpmean(2,jdim)*xpmean(2,jdim)
          xxp3(2,jdim)=tmpm1p3-xmean(2,jdim)*tmpmp3-3d0*xpmean(2,jdim)*
     .     tmpm1p2+3d0*xmean(2,jdim)*xpmean(2,jdim)*xp2(2,jdim)/
     .               Dfloat(Npart(2))+3d0*xpmean(2,jdim)*xpmean(2,jdim)*
     .                   xxp(2,jdim)/Dfloat(Npart(2))-3d0*xmean(2,jdim)*
     .                   xpmean(2,jdim)*xpmean(2,jdim)*xpmean(2,jdim)
          x3xp(2,jdim)=tmpm3p1-xpmean(2,jdim)*tmpm3-3d0*xmean(2,jdim)*
     .              tmpm2p1+3d0*xpmean(2,jdim)*xmean(2,jdim)*x2(2,jdim)/
     .     Dfloat(Npart(2))+3d0*xmean(2,jdim)*xmean(2,jdim)*xxp(2,jdim)/
     .                Dfloat(Npart(2))-3d0*xpmean(2,jdim)*xmean(2,jdim)*
     .                                       xmean(2,jdim)*xmean(2,jdim)
          x2xp2(2,jdim)=tmpm2p2-2d0*tmpm2p1*xpmean(2,jdim)-2d0*tmpm1p2*
     .     xmean(2,jdim)+x2(2,jdim)/Dfloat(Npart(2))*xpmean(2,jdim)*
     .     xpmean(2,jdim)+4d0*xmean(2,jdim)*xpmean(2,jdim)*xxp(2,jdim)/
     .     Dfloat(Npart(2))+xp2(2,jdim)/Dfloat(Npart(2))*xmean(2,jdim)*
     .     xmean(2,jdim)-3d0*xmean(2,jdim)*xmean(2,jdim)*xpmean(2,jdim)*
     .     xpmean(2,jdim)
c.......................................................................
          x2(2,jdim)=smom(1)/Dfloat(Npart(2))-
     .                                       xmean(2,jdim)*xmean(2,jdim)
          xxp(2,jdim)=smom(2)/Dfloat(Npart(2))-
     .                                      xmean(2,jdim)*xpmean(2,jdim)
          xp2(2,jdim)=smom(3)/Dfloat(Npart(2))-
     .                                     xpmean(2,jdim)*xpmean(2,jdim)
c.......................................................................
          etest=x2(2,jdim)*xp2(2,jdim)-xxp(2,jdim)*xxp(2,jdim)
c.......................................................................
          if (etest.gt.0d0) then
            emit(2,jdim)=Dsqrt(etest)
c.......................................................................
            beta(2,jdim)=x2(2,jdim)/emit(2,jdim)
            alpha(2,jdim)=xxp(2,jdim)/emit(2,jdim)
c.......................................................................
          else
            write(*,*) '*** Error: negative emittance.'
            emit(2,jdim)=-1d0
c.......................................................................
            beta(2,jdim)=0d0
            alpha(2,jdim)=0d0
c.......................................................................
          endif
c.......................................................................
          ftest=3d0*(x4(2,jdim)*xp4(2,jdim)+3d0*x2xp2(2,jdim)*
     .               x2xp2(2,jdim)-4d0*xxp3(2,jdim)*x3xp(2,jdim))
          if (ftest.gt.0d0) halop2(2,jdim)=dsqrt(ftest)/(2d0*etest)-2d0
c.......................................................................
          xaxe=0.5d0*(cmaxx-cminx)     !..Computes surface using ellipse
          xpaxe=0.5d0*(cmaxxp-cminxp)
          ellipse(2,jdim)=0.5d0*pi2*xaxe*xpaxe*Dett
          surf(2,jdim)=Convex_hull(isw,jdim,Npart(2),Ipoin)*Dett
c.......................................................................
          write(7+Iplane,'(1x,i3,1x,i7,7(1x,e12.6))') 
     .        0,Npart(2),emit(2,jdim),beta(2,jdim),alpha(2,jdim),
     .          ellipse(2,jdim),surf(2,jdim),halop1(2,jdim),
     .           halop2(2,jdim)
c.......................................................................
          write(10+Iplane,'(1x,i3,1x,i7,14(1x,e11.5))') 
     .        0,Npart(2),xmean(2,jdim),xpmean(2,jdim),
     .        x2(2,jdim),xxp(2,jdim),xp2(2,jdim),x3(2,jdim),
     .        x2xp(2,jdim),xxp2(2,jdim),xp3(2,jdim),x4(2,jdim),
     .        x3xp(2,jdim),x2xp2(2,jdim),xxp3(2,jdim),xp4(2,jdim)
c.......................................................................
          do ibin=1,Nbin+1                        !..distributions (i/o)
            write(6+Iplane,*) (Idist1D(ibin,jbox,jdim),jbox=1,2)
          enddo
          write(6+Iplane,*) (limitl(jbox,jdim),jbox=1,2)
          write(6+Iplane,*) (limitu(jbox,jdim),jbox=1,2)
c.......................................................................
          Id1Dmax(jdim)=100*(Idint(Id1Dmax(jdim)+
     .                          2d0*Dsqrt(Dfloat(Id1Dmax(jdim)))+1)/100)
c.......................................................................
          close(6+Iplane)
          close(7+Iplane)
          close(10+Iplane)
c.......................................................................
          iprojt=iproj                          !..stores value of iproj
          iproj=1                              !..changes value of iproj
          Call Init_graph(0)
          Call Opt_graph(0)
          Call Def_palette                     !..defines colour palette
c.......................................................................
          do jbox=1,2
c.......................................................................
            Call Histo1_book(Ianim,jbox,jdim)
c.......................................................................
            do ibx=1,Nbin+1
c.......................................................................
              xc=frame(1,Iplane)+Stp(Iplane)*(ibx+0.5)
              yc=Idist1D(ibx,jbox,jdim)
              Call hf1(Iden+jbox,xc,yc)
c.......................................................................
            enddo
c.......................................................................
            Call Histo1_plot(jbox,jdim)
c.......................................................................
          enddo          
c.......................................................................
          Call hbook1(Iden+3,' ',Nbin,real(frame(1,Iplane)),
     .                                real(frame(2,Iplane)),0.)
          Call htitle('Superposition of all profiles')
          Call hmaxim(Iden+3,real(Id1Dmax(jdim)))
          Call hplot(Iden+3,' ',' ',0)
          do jbox=1,2
            Call hplset('HCOL',real(jbox/6.*290.))     !..changes colour
            Call hplot(Iden+jbox,'S',' ',0)   !..superimposes histograms
            Call hdelet(Iden+jbox)                 !..deletes histograms
          enddo
          Call hdelet(Iden+3)
c.......................................................................
          iproj=iprojt                        !..restores value of iproj
          Call End_graph
          close(Iden)
c.......................................................................
          ext='ps '
          Call Gen_name(ext)
          Call system('mv -f fort.7 '//fnameo(4+jdim))
c.......................................................................
        enddo
c.......................................................................
      elseif (iflag.eq.1) then      !..final beam dist. with disc. parts
c.......................................................................
        isw=1   !..used to switch on surface computation via convex hull
c.......................................................................
        do jdim=1,Idim/2                          !..loop over H/V plane
          Iplane=1+2*(jdim-1)
          Ic=2-jdim
          Id1Dmax(jdim)=-1
c.......................................................................
          limitl(1,jdim)=(frame(1,Iplane)-frame(1,Iplane))*Stpi(Iplane)
     .                   +1
          limitu(1,jdim)=(frame(2,Iplane)-frame(1,Iplane))*Stpi(Iplane)
     .                   +1
          Npart(1)=Ntot
c.......................................................................
          ext='dis'
          Call Gen_name(ext)                  !..file with distributions
          open(6+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          ext='opt'
          Call Gen_name(ext)                 !..file with optical param.
          open(7+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          write(7+Iplane,'(1x,a3,1x,a7,7(a13))') 'Isl#','  Nb   ',
     .      '    Emit     ','    Beta     ','    Alpha    ',
     .      ' Ellipse app.',' Convex hull', ' h param.  ', ' H param.  '
          ext='mom'
          Call Gen_name(ext)                        !..file with moments
          open(10+Iplane,file=fnameo(Ic),form='formatted',
     .                                                 status='unknown')
          write(10+Iplane,'(1x,a3,1x,a7,14(a12))') 'Isl#','  Nb   ',
     .                                  (headmom(jdim,imom),imom=1,14)
c.......................................................................
          write(7+Iplane,'(1x,i3,1x,i7,7(1x,e12.6))') 
     .           0,Npart(1),
     .           emit(1,jdim),beta(1,jdim),alpha(1,jdim),
     .           ellipse(1,jdim),surf(1,jdim),halop1(1,jdim),
     .           halop2(1,jdim)
c.......................................................................
          write(10+Iplane,'(1x,i3,1x,i7,14(1x,e11.5))') 
     .           0,Npart(1),xmean(1,jdim),xpmean(1,jdim),x2(1,jdim),
     .           xxp(1,jdim),xp2(1,jdim),x3(1,jdim),x2xp(1,jdim),
     .           xxp2(1,jdim),xp3(1,jdim),x4(1,jdim),x3xp(1,jdim),
     .           x2xp2(1,jdim),xxp3(1,jdim),xp4(1,jdim)
c.......................................................................
          do jbox=1,Nboxes(1)
c.......................................................................
            bordc(1,jdim,jbox+1)=0.5d0*(bound(1,jbox+1)+
     .                                  bound(3,jbox+1))       !..centre
            bordc(2,jdim,jbox+1)=0.5d0*(bound(2,jbox+1)+
     .                                  bound(4,jbox+1))
c.......................................................................
            bordw(1,jdim,jbox+1)=0.5d0*(bound(3,jbox+1)- 
     .                                  bound(1,jbox+1))        !..width
            bordw(2,jdim,jbox+1)=0.5d0*(bound(4,jbox+1)-
     .                                  bound(2,jbox+1))
c.......................................................................
            limitl(jbox+1,jdim)=(Tmat(1,1,jdim)*bound(1,jbox+1)-
     .                         frame(1,Iplane))*Stpi(Iplane)+1
            limitu(jbox+1,jdim)=(Tmat(1,1,jdim)*bound(3,jbox+1)-
     .                         frame(1,Iplane))*Stpi(Iplane)+1
c.......................................................................
            Npart(jbox+1)=0                            !..initialisation
            xmean(jbox+1,jdim)=0d0 
            xpmean(jbox+1,jdim)=0d0
            x2(jbox+1,jdim)=0d0
            xxp(jbox+1,jdim)=0d0
            xp2(jbox+1,jdim)=0d0
            x3(jbox+1,jdim)=0d0          !..higher-order to compute halo
            x2xp(jbox+1,jdim)=0d0
            xxp2(jbox+1,jdim)=0d0 
            xp3(jbox+1,jdim)=0d0
            x4(jbox+1,jdim)=0d0
            x2xp2(jbox+1,jdim)=0d0
            x3xp(jbox+1,jdim)=0d0
            xxp3(jbox+1,jdim)=0d0
            xp4(jbox+1,jdim)=0d0
            do ii=1,Max_bin
              Idist1D(ii,jbox+1,jdim)=0
            enddo
c.......................................................................
            bord(1,jdim,jbox+1)=bordc(1,jdim,jbox+1)
     .                         -bordw(1,jdim,jbox+1)
            bord(2,jdim,jbox+1)=bordc(2,jdim,jbox+1)
     .                         -bordw(2,jdim,jbox+1)
            bord(3,jdim,jbox+1)=bordc(1,jdim,jbox+1)
     .                         +bordw(1,jdim,jbox+1)
            bord(4,jdim,jbox+1)=bordc(2,jdim,jbox+1)
     .                         +bordw(2,jdim,jbox+1)
c.......................................................................
            cminx=1d38                                 !..initialisation
            cminxp=1d38
            cmaxx=-1d38
            cmaxxp=-1d38
c.......................................................................
            do ipart=1,Ntot                       !..loop over particles
              Pcheck=(bord(1,jdim,jbox+1).le.x(ipart,Iplane)).and.
     .               (bord(3,jdim,jbox+1).gt.x(ipart,Iplane)).and.
     .               (bord(2,jdim,jbox+1).le.x(ipart,Iplane+1)).and.
     .               (bord(4,jdim,jbox+1).gt.x(ipart,Iplane+1))
              if (Pcheck) then
c.......................................................................
                if (Istab(ipart).eq.1) then
c.......................................................................
                  Npart(jbox+1)=Npart(jbox+1)+1
                  Ipoin(Npart(jbox+1))=ipart  !..defines part. positions
c.......................................................................
                  x2t=x(ipart,Iplane)*x(ipart,Iplane)
                  x3t=x2t*x(ipart,Iplane)
                  x4t=x2t*x2t
                  xp2t=x(ipart,Iplane+1)*x(ipart,Iplane+1)
                  xp3t=xp2t*x(ipart,Iplane+1)
                  xp4t=xp2t*xp2t
                  x2xpt=x2t*x(ipart,Iplane+1)
                  xxp2t=x(ipart,Iplane)*xp2t
                  x3xpt=x3t*x(ipart,Iplane+1)
                  xxp3t=x(ipart,Iplane)*xp3t
                  x2xp2t=x2t*xp2t
                  xmean(jbox+1,jdim)=xmean(jbox+1,jdim)+
     .                             x(ipart,Iplane)               !..mean
                  xpmean(jbox+1,jdim)=xpmean(jbox+1,jdim)+
     .                              x(ipart,Iplane+1)            !..mean
                  x2(jbox+1,jdim)=x2(jbox+1,jdim)+x2t
                  xxp(jbox+1,jdim)=xxp(jbox+1,jdim)+
     .                                 x(ipart,Iplane)*x(ipart,Iplane+1)
                  xp2(jbox+1,jdim)=xp2(jbox+1,jdim)+xp2t
                  x3(jbox+1,jdim)=x3(jbox+1,jdim)+x3t
                  x4(jbox+1,jdim)=x4(jbox+1,jdim)+x4t
                  xp3(jbox+1,jdim)=xp3(jbox+1,jdim)+xp3t
                  xp4(jbox+1,jdim)=xp4(jbox+1,jdim)+xp4t
                  x2xp(jbox+1,jdim)=x2xp(jbox+1,jdim)+x2xpt
                  xxp2(jbox+1,jdim)=xxp2(jbox+1,jdim)+xxp2t
                  xxp3(jbox+1,jdim)=xxp3(jbox+1,jdim)+xxp3t
                  x3xp(jbox+1,jdim)=x3xp(jbox+1,jdim)+x3xpt
                  x2xp2(jbox+1,jdim)=x2xp2(jbox+1,jdim)+x2xp2t
c.......................................................................
                  do kdim=1,Idim     !..copies co-ordinates into xtransf
                    xtransf(kdim)=x(ipart,kdim)
                  enddo
                  Call Transf                 !..transforms co-ordinates
                  x1coor=xtransf(Iplane)
c.......................................................................
                  ipointer=((x1coor-frame(1,Iplane))/Stepx(1))+1
                  Idist1D(ipointer,jbox+1,jdim)=
     .                                   Idist1D(ipointer,jbox+1,jdim)+1
                  Id1Dmax(jdim)=Max(Id1Dmax(jdim),
     .                                         Idist1D(ipointer,1,jdim))
c.......................................................................
                  cminx=Dmin1(cminx,x(ipart,Iplane))
                  cminxp=Dmin1(cminxp,x(ipart,Iplane+1))
                  cmaxx=Dmax1(cmaxx,x(ipart,Iplane))
                  cmaxxp=Dmax1(cmaxxp,x(ipart,Iplane+1))
c.......................................................................
                endif
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
            if (Npart(jbox+1).eq.0) then    !..checks for zero particles
              write(*,*) '*** Error: ',Npart(jbox+1),' particles in box'
     .                   ,Ibox(jbox)
              Idist1D(limitl(jbox+1,jdim),jbox+1,jdim)=1 !..fakes distr.
              goto 20
            endif
c.....................................................transforms moments
            avg(1)=Tmat(1,1,jdim)*xmean(jbox+1,jdim)+Tmat(1,2,jdim)*
     .             xpmean(jbox+1,jdim)
            avg(2)=Tmat(2,1,jdim)*xmean(jbox+1,jdim)+Tmat(2,2,jdim)*
     .             xpmean(jbox+1,jdim)
            xmean(jbox+1,jdim)=avg(1)/Dfloat(Npart(jbox+1))      !..mean
            xpmean(jbox+1,jdim)=avg(2)/Dfloat(Npart(jbox+1))     !..mean
c.......................................................................
            smom(1)=Tmom(1,1,jdim)*x2(jbox+1,jdim)+
     .              Tmom(1,2,jdim)*xxp(jbox+1,jdim)+
     .              Tmom(1,3,jdim)*xp2(jbox+1,jdim)
            smom(2)=Tmom(2,1,jdim)*x2(jbox+1,jdim)+
     .              Tmom(2,2,jdim)*xxp(jbox+1,jdim)+
     .              Tmom(2,3,jdim)*xp2(jbox+1,jdim)
            smom(3)=Tmom(3,1,jdim)*x2(jbox+1,jdim)+
     .              Tmom(3,2,jdim)*xxp(jbox+1,jdim)+
     .              Tmom(3,3,jdim)*xp2(jbox+1,jdim)
c.................higher-order moments are not transformed back!!!!!!!!!
            tmpm3=x3(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm4=x4(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            x3(jbox+1,jdim)=tmpm3-3d0*x2(jbox+1,jdim)/
     .                         Dfloat(Npart(jbox+1))*xmean(jbox+1,jdim)+
     .      2d0*xmean(jbox+1,jdim)*xmean(jbox+1,jdim)*xmean(jbox+1,jdim)
            x4(jbox+1,jdim)=tmpm4+6d0*x2(jbox+1,jdim)/
     .                         Dfloat(Npart(jbox+1))*xmean(jbox+1,jdim)*
     .              xmean(jbox+1,jdim)-4d0*tmpm3*xmean(jbox+1,jdim)-3d0*
     .         xmean(jbox+1,jdim)*xmean(jbox+1,jdim)*xmean(jbox+1,jdim)*
     .         xmean(jbox+1,jdim)
            dtest=x2(jbox+1,jdim)/Dfloat(Npart(jbox+1))-
     .                             xmean(jbox+1,jdim)*xmean(jbox+1,jdim)
            if (dtest.ne.0d0) halop1(jbox+1,jdim)=x4(jbox+1,jdim)/
     .                                                 (dtest*dtest)-2d0
c.......................................................................
            tmpmp3=xp3(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpmp4=xp4(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm2p1=x2xp(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm1p2=xxp2(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm1p3=xxp3(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm3p1=x3xp(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            tmpm2p2=x2xp2(jbox+1,jdim)/Dfloat(Npart(jbox+1))
            x2xp(jbox+1,jdim)=tmpm2p1-x2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xpmean(jbox+1,jdim)-2d0*
     .          xxp(jbox+1,jdim)/Dfloat(Npart(jbox+1))*
     .          xmean(jbox+1,jdim)+2d0*xmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)*xpmean(jbox+1,jdim)
            xxp2(jbox+1,jdim)=tmpm1p2-xp2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xmean(jbox+1,jdim)-2d0*
     .          xxp(jbox+1,jdim)/Dfloat(Npart(jbox+1))*
     .          xpmean(jbox+1,jdim)+2d0*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xmean(jbox+1,jdim)
            xp3(jbox+1,jdim)=tmpmp3-3d0*xp2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xpmean(jbox+1,jdim)+2d0*
     .          xpmean(jbox+1,jdim)*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)
            xp4(jbox+1,jdim)=tmpmp4+6d0*xp2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)-4d0*tmpmp3*xpmean(jbox+1,jdim)-3d0*
     .          xpmean(jbox+1,jdim)*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xpmean(jbox+1,jdim)
            xxp3(jbox+1,jdim)=tmpm1p3-xmean(jbox+1,jdim)*tmpmp3-3d0*
     .          xpmean(jbox+1,jdim)*tmpm1p2+3d0*xmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xp2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))+3d0*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xxp(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))-3d0*xmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)
            x3xp(jbox+1,jdim)=tmpm3p1-xpmean(jbox+1,jdim)*tmpm3-3d0*
     .          xmean(jbox+1,jdim)*tmpm2p1+3d0*xpmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)*x2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))+3d0*xmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)*xxp(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))-3d0*xpmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)*xmean(jbox+1,jdim)*xmean(jbox+1,jdim)
            x2xp2(jbox+1,jdim)=tmpm2p2-2d0*tmpm2p1*xpmean(jbox+1,jdim)-
     .          2d0*tmpm1p2*xmean(jbox+1,jdim)+x2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)+4d0*xmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)*xxp(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))+xp2(jbox+1,jdim)/
     .          Dfloat(Npart(jbox+1))*xmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)-3d0*xmean(jbox+1,jdim)*
     .          xmean(jbox+1,jdim)*xpmean(jbox+1,jdim)*
     .          xpmean(jbox+1,jdim)
c.......................................................................
            x2(jbox+1,jdim)=smom(1)/Dfloat(Npart(jbox+1))-
     .                    xmean(jbox+1,jdim)*xmean(jbox+1,jdim)
            xxp(jbox+1,jdim)=smom(2)/Dfloat(Npart(jbox+1))-
     .                     xmean(jbox+1,jdim)*xpmean(jbox+1,jdim)
            xp2(jbox+1,jdim)=smom(3)/Dfloat(Npart(jbox+1))-
     .                     xpmean(jbox+1,jdim)*xpmean(jbox+1,jdim)
c.......................................................................
            etest=x2(jbox+1,jdim)*xp2(jbox+1,jdim)-
     .                                 xxp(jbox+1,jdim)*xxp(jbox+1,jdim)
c.......................................................................
            if (etest.gt.0d0) then
              emit(jbox+1,jdim)=Dsqrt(etest)
c.......................................................................
              beta(jbox+1,jdim)=x2(jbox+1,jdim)/emit(jbox+1,jdim)
              alpha(jbox+1,jdim)=xxp(jbox+1,jdim)/emit(jbox+1,jdim)
c.......................................................................
            else
              write(*,*) '*** Error: negative emittance.'
              emit(jbox+1,jdim)=-1d0
c.......................................................................
              beta(jbox+1,jdim)=0d0
              alpha(jbox+1,jdim)=0d0
c.......................................................................
            endif
            ftest=3d0*(x4(jbox+1,jdim)*xp4(jbox+1,jdim)+3d0*
     .                 x2xp2(jbox+1,jdim)*x2xp2(jbox+1,jdim)-
     .                 4d0*xxp3(jbox+1,jdim)*x3xp(jbox+1,jdim))
            if (ftest.gt.0d0) halop2(jbox+1,jdim)=dsqrt(ftest)/
     .                                                   (2d0*etest)-2d0
c.......................................................................
            xaxe=0.5d0*(cmaxx-cminx)   !..Computes surface using ellipse
            xpaxe=0.5d0*(cmaxxp-cminxp)
            ellipse(jbox+1,jdim)=0.5d0*pi2*xaxe*xpaxe*Dett
            surf(jbox+1,jdim)=Convex_hull(isw,jdim,Npart(jbox+1),Ipoin)
     .                                                             *Dett
c.......................................................................
            write(7+Iplane,'(1x,i3,1x,i7,7(1x,e12.6))') 
     .      Ibox(jbox),Npart(jbox+1),
     .      emit(jbox+1,jdim),beta(jbox+1,jdim),
     .      alpha(jbox+1,jdim),ellipse(jbox+1,jdim),surf(jbox+1,jdim),
     .      halop1(jbox+1,jdim),halop2(jbox+1,jdim)
c.......................................................................
            write(10+Iplane,'(1x,i3,1x,i7,14(1x,e11.5))') 
     .      Ibox(jbox),Npart(jbox+1),
     .      xmean(jbox+1,jdim),xpmean(jbox+1,jdim),
     .      x2(jbox+1,jdim),xxp(jbox+1,jdim),xp2(jbox+1,jdim),
     .      x3(jbox+1,jdim),x2xp(jbox+1,jdim),xxp2(jbox+1,jdim),
     .      xp3(jbox+1,jdim),x4(jbox+1,jdim),x3xp(jbox+1,jdim),
     .      x2xp2(jbox+1,jdim),xxp3(jbox+1,jdim),xp4(jbox+1,jdim)
c.......................................................................
 20         continue
c.......................................................................
          enddo
c.......................................................................
          do ibin=1,Nbin+1                        !..distributions (i/o)
            write(6+Iplane,*) 
     .                    (Idist1D(ibin,jbox,jdim),jbox=1,Nboxes(1)+1)
          enddo
          write(6+Iplane,*) (limitl(jbox,jdim),jbox=1,Nboxes(1)+1)
          write(6+Iplane,*) (limitu(jbox,jdim),jbox=1,Nboxes(1)+1)
c.......................................................................
          Id1Dmax(jdim)=Idint(Id1Dmax(jdim)+
     .                               2d0*Dsqrt(Dfloat(Id1Dmax(jdim)))+1)
c.......................................................................
          close(6+Iplane)
          close(7+Iplane)
          close(10+Iplane)
c.......................................................................
          iprojt=iproj                          !..stores value of iproj
          iproj=1                              !..changes value of iproj
          Call Init_graph(0)
          Call Opt_graph(0)
          Call Def_palette                     !..defines colour palette
c.......................................................................
          do jbox=1,Nboxes(1)+1
c.......................................................................
            Call Histo1_book(Ianim,jbox,jdim)
c.......................................................................
            do ibx=1,Nbin+1
c.......................................................................
              xc=frame(1,Iplane)+Stp(Iplane)*(ibx+0.5)
              yc=Idist1D(ibx,jbox,jdim)
              Call hf1(Iden+jbox,xc,yc)
c.......................................................................
            enddo
c.......................................................................
            Call Histo1_plot(jbox,jdim)
c.......................................................................
          enddo          
c.......................................................................
          Call hbook1(Iden+Nboxes(1)+2,' ',Nbin,real(frame(1,Iplane)),
     .                                         real(frame(2,Iplane)),0.)
          Call htitle('Superposition of all profiles')
          Call hmaxim(Iden+Nboxes(1)+2,real(Id1Dmax(jdim)))
          Call hplot(Iden+Nboxes(1)+2,' ',' ',0)
          do jbox=1,Nboxes(1)+1
            Call hplset('HCOL',real(jbox/6.*290.))     !..changes colour
            Call hplot(Iden+jbox,'S',' ',0)   !..superimposes histograms
            Call hdelet(Iden+jbox)                 !..deletes histograms
          enddo
          Call hdelet(Iden+Nboxes(1)+2)
c.......................................................................
          iproj=iprojt                        !..restores value of iproj
          Call End_graph
          close(Iden)
c.......................................................................
          ext='ps '
          Call Gen_name(ext)
          Call system('mv -f fort.7 '//fnameo(4+jdim))
c.......................................................................
        enddo
          close(Iden)
c.......................................................................
      elseif (iflag.eq.2) then         !..histogram of beam distribution
c.......................................................................
        if (Ianim.ge.1) then       !..initialisation of graphics package
          Iout=Iout+1
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        command=script(0)(iscriptl(0):iscriptu(0))
     .          //'''2:0@1(0,-5cm)+1@1(0,1cm)'' -pa4 '
     .          //fnameo(Idim/2)(1:index(fnameo(Idim/2),' '))
     .          //'tmp.ps'
c.......................................................................
        Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')       !..sets zone def.
c.......................................................................
        if (iturn.le.Ninj+Icentre) then            !..no stability check
c.......................................................................
          Id2Dmax(iproj+1)=-1    !..variable for global max. bin content
c.......................................................................
          do jproj=1,iproj                      !..loop over projections
c.......................................................................
            Id2Dmax(jproj)=-1               !..finds maximum bin content
            Iplane=1+2*(jproj-1)
            do ipart=1,Ntot
              do kdim=1,Idim         !..copies co-ordinates into xtransf
                xtransf(kdim)=x(ipart,kdim)
              enddo
              Call Transf                     !..transforms co-ordinates
              x1coor=xtransf(icpr(Iplane))
              x2coor=xtransf(icpr(Iplane+1))
c....................................................starts computations
              ipx=(x1coor-frame(1,icpr(Iplane)))*Stpi(icpr(Iplane))+1
              ipy=(x2coor-frame(1,icpr(Iplane+1)))
     .                                        *Stpi(icpr(Iplane+1))+1
              Idist2D(ipx,ipy,jproj)=Idist2D(ipx,ipy,jproj)+1
              Id2Dmax(jproj)=Max(Id2Dmax(jproj),Idist2D(ipx,ipy,jproj))
c.......................................................................
            enddo
c.......................................................................
            Id2Dmax(iproj+1)=Max(Id2Dmax(iproj+1),Id2Dmax(jproj))
c.......................................................................
          enddo
c.......................................................................
          Id2Dmax(iproj+1)=Idint(Dfloat(Id2Dmax(iproj+1))+
     .               2d0*Dsqrt(Dfloat(Id2Dmax(iproj+1))))+1
c.......................................................................
        else                                          
c.......................................................................
          do jproj=1,iproj                      !..loop over projections
c.......................................................................
            do ibx=1,Nbin+1                            !..initialisation
              do iby=1,Nbin+1
                Idist2D(ibx,iby,jproj)=0
              enddo
            enddo
c.......................................................................
            Iplane=1+2*(jproj-1)
            do ipart=1,Ntot
              if (Istab(ipart).eq.1) then             !..stability check
                do kdim=1,Idim       !..copies co-ordinates into xtransf
                  xtransf(kdim)=x(ipart,kdim)
                enddo
                Call Transf                   !..transforms co-ordinates
                x1coor=xtransf(icpr(Iplane))
                x2coor=xtransf(icpr(Iplane+1))
c....................................................starts computations
                ipx=(x1coor-frame(1,icpr(Iplane)))*Stpi(icpr(Iplane))+1
                ipy=(x2coor-frame(1,icpr(Iplane+1)))
     .                                          *Stpi(icpr(Iplane+1))+1
                Idist2D(ipx,ipy,jproj)=Idist2D(ipx,ipy,jproj)+1
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
          enddo
c.......................................................................
        endif
c.......................................................................
        do jproj=1,iproj                        !..loop over projections
c.......................................................................
          Iplane=1+2*(jproj-1)
          Call Histo2_book(Ianim,iturn,jproj)
c.......................................................................
          do ibx=1,Nbin+1
c.......................................................................
            do iby=1,Nbin+1
c.......................................................................
              xc=frame(1,icpr(Iplane))+Stp(icpr(Iplane))*(ibx+0.5)
              yc=frame(1,icpr(Iplane+1))+Stp(icpr(Iplane+1))*(iby-0.5)
              zc=Idist2D(ibx,iby,jproj)
              Call hf2(Iden,xc,yc,zc)
c.......................................................................
            enddo
c.......................................................................
          enddo
c.......................................................................
          Call Histo2_plot(iturn,jproj)
c.......................................................................
        enddo
        Call hplzon(1,1,1,' ')                       !..resets zone def.
c.......................................................................
        if (Ianim.ge.1) then              !..terminates graphics package
          if (Ianim.ge.3) Call Tune_plot(iturn)
          Call End_graph
          close(Iden)
c.......................................................................
          if (Ianim.ge.3) then
            Call system(command)
            Call system('mv -f tmp.ps '//fnameo(Idim/2))
          endif
c.......................................................................
          Call system(script(Ianim)(iscriptl(Ianim):
     .               iscriptu(Ianim))//fnameo(Idim/2))
          Call system('rm -f '//fnameo(Idim/2))            !..deletes ps
        endif
c.......................................................................
      elseif (iflag.eq.3) then !..histo. of beam dist. + islands removal
c.......................................................................
        if (Ianim.ge.1) then       !..initialisation of graphics package
          Iout=Iout+1
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        command=script(0)(iscriptl(0):iscriptu(0))
     .          //'''2:0@1(0,-5cm)+1@1(0,1cm)'' -pa4 '
     .          //fnameo(Idim/2)(1:index(fnameo(Idim/2),' '))
     .          //'tmp.ps'
c.......................................................................
        Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')       !..sets zone def.
c.......................................................................
        if (iturn.eq.0) then                       !..no stability check
c.......................................................................
          Id2Dmax(iproj+1)=-1    !..variable for global max. bin content
c.......................................................................
          do jproj=1,iproj                      !..loop over projections
c.......................................................................
            Id2Dmax(jproj)=-1               !..finds maximum bin content
            Iplane=1+2*(jproj-1)
            do ipart=1,Ntot
              do kdim=1,Idim         !..copies co-ordinates into xtransf
                xtransf(kdim)=x(Iarray(ipart),kdim)
              enddo
              Call Transf                     !..transforms co-ordinates
              x1coor=xtransf(icpr(Iplane))
              x2coor=xtransf(icpr(Iplane+1))
c....................................................starts computations
              ipx=(x1coor-frame(1,icpr(Iplane)))*Stpi(icpr(Iplane))+1
              ipy=(x2coor-frame(1,icpr(Iplane+1)))
     .                                        *Stpi(icpr(Iplane+1))+1
              Idist2D(ipx,ipy,jproj)=Idist2D(ipx,ipy,jproj)+1
              Id2Dmax(jproj)=Max(Id2Dmax(jproj),Idist2D(ipx,ipy,jproj))
c.......................................................................
            enddo
c.......................................................................
            Id2Dmax(iproj+1)=Max(Id2Dmax(iproj+1),Id2Dmax(jproj))
c.......................................................................
          enddo
c.......................................................................
          Id2Dmax(iproj+1)=Idint(Dfloat(Id2Dmax(iproj+1))+
     .               2d0*Dsqrt(Dfloat(Id2Dmax(iproj+1))))+1
c.......................................................................
        else
c.......................................................................
          Iflag_proj=0                     !..flag for projection (x,xp)
c.......................................................................
          do jproj=1,iproj                      !..loop over projections
c.......................................................................
            do ibx=1,Nbin+1                            !..initialisation
              do iby=1,Nbin+1
                Idist2D(ibx,iby,jproj)=0
              enddo
            enddo
c.......................................................................
            Iplane=1+2*(jproj-1)
c.......................................................................
            if ((icpr(Iplane).eq.1).and.(icpr(Iplane+1).eq.2)) 
     .                                                  Iflag_proj=jproj 
            do ipart=1,Npmult
              if (Istab(Iarray(ipart)).eq.1) then     !..stability check
                do kdim=1,Idim       !..copies co-ordinates into xtransf
                  xtransf(kdim)=x(Iarray(ipart),kdim)
                enddo
                Call Transf                   !..transforms co-ordinates
                x1coor=xtransf(icpr(Iplane))
                x2coor=xtransf(icpr(Iplane+1))
c....................................................starts computations
                ipx=(x1coor-frame(1,icpr(Iplane)))*Stpi(icpr(Iplane))+1
                ipy=(x2coor-frame(1,icpr(Iplane+1)))
     .                                          *Stpi(icpr(Iplane+1))+1
                Idist2D(ipx,ipy,jproj)=Idist2D(ipx,ipy,jproj)+1
c.......................................................................
              endif
c.......................................................................
            enddo
c.......................................................................
          enddo
c.......................................................................
          if (Iflag_proj.eq.0) then !..projection (x,xp) is not done yet
            Iflag_proj=iproj+1
            do ipart=1,Npmult
              if (Istab(Iarray(ipart)).eq.1) then     !..stability check
                do kdim=1,Idim       !..copies co-ordinates into xtransf
                  xtransf(kdim)=x(Iarray(ipart),kdim)
                enddo
                Call Transf                   !..transforms co-ordinates
                x1coor=xtransf(1)
                x2coor=xtransf(2)
c....................................................starts computations
                ipx=(x1coor-frame(1,1))*Stpi(1)+1
                ipy=(x2coor-frame(1,2))*Stpi(2)+1
                Idist2D(ipx,ipy,iproj+1)=Idist2D(ipx,ipy,iproj+1)+1
c.......................................................................
              endif
c.......................................................................
            enddo
          endif
c.......................................................................
          write(*,*) '*** Starts border search'
          ipcx=(0d0-frame(1,1))*Stpi(1)+1      !..co-ordinates of centre
          ipcy=(0d0-frame(1,2))*Stpi(2)+1
          imaxx=(3*Sigman-frame(1,1))*Stpi(1)+1   !..estimate for search
          imaxy=(3*Sigman-frame(1,2))*Stpi(2)+1
          nmaxx=imaxx-ipcx
          nmaxy=imaxy-ipcy
          k12min=ipcy                !..values to find min y co-ordinate
          k22min=ipcy
c.......................................................................
          do ipx=0,nmaxx                   !..starts scan in x direction
c.......................................................................
            k1x=ipcx+ipx                           !..positive direction
c.......................................................................
            do ipy=0,nmaxy                 !..starts scan in y direction
c.......................................................................
              k1y=ipcy+ipy                         !..positive direction
c.......................................................................
              if (Idist2D(k1x,k1y,Iflag_proj).eq.0) then   !..zero-cross
                indk11=ipx
                itemp11(indk11)=k1y
                if (ipy.eq.0) goto 40
                goto 35
              endif
            enddo
 35         continue
          enddo
c.......................................................................
 40       do ipx=0,nmaxx                   !..starts scan in x direction
c.......................................................................
            k1x=ipcx+ipx                           !..positive direction
c.......................................................................
            do ipy=1,nmaxy                 !..starts scan in y direction
c.......................................................................
              k2y=ipcy-ipy                         !..negative direction
c.......................................................................
              if (Idist2D(k1x,k2y,Iflag_proj).eq.0) then   !..zero-cross
                indk12=ipx
                itemp12(indk12)=k2y
c.......................................................................
                if (k2y.lt.k12min) then          !..searches for minimum
                  isentinel12=indk12
                  k12min=k2y
                endif
c.......................................................................
                if (ipy.eq.1) goto 50
                goto 45
              endif
            enddo
 45         continue
          enddo
c.......................................................................
 50       do ipx=1,nmaxx                   !..starts scan in x direction
c.......................................................................
            k2x=ipcx-ipx                           !..negative direction
c.......................................................................
            do ipy=0,nmaxy                 !..starts scan in y direction
c.......................................................................
              k1y=ipcy+ipy                         !..positive direction
c.......................................................................
              if (Idist2D(k2x,k1y,Iflag_proj).eq.0) then   !..zero-cross
                indk21=ipx
                itemp21(indk21)=k1y
                if (ipy.eq.0) goto 60
                goto 55
              endif
            enddo
 55         continue
          enddo
c.......................................................................
 60       do ipx=1,nmaxx                   !..starts scan in x direction
c.......................................................................
            k2x=ipcx-ipx                           !..negative direction
c.......................................................................
            do ipy=1,nmaxy                 !..starts scan in y direction
c.......................................................................
              k2y=ipcy-ipy                         !..negative direction
c.......................................................................
              if (Idist2D(k2x,k2y,Iflag_proj).eq.0) then   !..zero-cross
                indk22=ipx
                itemp22(indk22)=k2y
c.......................................................................
                if (k2y.lt.k22min) then
                  isentinel22=indk22
                  k22min=k2y
                endif
c.......................................................................
                if (ipy.eq.1) goto 70
                goto 65
              endif
            enddo
 65         continue
          enddo
c.......................................................................
 70       Ncorners=0                     !..computes border of beam core
          if (k22min.le.k12min) then
            do i22=isentinel22,1,-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx-i22
              Ibord(Ncorners,2)=itemp22(i22)
            enddo
            do i12=0,indk12
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx+i12
              Ibord(Ncorners,2)=itemp12(i12)
            enddo
            do i11=indk11,0,-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx+i11
              Ibord(Ncorners,2)=itemp11(i11)
            enddo
            do i21=1,indk21
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx-i21
              Ibord(Ncorners,2)=itemp21(i21)
            enddo
            do i22=indk22,isentinel22+1,-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx-i22
              Ibord(Ncorners,2)=itemp22(i22)
            enddo
          else
            do i12=isentinel12,indk12
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx+i12
              Ibord(Ncorners,2)=itemp12(i12)
            enddo
            do i11=indk11,0,-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx+i11
              Ibord(Ncorners,2)=itemp11(i11)
            enddo
            do i21=1,indk21
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx-i21
              Ibord(Ncorners,2)=itemp21(i21)
            enddo
            do i22=indk22,1,-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx-i22
              Ibord(Ncorners,2)=itemp22(i22)
            enddo
            do i12=0,isentinel12-1
              Ncorners=Ncorners+1
              Ibord(Ncorners,1)=ipcx+i12
              Ibord(Ncorners,2)=itemp12(i12)
            enddo
          endif
c.......................................................................
          isw=2                        !..computes convex hull of border
          idummy=0
          dummy=Convex_hull(isw,idummy,idummy,Ipoin)
c.......................................................................
          Ibord(0,1)=Ibord(Ncorners,1)
          Ibord(0,2)=Ibord(Ncorners,2)
          Ibord(Ncorners+1,1)=Ibord(1,1)
          Ibord(Ncorners+1,2)=Ibord(1,2)
c.......................................................................
          Multpart=0
          xmu=0d0
          xsq=0d0
          do ipart=1,Npmult
            do kdim=1,Idim           !..copies co-ordinates into xtransf
              xtransf(kdim)=x(Iarray(ipart),kdim)
            enddo
            Call Transf                       !..transforms co-ordinates
            x1coor=xtransf(1)
            x2coor=xtransf(2)
c....................................................starts computations
            ipx=(x1coor-frame(1,1))*Stpi(1)+1
            ipy=(x2coor-frame(1,2))*Stpi(2)+1
            iptest(1)=ipx
            iptest(2)=ipy
            if (inside(iptest)) then
c.......................................................................
              Multpart=Multpart+1
              Iscratch(Multpart)=Iarray(ipart)
c.......................................................................
              xmu=xmu+x1coor
              xsq=xsq+x1coor*x1coor
c.......................................................................
            endif
          enddo
c.......................................................................
          if (Multpart.eq.0) then
            write(*,*) '*** Error: no particle left. Program stops.'
            stop
          endif
c.......................................................................
          Npmult=Multpart
          Sigman=Dsqrt(xsq/Dfloat(Npmult)-
     .                        (xmu/Dfloat(Npmult))*(xmu/Dfloat(Npmult)))
c.......................................................................
          do ipart=1,Npmult
            Iarray(ipart)=Iscratch(ipart)
          enddo
c.......................................................................
        endif
c.......................................................................
        do jproj=1,iproj                        !..loop over projections
c.......................................................................
          Iplane=1+2*(jproj-1)
          Call Histo2_book(Ianim,iturn,jproj)
c.......................................................................
          do ibx=1,Nbin+1
c.......................................................................
            do iby=1,Nbin+1
c.......................................................................
              xc=frame(1,icpr(Iplane))+Stp(icpr(Iplane))*(ibx+0.5)
              yc=frame(1,icpr(Iplane+1))+Stp(icpr(Iplane+1))*(iby-0.5)
              zc=Idist2D(ibx,iby,jproj)
              Call hf2(Iden,xc,yc,zc)
c.......................................................................
            enddo
c.......................................................................
          enddo
c.......................................................................
          Call Histo2_plot(iturn,jproj)
c.......................................................................
        enddo
        Call hplzon(1,1,1,' ')                       !..resets zone def.
c.......................................................................
        if (Ianim.ge.1) then              !..terminates graphics package
          if (Ianim.ge.3) Call Tune_plot(iturn)
          Call End_graph
          close(Iden)
c.......................................................................
          if (Ianim.ge.3) then
            Call system(command)
            Call system('mv -f tmp.ps '//fnameo(Idim/2))
          endif
c.......................................................................
          Call system(script(Ianim)(iscriptl(Ianim):
     .               iscriptu(Ianim))//fnameo(Idim/2))
          Call system('rm -f '//fnameo(Idim/2))            !..deletes ps
        endif
c.......................................................................
      elseif (iflag.eq.4) then         !..histogram of beam distribution
c.......................................................................
        if (Ianim.ge.1) then       !..initialisation of graphics package
          Iout=Iout+1
          ext='ps '
          Call Gen_name(ext)
          open(Iden,file=fnameo(Idim/2),form='formatted',
     .                                                 status='unknown')
          Call Init_graph(0)
          Call Opt_graph(Ianim)
        endif
c.......................................................................
        command=script(0)(iscriptl(0):iscriptu(0))
     .          //'''2:0@1(0,-5cm)+1@1(0,1cm)'' -pa4 '
     .          //fnameo(Idim/2)(1:index(fnameo(Idim/2),' '))
     .          //'tmp.ps'
c.......................................................................
        Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')       !..sets zone def.
c.......................................................................
        do jproj=1,iproj                        !..loop over projections
c.......................................................................
          do ibx=1,Nbin+1                              !..initialisation
            do iby=1,Nbin+1
              Idist2D(ibx,iby,jproj)=0
            enddo
          enddo
c.......................................................................
          Iplane=1+2*(jproj-1)
          do ipart=1,Npmult
            if (Istab(Iarray(ipart)).eq.1) then       !..stability check
              do kdim=1,Idim         !..copies co-ordinates into xtransf
                xtransf(kdim)=x(Iarray(ipart),kdim)
              enddo
              Call Transf                     !..transforms co-ordinates
              x1coor=xtransf(icpr(Iplane))
              x2coor=xtransf(icpr(Iplane+1))
c....................................................starts computations
              ipx=(x1coor-frame(1,icpr(Iplane)))*Stpi(icpr(Iplane))+1
              ipy=(x2coor-frame(1,icpr(Iplane+1)))
     .                                          *Stpi(icpr(Iplane+1))+1
              Idist2D(ipx,ipy,jproj)=Idist2D(ipx,ipy,jproj)+1
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
          Call Histo2_book(Ianim,iturn,jproj)
c.......................................................................
          do ibx=1,Nbin+1
c.......................................................................
            do iby=1,Nbin+1
c.......................................................................
              xc=frame(1,icpr(Iplane))+Stp(icpr(Iplane))*(ibx+0.5)
              yc=frame(1,icpr(Iplane+1))+Stp(icpr(Iplane+1))*(iby-0.5)
              zc=Idist2D(ibx,iby,jproj)
              Call hf2(Iden,xc,yc,zc)
c.......................................................................
            enddo
c.......................................................................
          enddo
c.......................................................................
          Call Histo2_plot(iturn,jproj)
c.......................................................................
        enddo
        Call hplzon(1,1,1,' ')                       !..resets zone def.
c.......................................................................
        if (Ianim.ge.1) then              !..terminates graphics package
          if (Ianim.ge.3) Call Tune_plot(iturn)
          Call End_graph
          close(Iden)
c.......................................................................
          if (Ianim.ge.3) then
            Call system(command)
            Call system('mv -f tmp.ps '//fnameo(Idim/2))
          endif
c.......................................................................
          Call system(script(Ianim)(iscriptl(Ianim):
     .               iscriptu(Ianim))//fnameo(Idim/2))
          Call system('rm -f '//fnameo(Idim/2))            !..deletes ps
        endif
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Transf
c.......................................................................
c.... Subroutine to transform co-ordinates from special C-S to standard
c.... C-S to physical co-ordinates.
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Character*29 script
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
c.......................................................................
      Dimension temp(Max_dim)
c.......................................................................
      do irow=1,2                                            !..row scan
c.......................................................................
        temp(irow)=0d0                                 !..initialisation
        temp(irow+2)=0d0                               !..initialisation
        do icol=1,2                                       !..column scan
c.......................................................................
          temp(irow)=temp(irow)+Tmat(irow,icol,1)*xtransf(icol)
          temp(irow+2)=temp(irow+2)+Tmat(irow,icol,2)*xtransf(icol+2)
c.......................................................................
        enddo
c.......................................................................
      enddo
c.......................................................................
      do irow=1,Max_dim
        xtransf(irow)=temp(irow)
      enddo
c.......................................................................
      return
c.......................................................................
      end

      Double Precision Function Convex_hull(isw,jdim,Np,Ipoin)
c.......................................................................
c.... Function to compute the surface (emittance) of each part using a
c.... convex hull algorithm (see R. Sedgewick "Algorithms" Second 
c.... Edition
c....
c.... isw      : switch used to compute the surface of a distribution of
c....            points (iswitch = 1), or to find the convex hull of the 
c....            beam core (iswitch = 2) for multi-extractions.
c.... jdim  = 1: Horizontal plane (iswitch = 1 only).
c.... jdim  = 2: Vertical plane (iswitch = 1 only).
c.... Np       : is the number of particles (isw = 1 only)
c....
c.... Author: M. Giovannozzi
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Dimension temp(2),scratch(Max_part,2)
      Dimension Ipoin(Max_part)
c.......................................................................
      Character*29 script
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
      Common/Multext/Ibord(0:Max_bord,2),Ixbordmax,Ncorners,
     .               Iarray(Max_part)
c.......................................................................
      if (isw.eq.1) then
c.......................................................................
        Iplane=1+2*(jdim-1)
        Imin=1
        xave=0d0
        xpave=0d0
        Fnorm=1d0/Dfloat(Np)
        do ipart=1,Np
          scratch(ipart,1)=x(Ipoin(ipart),Iplane)
          scratch(ipart,2)=x(Ipoin(ipart),Iplane+1)
          xave=xave+scratch(ipart,1)
          xpave=xpave+scratch(ipart,2)
        enddo
        xave=xave*Fnorm
        xpave=xpave*Fnorm
c.......................................................................
        do ipart=2,Np
c.......................................................................
          if (scratch(ipart,2).lt.scratch(Imin,2)) Imin=ipart
        enddo
        Mindex=0
        do icoord=1,2
          scratch(Np+1,icoord)=scratch(Imin,icoord)
        enddo
        A_min=0d0
c.......................................................................
 10     Mindex=Mindex+1
        do icoord=1,2
          temp(icoord)=scratch(Mindex,icoord)
          scratch(Mindex,icoord)=scratch(Imin,icoord)
          scratch(Imin,icoord)=temp(icoord)
        enddo
        Imin=Np+1
        Angle=A_min
        A_min=360d0
c.......................................................................
        do ipart=Mindex+1,Np+1
c.......................................................................
          Tangle=Theta(scratch(Mindex,1),scratch(Mindex,2),
     .                 scratch(ipart,1),scratch(ipart,2))
          if (Tangle.gt.Angle) then
            if (Tangle.lt.A_min) then
              Imin=ipart
              A_min=Tangle
            endif
          endif
        enddo
c.......................................................................
        diffx=Dabs(scratch(Imin,1)-scratch(Np+1,1))
        diffy=Dabs(scratch(Imin,2)-scratch(Np+1,2))
        if (diffx+diffy.gt.tol) goto 10
c.......................................................................
        surf=0d0                              !..starts area computation 
        scratch(Mindex+1,1)=scratch(1,1)          !..adds starting point
        scratch(Mindex+1,2)=scratch(1,2)
        do itri=1,Mindex
c.......................................................................
          aside=Dsqrt((scratch(itri+1,1)-scratch(itri,1))*
     .                (scratch(itri+1,1)-scratch(itri,1))+
     .                (scratch(itri+1,2)-scratch(itri,2))*
     .                (scratch(itri+1,2)-scratch(itri,2)))     !..side a
          bside=Dsqrt((scratch(itri,1)-xave)*(scratch(itri,1)-xave)+
     .              (scratch(itri,2)-xpave)*(scratch(itri,2)-xpave))!..b
          cside=Dsqrt((scratch(itri+1,1)-xave)*(scratch(itri+1,1)-xave)+
     .              (scratch(itri+1,2)-xpave)*(scratch(itri+1,2)-xpave))
          per=aside+bside+cside                             !..perimeter
          area2=per*(per-2d0*aside)*(per-2d0*bside)*(per-2d0*cside)
c.......................................................................
          if (area2.lt.0d0) then
            write(*,*) '*** Error: negative area. Program stops.'
            stop
          endif
c.......................................................................
          area=0.25d0*Dsqrt(area2)
          surf=surf+area
c.......................................................................
        enddo
c.......................................................................
        Convex_hull=surf
c.......................................................................
      else
c.......................................................................
        do ipart=1,Ncorners
          scratch(ipart,1)=Ibord(ipart,1)
          scratch(ipart,2)=Ibord(ipart,2)
        enddo
c.......................................................................
        Imin=1
        Mindex=0
        do icoord=1,2
          scratch(Ncorners+1,icoord)=scratch(Imin,icoord)
        enddo
        A_min=0d0
c.......................................................................
 20     Mindex=Mindex+1
        do icoord=1,2
          temp(icoord)=scratch(Mindex,icoord)
          scratch(Mindex,icoord)=scratch(Imin,icoord)
          scratch(Imin,icoord)=temp(icoord)
        enddo
        Imin=Ncorners+1
        Angle=A_min
        A_min=360d0
c.......................................................................
        do ipart=Mindex+1,Ncorners+1
c.......................................................................
          Tangle=Theta(scratch(Mindex,1),scratch(Mindex,2),
     .                 scratch(ipart,1),scratch(ipart,2))
          if (Tangle.gt.Angle) then
            if (Tangle.lt.A_min) then
              Imin=ipart
              A_min=Tangle
            endif
          elseif ((Tangle.eq.Angle).and.((ipart.ne.Ncorners+1).or.
     .                                   (Mindex.ne.1))) then
            Imin=ipart
            A_min=Angle
            Mindex=Mindex-1+1/Mindex !..playing with integer arithmetics
            goto 20
          endif
        enddo
c.......................................................................
        if (Imin.lt.Ncorners+1) goto 20
c.......................................................................
        Ncorners=Mindex             !..redefinition of number of corners 
c.......................................................................
        do ipart=1,Ncorners
          Ibord(ipart,1)=scratch(ipart,1)
          Ibord(ipart,2)=scratch(ipart,2)
          write(2,*) Ibord(ipart,1),Ibord(ipart,2)
        enddo
c.......................................................................
        Convex_hull=1d0
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Double Precision Function Theta(x1,xp1,x2,xp2)
c.......................................................................
c.... Computes angle between vectors at (x1,xp1) and (x2,xp2). It is 
c.... used for the convex hull.
c....
c.... Author: M. Giovannozzi - CERN

c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Deltax=x2-x1
      Deltay=xp2-xp1
c.......................................................................
      adx=Dabs(Deltax)
      ady=Dabs(Deltay)
c.......................................................................
      if ((Deltax.eq.0d0).and.(Deltay.eq.0d0)) then
c.......................................................................
        tvar=0d0
c.......................................................................
      else
c.......................................................................
        tvar=Deltay/(adx+ady)
c.......................................................................
        if (Deltax.lt.0d0) then
c.......................................................................
          tvar=2d0-tvar
c.......................................................................
        elseif (Deltay.lt.0d0) then
c.......................................................................
          tvar=4d0+tvar
c.......................................................................
        endif
c.......................................................................
      endif
c.......................................................................
      Theta=90d0*tvar
c.......................................................................
      return
c.......................................................................
      end

      Integer Function iccw(ip0,ip1,ip2)
c.......................................................................
c.... Function used to determine the intersection between two lines.
c.... ip0,ip1,ip2 are three integer-valued vectors (see the book by 
c.... R. Sedgewick "Algorithms" for more details).
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Integer ip0(2),ip1(2),ip2(2)
c.......................................................................
      idx1=ip1(1)-ip0(1)
      idx2=ip2(1)-ip0(1)
      idy1=ip1(2)-ip0(2)
      idy2=ip2(2)-ip0(2)
c.......................................................................
      if (idx1*idy2.gt.idy1*idx2) iccw=1
      if (idx1*idy2.lt.idy1*idx2) iccw=-1
c.......................................................................
      if (idx1*idy2.eq.idy1*idx2) then
c.......................................................................
        if ((idx1*idx2.lt.0).or.(idy1*idy2.lt.0)) then
          iccw=-1
        elseif ((idx1*idx1+idy1*idy1).ge.(idx2*idx2+idy2*idy2)) then
          iccw=0
        else
          iccw=1
        endif
c.......................................................................
      endif
c.......................................................................
      return
c.......................................................................
      end

      Logical Function intersect(ip0,ip1,ip2,ip3)
c.......................................................................
c.... This function determines whether two lines have a real 
c.... intersection (see the book by R. Sedgewick "Algorithms" for more 
c.... details and the web site 
c.http://condor.informatik.uni-oldenburg.de/~stueker/graphic/index.html)
c....
c.... ip0, ip1 : define the first line.
c.... ip2, ip3 : define the second line.
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Integer ip0(2),ip1(2),ip2(2),ip3(2)
c.......................................................................
      itest1=iccw(ip0,ip1,ip2)
      itest2=iccw(ip0,ip1,ip3)
      itest3=iccw(ip2,ip3,ip0)
      itest4=iccw(ip2,ip3,ip1)
c.......................................................................
      intersect=((itest1*itest2.lt.0).and.(itest3*itest4.lt.0)).or.
     .          ((itest1*itest2*itest3*itest4.eq.0))
c.......................................................................
      return
c.......................................................................
      end

      Logical Function inside(ipoint)
c.......................................................................
c.... This function determines whether a point is inside a polygon (see 
c.... the book by R. Sedgewick "Algorithms" for more details and the web 
c.... site 
c.http://condor.informatik.uni-oldenburg.de/~stueker/graphic/index.html)
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Integer ipoint(2),lt1(2),lt2(2),lp1(2),lp2(2),lv1(2),lv2(2),ip(2)
c.......................................................................
      Logical intersect
c.......................................................................
      Common/Multext/Ibord(0:Max_bord,2),Ixbordmax,Ncorners,
     .               Iarray(Max_part)
c.......................................................................
      icount=0
      j=0
c.......................................................................
      lt1(1)=ipoint(1)
      lt1(2)=ipoint(2)
      lt2(1)=Ixbordmax
      lt2(2)=ipoint(2)
c.......................................................................
      lv1(1)=ipoint(1)
      lv1(2)=ipoint(2)
      lv2(1)=ipoint(1)
      lv2(2)=ipoint(2)
c.......................................................................
      do i=1,Ncorners
        lp1(1)=Ibord(i,1)
        lp1(2)=Ibord(i,2)
c.......................................................................
        lp2(1)=Ibord(i-1,1)
        lp2(2)=Ibord(i-1,2)
c.......................................................................
        if (intersect(lv1,lv2,lp1,lp2)) then
          inside=.true.
          return
        endif
        lp2(1)=Ibord(i,1)
        lp2(2)=Ibord(i,2)
        if (.not.intersect(lp1,lp2,lt1,lt2)) then
          lp2(1)=Ibord(j,1)
          lp2(2)=Ibord(j,2)
          j=i
          if (intersect(lp1,lp2,lt1,lt2)) then
            icount=icount+1
          else
            ip(1)=Ibord(j,1)
            ip(2)=Ibord(j,2)
            iccw1=iccw(lt1,lt2,ip)
            ip(1)=Ibord(i,1)
            ip(2)=Ibord(i,2)
            iccw2=iccw(lt1,lt2,ip)
            if ((i.ne.j+1).and.(iccw1*iccw2.lt.1)) then
              icount=icount+1
            endif
          endif
        endif
c.......................................................................
      enddo
      ip(1)=Ibord(j,1)
      ip(2)=Ibord(j,2)
      iccw1=iccw(lt1,lt2,ip)
      ip(1)=Ibord(1,1)
      ip(2)=Ibord(1,2)
      iccw2=iccw(lt1,lt2,ip)
      if ((j.ne.Ncorners).and.(iccw1*iccw2.eq.1)) icount=icount-1         
c.......................................................................
      inside=(Mod(icount,2).eq.1)
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Gen_name(ext)
c.......................................................................
c.... Subroutine to merge file name in variable fnamei with plane name 
c.... (H/V) and extension.
c.... ext: string with extension
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*3 ext
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo,tempch(0:Max_pro)
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      do i=0,Max_pro
        fnameo(i)=' '
        tempch(i)=' '
        tempch(i)(1:30)=fnamei
      enddo
      ipos=index(tempch(1),' ')                          !..finds length
c.......................................................................
      write(tempch(3)(ipos:ipos+7),'(a5,a3)') 'hpar.',ext
      write(tempch(4)(ipos:ipos+8),'(a6,a3)') 'hvpar.',ext
      write(tempch(5)(ipos:ipos+7),'(a5,a3)') 'hpro.',ext
      write(tempch(6)(ipos:ipos+7),'(a5,a3)') 'vpro.',ext
      write(tempch(7)(ipos:ipos+8),'(a6,a3)') 'summa.',ext
      write(tempch(8)(ipos:ipos+8),'(a4,a3)') 'par.',ext
      write(tempch(9)(ipos:ipos+7),'(a5,a3)') 'dump.',ext
c.......................................................................
      if (Iout.lt.10) then
        write(tempch(0)(ipos:ipos+5),'(a1,i1,a1,a3)') 
     .                                      'v',Iout,'.',ext
        write(tempch(1)(ipos:ipos+5),'(a1,i1,a1,a3)') 
     .                                      'h',Iout,'.',ext
        write(tempch(2)(ipos:ipos+6),'(a2,i1,a1,a3)') 
     .                                      'hv',Iout,'.',ext
      elseif (Iout.lt.100) then
        write(tempch(0)(ipos:ipos+6),'(a1,i2,a1,a3)') 
     .                                      'v',Iout,'.',ext
        write(tempch(1)(ipos:ipos+6),'(a1,i2,a1,a3)') 
     .                                      'h',Iout,'.',ext
        write(tempch(2)(ipos:ipos+7),'(a2,i2,a1,a3)') 
     .                                      'hv',Iout,'.',ext
      elseif (Iout.lt.1000) then
        write(tempch(0)(ipos:ipos+7),'(a1,i3,a1,a3)') 
     .                                      'v',Iout,'.',ext
        write(tempch(1)(ipos:ipos+7),'(a1,i3,a1,a3)') 
     .                                      'h',Iout,'.',ext
        write(tempch(2)(ipos:ipos+8),'(a2,i3,a1,a3)') 
     .                                      'hv',Iout,'.',ext
      endif
c.......................................................................
      do i=0,Max_pro
        ipos=index(tempch(i),' ')-1
        fnameo(i)(1:ipos)=tempch(i)(1:ipos)
      enddo
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Init_graph(Ifile)
c.......................................................................
c.... Subroutine to initialise the graphics packages: HBOOK, HPLOT, 
c.... HIGZ.
c.... Ifile is a flag used to change the metafile unit
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      do iel=1,nplot                                   !..initialisation
        hmal(iel)=0.0
      enddo
c.......................................................................
      kwtype=0
      Call hlimit(nplot)                        !..limit for common pawc
      Call hplint(kwtype)                            !..workstation type
c.......................................................................
      Call igmeta(-Iden-Ifile,-111)                  !..defines metafile
c.......................................................................
      Return
c.......................................................................
      End


      Subroutine Opt_graph(isw)
c.......................................................................
c.... Subroutine to set the graphics options for the packages: HBOOK, 
c.... HPLOT, HIGZ.
c.... isw is a flag used to print the date and to open multiple files:
c.... isw = 0: date is printed
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      Call hplopt('NPTO',1)                                   !..options
      if (isw.eq.0) then
        Call hplopt('DATE',1)
      endif
      Call hplopt('NBOX',1)
c.......................................................................
      Call hplset('XSIZ',xgsiz(iproj))              !..resets parameters
      Call hplset('YSIZ',ygsiz(iproj))
c.......................................................................
      Call hplset('CFON',-130.)
      Call hplset('GFON',-130.)
      Call hplset('LFON',-130.)
      Call hplset('TFON',-130.)
      Call hplset('VFON',-130.)
c.......................................................................
      Call hplset('DATE',3.00)
      Call hplset('CSIZ',0.50*scft(iproj))
      Call hplset('ASIZ',0.60*scft(iproj))
      Call hplset('NDVX',-205.01)
      Call hplset('NDVY',-205.13)
      Call hplset('XVAL',0.20)
      Call hplset('XLAB',2.20)
      Call hplset('XMGL',3.00)
      Call hplset('XMGR',3.00)
      Call hplset('XWIN',2.00*scfw(iproj))
      Call hplset('VSIZ',0.50*scft(iproj))
      Call hplset('GSIZ',0.60*scft(iproj))
      Call hplset('YGTI',1.10)
      Call hplset('YVAL',0.20)
      Call hplset('YLAB',1.30)
      Call hplset('NCOL',255.)
      Call igset('LTYP',1.00)
      Call igset('LWID',2.00)
      Call igset('MSCF',2.00)
      Call igset('TXAL',0.00)
      Call igset('CHHE',0.05*real(bbox(4)-bbox(2)))
      Call igset('TXFP',-130.)
c.......................................................................
      xoff=3d0          !..defines normalisation transformation (HPLSOF)
      yoff=2d0
      scx=xgsiz(iproj)-6d0
      scy=ygsiz(iproj)-4d0
c.......................................................................
      Return
c.......................................................................
      End


      Subroutine Histo1_book(isw,iturn,jdim)
c.......................................................................
c.... Subroutine to book 1D histograms.
c.... isw is a flag used to print the title:
c.... isw = 0: Title is printed (for printout only)
c.... iturn is the turn number
c.... jdim is the plane (1 - H, 2 - V)
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
      Character*80 title
c.......................................................................
      Logical hexist
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      if (hexist(Iden+iturn)) Call hdelet(Iden+iturn)  !..deletes histo.
c.......................................................................
      if (isw.ne.1) then
        title='       Beam slice # '
        write(title(21:21),'(i1)') iturn-1
        Call htitle(title)
      endif
c.......................................................................
      Iplane=1+2*(jdim-1)
      Call hbook1(Iden+iturn,' ',Nbin,real(frame(1,Iplane)),
     .                                real(frame(2,Iplane)),0.)
      Call hminim(Iden+iturn,0.)
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Histo2_book(isw,iturn,jproj)
c.......................................................................
c.... Subroutine to book 2D histograms.
c.... isw is a flag used to print the title:
c.... isw = 0: Title is printed (for printout only)
c.... iturn is the turn number
c.... jproj is the number identifying the projection
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
      Character*80 title
c.......................................................................
      Logical hexist
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      if (hexist(Iden)) Call hdelet(Iden)  !..deletes existing histogram
c.......................................................................
      if ((isw.ne.1).and.(isw.ne.3)) then
        title='       Beam distribution after'
        write(title(32:42),'(i5,a6)') iturn,' turns'
        Call htitle(title)
      endif
c.......................................................................
      ipl=1+2*(jproj-1)
c.......................................................................
      Call hbook2(Iden,' ',Nbin,real(frame(1,icpr(ipl))),
     .                          real(frame(2,icpr(ipl))),
     .                     Nbin,real(frame(1,icpr(ipl+1))),
     .                          real(frame(2,icpr(ipl+1))),0.)
      Call hminim(Iden,0.)
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Def_palette
c.......................................................................
c.... Subroutine to define the colour palette
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal
c.......................................................................
      Common/pawc/hmal(nplot)
c.......................................................................
      Call setshd( -1,0.0,0.0,0.0)
      Call setshd(  8,1.0,1.0,1.0)
      Call setshd( 60,0.0,0.0,1.0)
      Call setshd(120,0.0,1.0,0.5)
      Call setshd(180,1.0,1.0,0.0)
      Call setshd(230,1.0,0.5,0.0)
      Call setshd(255,1.0,0.0,0.0)
      Call shade
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine setshd(idxi,ri,gi,bi)
c.......................................................................
c.... Subroutine to define colours for palette.
c.... From PAW FAQS
c....
c.......................................................................
      Parameter (nmx = 20)
c.......................................................................
      Common /shpt/ npt,idx(nmx),r(nmx),g(nmx),b(nmx)
c.......................................................................
      if (idxi.lt.0) then
        npt=0
        return
      endif
c.......................................................................
      npt=npt+1
c.......................................................................
      if (npt.gt.nmx) then
        write(*,*) 'Error: too many colours'
        return
      endif
c.......................................................................
      idx(npt)=idxi
      r(npt)=ri
      g(npt)=gi
      b(npt)=bi
c.......................................................................
      return
c.......................................................................
      end

      Subroutine shade
c.......................................................................
c.... Subroutine to interpolate colours for palette.
c.... From PAW FAQS
c....
c.......................................................................
      Parameter (nmx=20)
c.......................................................................
      Common /shpt/ npt,idx(nmx),r(nmx),g(nmx),b(nmx)
c.......................................................................
      if (npt.lt.2) then
        write(*,*) 'Error: at least two colours are needed'
        return
      endif
c.......................................................................
      do i=2,npt
        j=i-1
        i1=idx(j)
        i2=idx(i)
        r1=r(j)
        g1=g(j)
        b1=b(j)
        r2=r(i)
        g2=g(i)
        b2=b(i)
        n =i2-i1+1
c.......................................................................
        do ii=i1,i2
          scale=float(ii-i1)/(n-1)
          rs=(r2-r1)*scale+r1
          gs=(g2-g1)*scale+g1
          bs=(b2-b1)*scale+b1
          Call iscr(1,ii,rs,gs,bs)
        enddo
c.......................................................................
      enddo
c.......................................................................
      return
c.......................................................................
      end

      Subroutine Histo1_plot(iturn,jproj)
c.......................................................................
c.... Subroutine to plot 1D histograms.
c.... iturn is a control parameter. If 
c.... iturn > 1 then the histograms are plotted using a zoom
c.... jproj is the number identifying the projection
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      Iplane=1+2*(jproj-1)
c.......................................................................
      Call hmaxim(Iden+iturn,real(Id1Dmax(jproj))) !..defines max histo.
      if (iturn.ne.0) then
        Call hplot(Iden+iturn,' ',' ',0)    
      else
        Call hplzom(Iden+iturn,' ',limitl(iturn,jproj),
     .                             limitu(iturn,jproj))
      endif
      Call hplax(chtit(Iplane,Ico+1),'Particle distribution')  !..titles 
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Histo2_plot(iturn,jproj)
c.......................................................................
c.... Subroutine to plot 2D histograms.
c.... iturn is the turn number
c.... jproj is the number identifying the projection
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      ipl=1+2*(jproj-1)
      ieven=mod(jproj,2)
c.......................................................................
      Call hmaxim(Iden,real(Id2Dmax(iproj+1))) !..defines max histo val.
      Call Def_palette                         !..defines colour palette
      Call hplot(Iden,cmap(ieven),' ',0)  !..colour map displayed or not
      Call hplax(chtit(icpr(ipl),Ico+1),chtit(icpr(ipl+1),Ico+1)) !..tit
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Graph2_book(isw,iturn,jproj)
c.......................................................................
c.... Subroutine to prepare a plot.
c.... isw is a flag used to print the title:
c.... isw = 0: Title is printed (for printout only)
c.... iturn is the turn number
c.... jproj is the projection 
c.... 
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
      Character*80 title
c.......................................................................
      Logical hexist
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      if (hexist(Iden)) Call hdelet(Iden)  !..deletes existing histogram
c.......................................................................
      if ((isw.ne.1).and.(isw.ne.3)) then
        title='       Beam distribution after'
        write(title(32:42),'(i5,a6)') iturn,' turns'
        Call htitle(title)
      endif
c.......................................................................
      ipl=1+2*(jproj-1)
c.......................................................................
      Call hbook2(Iden,' ',1,real(frame(1,icpr(ipl))),
     .                       real(frame(2,icpr(ipl))),
     .                     1,real(frame(1,icpr(ipl+1))),
     .                       real(frame(2,icpr(ipl+1))),0.)
c.......................................................................
      Call hminim(Iden,0.)
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Graph_plot(jproj)
c.......................................................................
c.... Subroutine to plot the graph.
c.... jproj is the number identifying the projection
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      ipl=1+2*(jproj-1)
c.......................................................................
      Call hmaxim(Iden,1.)                    !..defines max histo value
      Call hplot(Iden,' ',' ',0)                           !..plot graph
      Call hplax(chtit(icpr(ipl),Ico+1),chtit(icpr(ipl+1),Ico+1)) !..tit
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine End_graph
c.......................................................................
c.... Subroutine to terminate the graphics packages: HBOOK, HPLOT, HIGZ
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal
c.......................................................................
      Common/pawc/hmal(nplot)
c.......................................................................
      Call igmeta(999,0)                  !..closes all active metafiles
      Call hplend                          !..ends the graphics packages
      Call igend
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Write_par
c.......................................................................
c.... Subroutine to write the input parameters and related quantities.
c.... It generates also a file with:
c.... - Survived/lost particles
c.... - tunes
c.... - nonlinear gradients
c.... as a function of turn number
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
      Real rxvec(Max_part),rxpvec(Max_part)
      Real ryvec(Max_part),rypvec(Max_part)
      Real Sexg,Octg,Tunx,Tuny,Laps1,Laps2,Laps3(11),
     .     Nsurv,Nlost,Nextr(11)
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Twiss/AlphaH(4),BetaH(4),AmuH(4),AlphaV(4),BetaV(4),AmuV(4)
     .            ,AmatSepOct(4,4),AmatOctSex(4,4),
     .             Tmat(2,2,2),Tmati(2,2,2),Dett,Detti,
     .             Tmom(3,3,2),Tmomi(3,3,2),Alambda,Alambdai,
     .             xtransf(Max_dim),emx,emy
      Common/Rot/csepoctH,ssepoctH,csepoctV,ssepoctV,AoffH,AoffV,
     .           coctsexH,soctsexH,coctsexV,soctsexV,FtuneH,FtuneV
      Common/Distribution/Radius,Sigma,Sigmac,Amean,Bmean,Radv,
     .                    Sigv,Sigvc,Ameav,Bmeav,Ameanc,Bmeanc,
     .                    Ameanvc,Bmeanvc,xsep,Sigman,Sigmah(Max_bord),
     .                    Akick(Max_bord),rxvec,rxpvec,
     .                    ryvec,rypvec,Namp,Nangles,Nampc,Nanglesc
      Common/Ripple/Amprq(Max_par),Freqrq(Max_par),Phaserq(Max_par),
     .              Ampro(Max_par),Freqro(Max_par),Phasero(Max_par),
     .              Nriplq,Nriplo
      Common/Coordinates/x(Max_part,Max_dim),Istab(Max_part)
      Common/Histogram/Stepx(Max_boxes),bound(4,Max_boxes),
     .                 frame(2,Max_dim),bord(4,Max_dim,Max_boxes),
     .                 Nbx,Nby,Nbin,Nboxes(Max_dim),Ibox(Max_boxes),
     .                 Id1Dmax(Max_dim),Id2Dmax(Max_pro2)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/Plot/Sexg(Max_dat),Octg(Max_dat),
     .            Tunx(Max_dat),Tuny(Max_dat),Laps1(Max_dat),
     .            Laps2(Max_dat),Nsurv(Max_dat),Nlost(Max_dat),
     .            Icount1,Icount2
c.......................................................................
      write(Iden+11,*) '***Summary of main input parameters ***'
      write(Iden+11,*) '***Turn number            : ',Nturn
      write(Iden+11,*) '***Print delay            : ',Nout
      write(Iden+11,*) '***Print offset           : ',Noff
      write(Iden+11,*) '***Phase space dim.       : ',Idim
      write(Iden+11,*) '***Switch for Extr.       : ',Iext
      if (Iext.eq.1) write(Iden+11,*) '***Septum position        : ',
     .                                                              xsep
c.......................................................................
      write(Iden+11,*) '***# of projections       : ',iproj
      write(Iden+11,*) '***Projections chosen     : ',
     .                                             (icpr(l),l=1,2*iproj)
c.......................................................................
      write(Iden+11,*) '***Initial hor. tune      : ',par(1)/pi2
      write(Iden+11,*) '***Intermediate hor. tune : ',par(2)/pi2
      write(Iden+11,*) '***Final hor. tune        : ',par(3)/pi2
      if (Idim.eq.4) then
        if (par(16)-par(14).ne.0d0) then
          write(Iden+11,*) '***Initial ver. tune      : ',par(14)/pi2
          write(Iden+11,*) '***Intermediate ver. tune : ',par(15)/pi2
          write(Iden+11,*) '***Final ver. tune        : ',par(16)/pi2
        else
          write(Iden+11,*) '***Vertical tune          : ',par(14)/pi2
        endif
      endif
      if (Icomp.ne.4) then
        write(Iden+11,*) '***Kappa factor           : ',par(4)
        write(Iden+11,*) '***Number of kicks        : ',ipar(8)
      endif
      if (Icomp.eq.5) then
        do imul=1,Mult
          write(Iden+11,*) '***Kick at extr. #     ',imul,' : ',
     .    Akick(imul)
        enddo
      endif
      if (Icomp.eq.4) then
        do imul=1,Mult
          write(Iden+11,*) '***Kappa at extr. #    ',imul,' : ',
     .    par(4)*Sigmah(Mult+1)/Sigmah(imul)
        enddo
      endif
      write(Iden+11,*) '***Twiss switch           : ',Itwiss
      if (Itwiss.eq.0) then
        if (Idim.eq.4) write(Iden+11,*) '***Beta ratio             : ',
     .                   par(11)
      else
        do ip=1,4
          if (Idim.eq.2) then
            write(Iden+11,*) '***Hor. beta       @ pos',ip,': ',
     .                                                         BetaH(ip)
            write(Iden+11,*) '***Hor. alpha      @ pos',ip,': ',
     .                                                        AlphaH(ip)
            write(Iden+11,*) '***Hor. phase adv. @ pos',ip,': ',AmuH(ip)
          else
            write(Iden+11,*) '***Hor. beta       @ pos',ip,': ',
     .                                                         BetaH(ip)
            write(Iden+11,*) '***Hor. alpha      @ pos',ip,': ',
     .                                                        AlphaH(ip)
            write(Iden+11,*) '***Hor. phase adv. @ pos',ip,': ',AmuH(ip)
            write(Iden+11,*) '***Ver. beta       @ pos',ip,': ',
     .                                                         BetaV(ip)
            write(Iden+11,*) '***Ver. alpha      @ pos',ip,': ',
     .                                                        AlphaV(ip)
            write(Iden+11,*) '***Ver. phase adv. @ pos',ip,': ',AmuV(ip)
          endif
        enddo
      endif
c.......................................................................
      write(Iden+11,*) '***Ramping time           : ',ipar(4)
      write(Iden+11,*) '***Time first tune-ramp   : ',ipar(1)-ipar(4)
      write(Iden+11,*) '***Exp. first tune-ramp   : ',par(7)
      write(Iden+11,*) '***Time intermediate tune : ',ipar(2)-ipar(1)
      write(Iden+11,*) '***Time second tune-ramp  : ',ipar(3)-ipar(2)
      write(Iden+11,*) '***Exp. second tune-ramp  : ',par(31)
      if (Icomp.eq.4) then
        write(Iden+11,*) '***Time back to init. tune: ',iramp
        if (iblowup.ne.0) 
     .      write(Iden+11,*) '***Time to apply blow-up  : ',iblowup
      endif
c.......................................................................
      write(Iden+11,*) '***Coordinate switch      : ',Ico
      if (Ico.ne.0) then
        write(Iden+11,*) '***Hor. phys. rms emit.   : ',emx
        write(Iden+11,*) '***Ver. phys. rms emit.   : ',emy
      endif
c.......................................................................
      if (Nriplq.ne.0) then
        write(Iden+11,*) '***Tune ripple specifications ***'
        do irip=1,Nriplq
          write(Iden+11,*) '***Amplitude              : ',Amprq(irip)
          write(Iden+11,*) '***Frequency              : ',Freqrq(irip)
          write(Iden+11,*) '***Phase                  : ',Phaserq(irip)
        enddo
      endif
c.......................................................................
      if (Nriplo.ne.0) then
        write(Iden+11,*) 
     .                 '***Nonlinear elements ripple specifications ***'
        do irip=1,Nriplo
          write(Iden+11,*) '***Amplitude              : ',Ampro(irip)
          write(Iden+11,*) '***Frequency              : ',Freqro(irip)
          write(Iden+11,*) '***Phase                  : ',Phasero(irip)
        enddo
      endif
c.......................................................................
      write(Iden+11,*) '***Computation switch     : ',Icomp
      if (Icomp.eq.4) then
        write(Iden+11,*) '***Number of mult. extr.  : ',Mult
        do imul=1,Mult
          write(Iden+11,*) '***Sigma at extr. #    ',imul,' : ',
     .                                                      Sigmah(imul)
        enddo
      endif
      if ((Icomp.eq.7).or.(Icomp.eq.8)) then
        write(Iden+11,*) '***Number of injections   : ',Ninj
        if (Icentre.ne.0) then
          if (Idist.eq.1) then
            write(Iden+11,*) '***Uniform distribution (centre beam) ***'
            write(Iden+11,*) '***Maximum hor. radius    : ',Sigmac
            write(Iden+11,*) '***Hor. mean              : ',Ameanc
            write(Iden+11,*) '***Ang. mean              : ',Bmeanc
            if (Idim.eq.4) then 
              write(Iden+11,*) '***Maximum ver. radius    : ',Sigvc
              write(Iden+11,*) '***Ver. mean              : ',Ameavc
              write(Iden+11,*) '***Ang. mean              : ',Bmeavc
            endif
          elseif (Idist.eq.2) then
            write(Iden+11,*) 
     .                      '***Gaussian distribution (centre beam) ***'
            write(Iden+11,*) '***Hor. mean              : ',Ameanc
            write(Iden+11,*) '***Ang. mean              : ',Bmeanc
            write(Iden+11,*) '***Hor. sigma             : ',Sigmac
            if (Idim.eq.4) then
              write(Iden+11,*) '***Ver. mean              : ',Ameavc
              write(Iden+11,*) '***Ang. mean              : ',Bmeavc
              write(Iden+11,*) '***Ver. sigma             : ',Sigvc
            endif
          endif
          write(Iden+11,*) '***Radial steps           : ',Nampc
          write(Iden+11,*) '***Angular steps          : ',Nanglesc
        endif
      endif
      if (Idist.eq.1) then
        write(Iden+11,*) '***Uniform distribution ***'
        write(Iden+11,*) '***Maximum hor. radius    : ',Radius
        write(Iden+11,*) '***Hor. mean              : ',Amean
        write(Iden+11,*) '***Ang. mean              : ',Bmean
        if (Idim.eq.4) then 
          write(Iden+11,*) '***Maximum ver. radius    : ',Radv
          write(Iden+11,*) '***Ver. mean              : ',Ameav
          write(Iden+11,*) '***Ang. mean              : ',Bmeav
        endif
        write(Iden+11,*) '***Radial steps           : ',Namp
        write(Iden+11,*) '***Angular steps          : ',Nangles
      elseif (Idist.eq.2) then
        write(Iden+11,*) '***Gaussian distribution ***'
        write(Iden+11,*) '***Hor. mean              : ',Amean
        write(Iden+11,*) '***Ang. mean              : ',Bmean
        write(Iden+11,*) '***Hor. sigma             : ',Sigma
        if (Idim.eq.4) then
          write(Iden+11,*) '***Ver. mean              : ',Ameav
          write(Iden+11,*) '***Ang. mean              : ',Bmeav
          write(Iden+11,*) '***Ver. sigma             : ',Sigv
        endif
      elseif (Idist.eq.3) then
        write(Iden+11,*) '***Hollow distribution ***'
        write(Iden+11,*) '***Hor. peak position     : ',Amean
        write(Iden+11,*) '***Ang. peak position     : ',Bmean
        write(Iden+11,*) '***Hor. sigma             : ',Sigma
        if (Idim.eq.4) then
          write(Iden+11,*) '***Ver. peak position     : ',Ameav
          write(Iden+11,*) '***Ang. peak position     : ',Bmeav
          write(Iden+11,*) '***Ver. sigma             : ',Sigv
        endif
      elseif (Idist.eq.4) then
        write(Iden+11,*) '***Hollow/Gaussian distribution ***'
        write(Iden+11,*) '***Hor. peak position     : ',Amean
        write(Iden+11,*) '***Ang. peak position     : ',Bmean
        write(Iden+11,*) '***Hor. sigma             : ',Sigma
        if (Idim.eq.4) then
          write(Iden+11,*) '***Ver. peak position     : ',Ameav
          write(Iden+11,*) '***Ang. peak position     : ',Bmeav
          write(Iden+11,*) '***Ver. sigma             : ',Sigv
        endif
      endif
      write(Iden+11,*) '***Radial steps           : ',Namp
      write(Iden+11,*) '***Angular steps          : ',Nangles
      write(Iden+11,*) '***Total init. conditions : ',Ntot
      if (Iext.eq.0) then
        Ilst=Nlost(Icount1)
        write(Iden+11,*) '***Lost particles         : ',Ilst
      endif
      if (Icomp.eq.8) then
        write(Iden+11,*) '***Exponent sext. curve   : ',par(36)
        write(Iden+11,*) '***Exponent oct. curve    : ',par(26)
      endif
c.......................................................................
      write(Iden+11,*) '***Bounding box           : ',
     .               '(',bbox(1),',',bbox(2),',',bbox(3),',',bbox(4),')'
      write(Iden+11,*) '***Grid size              : ',Nbx,' x',Nby
      write(Iden+11,*) '***Number of bins (histo) : ',Nbin
      write(Iden+11,*) '***Number of boxes        : ',Nboxes(1)
      write(Iden+11,*) '***Boxes id.              : ',
     .              (Ibox(jbox),jbox=1,Nboxes(1))
      write(Iden+11,*) '***Graphics format switch : ',Ianim
c.......................................................................
      if (Idim.eq.2) then !..writes numerical values of phys. parameters
c.......................................................................
        do idat=1,Icount1 
          write(Iden+12,*) Laps1(idat),Nsurv(idat),Nlost(idat),
     .                       Tunx(idat),Sexg(idat),Octg(idat)
        enddo
c.......................................................................
      else 
c.......................................................................
        do idat=1,Icount1 
          write(Iden+12,*) Laps1(idat),Nsurv(idat),Nlost(idat),
     .                       Tunx(idat),Tuny(idat),Sexg(idat),Octg(idat)
        enddo
c.......................................................................
      endif
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Plot_par(Ierr)
c.......................................................................
c.... Subroutine to plot the simulation parameters:
c.... 
c.... - intensity
c.... - nonlinear gradients
c.... - tunes
c.... - resonance net
c.... Ierr is an error flag
c....
c.... NB: some parameters are artificially modified to reuse some 
c....     routines without modifications.
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Dimension xlow(6),xhig(6),ylow(6),yhig(6),gcomp(4)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
      Real Sexg,Octg,Tunx,Tuny,Laps1,Laps2,Laps3(11),
     .     Nsurv,Nlost,Nextr(11)
      Real xskey(2),yskey(2),xokey(2),yokey(2)
      Real xskeyp,yskeyp,xokeyp,yokeyp
      Real Iextract(12),ext(12)
c.......................................................................
      Character*4 cmap
      Character*11 chtitx(6)
      Character*13 chtit
      Character*17 chtity(6)
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Logical hexist
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/Plot/Sexg(Max_dat),Octg(Max_dat),
     .            Tunx(Max_dat),Tuny(Max_dat),Laps1(Max_dat),
     .            Laps2(Max_dat),Nsurv(Max_dat),Nlost(Max_dat),
     .            Icount1,Icount2
      Common/pawc/hmal(nplot)
c.......................................................................
      Data Laps3/0.,1.,1.,2.,2.,3.,3.,4.,4.,5.,5./
c.......................................................................
      Ierr=0
      if (Icount2.gt.1000) then
        write(*,*) '*** Error: # of points exceeds HPLOT capability.'
        Ierr=1
        return
      endif
      if ((Icount1.lt.2).or.(Icount2.lt.2)) then
        write(*,*) '*** Error: # of points should be > 2.'
        Ierr=1
        return
      endif
c.......................................................................
      iprojt=iproj
      iproj=-2*(Idim-2)+7 !..warning: artificial redefinition of iproj!!
c.......................................................................
      if (Ntot.gt.50) then
        iresidue=Ntot/5           !..defines max int. for opt. labelling
        iexp=Dlog10(Dfloat(iresidue))
        amagn=10**iexp
        imant=iresidue/amagn
        Aintmax=5d0*amagn*(imant+1)
      else
        Aintmax=Dfloat(Ntot)
      endif
c.......................................................................
      Nturntot=Nturn+((Mult-1)*Iper+iramp+iblowup)*Icomp/4 !..mult. ext.
      if (Nturn.gt.50) then
        iresidue=(Nturntot-4*Iext)    !..defines max turn for opt. label
        iexp=Dlog10(Dfloat(iresidue))
        amagn=10**iexp
        imant=iresidue/amagn
        Antmax=amagn*(imant+0.5)
      else
        Antmax=Dfloat(Nturntot-4*Iext)
      endif
c.......................................................................
      amaxnux=Dmax1(par(1),par(3))
      aminnux=Dmin1(par(1),par(3))
      Dnux=(amaxnux-aminnux)/pi2                            !..hor. tune
      Dinux=1d0/Dnux
      iexp=Dlog10(Dinux)+1
      amagn=10**iexp
      imant=Dnux*amagn/5
      Dnux=5d0/amagn*(imant+1)
      centrex=Idint(5d-1*(aminnux+amaxnux)/pi2*amagn)/amagn
c.......................................................................
      if (par(16)-par(14).ne.0d0) then
        amaxnuy=Dmax1(par(14),par(16))
        aminnuy=Dmin1(par(14),par(16))
        Dnuy=(amaxnuy-aminnuy)/pi2                          !..ver. tune
        Dinuy=1d0/Dnuy
        iexp=Dlog10(Dinuy)+1
        amagn=10**iexp
        imant=Dnuy*amagn/5
        Dnuy=5d0/amagn*(imant+1)
        centrey=Idint(5d-1*(aminnuy+amaxnuy)/pi2*1d3)*1d-3
        tmin2=centre-Dnuy/2                                   
        tmax2=centre+Dnuy/2 
      else
        Dnuy=Dnux
        centrey=Idint(par(14)/pi2*1d3)*1d-3
      endif
      Dnu=Dmax1(Dnux,Dnuy)
      tmin1=centrex-Dnu/2                                   
      tmax1=centrex+Dnu/2 
      tmin2=centrey-Dnu/2                                   
      tmax2=centrey+Dnu/2 
c.......................................................................
      gcomp(1)=Sexg(1)
      gcomp(2)=Sexg(Icount2)
      gcomp(3)=Octg(1)
      gcomp(4)=Octg(Icount2)
      gmin=1e38
      gmax=-1e38
      do ii=1,4
        gmin=dmin1(gcomp(ii),gmin)
        gmax=dmax1(gcomp(ii),gmax)
      enddo
      gmin=gmin*(1.-0.1*gmin/Abs(gmin))
      gmax=gmax*(1.+0.1*gmax/Abs(gmax))
      Dgrad=gmax-gmin
      Dgrad=5d-1*(Idint(Dgrad/5*1d1)+1)
      centre=1d-1*Idint((gmin+gmax)/2*1d1)
      gmin=centre-Dgrad/2
      gmax=centre+Dgrad/2
c.......................................................................
      xskey(1)=8d-1*Antmax
      yskey(1)=5.5d-1*gmax
      xskey(2)=xskey(1)+Antmax/10
      yskey(2)=yskey(1)
      xokey(1)=8d-1*Antmax
      yokey(1)=3.d-1*gmax
      xokey(2)=xokey(1)+Antmax/10
      yokey(2)=yokey(1)
c.......................................................................
      do itur=1,5                     !..initialises extracted intensity
        Nextr(2*itur)=Nsurv(Icount1-5+itur)/Nsurv(Icount1-4)
        Nextr(2*itur-1)=Nextr(2*itur)
      enddo
      Nextr(11)=0.0
c.......................................................................
      if (Idim.eq.2) then        !..defines coordinates in cm for hplsof
c.......................................................................
        xskeyp=3+(xgsiz(iproj)-6)*(xskey(2)+Antmax/50)/Antmax
        yskeyp=(ygsiz(iproj)+4)/3+(ygsiz(iproj)-8)/3*
     .         (yskey(1)-5d-2*gmax-gmin)/Dgrad
c.......................................................................
        xokeyp=3+(xgsiz(iproj)-6)*(xokey(2)+Antmax/50)/Antmax
        yokeyp=(ygsiz(iproj)+4)/3+(ygsiz(iproj)-8)/3*
     .         (yokey(1)-5d-2*gmax-gmin)/Dgrad
c.......................................................................
      else
c.......................................................................
        xskeyp=(xgsiz(iproj)+4)/2+(xgsiz(iproj)-10)/2*
     .         (xskey(2)+Antmax/50)/Antmax
        yskeyp=(ygsiz(iproj)+2)/2+(ygsiz(iproj)-6)/2*
     .         (yskey(1)-5d-2*gmax-gmin)/Dgrad
c.......................................................................
        xokeyp=(xgsiz(iproj)+4)/2+(xgsiz(iproj)-10)/2*
     .         (xokey(2)+Antmax/50)/Antmax
        yokeyp=(ygsiz(iproj)+2)/2+(ygsiz(iproj)-6)/2*
     .         (yokey(1)-5d-2*gmax-gmin)/Dgrad
c.......................................................................
      endif
c.......................................................................
      xlow(1)=0d0                      !..graphical windows for plotting
      ylow(1)=0d0
      xhig(1)=Antmax
      yhig(1)=Aintmax   
      xlow(2)=0d0
      ylow(2)=0d0
      xhig(2)=6
      yhig(2)=1.2   
      xlow(3)=0d0
      ylow(3)=gmin
      xhig(3)=Antmax
      yhig(3)=gmax
      xlow(4)=0d0
      ylow(4)=tmin1
      xhig(4)=Antmax
      yhig(4)=tmax1   
      xlow(5)=0d0
      ylow(5)=tmin2
      xhig(5)=Antmax
      yhig(5)=tmax2   
c.......................................................................
      chtitx(1)='Turn number'
      chtity(1)='Particles        '
      chtitx(2)='Machine turn     '
      chtity(2)='Norm. particles  '
      chtitx(3)='Turn number'
      chtity(3)='Nonlinear grad.  '
      chtitx(4)='Turn number'
      chtity(4)='Hor. tune        '
      chtitx(5)='Turn number'
      chtity(5)='Ver. tune        '
      chtitx(6)='Hor. tune        '
      chtity(6)='Ver. tune        '
c.......................................................................
      Call Init_graph(10)
      Call Opt_graph(0)
c.......................................................................
      if (Iext.eq.0) then
c.......................................................................
        Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')        !..defines zones
c..............................................................Intensity
        if (hexist(Iden)) Call hdelet(Iden)   !..deletes existing histo.
        Call hbook2(Iden,' ',1,real(xlow(1)),real(xhig(1)),
     .                       1,real(ylow(1)),real(yhig(1)),0.)
        Call hminim(Iden,0.)
        Call hmaxim(Iden,1.)                    
c.......................................................................
        Call hplot(Iden,' ',' ',0)                         !..plot graph
        Call hplax(chtitx(1),chtity(1))                !..titles on axis
c.......................................................................
        Call hpline(Laps1,Nsurv,Icount1,' ')               !..plots line
c.......................................................................
      else
c.......................................................................
        if (Idim.eq.2) then                   !..splits one of the zones
          Call hplzon(2,3,1,' ')
        else
          Call hplzon(2,4,1,' ')
        endif
c.......................................................................
        if (hexist(Iden)) Call hdelet(Iden)   !..deletes existing histo.
        Call hbook2(Iden,' ',1,real(xlow(1)),real(xhig(1)),
     .                       1,real(ylow(1)),real(yhig(1)),0.)
        Call hminim(Iden,0.)
        Call hmaxim(Iden,1.)                    
c.......................................................................
        Call hplot(Iden,' ',' ',0)                         !..plot graph
        Call hplax(chtitx(1),chtity(1))                !..titles on axis
c.......................................................................
        Call hpline(Laps1,Nsurv,Icount1,' ')               !..plots line
c.......................................................................
        if (Idim.eq.4) Call hplzon(2,4,3,'S')           !..sets position
        Call hplset('NDVX',-206.01)                    !..defines labels
        Call hplset('NDVY',-206.13)
        if (hexist(Iden)) Call hdelet(Iden)   !..deletes existing histo.
        Call hbook2(Iden,' ',1,real(xlow(2)),real(xhig(2)),
     .                       1,real(ylow(2)),real(yhig(2)),0.)
        Call hminim(Iden,0.)
        Call hmaxim(Iden,1.)                    
c.......................................................................
        Call hplot(Iden,' ',' ',0)                         !..plot graph
        Call hplax(chtitx(2),chtity(2))                !..titles on axis
c.......................................................................
        Call hpline(Laps3,Nextr,11,' ')                    !..plots line
        Call hplset('NDVX',-205.01)          !..resets labels definition
        Call hplset('NDVY',-205.13)
c.......................................................................
      endif
c.......................................................................
      if (Iext.eq.1) Call hplzon(Nxz(iproj),Nyz(iproj),2,'S') !..resets
c....................................................Nonlinear gradients
      if (hexist(Iden)) Call hdelet(Iden) !..deletes existing histogram
c.......................................................................
      Call hbook2(Iden,' ',1,real(xlow(3)),real(xhig(3)),
     .                     1,real(ylow(3)),real(yhig(3)),0.)
      Call hminim(Iden,0.)
      Call hmaxim(Iden,1.)
c.......................................................................
      Call hplot(Iden,' ',' ',0)                           !..plot graph
      Call hplax(chtitx(3),chtity(3))                  !..titles on axis
c.......................................................................
      Call hpline(Laps2,Sexg,Icount2,' ')                  !..plots line
      Call hpline(xskey,yskey,2,' ')
      Call hplsof(xskeyp,yskeyp,'K?2!',0.30,0.,0.,-1)
c.......................................................................
      Call igset('LTYP',2.00)                       !..changes line type
      Call hpline(Laps2,Octg,Icount2,' ')                  !..plots line
      Call hpline(xokey,yokey,2,' ')
      Call hplsof(xokeyp,yokeyp,'K?3!',0.30,0.,0.,-1)
c..............................................................Hor. tune
      if (hexist(Iden)) Call hdelet(Iden) !..deletes existing histogram
c.......................................................................
      Call hbook2(Iden,' ',1,real(xlow(4)),real(xhig(4)),
     .                     1,real(ylow(4)),real(yhig(4)),0.)
      Call hminim(Iden,0.)
      Call hmaxim(Iden,1.)
c.......................................................................
      Call hplot(Iden,' ',' ',0)                           !..plot graph
      Call hplax(chtitx(4),chtity(4))                  !..titles on axis
c.......................................................................
      Call igset('LTYP',1.00)                       !..changes line type
      Call hpline(Laps2,Tunx,Icount2,' ')                  !..plots line
c.......................................................................
      if (Idim.eq.4) then
c..............................................................Ver. tune
        if (hexist(Iden)) Call hdelet(Iden)!..deletes existing histogram
c.......................................................................
        Call hbook2(Iden,' ',1,real(xlow(5)),real(xhig(5)),
     .                       1,real(ylow(5)),real(yhig(5)),0.)
        Call hminim(Iden,0.)
        Call hmaxim(Iden,1.)
c.......................................................................
        Call hplot(Iden,' ',' ',0)                         !..plot graph
        Call hplax(chtitx(5),chtity(5))                !..titles on axis
c.......................................................................
        Call hpline(Laps2,Tuny,Icount2,' ')                !..plots line
c...........................................................Tune diagram
        if (hexist(Iden)) Call hdelet(Iden)!..deletes existing histogram
c.......................................................................
        iproj=1          !..warning: artificial redefinition of iproj!!
        Call Opt_graph(0)
        Call hplzon(1,1,1,' ')                 !..resets zone definition
c.......................................................................
        xlow(6)=1d-2*Idint(tmin1*1d2)                       !..hor. tune
        xhig(6)=1d-2*Idint(tmax1*1d2+1d0)
        ylow(6)=1d-2*Idint(tmin2*1d2)                       !..ver. tune
        yhig(6)=1d-2*Idint(tmax2*1d2+1d0)
c.......................................................................
        Call hbook2(Iden,' ',1,real(xlow(6)),real(xhig(6)),
     .                       1,real(ylow(6)),real(yhig(6)),0.)
        Call hminim(Iden,0.)
        Call hmaxim(Iden,1.)
c.......................................................................
        Call hplot(Iden,' ',' ',0)                         !..plot graph
        Call hplax(chtitx(6),chtity(6))                !..titles on axis
c.......................................................................
        Call igset('LTYP',2.00)                     !..changes line type
        Call Tunediag(xlow(6),xhig(6),ylow(6),yhig(6)) !..draws tunediag
c.......................................................................
        Call igset('LTYP',1.00)                     !..changes line type
        Call hpline(Tunx,Tuny,Icount2,' ')                 !..plots line
c.......................................................................
      endif
c.......................................................................
      iproj=iprojt                    !..restores initial value of iproj
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Tune_plot(iturn)
c.......................................................................
c.... Subroutine to plot the horizontal tune. It is used to add the tune
c.... to the evolution of the beam distribution. It is meant for movies.
c....
c.... NB: some parameters are artificially modified to reuse some 
c....     routines without modifications.
c....
c.... Author: M. Giovannozzi - CERN
c....  
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real hmal,xgsiz,ygsiz,scft,scfw
      Real Tunx(Max_dat),Turn(Max_dat)
c.......................................................................
      Character*4 cmap
      Character*11 chtitxx
      Character*13 chtit
      Character*17 chtityy
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Logical hexist
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
      Common/pawc/hmal(nplot)
c.......................................................................
      if (Nturn.gt.50) then
        iresidue=(Nturn-4*Iext)/5    !..defines max turn for opt. label
        iexp=Dlog10(Dfloat(iresidue))
        amagn=10**iexp
        imant=iresidue/amagn
        Antmax=5d0*amagn*(imant+1)
      else
        Antmax=Dfloat(Nturn-4*Iext)
      endif
c.......................................................................
      amaxnux=Dmax1(par(1),par(3))
      aminnux=Dmin1(par(1),par(3))
      Dnux=(amaxnux-aminnux)/pi2  !..hor. tune
      Dinux=1d0/Dnux
      iexp=Dlog10(Dinux)+1
      amagn=10**iexp
      imant=Dnux*amagn/5
      Dnux=5d0/amagn*(imant+1)
      centre=Idint(5d-1*(aminnux+amaxnux)/pi2*1d3)*1d-3
      tmin1=centre-Dnux/2                                   
      tmax1=centre+Dnux/2 
c.......................................................................
      xlow=0d0                      !..graphical windows for plotting
      ylow=tmin1
      xhig=Antmax
      yhig=tmax1   
c.......................................................................
      chtitxx='Turn number'
      chtityy='Hor. tune        '
c.......................................................................
      iprojt=iproj
      iproj=1             !..warning: artificial redefinition of iproj!!
      Call Opt_graph(1)
      iproj=7             !..warning: artificial redefinition of iproj!!
c.......................................................................
      Call hplzon(Nxz(iproj),Nyz(iproj),1,' ')          !..defines zones
c..............................................................Hor. tune
      if (hexist(Iden)) Call hdelet(Iden) !..deletes existing histogram
c.......................................................................
      Call htitle(' ')
      Call hbook2(Iden,' ',1,real(xlow),real(xhig),
     .                     1,real(ylow),real(yhig),0.)
      Call hminim(Iden,0.)
      Call hmaxim(Iden,1.)
c.......................................................................
      Call hplot(Iden,' ',' ',0)                           !..plot graph
      Call hplax(chtitxx,chtityy)                      !..titles on axis
c.......................................................................
      Call igset('LTYP',1.00)                       !..changes line type
      Call igset('LWID',8.00) !..changes line thickness
      Call igset('PLCI',2.00) !..changes line colour
c.......................................................................
      Tunx(1)=par(1)/pi2
      Turn(1)=0
      Istep=Nturn/1000
      Icount=1
      do kturn=0,iturn,Istep                          !..loop over turns
c.......................................................................
        if (kturn.le.ipar(4)) then                          !..tune ramp
          freq=par(1)
        elseif ((kturn.gt.ipar(4)).and.(kturn.le.ipar(1))) then
          freq=par(2)+par(5)*Dexp(par(7)*Dlog(Dfloat(ipar(1)-kturn)))
        elseif ((kturn.gt.ipar(1)).and.(kturn.le.ipar(2))) then
          freq=par(2)
        elseif ((kturn.gt.ipar(2)).and.(kturn.le.ipar(3))) then
          freq=par(2)+par(6)*Dexp(par(31)*Dlog(Dfloat(kturn-ipar(2))))
        elseif ((kturn.gt.ipar(3)).and.(kturn.le.ipar(5))) then
          freq=par(3)
        elseif ((kturn.gt.ipar(5)).and.(kturn.le.ipar(6))) then
          freq=par(3)+par(10)*(kturn-ipar(5))
        elseif ((kturn.gt.ipar(6)).and.(kturn.le.ipar(7))) then
          freq=par(1)
        endif
c.......................................................................
        Icount=Icount+1
c.......................................................................
        Tunx(Icount)=freq/pi2
        Turn(Icount)=kturn
c.......................................................................
      enddo
c.......................................................................
      if (Icount.ge.1002) Icount=1001
      Call hpline(Turn,Tunx,Icount,' ')                    !..plots line
c.......................................................................
      iproj=iprojt                    !..restores initial value of iproj
c.......................................................................
      Return
c.......................................................................
      End

      Double Precision Function Transfx(xmin,xmax,x)
c.......................................................................
c.... Trivial function to define a linear transformation of x. It is 
c.... used to generate coordinates in cm for HPLSOF
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      scalex=scx/(xmax-xmin)
      Transfx=xoff+scalex*(x-xmin)
c.......................................................................
      Return
c.......................................................................
      End
 
      Double Precision Function Transfy(ymin,ymax,y)
c.......................................................................
c.... Trivial function to define a linear transformation of y. It is 
c.... used to generate coordinates in cm for HPLSOF
c....
c.... Author: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw
c.......................................................................
      Character*4 cmap
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      scaley=scy/(ymax-ymin)
      Transfy=xoff+scaley*(y-ymin)
c.......................................................................
      Return
c.......................................................................
      End

      Subroutine Tunediag(xmin,xmax,ymin,ymax)
c.......................................................................
c.... Creates two files in order to draw the resonant lines in the tune 
c.... diagram. All the resonant lines from order Iorl to order 
c.... Ioru and drawn. xmin, xmax, ymin and ymax are the coordinates 
c.... of the square in the tune where one wants to plot the diagram
c....
c.... Author: E. Todesco - INFN of Bologna
c....
c.... This routine has been modified in order to fit the needs of our
c.... analysis
c....
c.... Editor: M. Giovannozzi - CERN
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      Real xgsiz,ygsiz,scft,scfw,resx(Max_dat),resy(Max_dat),xcoor,ycoor
c.......................................................................
      Character*4 cmap
      Character*5 texts
      Character*6 textl
      Character*13 chtit
      Character*29 script
      Character*30 fnamei
      Character*40 fnameo
c.......................................................................
      Common/Para/pi2,par(Max_par),bbox(4),ipar(Max_par),Nturn,Nout,
     .            Noff,Icomp,Idist,Ntot,Npmult,Idim,Iext,Iorl,Ioru,Mult,
     .            Idump,Itwiss,Ico,iramp,iblowup,Iper,Ninj,Icentre,
     .            iscriptl(0:4),iscriptu(0:4),script(0:4)
      Common/Graph/xgsiz(Max_pro),ygsiz(Max_pro),scft(Max_pro),
     .             scfw(Max_pro),Stp(Max_dim),Stpi(Max_dim),
     .             Ianim,Iout,Iden,Nxz(Max_pro),Nyz(Max_pro),
     .             iproj,icpr(Max_pro2),scx,scy,xoff,yoff,
     .             limitl(Max_boxes,Max_dim),limitu(Max_boxes,Max_dim),
     .             fnamei,fnameo(0:Max_pro),
     .             chtit(Max_dim,Max_dim),cmap(0:Max_dim)
c.......................................................................
      epsilon=(xmax-xmin)/80
c.......................................................................
      resx(1)=xmin
      resy(1)=ymin
      resx(2)=xmin
      resy(2)=ymax
      resx(3)=xmax
      resy(3)=ymax
      resx(4)=xmax
      resy(4)=ymin
      resx(5)=xmin
      resy(5)=ymin
      icont=5
c.......................................................................
      do iq=Iorl,Ioru
        do ip=1,iq-1
          if (Dfloat(ip)/iq.gt.xmin.and.Dfloat(ip)/iq.lt.xmax) then
            Call simpl1(ip,iq,ipp,iqq)
            if (iqq.lt.Iorl.or.iqq.eq.iq) then
              resx(icont+1)=Dfloat(ip)/iq
              resy(icont+1)=ymin
              resx(icont+2)=Dfloat(ip)/iq
              resy(icont+2)=ymax
              resx(icont+3)=xmin
              resy(icont+3)=ymax
              resx(icont+4)=xmin
              resy(icont+4)=ymin
              write(texts(1:5),'(a1,i1,a1,i1,a1)') '(',iq,',',0,')'
              xcoor=real(Transfx(xmin,xmax,Dfloat(ip)/iq-epsilon))
              ycoor=real(Transfy(ymin,ymax,ymin-4.5*epsilon))
              Call hplsof(xcoor,ycoor,texts,0.24,90.,0.,-1)
              icont=icont+4
            endif
          endif
        enddo
      enddo
c.......................................................................
      do iq=Iorl,Ioru
        do ip=1,iq-1
          if (Dfloat(ip)/iq.gt.ymin.and.Dfloat(ip)/iq.lt.ymax) then
            Call simpl1(ip,iq,ipp,iqq)
            if (iqq.lt.Iorl.or.iqq.eq.iq) then
              resx(icont+1)=xmin
              resy(icont+1)=Dfloat(ip)/iq
              resx(icont+2)=xmax
              resy(icont+2)=Dfloat(ip)/iq
              resx(icont+3)=xmax
              resy(icont+3)=ymin
              resx(icont+4)=xmin
              resy(icont+4)=ymin
              write(texts(1:5),'(a1,i1,a1,i1,a1)') '(',0,',',iq,')'
              xcoor=real(Transfx(xmin,xmax,xmin+epsilon))
              ycoor=real(Transfy(ymin,ymax,Dfloat(ip)/iq-4.5*epsilon))
              Call hplsof(xcoor,ycoor,texts,0.22,0.,0.,-1)
              icont=icont+4
            endif
          endif
        enddo
      enddo
c.......................................................................
      do i=Iorl,Ioru
c.......................................................................
        do ip=1,i-1
          resx(icont+1)=xmin
          resy(icont+1)=ymax
          icont=icont+1
          iq=i-ip
c.......................................................................
          do l=1,ip+iq-1
            Call simpl2(ip,iq,l,ipp,iqq,ll)
c.......................................................................
            if (ipp+iqq.lt.Iorl.or.iqq.eq.iq) then
c.......................................................................
              iang=-atan2(Dfloat(ip),Dfloat(iq))/pi2*360
              x=(l-iq*ymax)/ip
              y=(l-ip*xmin)/iq
c.......................................................................
              if (x.ge.xmin.and.x.le.xmax) then
c.......................................................................
                write(texts(1:5),'(a1,i1,a1,i1,a1)') '(',ip,',',iq,')'
                xcoor=real(
     .          Transfx(xmin,xmax,x+2*epsilon*Dsin(iang/360.*pi2)))
                ycoor=real(
     .          Transfx(ymin,ymax,ymax-2*epsilon*Dcos(iang/360.*pi2)-5*
     .                  epsilon))
                Call hplsof(xcoor,ycoor,texts,0.22,real(iang),0.,-1)
                resx(icont+1)=x
                resy(icont+1)=ymax
c.......................................................................
              elseif (y.ge.ymin.and.y.le.ymax) then
c.......................................................................
                write(texts(1:5),'(a1,i1,a1,i1,a1)') '(',ip,',',iq,')'
                xcoor=real(
     .              Transfx(xmin,xmax,xmin-epsilon*Dsin(iang/360.*pi2)))
                ycoor=real(
     .              Transfx(ymin,ymax,y+epsilon*Dcos(iang/360.*pi2)-5*
     .                      epsilon))
                Call hplsof(xcoor,ycoor,texts,0.22,real(iang),0.,-1)
                resx(icont+1)=xmin
                resy(icont+1)=y
c.......................................................................
              endif
              x=(l-iq*ymin)/ip
              y=(l-ip*xmax)/iq
              if (x.ge.xmin.and.x.le.xmax) then
c.......................................................................
                resx(icont+2)=x
                resy(icont+2)=ymin
                resx(icont+3)=xmin
                resy(icont+3)=ymin
                resx(icont+4)=xmin
                resy(icont+4)=ymax
                icont=icont+4
c.......................................................................
              elseif (y.ge.ymin.and.y.le.ymax) then
                resx(icont+2)=xmax
                resy(icont+2)=y
                resx(icont+3)=xmax
                resy(icont+3)=ymax
                resx(icont+4)=xmin
                resy(icont+4)=ymax
                icont=icont+4
c.......................................................................
              endif
c.......................................................................
            endif
          enddo
c.......................................................................
          iq=ip-i
          resx(icont+1)=xmin
          resy(icont+1)=ymin
          icont=icont+1
c.......................................................................
          do l=iq+1,ip-1
            Call simpl2(ip,iq,l,ipp,iqq,ll)
c.......................................................................
            if (ipp-iqq.lt.Iorl.or.iqq.eq.iq) then
c.......................................................................
              iang=atan2(Dfloat(ip),-Dfloat(iq))/pi2*360
              x=(l-iq*ymin)/ip
              y=(l-ip*xmin)/iq
c.......................................................................
              if (x.ge.xmin.and.x.le.xmax) then
c.......................................................................
                write(textl(1:6),'(a1,i1,a1,i2,a1)') '(',ip,',',iq,')'
                xcoor=real(
     .             Transfx(xmin,xmax,x-epsilon*Dsin(iang/360.*pi2)))
                ycoor=real(
     .             Transfy(ymin,ymax,ymin+epsilon*Dcos(iang/360.*pi2)-5*
     .                     epsilon))
                Call hplsof(xcoor,ycoor,textl,0.22,real(iang),0.,-1)
                resx(icont+1)=x
                resy(icont+1)=ymin
c.......................................................................
              elseif (y.ge.ymin.and.y.le.ymax) then
c.......................................................................
                write(textl(1:6),'(a1,i1,a1,i2,a1)') '(',ip,',',iq,')'
                xcoor=real(
     .             Transfx(xmin,xmax,xmin+2*epsilon*sin(iang/360.*pi2)))
                ycoor=real(
     .             Transfy(ymin,ymax,y-2*epsilon*cos(iang/360.*pi2)-5*
     .                     epsilon))
                Call hplsof(xcoor,ycoor,textl,0.22,real(iang),0.,-1)
                resx(icont+1)=xmin
                resy(icont+1)=y
c.......................................................................
              endif
c.......................................................................
              x=(l-iq*ymax)/ip
              y=(l-ip*xmax)/iq
c.......................................................................
              if (x.ge.xmin.and.x.le.xmax) then
c.......................................................................
                resx(icont+2)=x
                resy(icont+2)=ymax
                resx(icont+3)=xmin
                resy(icont+3)=ymax
                resx(icont+4)=xmin
                resy(icont+4)=ymin
                icont=icont+4
c.......................................................................
              elseif (y.ge.ymin.and.y.le.ymax) then
c.......................................................................
                resx(icont+2)=xmax
                resy(icont+2)=y
                resx(icont+3)=xmax
                resy(icont+3)=ymin
                resx(icont+4)=xmin
                resy(icont+4)=ymin
                icont=icont+4
c.......................................................................
              endif
c.......................................................................
            endif
c.......................................................................
          enddo
c.......................................................................
        enddo
c.......................................................................
      enddo
c.......................................................................
      Call hpline(resx,resy,icont,' ')      !..plots the resonance lines
c.......................................................................
      END

      Subroutine Simpl1(ip,iq,ipp,iqq)
c.......................................................................
c.... Given two integer numbers ip and iq, writes in ipp and iqq
c....  1)  ip and iq, if have no common divisor
c....  2)  ip/i and iq/i, where i is the smallest common divisor
c....      (greater than one)
c....
c.... Author: E. Todesco - INFN of Bologna
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      ipp=ip
      iqq=iq
      do i=2,ip
        iscr1=ip/i
        iscr2=iq/i
        if (iscr1*i.eq.ip.and.iscr2*i.eq.iq) then
          ipp=iscr1
          iqq=iscr2
          goto 20
        endif
      enddo
 20   continue
c.......................................................................
      end
 
      Subroutine Simpl2(ip,iq,l,ipp,iqq,ll)
c.......................................................................
c.... Given three integer numbers ip,iq and l, writes in ipp, iqq and ll
c....
c....  1)  ip, iq and ll, if they have no common divisor
c....  2)  ip/i, iq/i and ll/i, where i is the smallest common
c....      divisor (greater  than one)
c....
c.... Author: E. Todesco - INFN of Bologna
c....
c.......................................................................
      Implicit Double Precision (A-H,O-Z)
      Implicit Integer (I-N)
c.......................................................................
      Parameter(Max_part=8700000,Max_par=100,Max_bord=10000,
     .          Max_boxes=20,Max_bin=500,Max_dim=4,Max_pro=9,
     .          Max_pro2=2*Max_pro,Max_box2=2*Max_boxes,Max_dat=2000,
     .          Max_box4=2*Max_box2,Max_e1D=Max_bin*Max_box2,
     .          Max_e2D=Max_pro*Max_bin*Max_bin,tol=1d-08,nplot=100000)
c.......................................................................
      ipp=ip
      iqq=iq
      ll=l
      do i=2,ip
        iscr1=ip/i
        iscr2=iq/i
        iscr3=l/i
        if (iscr1*i.eq.ip.and.iscr2*i.eq.iq.and.iscr3*i.eq.l) then
          ipp=iscr1
          iqq=iscr2
          ll=iscr3
          goto 20
        endif
      enddo
 20   continue
c.......................................................................
      end

*
* $Id: ranlux.F,v 1.2 1997/09/22 13:45:47 mclareni Exp $
*
* $Log: ranlux.F,v $
* Revision 1.2  1997/09/22 13:45:47  mclareni
* Correct error in initializing RANLUX by using RLUXIN with the output of
* RLUXUT from a previous run.
*
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
C#include "pilot.h"
      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END

*
* $Id: rnormx.F,v 1.1.1.1 1996/04/01 15:02:55 mclareni Exp $
*
* $Log: rnormx.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
C#include "pilot.h"
      SUBROUTINE RNORMX(DEVIAS,NDEV,ROUTIN)
C        Generator of a vector of independent Gaussian-distributed 
C        (pseudo-)random numbers, of mean zero and variance one,
C        making use of a uniform pseudo-random generator (RANMAR).
C        The algorithm for converting uniform numbers to Gaussian
C        is that of "Ratio of Uniforms with Quadratic Bounds."  The
C        method is in principle exact (apart from rounding errors),
C        and is based on the variant published by Joseph Leva in
C        ACM TOMS vol. 18(1992), page 449 for the method and 454 for
C        the Fortran algorithm (ACM No. 712).
C        It requires at least 2 and on average 2.74 uniform deviates
C        per Gaussian (normal) deviate.
C   WARNING -- The uniform generator should not produce exact zeroes,
C   since the pair (0.0, 0.5) provokes a floating point exception.
      SAVE  S, T, A, B, R1, R2
      DIMENSION U(2), DEVIAS(*)
      EXTERNAL ROUTIN
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2/ 0.27597, 0.27846/
C         generate pair of uniform deviates
      DO 200 IDEV = 1, NDEV
   50 CALL ROUTIN(U,2)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C           accept P if inside inner ellipse
      IF (Q .LT. R1)  GO TO 100
C           reject P if outside outer ellipse
      IF (Q .GT. R2)  GO TO 50
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)
  200 DEVIAS(IDEV) = DEVIAT
      RETURN
      END

