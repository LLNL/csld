      real*8 function derfc (x)
c  adapted from july 1977 edition, w. fullerton, c3,
c  los alamos scientific lab.  supplied by mark.
      implicit real*8 (a-h,o-z)
      dimension erfcs(13), erfccs(24), erc2cs(23), d1mach(3)
      real*8 eta
      external dcsevl, initds
      save
c  d1mach(1)=smallest non-zero number.
c  d1mach(3)=smallest number which can be added to 1
      data d1mach /1.d-30, 0.d0, 1.d-15 /

      data erf cs( 1) /   -.0490461212 34691808d0 /
      data erf cs( 2) /   -.1422612051 0371364d0 /
      data erf cs( 3) /    .0100355821 87599796d0 /
      data erf cs( 4) /   -.0005768764 69976748d0 /
      data erf cs( 5) /    .0000274199 31252196d0 /
      data erf cs( 6) /   -.0000011043 17550734d0 /
      data erf cs( 7) /    .0000000384 88755420d0 /
      data erf cs( 8) /   -.0000000011 80858253d0 /
      data erf cs( 9) /    .0000000000 32334215d0 /
      data erf cs(10) /   -.0000000000 00799101d0 /
      data erf cs(11) /    .0000000000 00017990d0 /
      data erf cs(12) /   -.0000000000 00000371d0 /
      data erf cs(13) /    .0000000000 00000007d0 /
c
      data erc2cs( 1) /   -.0696013466 02309501d0 /
      data erc2cs( 2) /   -.0411013393 62620893d0 /
      data erc2cs( 3) /    .0039144958 66689626d0 /
      data erc2cs( 4) /   -.0004906395 65054897d0 /
      data erc2cs( 5) /    .0000715747 90013770d0 /
      data erc2cs( 6) /   -.0000115307 16341312d0 /
      data erc2cs( 7) /    .0000019946 70590201d0 /
      data erc2cs( 8) /   -.0000003642 66647159d0 /
      data erc2cs( 9) /    .0000000694 43726100d0 /
      data erc2cs(10) /   -.0000000137 12209021d0 /
      data erc2cs(11) /    .0000000027 88389661d0 /
      data erc2cs(12) /   -.0000000005 81416472d0 /
      data erc2cs(13) /    .0000000001 23892049d0 /
      data erc2cs(14) /   -.0000000000 26906391d0 /
      data erc2cs(15) /    .0000000000 05942614d0 /
      data erc2cs(16) /   -.0000000000 01332386d0 /
      data erc2cs(17) /    .0000000000 00302804d0 /
      data erc2cs(18) /   -.0000000000 00069666d0 /
      data erc2cs(19) /    .0000000000 00016208d0 /
      data erc2cs(20) /   -.0000000000 00003809d0 /
      data erc2cs(21) /    .0000000000 00000904d0 /
      data erc2cs(22) /   -.0000000000 00000216d0 /
      data erc2cs(23) /    .0000000000 00000052d0 /
c
      data erfccs( 1) /   0.0715179310 202925d0 /
      data erfccs( 2) /   -.0265324343 37606719d0 /
      data erfccs( 3) /    .0017111539 77920853d0 /
      data erfccs( 4) /   -.0001637516 63458512d0 /
      data erfccs( 5) /    .0000198712 93500549d0 /
      data erfccs( 6) /   -.0000028437 12412769d0 /
      data erfccs( 7) /    .0000004606 16130901d0 /
      data erfccs( 8) /   -.0000000822 77530261d0 /
      data erfccs( 9) /    .0000000159 21418724d0 /
      data erfccs(10) /   -.0000000032 95071356d0 /
      data erfccs(11) /    .0000000007 22343973d0 /
      data erfccs(12) /   -.0000000001 66485584d0 /
      data erfccs(13) /    .0000000000 40103931d0 /
      data erfccs(14) /   -.0000000000 10048164d0 /
      data erfccs(15) /    .0000000000 02608272d0 /
      data erfccs(16) /   -.0000000000 00699105d0 /
      data erfccs(17) /    .0000000000 00192946d0 /
      data erfccs(18) /   -.0000000000 00054704d0 /
      data erfccs(19) /    .0000000000 00015901d0 /
      data erfccs(20) /   -.0000000000 00004729d0 /
      data erfccs(21) /    .0000000000 00001432d0 /
      data erfccs(22) /   -.0000000000 00000439d0 /
      data erfccs(23) /    .0000000000 00000138d0 /
      data erfccs(24) /   -.0000000000 00000048d0 /
c
      data sqrtpi /1.772453850 9055160d0/
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0.d0/
c
      if (nterf .ne. 0d0) goto 10
      eta = 0.1*sngl(d1mach(3))
      nterf = initds (erfcs, 13, eta)
      nterfc = initds (erfccs, 24, eta)
      nterc2 = initds (erc2cs, 23, eta)
      xsml = -dsqrt (-dlog(sqrtpi*d1mach(3)))
      xmax = dsqrt (-dlog(sqrtpi*d1mach(1)))
      xmax = xmax - 0.5d0*dlog(xmax)/xmax - 0.01d0
      sqeps = dsqrt (2.0d0*d1mach(3))
   10 if (x .gt. xsml) goto 20
      derfc = 2.d0
      return
   20 if (x .gt. xmax) goto 40
      y = dabs(x)
      if (y .gt. 1.0d0) goto 30
      if (y .lt. sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y .ge. sqeps) derfc = 1.0d0 -
     .  x*(1.0d0 + dcsevl (2.d0*x*x-1.d0, erfcs, nterf) )
      return
   30 y = y*y
      if (y .le. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl ((8.d0/y-5.d0)/3.d0, erc2cs, nterc2) )
      if (y .gt. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl (8.d0/y-1.d0, erfccs, nterfc) )
      if (x .lt. 0.d0) derfc = 2.0d0 - derfc
      return
  40  continue
c|    write(6,*) 'derfc: underflow'
      derfc = 0.d0
      return
      end
c ----------- function initds ------------------
      integer function initds (dos, nos, eta)
      implicit real*8 (a-h,o-z)
c- initialize things for chebychev series
      real*8 dos(nos)
      real*8 eta, err
      err = 0.
      do  10  ii = 1, nos
        i = nos + 1 - ii
        err = err + dabs(dos(i))
        if (err .gt. eta) goto 20
   10 continue
   20 continue
c     if (i .eq. nos) write(6,*) 'initds: eta may be too small'
      initds = i
      return
      end
c ----------- function dcsevl ------------------
      real*8 function dcsevl (x, a, n)
      implicit real*8 (a-h,o-z)
c- evaluate the n-term chebyshev series a at x.
      real*8 a(n), x, twox, b0, b1, b2
      if (dabs(x) .gt. 1.1d0) write(6,*) 'dcsevl:  x outside (-1,1)'
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do  10  i = 1, n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
   10 continue
      dcsevl = .5d0*(b0-b2)
      return
      end
c ---------- derf ---------------------
      real*8 function derf (x)
      implicit real*8 (a-h,p-z), integer(o)
      derf=1.d0-derfc(x)
      return
      end
