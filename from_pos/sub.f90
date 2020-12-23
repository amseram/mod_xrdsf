!######################################################
!# Original XRDSF url: https://github.com/rraffiu/XRDSF
!# Author. : Amsera
!# Date.   : 21st Dec. 2020
!######################################################
Module xrdsf_lib
Contains
  Function Calc_intents(msd, natm, atm_loc, curg, prma, prmb, prmc) Result (sq)
    Implicit None
    Integer, Intent(in) :: natm
    Real(8), Intent(in) :: atm_loc(natm,3), prma(natm,4), prmb(natm,4), prmc(natm), curg(3), msd
    Integer    i
    Real(8)    nrmg, sq(2), tsq, dwf, proj(natm)
    Complex(8) v, prow
    v    = cmplx(0.,1.)
    Call Calc_proj(natm,atm_loc,curg,proj)
    nrmg = Norm(curg)
    dwf  = Exp(-1.d0*msd*nrmg*nrmg/3.d0)
    prow = 0.d0
    Do i = 1, natm
      prow = prow + Calc_eaff(prma(i,:),prmb(i,:),prmc(i),nrmg)*Exp(proj(i)*v)
    Enddo
    tsq  = Abs(prow*dwf)
    tsq  = tsq * tsq
    sq(1)= nrmg
    sq(2)= tsq
    Return 
  End function Calc_intents
  
  Subroutine Calc_proj(natm,aloc,curg,proj) 
    Implicit None
    Integer, Intent(in) :: natm
    Real(8), Intent(in) :: aloc(natm,3), curg(3)
    Real(8), Intent(inout) :: proj(natm)
    Integer i
    Do i = 1, natm
      proj(i) = Dot_product(aloc(i,:),curg)
    Enddo
    Return
  End subroutine Calc_proj
  
  Function Calc_eaff(a,b,c,g) Result(f)
    Implicit None
    Real(8), Intent(in) :: a(4), b(4), c, g
    Real(8) tmp,f,pi 
    Integer i
    pi  = 3.1415926d0
    tmp = -1.d0*((0.25*g/pi)**2.d0)
    f   = c
    Do i = 1, 4
      f = f + a(i)*Exp(b(i)*tmp) 
    Enddo
    Return
  End function Calc_eaff
  
  Function Norm(v) Result(l)
    Implicit None
    Real(8), Intent(in) :: v(3)
    Real(8) l
    l = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
    l = Sqrt(l)
    Return
  End function Norm

  Function Get_g(rvec,indx) Result(cg)
    Implicit None
    Integer, Intent(in)  :: indx(3)
    Real(8), Intent(in)  :: rvec(3,3)
    Real(8) cg(3)
    cg = matmul(indx,rvec) 
    Return
  End function Get_g
End Module xrdsf_lib

Subroutine f_calc_ints(m,n,l,g,a,b,c,sq)
  Use xrdsf_lib
  Integer, Intent(in)  :: n
  Real(8), Intent(in)  :: m, l(n,3), g(3), a(n,4), b(n,4), c(n) 
  Real(8), Intent(out) :: sq(2)
!f2py intent(in)  :: m,m,l,g,a,b,c
!f2py intent(out) :: sq
  sq = Calc_intents(m,n,l,g,a,b,c)
  Return
End subroutine f_calc_ints

Subroutine f_loop(msd,mh,mk,ml,ng,rvec,na,aloc,prma,prmb,prmc,sq) 
  Use xrdsf_lib
  Integer, Intent(in)   :: mk, mh, ml, ng, na
  Real(8), Intent(in)   :: rvec(3,3),aloc(na,3),prma(na,4),prmb(na,4),prmc(na),msd
  Real(8), Intent(out)  :: sq(ng,5)
  Integer  h,k,l,ig
  Real(8)  cg(3)
!f2py intent(out) :: sq
  ig = 0
  Do h = -mh, mh
    Do k = -mk, mk
      Do l = -ml, ml
        ig = ig + 1
        sq(ig,1:3) = (/h,k,l/)
        cg = Get_g(rvec,(/h,k,l/))
        sq(ig,4:5) = Calc_intents(msd,na,aloc,cg,prma,prmb,prmc) 
      Enddo
    Enddo
  Enddo
  Return
End Subroutine f_loop
