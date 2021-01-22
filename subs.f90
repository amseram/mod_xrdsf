Module VASP_lib
  Implicit None
  Private 

  Type Atom_info
  Character(2), Allocatable :: Atm_typ(:)
  Integer                   :: Atm_tot
  Integer,      Allocatable :: Atm_num(:)
  Integer,      Allocatable :: Atm_rnk(:)
  Real(8),      Allocatable :: Atm_loc(:,:) 
  End type Atom_info
  Type Cell_info
  Real(8) :: Cel_rat, Cel_vol
  Real(8) :: Cel_vec(3,3)
  Real(8) :: Cel_rvc(3,3)
  Real(8) :: Cel_len(3)
  Real(8) :: Cel_rln(3)
  Logical :: Isdirec = .True.
  Logical :: Isconst = .False.
  End type Cell_info
  Type Pos_info
  Type(Atom_info) :: atminf
  Type(Cell_info) :: celinf
  Character(8)    :: Fname
  Character(64)   :: Pname
  Integer         :: Fobj
  End type Pos_info
Contains
  Function Read_POSCAR (Fname) Result (Posinf)
    Implicit None
    Character(8), Intent(in  ) :: Fname
    Type(Pos_info)             :: Posinf
    Posinf%Fname = Fname
    Posinf%Fobj  = 99999
    Open (Posinf%Fobj, File=Fname, Status='Unknown')
    Read (Posinf%Fobj, '(a64)') Posinf%Pname
    Call Get_Cell_Info(Posinf)
    Call Get_Atom_Info(Posinf)
    Close (Posinf%Fobj)
    Return
  End function Read_POSCAR 

  Subroutine Get_Cell_Info (Pos)
    Type(Pos_info), Intent(inout) :: Pos 
    Integer i,j
    Character(64) tmp
    Read (Pos%Fobj, *) Pos%celinf%Cel_rat
    Do i = 1, 3
      Read (Pos%Fobj, *) (Pos%celinf%Cel_vec(j,i), j=1,3)
    Enddo
    Read (Pos%Fobj, '(a64)') tmp
    Backspace (Pos%Fobj)
    j = 0 
    Do i = 2, 64
      If (tmp(i:i).ne." " .and. tmp(i-1:i-1) .eq. " ") Then
        j = j + 1
      Endif
    Enddo
    Allocate(Pos%atminf%Atm_typ(j))
    Allocate(Pos%atminf%Atm_num(j))
    Read (Pos%Fobj, *) (Pos%atminf%Atm_typ(i), i =1, j)
    Read (Pos%Fobj, *) (Pos%atminf%Atm_num(i), i =1, j)
    Pos%atminf%Atm_tot = Sum(Pos%atminf%Atm_num)
    Allocate(Pos%atminf%Atm_loc(3,Pos%atminf%Atm_tot))
    Read (Pos%Fobj, '(a64)') tmp
    If (tmp(1:1) .eq. "S" .or. tmp(1:1) .eq. "s") Then
        Pos%celinf%Isconst = .True.
        Read (Pos%Fobj, '(a64)') tmp
    EndIf 
    If (tmp(1:1) .eq. "C" .or. tmp(1:1) .eq. "c") Then
        Pos%celinf%Isdirec = .False.
    Endif
    Return
  End subroutine Get_Cell_Info
  Subroutine Get_Atom_Info(Pos)
    Type(Pos_info), Intent(inout) :: Pos 
    Return
  End subroutine Get_Atom_Info
End Module VASP_lib
