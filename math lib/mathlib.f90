Module Mathlib
  contains
  Recursive Function det(A,col,row) Result(D)
    Integer row, col, idxRows(row-1), idxCOLS(col-1)
    Real(8) A(col,row), D
    Integer m, n, k, c, f
    m=row
    n=col
    If (m>1) Then
      D = 0;
      Do k=1,m
        idxCOLS = sitdiff([1:m],m,k)
        idxROWS = [2:n]
        c = col-1
        f = row-1
        D = D + (-1)**(1+k) * A(k,1) * det(A(idxCOLS,idxROWS),c,f)
      EndDo
    Else
      D = A(1,1)
    EndIf
  End Function det

  Function sitdiff(A,n,k) Result (B)
    Integer n,d,k
    Integer A(:),B(n-1)
    If(k .eq. 1) Then
      B = A(k+1:n)
    ElseIf(k .eq. n) Then
      B = A(1:k-1)
    Else
      B = A((/1:k-1,k+1:n/))
    EndIf
  End Function sitdiff

  Function inv(A,row) Result (inv_A)
    Integer row,m(row-1),n(row-1)
    Real(8) A(row,row),inv_A(row,row)
    Do i=1,row
      Do j=1,row
        m = sitdiff([1:row],row,i)
        n = sitdiff([1:row],row,j)
        inv_A(i,j) = (-1)**(i+j)*det(A(m,n),row-1,row-1)
      Enddo
    Enddo
    inv_A = inv_A/det(A,row,row)
    inv_A = tran(inv_A,row,row)
  End function inv

  Function arange(n) Result(a)
    Integer i, n, a(n)
    Do i = 1, n
      a(i) = i
    Enddo 
  End function arange

  Function tran(A,col,row) Result(t_A)
    Integer row,col
    Real(8) t_A(row,col),A(col,row)
    Do i=1,col
      Do j=1,row
        t_A(j,i)=A(i,j)
      Enddo
    Enddo
  End function tran

  Subroutine output_mat(A,row,col)
    Integer row,col
    Real(8) A(row,col)
    Character(12) ftmp 
    Write(ftmp,'(I3,A9)') col,"(F16.8,x)" 
    Do i = 1,row
       Write(6,ftmp) A(:,i)
    Enddo
  End subroutine output_mat

End module Mathlib 
