MOD_XRDSF
Requirements:
  Argparse;
  Pandas;
  Matplotlib;
Compile:
  f2py -c -m xrdsflib  sub.f90
Usage:
  ./modified_xrdsf.py -h : check help info
  ./modified_xrdsf.py test : plot test.vasp sq
  ./modified_xrdsf.py test -vb 0 : plot test.vasp sq without any verbs
  ./modified_xrdsf.py test -is_plot F : output test.vasp sq in text mode
