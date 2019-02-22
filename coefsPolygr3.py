from sympy import *

L  = symbols('L')
w1, w2  = symbols('w1 w2')
t1, t2  = symbols('t1 t2')
x  = symbols('x')
E  = symbols('E')
I  = symbols('I')

latexflag = 0

pprint("===================================")

# convencion reddy
#~ M = Matrix( [ [ 0, 0, 0, 1], [0, 0, -1, 0 ], [L**3, L**2, L, 1 ], [ -3*L**2, -2*L, -1,0 ] ] )

# convencion: w1, t1, w2, t2 desplazamientos hacia arriba y giros antihorario.
# matriz de condiciones para calculo de funciones de interpolacion
# v(0) theta(0) v(L) theta(L):
Matcond = Matrix( [ [ 0, 0, 0, 1], [0, 0, 1, 0 ], [L**3, L**2, L, 1 ], [ 3*L**2, 2*L, 1,0 ] ] )
# ~ Matcond = Matrix( [ [ 0, 0, 0, 1], [0, 0, -1, 0 ], [L**3, L**2, L, 1 ], [ -3*L**2, -2*L, -1,-0 ] ] )

Matcondinv = Matcond.inv()

# vector de grados de libertad
U = Matrix( [ [ w1 ], [t1 ], [w2 ], [t2 ] ] )

pprint("  Coeficientes despejados: ")
pprint(Matcondinv*U)
pprint("----------------------------------")
if (latexflag>0):
  pprint("  Codigo LaTeX de coeficientes despejados: ")
  pprint( latex( Matcondinv*U ) )
  pprint("----------------------------------")


input("  Presione una tecla para continuar.")

# vector de monomios
X = Matrix( [ [ x**3, x**2, x, 1 ] ] )

# funciona de flecha
w = X * Matcondinv * U

pprint("  Expresion de la flecha: ")
pprint( simplify(factor( w,x)) )
pprint("----------------------------------")
if (latexflag>0):
  pprint("  Expresion LaTeX de la flecha: ")
  pprint( latex(factor( w,x)) )
  pprint("----------------------------------")
  pprint( latex(w) )
  pprint("----------------------------------")

input("  Presione una tecla para continuar.")


# para ploteo
# ~ f1 = w.subs( [ (L,1), (w1, 1), (t1, 0), (w2, 0), (t2, 0) ] )
# ~ f2 = w.subs( [ (L,1), (w1, 0), (t1, 1), (w2, 0), (t2, 0) ] )
# ~ f3 = w.subs( [ (L,1), (w1, 0), (t1, 0), (w2, 1), (t2, 0) ] )
# ~ f4 = w.subs( [ (L,1), (w1, 0), (t1, 0), (w2, 0), (t2, 1) ] )

#~ plot(  (x,0,1))
# ~ plot( f1[0],f3[0], (x,0,1))

# ~ plot( f2[0],f4[0], (x,0,1))
#~ plot(  (x,0,1))
#~ plot( f4[0], (x,0,1))


# funciones de flecha en funcion de L
fL1 = w.subs( [  (w1, 1), (t1, 0), (w2, 0), (t2, 0) ] )
fL2 = w.subs( [  (w1, 0), (t1, 1), (w2, 0), (t2, 0) ] )
fL3 = w.subs( [  (w1, 0), (t1, 0), (w2, 1), (t2, 0) ] )
fL4 = w.subs( [  (w1, 0), (t1, 0), (w2, 0), (t2, 1) ] )


ddwdx = diff( diff(w,x),x)
pprint("  Curvatura: ")
pprint( ddwdx )
pprint("----------------------------------")

ddf1dx = ddwdx.subs([(w1, 1), (t1, 0), (w2, 0), (t2, 0)])
ddf2dx = ddwdx.subs([(w1, 0), (t1, 1), (w2, 0), (t2, 0)])
ddf3dx = ddwdx.subs([(w1, 0), (t1, 0), (w2, 1), (t2, 0)])
ddf4dx = ddwdx.subs([(w1, 0), (t1, 0), (w2, 0), (t2, 1)])

pprint("  Expresiones de curvaturas nodales: ")
pprint( ddf1dx)
pprint( ddf2dx)
pprint( ddf3dx)
pprint( ddf4dx)
pprint("----------------------------------")

input("  Presione una tecla para continuar...")


V1 = integrate(fL1, (x, 0, L))
M1 = integrate(fL2, (x, 0, L))
V2 = integrate(fL3, (x, 0, L))
M2 = integrate(fL4, (x, 0, L))


if (latexflag >0):
  pprint("  Expresion LaTeX de momentos nodales =====")
  pprint( latex(V1))
  pprint( latex(M1))
  pprint( latex(V2))
  pprint( latex(M2))


fsmat = Matrix( [ ddf1dx, ddf2dx, ddf3dx, ddf4dx ] )

qsmat = Matrix( [ fL1, fL2, fL3, fL4 ] )


K11 = E*I*integrate(ddf1dx*ddf1dx, (x, 0, L))
pprint( K11 )

K12 = E*I*integrate(ddf1dx*ddf2dx, (x, 0, L))
pprint( K12 )

K = E*I*integrate( simplify( fsmat * fsmat.T ) , (x,0,L))

Nq =integrate( simplify( qsmat ) , (x,0,L))
 
pprint("  Matriz de rigidez")

pprint( K )
if (latexflag >0):
  pprint( latex( K) )

input("  Presione una tecla para continuar...")


# ~ pprint( K.row(1) )
t1art = solve( K.row(1)*U, t1)[t1]
# ~ pprint( K.row(3) )
t2art = solve( K.row(3)*U, t2)[t2]

pprint("  Giro en nodo inicial con articulacion en ese nodo:")
pprint( t1art )

pprint("  Giro en nodo final con articulacion en ese nodo:")
pprint( t2art )

aux = t1art

vectorfacts1 = zeros(4,1)
vectorfacts2 = zeros(4,1)
for i in range(4):
  vectorfacts1[i] = diff( t1art, U[i] )
  vectorfacts2[i] = diff( t2art, U[i] )

Kart1 = zeros(4)
Kart2 = zeros(4)

for i in range(4):
  if (i!=1):
    Kart1[:,i] = simplify( expand( K.col(i) + K.col(1) * vectorfacts1[i] ) )
  if (i!=3):
    Kart2[:,i] = simplify( expand( K.col(i) + K.col(3) * vectorfacts2[i] ) )

pprint(" vector facts1: " )
pprint( vectorfacts1 )
pprint(" vector facts2: " )
pprint( vectorfacts2 )

input("  Presione una tecla para ver matrices...")

pprint( simplify( cancel( expand( Kart1) ) ) )
pprint( simplify( cancel( expand( Kart2) ) ) )

if (latexflag>0):
  pprint( latex( Nq) )


