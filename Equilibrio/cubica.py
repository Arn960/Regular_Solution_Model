# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 19:33:39 2021

@author: damv_

Ecuaciones cubicas de estado
"""
from scipy.special import gamma
import numpy as np
"""****************************************************************************************************************************
Indices para ecuaciones cubicas de estado;
    0 - Peng Robinson
    1 - Soave Riedlich Kwon
    2 - Van der Waals
    
R = 83.14472 #Barcm3/molK
En la rutina PseudoComp, se generan los pseudocomponentes a partir de la fracción de nC7+. Min se toma 
de la literatura (Tablas de SCN) para la masa molar de nC7 como fraccion SCN.

Metodo para generar pseudocomponentes:
    Distribucion: Gamma
    Metodo de integracion: Gauss - Laguerre
    
    Lista de datos requeridos:
        zpseud = Composicion de fraccion C7+
        Mprom = Masa molar promedio de fraccion C7+
        n = numero de pseudocomponentes, 4 es un buen valor para no forzar la PC.
        alfa = parametro de forma, un valor de entre 1 y 2 es buen estimado para crudos.
        Mmin = Masa molar minima para la fraccion C7+, ya que toma la masa molar del nC7 (SCN).

****************************************************************************************************************************"""

def PseudoComp( zpseud , Mprom , n  , alfa , Mmin = 90): #Entradas escalares
    z = np.array([ ])
    M = np.array([ ])
    X,w = np.polynomial.laguerre.laggauss( n ) #Puntos y pesos para Gauss Laguerre
    beta = ( 2.5*Mprom-Mmin )/X[ -1 ]
    exs = alfa*beta/( Mprom-Mmin )
    delta = np.exp( exs-1 )
    f = lambda j : ( X[ j ]**( alfa-1 ) )*( exs**alfa )/( gamma( alfa )*delta**X[ j ] )
    for i in range( n ):
        zi = zpseud*( w[ i ]*f( i ) )
        Mi = Mmin+beta*X[ i ]
        z = np.append( z , zi )
        M = np.append( M , Mi )
    cz = zpseud/sum( z ) #Correccion de z
    z = z*cz
    cM = Mprom*sum( z )/sum( M*z ) #Correccion de M
    M = M*cM
    return np.array([ z , M ])

"""****************************************************************************************************************************
Entradas en las ecuaciones:
    Las propiedades de substancias como Pc, Tc, w, etc,
    son np.array. El indice corresponde a una substancia especifica.

Indices para ecuaciones cubicas de estado;
    0 - Peng Robinson -> RECOMENDADA
    1 - Soave Riedlich Kwon
    2 - Van der Waals

Rutina para funcion alfa( Tr ):
    a) Alfa de Stryjek y Vera con parametro polar q - 1986 ( alfa_SV )

Se usan las reglas de mezclado de VdW.
Oma = np.array([ 0.457235530 , 0.42748 , 27/64 ])    Omb = np.array([ 0.077796074 , 0.08664 , 1/8 ])
    A = Oma[Ec]*alfa( Tr )*Pr/Tc**2
    B = Omb[Ec]*Pr/Tc
    Bm = B_i*x_i
    A_ij = ( 1 - k_ij)*( A_i*A_j )**0.5
    Am = x_i*x_j*A_ij
    A_i´ = x_j*A_ij
    
    Rutinas para los distintos parametros:
        a) Parametros B individuales ( B_i )
        b) B de mezcla, calculando Bi ( Bm_0 )
        c) Parametros A individuales ( A_i )
        d) A de mezcla, calculando Ai ( Am_0 )
        e) A de mezcla, otorgando Ai ( Am_1 )
        f) A´ de un componente en mezcla ( Aprim_i )

****************************************************************************************************************************"""

def B_i( P ,T , Pc, Tc, Ec = 0 ):
    Omb = np.array([ 0.077796074 , 0.08664 , 1/8 ])
    Pr = P/Pc
    Tr = T/Tc
    B = Pr*Omb[ Ec ]/Tr
    return B

def Bm_0( P ,T , Pc, Tc, x , Ec = 0 ):
    return np.dot( x , B_i( P ,T , Pc, Tc, Ec) )

#Alfa de Stryjek y Vera, 1986
def alfa_SV( T , Tc , w , q , Ec = 0 ):
    r = np.array([ [ 0.378893, 1.4897153, -0.17131848, 0.0196554 ] , [ 0.48508, 1.55191, -0.15613, 0 ] , [0,0,0,0] ])
    alfa = np.zeros([ Tc.size ])
    mi = lambda n: r[ Ec,0 ] + r[ Ec,1 ]*w[ n ] + r[ Ec,2 ]*np.power(w[n],2) + r[ Ec,3 ]*np.power(w[n],3)
    ci = lambda n: 1 + mi( n )/2 +0.3*q[ n ]
    for i in range( Tc.size ):
        if T/Tc[ i ] <= 1:
            alfa[ i ] = np.power( 1 + mi( i )*( 1 - np.power( T/Tc[ i ], 0.5 ) - q[ i ]*( 1 - T/Tc[ i ] )*( 0.7 - T/Tc[ i ] ) ), 2 )
        elif T/Tc[ i ] > 1:
            alfa[ i ] = np.exp( 2*( ci( i ) - 1 )*( 1 - np.power( T/Tc[ i ], ci(i) ) )/ci( i ) )
    return alfa

def A_i( P , T , Pc , Tc , w , q , Ec = 0 ):
    Oma = np.array([ 0.457235530 , 0.42748 , 27/64 ])
    # A = ( Oma[Ec]*P/np.power(T,2) )*alfa_SV( T , Tc , w , q , Ec )*np.power( Tc,2 )/Pc
    A = Oma[Ec]*(P/Pc)*np.power(Tc/T, 2)*alfa_SV( T , Tc , w , q , Ec )
    return A

def Am_0( P , T , Pc , Tc , w , q , x , kij ,  Ec = 0 ):
    Ai = np.power( A_i( P , T , Pc , Tc , w , q , Ec ) , 0.5 )
    Aij = np.outer( Ai , Ai )
    Kij = 1 - kij
    Aij = Aij*Kij
    # print( Aij )
    return np.dot( x.T , np.dot( x.T , Aij ) )

def Am_1( x , Ai, kij ):
    Ai = np.power( Ai, 0.5 )
    Aij = np.outer(Ai, Ai)*(1 - kij)
    return np.dot( x.T, np.dot( x.T , Aij ) )

def Aprim_i( P , T , Pc , Tc , w , q , x , kij ,  Ec = 0 ):
    Ai = np.power( A_i( P , T , Pc , Tc , w , q , Ec ) , 0.5 )
    Aij = np.outer( Ai , Ai )
    Kij = 1 - kij
    Aij = Aij*Kij
    return 2*np.dot( x.T , Aij )

def Aprim_i_1( x, Ai, kij ):
    Ai = np.power( Ai, 0.5 )
    Aij = np.outer( Ai, Ai )*(1 - kij)
    return 2*np.dot( x.T, Aij )

"""****************************************************************************************************************************
Ecuaciones Cubicas de Estado:
    0 - Peng Robinson
    1 - Soave Riedlich Kwon
    2 - Van der Waals

Rutina para solucion de Ec cubica, dando T,P y x ( z_0 )
Rutina para solucion de Ec cubica, dando Am y Bm ( z_1 )

    u = np.array([ 2 , 1 , 0 ])   v = np.array([ -1 , 0, 0 ])
 Polinomio: z**3 + (  (u[Ec] - 1)*B  - 1 )*z**2  ...
  ... + ( A - u[Ec]*B +(v[Ec] - u[Ec])*B**2 )*z - ( A*B + v[Ec]*( B**2 + B**3 ) )

Rutina necesaria para calcular fugacidad. f_i = P*phi_i*x_i
Rutinas para calcular coeficientes de fugacidad "phi":
        a) Coeficiente de fugacidad de un componente ( Phi_i )
        b) Lista de coeficientes de fugacidad para todos los componentes ( phi_ )
****************************************************************************************************************************"""

def z_0( P , T , Pc , Tc , w , q , x , kij , Fase , Ec = 0 ):
    Fase = str( Fase )
    u = np.array([ 2 , 1 , 0 ])
    v = np.array([ -1 , 0, 0 ])
    A = Am_0( P , T , Pc , Tc , w , q , x , kij ,  Ec )
    B = Bm_0( P ,T , Pc, Tc, x , Ec )
    alf = 1 - (u[Ec] - 1)*B
    beta = A - u[Ec]*B +(v[Ec] - u[Ec])*np.power(B,2)
    gamma = A*B + v[Ec]*( np.power(B,2) + np.power(B,3) )
    if Fase == "V":
        z = np.roots([ 1 , -alf , beta , -gamma ])[ 0 ]
    elif Fase == "L":
        z = np.roots([ 1 , -alf , beta , -gamma ])[ 2 ]
    else:
        z = np.roots([ 1 , -alf , beta , -gamma ])
    return z

def z_1( Am , Bm , Fase , Ec = 0 ):
    Fase = str( Fase )
    u = np.array([ 2 , 1 , 0 ])
    v = np.array([ -1 , 0, 0 ])
    alf = 1 - (u[Ec] - 1)*Bm
    beta = Am - u[Ec]*Bm +(v[Ec] - u[Ec])*Bm**2
    gamma = Am*Bm + v[Ec]*( Bm**2 + Bm**3 )
    if Fase == "V":
        z = np.roots([ 1 , -alf , beta , -gamma ])[ 0 ]
    elif Fase == "L":
       z =  np.roots([ 1 , -alf , beta , -gamma ])[ 2 ]
    else:
        z = np.roots([ 1 , -alf , beta , -gamma ])
    return z

# def z_0( P , T , Pc , Tc , w , q , x , kij , Fase , Ec = 0 ):
#     Fase = str( Fase )
#     u = np.array([ 2 , 1 , 0 ])
#     v = np.array([ -1 , 0, 0 ])
#     A = Am_0( P , T , Pc , Tc , w , q , x , kij ,  Ec )
#     B = Bm_0( P ,T , Pc, Tc, x , Ec )
#     alf = 1 - (u[Ec] - 1)*B
#     beta = A - u[Ec]*B +(v[Ec] - u[Ec])*np.power(B,2)
#     gamma = A*B + v[Ec]*( np.power(B,2) + np.power(B,3) )
#     C = 3*beta - np.power( alf, 2 )
#     D = -np.power( alf, 3 ) + 4.5*alf*beta - 13.5*gamma
#     Q = np.power( C, 3 ) + np.power( D, 2 )
#     if Q > 0:
#         D1 = np.power( abs(-D + Q**0.5) )
#         z = (alf + ( ))

# def phi_i( P , T , Pc , Tc , w , q , x , kij , n , Fase , Ec = 0 ):
#     u = np.array([ 2 , 1 , 0 ])
#     v = np.array([ -1 , 0, 0 ])
#     Am = Am_0( P , T , Pc , Tc , w , q , x , kij , Ec )
#     Bi = B_i( P , T , Pc , Tc , Ec )
#     Bm =  np.dot( x , Bi )
#     A_prim = Aprim_i( P , T , Pc , Tc , w , q , x , kij , Ec )
#     z = z_1( Am , Bm , Fase , Ec )
#     # print(z)
#     delta = ( u[Ec]**2 - 4*v[Ec] )**0.5
#     L = np.log( ( 2*z + Bm*(u[Ec] + delta) )/( 2*z + Bm*( u[Ec] - delta) ) )/delta
#     ln_Phi = -np.log( z - Bm ) - ( z - 1 )*Bi[ n ]/Bm + ( Am/Bm )*( Bi[ n ]/Bm - A_prim[ n ]/ Am )*L
#     return np.exp( ln_Phi )

def phi_i( P , T , Pc , Tc , w , q , x , kij , n , Bi , Am , Bm , A_prim ,Fase , Ec = 0 ):
    u = np.array([ 2 , 1 , 0 ])
    v = np.array([ -1 , 0, 0 ])
    z = z_1( Am , Bm , Fase , Ec )
    # print(z)
    delta = np.power( np.power(u[Ec], 2) - 4*v[Ec] , 0.5 )
    L = np.log( ( 2*z + Bm*( u[Ec] + delta) )/( 2*z + Bm*( u[Ec] - delta ) ) )/delta
    ln_Phi = -np.log( z - Bm ) + ( z - 1 )*Bi[ n ]/Bm + ( Am/Bm )*( Bi[ n ]/Bm - A_prim[ n ]/ Am )*L
    return np.exp( ln_Phi )

# def phi_( P , T , Pc , Tc , w , q , x , kij , Fase , Ec = 0 ):
#     Phi = np.zeros( x.size )
#     for i in range( 0 , x.size ):
#         Phi[ i ] = phi_i( P , T , Pc , Tc , w , q , x , kij , i , Fase , Ec )
#     return Phi

def phi_( P , T , Pc , Tc , w , q , x , kij , Fase , Ec = 0 ):
    Ai = A_i( P , T , Pc , Tc , w , q , Ec )
    Am = Am_1( x , Ai, kij )
    Bi = B_i( P , T , Pc , Tc , Ec )
    Bm =  np.dot( x , Bi )
    A_prim = Aprim_i_1( x, Ai, kij )
    Phi = np.array([])
    for i in range( x.size ):
        Phi = np.append( Phi, phi_i( P , T , Pc , Tc , w , q , x , kij , i , Bi , Am , Bm, A_prim ,Fase , Ec ) )
    return Phi

"""****************************************************************************************************************************
Rutina para Presion de saturacion ( Psat_i )

Rutinas para los siguientes problemas de equilibrio fases ideales:
        a) Presion de Burbuja ( Pbur_id )
        b) Temperatura de rocio ( Troc_id )
****************************************************************************************************************************"""

def Psat_i( T , Pc , Tc , w ):
    log_10Pr = ( 7/3 )*( 1 + w )*( 1 - Tc/T )
    return Pc*np.power( 10, log_10Pr )

def Pbur_id( T , Pc , Tc , w , x ):
    return np.dot( Psat_i( T , Pc , Tc , w ) , x )

def Troc_id(  P , Pc , Tc , w , y , T = 350 , delta = 1e-8 , tol = 1e-6 , Nmax = 20 ):
    ferror = lambda T : np.log( P*sum( y/Psat_i( T , Pc , Tc , w ) ) )
    dferror = lambda T: ( ferror( T + delta ) - ferror( T ) )/( 1/(T +delta) - 1/T )
    N = 0
    while N < Nmax:
        N = N + 1
        if abs( ferror( T ) ) < tol:
            break
        T = 1/T - ferror( T )/dferror( T )
        T = 1/T
    # print("N. de iteraciones (Rutina Troc ideal):",N,"\n")
    return T

"""****************************************************************************************************************************
Ecuaciones Cubicas de Estado:
    0 - Peng Robinson
    1 - Soave Riedlich Kwon
    2 - Van der Waals

Rutinas para los siguientes problemas de equilibrio:
        a) Presion de Burbuja, genera estimados iniciales ( Pbur )
        b) Temperatura de rocio, genera estimados iniciales ( Troc )
****************************************************************************************************************************"""


def Pbur( T , Pc , Tc , w , q , z , kij , delta = 1e-10 , tol = 1e-6 , Nmax = 20 , Ec = 0 ):
    Psat = Psat_i( T , Pc , Tc , w )
    Pbur_r = np.dot( Psat , z )
    y_i = Psat*z/Pbur_r
    del( Psat )
    ferror = lambda P : sum( z*phi_( P , T , Pc , Tc , w , q , z , kij , "L" )/phi_( P , T , Pc , Tc , w , q , y_i , kij , "V" ) ) - 1
    dferror = lambda P : ( ferror( P + delta ) - ferror( P ) )/( 1/(P + delta) - 1/P )
    N = 0
    while N < Nmax:
        N = N + 1
        if abs( ferror( Pbur_r ) ) < tol:
            break
        Pbur_rr = 1/Pbur_r - ferror( Pbur_r )/dferror( Pbur_r )
        Pbur_r = 1/Pbur_rr
        y_i = z*phi_( Pbur_r , T , Pc , Tc , w , q , z , kij , "L" )/phi_( Pbur_r , T , Pc , Tc , w , q , y_i , kij , "V" )
        y_i = y_i/sum( y_i )
    # print( "N. de iteraciones (Rutina Pbur real):", N , "\n" )
    return  np.append( np.array([ Pbur_r ]) , y_i )

def Troc( P , Pc , Tc , w , q , y , kij , delta = 1e-8 , tol = 1e-6 , Nmax = 20 , Ec = 0 ):
    Troc_r = Troc_id(  P , Pc , Tc , w , y )
    x_i = P*y/Psat_i( Troc_r , Pc, Tc, w )
    ferror = lambda T: np.log( sum( y*phi_( P , T , Pc , Tc , w , q , y , kij , "V" )/phi_( P , T , Pc , Tc , w , q , x_i , kij , "L" ) ) )
    dferror = lambda T: ( ferror(T + delta) - ferror( T ))/( 1/(T + delta) - 1/T )
    N = 0
    while N < Nmax:
        N = N + 1
        if abs( ferror( Troc_r ) ) < tol:
            break
        Troc_r = 1/Troc_r - ferror( Troc_r )/dferror( Troc_r )
        Troc_r = 1/Troc_r
        x_i = y*phi_( P , T , Pc , Tc , w , q , y , kij , "V" )/phi_( P , T , Pc , Tc , w , q , x_i , kij , "L" )
        x_i = x_i/sum( x_i )
    # print( "N.de iteraciones (Rutina Troc real):",N,"\n")
    return np.append( np.array([ Troc_r ]) , x_i )

#Pruebas: nC3 y nC7 
if __name__=="__main__":
    P = 1.01013 #Bar
    T = 298.15 #K
    Pc = np.array([ 42.4953 , 27.4084 ])
    Tc = np.array([ 369.82 , 540.14 ])
    w = np.array([ 0.15416, 0.35 ])
    q = np.array([ -0.03136 , -0.02325 ])
    kij = np.array([ [0 , 0.0067],[ 0.0067 , 0] ])
    z = np.array([ 0.6 , 0.4 ])
    B = B_i( P , T , Pc , Tc )
    Bm = Bm_0( P ,T , Pc, Tc, z )
    A = A_i( P , T , Pc , Tc , w , q )
    alfa = alfa_SV( T , Tc , w , q )
    Am = Am_0( P , T , Pc , Tc , w , q , z , kij )
    A_prim = Aprim_i( P , T , Pc , Tc , w , q , z , kij )
    Z = z_1( Am, Bm, "Prueba")
    phi_nC3_L = phi_i( P , T , Pc , Tc , w , q , z , kij , 0 , B, Am, Bm , A_prim, "L" )
    phi_nC3_V = phi_i( P , T , Pc , Tc , w , q , z , kij , 0 , B, Am, Bm , A_prim, "V" )
    phi_nC7_L = phi_i( P , T , Pc , Tc , w , q , z , kij , 1 , B, Am, Bm , A_prim, "L" )
    phi_nC7_V = phi_i( P , T , Pc , Tc , w , q , z , kij , 1 , B, Am, Bm , A_prim, "V" )
    phi = phi_( P , T , Pc , Tc , w , q , z , kij , "L" )
    psat = Psat_i( T, Pc, Tc, w )
    Pb_id = Pbur_id( T , Pc , Tc , w , z )
    Pbur_r = Pbur( T , Pc , Tc , w , q , z , kij )
    Tr_id = Troc_id(  P , Pc , Tc , w , z )
    Troc_r = Troc( P , Pc , Tc , w , q , z , kij )
    print("R E S U L T A D O S   D E   P R U E B A S :")
    print("T:",T,"K  P:",P,"Bar \n")
    # print("___________________________________________________________________________________________","\n")
    # print("Bi:",B,"\n")
    # print("Bm:",Bm,"\n")
    # print("alfa_i:",alfa,"\n")
    # print("Ai:",A,"\n")
    # print("Am:",Am,"\n")
    # print("A_prim:",A_prim,"\n")
    # print("z:",z,"\n")
    print("phi_nC3 L:",phi_nC3_L,"\n")
    print("phi_nC3 V:",phi_nC3_V,"\n")
    print("phi_nC7 L:",phi_nC7_L,"\n")
    print("phi_nC7 V:",phi_nC7_V,"\n")
    # print("phi:",phi,"\n")
    print("psat:",psat,"Bar \n")
    print("Pbur_id:",Pb_id,"Bar \n")
    print("Pbur_real:",Pbur_r[0],"Bar \n")
    # print("yi_id:",psat*z/sum( psat*z ),"[nC3 , nC7] \n")
    # print("yi_real:",Pbur_r[ 1 : z.size + 1 ],"[nC3 , nC7] \n")
    print("Troc_id:",Tr_id,"K \n")
    print("Troc_real:",Troc_r[0],"K \n")
    # print("xi_id:",P*z/Psat_i( Tr_id , Pc, Tc, w ),"[nC3 , nC7] \n")
    # print("xi_real:",Troc_r[ 1 : z.size + 1 ],"[nC3 , nC7] \n")