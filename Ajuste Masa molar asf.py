# -*- coding: utf-8 -*-
"""****************************************************************************************************************************
Created on Sat Feb 27 11:38:23 2021

@author: damv_

Indonesia
w = np.array([0.232, 0.339, 0.382, 0.047, 0])
w_nC5 = 0.65
Ajuste de Masa molar = 2490.6117315383026 g/mol    en 4 iteraciones

Cold Lake
w = np.array([0.194, 0.381, 0.267, 0.155, 0.3/100])
w_nC7* = 0.515
Ajuste de Masa molar* = 3248.2131162411993 g/mol   en 4 iteraciones
w_nC7 = 0.490
Ajuste de Masa molar = 3359.0696375286734 g/mol   en 5 iteraciones
w_nC7 = 0.480 ,  n = 7  <-- Mejor ajuste con datos exp
Ajuste de Masa molar = 3431.6477618206495 g/mol   en 5 iteraciones

Venezuela 2 c:
w = np.array([0.205, 0.380, 0.196, 0.218, 0.1/100])
w_nC7 = 0.46
Ajuste de Masa molar = 3456.8069233264805 g/mol   en 5 iteraciones
w_nC7 = 0.44
Ajuste de Masa molar = 3560.7996045373425 g/mol   en 5 iteraciones

Lloydminster
w = np.array([0.231, 0.417, 0.195, 0.153, 0.4/100])
w_nC7* = 0.487
Ajuste de Masa molar* = 3316.0597650938685 g/mol   en 4 iteraciones
w_nC7 = 0.470
Ajuste de Masa molar = 3392.472440187809 g/mol en 5 iteraciones
w_nC7 = 0.46 <- Mejor ajuste 
Ajuste de Masa molar = 3482.818419746553 g/mol en 5 iteraciones


LLOYDMINSTER + COLD LAKE (PROPORCION DESCONOCIDA)
w = np.array([0.290, 0.422, 0.158, 0.130, 0.0])
w_nC7 = 0.4665 , n = 8 , Masf = 3287.2249486131836 g/mol

Golfo de Mx
w = np.array([0.503, 0.305, 0.146, 0.040, 0.006])
w_nC7 = 0.4188:
    n = 5, Masf = 2853.9742651039096
    n = 6, Masf = 2861.7467201202135 <- Bueno
    n = 7, Masf = 2866.95625873616
    n = 8, Masf = 2871.126870517777
    n = 10, Masf = 2877.5285555868336
****************************************************************************************************************************"""
import numpy as np
import Equilibrio.equilibrio as eq

#DATOS
T = 23 + 273.15 #K
P = 1.01325 #Bar
n = 9 #num. de pseudocomponentes en que se parte el asfalteno
#w_nC7 es el punto de floculacion
w_nC7 = 0.46
#w = np.array([Saturado,Aromatico,Resina,Asfalteno,Solidos])
w = np.array([0.231, 0.417, 0.195, 0.153, 0.4/100]) #Fraccion masa
#Eliminar fraccion de solidos inorganicos (no participan en el equilibrio)
w = eq.Recomponer( w ) #Eliminar fraccion de solidos y reajuste de fraccion masa
M = np.array([460.0, 521.86, 1042.86, 0]) #Masa molar g/mol
#Parametros para partir la fraccion de asfaltenos con distribucion gamma NO MODIFICAR
Mmin_asf = 1800 #g/mol / no modificar
beta = 3.5 #Parametro de forma distribucion Gamma / no modificar

"""****************************************************************************************************************************
Corregir las fracciones masa, eliminando los solidos y reajustando las
fracciones. Posteriormente se calculan las fracciones mol de acuerdo a la
composicion masica. Generar las composiciones molares

Solvente: nC7 algunos casos (Indonesia) se usa nC5
Se añade al final de las listas np.array 

###############################################################################
Propiedades criticas de fracciones SAR promedio y solvente utilizado nC5

PROPIEDAD = np.array([SATURADO, AROMATICO, RESINA, SOLVENTE])
El parametro de solubilidad se calcula ya que depende de la temperatura. 
De igual forma, los volumenes molares se calculan "in situ", por su dependencia con T.

Propiedades nC5: M = 72.15 g/mol, Pc = 33.70 Bar, Tc = 469.70 K, V* = 314.8268967 cm3/mol, w = 0.252 , d = 22.222 - 0.0264*T MPa**0.5
Propiedades nC7: M = 100.204 g/mol, Pc = 27.40 Bar, Tc = 540 K, V* = 434.4297835197149 cm3/mol, w = 0.35 , d = 22.121 - 0.0232*T MPa**0.5
Propiedades de solventes obtenidas de B. E. Poling, J. M. Prausnitz y J. P. O'Connell, The Properties of Gases and Liquids, quinta ed., Nueva York: McGraw-Hill, 2000
****************************************************************************************************************************"""
#SOLVENTE: nC5
# M = np.append( M, 72.15 ) #g/mol
# Pc = np.array([5.91809024, 12.86868178, 7.40880462, 33.70]) #Bar
# Tc = np.array([856.50785811, 1481.55069013, 2052.45144954, 469.70]) #K
# omega = np.array([1.36549975, 0.90174793, 1.20046341, 0.252]) #Factor acentrico
# va = np.array([2420.4837844,   2233.96976857,  4995.19050729, 314.8268967]) #Parametro correlacions Costald
# d = np.append( eq.del_SARA( T , M[0 : 3] ) , 22.222 - 0.0264*T )
# vsat = eq.Vs_Costald( T , va , omega , Tc )

# #SOLVENTE: nC7
M = np.append( M, 100.204 ) #g/mol
Pc = np.array([5.91809024, 12.86868178, 7.40880462, 27.40]) #Bar
Tc = np.array([856.50785811, 1481.55069013, 2052.45144954, 540.0]) #K
omega = np.array([1.36549975, 0.90174793, 1.20046341, 0.35]) #Factor acentrico
va = np.array([2420.4837844,   2233.96976857,  4995.19050729, 434.4297835197149]) #Parametro correlacions Costald
d = np.append( eq.del_SARA( T , M[0 : 3] ) , 22.121 - 0.0232*T )
vsat = eq.Vs_Costald( T , va , omega , Tc )

#La mejor ferror es:
def ferror( z, ki ): #Nombrada funcion error 2 en tesis
    x_h = z[2 : n + 3]*ki
    return np.log( sum( x_h ) )

#Calcula los parametros para asfaltenos dando una masa molar y 
def Rec_Eq( Mprom ):
    M_ = M.copy()
    M_[ 3 ] = Mprom
    w_ = w.copy()
    w_ = np.append( ( 1 - w_nC7 )*w_ , w_nC7 )
    z = ( w_/M_ )/sum( w_/M_ )
    zasf = z[ -2 ]
    z = np.delete( z , 3 )
    z = np.insert( z , 3 , eq.PseudoComp( zasf , Mprom , n , Mmin_asf , beta)[ 0 ] )
    M_ = np.insert( M_ , 3 , eq.PseudoComp( zasf , Mprom , n , Mmin_asf , beta)[ 1 ] )
    M_ = np.append( M_ , M[ -1] )
    #Las correlaciones son funcion de M, es necesario calcular todas las propiedades necesarias
    Pc_ = Pc.copy()
    Pc_ = np.insert( Pc_, 3, eq.Pc_Asf( M_ , n ) )
    Tc_ = Tc.copy()
    Tc_ = np.insert( Tc_, 3, eq.Tc_Asf( M_ , n ) )
    omega_ = omega.copy()
    omega_ = np.insert( omega_, 3, eq.Omega_Asf( M_ , n ) )
    va_ = va.copy()
    va_ = np.insert( va_, 3, eq.V_a_asf( Tc_ , Pc_ , omega_ , n ) )
    vsat_ = eq.Vs_Costald( T , va_ , omega_ , Tc_ )
    vsat_ = eq.V_Comp( P , T , vsat_ , Pc_ , Tc_ , omega_ )
    d_ = d.copy()
    d_ = np.insert( d_, 3, eq.del_asf( T , M_ , n ) )
    xl = z.copy()
    vml_ = eq.Vm_Costald( T, va_, Tc_, omega_, xl )
    # Estimado inicial 8% resinas y 92% asfaltenos
    xh_ = np.append( 0.08 , 0.92*z[ 3 : n + 3 ]/sum( z[ 3 : n + 3 ] ) ) 
    xh_ = xh_/sum( xh_ )
    vmh_ = eq.Vm_Costald( T, va_[ 2 : n + 3 ], Tc_[ 2 : n + 3 ] , omega_[ 2 : n + 3 ] , xh_ )
    # print("vmh:",vmh,"\n")
    #Se corrigen las composiciones en 3 iteraciones. En el estudio de sensibilidad de estimados
    #iniciales se encontro que la convergencia se cumple en la mayoria de los sistemas en dos iteraciones
    xh_ = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_) )
    vmh_ = eq.Vm_Costald( T, va_[ 2 : n + 3 ], Tc_[ 2 : n + 3 ] , omega_[ 2 : n + 3 ] , xh_ )
    xh_ = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_) )
    vmh_ = eq.Vm_Costald( T, va_[ 2 : n + 3 ], Tc_[ 2 : n + 3 ] , omega_[ 2 : n + 3 ] , xh_ )
    xh_ = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_) )
    vmh_ = eq.Vm_Costald( T, va_[ 2 : n + 3 ], Tc_[ 2 : n + 3 ] , omega_[ 2 : n + 3 ] , xh_ )
    # print("Error:",ferror( xl, eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_) ),"\t")
    return  ferror( xl, eq.K_ihl( T, vsat_ , vml_, vmh_, xl, xh_, d_) )

#Aproximacion de la derivada de la funcion error, se usa metodo de la secante
def dferror_dM( Mprom , delta = 1e-3 ):
    f_1 = Rec_Eq( Mprom + delta )
    f_2 = Rec_Eq( Mprom )
    return ( f_1 - f_2 )/delta

#Estimado inicial y tolerancia
M_ajustada = 2865
tol = 1e-8

print("#### Ajuste final de Masa molar de asfaltenos ####\n")
for i in range( 10 ):
    if abs( Rec_Eq( M_ajustada ) ) < tol:
        # print("################# CONVERGENCIA ALCANZADA ################# \n")
        print( "Funcion Error < " , tol, ", Iteraciones Necesarias:",i )
        # print("ERROR:", Rec_Eq( M_ajustada ) ,"\n")
        break
    # print("############### i =",i + 1,"###############")
    M_ajustada = M_ajustada - Rec_Eq( M_ajustada )/dferror_dM( M_ajustada )

print("\t RESULTADOS FINALES:")
M[ 3 ] = M_ajustada
print("\t Masf ajustada:",M_ajustada,"g/mol \t")
print("Fraccion solvente:",w_nC7,"\t N. de particiones asfalteno:",n)
print("############################## FIN ##############################")