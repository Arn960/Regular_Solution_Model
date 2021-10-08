# -*- coding: utf-8 -*-
"""****************************************************************************************************************************
Created on Wed Feb 10 18:59:20 2021

@author: damv_

El crudo de Indonesia es el unico que utiliza nC5 como solvente en la prueba de floculacion.

Indonesia
w = np.array([0.232, 0.339, 0.382, 0.047, 0])
w_nC7 = 0.65 #El solvente usado, realmente es nC5
M_asf = 2487.1212963954586 #g/mol

Venezuela 2
w = np.array([0.205, 0.380, 0.196, 0.218, 0.1])
w_nC7 = 0.46
Masf = 3456.8069 #g/mol
****************************************************************************************************************************"""
import numpy as np
import Equilibrio.equilibrio as eq
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import PercentFormatter

#DATOS
T = 23 + 273.15 #K
P = 1.01325 #Bar
n = 5 #num. de pseudocomponentes en que se parte el asfalteno
w_nC7 = 0.65

"""###########[ Sat ,Arom ,Resin ,Asfal ,Solid2 ]###########################################################################"""
w = np.array([0.232, 0.339, 0.382, 0.047, 0]) #Fraccion masa
M = np.array([460.0, 521.86, 1042.86, 2487.1212963954586 ]) #Masa molar g/mol

"""****************************************************************************************************************************
Corregir las fracciones masa, eliminando los solidos y reajustando las
fracciones. Posteriormente se calculan las fracciones mol de acuerdo a la
composicion masica. Generar las composiciones molares

Solvente: nC7 algunos casos (Indonesia) se usa nC5
Se añade al final de las listas np.array 
****************************************************************************************************************************"""
w = eq.Recomponer( w ) #Elimina solidos
w = np.append( ( 1 - w_nC7 )*w , w_nC7 ) #Añade el solvente
#Caso nC5
M = np.append( M, 72.15 ) #g/mol
#Caso nC7
# M = np.append( M, 100.204 )

z = ( w/M )/sum( w/M ) #Pasamos a fraccion mol
# del( w_nC7 )
del( w ) #liberamos ram
beta = 3.5
"""****************************************************************************************************************************
Generación de pseudocomponentes, bajo los parametros:
    M prom: M de Asfaltenos.
    M min*: 1800 g/mol, monomero de asfalteno.
    alfa*: 3.5 
    n: 5
* K.Akbarzadeh et al./ Fluid Phase Equilibria 232 (2005)
Los valores de alfa y Mmin pueden modificarse, por defecto se usan los valores
ya mencionados.
****************************************************************************************************************************"""
Mprom = M[ -2 ] #Masa molar promedio de asfaltenos
zasf = z[ -2 ]  #Fracción de asfaltenos
z = np.delete( z , 3 )
M = np.delete( M , 3 )
z = np.insert( z , 3 , eq.PseudoComp( zasf , Mprom , n , 1800 , beta)[ 0 ] )
M = np.insert( M , 3 , eq.PseudoComp( zasf , Mprom , n , 1800 , beta)[ 1 ] )

"""****************************************************************************************************************************
Propiedades criticas y parametros para fracciones SARA
Propiedades criticas: Pc, Tc y Vc (esta última para matriz Kij)
Paramertos: Factor acéntrico y V* para correlacion HBT (volumen molar)
Propiedades nC5: M = 72.15 g/mol, Pc = 33.70 Bar, Tc = 469.70 K, Vc = 331 cm3/mol, w = 0.252, d = 22.222 - 0.0264*T MPa**0.5
Propiedades nC7: M = 100.204 g/mol, Pc = 27.40 Bar, Tc = 540 K, Vc = 428 cm3/mol, w = 0.35, d = 22.121 - 0.0232*T MPa**0.5
****************************************************************************************************************************"""
#nC5 Datos de disolvente
Pc = np.append( eq.Pc_SARA( M[0 : n + 3] ) , 33.70  ) #Bar
Tc = np.append( eq.Tc_SARA( M[0 : n + 3] ) , 469.70 ) #K
omega = np.append( eq.Omega_SARA( M[0 : n + 3] ) , 0.252 )#factor acentrico
Va = eq.V_a( Tc[0 : n + 3] , Pc[0 : n + 3] , omega[0 : n + 3] ) #cm3/mol parametro correlacion HBT
Va = np.append( Va , eq.V_a_solv( Tc[ -1 ] , Pc[ -1 ] , M[ -1 ] , omega[ -1 ] ) ) #cm3/mol nC5
d = np.append( eq.del_SARA( T , M[0 : n + 3] ) , 22.222 - 0.0264*T )#MPa**0.5

#nC7 Datos de disolvente
# Pc = np.append( eq.Pc_SARA( M[0 : n + 3] ) , 27.40  ) #Bar
# Tc = np.append( eq.Tc_SARA( M[0 : n + 3] ) , 540.0 ) #K
# omega = np.append( eq.Omega_SARA( M[0 : n + 3] ) , 0.35 )#factor acentrico
# Va = eq.V_a( Tc[0 : n + 3] , Pc[0 : n + 3] , omega[0 : n + 3] ) #cm3/mol parametro correlacion HBT
# Va = np.append( Va , eq.V_a_solv( Tc[ -1 ] , Pc[ -1 ] , M[ -1 ] , omega[ -1 ] ) ) #cm3/mol nC7
# d = np.append( eq.del_SARA( T , M[0 : n + 3] ) , 22.121 - 0.0232*T )#MPa**0.5

"""****************************************************************************************************************************
Volumen molar de líquido saturado, he observado que la diferencia entre el 
liquido subenfriado y saturado son minimas para las fracciones SARA.
vsat [=] cm^3/mol                densidad = M/v

Ademas se generan estimados iniciales.
xh = fraccion en fase pesada (solo hay resinas y asfaltenos)
xl = fraccion en fase ligera (se encuentran todos los componentes)
El ultimo elemento de xl corresponde al nC7 añadido

Contador composiciones de Resina y todos los asfaltenos: [2 : n + 3]
Contador todas las composiciones SARA: [0 : n + 3]
Contador para nC7: [ -1 ] y [ n + 3 ]
Contador Resinas: [ 2 ]
ESTIMADOS INICIALES 
****************************************************************************************************************************"""
vsat = eq.Vs_Costald( T , Va , omega , Tc )
#V_comp traslada a condiciones de L.subenfriado 
vsat = eq.V_Comp( P , T , vsat , Pc , Tc , omega )

# Estimados iniciales
xh = np.random.random( n + 1 )
# xh = np.random.random(n)
# xh = np.append(0.1109 , 0.8891*xh/sum( xh ) )
xh = xh/sum( xh )
xl = z.copy(  )

vml = eq.Vm_Costald( T, Va, Tc, omega, xl )
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )
xh_estimados = xh.copy()

"""
#Recalcular xh con metodo de Michelsen
xh = eq.Yi_Mich( T ,vsat , vmh, xh , d , z )/sum( eq.Yi_Mich( T ,vsat , vmh, xh , d , z ) )
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )

"""
# print("___________________________________________________________________________________________","\n")

#Recalcular con xh_i = K_i*xl_i
xh = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d) )
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )

xh_inter = xh.copy()

xh = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d) )
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )


"""****************************************************************************************************************************
Calculo de equilibrio:
Ki = xi_h / xi_l
Ferr = sum( zi*(Ki - 1) ) -> 0 
dFerr/dxnC7 se aproxima con diferenciales, delta = 1e-7
tol = 1e-6

def recal_xh( z ):
    return z[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, z, xh, d)/sum( z[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, z, xh, d) )
****************************************************************************************************************************"""

def ferror( z, ki ):
    f_error = np.zeros( z.size )
    f_error[0 : 2] = z[0 : 2]
    f_error[2 : n + 3] = z[2 : n + 3]*(1 - ki)
    f_error[ -1 ] = z[ -1 ]
    return sum( f_error )

ki = eq.K_ihl( T, vsat , vml, vmh, xl, xh, d )
xh_recalculados = xh.copy()
f_error = np.array([ ferror( z, ki ) ])
xl = z.copy(  )

print("########################### INICIO MONTECARLO ###########################")
N = 5000
for i in range( N ):
    xh = np.random.random( n + 1 )
    xh = xh/sum( xh )
    xh_estimados = np.vstack(( xh_estimados, xh ))
    #Correccion de xh
    #Correccion, ciclo 1
    vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )
    xh = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d) )
    xh_inter = np.vstack((xh_inter, xh))
    #Correccion, ciclo 2
    vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )
    xh = xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d)/sum( xl[ 2 : n + 3 ]*eq.K_ihl( T, vsat , vml, vmh, xl, xh, d) )
    #Guardar recalculados
    xh_recalculados = np.vstack(( xh_recalculados, xh ))
    vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh )
    ki = eq.K_ihl( T, vsat , vml, vmh, xl, xh, d )
    f_error = np.append( f_error , ferror( z, ki ) )

print("\t RESULTADOS \n")
ferror_min = np.amin( abs( f_error ) )
prom = np.mean( abs(f_error) )
xh_r_prom = np.mean( xh_recalculados , axis = 0 )
xh_inter_prom = np.mean( xh_inter , axis = 0 )
minimos = np.where( abs(f_error) == np.amin( abs(f_error) ) )
print("Fraccion masa de solvente:",w_nC7,"\n")
print("N. de pseudo componenetes de asfaltenos:",n,"\n")
print("Ferror Min:",ferror_min,"\n")
print("Ferror Max:",np.amax( abs(f_error) ),"\n")
print("% dif min y promedio:",100*abs( ( prom - ferror_min)/ferror_min ),"\n")
print("Varianza del Error:",np.var( f_error ),"\n")
print("Ubicacion del minimo:",minimos,"\n")
print("composicion estimada en el minimo:",xh_estimados[int( minimos[0] )],"\n")
print("Asfaltenos en fase pesada (xh en el min):",1 - xh_estimados[minimos[0],0],"\n")
print("xh en mín:",xh_recalculados[minimos[0]],"\n")
print("##### COMPOSICION FASE PESADA CORREGIDA ########## RESULTADOS ##### \n")
print("##### UNA ITERACION #### \n")
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh_inter_prom )
ki = eq.K_ihl( T, vsat , vml, vmh, xl, xh_inter_prom, d )
#Resultados de una y dos iteraciones
print("xh promedios:",xh_inter_prom,"\n ERROR xhprom:",ferror( z, ki ),"\n")
print("xh varianza:",np.var( xh_inter , axis = 0 ),"\n")
print("xh D.Est.:",np.power(np.var( xh_inter , axis = 0 ),0.5),"\n")
vmh = eq.Vm_Costald( T, Va[ 2 : n + 3 ], Tc[ 2 : n + 3 ] , omega[ 2 : n + 3 ] , xh_r_prom )
ki = eq.K_ihl( T, vsat , vml, vmh, xl, xh_r_prom, d )
print("##### DOS ITERACIONES #### \n")
print("xh promedios:",xh_r_prom,"\n ERROR xhprom:",ferror( z, ki ),"\n")
print("xh varianza:",np.var( xh_recalculados , axis = 0 ),"\n")
print("xh D.Est.:",np.power(np.var( xh_recalculados , axis = 0 ),0.5),"\n")

print("############################## FIN ##############################")

#Extraer figura y ejes
#gfiguras acomodadas en 2 filas de 3 columnas
fig, ( (ax1, ax2, ax3 ), (ax4, ax5, ax6) ) = plt.subplots(nrows=2, ncols=3)
fig.subplots_adjust(hspace=0.3) #Separacion entre figuras
plt.rc('axes', titlesize=20)

#Histogramas (Ejes de la figura general)
N_bins = 13
font_size = 16
cod_color_1 = mpatches.Patch(color="tab:blue", label="Una iteración")
cod_color_2 = mpatches.Patch(color="tab:orange", label="Dos iteraciones")
fig.legend(handles=[cod_color_1,cod_color_2], loc=1, fontsize="xx-large")

#Eje 1: Resinas
ax1.hist(xh_inter[:,0], bins=N_bins) #Datos y num. de barras. 1 iter
ax1.hist(xh_recalculados[:,0], bins=N_bins) #Datos y num. de barras. 2 iter
ax1.set_title("Resinas", fontsize=18) #Titulo
ax1.set_xlabel("Fracción mol",fontsize = font_size) #Etiquetas eje x
ax1.set_ylabel("Frecuencia",fontsize = font_size) #Etiquetas eje y
ax1.yaxis.set_major_formatter(PercentFormatter(xmax=N)) #Cambiar a frecuencia relativa
ax1.tick_params(axis='x', labelsize=12) #Tamaño de numeros en eje x
ax1.tick_params(axis='y', labelsize=12) #Tamaño de numeros en eje y

#Eje 2: Asfalenos 1
ax2.hist(xh_inter[:,1], bins=N_bins)
ax2.hist(xh_recalculados[:,1], bins=N_bins)
ax2.set_title("Asfaltenos 1", fontsize=18)
ax2.set_xlabel("Fracción mol",fontsize = font_size)
ax2.yaxis.set_major_formatter(PercentFormatter(xmax=N))
ax2.tick_params(axis='x', labelsize=12)
ax2.tick_params(axis='y', labelsize=12)

#Eje 3: Asfaltenos 2
ax3.hist(xh_inter[:,2], bins=N_bins)
ax3.hist(xh_recalculados[:,2], bins=N_bins)
ax3.set_title("Asfaltenos 2", fontsize=18)
ax3.set_xlabel("Fracción mol",fontsize = font_size)
ax3.yaxis.set_major_formatter(PercentFormatter(xmax=N))
ax3.tick_params(axis='x', labelsize=12)
ax3.tick_params(axis='y', labelsize=12)

#Eje 4: Asfaltenos 3
ax4.hist(xh_inter[:,3], bins=N_bins)
ax4.hist(xh_recalculados[:,3], bins=N_bins)
ax4.set_title("Asfaltenos 3", fontsize=18)
ax4.set_xlabel("Fracción mol",fontsize = font_size)
ax4.set_ylabel("Frecuencia",fontsize = font_size)
ax4.yaxis.set_major_formatter(PercentFormatter(xmax=N))
ax4.tick_params(axis='x', labelsize=12)
ax4.tick_params(axis='y', labelsize=12)

#Eje 5: Asfaltenos 4
ax5.hist(xh_inter[:,4], bins=N_bins)
ax5.hist(xh_recalculados[:,4], bins=N_bins)
ax5.set_title("Asfaltenos 4", fontsize=18)
ax5.set_xlabel("Fracción mol",fontsize = font_size)
ax5.yaxis.set_major_formatter(PercentFormatter(xmax=N))
ax5.tick_params(axis='x', labelsize=12)
ax5.tick_params(axis='y', labelsize=12)

#Eje 6: Asfaltenos 5
ax6.hist(xh_inter[:,5], bins=N_bins)
ax6.hist(xh_recalculados[:,5], bins=N_bins)
ax6.set_title("Asfaltenos 5", fontsize=18)
ax6.set_xlabel("Fracción mol",fontsize = font_size)
ax6.yaxis.set_major_formatter(PercentFormatter(xmax=N))
ax6.tick_params(axis='x', labelsize=12)
ax6.tick_params(axis='y', labelsize=12)

"""
#Diagramas violin
#Eje 1: Resinas
ax1.violinplot(xh_recalculados[:,0], showmedians=True)
ax1.set_xticks([1])
# ax1.set_xticklabels([r"$x_{resina}$"], fontsize=18)
ax1.set_title(r"$x_{resina}$")
ax1.set_ylabel('Fracción mol')

#Eje 2: Asfalteno partición 1
ax2.violinplot(xh_recalculados[:,1], showmedians=True)
ax2.set_xticks([1])
# ax2.set_xticklabels([r"$x_{asfalteno 1}$"], fontsize=18)
ax2.set_title(r"$x_{asfalteno 1}$")

#Eje 3: Asfalteno partición 2

ax3.violinplot(xh_recalculados[:,2], showmedians=True)
ax3.set_xticks([1])
# ax3.set_xticklabels([r"$x_{asfalteno 2}$"], fontsize=18)
ax3.set_title(r"$x_{asfalteno 2}$")

#Eje 4: Asfalteno partición 3
ax4.violinplot(xh_recalculados[:,3], showmedians=True)
ax4.set_xticks([1])
# ax4.set_xticklabels([r"$x_{asfalteno 3}$"], fontsize=18)
ax4.set_title(r"$x_{asfalteno 3}$")
ax4.set_ylabel('Fracción mol')

#Eje 5: Asfalteno partición 4
ax5.violinplot(xh_recalculados[:,4], showmedians=True)
ax5.set_xticks([1])
ax5.set_xticklabels([r"$x_{asfalteno 4}$"], fontsize=18)
# ax5.set_title(r"$x_{asfalteno 4}$")

#Eje 6: Asfalteno partición 5
ax6.violinplot(xh_recalculados[:,5], showmedians=True)
ax6.set_xticks([1])
ax6.set_xticklabels([r"$x_{asfalteno 5}$"], fontsize=18)
# ax6.set_title(r"$x_{asfalteno 5}$")

fig.suptitle("Distribución de composiciones recalculadas a partir de composiciones aleatorias", fontsize=25)
plt.show()
"""

#PROM recalculando con Kij
#PROM n = 4 , beta = 3.5 - 0.24725672857617198
#PROM n = 4 , beta = 3.7 - 0.23999448012870475
#PROM n = 4 , beta = 3.9 - 0.2342588848166447
#PROM n = 4 , beta = 4.5 - 0.22268288100759276
#PROM n = 4 , beta = 5.4 - 0.21293248519866510
#PROM n = 4 , beta = 6.0 - 0.2087131398056176
#PROM n = 4 , beta = 8.0 - 0.20054034435944137

#Datos segunda fase: Estimado con ferror minimo
# z_asf = 0.69448043 ,  ferror = 1.4271320791792164e-09
# z_asf = 0.8529445 ,  ferror = 3.9659708850336983e-10
# z_asf = 0.87119182 ,  ferror = 2.2906854191262482e-10
# z_asf = 0.73124096 ,  ferror = 1.0634586544711055e-09
# z_asf = 0.88051667 ,  ferror = 4.210317317010492e-09 
# z_asf = 0.81551756 ,  ferror = 2.8330274792764953e-09
# z_asf = 0.98378219 ,  ferror = 2.145998601044141e-09
# z_asf = 0.82764264 ,  ferror = 5.020034388181216e-10 
# z_asf = 0.82764264 ,  ferror = 5.291876936652784e-10 

