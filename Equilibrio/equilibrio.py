"""****************************************************************************
@author: damv_         Hecho con amor, no olvides eso
Equilibrio para proyecto Calculos Equilibrio
Funciones y algoritmos necesarios para calculo de equilibrio. Libreria personal.

VALORES CONSTANTE R
R = 83.14472 Bar cm3 / mol K
R = 8.314472 MPa cm3 / mol K

****************************************************************************"""
from scipy.special import gamma
import numpy as np

"""****************************************************************************
Distribución gamma, alfa para sistemas polimericos = 3.5, probado y aprobado
zpseud = fracción de interés, Mprom = M.molar de fracción, Mmin = Masa molar mínima,
n = número de pseudo componentes, alfa = Parámetro de forma para distribución Gamma
****************************************************************************"""
#PseudoComp parte un pseudocomponente en varias fracciones siguiendo el metodo de la funcion
#Gamma reportada por Whitson. Mmin y alfa fijadas estan fijas para asfaltenos. 
def PseudoComp( zpseud , Mprom , n , Mmin = 1800 , alfa = 3.5 ): #Entradas escalares
    z = np.array([])
    M = np.array([])
    X,W = np.polynomial.laguerre.laggauss( n ) #Puntos y pesos para Gauss Laguerre
    beta = ( 2.5*Mprom-Mmin )/X[ -1 ]
    exs = alfa*beta/(Mprom-Mmin)
    delta = np.exp(exs-1)
    f = lambda j : ( X[ j ]**( alfa-1 ) )*( exs**alfa )/( gamma( alfa )*delta**X[ j ] )
    for i in range( n ):
        zi = zpseud*( W[ i ]*f( i ) )
        Mi = Mmin+beta*X[ i ]
        z = np.append( z , zi )
        M = np.append( M , Mi )
    cz = zpseud/sum( z ) #Correccion de z
    z = z*cz
    cM = Mprom*sum( z )/sum( M*z ) #Correccion de M
    M = M*cM
    return np.array([ z , M ])

#Recomponer, en caso de que se reporte la fraccion de solidos (arcilla, arena, etc.)
def Recomponer(w): #Entrada np.array
    if w.size > 4 :
        w = w[ :4 ]
        w = w/sum(w)
    return w

# Blend: FUSIONA COMPONENTES SAR, se usa con propiedades promedio preferible para el
# modelo de solucion regular "SR"
def Blend_SR(w_1,w_2,Porcentaje_masa_crudo_1): #w_1 y w_2 son np.array, el % es escalar
    w = np.zeros( 3 )
    for i in range(0,3):
        w[i] = Porcentaje_masa_crudo_1*w_1[i] + (1 - Porcentaje_masa_crudo_1)*w_2[i]
    w_1 = np.delete(np.flip(w_1,0),[-1,-2,-3])
    w_2 = np.delete(np.flip(w_2,0),[-1,-2,-3])
    w = np.append(np.append(w,w_1),w_2)
    return w / sum( w )

#NO USADO
def Blend_Oil( w_1, w_2, Fraccion_masa_crudo_1 ):
    w = Fraccion_masa_crudo_1*w_1[0 : 3] + (1 - Fraccion_masa_crudo_1)*w_2[0 : 3]
    w = np.append( w , Fraccion_masa_crudo_1*w_1[ -1 ] )
    w = np.append( w , (1 - Fraccion_masa_crudo_1)*w_2[ -1 ] )
    return w


# Blend_VL: NO FUSIONA COMPONENTES SAR

def Blend_VL(w_1,w_2,Porcentaje_masa_crudo_1): #w_1 y w_2 son np.array, el % es escalar
    w_1 = Porcentaje_masa_crudo_1*w_1
    w_2 = (1 - Porcentaje_masa_crudo_1)*w_2
    w = np.append( w_1 , w_2 )
    return w/sum(w)

"""
#Método de prueba para saturados, queda pendiente de probar

def Propcrit_S(M): #Saturados y aromaticos, método de 
    g=-0.567944+(M-66)**0.0590526 #gravedad especifica, Soreide
    Tb=1928.3-1.695e5*(M**-0.03522)*(g**3.266)*np.exp(-4.922e-3*M-4.7685*g+3.462e-3*M*g) #°R, Soreide
    Tc=341.7+811*g+(0.4244+0.1174*g)*Tb+(0.4669-3.2623*g)*1e5/Tb #°R, Kesler - Lee
    Pc=8.3634-0.0566/g-(0.24244+2.2898/g+0.11857/(g**2))*1e-3*Tb+(1.4685+3.648/g+0.47227/(g**2))*1e-7*Tb**2-(0.42019+1.6977/(g**2))*1e-10*Tb**3 #lnPc, Kesler - Lee
    Pc=np.exp(Pc) #Psia
    Kw=(Tb**(1/3))/g #Factor K de Watson
    A=[-5.92714,6.09648,1.28862,-0.16934,15.2518,-15.6875,-13.4721,0.43577] #Parametros de correlacion de Lee - Kesler
    Tbr=Tb/Tc
    if Tbr <0.8: #Correlacion de Lee Kesler para el factor acentrico
        w=(-np.log(Pc/14.7)+A[0]+A[1]/Tbr+A[2]*np.log(Tbr)+A[3]*Tbr**6)/(A[4]+A[5]/Tbr+A[6]*np.log(Tbr)+A[7]*Tbr**6)
    else:
        w=-7.904+0.1352*Kw-0.007465*Kw**2+8.359*Tbr+(1.408-0.01063*Kw)/Tbr
    Pc=Pc/14.5038 #Bar
    Tc=((Tc-491.67)/1.8)+273.15 #Kelvin
    return{"Pc":Pc,"Tc":Tc,"w":w}
"""

#Propiedades criticas
#presion critica SARA
def Pc_SARA(M): #Bar
    m = M.size
    Pc = np.zeros(m)
    #Saturados, Riazi
    Pc[0] = np.exp(4.65757-0.13426*M[0]**0.5)
    #Aromaticos
    Pc[1] = 1891.4*M[1]**-0.7975
    #Resinas
    Pc[2] = 1891.4*M[2]**-0.7975
    for i in range(3,m):
        Pc[i] = 1891.4*M[i]**-0.7975
    return Pc

#Presion critica 
def Pc_Asf( M , n ): #Usar cuando se incluye nC-7
    Pc = np.zeros( n )
    for i in range( 0 , n  ):
        Pc[ i ] = 1891.4*M[ i + 3 ]**-0.7975
    return Pc

#Temperatura critica SARA
def Tc_SARA(M): #Kelvin
    m = M.size
    Tc = np.zeros(m)
    #Saturados Riazi
    Tb = 1070 - np.exp( 6.98291 - 0.02013*M[0]**(2/3) )
    Tbr = 1.15 - np.exp( -0.41966 - 0.02436*M[0]**0.58)
    Tc[0] = Tb/Tbr
    #Aromaticos
    Tc[1] = 77.856*M[1]**0.4708
    #Resinas
    Tc[2] = 77.856*M[2]**0.4708
    #Asfaltenos
    for i in range(3,m):
        Tc[i] = 77.856*M[i]**0.4708
    return Tc

#Temperatura critica asfaltenos 
def Tc_Asf( M , n ): #Usar cuando se incluye nC-7
    Tc = np.zeros( n )
    for i in range( 0 , n  ):
        Tc[ i ] = 77.856*M[ i + 3 ]**0.4708
    return Tc

#Volumen molar de fracciones SARA - NO USADO
def Vc_SARA(M): #cm3/mol
    m = M.size
    Vc = np.zeros(m)
    #Saturados
    densidad_crit = 0.26 - np.exp( -3.50532 - 1.5e-6*M[0]**2.38)#g/cm3
    Vc[0] = M[0]/densidad_crit
    #Aromaticos
    Vc[1] = 2.4988*M[1]+116.8879
    #Resinas
    Vc[2] = 2.4988*M[2]+116.8879
    #Asfaltenos
    for i in range(3,m):
        Vc[i] = 2.4988*M[i]+116.8879
    return(Vc)

#Factor acentrico "w" para fracciones SARA
def Omega_SARA( M ):
    m = M.size
    omega = np.zeros(m)
    cf = np.array([ 0.8098 , 0.7910 , 0.7940 ])
    om = lambda M,n: cf[ n ]*(0.5837*np.log( M )-2.5389)
    #Saturados
    omega[0] = np.exp( -3.06826 + 1.04987*M[ 0 ]**0.2 ) - 0.3
    #Aromaticos
    omega[1] = om( M[ 1 ] , 0 )
    #Resinas
    omega[2] = om( M[ 2 ] , 1 )
    #Asfaltenos
    for i in range( 3 , m ):
        omega[i] = om( M[ i ] , 2 )
    return omega

#Factor acentrico "w"
def Omega_Asf( M , n ):
    w = np.zeros( n )
    for i in range( 0 , n ):
        w[ i ] = 0.7940*( 0.5837*np.log( M[ i + 3 ] ) - 2.5389 )
    return w
#Matriz de parametros de interaccion binaria Kij - NO USADO

def Kij(v): #usar con np.array
    n=v.size
    K=np.zeros((n,n))
    for i in range(0,n-1):
        for j in range(i+1,n):
            K[i,j]=1-8*((v[i]*v[j])**.5)/((v[i]**(1/3))+(v[j]**(1/3)))**3
            K[j,i]=K[i,j]
    return K

"""
Correlaciones para volumen molar, usar np.array para Tc, Pc,vc, w y fracciones x
"""

#Correlacion de Costald, volumen molar en cm3/mol
#Reglas de mezclado Costald

def V_a( Tc , Pc , w ): #entradas como np.array cm3/mol
    m = Tc.size #Contador para
    va = np.zeros( m )
    # va = RTc/Pc*( alfa + beta*w + gam*w**2 )
    #Constantes ajustadas para fracciones SARA
    alf = np.array([ 0.332482 , 1.45477 , -0.151344 , 1.04739 ]) 
    bet = np.array([ -0.163673 , -2.7688 , 0.660325 , -0.835364 ])
    gam = np.array([ 0.0494277 , 1.56843 , -0.294554 , 0.205254 ])
    #Saturado
    va[0] = ( 83.14472*Tc[0]/Pc[0] )*( alf[0]+bet[0]*w[0]+gam[0]*w[0]**2 )
    #Aromaticos
    va[1] = ( 83.14472*Tc[1]/Pc[1] )*( alf[1]+bet[1]*w[1]+gam[1]*w[1]**2 )
    #Resinas
    va[2] = ( 83.14472*Tc[2]/Pc[2] )*( alf[2]+bet[2]*w[2]+gam[2]*w[2]**2 )
    #Asfaltenos
    for i in range( 3 , m ):
        va[i] = ( 83.14472*Tc[i]/Pc[i] )*( alf[3]+bet[3]*w[i]+gam[3]*w[i]**2 )
    return va

def V_a_asf( Tc , Pc , w , n ):
    Va_asf = np.zeros( n )
    for i in range( 0 , n ):
        Va_asf[ i ] = ( 83.14472*Tc[i + 3]/Pc[i + 3] )*( 1.04739 - 0.835364*w[i + 3] + 0.205254*w[i + 3]**2 )
    return Va_asf

def V_a_solv( Tc , Pc , M , w ):
    ai=np.array([ 0.2905331 , -0.08057958 , 0.02276965] )
    return ( 83.14472*Tc/Pc )*( ai[0] + ai[1]*w + ai[2]*w**2 )

#Reglas de mezclado Costald, usar np.array en todas las reglas de mezclado

def Vcm_Costald(x,vc): #usar con np.array
    vi_a=np.power(vc,2/3)
    vi_b=np.power(vc,1/3)
    return 0.25*(np.vdot(x,vc)+3*np.dot(x,vi_a)*np.dot(x,vi_b))

def Tcm_Costald(x,vc,Tc,vm): #vm es escalar, x,vc, Tc son np.array
    z=np.power(vc*Tc,0.5)
    return(np.dot(x,z)**2)/vm

def w_m( x , w ): #Todas las entradas son np.array
    return np.dot(x,w)

"""
Las propiedades de mezcla se determinan con Vcm, Tcm y wm.
V_Costald cm3/mol, sirve para mezclas o componentes simples.
va es vaster o v* por como aparece en la correlación
"""
#volumen molar en punto de saturacion - Correlacion COSTALD para un componente
def V_Costald( T , va , w , Tc ): #uno o varios componentes, para mezcla usar las reglas 
    a = np.array([-1.52816 , 1.43907 , -0.81446 , 0.190454 , -0.296123 , 0.386914 , -0.0427258 , -0.0480645])
    Tr = T/Tc
    Tau = 1 - Tr
    vo = lambda Tau: 1 + a[ 0 ]*Tau**(1/3) + a[ 1 ]*Tau**(2/3) + a[ 2 ]*Tau + a[ 3 ]*Tau**(4/3)
    vd = lambda Tr: ( a[ 4 ] + a[ 5 ]*Tr + a[ 6 ]*Tr**2 + a[ 7 ]*Tr**3 )/( Tr - 1.00001 )
    # print("vo_sat", vo( Tau ) )
    # print("vd_sat",vd( Tr ) )
    return va*vo( Tau )*( 1 - w*vd( Tr ) )

#Lista de Vsat usando Costald 
def Vs_Costald( T , va , w , Tc ): #Solo T es escalar, los demas son np.array
    m = Tc.size
    vs = np.zeros( m )
    for i in range( m ):
        vs[ i ] = V_Costald( T , va[ i ] , w[ i ] , Tc[ i ] )
    return vs

#Volumen molar de una mezcla usando la correlacion COSTALD
def Vm_Costald( T, Va, Tc, w, x ): #Solo T es escalar, el resto np.array
    Vm = Vcm_Costald( x , Va )
    Tcm = Tcm_Costald( x , Va , Tc , Vm )
    wm = w_m( x , w )
    return V_Costald(T, Vm, wm, Tcm)

"""
Correlación de Li, es unicamente para mezclas - NO USADO
"""
#Reglas de Mezclado Li

#Phi tambien sirve para el modelo de solucion regular
#Phi para uno de los componentes, especificar con n
def Phi_i( x , v , n ): #x y v np.array, n es entero
    return x[ n ]*v[ n ]/np.dot( x , v )

def Phi( x , v ):
    return x*v/np.dot( x, v )

#Todas las entradas son np.array, el último es un array en forma de matriz
def Tcm_Li(x,Vc,Tc,kij):
    n=Tc.size
    Tcbin=np.zeros((n,n))
    phi=np.zeros(n)
    f=lambda i,j: (1-kij[i,j])*(Tc[i]*Tc[j])**0.5
    for i in range(n):
        Tcbin[i,i]=Tc[i]
        phi[i]=Phi_i(x,Vc,i)
        for j in range(i,n):
            Tcbin[i,j]=f(i,j)
            Tcbin[j,i]=Tcbin[i,j]
    phi=phi.T
    Tcm=np.dot(np.dot(Tcbin,phi).T,phi)
    return Tcm

#Correlacion de Li - NO USADO
def V_Li(T,x,Tc,Pc,w,Tcm): #Tc y Pc son np.array
    n=w.size
    zra=0.29056-0.08775*w
    zram=np.dot(x,zra)
    exponente=1+(1-T/Tcm)**0.2857
    f= lambda i: x[i]*Tc[i]/Pc[i]
    suma=0
    for i in range(n):
        suma=f(i)+suma
    return 83.14472*suma*zram**exponente

#Volumen subenfriado metodo de Chang y Zhao.
#V. Componentes puros.
def v_sub(P, T, Vs, Pc, Tc, w):
    Tr = T/Tc
    a = [-170.335, -28.578, 124.809, -55.5393, 130.01]
    b = [0.164813, -0.0914427]
    A = a[ 0 ] + a[ 1 ]*Tr + a[ 2 ]*Tr**3 + a[ 3 ]*Tr**6 + a[ 4 ]/Tr
    B = b[ 0 ] + w*b[ 1 ]
    C = np.exp( ( 1.00588 - Tr )**B )
    Pv = Pc*10**( ( 7/3 )*( 1 + w )*( 1 - 1/Tr ) )
    return Vs*( A*Pc + C*( P - Pv ) )/( A*Pc + np.exp( 1 )*( P - Pv ) )

def V_Comp( P, T, Vs, Pc, Tc, w):
    m = Vs.size
    V = np.array([])
    for i in range( 0 , m ):
        V = np.append( V , v_sub( P, T, Vs[i], Pc[i], Tc[i], w[i]) )
    return V


"""
modelos termodinamicos
Parametros de solubilidad para el modelo de solucion regular. Las correlaciones
son las utilizadas por Akbarzadeh. Dichas correlaciones se dan en MPa**0.5, por
lo que se convierte a Bar al multiplicar por *(10**.5)

Suposiciones: Solo encontramos resinas y asfaltenos en la fraccion pesada, por 
lo que el coeficiente de reparto "Ki" solo se programa para esos componentes.
"""
#Parametros de solubilidad SARA
def del_SARA(T,M): #M es Numpy array, parametro de sol. [=] MPa**0.5
    delt = np.zeros( M.size )
    #Saturados
    delt[0] = ( 23.01893 - 0.0222*T )
    #Aromaticos
    delt[1] = ( 26.333 - 0.0204*T )
    # delt[2] = ( 14821/150 - (4/15)*T ) #Opcion para resinas, se prefiere calcular de otra forma
    A = 0.579 - ( 7.5e-4 )*T
    v = lambda n: 1.493*M[ n ]**0.936
    #Resinas y asfaltenos
    for i in range( 2 , len( M ) ):
        delt[ i ] = ( 1e3*A*M[ i ]/v( i ) )**0.5
    return delt

#Parametros de solubilidad solo para resinas y asfaltenos
def del_asf( T , M , n ):
    delt_asf = np.zeros( n )
    A = (0.579 - (7.5e-4)*T)*1e3
    v = lambda n: 1.493*M[ n ]**0.936
    for i in range( 0 , n ):
        delt_asf[ i ] = ( A*M[i + 3]/v(i + 3) )**0.5
    return delt_asf

#ln(Gam_i) con modelo de solucion regular, para un componente en fase pesada
def ln_Gam_i( T ,vs , vm_h, x_h , d , n ): #np.array, menos T (float) y n (int)
    phi_h = Phi( x_h, vs[ 2 : x_h.size + 2 ] ) #El contador excluye al solvente, saturados y aromaticos
    dm_h = np.dot( d[ 2 : x_h.size + 2 ], phi_h )
    # print(vs[ n ])
    ln_G = 1 + np.log( vs[ n ]/vm_h ) - vs[ n ]/vm_h + ( vs[ n ]/( 8.314472*T ) )*( d[ n ] - dm_h)**2
    return ln_G

#ln(G_i) (NOTA: NO ES energia de Gibbs) con modelo de solucion regular, para fase ligera
def ln_G_i( T ,vs , vm_h, x , d , n ):
    phi = Phi( x , vs )
    dm_h = np.dot( d , phi )
    # print(vs[ n ])
    ln_g = 1 + np.log( vs[ n ]/vm_h ) - vs[ n ]/vm_h + ( vs[ n ]/( 8.314472*T ) )*( d[ n ] - dm_h)**2
    return ln_g

#ln(Gam_i) para todos los componentes
def ln_Gam( T , vs , vm_h , x_h , d ):
    n = x_h.size
    G = np.array([])
    for i in range( 2 , n + 2 ):
        G = np.append( G, ln_Gam_i(T, vs, vm_h, x_h, d, i) )
    return G

#ln(G_i) (NOTA: NO ES energia de Gibbs) para todos los componentes
def ln_G( T , vs , vm_h , x_h , d ):
    G = np.array([ ])
    for i in range(0 , x_h.size):
        G = np.append( G, ln_G_i( T ,vs , vm_h, x_h , d , i ) )
    return G

#Coeficiente de reparto individual
def K_i( T ,vs , vm_l, vm_h , x_l , x_h , d , n): #vs d y x son array, vm se calcula antes
    m = x_h.size
    phi_l = Phi( x_l , vs )
    phi_h = Phi( x_h, vs[ 2 : m + 2 ] )
    dm_l = np.dot( d , phi_l )
    dm_h = np.dot( d[ 2 : m + 2 ], phi_h )
    ln_Ki = lambda n: np.log( vm_h/vm_l )+vs[ n ]*( 1/vm_h - 1/vm_l )+( vs[ n ]/( 8.314472*T ) )*( (d[n] - dm_l)**2 - (d[n] - dm_h)**2 )
    return np.exp( ln_Ki( n ) )

#Coeficiente de reparto para todos los componentes
def K_ihl( T, vs, vm_l, vm_h, x_l, x_h, d):
    n = x_h.size
    Ki = np.array([ ])
    for i in range( 2 , n + 2 ):
        Ki = np.append( Ki, K_i( T ,vs , vm_l, vm_h , x_l , x_h , d , i) )
    return Ki

#Yi Michelsen componente individual
def Yi_Mich_i( T ,vs , vm_h, x_h , d , z , n ):
    ln_Gi =  ln_Gam_i( T ,vs , vm_h, x_h , d , n )
    ln_G_z =  ln_Gam_i( T ,vs , vm_h, z[ 2 : x_h.size + 2 ] , d , n )
    # print(ln_G_z)
    return np.exp( np.log( z[ n ] ) + ln_G_z - ln_Gi )

#Yi Michelsen para todos los componentes
def Yi_Mich( T ,vs , vm_h, x_h , d , z ):
    n = x_h.size
    Y = np.array([ ])
    for i in range( 2 , n + 2 ):
        Y = np.append( Y , Yi_Mich_i( T ,vs , vm_h, x_h , d , z , i ) )
    return Y

#Pendiente: de K_i, programar para que calcule todos los coeficientes de reparto y
#los guarde en un array.
    
def ferror_floc(z,k): #z y k son np. array
    return np.dot(z,k-1)

def x_h(z,k,Teta):
    xi = z*k/(1+Teta*(k-1))
    return xi/sum(xi)

def x_l(z,k,Teta):
    xi = z/(1+Teta*(k-1))
    return xi/sum(xi)

#Pruebas
if __name__=="__main__":
    #Prueba de distribucion Gamma
    A=PseudoComp(0.60447455952728, 412.3, 10,10) #z_pseudo, Mprom, Mmin, n_puntos, alfa estandar es 3.5
    print("z fracciones",A[0])
    # x=np.array([0.3,0.7])
    # Tc=np.array([450.0,380.7])
    # Pc=np.array([33.7,44.1])
    # Vc=np.array([236.8,321.0])
    # k=np.array([0.001,0.2]) #Coeficiente de reparto
    # w=np.array([0.023,0.014])
    # Vm=Vcm_Costald(x, Vc)
    # Tc_Costald=Tcm_Costald(x, Vc, Tc, Vm)
    # Tc_Li=Tcm_Li(x, Vc, Tc, Kij(Vc))
    # # vsat_Costald = V_Costald(350, 269.77, w_m, Tc_Costald)
    # vsat_Li = V_Li(350,x,Tc,Pc,w,Tc_Li)
    # print("Vm_Costald:",Vc)
    # print("Tcm_Costald:",Tc_Costald,"\t","Tcm_Li:",Tc_Li)
    # # print("v Costald:",vsat_Costald,"\t",vsat_Li)
    # # print("Kij:",Kij(Vc))
    # print("xh:",x_h(x,k,0.4),"\n","xl:",x_l(x,k,0.4))
    # print("suma xh:",sum(x_h(x,k,0.4)),"\t","suma xl:",sum(x_l(x,k,0.4)))
    

"""
Resultados esperados:
[1.00759128e-05 1.24747424e-03 1.33545576e-02 5.26786256e-02
 1.13165842e-01 1.54018402e-01 1.40916773e-01 8.71108667e-02
 3.45529782e-02 7.41896455e-03]

"""