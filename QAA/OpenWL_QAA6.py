# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 22:26:49 2020
Implemenation of Quasi-Analytical Algorithm (QAA)
http://www.ioccg.org/groups/Software_OCA/QAA_v6_2014209.pdf
@author: ZengC
"""
import numpy as np

def rrs_from_Rrs(Rrs,a=0.52,b=1.7):
    """
    estimate rrs from Rrs,
    Input:
        Rrs as N*1 array
        a,b are optional parameters
    Output:
        rrs as N*1 array
    """
    rrs=Rrs/(a+ np.multiply(Rrs,b))
    return rrs

def Rrs_from_rrs(rrs,a=0.52,b=1.7):
    """
    estimate Rrs from rrs,
    Input:
        rrs as N*1 array
        a,b are optional parameters
    Output:
        Rrs as N*1 array
    """
    #rrs=Rrs/(a+ np.multiply(Rrs,b))
    Rrs= a*rrs/(1-np.multiply(rrs,b))
    return Rrs

def u_from_rrs(rrs,g0=0.08945,g1=0.1247):
    """
    estimate u from rrs,
    Input:
        rrs as N*1 array
        g0,g1 are optional parameters
    Output:
        u as N*1 array
    """    
    u= ( -g0 + np.sqrt(g0**2+np.multiply(rrs,4*g1)) ) /(2*g1)
    return u

def a_lambd0_from_rrs(wavelength,rrs,lambd=555,a_w_lambd0=None,h55x=[-1.1459,-1.36583,-0.469266],A=5):  #QQA_v5
    """
    calculate a_lambda at lambda 0 using rrs approximation, for Rrs(670) < 0.0015 sr-1 case (QAA5)
    Input:
        wavelength: N*1  of wavelengths
        rrs: N*1 of bands
        h: 3*1 parameters
        lambd: the lambda 0 wavelength, scalar value in nm
        a_w_lambd0: scalar value for a_w at lambda0
    Output:
        a_lambda at lambda0,  scalar value
    """
    if a_w_lambd0 is None:
        a_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        a_w=a_w[:,[0,1]]
        a_w_lambd0=np.interp(lambd,a_w[:,0], a_w[:,1])
    key_rrs=np.interp([443,490,lambd,670], wavelength, rrs)
    chi=np.log10( (key_rrs[0]+key_rrs[1])/(key_rrs[2]+A*key_rrs[3]/key_rrs[1]*key_rrs[3]) )
    a_lambd0=a_w_lambd0 + np.power(10,h55x[0]+h55x[1]*chi+h55x[2]*chi**2)  #absoprtoin at 55x nm.
    #print("a_w_lambd0: {}, key_rrs:{}, chi:{} chi_power: {}".format(a_w_lambd0,key_rrs,chi,h[0]+h[1]*chi+h[2]*chi**2))
    #print("a_lambd0: ",a_lambd0)
    return a_lambd0

def a_lambd0_from_Rrs(wavelength,Rrs,lambd=670,a_w_lambd0=None,h66x=[0.39,1.14]):  #QQA_v6 others
    """
    calculate a_lambda at lambda 0 using rrs approximation, for Rrs(670) >= 0.0015 sr-1 case (QAA6 other)
    Input:
        wavelength: N*1  of wavelengths
        Rrs: N*1 of bands
        lambd: the lambda 0 wavelength, scalar value in nm
        a_w_lambd0: scalar value for a_w at lambda0
    Output:
        a_lambda at lambda0,  scalar value
    """
    if a_w_lambd0 is None:
        a_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        a_w=a_w[:,[0,1]]
        a_w_lambd0=np.interp(lambd,a_w[:,0], a_w[:,1])
    key_Rrs=np.interp([443,490,lambd,670], wavelength, Rrs)
    a_lambd0=a_w_lambd0 + h66x[0]* np.power(key_Rrs[2]/(key_Rrs[0]+key_Rrs[1]) ,h66x[1])  #absoprtoin at 55x nm.
    #print("key_Rrs:{},a_lambd0:{}".format(key_Rrs,a_lambd0))
    return a_lambd0

def bb_p_lambd0_from_u_and_a(wavelength,u,a_lambd0,lambd,bb_w_lambd0=None):
    if bb_w_lambd0 is None:
        bb_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        bb_w=bb_w[:,[0,2]]/2   # convert b to bb, water as half for backscatter
        bb_w_lambd0=np.interp(lambd, bb_w[:,0], bb_w[:,1])
    u0=np.interp(lambd, wavelength, u)
    bb_p_lambd0=u0*a_lambd0/(1-u0) - bb_w_lambd0
    #print('bb_p_lambd0:',bb_p_lambd0)
    return bb_p_lambd0

def bbp_labmda0_from_rrs(wavelength,rrs,Rrs=None,a_w=None,bb_w=None,\
                         params={'h55x':[-1.1459,-1.36583,-0.46927],'A':5,'h66x':[0.39,1.14]}): #set the initial values of abs and bbp
    """
    calculate a_lambda at lambda 0 using rrs approximation, for Rrs(670) >= 0.0015 sr-1 case (QAA6 other)
    Input:
        wavelength: N*1  of wavelengths
        rrs , Rrs: N*1 of bands
        a_w, bb_w: water abs and bb,in format N*2
        params: the parameters that can be adjusted; params as a dict of h, A, 
    """
    if (Rrs is None) and (rrs is None): 
        print('----------no valid input of Rrs or rrs, will exit-----')
        return None
    if Rrs is None: Rrs=Rrs_from_rrs(rrs)
    if rrs is None: rrs=rrs_from_Rrs(Rrs)
    u=u_from_rrs(rrs)
    
    if np.interp(670, wavelength, Rrs)<0.0015: 
#    if False:
        lambd=555
    else: 
        lambd=670
    
    if a_w is None:         a_w_lambd0=None
    else:         a_w_lambd0=np.interp(lambd,a_w[:,0], a_w[:,1])
    
    if bb_w is None:        bb_w_lambd0=None
    else:        bb_w_lambd0=np.interp(lambd,bb_w[:,0], bb_w[:,1])
    
    #print('a_w_lambd0:{},bb_w_lambd0:{}'.format(a_w_lambd0,bb_w_lambd0))    
    
    if np.interp(670, wavelength, Rrs)<0.0015 : # QAA v5 case
#    if False:
        lambd=555  #use 555nm as the 55x band
        
        #h=[-1.1459,-1.36583,-0.46927]  #according to QAA v6 spreedsheet
        #A=5
#        a_lambd0=a_lambd0_from_rrs(wavelength,Rrs,h,lambd,a_w_lambd0=a_w_lambd0)
        a_lambd0=a_lambd0_from_rrs(wavelength,rrs,lambd,a_w_lambd0=a_w_lambd0,h55x=params['h55x'],A=params['A'])
        
        bb_p_lambd0=bb_p_lambd0_from_u_and_a(wavelength,u,a_lambd0,lambd,bb_w_lambd0=bb_w_lambd0)
    else:   # QAA v6 other cases
        lambd=670
        a_lambd0=a_lambd0_from_Rrs(wavelength,Rrs,lambd,a_w_lambd0=a_w_lambd0,h66x=params['h66x'])
        bb_p_lambd0=bb_p_lambd0_from_u_and_a(wavelength,u,a_lambd0,lambd,bb_w_lambd0=bb_w_lambd0)
    #print('bb_p_lambd0:{},a_lambd0:{}'.format(bb_p_lambd0,a_lambd0))    
    
    return bb_p_lambd0,lambd 
    
def bb_p_lambda(wavelength,bb_p_lambd0,rrs,lambd0=555,eta_coeff=[2.0,1.2,-0.9]):  #tune params with data if possible
    """
    calculate bb_p_lambda over a given wavelength list
    Input:
        wavelength: N*1  of wavelengths
        bb_p_lambd0: a scalar value to anchor the magnitude of the curve
        rrs: N*1 of bands
        eta_coeff: the 3 coefficients for eta
    """
    bands=[443,555]
#    bands=[440,555]   
    
    rrs_key=np.interp(bands, wavelength, rrs)
    eta=eta_coeff[0] * (1 - eta_coeff[1] *np.exp(eta_coeff[2]*np.divide(rrs_key[0],rrs_key[1])) )
    bb_p_lambda=list()
    for x in wavelength:
        bb_p=bb_p_lambd0*np.power(lambd0/x,eta)
        bb_p_lambda.append(bb_p)
    return bb_p_lambda


def a_from_u_and_b(wavelength,u,bb_p,bb_w=None): #u=bb/(a+bb)
    """
    calculate a_lambda over a given wavelength list
    Input:
        wavelength: N*1  of wavelengths
        u,bb_p,bb_w: N*1 list of u and bb_p over the wavelength list 
    """
    if bb_w is None:
        bb_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        bb_w=bb_w[:,[0,2]]
    bb_w=np.interp(wavelength, bb_w[:,0], bb_w[:,1])
    a_lambd=np.multiply(1-u, np.divide((bb_w+bb_p), u))
    return a_lambd

def a_coefs_from_rrs(wavelength,rrs,A=[0.74,0.2,0.8],B=[0.015,0.002,0.6],xi_bands=[442.5, 415.5]):
    """
    calculate coefficients for a components 
    Input:
        wavelength: N*1  of wavelengths
        rrs: N*1 of bands
        A: coefficients for zeta
        B: cofficients for xi
        xi_bands: the 2 bands used for xi 
    """
    a_coeff_bands=[443,555]
    rrs_key=np.interp(a_coeff_bands, wavelength, rrs)
    zeta=A[0]+ A[1]/ (A[2] + np.divide(rrs_key[0], rrs_key[1]))
    S=   B[0]+ B[1]/ (B[2] + np.divide(rrs_key[0], rrs_key[1]))
    xi=np.exp((xi_bands[0]- xi_bands[1])*S)
    
    return zeta, xi,S

def a_dg_lambd(wavelength,a_dg_lambd_0,S,lambd0=443):
    a_dg_lambd=list()
    for x in wavelength:
        temp=a_dg_lambd_0 * np.exp(-1*S*(x -lambd0))
        a_dg_lambd.append(temp)
    return a_dg_lambd

def a_phyto_from_a(wavelength,a,zeta, xi,S,a_w=None): #seperate a_ph from a and aw
    if a_w is None:
        a_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        a_w=a_w[:,[0,1]]   
    
    bands=[412,443]
    a_w_key=np.interp(bands, a_w[:,0], a_w[:,1])
    a_w_lambd=np.interp(wavelength, a_w[:,0], a_w[:,1])
    a_key=np.interp(bands, wavelength, a)
    
    a_g_443=np.divide(a_key[0]-np.multiply(a_key[1],zeta) , xi-zeta ) - \
            np.divide(a_w_key[0]-np.multiply(a_w_key[1],zeta) , xi-zeta ) 
    
    a_dg= a_dg_lambd(wavelength,a_g_443,S)
    #print('a_g_443:{},a_dg:{}'.format(a_g_443,a_dg))
    
    a_ph=np.subtract(a-a_dg,a_w_lambd)
    a_ph[a_ph<0] = 0  #set negative aph as 0
    
    return a_ph

def QAA6(params,rrs,wavelength,Rrs=None, a_w=None,bb_w=None):
    """
    the main function for QQA version 6:
    INput:
        wavelength: the wavelength of the input Rrs or rrs bands
        rrs,Rrs: the water reflectance, in shape of N*1; if rrs is given, then Rrs is ignored
        a_w, bb_w: water abs and backscter, in shapef of N*2; first column as wavelength, 2nd column value (1/m)
    output:
        a_ph:  the absorption of phytoplankton at input bands, shape of N*1;  set 0 for failed values; unit: (1/m)
    """
    if a_w is None:  #load the default water abs and bb
        ab_w=np.loadtxt('../IOPfiles/waterIOP_sea_water.txt',skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        bb_w=ab_w[:,[0,2]] /2  #convert b to bb for water, 
        a_w=ab_w[:,[0,1]] 
    #print('aw:',np.interp(wavelength, a_w[:,0],a_w[:,1]))
    
    if rrs is None and Rrs is None: 
        print('------------please provide valid data to run , exit--------------')
        return None
    
    if rrs is None :
        rrs=rrs_from_Rrs(Rrs)
    u=u_from_rrs(rrs)
    
    ##seperate the parameter (input is 1d) to meaningful sub parameters used in each step
#    params={'h55x':[-1.1459,-1.36583,-0.46927],'A':5,'h66x':[0.39,1.14],\
#            'eta_coeff':[2.0,1.2,-0.9],\
#            'zeta_coeff':[0.74,0.2,0.8],'xi_coeff':[0.015,0.002,0.6]}
#
#    params=[-1.1459,-1.36583,-0.46927,5,0.39,1.14,2.0,1.2,-0.9,0.74,0.2,0.8,0.015,0.002,0.6]
    param1={'h55x':params[0:3],'A':params[3],'h66x':params[4:6]}
    eta_coeff=params[6:9]
    zeta_coeff=params[9:12]
    xi_coeff=params[12:15]
    
    bb_p_lambd0,lambd0_4_bbp=bbp_labmda0_from_rrs(wavelength,rrs,a_w=a_w,bb_w=bb_w,params=param1)  #get bb0
    
    bb_p=bb_p_lambda(wavelength,bb_p_lambd0,rrs,lambd0_4_bbp,eta_coeff=eta_coeff)
    
    a=a_from_u_and_b(wavelength,u,bb_p,bb_w=bb_w) 
    
    zeta, xi,S=a_coefs_from_rrs(wavelength,rrs,A=zeta_coeff,B=xi_coeff)  #A=[0.74,0.2,0.8],B=[0.015,0.002,0.6],xi_bands=[442.5, 415.5]
    
    a_ph=a_phyto_from_a(wavelength,a,zeta, xi,S, a_w=a_w)
    
    return a_ph

#    params={'h55x':[-1.1459,-1.36583,-0.46927],'A':5,'h66x':[0.39,1.14],\
#            'eta_coeff':[2.0,1.2,-0.9],\
#            'zeta_coeff':[0.74,0.2,0.8],'xi_coeff':[0.015,0.002,0.6]}

##if run as a main program,[rather than import as a model]
if __name__ == "__main__":
    ##parse arguments
    import os,sys,argparse
    parser = argparse.ArgumentParser(description='Implementation of QAA version 6')
    parser.add_argument('-N','--wavelength', metavar='', type=str, 
                        default='412,443,489,510,555,670',
                        help='the wavelength of the input Rrs or rrs bands. default:[412,443,489,510,555,670]')
    parser.add_argument('-U', '--rrs', metavar='', type=str,
                        default=None,
                        help='an sample of under water remote sensing reflectance, dimensions:N*1')
    parser.add_argument('-R', '--Rrs', metavar='', type=str,
                        default=None,
                        help='an sample of above water remote sensing reflectance, dimensions:N*1; \
                        if rrs is given, then Rrs is ignored')
    parser.add_argument('-F', '--file', metavar='', type=str,
                        default=None,
                        help='a csv file of *above* water remote sensing reflectance, if --file given, --rrs/--Rrs are ignored ; \
                        dimenion:N*M; where there is N samples, with M bands, each row is a sample;\
                        example data: Rrs_QAA_test.csv')                        
    parser.add_argument('-W', '--waterIOP', metavar='', 
                        default=None,
                        help='water absorption curve, dimensions: N*3: 1st column as wavelength, 2nd column abs value (1/m), 3rd for scatter, consistent with Hydrolight;\
                        if not provided, use water abs and bb from ../IOPfiles/waterIOP_sea_water.txt')
    parser.add_argument('--h55x', metavar='', default=[-1.1459,-1.36583,-0.46927], type=float,nargs=3, help='parameters in QAA; step 2 left IF, h0,h1,h2')       
    parser.add_argument('--A', metavar='', default=5, type=float, help='parameters in QAA; step 2 left IF, coefficient in the log(...)') 
    parser.add_argument('--h66x', metavar='', default=[0.39,1.14], type=float,nargs=2, help='parameters in QAA; step 2 Right ELSE') 
    parser.add_argument('--eta', metavar='', default=[2.0,1.2,-0.9], type=float,nargs=3, help='parameters in QAA; step 4') 
    parser.add_argument('--zeta', metavar='', default=[0.74,0.2,0.8], type=float,nargs=3, help='parameters in QAA; step 7A') 
    parser.add_argument('--xi', metavar='', default=[0.015,0.002,0.6], type=float,nargs=3, help='parameters in QAA; step 7B') 
    parser.add_argument('-S', '--save', action='store_true', help='flag to save the a_ph as .csv file, default is False')
    
    #args = parser.parse_args(['-F','Rrs_QAA_test.csv'])
    args = parser.parse_args()
    
    ##determine the input Rrs/rrs
    wavelength=[float(x) for x in args.wavelength.split(',')]
    if not args.rrs is None: 
        rrs=[float(x) for x in args.rrs.split(',')]
    else: 
        rrs=None
    if not args.Rrs is None: 
        Rrs=[float(x) for x in args.Rrs.split(',')]
    else:
        Rrs=None
    if (not (args.file is None)) and os.path.exists(args.file):
        Rrs=np.loadtxt(args.file)  #N*M where there is N samples, with M bands, each row is a sample
        rrs=None #overwrite the rrs input
    if (not (Rrs is None) and len(wavelength) != np.transpose(np.array(Rrs)).shape[0] ) or (not (rrs is None) and len(wavelength) != len(rrs)):
        print("=====incosistent dimensions of the input: wavelength and rrs/Rrs, program will exit =====")
        sys.exit()
    
    ##load water IOP if needed
    if (not (args.waterIOP is None)) and os.path.exists(args.waterIOP):
        ab_w=np.loadtxt(args.waterIOP,skiprows=10,comments='-1')  #3 columns: [wavelen(nm), a(1/m), b(1/m)]
        bb_w=ab_w[:,[0,2]] *0.5  #convert b to bb for water, 
        a_w=ab_w[:,[0,1]]
    else:
        a_w=None
        bb_w=None

    #default QAA params, for their meaning, pease check [QAA document](https://www.ioccg.org/groups/Software_OCA/QAA_v6_2014209.pdf)
    #params=[-1.1459,-1.36583,-0.46927,5,0.39,1.14,2.0,1.2,-0.9,0.74,0.2,0.8,0.015,0.002,0.6]
    ##combine the parameters (input is 1d) to pass into QAA
    params=[*args.h55x, args.A, *args.h66x, *args.eta, *args.zeta, *args.xi] 
    
    ##debug
    h55x,A,h66x=params[0:3],params[3],params[4:6]
    eta_coeff,zeta_coeff,xi_coeff=params[6:9],params[9:12],params[12:15]
    #print("===input params:h55x :{},A:{},h66x:{},\n eta:{},zeta:{},xi :{}".format(h55x,A,h66x,eta_coeff,zeta_coeff,xi_coeff))
    
    ##currently not support matrix operation, need to pass as vectors per sample
    if np.array(Rrs).ndim>1:  #multiple samples
        nsample=Rrs.shape[0]
        a_ph=np.empty([nsample, len(wavelength)])
        for idx in range(nsample):
            a_ph[idx,:]=QAA6(params,rrs,wavelength,Rrs=Rrs[idx,:], a_w=a_w,bb_w=bb_w)
    else: #single sample
        a_ph=QAA6(params,rrs,wavelength,Rrs=Rrs, a_w=a_w,bb_w=bb_w)
    
    if args.save:
        f_aph=os.path.join(os.path.dirname(args.file),'aph_'+os.path.basename(args.file))
        np.savetxt(f_aph, a_ph, delimiter=",")
        print("===QAA completed,the absorption of phytoplankton is written to file: {} ===".format(f_aph))
    else:
        print("===QAA completed,the absorption of phytoplankton at input bands [N*1,set 0 for failed values; unit: (1/m)]: \n {} ===".format(a_ph))
