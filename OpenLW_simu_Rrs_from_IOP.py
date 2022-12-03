
import os,sys, argparse
import numpy as np
import pandas as pd
from scipy import interpolate   

def load_iop_file(str_path,header_lines=10,skip_footer=1,wavelength=range(400,851,5)):
    ##--------------LOAD water quality parameter and resample if provided. -------------
    if not os.path.exists(str_path):
       print(f'input IOP file {str_path} does not exist, will exit')
       return None
    IOP=np.genfromtxt(str_path,skip_header=header_lines,skip_footer=skip_footer)
    
    if isinstance(wavelength, type(None)):  #not resampled
        return IOP[:,1:]
    else:    
        f = interpolate.interp1d(IOP[:,0], np.transpose(IOP[:,1:]), fill_value='extrapolate')
        return np.squeeze(np.transpose(f(wavelength)))  #need to convert the shape
 
def OpenLW_simu_Rrs_from_IOP(v_chl,v_mspm ,v_cdom=[0.994],sDir_IOP='./IOPfiles',IOPname='CCIW',
                             wavelength=range(400,851,5),bAB=False,f=0.319,bLee=False,str_path_waterIOP='waterIOP_sea_water.txt',
                             header_lines=10,water_bb2b_ratio=0.5,chl_bb2b_ratio=0.01,mspm_bb2b_ratio=0.02,
                             theta_sun=0,theta_view=0,wind_speed=5,eta=0.52,gamma=1.6):
    """
    ## this is the implementation of water-leaving raidance (Rrs) simulation as used in: 
    ## 		Zeng C, Binding C. The effect of mineral sediments on satellite chlorophyll-a retrievals from line-height algorithms using red and near-infrared bands. Remote Sensing. 2019;11(19). doi:10.3390/rs11192306
    ## 		Author: Chui Zeng;  chqzeng@gmail.com
    ## 
    ## It uses a very simplified light and water interactcion empricical model:  Rrs = f/PI *  bb / (a + bb) according to:
    ## 	Albert, A. and C. D. Mobley (2003). "An analytical model for subsurface irradiance and remote sensing reflectance in deep and shallow case-2 waters." Optics Express 11(22): 2873-2890;
    ## for a more strict simulation, please use HydroLight/Ecolight [Mobley, et al,https://www.sequoiasci.com/product/hydrolight/]
    ## 
    ## Input arguments:
    ## 		v_chl,v_cdom,v_mspm: the value list of chlorophyll-a[mg/m^3], cdom[1/m], and mineral content/mspm[g/m^3]; default cdom is: 0.994 (1/m)
    ## 		sDir_IOP:  the directory where the IOP files are stored; default is the "./IOPfiles"
    ## 		IOPname  : the name ID key for a provided IOP set. corresponding to the files in: `./IOPfiles/` folder; default IOP model:" "
    ## 		wavelength: the list of wavelength to simulate the output Rrs; default: [400:5:850] nm
    ## 		bAB: flag to use Chl absorption with model: abs(chl) = A x chl^(1-B) model, otherwise use simpler: abs(chl)= a* x chl   model
    ##        bSave: save the simulation result in current directory, in format: 'OpenLW_simu_Rrs_from_IOP_[%s]_chl[%d]_mspm[%d]_cdom[%f].csv'
    ## 		f:  the emperical value f in formula: Rrs = f *  bb / (a + bb);  default use 0.319 for inland water after Jerome et al. (1988). 
    ## 		waterIOP: the file name of water IOP, in the './IOPfiles/'	folder, defualt is: 'waterIOP_sea_water.txt'  for clear fresh water
    ## 
    ## Output arguments: [Rrs_table,Rrs]
    ## 		Rrs: the water leaving radiance (sr-1); dimensions: N x M; where
    ##             N=length(v_chl)*length(v_mspm)*length(v_cdom),  M=length(wavelength)
    ##        Rrs_table: N x (M+3); the table format of Rrs, while add 3 columns at
    ##             beginning for: chl, mspm, and cdom
    """
    ##load water IOP
    #IOP_water=np.loadtxt(os.path.join(sDir_IOP,str_path_waterIOP),skiprows=header_lines,comments='-1')  #risky using that comment=??
    if isinstance(str_path_waterIOP, type(None)) :
        str_path_waterIOP='waterIOP_sea_water.txt'
    if not '/' in str_path_waterIOP: #only file name given, try to predict the path
        if os.path.exists(os.path.join('./',str_path_waterIOP)): 
            pass
        elif os.path.exists(os.path.join(sDir_IOP,str_path_waterIOP)): 
            str_path_waterIOP=os.path.join(sDir_IOP,str_path_waterIOP)
        else:
            print(f'failed to intepret water iop path: {str_path_waterIOP}, will exit')
            return None,None
    IOP_W=load_iop_file(str_path_waterIOP,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_W is None: 
        print(f"=====failed to load IOP file : {str_path_waterIOP}, will exit")
        return None,None        
    IOP_W_a,IOP_W_b=IOP_W[:,0],IOP_W[:,1]
    IOP_W_bb=IOP_W_b * water_bb2b_ratio;
    
    ##load mass-specific absoprtion of chl, cdom, and mspm
    str_path_abs_chl =os.path.join(sDir_IOP,f"{IOPname}_abs_coeff_chl{'_AB'if bAB else ''}.txt")
    str_path_abs_cdom=os.path.join(sDir_IOP,f'{IOPname}_abs_coeff_CDOM.txt')
    str_path_abs_mspm=os.path.join(sDir_IOP,f'{IOPname}_abs_coeff_mineral.txt')
    if not  os.path.exists(str_path_abs_chl): 
        print(f"=====IOP file not exist: {str_path_abs_chl}, will exit")
        return None,None
    if not  os.path.exists(str_path_abs_cdom): 
        print(f"=====IOP file not exist: {str_path_abs_cdom}, will exit")
        return None,None
    if not  os.path.exists(str_path_abs_mspm): 
        print(f"=====IOP file not exist: {str_path_abs_mspm}, will exit")
        return None,None
    
    ##--------------1. chl -------------
    IOP_CHL_a=load_iop_file(str_path_abs_chl,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_CHL_a is None: 
        print(f"=====failed to load IOP file : {str_path_abs_chl}, will exit")
        return None,None
    #IOP_temp[IOP_temp[:,0]>720,2]=0;  # set the B value for 720nm -800nm to zero; (according to Mobley's book) 
    ##--------------2. cdom -------------
    IOP_CDOM_a=load_iop_file(str_path_abs_cdom,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_CDOM_a is None: 
        print(f"=====failed to load IOP file : {str_path_abs_cdom}, will exit")
        return None,None
    ##--------------3. mspm -------------
    IOP_MSPM_a=load_iop_file(str_path_abs_mspm,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_MSPM_a is None: 
        print(f"=====failed to load IOP file : {str_path_abs_mspm}, will exit")
        return None,None
    
    #load mass-specific scatter properties of chl and mspm
    str_path_sct_chl=os.path.join(sDir_IOP,f'{IOPname}_scatter_coeff_chl.txt')
    str_path_sct_mspm=os.path.join(sDir_IOP,f'{IOPname}_scatter_coeff_mineral.txt')  
    #read chl mass-specific scatter
    IOP_CHL_b=load_iop_file(str_path_sct_chl,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_CHL_b is None: 
        print(f"=====failed to load IOP file : {str_path_sct_chl}, will exit")
        return None,None
    IOP_MSPM_b=load_iop_file(str_path_sct_mspm,header_lines=header_lines,skip_footer=1,wavelength=wavelength)
    if IOP_MSPM_b is None: 
        print(f"=====failed to load IOP file : {str_path_sct_mspm}, will exit")
        return None,None
	#bb/b ratio : from scatter to back-scatter
    IOP_CHL_bb=IOP_CHL_b*chl_bb2b_ratio;   
    IOP_MSPM_bb=IOP_MSPM_b*mspm_bb2b_ratio;   

    #mix the chl, mspm,cdom to different combinations  
    XX, YY,ZZ = np.meshgrid(v_chl,v_mspm,v_cdom)
    IOP_list = np.stack([XX.flatten(), YY.flatten(), ZZ.flatten()]) #the list of IOP in shape: 3xN; where N is the sample number
    #N=IOP_list.shape[-1]
    
    ##the simple cacluation of the Rrs with Rrs <--: bb/(a+bb) according to: 
    ##Gege, P. (2017). Chapter 2 - Radiative Transfer Theory for Inland Waters. Bio-optical Modeling and Remote Sensing of Inland Waters. D. R. Mishra, I. Ogashawara and A. A. Gitelson, Elsevier: 25-67.    
    #a_total=[chl]*a_chl + [cdom] *a_cdom + [mspm] *a_mspm +aw
    #bb_total=[chl] *bb_chl + [mspm] * bb_mspm +bbw 
    #Rrs = f .*  bb_total ./ (a_total+bb_total) ./pi;
    if bAB and IOP_CHL_a.ndim==2:  
        #need to calcuate a_chl_total= A.* [chl].^(-1*B), where A,B=IOP_CHL_a[:,0],IOP_CHL_a[:,1]
        ##reference: Bricaud, A., Babin, M., Morel, A., Claustre, H., 1995. Variability in the chlorophyll specific absorption coefficients of natural phytoplankton: analysis and parameterization. J. Geophys. Res. 100, 13,321-13,332.
        a_chl_total=IOP_CHL_a[:,0]*np.power(IOP_list[0,:].reshape([-1,1]), -1*IOP_CHL_a[:,1].reshape([1,-1]))
    else:   #calculate a_chl_total= A * [chl]
        a_chl_total=np.matmul(IOP_list[0,:].reshape([-1,1]),IOP_CHL_a.reshape([1,-1]))   
    a_total=a_chl_total+np.matmul(IOP_list[1,:].reshape([-1,1]),IOP_MSPM_a.reshape([1,-1]))+\
    np.matmul(IOP_list[2,:].reshape([-1,1]),IOP_CDOM_a.reshape([1,-1]))+IOP_W_a
    
    bb_total=np.matmul(IOP_list[0,:].reshape([-1,1]),IOP_CHL_bb.reshape([1,-1]))+\
    np.matmul(IOP_list[1,:].reshape([-1,1]),IOP_MSPM_bb.reshape([1,-1]))+IOP_W_bb    
    
    u=bb_total/(a_total+bb_total)
    
    
    ##for optical deep water: r_rs(λ)= f_rs(λ) * u(λ)
    if not isinstance(f, type(None)):  #user assigned f value
        f_rs= f
    elif bLee:
        ##method 1: Lee, Z.-P., Carder, K.L., Mobley, C.D., Steward, R.G., Patch, J.S., 1999. Hyperspectral remote sensing for shallow waters: 2. Deriving bottom depths and water properties by optimization. Appl. Opt. 38, 3831-3843.
        f_rs= 0.084 +0.017*u
    else :
        #method 2: Albert, A., 2004. Inversion technique for optical remote sensing in shallow water. Ph.D. Dissertation. Universita¨t Hamburg, Hamburg, Germany, p. 188.
        #wind_speed: is the wind speed in units of (m/s)
        #theta_sun: the sun zenith angle in water
        #theta_view: the viewing zenith angle in water
        f_rs=0.0512*(1+ 4.6659*u + 7.8387*u**2 + 5.4571*u**3)*(1+ 0.1098/np.cos(theta_sun/180*np.pi))*(1+ 0.4021/np.cos(theta_view/180*np.pi))*(1-0.0044*wind_speed)
    rrs=f_rs * u
    
    ##from r_rs(λ) (under water ) to Rrs(λ) (above water)
    #Lee, Z.-P., Carder, K.L., Mobley, C.D., Steward, R.G., Patch, J.S., 1998. Hyperspectral remote sensing for shallow waters: 1. A semi-analytical model. Appl. Opt. 37, 6329-6338.
    ##eta: the water-to-air radiance divergence factor
    ##gamma: accounts for the effects of internal reflection from water to air
    ##Rrs_surf: the reflections at the water surface, ignore here temporarily
    #eta,gamma=0.52, 1.6
    Rrs_surf=0
    Rrs=eta*rrs / (1-gamma*rrs) +Rrs_surf
    
    ##convert numpy array to pandas table as output
    cols=['CHL','MSPM','CDOM']+[f'{w}nm' for w in wavelength]
    rrs=pd.DataFrame(np.concatenate((np.transpose(IOP_list),rrs), axis=1),columns=cols)
    Rrs=pd.DataFrame(np.concatenate((np.transpose(IOP_list),Rrs), axis=1),columns=cols)
    #Rrs = (f *  bb_total / (a_total+bb_total) ) /np.pi
    return Rrs,rrs

#auxilary functions for the input argument intepretation
def str_to_value_list(v):
    if isinstance(v, type(None)):
        return None
    elif v.strip().startswith('[') and v.strip().endswith(']'):
        val=(v.strip())[1:-1]  #get the list of [val1,val2,....]
        if ',' in val:   #list split by ','
            return [float(x) for x in val.split(',')]
        elif ' ' in val:  #list split by whitespace
            return [float(x) for x in val.split(' ')]
        else:  #only one value
            return [float(val)]
    elif v.strip().count(':')==2 :  #interpret as a min_bound:interval:max_bound
        try: val=[float(x) for x in v.strip().split(':')]
        except: return None
        if len(val)==3:  return np.arange(val[0],val[2],val[1])  
        else: return None
    else:
        return None

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
        
##if run as a main program,[rather than import as a model]
if __name__ == "__main__":
    ##parse arguments
    
    parser = argparse.ArgumentParser(description='Implementat an approximate simulation of Hydrolight for water leaveing remote sensing reflectance using water IOP')
    parser.add_argument('--chl', metavar='', type=str_to_value_list, default=[1,2,3,4,5,10,20,50,100], help='the CHL [mg/m3] value range to be simulated in format:  min:interval:max , or [chl1,chl2, chl3,...]')
    parser.add_argument('--mspm', metavar='', type=str_to_value_list, default=[1,2,3,4,5,10,20], help='the MSPM [g/m3] value range to be simulated in format:  min:interval:max , or [mspm1,mspm2, mspm3,...]')
    parser.add_argument('--cdom', metavar='', type=str_to_value_list, default=[0.994], help='the CDOM [1/m] value range to be simulated in format:  min:interval:max , or [cdom1,cdom2, cdom3,...]')
    
    ##IOP files
    parser.add_argument('--sIOP', metavar='', type=str, default='./IOPfiles',help='the path to the direotory of IOP files. default: ./IOPfiles')  #sDir_IOP
    parser.add_argument('--IOPname', metavar='', type=str, default='LELW',help='the name of defualt IOP profile to be used. need to have IOP abs and scatter files stored in folder: sDir_IOP/[IOPname]_[abs|scatter]_coeff_[CDOM|chl|mineral].txt. default: LELW, for Lake Erie and Lake Winnipege averaged') 
    parser.add_argument('-N','--wavelength', metavar='', type=str_to_value_list, default=range(400,851,5),
                        help='the wavelength of the output Rrs or rrs bands. in format:  min:interval:max , or [wav1,wav2, wav3,...].  default:[400,405,410,...,850]')
    parser.add_argument('--bAB', metavar='',type=str2bool, nargs='?',const=False, default=False, help="whether use the abs_chl=A*[chl]^(-B) model? if yes, will need a ###_chl_AB.txt for absorption profile with two abs columns. default False, \n for True use:'yes', 'true', 't', 'y', '1' \n ; for False, use: 'no', 'false', 'f', 'n', '0', same for the other bool parameters")
    
    ##parameters
    parser.add_argument('-f', metavar='', default=None, type=float, help='parameters f for: f/Q * bb/(a+bb); e.g.,0.319. if set, then will ignore further Lee1998 model or Albert2004 models') 
    parser.add_argument('--bLee', metavar='',type=str2bool, nargs='?',const=False, default=False, help="whether use the Lee et al 1999 model for the rrs=f(u). default False, to use the Ablert2004 model; \n for True use:'yes', 'true', 't', 'y', '1' \n ; for False, use: 'no', 'false', 'f', 'n', '0'")

    #str_path_waterIOP
    parser.add_argument('-W', '--waterIOP', metavar='', default=None,
                        help='water absorption curve, dimensions: N*3: 1st column as wavelength, 2nd column abs value (1/m), 3rd for scatter, consistent with Hydrolight;\
                        if not provided, use water abs and bb from ../IOPfiles/waterIOP_sea_water.txt')
    parser.add_argument('--header_lines', metavar='', default=10, type=float, help='number of lines as header for the IOP profile files, default is 10') 
    parser.add_argument('--water_bb2b_ratio', metavar='', default=0.5, type=float, help='the water_bb2b_ratio: the ratio of backscater ratio to scatter for water, default is 0.5') 
    parser.add_argument('--chl_bb2b_ratio', metavar='', default=0.01, type=float, help='the chl_bb2b_ratio: the ratio of backscater ratio to scatter for chl, default is 0.01') 
    parser.add_argument('--mspm_bb2b_ratio', metavar='', default=0.02, type=float, help='the mspm_bb2b_ratio: the ratio of backscater ratio to scatter for mspm, default is 0.02') 
    
    
    ##observe geometry: theta_sun=0,theta_view=0,wind_speed=5,eta=0.52,gamma=1.6
    #wind_speed: is the wind speed in units of (m/s)
    #theta_sun: the sun zenith angle in water
    #theta_view: the viewing zenith angle in water
    parser.add_argument('--theta_sun', metavar='', default=0,choices=range(0,89), type=int, help='the theta_sun: parameter used in Albert 2004, the sun zenith angle in water in degree 0-90, default is 0, nadir') 
    parser.add_argument('--theta_view', metavar='', default=0,choices=range(0,89), type=int, help='the theta_view: parameter used in Albert 2004, the view zenith angle in water in degree 0-90, default is 0, nadir') 
    parser.add_argument('--wind_speed', metavar='', default=5, type=float, help='the the wind speed in units of (m/s), default is: 5') 
    
    
    ##from r_rs(λ) (under water ) to Rrs(λ) (above water)
    #Lee, Z.-P., Carder, K.L., Mobley, C.D., Steward, R.G., Patch, J.S., 1998. Hyperspectral remote sensing for shallow waters: 1. A semi-analytical model. Appl. Opt. 37, 6329-6338.
    ##eta: the water-to-air radiance divergence factor
    ##gamma: accounts for the effects of internal reflection from water to air
    ##Rrs_surf: the reflections at the water surface, ignore here temporarily
    #eta,gamma=0.52, 1.6
    parser.add_argument('--eta', metavar='', default=0.52, type=float, help='parameter for rrs--> Rrs in Lee et al 1998. eta: the water-to-air radiance divergence factor, default is: 0.52') 
    parser.add_argument('--gamma', metavar='', default=1.6, type=float, help='parameter for rrs--> Rrs in Lee et al 1998. gamma: accounts for the effects of internal reflection from water to air, default is: 1.6')                    
    
    #export result as file
    parser.add_argument('-S', '--save',metavar='',  type=str2bool, nargs='?',const=True, default=True, help='flag to save the simulation as two [_Rrs.csv, _rrs.csv] temp files, default is False')
    parser.add_argument('--out_prefix',metavar='',  default=None, help='the output csv files prefix, can be a full path. default is no prefix, will rewrite the two files in current folder [_Rrs.csv, _rrs.csv')
    
    #args = parser.parse_args(['-F','Rrs_QAA_test.csv'])
    args = parser.parse_args()
    
    print(f"==== input argument inspection:v_chl {args.chl},v_mspm={args.mspm} ,v_cdom={args.cdom},sDir_IOP={args.sIOP},IOPname={args.IOPname}, \
                                 wavelength={args.wavelength},bAB={args.bAB},f={args.f},bLee={args.bLee}, str_path_waterIOP={args.waterIOP},header_lines={args.header_lines},\
                                 water_bb2b_ratio={args.water_bb2b_ratio},chl_bb2b_ratio={args.chl_bb2b_ratio},mspm_bb2b_ratio={args.mspm_bb2b_ratio},\
                                 theta_sun={args.theta_sun},theta_view={args.theta_view},wind_speed={args.wind_speed},eta={args.eta},gamma={args.gamma} ")
                            
    Rrs,rrs=OpenLW_simu_Rrs_from_IOP(args.chl,args.mspm ,v_cdom=args.cdom,sDir_IOP=args.sIOP,IOPname=args.IOPname,
                                 wavelength=args.wavelength,bAB=args.bAB,f=args.f,bLee=args.bLee, str_path_waterIOP=args.waterIOP,header_lines=args.header_lines,
                                 water_bb2b_ratio=args.water_bb2b_ratio,chl_bb2b_ratio=args.chl_bb2b_ratio,mspm_bb2b_ratio=args.mspm_bb2b_ratio,
                                 theta_sun=args.theta_sun,theta_view=args.theta_view,wind_speed=args.wind_speed,eta=args.eta,gamma=args.gamma)    
    if (Rrs is None) or (rrs is None): 
        print("=====failed to run the simulation, no result is generated, will exit")
        sys.exit(0)  
    else:
        print("===== successfully finish the Rrs simulation ===========")
    
    ##export the simulation result as csv
    if  args.save:
        if not isinstance(args.out_prefix,type(None)):  
            prefix=args.out_prefix
        else: 
            prefix=r'./'
        Rrs.to_csv(prefix+r'_Rrs.csv',index=False)
        rrs.to_csv(prefix+r'_rrs_.csv',index=False)  #bug of to_csv:  Windows path case insensitive. need to use a different file name; https://github.com/pandas-dev/pandas/issues/19126
        
        print(f"===== written simulation result to: {prefix+'[_Rrs , _rrs_].csv'} ===========")
