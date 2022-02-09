function [ Rrs ] = Simplified_a_bb_2_Rrs_fit_IOP_distr(str_name_key,v_chl,v_cdom,v_mineral, wavelength,bAB,f)
%this is the implementation of a simplified model to exam the IOPs.
%  Rrs = f *  bb / (a + bb);;   f = 0.319 after Jerome et al. (1988).
    if nargin < 6
        bAB=false;  %whether use the a*chl= A* Chl^B model (default is a*chl= constant(lampta) )
    end
    if nargin < 7 || ~exist('f','var') || isempty(f)  %if not defined f
        f=0.319;  %default f = 0.319 after Jerome et al. (1988)
    end
    str_IOP_avg='./IOPfiles/'  ; %'./HE52_IOPs_avg/'; %'./HE52_IOPs_avg/';  %average IOP generated from all samples,  Absorption_regression

    % ----- these two scatter function are not used currently, because direclty read from files
    %purewater_scatter=@(lampta)16.06*1.08e-4*(500./lampta).^4.32;  %purewater scatter upon wavelength
    %purewater_phase=@(phi) 0.06225*(1+(cos(phi)).^2);  %the scatter phase function
    IOP_water=readIOP([str_IOP_avg,'waterIOP_sea_water.txt'],3);
    %IOP_water=readIOP([str_IOP_avg,'Caren_H2OabClearNat.txt'],3);
    
    %read the pure water absorption and scattering
%     fileID = fopen([str_IOP_avg,'H2OabClearNat.txt']);
%     IOP_water=textscan(fileID,'%f %f %f','CollectOutput',true, 'HeaderLines',10) ; %skip the 10 lines as header, and read the rest as 
%     IOP_water=IOP_water{1,1};  %the pure water IOP
%     fclose(fileID);
    
    %set the pure water IOPs
    IOP_W_a=interp1(IOP_water(:,1),IOP_water(:,2),wavelength,'linear','extrap');
    IOP_W_b=interp1(IOP_water(:,1),IOP_water(:,3),wavelength,'linear','extrap');
    IOP_W_bb=IOP_W_b * 0.5;  %currently set the bb/b ratio as: 0.5 for pure water
    
    %str_abs_path=strcat(str_IOP_avg,str_name_key,'_abs_coeff_', {'chl','CDOM','mineral'}, '.txt');
    if bAB
        str_abs_path{1}=[str_IOP_avg,str_name_key,'_abs_coeff_chl_AB.txt'];  %using the abs(chl) = A x chl^(1-B) model
    else
        str_abs_path{1}=[str_IOP_avg,str_name_key,'_abs_coeff_chl.txt'];     % using the abs(chl)= a* x chl   model
    end
    str_abs_path{2}=[str_IOP_avg,str_name_key,'_abs_coeff_CDOM.txt'];
    str_abs_path{3}=[str_IOP_avg,str_name_key,'_abs_coeff_mineral.txt'];
    
    if bAB
        IOP_temp=readIOP(str_abs_path{1},3);  %read chl mass-specific absoprtion
        %IOP_temp(IOP_temp(:,1)>720,3)=0;  % set the B value for 720nm -800nm to zero; (according to Mobley's book) 
        IOP_CHL_AB=interp1(IOP_temp(:,1),IOP_temp(:,[2 3]),wavelength,'linear','extrap');
    else
        IOP_temp=readIOP(str_abs_path{1});  %read chl mass-specific absoprtion
        IOP_CHL_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    end
    IOP_temp=readIOP(str_abs_path{2});  %read CDOM mass-specific absoprtion
    IOP_CDOM_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    IOP_temp=readIOP(str_abs_path{3});  %read mineral mass-specific absoprtion
    IOP_MSPM_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    
    %read the scatter : only use the average coeffs currently
    str_sct_path=strcat(str_IOP_avg,str_name_key, '_scatter_coeff_',{'chl.txt', 'mineral.txt'});
    %str_sct_path=strcat(str_IOP_avg, {'LELW_scatter_coeff_chl.txt', 'LELW_scatter_coeff_mineral.txt'});
    %str_sct_path=strcat('D:\VboxShare\EC_work_2017\Absorption_regression\HE52_IOPs_avg\','Caren_scatter_coeff_', {'chl','mineral'}, '.txt');
    
    IOP_temp=readIOP(str_sct_path{1});  %read chl mass-specific absoprtion
    IOP_CHL_b=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    %opts = optimset('Display','off');
    %b_over_wave=@(a,x)a(1)*x+a(2) ; %fit a basica linear equation:  b = k * wavelength + b
    %coef=lsqcurvefit(b_over_wave,[-0.01 0.1],IOP_temp(:,1),IOP_temp(:,2),[-10 0],[0 10],opts); %lsqcurvefit(fun,x0_init,xdata,ydata,lowerbound,upperbound)
    %IOP_CHL_b=b_over_wave(coef,wavelength);
    
    IOP_temp=readIOP(str_sct_path{2});  %read mineral mass-specific absoprtion
    IOP_MSPM_b=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    %coef=lsqcurvefit(b_over_wave,[-0.01 0.1],IOP_temp(:,1),IOP_temp(:,2),[-10 0],[0 10]); %lsqcurvefit(fun,x0_init,xdata,ydata,lowerbound,upperbound)
    %IOP_MSPM_b=b_over_wave(coef,wavelength);
    
    
    IOP_CHL_bb=IOP_CHL_b*0.01;   %bb/b ratio is 0.02
    IOP_MSPM_bb=IOP_MSPM_b*0.02; %bb/b default is 0.02
    
    %do the calculation
    num_chl=length(v_chl);
    Rrs=zeros(num_chl,length(wavelength));
    calc_absorp=@(chl,A,B)(A.*chl.^(-1*B));  %define model of absorption calcuation: abs(chl) = A x chl^(1-B)
    for idx=1:num_chl  %iterate each chl value
        if bAB
            IOP_CHL_a=calc_absorp(v_chl(idx),IOP_CHL_AB(:,1),IOP_CHL_AB(:,2));  %calcuate the absroption of Chl dependent on chl
            IOP_CHL_a(IOP_CHL_a<0)=0;
            IOP_CHL_a=IOP_CHL_a';
        end
        Rrs(idx,:)=simplified_Rrs(v_chl(idx),v_cdom,v_mineral(idx),IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f) ;
    end
end

function Rrs=simplified_Rrs(v_chl,v_cdom,v_mineral,IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f)
    a_total= IOP_W_a + v_chl.*IOP_CHL_a +v_cdom.*IOP_CDOM_a + v_mineral.*IOP_MSPM_a;
    bb_total=IOP_W_bb + v_chl.* IOP_CHL_bb + v_mineral.*IOP_MSPM_bb;
    
    Rrs = f .*  bb_total ./ (a_total+bb_total) ./pi;  %conversion from LU/ED to Rrs, as a ratio of Pi
end
% function Rrs=simplified_Rrs2(IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f)%call the func without concetnration
%     a_total= IOP_W_a + IOP_CHL_a +IOP_CDOM_a + IOP_MSPM_a;
%     bb_total=IOP_W_bb + IOP_CHL_bb + IOP_MSPM_bb;
%     
%     Rrs = f .*  bb_total ./ (a_total+bb_total);
% end
function IOP= readIOP(str_in_file, columns)  %read the IOP from file
    if ~exist('columns','var') || isempty(columns)
        columns=2; %deault columns is 2 in the input IOP files (first column wavelength, 2nd as value);
    end
    
    fileID = fopen(str_in_file);
    if columns ==2
        IOP=textscan(fileID,'%f %f','CollectOutput',true, 'HeaderLines',10,'CommentStyle','-1.0') ;  %rows '-1.0 means eand'
    else
        IOP=textscan(fileID,'%f %f %f','CollectOutput',true, 'HeaderLines',10,'CommentStyle','-1.0') ; %skip the 10 lines as header, and read the rest as 
    end
    fclose(fileID);
    
    IOP=IOP{1,1};  %compress the IOP
end