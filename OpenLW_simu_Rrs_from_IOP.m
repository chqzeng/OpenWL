function [ Rrs_table, Rrs ] = OpenLW_simu_Rrs_from_IOP(v_chl,v_mspm, varargin ) 
%this is the implementation of water-leaving raidance (Rrs) simulation as used in: 
%		Zeng C, Binding C. The effect of mineral sediments on satellite chlorophyll-a retrievals from line-height algorithms using red and near-infrared bands. Remote Sensing. 2019;11(19). doi:10.3390/rs11192306
%		Author: Chui Zeng;  chqzeng@gmail.com
%
%It uses a very simplified light and water interactcion empricical model:  Rrs = f *  bb / (a + bb);
%for a more strict simulation, please use HydroLight/Ecolight [Mobley, et al,https://www.sequoiasci.com/product/hydrolight/]
%
%Input arguments:
%		v_chl,v_cdom,v_mspm: the value list of chlorophyll-a[mg/m^3], cdom[1/m], and mineral content/mspm[g/m^3]; default cdom is: 0.994 (1/m)
%		sDir_IOP:  the directory where the IOP files are stored; default is the "./IOPfiles"
%		IOPname  : the name ID key for a provided IOP set. corresponding to the files in: `./IOPfiles/` folder; default IOP model:" "
%		wavelength: the list of wavelength to simulate the output Rrs; default: [400:5:850] nm
%		bAB: flag to use Chl absorption with model: abs(chl) = A x chl^(1-B) model, otherwise use simpler: abs(chl)= a* x chl   model
%       bSave: save the simulation result in current directory, in format: 'OpenLW_simu_Rrs_from_IOP_[%s]_chl[%d]_mspm[%d]_cdom[%f].csv'
%		f:  the emperical value f in formula: Rrs = f *  bb / (a + bb);  default use 0.319 after Jerome et al. (1988). 
%		waterIOP: the file name of water IOP, in the './IOPfiles/'	folder, defualt is: 'waterIOP_sea_water.txt'  for clear fresh water
%
%Output arguments: [Rrs_table,Rrs]
%		Rrs: the water leaving radiance (sr-1); dimensions: N x M; where
%            N=length(v_chl)*length(v_mspm)*length(v_cdom),  M=length(wavelength)
%       Rrs_table: N x (M+3); the table format of Rrs, while add 3 columns at
%            beginning for: chl, mspm, and cdom

	%setup the input argument parser
	[v_cdom,sDir_IOP,IOPname,wavelength,bAB,f,waterIOP]=deal([0.994],'./IOPfiles','LELW',400:5:850,false, 0.319,'waterIOP_sea_water.txt');
	p = inputParser;
	vld = @(x) isnumeric(x) && isvector(x) && all(x > 0); %a valid Scalar Positive vector
	vnum=@(x) isnumeric(x) && isscalar(x) && (x > 0);%a valid Scalar Positive Number
	addRequired(p,'v_chl',vld);
	addRequired(p,'v_mspm',vld);
	addOptional (p,'v_cdom',v_cdom,vld);
    addParameter(p,'sDir_IOP',sDir_IOP,@ischar);
	addParameter(p,'IOPname',IOPname,@ischar);
	addParameter(p,'wavelength',wavelength,vld);
	addParameter(p,'bAB',bAB,@islogical);
    addParameter(p,'bSave',true,@islogical);
	addParameter(p,'f',f,vnum);
	addParameter(p,'waterIOP',waterIOP,@ischar);
	addParameter(p,'water_bb2b_ratio',0.5,vnum);
	addParameter(p,'chl_bb2b_ratio',0.01,vnum);
	addParameter(p,'mspm_bb2b_ratio',0.02,vnum);
    parse(p,v_chl,v_mspm,varargin{:});
    
	%update the input arguments
	sDir_IOP=p.Results.sDir_IOP;
    IOPname=p.Results.IOPname;
    wavelength=p.Results.wavelength;
	
    Rrs=[];
	%load pure water absroption and scatter
	IOP_water=readIOP(fullfile(sDir_IOP,p.Results.waterIOP),3);
	if isempty(IOP_water); return;  end;
    % ----- these two scatter function are not used currently, because direclty read from files
    %purewater_scatter=@(lampta)16.06*1.08e-4*(500./lampta).^4.32;  %purewater scatter upon wavelength
    %purewater_phase=@(phi) 0.06225*(1+(cos(phi)).^2);  %the scatter phase function
    IOP_W_a=interp1(IOP_water(:,1),IOP_water(:,2),wavelength,'linear','extrap');
    IOP_W_b=interp1(IOP_water(:,1),IOP_water(:,3),wavelength,'linear','extrap');
    IOP_W_bb=IOP_W_b * p.Results.water_bb2b_ratio;  %currently set the bb/b ratio as: 0.5 for pure water
    
	%configure the correct absorption IOP files
    if p.Results.bAB %using the abs(chl) = A x chl^(1-B) model
        str_abs_path{1}=fullfile(sDir_IOP,[IOPname,'_abs_coeff_chl_AB.txt']);  
    else   % using the abs(chl)= a* x chl   model
        str_abs_path{1}=fullfile(sDir_IOP,[IOPname,'_abs_coeff_chl.txt']);    
    end
    str_abs_path{2}=fullfile(sDir_IOP,[IOPname,'_abs_coeff_CDOM.txt']);
    str_abs_path{3}=fullfile(sDir_IOP,[IOPname,'_abs_coeff_mineral.txt']);
    
	%load mass-specific absoprtion of chl, cdom, and mspm
    if p.Results.bAB
        IOP_temp=readIOP(str_abs_path{1},3);  %read chl 
		if isempty(IOP_temp); return;    end;
        %IOP_temp(IOP_temp(:,1)>720,3)=0;  % set the B value for 720nm -800nm to zero; (according to Mobley's book) 
        IOP_CHL_AB=interp1(IOP_temp(:,1),IOP_temp(:,[2 3]),wavelength,'linear','extrap');
    else
        IOP_temp=readIOP(str_abs_path{1});  %read chl 
		if isempty(IOP_temp); return;   end;
        IOP_CHL_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    end
    IOP_temp=readIOP(str_abs_path{2});  %read CDOM 
	if isempty(IOP_temp); return;   end;
    IOP_CDOM_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    IOP_temp=readIOP(str_abs_path{3});  %read mspm 
	if isempty(IOP_temp); return;   end;
    IOP_MSPM_a=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    
	%load mass-specific scatter properties of chl and mspm
    str_sct_path=fullfile(sDir_IOP,strcat(IOPname, '_scatter_coeff_',{'chl.txt', 'mineral.txt'}));    
    IOP_temp=readIOP(str_sct_path{1});  %read chl mass-specific scatter
	if isempty(IOP_temp); return;   end;
    IOP_CHL_b=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');    
    IOP_temp=readIOP(str_sct_path{2});  %read mineral mass-specific scatter
	if isempty(IOP_temp); return;   end;
    IOP_MSPM_b=interp1(IOP_temp(:,1),IOP_temp(:,2),wavelength,'linear','extrap');
    %coef=lsqcurvefit(b_over_wave,[-0.01 0.1],IOP_temp(:,1),IOP_temp(:,2),[-10 0],[0 10]); %lsqcurvefit(fun,x0_init,xdata,ydata,lowerbound,upperbound)
    %IOP_MSPM_b=b_over_wave(coef,wavelength); %fit a curve?
	
	%bb/b ratio : from scatter to back-scatter
    IOP_CHL_bb=IOP_CHL_b*p.Results.chl_bb2b_ratio;   
    IOP_MSPM_bb=IOP_MSPM_b*p.Results.mspm_bb2b_ratio; 
    
    %mix the chl, mspm,cdom to different combinations
	[an, bn, cn] = ndgrid(p.Results.v_chl,p.Results.v_mspm,p.Results.v_cdom); 
	IOP_list = [an(:), bn(:), cn(:)];
    num_chl=length(p.Results.v_chl)*length(p.Results.v_mspm)*length(p.Results.v_cdom); %length(IOP_list)
    Rrs=zeros(num_chl,length(wavelength));
    calc_absorp=@(chl,A,B)(A.*chl.^(-1*B));  %define model of absorption calcuation: abs(chl) = A x chl^(1-B)
	
	%parpool(8)  %parallel the porcessing
    %parfor idx=1:num_chl  
    if p.Results.bAB  
        for idx=1:length(v_chl) %need iteration for each [chl], as different IOP_chl_a with varied on chl 
                IOP_CHL_a=calc_absorp(v_chl(idx),IOP_CHL_AB(:,1),IOP_CHL_AB(:,2));  %calcuate the absroption of Chl dependent on chl
                IOP_CHL_a(IOP_CHL_a<0)=0;
                IOP_CHL_a=IOP_CHL_a';
            %matrix operation over wavelength
            Rrs(idx:length(v_chl):end,:)=simplified_Rrs(IOP_list(idx:length(v_chl):end,:),IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f);
        end
    else 
        Rrs=simplified_Rrs(IOP_list,IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f);
    end
    
    Rrs_table=[IOP_list,Rrs];  %add extra info to the Rrs
    Rrs_table=array2table(Rrs_table,'VariableNames',[{'Chl_mg_m3','MSPM_g_m3','CDOM_1_m'},strcat({'b'},num2str(wavelength'))']);
    if p.Results.bSave
        %csvwrite('test.csv',Rrs);
        out_file=sprintf('OpenLW_simu_Rrs_from_IOP_%s_chl%.0f_mspm%.0f_cdom%.3f.csv',p.Results.IOPname, max(p.Results.v_chl),max(p.Results.v_mspm),max(p.Results.v_cdom));  
        writetable(Rrs_table, out_file, 'WriteVariableNames', true) 
		disp(['====finish Rrs simulation, result written to:',out_file,'===='])
    end
end

function Rrs=simplified_Rrs(IOP_list,IOP_W_a,IOP_CHL_a,IOP_CDOM_a,IOP_MSPM_a,IOP_W_bb,IOP_CHL_bb,IOP_MSPM_bb,f)
% the core function of  Rrs calcuation
    num_chl=size(IOP_list,1);
    a_total= repmat(IOP_W_a,num_chl,1) + IOP_list(:,1)*IOP_CHL_a  + IOP_list(:,2)*IOP_MSPM_a + +IOP_list(:,3)*IOP_CDOM_a;
    bb_total=repmat(IOP_W_bb,num_chl,1)+ IOP_list(:,1)*IOP_CHL_bb + IOP_list(:,2)*IOP_MSPM_bb;
    Rrs = f .*  bb_total ./ (a_total+bb_total) ./pi;  %conversion from LU/ED to Rrs, as a ratio of Pi
end

function IOP= readIOP(str_in_file, columns)  %read the IOP from file
%load a IOP file
    IOP=[];
    if exist(str_in_file, 'file') ~= 2 % File not exists.
        disp(['the IOP file: ',str_in_file, 'does not exist! program will exit'])
        return ;
    end
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