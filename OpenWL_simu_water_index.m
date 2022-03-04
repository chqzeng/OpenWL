%%--------------------------------------------------------------------------------------------------------------------------
%this script applies the Rrs simulation with given water IOPs, over a specific sensor
%it then computes line-height water index, and plot the water index vs pre-defined chl concnetration.
%	Input arguments:
%				v_chl,v_mspm,v_cdom:  vectors of IOPs. 
%                    defualt: 5:5:100  2:2:20  0.994, respectively
%				wavelength: 			wavelength of simulation, default  400:800 at 1nm resolution
%				bPlot: whether plot the result of water idx after simulation, default is true
%				bLog: whether plot with log(Chl)
%				bSave: whether save the integrated bands and water index
%				sensor: sensor selection from list: ['HYPER','OLCI','MSI','MERIS','OCI','VIIRS','OLI','MODIS_Aqua','MODIS_Terra'] for simulation
%				waterIDX: line height water index selection: ['MCI','FLH','CI','user']; if choose 'user' model, then the [left, center,right] bands need to set for "designed_bands"
%				designed_bands: user-defined 3 bands for water index calculation, not used if a waterIDX predefined model is chosen. default is [681 708 753].
%				Band_diff_tor: tolerance of index designed band and satellite band difference, default 20nm
%--------------------------------------------------------------------------------------------------------------------------
function Simu_water_index(varargin)
p = inputParser;
vld = @(x) isnumeric(x) && isvector(x) && all(x > 0); %a valid Scalar Positive vector
vnum=@(x) isnumeric(x) && isscalar(x) && (x > 0);%a valid Scalar Positive Number
addParameter(p,'v_chl',5:5:100,vld);  %concnetraion of chl, mspm, and cdom, same as simulation function
addParameter(p,'v_mspm',2:2:20,vld);
addParameter(p,'v_cdom',0.994,vld);
addParameter(p,'wavelength',400:800,vld);  %wavelength of simulation
addParameter(p,'bPlot',true,@islogical);     %whether plot the result of water idx after simulation, default is true
addParameter(p,'bSave',true,@islogical);     %whether save the integrated bands and water index
addParameter(p,'sensor','OLCI',@(x) ischar(x) &&  any(validatestring(x,{'HYPER','OLCI','MSI','MERIS','OCI','VIIRS','OLI','MODIS_Aqua','MODIS_Terra'})));  %selection of sensor for simulation
addParameter(p,'waterIDX','MCI',@(x) ischar(x) && any(validatestring(x,{'MCI','FLH','CI','user'})));    %selection of line-height index for calculation
addParameter(p,'designed_bands',[681 708 753],vld);  %user-defined 3 bands for water index calculation, will overwrite the `waterIDX` if set
addParameter(p,'Band_diff_tor',20,vnum);   %tolerance of index designed band and satellite band difference, default 20nm
addParameter(p,'bLog',false,@islogical);   %whether plot with log(Chl)
parse(p,varargin{:});
wavelength=p.Results.wavelength;  %400:800;  %1nm resolution, too much for general app. except here for sensor band integration

% load band description of existing sensors: ['OLCI','MSI','MERIS','OCI','VIIRS','OLI','MODIS_Aqua','MODIS_Terra'], add new sensor to this file if needed
descr=load('Sensor_RSR.mat');    
 
%run a simulation
% OpenWL_simu_Rrs_from_IOP options, e.g.: ['waterIOP','waterIOP_SmithandBaker.txt','IOPname','LELW','wavelength',400:100:800,'bAB',false,'bSave',true]
[Simu_Rrs,simu_spectra]=OpenWL_simu_Rrs_from_IOP(p.Results.v_chl,p.Results.v_mspm,p.Results.v_cdom,'wavelength',wavelength,'bSave',false); 
Chl=Simu_Rrs.Chl_mg_m3;

%choose the sensor description 
try
	if ~strcmp(p.Results.sensor, 'HYPER')
		descr_sensor=descr.(sprintf('descr_%s', p.Results.sensor));
	else
		descr_sensor.nominal_centre_wavelength=wavelength'; %build a temp sensor model
	end
catch
	fprintf('====failed to load the sensor band description: descr_%s. program will exit', p.Results.sensor);
	return;
end
%choose the water index bands
if ~strcmp(p.Results.waterIDX, 'user')  %not user-defined water index, then select pre-defined index
	if strcmp(p.Results.waterIDX,'MCI'); designed_bands=[681 708 753]; end
	if strcmp(p.Results.waterIDX,'CI'); designed_bands=[665, 681, 708]; end
	if strcmp(p.Results.waterIDX,'FLH'); designed_bands=[665, 674, 753]; end
	
	%add more choice of water index here, and update the input parameter: addParameter(p,'waterIDX',...)
else
	designed_bands=p.Results.designed_bands;
end

%%--------------simulate varied senosrs ---------------------
Band_diff_tor=p.Results.Band_diff_tor;
if ~strcmp(p.Results.sensor, 'HYPER')
	%generate the band integration for a given sensor
	inter_bands=Band_integration(descr_sensor,wavelength, simu_spectra);
	%calculate water index from the satellite bands after simu integration
else  %HYPER spectral 
	inter_bands=simu_spectra;  %HYPER spectral modelling, no need to integrate to sensor bands
	Band_diff_tor=2;
	%waterIDX = WaterIndex_from_spectra(designed_bands,descr_sensor,simu_spectra,2);
end	

waterIDX = WaterIndex_from_spectra(designed_bands,descr_sensor,inter_bands,Band_diff_tor); 	
if isempty(waterIDX); return ; end;  %if failed to calcuate the water index
if strcmp(p.Results.waterIDX,'CI'); waterIDX=-1*waterIDX; end;  %CI needs to reverse the symbol

%save the simulation result
if p.Results.bSave
    Rrs_table=array2table([inter_bands,waterIDX],'VariableNames',[cellstr(num2str(descr_sensor.nominal_centre_wavelength,'b%3.f'))',{p.Results.waterIDX}]);
	Rrs_table=[Simu_Rrs(:,1:3),Rrs_table];
    out_file=sprintf('./TestData/OpenWL_simu_waterIDX_%s_of_sensor_%s.csv',p.Results.waterIDX,p.Results.sensor);  
    writetable(Rrs_table, out_file, 'WriteVariableNames', true) 
	disp(['====finish sensor water index simulation, result written to:',out_file,'===='])
end

if p.Results.bPlot
	figure; 
	%color=zeros(length(Chl),3);
	%color(:,1)=Simu_Rrs.MSPM_g_m3 ./max(Simu_Rrs.MSPM_g_m3);
	color=Simu_Rrs.MSPM_g_m3;
    if p.Results.bLog
        scatter(log(Chl),waterIDX,50,color,'filled');   %log(Chl)  pointsize:50
        xlim([prctile(log(Chl),0.5) prctile(log(Chl),99.5)]); xlabel('Predefined log Chl (\mug/L)');
    else
        scatter(Chl,waterIDX,50,color,'filled');   %log(Chl)  pointsize:50
        xlabel('Predefined Chl (\mug/L)');
    end
	colormap(jet);chb = colorbar(); caxis([min(Simu_Rrs.MSPM_g_m3),max(Simu_Rrs.MSPM_g_m3)]);ylabel(chb, 'MSPM(g/m^3)');
	ylabel(sprintf('Simulated %s (Sr^{-1})',p.Results.waterIDX)); 
	
	grid on; 
	title(sprintf('    waterIDX [%s] simulated for sensor %s',p.Results.waterIDX,p.Results.sensor))
end
end
