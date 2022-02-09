function [ IDX ,IDX_slp ] = WaterIndex_from_spectra(index_bands,descr,simu_spectra,bandThreshold,bInterpolate)
%this function is used to define water index based on given spectra
%   Detailed explanation goes here
%build a Chl-index :  "MCI	R(2) - R(1) - (R(3) - R(1)).((2 - 1)/(3 - 1))	681, 708, 753"
% this function is compatible for 3-band line height indexes, since they share the same
% format of equation, but using different bands.
if nargin <4 || isempty(bandThreshold) 
    bandThreshold=10; %the threshold of difference between designed water index band, and the given sensor bands (nm)
    %bInterpolate=false; 
elseif nargin <5 || isempty(bInterpolate)
    bInterpolate=false;  %using interplotted Rrs from satellite bands, to calculate water index.
end

%generate water index
if bInterpolate
    val_bands=interp1(descr.nominal_centre_wavelength,simu_spectra',index_bands,'linear','extrap');
else
    [closest_bands,dist]=knnsearch(descr.nominal_centre_wavelength,index_bands');
    if sum( abs(dist)>bandThreshold ) || ... %there is large difference between water index bands and satellite badn central
            length(closest_bands)>length(unique(closest_bands))  % water index bands map to same satellite bands.
        IDX=[];   %if one sample is wrong, means this sensor is inappropirate for the given water index.
        IDX_slp=[];
        disp('==== failed to calcuate water index: distance(sensor band, water index band)>threshold, i.e., bands missing in sensor for the index ===');
        return;
    else
        val_bands= simu_spectra(:,closest_bands); %fetch the satellite Rrs at its bands directly (without interpolation)
    end
end
%calculate MCI
IDX=val_bands(:,2)-val_bands(:,1)-(val_bands(:,3)-val_bands(:,1))*(index_bands(2)-index_bands(1))./(index_bands(3)-index_bands(1));
IDX_slp=(val_bands(:,3)-val_bands(:,1))./(index_bands(3)-index_bands(1));

end

