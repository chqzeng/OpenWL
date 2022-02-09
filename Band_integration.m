function [ inter_Samples ] = Band_integration(descr, wavelength, spectra)
%spectra: given an input Rrs cruve group, as table of size(M,N) , where M is the
%number of samples, and N is the number of bands for each sample
%wavelength:  is a N*1 size vector of wavelength
%descr: the band description structure, at least including 'RSR' and 'RSR_wv' fields

    %descr=OLCI_descr; strtitle='RSR OLCI';
    nBands=length(descr.RSR(1,:));
    nSamples=length(spectra(:,1));
    inter_Samples=zeros(nSamples,nBands);  
    parfor id_smp=1:nSamples  %iterate for each sample
        T=1e-3;  %the minimum threshold of valid RSR to be condisered
        spectrum=spectra(id_smp,:);
        %inter_bands=zeros(nBands,1); %the bands of a given sensor, based on integration with RSR function
        for idx=1:nBands  %length(descr.RSR(1,:))
            if length(descr.RSR_wv(1,:))==1; %all of the bands share the common wavelength list
                wave_idx=1; 
            else
                wave_idx=idx; 
            end
            flt=descr.RSR(:,idx)>T & descr.RSR(:,idx) <= 1;
            valid_RSR=[descr.RSR_wv(flt,wave_idx) descr.RSR(flt,idx)];
            spec_interp=interp1(wavelength,spectrum,valid_RSR(:,1),'nearest', 'extrap'); %linear  nearest
            %apply the integration
            inter_Samples(id_smp,idx)=sum(spec_interp.*valid_RSR(:,2))/sum(valid_RSR(:,2));

            if mod(id_smp,1e5)==0
                disp(['Progress: ',num2str(id_smp),'/',num2str(nSamples)]);
            end
        end  
     end
 
end

