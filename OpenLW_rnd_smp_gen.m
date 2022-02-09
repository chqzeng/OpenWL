% STEP 0:  joint [Chl, MSPM] distribution based on in-situ data is prepared     

function Joint_rnd_pts=OpenLW_rnd_smp_gen(varargin)	 
	 %------generate random joint distribution of water IOPs [chl, mspm, cdom] ----------------
	 [N,min_chl,gap_chl,max_chl,min_mspm,gap_mspm,max_mspm,min_cdom,gap_cdom,max_cdom]=deal(100,0,0.05,150,0,0.1,100,0.994,0.994,0.994);
	%default: N=100; % #number of random IOP samples
	%default: max_chl=150;  %range of Chlorophyll-a (mg/m^3): [min_chl:gap_chl:max_chl]
    %default: max_mspm=100; %range of sediments (g/m^3): [min_mspm:gap_mspm:max_mspm]
	%default: max_cdom=0.994;%CDOM is not tested in this study, an average from in-situ data is used. can create grid just as the other two params
	
	%grid_chl=0.1:0.05:150;  %range of Chlorophyll-a (mg/m^3)
    %grid_mspm=0.1:0.1:100; %range of sediments (g/m^3)
    %N_type = {'chl','mspm','cdom'};

    p = inputParser;
    vld = @(x) isnumeric(x) && isscalar(x) && (x > 0);  %a valid Scalar Positive Number
    %addRequired(p,'N',validScalarPosNum);
    addParameter(p,'N',N,vld);
    addParameter(p,'cdom',cdom,vld);
    addParameter(p,'max_chl',max_chl,vld);
    addParameter(p,'max_mspm',max_mspm,vld);
    addParameter(p,'min_chl',min_chl,vld);
    addParameter(p,'min_mspm',min_mspm,vld);
    addParameter(p,'gap_chl',gap_chl,vld);
    addParameter(p,'gap_mspm',gap_mspm,vld);
    parse(p,varargin{:});

	grid_chl=p.Results.min_chl : p.Results.gap_chl : p.Results.max_chl;
	grid_mspm=p.Results.min_mspm : p.Results.gap_mspm : p.Results.max_mspm;
	grid_cdom=p.Results.min_cdom : p.Results.gap_cdom : p.Results.max_cdom;
    [gridx,gridy,gridz]=meshgrid(grid_chl,grid_mspm,grid_cdom);
	
    sf = fit(JointDist_manual(:,1:2),JointDist_manual(:,3),'linear'); %'cubicinterp' 'linearinterp' %interplation of the histrogram field
    Join_grid=reshape(sf(reshape(gridx,[],1), reshape(gridy,[],1)),length(grid_mspm),length(grid_chl));
    
    [Joint_rnd_pts(:,1), Joint_rnd_pts(:,2)]=Generate_rnd_IOP_params(grid_chl,grid_mspm,Join_grid,N);
end