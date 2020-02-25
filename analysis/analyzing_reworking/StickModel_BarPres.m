% StickModel_BarPres
%
%   Run avulsion model using rectangular channels and simple floodplain filling rules, then calculate channel preservation.
%   Preservation is calculated for all channel elements, and requires the Polygon Clip package on the Matlab exchange. Model output is saved as a .mat file in the last cell
%
%   preservation code written by EPC March 2017
%
%   light editing and wrapping by EAB Nov 2019

%% %%%%%% COMMONLY MODIFIED INPUT VARIABLES

function [mod_str] = StickModel_BarPres(in_str)

nits = in_str.nits;
mw = in_str.mw;
avul = in_str.avul;
ch_z = in_str.ch_z;
ch_w_i = in_str.ch_w_i;
IR = in_str.IR;
randint = in_str.randint;
vertagg_rate = in_str.vertagg_rate;
fp_lat = in_str.fp_lat;
walk_mag = in_str.walk_mag;
wwide = in_str.wwide;
wthick = in_str.wthick;
plot_bool = in_str.plot_bool;
latmob = in_str.latmob;
fp_filling = in_str.fp_filling;

rel_window = in_str.rel_window_mult;

% nits = 1000; %% number of iterations

init_z = 3;

% model and channel element dimensions
% mw = 100; %%%% model width multiplier [I've been using constant @ 50x channel width]
% ch_z=1; %sand body height (generic length units)
ch_w=nan(nits,1); % pre-allocate channel width vector
ch_w(1)=ch_w_i; %initial sand body width (generic length units) MUST BE AN ODD NUMBER=

% incision, floodplain aggradation, and avulsion parameters
% IR = 0.65; %%%% incision rate multiplier [actual incision at the base of the channel element every timestep is IR times channel thickness]
% randint = 1; %%% interval within which there is 1 random timestep
%%% for COMPENSATIONAL runs, set this to 10 (i.e. one random
%%% timestep every 10 timesteps, the rest are compensational)

% vertagg_rate = 0.25; %%% floodplain agg rate at each side of channel [depositional width can vary too; see next parameter]
% fp_lat = 4; %% the width range where fp is deposited on each side of channel, multipied by channel width
% lambda = 0.005; % exponential decay exponent for floodplain aggradation

%%% window dimensions (how much of the model do you want to sample for
%%% analysis? samples a window from the center to avoid edge effects
deswidth = wwide; % desired window width (scalar multiple of ch_w)
desthick = wthick; % desired window thickness (scalar multiple of ch_z)

% MODEL DOMAIN, TIMESTEPS, & MAXIMUM TOPOGRAPHIC RELIEF LIMIT
model_w = (ch_w(1)*mw); % width of model domain (must be larger than channel body width)
model_idx = 1:model_w; % index values across model domain
max_relief = 1*ch_z; % maximum permissible topographic relief, as a multiple of sand body height
% walk_mag=4; % multiple of the channel width that the random walk is allowed to move within


% AVULSION PATTERN
avuls_type = avul*ones(1,nits); % vector with flags for clustered (0), compensation (1), or random avulsions (2)

% %%%% OPTION 1: FIXED INTERVAL BETWEEN RANDOM/COMP
% % Random timesteps
% %randtime = 0; % number of steps between random avulsions
% randints = randtime:randtime:nits; % vector of timesteps when random avulsions will occur
% avuls_type(randints) = 2;

% Compensation timesteps
% %comptime = 1; % number of steps between compensational avlusions (compensation timescale)
% compints = comptime:comptime:nits; %vector of timesteps when compensation will occur
% avuls_type(compints) = 1;

%%%% OPTION 2: RANDOMIZED PLACEMENT OF AVULSION PATTERN
if randint > 0
	%     avuls_type = nan;
	randplac = randi(randint,[1 floor(nits/randint)]);
	%     randplac = nan(1,nits/randint);
	%     for ix =1:length(randplac)
	%         randplac(ix) = randi(randint,1);
	%     end
	for kx = 1:length(randplac)
		ff = (kx-1)*randint;
		avuls_type(ff + randplac(kx)) = 2;
	end
	
	%     for kx = 1:length(randplac)
	%         at = ones(1,randint); at(randplac(kx)) = 2;
	%         avuls_type = [avuls_type at];
	%     end
	%     avuls_type(1) = [];
end


% CHANNEL INCISION
% option 1: constant incision rate
incision_rate = (ch_z*(IR))*ones(nits,1); % currently defined as half of the sand body height

% LATERAL MOBILITY PATTERN
% lat_mob = zeros(1,nits); % pre-allocate a lateral mobility vector
% latmob=0;
if latmob ~= 0
	lat_dist = -latmob*ch_w(1):latmob*ch_w(1);  % the possible range of mobility distances, including left (neg.) and right (pos) directions, as a function of channel width
else
	lat_dist = nan;
end


% MODEL OUTPUT
Topo = nan(nits,model_w);       % Topo profile storage matrix (each row is topography at a given timestep)
Topo(1,:) = init_z;                  % initial condition for elevation is two
Topo(1,1:(ch_w(1)-1)/2) = init_z+ch_z;       % set buffer zone on lefthand side
Topo(1,model_w-(ch_w(1)-1)/2:end) = init_z+ch_z;       % set buffer zone on righthand side
% to add white noise to initial topography:
Tnoise = 0.1*vertagg_rate.*randn(1,model_w);     % generate vector of random noise to add
Topo(1,:) = Topo(1,:) + Tnoise;
%
% Topo_reference = NaN(nits,model_w);
% Topo_reference(1,:) = 0;
ch_idx = NaN(nits+1,1);         % vector for storing channel index location through time
ch_idx(1) = randi(model_w,1,1);  % choose a random initial channel location
chan_pos = nan(nits,4);         % channel position matrix
ch_center = nan(nits,2);          % channel centroids matrix
% chXY = NaN(nits+1,2);           % vector for storing channel centroid coordinates for SPP analysis
relief = NaN(nits,1);           % vector to store maximum relief at each timestep
fpagg_area = nan(nits,1);       % vector to store floodplain area deposited each timestep
fpagg_area_window = nan(nits,1);  % vector to store floodplain area deposited each timestep in evaluation window
% Randomly choose the lateral mobility of the channel for this timestep
% if isnan(lat_dist)==0
lat_mob = randsample(lat_dist,nits,true);
% end

% keyboard

%% Stratigraphic model
% figure;

for t = 2:nits
	
	ch_w(t) = ch_w(1)+abs(lat_mob(t)); % set a new channel width based on lateral mobility
	
	redge = model_w - ch_w(t) - 1;
	ledge = ch_w(t) + 1;
	%     keyboard
  %%% CALCULATE MOVING AVERAGE OF TOPOGRAPHY TO HELP COMP TIMESTEP
  %%% AVOID SPOTS RIGHT NEXT TO CHANNELS
  smoothtopo = movmean(Topo(t-1,:),(ch_w(1)-1) + latmob);
  %smoothtopo(1:ledge) = smoothtopo(1:ledge) + ch_z + relief(t-1);
  % smoothtopo(redge:end) = smoothtopo(redge:end) + ch_z + relief(t-1);
    
	%%%%%%%%%%%%%%%%% COMPENSTION TIMESTEPS
	if avuls_type(t-1) == 1
		
		%%% CALCULATE MOVING AVERAGE OF TOPOGRAPHY TO HELP COMP TIMESTEP
		%%% AVOID SPOTS RIGHT NEXT TO CHANNELS
		% smoothtopo = movmean(Topo(t-1,:),(ch_w(1)-1) + latmob);
		%smoothtopo(1:ledge) = smoothtopo(1:ledge) + ch_z + relief(t-1);
		% smoothtopo(redge:end) = smoothtopo(redge:end) + ch_z + relief(t-1);
		
		sortsmooth = sort(smoothtopo);
		
		%mins_idx = find(Topo(t-1,1:model_w-ch_w(1)) == min(Topo(t-1,1:model_w-ch_w(1)))); % find where the lowest points are across the model domain
		[~, mins_idx] = ismember(sortsmooth(1:100), smoothtopo); % find where the lowest points are across the model domain
		% clear smoothtopo
		
		% if length(mins_idx)>1
		%     disp(length(mins_idx));
		% end
		r_compind = randi(length(mins_idx)); % choose a random lowest point OR
		%[k,r_compind] = min(abs(mins_idx-ch_idx)); % choose the closest lowest point OR
		
		%         if mod(t,10)==0
		%         [k,r_compind] = max(abs(mins_idx-ch_idx(t-1))); % choose the furthest lowest point
		%         end
		%
		ch_idx(t) = mins_idx(r_compind); % new channel location for the current timestep
		
		% if ch_idx(t) > redge   % if the new channel index is too far to the right
		% 	ch_idx(t) = ledge + abs(model_w - ch_idx(t));
		% elseif ch_idx(t) < ledge  % if the new channel index is too far to the left
		% 	ch_idx(t) = redge - abs(model_idx(1) - ch_idx(t) - 1); % move the index to the right side of the model domain
		% end
		
		% if ch_idx(t) < ch_w(t)/2
		% 	keyboard
		% end
		
		%%%%%%%%%%%%%%%%% RANDOM AVLUSION CHANNEL LOCATIONS
		
	elseif avuls_type(t-1) == 2
		
		% Make sure new channel location doesn't exceed max relief at any
		% point along the width of the channel
		% WITH ELEVATION GOVERNER:
		% ch_idx_OK = find((Topo(t-1,1:model_w-ch_w(1)) + ch_z - incision_rate(t)) - min(Topo(t-1,1:model_w-ch_w(1))+vertagg_rate) < max_relief); %index values for locations where aggradation of ch_z won't exceed threshold
		topo_inbounds = smoothtopo(ledge:redge);
    topo_plus_agg = min(smoothtopo(ledge:redge)+vertagg_rate);
    potential_topo = topo_inbounds + ch_z - incision_rate(t) - topo_plus_agg;
    ch_idx_OK = find(potential_topo < max_relief); %index values for locations where aggradation of ch_z won't exceed threshold
		ind_buff = find(ch_idx_OK > ledge & ch_idx_OK < redge);
		ch_idx_OK = ch_idx_OK(ind_buff);
		ch_draw_OK = randi(length(ch_idx_OK)); %randomly select new channel location from model domain
		
		ch_idx(t) = ch_idx_OK(ch_draw_OK); %new ch_idx in model domain
		
		% adjust channel index location for lateral mobility direction
		ch_idx(t) = ch_idx(t) + floor(lat_mob(t)/2);
		
		% if ch_idx(t) > redge   % if the new channel index is too far to the right
		% 	ch_idx(t) = ledge + abs(model_w - ch_idx(t));
		% elseif ch_idx(t) < ledge  % if the new channel index is too far to the left
		% 	ch_idx(t) = redge - abs(model_idx(1) - ch_idx(t) - 1); % move the index to the right side of the model domain
		% end
		% CLUSTERED TIMESTEPS: RANDOM-WALK
		
	elseif avuls_type(t-1) == 0
		
		% Randomly choose how far to move the channel w/ random walk (as a function
		% of the DEFAULT channel width, ch_w); could be changed to be a
		% function of the channel width at that timestep, if desired
		r_possible = -(walk_mag)*(ch_w(1)-1):(walk_mag)*(ch_w(1)-1); % neg. is left, pos. is right
		rwalk = randsample(r_possible,1);
		ch_idx(t) = ch_idx(t-1) + rwalk;
		
		% adjust channel index location for lateral mobility direction
		ch_idx(t) = ch_idx(t) + floor(lat_mob(t)/2);
		
	end
	
	if ch_idx(t) > redge   % if the new channel index is too far to the right
		ch_idx(t) = ledge + abs(model_w - ch_idx(t));
	elseif ch_idx(t) < ledge  % if the new channel index is too far to the left
		ch_idx(t) = redge - abs(model_idx(1) - ch_idx(t) - 1); % move the index to the right side of the model domain
	end
	
	%%%%%%%%%%%%%%%%%%%% DETERMINE FLOODPLAIN AND CHANNEL DEPOSITION LOCATIONS
	
	%%%% determine channel indexes
	edge_buff = floor((ch_w(t)-1)/2);
	% S_left = :ch_idx(t)-1;
	% S_right = ch_idx(t):;
	S_range = ch_idx(t)-edge_buff:ch_idx(t)+edge_buff;
	
	S_range(S_range <= 0) = [];
	S_range(S_range >= model_w) = [];
	
	%     if S_range(1) > 5460
	%         keyboard
	%     end
	
	%%%% determine floodplain indexes
	switch fp_filling
		case 'everywhere' %%%%% OPTION 1: FLOODPLAIN EVERYWHERE BESIDES CHANNEL
			fp_range = setdiff(model_idx,S_range);  %index values of floodplain locations
		case 'local' %%%%% OPTION 2: FLOODPLAIN ONLY WITHIN X CHANNEL WIDTHS ON EITHER SIDE
			%fp_lat = 2;
			[total_range] = edgewrap(t,fp_lat,S_range,ch_w,model_idx,model_w);
			fp_range = setdiff(total_range,S_range);
	end
	
	%%%%%%%%%%%%%%%%%%% SET NEW TOPOGRAPHY
	%%%%%% FROM MSB MODEL
	% set new channel elevation
	Topo(t,ch_idx(t)) = smoothtopo(ch_idx(t)) + ch_z - incision_rate(t); %calculate at channel index location
	Topo(t,S_range)=Topo(t,ch_idx(t)); %set elevation for entire channel element based on channel index height
	
	%     chan_elev = Topo(t,ch_idx(t));
	switch fp_filling
		case 'local'
			fl_decay_x = 1:ceil(length(fp_range)/2);
			
			thresh = 0.1 * vertagg_rate;
			
			lambda = log(thresh * ch_z/vertagg_rate)/(length(fp_range)/2);
			
			fp_agg_exp_r = (vertagg_rate) * exp(lambda * fl_decay_x);
			fp_agg_exp_l = flip(fp_agg_exp_r);
			
			fp_agg_exp = [fp_agg_exp_l fp_agg_exp_r];
			
			fp_extent = length(fp_range);
			
			if mod(fp_extent, 2) ~= 0
				fp_agg_exp(end) = [];
			end
			
			lowerbound = (S_range(1) - fp_lat * ch_w(t));
			upperbound = (S_range(end) + fp_lat * ch_w(t));
			
% 			fp_agg_vec = zeros(model_w,1);
			
			if lowerbound < model_idx(1) 
				over = model_idx(1) - lowerbound;
				ledge = over;
				okbit = fp_agg_exp(ledge+1:fp_extent);
				overbit = fp_agg_exp(1:ledge);
				fp_agg_exp = [okbit overbit];
			elseif upperbound > model_w
				over = upperbound - model_w;
				redge = fp_extent - over;
				overbit = fp_agg_exp(redge+1:fp_extent);
				okbit = fp_agg_exp(1:redge);
				fp_agg_exp = [overbit okbit];
			end
						
			Topo(t,fp_range) = smoothtopo(fp_range) + fp_agg_exp;  %update fp locations by adding overbank aggradation rate
		case 'everywhere'
			Tnoise = 0.01 * vertagg_rate .* randn(1,length(fp_range));     % generate vector of random noise to add
			fp_agg_exp = vertagg_rate + Tnoise;
			Topo(t,fp_range) = Topo(t-1,fp_range) + fp_agg_exp;  %update fp locations by adding overbank aggradation rate
	end% update floodplain elevation
	
	%     Topo(t,fp_range) = Topo(t-1,fp_range)+fp_agg_exp;
	
	% does there need to be the exp. decay with lambda here?
	
	Topo(t,isnan(Topo(t,:))) = Topo(t-1,isnan(Topo(t,:)));
	
	% check to see if fp elevation exceeds max elevation
	if vertagg_rate>0
		max_elev = Topo(t,ch_idx(t));   %maximum aggradation elevation, the height of the active channel element
		toohigh_idx = fp_range(find(Topo(t,fp_range)>max_elev));  % find locations where fp height exceeds max
		
		for m=1:length(toohigh_idx)
			if Topo(t,toohigh_idx(m))-max_elev <= vertagg_rate
				Topo(t,toohigh_idx(m))= max_elev;
			else
				Topo(t,toohigh_idx(m)) = Topo(t-1,toohigh_idx(m));
			end
		end
		
		Topo(t,S_range)=Topo(t,ch_idx(t)); %update channel elevations
		
	end
	
	
	% calculate floodplain aggradation area for this timestep
	fpagg_area(t) = sum(Topo(t,fp_range)-Topo(t-1,fp_range));
	
% 	calculate floodplain aggradation in window 
	xind_low = floor(model_w * deswidth/2);
	xind_high = model_w - floor(model_w * deswidth/2);
	fp_range_window = fp_range(fp_range > xind_low & fp_range < xind_high);
	
	fpagg_area_window(t) = sum(Topo(t,fp_range_window)-Topo(t-1,fp_range_window));

	
	
	% update non-channel or fp cells
	% %           Topo(t,setdiff(model_idx((ch_w-1)/2:end),total_range)) = Topo(t-1,setdiff(model_idx((ch_w-1)/2:end),total_range));
	
	%     Topo(t,1:(ch_w-1)/2) = max(Topo(t,:)); % set buffer zone
	
	
	
	%%%%%%%%%%%%%%%%%%%% STORE CH AND DEP ELEMENT LOCATIONS
	% channel center
	ch_center(t,1)=ch_idx(t);
	ch_center(t,2)=(Topo(t,ch_idx(t))-(ch_z/2));
	
	% Channel position matrix for use with matlab fcn rectint
	chan_pos(t,1)=S_range(1);
	
	% FOR TOPO AGG OPTION 1:
	chan_pos(t,2)=Topo(t-1,ch_idx(t))-incision_rate(t);
	chan_pos(t,3)=ch_w(t)-1;
	chan_pos(t,4)=ch_z;
	
	% max relief calculation
	relief(t) = max(Topo(t,:))-min(Topo(t,:));
	
	% local relief calculation
	relief_window = rel_window; %how far to look on each side of the channel to calculate the relief
	relrange=edgewrap(t,relief_window,S_range,ch_w,model_idx,model_w);
	LR(t) = max(Topo(t,relrange))-min(Topo(t,relrange));
	clear relrange
	
	
	
	% % make the next avulsion compensational if the relief is beyond a threshold
	if t < nits
		if LR(t) >= max_relief
			avuls_type(t) = 1;
		end
	end
	
	% make the next avulsion random if the avulsion type is compensational and
	% there is no relief on the Topo surface
	if avuls_type(1) == 1 && t > 1
		if length(unique(Topo(t, 1:(end - ch_w(1))))) == 1
			avuls_type(t) = 2;
		end
	end
	
	
% 	keyboard 
	
	ledgetaper = floor(0.05 * model_w);
	redgetaper = floor(model_w - 0.05 * model_w);

	close_mean = abs(Topo(t,:) - nanmean(Topo(t,ledgetaper:redgetaper)));

	ledgetaper = find(close_mean < 0.02, 1, 'first');
	redgetaper = find(close_mean < 0.02, 1, 'last');

	lmean = nanmean(Topo(t,ledgetaper:ledgetaper + (5 * ch_w(1))));
	rmean = nanmean(Topo(t,redgetaper-(5 * ch_w(1)):redgetaper));

	edge_decay_xL = 1:ledgetaper;
	edge_decay_xR = 1:model_w - redgetaper;

	int_max = max(Topo(t,ledgetaper:redgetaper));
	int_mean = mean(Topo(t,ledgetaper:redgetaper));
	edgebound = int_max - int_mean;
	
	thresh = 0.1 * vertagg_rate;
	
	lambdaL = log(thresh./edgebound)./(length(edge_decay_xL));
	lambdaR = log(thresh./edgebound)./(length(edge_decay_xR));

	edge_agg_exp_l = (edgebound) .* exp(lambdaL .* edge_decay_xL);
	edge_agg_exp_r = (edgebound) .* exp(lambdaR .* edge_decay_xR);

	Topo(t, 1:ledgetaper) = edge_agg_exp_l + int_mean;
	Topo(t, redgetaper + 1:model_w) = flip(edge_agg_exp_r) + int_mean;
	
	Tnoise = 0.01 * vertagg_rate .* randn(1, model_w);     % generate vector of random noise to add
	Topo(t, :) = Topo(t, :) + Tnoise;
	
end

% keyboard

%% Calculate channel preservation
id = 'MATLAB:polyshape:repairedBySimplify';
warning('off', id)
id = 'MATLAB:polyshape:boolOperationFailed';
warning('off', id)

%
% Preallocate matrix to store percent preservation of every channel body at
% every timestep
p_pres = nan(nits,nits); % rows are timesteps, columns channel objects
p_erod = nan(nits,nits);
final_pres = nan(1,length(2:(size(p_pres,1))));
ch_touch = zeros(1,nits); % vector to store how many other channel objects touch this one
lat_erod = nan(1,length(2:(size(p_pres,1))));
vert_erod = 0*ones(1,length(2:(size(p_pres,1))));
vert_pres = ch_z*ones(1,length(2:(size(p_pres,1))));
fvp = ones(1,nits-1)*100; % percent of channel element width that is fully preserved vertically
et_mat = nan(nits,nits); % save all of the erosion times

% Calculate area of the rectangle at each timestep
ch_area=(ch_w-1)*ch_z; % channel width times channel height

% Determine overlap of rectangles between each timestep
for tt = 2:(size(p_pres,1)-1)   % time counter
	
	for ii = 1:(size(p_pres,2)-tt)  % channel object (rect) counter
		
		% when a subsequent rectangle erodes the entire thickness
		if chan_pos(tt+ii,2)<chan_pos(tt,2) && rectint(chan_pos(tt,:),chan_pos(tt+ii,:))~=0
			
			overlap_area = (length((intersect(chan_pos(tt,1):(chan_pos(tt,1)+chan_pos(tt,3)),chan_pos(tt+ii,1):(chan_pos(tt+ii,1)+chan_pos(tt+ii,3)))))-1) * ch_z;
			lat_erod(tt) = overlap_area/ch_z;
			vert_erod(tt) = ch_z;
			
			% when eroding channel does not erode below first channel
		else
			overlap_area = rectint(chan_pos(tt,:),chan_pos(tt+ii,:));
			if overlap_area > 0
				lat_erod(tt) = length((intersect(chan_pos(tt,1):(chan_pos(tt,1)+chan_pos(tt,3)),chan_pos(tt+ii,1):(chan_pos(tt+ii,1)+chan_pos(tt+ii,3)))))-1;
				vert_erod(tt) = overlap_area / lat_erod(tt);
			end
		end
		
		p_pres(ii,tt) = 100 - (overlap_area/ch_area(tt))*100; % percent of rectangle preserved at each timestep
		p_erod(ii,tt) = (overlap_area/ch_area(tt))*100; % percent eroded per timestep
	end
	
	% record when the erosion happened
	et = find(p_pres(:,tt)<100);
	eros_time = tt + et ;
	if ~isempty(eros_time>=1)
		et_mat(1:length(eros_time),tt)=eros_time;
	end
	toterod = sum(p_erod(et,tt)); % total erosion at that timestep
	ch_touch(tt)=length(et);
	
	
	% preallocate variables
	ch_int = nan(length(eros_time),4); % store intersection between eroding channel and original channel
	ov_ar = zeros(length(eros_time),1); % actual area eroded per timestep
	
	% find intersection of rectangles for each time
	for e = 1:length(eros_time)
		
		% Find the horizontal overlap of eroding channel with original
		ch_int(e,1) = min(intersect(chan_pos(tt,1):(chan_pos(tt,1)+chan_pos(tt,3)),chan_pos(eros_time(e),1):(chan_pos(eros_time(e),1)+chan_pos(eros_time(e),3))));
		ch_int(e,3) = (max(intersect(chan_pos(tt,1):(chan_pos(tt,1)+chan_pos(tt,3)),chan_pos(eros_time(e),1):(chan_pos(eros_time(e),1)+chan_pos(eros_time(e),3))))) - ch_int(e,1);
		
		% Find the vertical overlap of eroding channel with original
		if chan_pos(tt,2)>chan_pos(eros_time(e),2) % if overlying channel totally erodes
			ch_int(e,2) = chan_pos(tt,2); % channel totally intersects first channel
			ch_int(e,4) = chan_pos(tt,4);
		else
			ch_int(e,2) = chan_pos(eros_time(e),2); %if overlying channel partially erodes
			ch_int(e,4) = (chan_pos(tt,2)+ch_z) - chan_pos(eros_time(e),2);
		end
		
		%calculate the total area of erosion for each eroding channel
		ch_int(e,5) = ch_int(e,3) * ch_int(e,4);
		
		% find polygon vertices of each eroding channel
		xerod = ch_int(e,1):ch_int(e,1)+ch_int(e,3);
		yerod = linspace(ch_int(e,2),(ch_int(e,2)+ch_int(e,4)),length(xerod));
		Erod(e).x = [xerod fliplr(xerod)];
		Erod(e).y = [ch_int(e,2)*ones(1,length(xerod)) (ch_int(e,2)+ch_int(e,4))*ones(1,length(xerod))];
		Erod(e).hole = 0;
	end
	
	%%%%% CALCULATE PRESERVED AREA AND VERTICAL PRESERVATION
	if min(p_pres(:,tt)) == 100 %%%% if no channels erode
		final_pres(tt) = 100;
		
	elseif length(eros_time) == 1 %%%% if only one channel is eroding
		
		final_pres(tt) = min(p_pres(:,tt));
		
		% find percent of channel width that is fully vertically preserved
		fvp(tt) = (ch_w(tt)-1 - ch_int(e,3))/(ch_w(tt)-1)*100;
		
		% find max vertical preservation
		if fvp(tt) == 0
			vert_pres(tt) = ch_z - ch_int(4);
		end
		
		
	elseif length(eros_time) > 1 %%%% if multiple channels erode
		
		% find the unification polygon using PolygonClip
		for k=1:length(eros_time)-1
			if k==1
				clear place
				p1 = polyshape(Erod(k).x, Erod(k).y);
				p2 = polyshape(Erod(k+1).x, Erod(k+1).y);
				place = union(p1,p2);
			else
				p3 = polyshape(Erod(k+1).x, Erod(k+1).y);
				place = union(place,p3);
				%                 [Uni(1:size(place,2)).x] = place(1:end).x;
				%                 [Uni(1:size(place,2)).y] = place(1:end).y;
			end
			%             xx = place.Vertices(:,1);
			%             yy = place.Vertices(:,2);
			%             Uni(1:place.NumRegions).x = xx(~isnan(xx));
			%             [Uni(1:place.NumRegions).y] = yy(~isnan(yy));
			%             if size(Uni,2)>size(place,2)
			%                 Uni((size(place,2)+1):end)=[];
			%             end
		end
		
		polyout = regions(place);
		
		for ii=1:place.NumRegions
			Uni(ii).x = polyout(ii).Vertices(:,1);
			Uni(ii).y = polyout(ii).Vertices(:,2);
		end
		
		for b=1:size(Uni,2)
			
			%%% find maximum erosion for each unification polygon
			e_area(b) = polyarea(Uni(b).x,Uni(b).y);
			l_erod(b) = max(Uni(b).x) - min(Uni(b).x);
			v_erod(b) = max(Uni(b).y) - min(Uni(b).y);
			
			%%% find maximum vertical preservation for each
			%%% unification polygon
			
			
			%%% find corner vertices in unification polygon only
			
			%             ind = 2:length(Uni(b).x)-1;
			%             cornerind = nan(1,length(ind));
			%
			%             for d = 1:length(ind)
			%                 if Uni(b).x(ind(d)-1) ~= Uni(b).x(ind(d)) || Uni(b).x(ind(d)+1) ~= Uni(b).x(ind(d))
			%                     cornerind(ind(d)) = ind(d);
			%                 else
			%                     cornerind(ind(d)) = nan;
			%                 end
			%             end
			%
			%             testingy = Uni(b).y(cornerind(isnan(cornerind)==0));
			
			%             sorty = sort(testingy);
			%             unsorty = flipud(unique(sorty));
			
			
			%             v_pres(b) = ch_z - (unsorty(1) - unsorty(2));
			v_pres(b) = ch_z - v_erod(b);
			
		end
		
		final_pres(tt) = 100 - (sum(e_area)/ch_area(tt))*100;
		lat_erod(tt) = sum(l_erod); % how much of the element width is removed
		vert_erod(tt) = max(v_erod);
		
		if lat_erod(tt) == (ch_w(tt)-1)
			vert_pres(tt) = max(v_pres); %%%% maximum element thickness preserved
		end
		%
		%                 if isempty(vertp)==0
		%                 fvp(tt) = (length(vertp)+1)/(ch_w(tt)-1)*100;
		%                 else fvp(tt) = 0;
		%                 end
		
		clear Uni e_area Erod xerod yerod l_erod v_erod v_pres ind cornerind
		
		
	end
	
end

final_erosion = 100 - final_pres; % in percent removed per channel

% Find the location of partially eroded channels
erod_loc = find(final_pres ~= 100);

% Flag whether the vertical erosion is more than half the thickness at any
% point
vflag = zeros(1,length(2:(size(p_pres,1))));
vflag(vert_erod>(ch_z/2)) = 1;

xind_low = floor(model_w * deswidth/2);
xind_high = model_w - floor(model_w * deswidth/2);

ch_in_window = ch_idx > xind_low & ch_idx < xind_high;



% keyboard

%%%%%%%% CALCULATE NET AGGRADATION RATE
nagg = nanmean(Topo(end,:))/nits;

%%%%%%%%%% CALCULATE FINAL REWORKING

total_topo = nansum(Topo(end,xind_low:xind_high) - init_z);

total_fp = nansum(fpagg_area_window(ch_in_window));

total_ch = nansum(ch_area(ch_in_window));

reworking = total_topo/(total_fp + total_ch);



% % Bulk statistics on amount of preservation (final amounts)
% pres_stats_labels = [{'Mean'},{'Median'},{'Mode'},{'Min'},{'Max'},{'StDev'},{'% eroded'}]; % easy reference for pres_stats vector labels
% pres_stats(1) = mean(final_pres(~isnan(final_pres)));
% pres_stats(2) = median(final_pres(~isnan(final_pres)));
% pres_stats(3) = mode(final_pres(~isnan(final_pres)));
% pres_stats(4) = min(final_pres(~isnan(final_pres)));
% pres_stats(5) = max(final_pres(~isnan(final_pres)));
% pres_stats(6) = std(final_pres(~isnan(final_pres)));
% pres_stats(7) = (length(erod_loc)/length(final_pres(~isnan(final_pres))))*100;
%
% %%%% HISTOGRAM OF PRESERVATION RESULTS
% % endpres = find(final_pres(1:end-1)~=100,1,'last'); %don't consider final timesteps where channel is always preserved
% % bins=10:10:100;
% % preshist = hist(final_pres(1:endpres),bins);
% % figure; bar((preshist/length(final_pres(1:endpres)))*100);
%
% % To consider all timesteps:
% bins=0:10:100;
% preshist = hist(final_pres,bins);
% %figure; bar((preshist/length(final_pres(250:end)))*100);
% % set(gca,'XTickLabel',[0 10 20 30 40 50 60 70 80 90 100]);
% % xlabel('Preserved amount of channel (percent)')
% % ylabel('Percent of total channels')
%
% % Calculate percent of fully preserved channels at each timestep:
% cumper_EPC; % this uses final_pres and creates cumper100 and cumper0
%
% %

%% Analysis of vertical element preservation for bar preservation paper
%%%%%%%%% CLIP DATASET TO DESIRED PANEL WINDOW SIZE
%%% select data within central window

xind_low = floor(model_w * deswidth/2);
xind_high = model_w - floor(model_w * deswidth/2);

thickness = max(Topo(:)) - min(Topo(:));

yind_low = min(Topo(:)) + (thickness * desthick/2);
yind_high = max(Topo(:)) - thickness * desthick/2;

w_ind = find(chan_pos(1:end-1,1)>(xind_low) & chan_pos(1:end-1,1)<(xind_high));

w_ind(isnan(w_ind)==1) = [];
v_ind = find(chan_pos(1:end-1,2) > (yind_low) & chan_pos(1:end-1,2) < (yind_high));
v_ind(isnan(v_ind)==1) = [];

% issue error warning if there isn't enough thickness
if yind_high > max(Topo(:))
	disp('WARNING: WINDOW USES FULL MODEL THICKNESS');
end

ind1 = intersect(w_ind,v_ind); % find index of elements within the window
length_ind1 = length(ind1); % export number of elements within window

if isempty(ind1)
	error('Evaluation window is too small. No channels captured in window.')
end

%%%%% FULL PRESERVATION (PERCENT): HALF OF ELEMENT AREA PRESERVED AND MAX
%%%%% VERTICAL PRESERVATION > 90%
fullpres = length(find(vert_pres(ind1) > 0.99*ch_z & final_pres(ind1)>=50))/length(ind1)*100;
%fullpres = length(find(fvp(ind1) >= 25 & final_pres(ind1)>=50))/length(ind1)*100; % get PERCENT of elements with full vertical pres > 25% and preserved area >=50%
fullpres_ind = ind1(find(vert_pres(ind1) > 0.99*ch_z & final_pres(ind1)>=50));

%%%%% SIGNFICIANT PRESERVATION (PERCENT)
sigpres1 = length(find(vert_pres(ind1) <= 0.99*ch_z & vert_pres(ind1) > 0.75*ch_z & (final_pres(ind1)>=50)))/length(ind1)*100;
% max vert erosion >75% & preserved area > 50%
sigpres2 = length(find(vert_pres(ind1) >= 0.75*ch_z & final_pres(ind1)<50)) / length(ind1)*100;
% max vert erosion > 75% ch_z, but preserved area < 50%

sgpres1_ind = ind1(find(vert_pres(ind1) <= 0.99*ch_z & vert_pres(ind1) > 0.75*ch_z & (final_pres(ind1)>=50)));
sgpres2_ind = ind1(find(vert_pres(ind1) >= 0.75*ch_z & final_pres(ind1)<50));

% TRUNCATED (PERCENT)
truncated = length(find(vert_pres(ind1) < 0.75*ch_z))/length(ind1)*100;
truncated_ind = ind1(find(vert_pres(ind1) < 0.75*ch_z));

% PLOT PIE CHART OF RESULTS
% labels = {'full','sig1','sig2','trunc'};
% figure; pie([fullpres sigpres1 sigpres2 truncated],labels)

class = cell(size(vert_pres'));
class(fullpres_ind) = {'full'};
class(sgpres1_ind) = {'sig1'};
class(sgpres2_ind) = {'sig2'};
class(truncated_ind) = {'trunc'};

pres_tab = table(vert_pres(ind1)', final_pres(ind1)'/100, class(ind1));


% MAKE PRES INTO AN OUTPUT VECTOR
pres_percents = [fullpres sigpres1 sigpres2 truncated];
% keyboard
%% Calculate net-to-gross
figure('Visible', 'off');
hold on
for i=1:length(ind1)
	rectangle('Position',chan_pos(ind1(i),:),'FaceColor','k')
end
set(gca,'xlim',[min(chan_pos(ind1,1))+ch_w(1) max(chan_pos(ind1,1))],'Visible', 'off')
set(gca,'ylim',[min(chan_pos(ind1,2)) max(chan_pos(ind1,2))],'Visible', 'off')

axis off
F=getframe;
bw = im2bw(F.cdata);
n2g = length(find(bw==0))/(size(bw,1)*size(bw,2));
clf
% n2g = nan;

sand_in = nansum(ch_area(ind1));

input_n2g = sand_in / (sand_in + nansum(fpagg_area(ind1)));

% keyboard

% % %
mod_str = in_str;
mod_str.avuls_type = avuls_type ;
mod_str.chan_pos = chan_pos ;
mod_str.pres_percents = pres_percents ;
mod_str.n2g = n2g ;
mod_str.pres_tab = pres_tab;
mod_str.length_ind1 = length_ind1 ;
mod_str.fpagg_area = fpagg_area ;
mod_str.reworking = reworking ;
mod_str.nagg = nagg ;
if in_str.store_topo
	mod_str.Topo = Topo;
end
mod_str.model_w = model_w;
mod_str.xind_low = xind_low;
mod_str.xind_high = xind_high;
mod_str.yind_low = yind_low;
mod_str.yind_high = yind_high;
mod_str.relief = relief;
mod_str.input_n2g = input_n2g;

if plot_bool && in_str.store_topo
	plot_strat_str(mod_str, true, 25)
end

end
