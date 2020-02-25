%%%%%%%%%%%%%%%%%%%%%%%%%% edgewrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% last updated 03/28/2013
% by Ellen P. Chamberlin, based on code from E.A. Hajek
%
%
% This Matlab code is a function file embedded in StickModel_EPC that
% addresses the edge effects in the model.
%
%

function [total_range] = edgewrap(t,fp_lat,S_range,ch_w,model_idx,model_w)

lowerbound = (S_range(1) - fp_lat * ch_w(t));
upperbound = (S_range(end) + fp_lat * ch_w(t));

total_range = lowerbound:upperbound;

%make sure to stay inside model domain (this is w/o wrap around)

if lowerbound < model_idx(1)
	
	% keyboard

	% total_range = S_range(end) : upperbound;
	%
	% if lowerbound >= model_idx(1)
	%
	% 	total_range = [total_range (lowerbound : S_range(1))];
	%
	% else
	
	over = sum(total_range < model_idx(1));
	% lwrap = model_idx(1) : min(total_range) - 1;
	lwrap = model_idx(1) : upperbound;
	rwrap = (model_w - over) : model_w;
	% total_range = [(ledgefp) total_range (redgefp)];
	
	clear total_range
	total_range = [lwrap rwrap];
	
% end

end

if upperbound > model_w
	
% 	keyboard

	over = sum(total_range > model_w);
	
	lwrap = model_idx(1) : model_idx(over);
	rwrap = total_range(1 : (end-over));
	
	clear total_range
	
	total_range = [lwrap rwrap];
	
end
