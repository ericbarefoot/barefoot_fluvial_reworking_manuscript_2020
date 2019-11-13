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

function [total_range]=edgewrap(t,fp_lat,S_range,ch_w,model_idx,model_w)

        total_range = (min(S_range)-fp_lat*ch_w(t)):max(S_range)+fp_lat*ch_w(t);
            %make sure to stay inside model domain (this is w/o wrap around)
            if min(total_range)<model_idx(1)
                clear total_range
                total_range = S_range(end):S_range(end)+fp_lat*ch_w(t);
                if S_range(1)-fp_lat*ch_w >= model_idx(1)
                    total_range = [total_range S_range(1)-fp_lat*ch_w(t):S_range(1)];
                else 
                    over = abs(S_range(1)-fp_lat*ch_w(t));
                    total_range = [model_idx(1):min(total_range)-1 total_range model_w-over:model_w];
                end
            end
            if max(total_range)>model_w
                over = length(find(total_range>model_w));
                total_range = [model_idx(1):model_idx(max(over)) total_range(1:end-over)];
            end