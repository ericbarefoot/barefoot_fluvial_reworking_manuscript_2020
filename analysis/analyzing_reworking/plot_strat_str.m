
function plot_strat_str(mod_str)

nits = mod_str.nits;
chan_pos = mod_str.chan_pos;
Topo = mod_str.Topo;
model_w = mod_str.model_w;
%%%%%%%%%%%%%%%%%%%%% PLOTTING
% FIGURE 1 - channel body rectangles and topography
%     %%% OPTION 2: Color code rectangles by timestep
figure(1)
% plot(Topo','color','k');
hold on;

% cc=hsv(nits);
cc=zeros(nits,3);
% keyboard
for ii = 1:nits
    if ~any(isnan(chan_pos(ii,:)))
       rectangle('Position',chan_pos(ii,:),'FaceColor',cc(ii,:)); % this plots the rectangular channel object on a figure     
    end
end

% % %            rectangle('Position',chan_pos(t,:),'FaceColor','k')
%             pause(.1)
% %
% %     %%% Plot topographic surfaces
% % % %
%            pause(1)
%
%     %%% Plot channel body centroids
% %         plot(ch_center(t,1),ch_center(t,2),'*')
%
%%Set axis properties
%axis equal
set(gca,'xlim',[0,model_w],'ylim',[0,(max(max(Topo))+1)]); %[0,100]); %[0,(max(max(Topo))+1)])  %,'ylim',[0,nits/(model_w/4)])
title('Synthetic stratigraphy'); xlabel('Model width'); ylabel('vertical aggradation (time)');

end
