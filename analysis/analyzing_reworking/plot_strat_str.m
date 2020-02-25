
function plot_strat_str(mod_str)

nits = mod_str.nits;
chan_pos = mod_str.chan_pos;
model_w = mod_str.model_w;
%%%%%%%%%%%%%%%%%%%%% PLOTTING
% FIGURE 1 - channel body rectangles and topography
%     %%% OPTION 2: Color code rectangles by timestep
figure(1)
cla

% hax=axes; 

topolines = false;

if topolines
	Topo = mod_str.Topo;
	n = 25;
    [top, ~] = size(Topo);
    tplot = floor(linspace(1,top,n));
    plot(Topo(tplot,:)','color',[0.75 0.75 0.75]);
end
hold on;

% cc=hsv(nits);
cc=[zeros(nits,3), 1.0*ones(nits,1)];
% keyboard
for ii = 1:nits
    if ~any(isnan(chan_pos(ii,:)))
       rectangle('Position',chan_pos(ii,:),'FaceColor',cc(ii,:)); % this plots the rectangular channel object on a figure     
    end
end

set(gca,'xlim',[0,model_w],'ylim',[0,(max(chan_pos(:,2))*1.2)]); %[0,100]); %[0,(max(max(Topo))+1)])  %,'ylim',[0,nits/(model_w/4)])

% line([mod_str.xind_low mod_str.xind_low],get(gca,'YLim'),'Color',[1 0 0])
% line([mod_str.xind_high mod_str.xind_high],get(gca,'YLim'),'Color',[1 0 0])
% line(get(gca,'XLim'),[mod_str.yind_low mod_str.yind_low],'Color',[1 0 0])
% line(get(gca,'XLim'),[mod_str.yind_high mod_str.yind_high],'Color',[1 0 0])

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
title('Synthetic stratigraphy'); xlabel('Model width'); ylabel('vertical aggradation (time)');

% hold off

end
