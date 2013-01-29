function raster(t_trunc,r_est,dt)
figure(4)
%subplot(2,1,1)
Y0=[ones(1,size(t_trunc,1));zeros(1,size(t_trunc,1))];
Y=Y0;

X = zeros([2*size(t_trunc,2) size(t_trunc,1)]);
for i=1:size(t_trunc,2)
    X((2*i-1):(2*i),:) = [t_trunc(:,i)';t_trunc(:,i)'];
end

for i=1:(size(t_trunc,2)-1)
    Y=[Y;Y0+i*ones(2,size(t_trunc,1))];
end

%axis([0 max(max(t_trunc))+1 -1 2])
for i=1:size(t_trunc,2)
    plot(X((2*i-1):(2*i),:),Y((2*i-1):(2*i),:),'k-')
    hold on
end
set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'YTick', []) % don't draw y-axis ticks
%set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
set(gca,'Color',get(gcf,'Color')) % match figure background
set(gca,'YColor',get(gcf,'Color')) % hide the y axis
axis('tight')
box off

% subplot(2,1,2)
% end_time = floor(max(max(t_trunc))/dt);
% plot((1:end_time)*dt,1000*r_est(1:end_time))
% set(gca,'Color',get(gcf,'Color')) % match figure background
% axis('tight')

box off
