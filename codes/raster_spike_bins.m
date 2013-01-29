function raster_spike_bins(n,bin_size)
figure(5)
%subplot(2,1,1)
n(n>1)=1;
X = [1:size(n,1);1:size(n,1)]*bin_size;
for i=1:size(n,2)
    Y = [i*n(:,i)';(i-1)*n(:,i)'];
    plot(X,Y,'k-')
    hold on
end
hold off

set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'YTick', []) % don't draw y-axis ticks
set(gca,'Color',get(gcf,'Color')) % match figure background
set(gca,'YColor',get(gcf,'Color')) % hide the y axis
axis('tight')
box off

% subplot(2,1,2)
% end_time = max(max(t_trunc))/dt;
% plot((1:end_time)*dt,1000*r_est(1:end_time))
% set(gca,'Color',get(gcf,'Color')) % match figure background
% axis('tight')

%box off