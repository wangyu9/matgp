function [] = bold_plot(u)

h = figure('Position',[100,100,1600,1200],'PaperPositionMode', 'auto');
%h = figure('PaperPositionMode', 'auto');
set(gcf,'papersize',[16 12])

set(gca,'fontsize',20,'linewidth',2)

hold on;

if size(u,1)==1
   u = u'; 
end

c = size(u,2);
for i=1:c
    plot(0:size(u,1)-1,u(:,i),'LineWidth',3);
end

%ylim([0,35])

%h_legend = legend('','');
%set(h_legend,'FontSize',25,'Location','best');

box on;

% make the border tight to save space
% https://www.mathworks.com/matlabcentral/newsreader/view_thread/156326
if 0
    % tight
    set(gca,'LooseInset',get(gca,'TightInset'))
end
% https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
