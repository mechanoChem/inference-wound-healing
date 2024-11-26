% clear; clc; close all;
set_plot_defaults_export
basisNamesWF = {'$-\nabla W \cdot \nabla C$','$-C*\nabla W \cdot \nabla C$', '$-C^{2}*\nabla W \cdot \nabla C$', ...
    '$C*\nabla W \cdot v$', '$C \nabla W \cdot v * C$', '$C* \nabla W \cdot v*C^{2}$', ...
    '$C*W$', '$C^{2}*W$'};
% basisNames = {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
%     '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
%     '$C$', '$C^{2}$', 'tracking'};
basisNames =  {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
    '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
    '$C$', '$C^{2}$'};
figLoc = 'figures_datameeting\';

initConcList = {'10000', '12000', '14000', '16000', '18000', '20000'};
fig = figure()

timeLabels2 = {'  0h Experiment', '12h Experiment', '24h Experiment', '36h Experiment', '48h Experiment','  0h Jin et.al.', '12h Jin et.al.', '24h Jin et.al.', '36h Jin et.al.', '48h Jin et.al.','  0h PDE+Optim', '12h PDE+Optim', '24h PDE+Optim', '36h PDE+Optim', '48h PDE+Optim'};
% timeLabels2 = {'0 hours', '12 hours', '24 hours', '36 hours', '48 hours'};

for i=1:length(initConcList)
    
    initConc = initConcList{i};
    
    % h5Loc = strcat('\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\Adjoint_1D_Time_Independant\initCond', num2str(initConc), '\step6\density.h5');
    h5Loc = strcat('../results/forward_solution/Adjoint_1D_Time_Independant/initCond', num2str(initConc), '/step6/density.h5');
    
    [mdlDensity, mdlMesh] = importFenicsModelDensity1D(h5Loc, {'0', '1', '2', '3'});
    
    % exptLoc = strcat('\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\PreProcess\density', num2str(initConc), '\density_1D_3_3_rolling_win3_refine4.h5');
    exptLoc = strcat('../results/PreProcess/density', num2str(initConc), '/density_1D_3_3_rolling_win1_refine4.h5');
    
    [exptDensity, exptMesh] = importFenicsModelDensity1D(exptLoc, {'0', '1','2', '3','4'});


    nTimes = 5;
    %% Plot comparison of VSI and Adjoint
     % adjFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\Adjoint_1D_Time_Independant\initCond', initConc];
    % vsiFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\VSI_1D_Time_Independant\initCond', initConc];
    % jinFwdDir = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\forward_solution\UserDefined\initCond', initConc, '\step0\'];
    adjFwdDir = ['../results/forward_solution/Adjoint_1D_Time_Independant/initCond', num2str(initConc)];
    vsiFwdDir = ['../results/forward_solution/VSI_1D_Time_Independant/initCond', num2str(initConc)];
    jinFwdDir = ['../results/forward_solution/UserDefined/initCond', num2str(initConc), '/step0/'];    
        
    toImport = 6;
    
    % adjName = sprintf('%s\\step%i\\%s',adjFwdDir, toImport, 'density.h5');
    adjName = sprintf('%s/step%i/%s',adjFwdDir, toImport, 'density.h5');
    
    [adjDense, adjMesh] = importFenicsModelDensity1D(adjName, {'0', '1', '2','3'});
    adjDenseAll = adjDense;
    % vsiName = sprintf('%s\\step%i\\%s',vsiFwdDir, toImport, 'density.h5');
    vsiName = sprintf('%s/step%i/%s',vsiFwdDir, toImport, 'density.h5');
    
    [vsiDense, vsiMesh] = importFenicsModelDensity1D(vsiName, {'0', '1', '2','3'});
    vsiDenseAll = vsiDense;
    jinName = sprintf('%s%s', jinFwdDir, 'density.h5');
    [jinDense, jinMesh] = importFenicsModelDensity1D(jinName, {'0', '1', '2', '3'});
    jinDenseAll = jinDense
    
    % exptLoc = strcat('\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_rerun10_17_22\results\PreProcess\density', num2str(initConc), '\density_1D_3_3_rolling_win3_refine4.h5');
    exptLoc = strcat('../results/PreProcess/density', num2str(initConc), '/density_1D_3_3_rolling_win1_refine4.h5');
    
    [exptDensity, exptMesh] = importFenicsModelDensity1D(exptLoc, {'0', '1','2', '3','4'});
    
    
    %% Plot bar graph and model results
    %Bar graphs need to be formatted as [nBasis x 2];
    
    % figure()
    % tl = tiledlayout(1,3)
    % nexttile()
    subplot(4,2,i)
    
    timeCmap = lines(nTimes);
    
    timeColors = [[0 0 0]; lines(nTimes-1)];
    timeColors = timeCmap;
    
    plot(exptMesh, exptDensity(1,:),'LineStyle','none','Marker', 'square','Color', timeColors(1,:),'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor', timeColors(1,:))
    hold on
    % plot(exptMesh, exptDensity(1,:), 'Color', timeColors(1,:), 'LineStyle', '-','LineWidth',2.0);
    % grid on
    for ii = 1:nTimes-1
        hold on
        disp(ii)
        plot(exptMesh, exptDensity(ii+1,:),'LineStyle','none','Marker', 'square','Color', timeColors(ii+1,:),'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor', timeColors(ii+1,:))
        % h = plot(exptMesh, jinDenseAll(ii, :), 'Color', timeColors(ii+1,:), 'LineStyle', '-','LineWidth',2.0);
        % h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    % plot(exptMesh, exptDensity(1,:),'LineStyle','none','Marker', 'o','Color', timeColors(1,:),'LineWidth',1.2,'MarkerSize',2,'MarkerFaceColor', timeColors(1,:))
    hold on
    plot(exptMesh, exptDensity(1,:), 'Color', timeColors(1,:), 'LineStyle', '--','LineWidth',2.0);
    grid on
    for ii = 1:nTimes-1
        hold on
        disp(ii)
        % plot(exptMesh, exptDensity(ii+1,:),'LineStyle','none','Marker', 'o','Color', timeColors(ii+1,:),'LineWidth',1.2,'MarkerSize',2,'MarkerFaceColor', timeColors(ii+1,:))
        h = plot(exptMesh, jinDenseAll(ii, :), 'Color', timeColors(ii+1,:), 'LineStyle', '--','LineWidth',2.0);
    end
    hold on
    plot(exptMesh, exptDensity(1,:), 'Color', timeColors(1,:), 'LineStyle', '-','LineWidth',2.0);
    grid on
    for ii = 1:nTimes-1
        hold on
        disp(ii)
        h = plot(exptMesh, adjDenseAll(ii, :), 'Color', timeColors(ii+1,:), 'LineStyle', '-','LineWidth',2.0);
    end

    if i == 1
        leg = legend(timeLabels2, 'Orientation','Vertical','Location','EastOutside','NumColumns',3, FontSize=14);
    end 

    ylim([0 2.5e-3]);    
    title(strcat('Initial density: ',initConcList{i},' cells'),FontSize=14);
    % han=axes(fig,'visible','off'); 
    % han.XLabel.Visible='on';
    % han.YLabel.Visible='on';
    % xlabel(han,'x-Dimension (\mu m)',FontSize=14);
    % ylabel(han,'Concentration (cells/\mum^2)',FontSize=14);
    xlabel('x-Dimension (\mu m)',FontSize=14);
    ylabel({'Concentration','(cells/\mum^2)'},FontSize=14);
     



    
    
    % pause
    % saveas(gcf, [figLoc '*jCompareStep' num2str(basisIdx) '.tif'])
    % saveas(gcf, [figLoc 'VSIAdjCompareStep' num2str(basisIdx) '.jpg'])
end

% Create dummy subplot for legend
hLegend = subplot(4,2,7.5);
posLegend = get(hLegend,'Position');
axis(hLegend,'off');
set(leg,'Position',posLegend);

set(gcf, 'Position',[659    91   833   702]);
figFolder = './figures/';
exportgraphics(gcf, [figFolder 'jin_result.pdf'], 'ContentType','vector');
