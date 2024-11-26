function visualizeFS_all()
clear; clc; %set_plot_defaults;

well_list = {'A01', 'A02', 'A03', 'A04', 'A05', 'A06', ...
'B01', 'B02', 'B03', 'B04', 'B05', 'B06', ...  
'C01', 'C02', 'C03', 'C04', 'C05', 'C06', ...
'D01', 'D02', 'D03', 'D04', 'D05', 'D06', ...  
};

for i=1:24
    plot_well(well_list{i})
end
%Code to visualize forward solutions:
%{
1. Locate data
2. Import and format H5/XDMF files
3. Visualize and compare experimental data and forward solution
4. Quantitatively compare expt and fwd solution (calculate loss externally)
%}

end
function plot_well(well)


% well = 'A02';
close all
stepVsi = 5;
stepAdj = stepVsi ;
figSave = './figuresChapter/FS/';
vsiFwdPath = '../results/forward_solution/VSI_2D_Time_Independent';
vsiFwdLoc = sprintf('%s/well%s/step%i/',vsiFwdPath, well, stepVsi);
vsiFwdDensH5 = [vsiFwdLoc 'density.h5'];
vsiFwdDensXD = [vsiFwdLoc 'density.xdmf'];

exptPath = '../results/PreProcess';
exptLoc = sprintf('%s/%s/', exptPath, well);
exptDensH5 = [exptLoc 'density_3_3_rolling_win5.h5'];
exptCellDensh5 = [exptLoc 'cell_density_3_3_rolling_win5.h5'];
exptXD = [exptLoc 'density_3_3_rolling_win5.xdmf'];

adjFwdPath = '../results/forward_solution/Adjoint_2D_Time_Independent';
adjFwdLoc = sprintf('%s/well%s/step%i/',adjFwdPath, well, stepAdj);
adjFwdDensH5 = [adjFwdLoc 'density.h5'];
adjFwdDensXD = [adjFwdLoc 'density.xdmf'];

nTimes = 72;
for ii= 1:nTimes
    timeFields{ii} = num2str(ii-1);
end
totalTime = 24;
timeVec = linspace(0,totalTime, nTimes);
% [mdlDensity, mdlMesh] = importFenicsModelDensity2D(vsiFwdDensH5, timeFields)

info = h5info(vsiFwdDensH5);
geo = h5read(vsiFwdDensH5, '/Mesh/0/mesh/geometry');
topo = h5read(vsiFwdDensH5, '/Mesh/0/mesh/topology');
dGrid = (geo(1,2)-geo(1,1));
geoIdx = geo./dGrid;
for ii = 1:nTimes
    dens(ii,:) = h5read(vsiFwdDensH5, sprintf('/VisualisationVector/%s', timeFields{ii}));
    densExpt(ii,:) = h5read(exptDensH5, sprintf('/VisualisationVector/%s', timeFields{ii}));
    densAdj(ii,:) = h5read(adjFwdDensH5, sprintf('/VisualisationVector/%s', timeFields{ii}));
    for jj = 1:length(geoIdx)
        denseGrid(geoIdx(1,jj)+1, geoIdx(2,jj)+1,ii) = dens(ii,jj);
        denseGridExpt(geoIdx(1,jj)+1, geoIdx(2,jj)+1,ii) = densExpt(ii,jj);
        denseGridAdj(geoIdx(1,jj)+1, geoIdx(2,jj)+1,ii) = densAdj(ii,jj);
    end
end
zLim = max([max(dens(:)) max(densExpt(:))]);
dim1rng = (min(geo(1,:))):(geo(1,2)-geo(1,1)):(max(geo(1,:)));
dim2rng = (min(geo(2,:))):(geo(1,2)-geo(1,1)):(max(geo(2,:)));
[dim1Grid, dim2Grid] = meshgrid(dim1rng, dim2rng);

figure('Position', [100 100 900 400])

set(gcf, 'color', 'w')
timeIdx = [1 36 72];
tl = tiledlayout(3,3)
colTitles = {'Experimental Data', 'VSI solution', 'Adjoint Solution'};
rowLabels = {'T = 0 hrs.', 'T = 12 hrs.', 'T = 24 hrs.'};
for ii = 1:3

    nexttile()
    surf(dim1Grid', dim2Grid', denseGridExpt(:,:,timeIdx(ii)));
    view(2)
    clim([0 zLim])
    xlim([0 max(dim1rng)])
    ylim([0 max(dim2rng)])


    if ii == 1
        title(colTitles{1})
    end
    ylabel(rowLabels{ii})
    

    nexttile()
    surf(dim1Grid', dim2Grid', denseGrid(:,:,timeIdx(ii)));
    view(2)
    clim([0 zLim])
        xlim([0 max(dim1rng)])
    ylim([0 max(dim2rng)])
    if ii == 1
        title(colTitles{2})
    end


    h = nexttile()
    surf(dim1Grid', dim2Grid', denseGridAdj(:,:,timeIdx(ii)));
    view(2)
    clim([0 zLim])
        xlim([0 max(dim1rng)])
    ylim([0 max(dim2rng)])
    if ii == 1
        title(colTitles{3})
    end
end
cbh = colorbar(h);
cbh.Layout.Tile = 'east'
ylabel(cbh, 'Cell Density (cells/um^2)')
% title(tl, sprintf('Well: %s Step: %i  Expt/VSI Fwd/Adj Fwd', well, stepAdj), 'FontSize', 24)
% ylabel(tl, sprintf('Top to bottom: T = %.1f, T = %.1f, T = %.1f', timeVec(timeIdx(1)), timeVec(timeIdx(2)), timeVec(timeIdx(3))), ...
%     'FontSize', 24)
exportgraphics(gcf, sprintf('%sFS_%s_step%i.jpg', figSave, well, stepVsi), 'Resolution',1500);
exportgraphics(gcf, sprintf('%sFS_%s_step%i_nonVector.pdf', figSave, well, stepVsi), 'ContentType','image');



%% Plot adjoint only
% figure()
% figure(Position=[100 100 1400 1000])
% set(gcf, 'color', 'w')
% timeIdx = [1 36 72];
% tl = tiledlayout(3,2)
% colTitles = {'Experimental Data', 'Adjoint Solution'};
% rowLabels = {'T = 0 hrs.', 'T = 12 hrs.', 'T = 24 hrs.'};
% for ii = 1:3
% 
%     nexttile()
%     surf(dim1Grid', dim2Grid', denseGridExpt(:,:,timeIdx(ii)));
%     view(2)
%     clim([0 zLim])
%     xlim([0 max(dim1rng)])
%     ylim([0 max(dim2rng)])
% 
% 
%     if ii == 1
%         title(colTitles{1})
%     end
%     ylabel(rowLabels{ii})
% 
% 
%     h = nexttile()
%     surf(dim1Grid', dim2Grid', denseGridAdj(:,:,timeIdx(ii)));
%     view(2)
%     clim([0 zLim])
%         xlim([0 max(dim1rng)])
%     ylim([0 max(dim2rng)])
%     if ii == 1
%         title(colTitles{2})
%     end
% end
% cbh = colorbar(h);
% cbh.Layout.Tile = 'east'
% ylabel(cbh, 'Cell Density (cells/um^2)')
% % title(tl, sprintf('Well: %s Step: %i  Expt/VSI Fwd/Adj Fwd', well, stepAdj), 'FontSize', 24)
% % ylabel(tl, sprintf('Top to bottom: T = %.1f, T = %.1f, T = %.1f', timeVec(timeIdx(1)), timeVec(timeIdx(2)), timeVec(timeIdx(3))), ...
% %     'FontSize', 24)
% 
% exportgraphics(gcf, sprintf('%sexptAdjCompare_%s_step%i.jpg', figSave, well, stepVsi), 'Resolution',1500);
% exportgraphics(gcf, sprintf('%sexptAdjCompare_%s_step%i_nonVector.pdf', figSave, well, stepVsi), 'ContentType','image');
end
