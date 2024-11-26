clear; clc; close all;


%{
%% plot gamma matrix
for ii = 1:nDensity
    gammaFolder = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\VSI_gamma_matrix\Physics_Based_Time_Independent_1D\density' num2str(densityVec(ii)) '\' ];

    % gamma_Group_3_3_rolling_win3_F200000_refine4
    gammaFile = 'gamma_history_Group_3_3_rolling_win3_F200000_refine4.dat';
    lossFile = 'loss_Group_3_3_rolling_win3_F200000_refine4.dat';
    gammaFinalFile = 'gamma_Group_3_3_rolling_win3_F200000_refine4.dat';
    gammaAll(:,:,ii) = readmatrix(strcat(gammaFolder, gammaFile));
    gammaFinal(:,:,ii) = readmatrix(strcat(gammaFolder, gammaFinalFile));
    loss(:,ii) = readmatrix(strcat(gammaFolder, lossFile), 'Delimiter', ' ');
end
gammaFolderSimul = '\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\VSI_gamma_matrix\Physics_Based_Time_Independent_1D\densitySimul\';
gammaFileSimul = 'gamma_history_Group_3_3_rolling_win3_F200000_refine4.dat';
lossFileSimul = 'loss_Group_3_3_rolling_win3_F200000_refine4.dat';
gammaFinalFileSimul = 'gamma_Group_3_3_rolling_win3_F200000_refine4.dat';
gammaAllSimul = readmatrix(strcat(gammaFolderSimul, gammaFileSimul));
gammaFinalSimul = readmatrix(strcat(gammaFolderSimul, gammaFinalFileSimul));
lossSimul = readmatrix(strcat(gammaFolderSimul, lossFileSimul), 'Delimiter', ' ');

nIter = length(loss);
%}


basisNames = {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
    '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
    '$C$', '$C^{2}$', 'tracking'};

figFolder = './figuresChapter/';

gammaPath = '../results/VSI_gamma_matrix/';
wellGroups = {'A01_to_D01', 'A02_to_D02', 'A03_to_D03', 'A04_to_D04', 'A05_to_D05', 'A06_to_D06'};
groupNames =  {'10 \muM Tram.', '5 \muM Tram.', '1 \muM Tram.', '500 nM Tram.', '100 nM Tram.', 'NT'};
nGroups = length(wellGroups);
gammaFile = 'gamma_Group_3_3_rolling_win5_F200000.dat';
gammaHistoryFile = 'gamma_history_Group_3_3_rolling_win5_F200000.dat';
lossFile = 'loss_Group_3_3_rolling_win5_F200000.dat';
for ii = 1:nGroups
    gammaFolder = sprintf('group%i_%s/Physics_Based_Time_Independent/', ii, wellGroups{ii});
    gammaFinalAll(:,:,ii) = readmatrix([gammaPath gammaFolder gammaFile]);
    gammaHistAll(:,:,ii) = readmatrix([gammaPath gammaFolder gammaHistoryFile]);
    lossAll(:,ii) = importLossFile([gammaPath gammaFolder lossFile])
end
gammaPathAdj = '../results/Adjoint_gamma_matrix/';
adjSteps = 5:7;
for ii = 1:nGroups
    for jj = 1:length(adjSteps)
        step = adjSteps(jj);
        adjGroup = sprintf('group%d_%s_step%d/', ii, wellGroups{ii},step);
        adjFile = sprintf('%s%s/Physics_Based_Time_Independent/%s', gammaPathAdj, adjGroup, gammaFile);
        adjMatrix = readmatrix(adjFile);
        adjAll(:,jj,ii) = adjMatrix(1,:);
    end
end


nTimes = 72;
tVec = linspace(0, 24, nTimes);
nBasis = 8; 
bVec = 1:nBasis;

%% Plot VSI progress
figure('Position', [100 100 900 400])
% colors = {'r','b', 'r', 'b', ' '};
% linetype = {'-', '-', '--', '--'};
tl = tiledlayout(3,3)
for ii = 1:nBasis
    nexttile
    for jj = 1:nGroups
        plot(bVec, gammaHistAll(ii,:,jj) )
%         plot(bVec, gammaHistAll(ii,:,jj), "Color",colors{jj}, "LineStyle",linetype{jj})
        title(basisNames{ii}, 'Interpreter','latex')
        hold on
    end

end
xlabel(tl, 'VSI iteration', 'FontSize', 24)
ylabel(tl, 'Basis magnitude', 'FontSize', 24)
legend(groupNames)

exportgraphics(gcf, [figFolder 'plotResults1.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults1.jpg'], 'Resolution', 1500);

%% Compare VSI and Adjoint
figure('Position', [100 100 900 400])
tl = tiledlayout(2,3)
barLbl = basisNames([1 7 8])
barLbl = categorical(barLbl);
barLbl = reordercats(barLbl, basisNames([1 7 8]));
for ii = 1:nGroups
    nexttile
    bar(barLbl, [gammaHistAll([1 7 8], 6,ii) adjAll([1 7 8],1, ii)])
    title(groupNames{ii})

    set(gca, 'TickLabelInterpreter','latex');
%     if ii<6
%         nexttile
%     end

end
ylabel(tl,'Basis Magnitude', 'FontSize', 24);

exportgraphics(gcf, [figFolder 'plotResults2.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults2.jpg'], 'Resolution', 1500);


%%
figure('Position', [100 100 900 400])

tl = tiledlayout(2,3)
for ii = 1:nGroups
    nexttile
    basisRatio = adjAll([1 7 8], 1, ii)./gammaHistAll([1 7 8], 6, ii);
    bar(barLbl, basisRatio);
    title(groupNames{ii})

    set(gca, 'TickLabelInterpreter','latex');
%     if ii<6
%         nexttile
%     end

end
ylabel(tl, 'Adj/VSI ratio', 'FontSize', 24);
   

exportgraphics(gcf, [figFolder 'plotResults3.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults3.jpg'], 'Resolution', 1500);

%%
figure('Position', [100 100 900 400])
tl = tiledlayout(2,3)
barLbl = basisNames([7])
barLbl = categorical(barLbl);
barLbl = reordercats(barLbl, basisNames([7]));
for ii = 1:nGroups
    nexttile
    bar(barLbl, [gammaHistAll([7], 6,ii) adjAll([7],1, ii)])
    title(groupNames{ii})
    set(gca, 'TickLabelInterpreter','latex');

end
ylabel(tl, 'Basis Magnitude', 'FontSize', 24);
exportgraphics(gcf, [figFolder 'plotResults4.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults4.jpg'], 'Resolution', 1500);

%% Compare final bases b/w conditions
figure('Position', [100 100 900 400])
tl = tiledlayout(1,3);
barLbl = groupNames;
barLbl = categorical(barLbl);
barLbl = reordercats(barLbl, groupNames);
basesToPlot = [1 7 8];
for ii = 1:length(basesToPlot)
    nexttile()
    bar(barLbl, squeeze(adjAll(basesToPlot(ii), 1, :)))
    title(basisNames{ii}, 'Interpreter','latex')
end
exportgraphics(gcf, [figFolder 'plotResults5.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults5.jpg'], 'Resolution', 1500);

%%
figure('Position', [100 100 900 400])
cellConcVec = linspace(0, 0.003, 100);
figure()
for ii = 1:nGroups
    rxnVec = cellConcVec.*adjAll(7, 1, ii) + (cellConcVec.^2).*adjAll(8, 1, ii);
    plot(cellConcVec, rxnVec)
    hold on
end
legend(groupNames)
xlabel('Concentration (cells/um^2)')
ylabel("Reaction Term")
exportgraphics(gcf, [figFolder 'plotResults6.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'plotResults6.jpg'], 'Resolution', 1500);
