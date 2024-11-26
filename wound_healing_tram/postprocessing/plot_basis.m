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
groupNames =  {'10 \muM Tram.', '5 \muM Tram.', '1 \muM Tram.', '500 nM Tram.', '100 nM Tram.', 'NT'}
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


%% Compare final bases b/w conditions
% figure()
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

exportgraphics(gcf, [figFolder 'basisPlot1.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'basisPlot1.jpg'], 'Resolution', 1500);

%% Plot both bases in one figure
% figure()
figure('Position', [100 100 900 400])
tl = tiledlayout(1,2);
cellConcVec = linspace(0, 0.003, 100);
colors = {'g','k','m','c','r', 'b'};
nexttile
nGroups_start = 5;
nGroups_end = 6;
for ii = nGroups_start :nGroups_end
    difVec = adjAll(1,1,ii) + cellConcVec.*adjAll(2, 1, ii) + (cellConcVec.^2).*adjAll(3, 1, ii);
    plot(cellConcVec, difVec, colors{ii})
    hold on
end
ylim([0 60])
ylabel('Diffusivity (\mum^2/hr)')

nexttile
for ii = nGroups_start :nGroups_end
    rxnVec = cellConcVec.*adjAll(7, 1, ii) + (cellConcVec.^2).*adjAll(8, 1, ii);
    plot(cellConcVec, rxnVec, colors{ii})
    hold on
end
legend(groupNames(nGroups_start :nGroups_end))
xlabel(tl, 'Concentration (cells/um^2)')
ylabel("Reaction Term (1/hr)")

exportgraphics(gcf, [figFolder 'basisPlot2.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'basisPlot2.jpg'], 'Resolution', 1500);


%% Plot loss and both bases in one figure
figure('Position', [100 100 900 400])
tl = tiledlayout(1,3);
cellConcVec = linspace(0, 0.003, 100);
colors = {'g','k','m','c','r', 'b'};
nexttile

nNonZeroBases = squeeze(sum(gammaHistAll~=0));
for ii = nGroups_start :nGroups_end
    plot(nNonZeroBases(:,ii), lossAll(:,ii), colors{ii})
    hold on
end
xlabel('Number of terms')
ylabel('Loss')
xlim([1 8])
legend(groupNames(nGroups_start :nGroups_end))

nexttile

for ii = nGroups_start :nGroups_end
    difVec = adjAll(1,1,ii) + cellConcVec.*adjAll(2, 1, ii) + (cellConcVec.^2).*adjAll(3, 1, ii);
    plot(cellConcVec, difVec, colors{ii})
    hold on
end
ylim([0 60])
ylabel('Diffusivity (\mum^2/hr)')

nexttile
for ii = nGroups_start :nGroups_end
    rxnVec = cellConcVec.*adjAll(7, 1, ii) + (cellConcVec.^2).*adjAll(8, 1, ii);
    plot(cellConcVec, rxnVec, colors{ii})
    hold on
end
xlabel(tl, 'Concentration (cells/um^2)')
ylabel("Reaction Term (1/hr)")

exportgraphics(gcf, [figFolder 'basisPlot3.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'basisPlot3.jpg'], 'Resolution', 1500);
