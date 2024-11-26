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
groupNames = {'10 \muM Tram.', '5 \muM Tram.', '1 \muM Tram.', '500 nM Tram.', '100 nM Tram.', 'NT'};
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

%% Plot loss
nNonZeroBases = squeeze(sum(gammaHistAll~=0));
figure()
plot(nNonZeroBases(:,1:6), lossAll(:,1:6))
xlabel('Number of terms')
ylabel('Loss')
legend(groupNames(1:6))

exportgraphics(gcf, [figFolder 'lossPlot1.pdf'], 'ContentType','vector');
exportgraphics(gcf, [figFolder 'lossPlot1.jpg'], 'Resolution', 1500);