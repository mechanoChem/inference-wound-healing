clear; clc; close all;
set_plot_defaults_export
%% Import data
% basisFolder = '\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\basis\Physics_Based\basis_1D_3_3_rolling_win3\';
basisFolder = '../results/basis/Physics_Based/basis_1D_3_3_rolling_win1/';

basisFileParts = 'basis_step_%i_refine4.dat';
% basisFile = '\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\basis\Physics_Based\basis_1D_3_3_rolling_win3\basis_step_0_refine4.dat';
basisNamesWF = {'$-\nabla W \cdot \nabla C$','$-C*\nabla W \cdot \nabla C$', '$-C^{2}*\nabla W \cdot \nabla C$', ...
    '$C*\nabla W \cdot v$', '$C \nabla W \cdot v * C$', '$C* \nabla W \cdot v*C^{2}$', ...
    '$C*W$', '$C^{2}*W$', 'tracking'};

basisNames = {'$-\nabla^2 C$','$-C*\nabla^2 C$', '$-C^{2}*\nabla^2 C$', ...
    '$-v \cdot \nabla C$', '$-C* v \cdot \nabla C$', '$-C^2* v \cdot \nabla C$', ...
    '$C$', '$C^{2}$', 'tracking'};

figFolder = './figures_chapter/';
figFolder = './figures/';


origDataLoc = '../jtbDensityData.mat'
origData = load(origDataLoc)
xData = origData.data.pos; %in microns
nTimes = 5;
nDensity = 6;
densityVec = 20000:-2000:10000;
tVec = 0:12:48;


for tt = 1:nTimes
    for dd = 1:nDensity
        % basisFolder = strcat('\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\basis\Physics_Based\density',num2str(densityVec(dd)),'\basis_1D_3_3_rolling_win3\');
        basisFolder = strcat('../results/basis/Physics_Based/density',num2str(densityVec(dd)),'/basis_1D_3_3_rolling_win1/');

        basisFile = sprintf('%sbasis_step_%i_refine4.dat', basisFolder, tt-1);
        basis(:,:,tt,dd) = readmatrix(basisFile);
    end
end
denseCmap = flipud(lines(nDensity));
timeCmap = lines(nTimes);
    nBasis = size(basis, 2);




%% plot gamma matrix
for ii = 1:nDensity
    % gammaFolder = ['\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\VSI_gamma_matrix\Physics_Based_Time_Independent_1D\density' num2str(densityVec(ii)) '\' ];
    gammaFolder = ['../results/VSI_gamma_matrix/Physics_Based_Time_Independent_1D/density' num2str(densityVec(ii)) '/' ];

    % gamma_Group_3_3_rolling_win3_F200000_refine4
    gammaFile = 'gamma_history_Group_3_3_rolling_win1_F200000_refine4.dat';
    lossFile = 'loss_Group_3_3_rolling_win1_F200000_refine4.dat';
    gammaFinalFile = 'gamma_Group_3_3_rolling_win1_F200000_refine4.dat';
    gammaAll(:,:,ii) = readmatrix(strcat(gammaFolder, gammaFile));
    gammaFinal(:,:,ii) = readmatrix(strcat(gammaFolder, gammaFinalFile));
    loss(:,ii) = readmatrix(strcat(gammaFolder, lossFile), 'Delimiter', ' ');
end
% gammaFolderSimul = '\\wsl$\Ubuntu-20.04\home\pkinn\vsiTestWSL\Cell_Migration_density_checkbasisgen\results\VSI_gamma_matrix\Physics_Based_Time_Independent_1D\densitySimul\';
% gammaFolderSimul = '../old_results/VSI_gamma_matrix/Physics_Based_Time_Independent_1D/densitySimul/';
% 
% gammaFileSimul = 'gamma_history_Group_3_3_rolling_win3_F200000_refine4.dat';
% lossFileSimul = 'loss_Group_3_3_rolling_win3_F200000_refine4.dat';
% gammaFinalFileSimul = 'gamma_Group_3_3_rolling_win3_F200000_refine4.dat';
% gammaAllSimul = readmatrix(strcat(gammaFolderSimul, gammaFileSimul));
% gammaFinalSimul = readmatrix(strcat(gammaFolderSimul, gammaFinalFileSimul));
% lossSimul = readmatrix(strcat(gammaFolderSimul, lossFileSimul), 'Delimiter', ' ');

nIter = length(loss);
figure()
nTerms = sum(gammaAll(:,1:nIter) >0, 1);
timeColors = lines(nDensity);
for jj = 1:nDensity

    legText{jj} = sprintf('Initial density: %i cells', densityVec(jj))
    h = semilogy(nTerms, loss(:,jj), '-s','Color', timeColors(jj,:),'LineWidth',2,'MarkerSize',4,'MarkerFaceColor', timeColors(jj,:))
    hold on
end
% figure()
% nTerms = sum(gammaAll(:,1:nIter) >0, 1);
% timeColors = lines(nDensity);
% h = semilogy(nTerms, loss, '-s','Color', timeColors,'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor', timeColors)
% set(h, {'Color'}, num2cell(denseCmap,2))
xlabel('Number of Terms',FontSize=14);
ylabel('$\min_{\theta} ||\mathbf{R}(\mathbf{C^{h}}; \mathbf{\theta}) ||_2^2 $',FontSize=14, Interpreter='latex');
ylim([1e-8,1e-5])
grid on
legend(legText,FontSize=14);
figFolder = './figures/';
exportgraphics(gcf, [figFolder 'AllLossFunctions.pdf'], 'ContentType','vector');
% exportgraphics(gcf, [figFolder 'AllLossFunctions.jpg'], 'Resolution', 1500);
