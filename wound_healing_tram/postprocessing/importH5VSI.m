function [denseGrid, dim1Grid, dim2Grid] = importH5VSI(h5name, h5loc, nTimes)
%importh5VSI import h5 files with density from vsi code
%   Can be used to import experimental density or fwd solution h5 files.  

h5fn = [h5loc h5name];
geo = h5read(h5fn, '/Mesh/0/mesh/geometry');
topo = h5read(h5fn, '/Mesh/0/mesh/topology');
dGrid = geo(1,2)-geo(1,1);
geoIdx = geo./dGrid;

for ii= 1:nTimes
    timeFields{ii} = num2str(ii-1);
end

for ii = 1:nTimes
    dens(ii,:) = h5read(h5fn, sprintf('/VisualisationVector/%s', timeFields{ii}));
    for jj = 1:length(geoIdx)
        denseGrid(geoIdx(1,jj)+1, geoIdx(2,jj)+1, ii) = dens(ii,jj);
    end

end
dim1rng = (min(geo(1,:))):(geo(1,2)-geo(1,1)):(max(geo(1,:)));
dim2rng = (min(geo(2,:))):(geo(1,2)-geo(1,1)):(max(geo(2,:)));
[dim1Grid, dim2Grid] = meshgrid(dim1rng, dim2rng);


end