function [density,mesh] = importFenicsExptDensity1D(filePath, timeFields)
%importFenicsDensity1D(filepath): Get density and mesh characteristics from
%fenics file

nTimes = length(timeFields);
density = [];
for ii = 1:nTimes
    dc = h5read(filePath,strcat('/density_', timeFields{ii}));
    density = [density; dc];

end

meshAll = h5read(filePath,strcat('/Mesh/0/mesh/geometry'));
mesh = meshAll(1,:);
end