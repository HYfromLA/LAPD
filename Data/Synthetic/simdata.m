function [Data, LabelsGT] = simdata(DATAopts)

%% 1d Examples

if strcmp(DATAopts.shape, 'two lines')
    [Data, LabelsGT] = shape_2lines(DATAopts);    
end

if strcmp(DATAopts.shape, 'dollar sign')
    [Data, LabelsGT] = shape_dollarsign(DATAopts);    
end

if strcmp(DATAopts.shape, 'olympic rings')
    [Data, LabelsGT] = shape_olympicrings(DATAopts);    
end

if strcmp(DATAopts.shape, 'three curves')
    [Data, LabelsGT] = shape_3curves(DATAopts);    
end

if strcmp(DATAopts.shape, 'rose circles')
    [Data, LabelsGT] = shape_rosecircles(DATAopts);    
end

%% 2d Examples

if strcmp(DATAopts.shape, 'two planes')
    [Data, LabelsGT] = shape_2planes(DATAopts);    
end

if strcmp(DATAopts.shape, 'two 2spheres')
    [Data, LabelsGT] = shape_2twospheres(DATAopts);    
end

if strcmp(DATAopts.shape, 'two triangles') 
    [Data, LabelsGT] = shape_2triangles(DATAopts);    
end

if strcmp(DATAopts.shape, 'three planes') 
    [Data, LabelsGT] = shape_3planes(DATAopts);    
end

if strcmp(DATAopts.shape, 'cone plane') 
    [Data, LabelsGT] = shape_coneplane(DATAopts);    
end

if strcmp(DATAopts.shape, 'swiss roll')  
    [Data, LabelsGT] = shape_swissroll(DATAopts);    
end

%% 3d Examples

if strcmp(DATAopts.shape, 'two cuboids')
    [Data, LabelsGT] = shape_2cuboids(DATAopts);    
end

if strcmp(DATAopts.shape, 'two 3spheres')
    [Data, LabelsGT] = shape_2threespheres(DATAopts);    
end

%% 4d examples
if strcmp(DATAopts.shape, 'two tesseracts')
    [Data, LabelsGT] = shape_2tesseracts(DATAopts);    
end

if strcmp(DATAopts.shape, 'two 4spheres')
    [Data, LabelsGT] = shape_2fourspheres(DATAopts);    
end

%% 5d examples
if strcmp(DATAopts.shape, 'two penteracts')
    [Data, LabelsGT] = shape_2penteracts(DATAopts);    
end

if strcmp(DATAopts.shape, 'two 5spheres')
    [Data, LabelsGT] = shape_2fivespheres(DATAopts);    
end




end
