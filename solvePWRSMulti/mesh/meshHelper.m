function [ cleanedMesh ] = meshHelper( meshPath, varargin )
%MESHHELPER Performs common operations on .mat files containing meshes
%   Given a path to a .mat file and some options, it loads a mesh and
%   texture from that file and performs operations on it. The .mat file
%   must contain variables Mesh and imgLRGB. Mesh is a struct specifying
%   mesh faces with 3D world vertices and 2D mesh vertices that correspond
%   to the image in imgLRGB. Parameters then control which operations this
%   function performs. This function can use the Wavefront OBJ toolbox to
%   save to a common mesh file format (wavefront wobj) which allows it to
%   be opened in other applications, such as meshlab
%
%   Toolbox:
%   https://www.mathworks.com/matlabcentral/fileexchange/27982-wavefront-obj-toolbox
%   Meshlab viewing tool:
%   http://www.meshlab.net/
%
%   Parameters (specifed after path as 'parametername', parameter value)
%
%   plotWPatch        Whether to create a figure and draw a 3D, untextured,
%                     representation of the mesh. 
%                     Default: 1
%   saveWWOBJToolbox  Whether to save the mesh as a wavefront wobj. This
%                     create 3 files alognside the input .mat file: .obj
%                     file specifying triangles in space and on a texture
%                     image, the texture image as a .bmp, and a .mtl
%                     parameter file.
%                     Default: 1
%   maxCoord          Filter value for excluding triangles with points that
%                     are too far away. Helpful for viewing mesh.
%                     Default: 40
%
% See also run_pwrs_red.m, write_wobj, reconstruct3dMesh.m

p = inputParser;
valBool = @(x)validateAttributes(x,{'numeric'},{'finite','scalar'});
addParameter(p, 'plotWPatch',       1, valBool);
addParameter(p, 'saveWWOBJToolbox', 1, valBool);
addParameter(p, 'maxCoord',        40, ... Whatever units in .mat file
    @(x)validateAttributes(x,{'numeric'},{'finite','positive','scalar'}));
parse(p,varargin{:});    

[inputDir,inputName,~] = fileparts(meshPath);
load(meshPath);

% Clean mesh
% Currently only removes triangles with at least one far coordinate could
% consider removing small triangles, hidden triangles, etc.
validVertices = all(abs(Mesh.vertices) < p.Results.maxCoord, 2);
validFaces = all(reshape(validVertices, 4, size(validVertices, 1)/4)', 2);
validVertices = Mesh.Findices(validFaces, :)';
cleanedMesh.vertices = Mesh.vertices(validVertices(:), :);
cleanedMesh.vertices_texture = Mesh.vertices_texture(validVertices(:), :);
cleanedMesh.Findices = reshape((1:size(cleanedMesh.vertices,1)),...
    4, size(cleanedMesh.vertices, 1)/4)';

% Export Mesh
if p.Results.saveWWOBJToolbox
    if exist('write_wobj','file')
        outputPath = fullfile(inputDir,inputName);    % With extension stripped
        cleanedMesh = saveWOBJToolbox(outputPath,cleanedMesh,imgLRGB);
    else
        warning(['Failed to save, write_wobj not found, get from:'...
            'https://www.mathworks.com/matlabcentral/fileexchange/27982-wavefront-obj-toolbox']);
    end
end

% Plot Mesh
if p.Results.plotWPatch
    figure; patch('Faces',cleanedMesh.Findices,'Vertices',cleanedMesh.vertices);
    title(inputName);
end
end

function cleanedMesh = saveWOBJToolbox(outputPath,cleanedMesh,img)
% Saves mesh using example paramters from WOBJ_toolbox

[outputDir,outputName,~] = fileparts(outputPath);
cleanedMesh.objects(1).type ='g';
cleanedMesh.objects(2).type ='usemtl';
cleanedMesh.objects(3).type ='f';
material(1).type='newmtl';
material(2).type='Ka';     material(2).data=[0.8 0.8 0.8];
material(3).type='Kd';     material(3).data=[0.8 0.8 0.8];
material(4).type='Ks';     material(4).data=[1 1 1];
material(5).type='illum';  material(5).data=1;
material(6).type='Ns';     material(6).data=27;
material(7).type='map_Kd';

% Not very space efficient, could use unique to clean up significantly
cleanedMesh.objects(3).data.vertices = cleanedMesh.Findices;
cleanedMesh.objects(3).data.texture = cleanedMesh.Findices;

cleanedMesh.objects(1).data = outputName;
cleanedMesh.objects(2).data = [outputName,'_tex'];
material(7).data = [outputName,'.bmp'];
material(1).data = [outputName,'_tex'];
cleanedMesh.material = material;

curDir = cd(outputDir);
% needs to save to current directory
write_wobj(cleanedMesh,[outputName,'.obj']);
imwrite(img,[outputName,'.bmp']);
cd(curDir);
end
