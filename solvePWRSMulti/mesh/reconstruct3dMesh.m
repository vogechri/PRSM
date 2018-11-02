function Mesh =  reconstruct3dMesh(cam, N_res, Seg)
% RECONSTRUCT3DMESH Compute 2D/3D triangulation of the PRSM segmentation
%   [Mesh,ref,cam] =  reconstruc3dMesh( ref, cam, N_res, Seg ) returns a
%   struct, Mesh, that contains a 2D/3D triangulation of the Piecewise
%   Rigid Planar representation of the scene in the image corresponding to
%   Seg and N_res. cam is a PRSM camera struct, containing the camera
%   intrinsics. Seg is a PRSM segmentation struct containing a label image,
%   a int image where every pixel is a plane assignment label. N_res is a
%   list of normals specifying candidate planes in the scene flow estimate.
%
%   Example (see more in meshHelper.m)
%   -------
%     %  call in pwrsf_v4.m
%     Mesh =  reconstruc3dMesh(cam(1), N_lin, SegNew);
%
%     % For display, filter out far away planes
%     validVertices = all(abs(Mesh.vertices)<40,2);
%     validFaces = all(reshape(validVertices,4,size(validVertices,1)/4)',2);
%     validVertices = Mesh.Findices(validFaces,:)';
%     Mesh.vertices = Mesh.vertices(validVertices(:),:);
%     Mesh.vertices_texture = Mesh.vertices_texture(validVertices(:),:);
%     Findices = reshape((1:size(Mesh.vertices,1)),4,size(Mesh.vertices,1)/4)';
%
%     % display
%     figure; patch('Faces',Findices,'Vertices',Mesh.vertices);
%     title('Mesh without color');
%
%
% See also pwrsf_v4, reconstruc2dFlowHom, convert2DwPlaneto3D,
% meshHelper.m, write_wobj
V3D = {}; V2D = {};
klinv = cam.Kl^-1;
% Triangulate in image every segment
% Seg.Img contains zero indexed segment Ids. To access correct index in
% normal vector, increment
[segBoundarySubscripts, segmentIDs] = findTrianglesSubPixel(Seg.Img+1);
% Convert image triangles to 3D mesh triangles
for i = 1:numel(segmentIDs)
    segmentID = segmentIDs(i);
    triangleCollection = segBoundarySubscripts{segmentID};
    % loop over triangles with ID and compute their 3D position
    for j = 1:numel(triangleCollection)
        poly2 =  triangleCollection{j};
        poly3 = convert2DwPlaneto3D(N_res(:,segmentID),poly2(:,[2,1]),klinv);
        V3D{end+1,1} = poly3; V2D{end+1,1} = poly2;
    end
end
% Convert from linkedlist-like-structure to array-like-structure
V2D = cell2mat(V2D);  V3D = cell2mat(V3D);

% scale to [0,1] and swap x and y
V2D = (V2D(:,[2,1])) ./ [ size(Seg.Img,2), size(Seg.Img,1)];
V2D(:,2) = 1-V2D(:,2);

Mesh.Findices = reshape((1:size(V2D,1)),4,size(V2D,1)/4)';
Mesh.vertices = V3D;
Mesh.vertices_texture = V2D;
