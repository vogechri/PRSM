function [segBoundarySubscripts,segmentIDs] = ...
    findTrianglesSubPixel(inputSegmentationImage)
% for a label/segmentation image, finds all edges of segments, then
% triangulates the edges. Subpixel boundary computation inspired by
% <http://scikit-image.org/docs/dev/api/skimage.segmentation.html#skimage.segmentation.find_boundaries scikit-image>
% returns boundaries as cell array associated with segmentID in
% correspponding index in segmentIDs
%
% See also reconstruct3dMesh


% Add interstitial rows/cols that will contain segment boundaries
bigSeg = -ones(2*size(inputSegmentationImage)+1);
[xinSeg,yinSeg] = meshgrid(2:2:(size(inputSegmentationImage,2)*2), 2:2:(size(inputSegmentationImage,1)*2));
indicesInSeg = sub2ind(size(bigSeg),yinSeg(:),xinSeg(:));
bigSeg(indicesInSeg) = inputSegmentationImage(:);
% setup
segmentIDs = unique(inputSegmentationImage(:));
dilateStrel = strel('arbitrary',ones(3));
segBoundarySubscripts = cell(max(segmentIDs),1);
% matlab only supports binary contour extraction for segmentations, so loop
for segmentID = segmentIDs'
    tmpMask = bigSeg == segmentID;
    % grow mask into interstitial rows, allowing for subpixel countour
    tmpMask = imdilate(tmpMask, dilateStrel);
    tmpMaskHoles = imfill(tmpMask,8,'holes') & ~tmpMask;
    % grow holes/shrink mask into intersititial rows for holes
    tmpMaskHoles = imdilate(tmpMaskHoles, dilateStrel);
    tmpMask = tmpMask & ~tmpMaskHoles;
    [tmpBoundaryCell,~,~,familyTree] = bwboundaries(tmpMask,8);
    % Remove redundant points and convert to original image coordinates
    tmpBoundaryCell = simplify_polygon(tmpBoundaryCell);
    boundaryTriangles = triangulationWrapper(tmpBoundaryCell, familyTree);
    segBoundarySubscripts{segmentID} = boundaryTriangles;
end
end

function boundaryTriangles = triangulationWrapper(boundaryShapes, familyTree)
% combines edges and holes for one plane, projected to some number of
% segments for constrained delauney, looped above to cover all planes
boundaryTriangles = {};
for i = 1:numel(boundaryShapes)
    if ~any(familyTree(i,:))
        % loop over parents, then add their children (holes) as additional
        % constraints
        coords = boundaryShapes{i}(1:(end-1),:);
        last = size(coords,1);
        consts = [(1:(last-1))',(2:last)';last,1];
        for j = 1:numel(boundaryShapes)
            if familyTree(j,i)
                % hole found, add as additional constraint on parent
                start = size(coords,1)+1; 
                last = start + size(boundaryShapes{j},1)-2;
                coords = [coords; boundaryShapes{j}(1:(end-1),:)];
                consts = [consts; (start:(last-1))',((start+1):last)';last,start];            
            end
        end
        % remove duplicates and remap constraints to match
        [coords,~,iCoords] = unique(coords,'rows');
        consts2 = iCoords(consts);
        DT = delaunayTriangulation(coords(:,1),coords(:,2),consts2);
        DT = DT(isInterior(DT),:);   % dims: [triangle, vertex, singleton]
        DT = DT(:,[1,2,3,1]); % Complete the triangle
        DT = permute(DT, [2, 3, 1]); % dims: [vertex, singleton, triangle]
        shapeMat = [coords(sub2ind(size(coords),DT,ones(size(DT)))),...
            coords(sub2ind(size(coords),DT,2*ones(size(DT))))];   % dims: [vertex, coord(x,y), triangle]
        % permute changes from [1,1,n] to [n,(1),(1)] size from num2cell
        shapeCell = permute(num2cell(shapeMat,[1,2]),[3,1,2]);
        boundaryTriangles = cat(1,boundaryTriangles,shapeCell);
    end
end
end

function Pc = simplify_polygon(Pc)
% taken from
%http://blogs.mathworks.com/steve/2012/08/28/wrapping-up-the-analysis-of-cody-solutions/
for i = 1:numel(Pc)
    P = Pc{i}/2 - 0.5; % convert coordinate systems
try
    diag(sum(abs(diff(P)),2)) \ diff(P);
    P(any(ans - circshift(ans,1),2),:);
    P = vertcat(ans, ans(1,:));
end
    Pc{i} = P;
end

end