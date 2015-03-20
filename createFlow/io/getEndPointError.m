function [err, errImg] = getEndPointError(flow, flowGT)

[N, M, D] = size(flow);

test1 = (flow(:,:,1) - flowGT(:,:,1)).^2;
test2 = (flow(:,:,2) - flowGT(:,:,2)).^2;

errImg = sqrt( test1 + test2 );
if (size(flowGT,3) ~= 3)

  undef = ( (abs(flowGT(:,:,2)) > M ) | (abs(flowGT(:,:,1)) > N ));
  sumPix = N*M - sum(undef(:));
  errImg = errImg .* (1-undef);
  err = sum(sum( errImg )) / sumPix;

else
  errImg = errImg .* (flowGT(:,:,3)~=0);
  err = sum(sum( errImg(flowGT(:,:,3)~=0) )) / (sum(sum(flowGT(:,:,3)~=0)));
end
