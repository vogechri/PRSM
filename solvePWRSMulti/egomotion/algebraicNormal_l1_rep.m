%%% fit normals inot depth map -- IRLS
function [normal, l1Error, avS] = algebraicNormal_l1_rep (cam, p, q, repeats )

%e=  cam.Kl*cam.Tr;
e=  cam.Kr*cam.Tr;
A = cam.Kr*cam.Rr;

Ap = A*p;

b1 = (bsxfun(@times, e,p(1,:)));
b2 = (bsxfun(@times, e,p(2,:)));
b3 = (bsxfun(@times, e,p(3,:)));
B = -cat(2, b1(:), b2(:), b3(:));

b = q-Ap;b=b(:);


for i=1:repeats
  if exist('diagOut','var')
    B = bsxfun(@times, B, diagOut);
    b = (diagOut.*b);
  end

  v = (B'*B)\(B'*b(:));

  normal = v;
  l1Error = abs(B*v - b);
  diagOut = 1./max(0.1, l1Error); % micro gain
  
end

avS     = sum(l1Error) / size(p, 2);
return;