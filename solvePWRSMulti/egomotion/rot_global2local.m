% global: R x +t encode as R(x-c)+c+t = Rx -Rc+c+t, so R|Rc-c+t
function Rloc = rot_global2local(R, centers)

Rloc = R;

for i = 1:size(R, 3)
  Rloc(1:3,4,i) = R(1:3,4,i) + squeeze(R(1:3,1:3,i)) * centers(:,i) - centers(:,i);
end
