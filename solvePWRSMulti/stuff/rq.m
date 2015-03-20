%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% computes the rq decomposition of the matrix P. Returns the matrices R and
%% Q, both \in \mathbf{R^{nxm}} with P \in \mathbf{R^[nxm}}. Here 
%% R is an upper triangle matrix and Q an unitary matrix so that P = R*Q.
%%%%%%%
%% P is a matrix P \in \mathbf{R^[nxm}}
%% R is a matrix R \in \mathbf{R^[nxm}}, R is an upper triangle matrix
%% Q is a matrix Q \in \mathbf{R^[nxm}}, Q is an unitary matrix
function [ R,Q ] =  rq ( P )
% some algebraic witchcraft to get RQ-decomposition
% from matlab built-in QR-decomposition

[Qh Rh] = qr(flipud(P).');
Rh = flipud(Rh.');
R = Rh(:,end:-1:1);
Qh = Qh.';
Q = Qh(end:-1:1,:);  