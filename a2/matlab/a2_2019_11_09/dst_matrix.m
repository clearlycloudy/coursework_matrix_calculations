function c = dst_matrix(n)
%DSTMTX Discrete sine transform matrix.
%   D = DSTMTX(N) returns the N-by-N DST transform matrix.  D*A
%   is the DST of the columns of A and D'*A is the inverse DST of
%   the columns of A (when A is N-by-N).
%
%   If A is square, the two-dimensional DST of A can be computed
%   as D*A*D'. This computation is sometimes faster than using
%   DST2, especially if you are computing large number of small
%   DST's, because D needs to be determined only once.
%
%   Class Support
%   -------------
%   N is an integer scalar of class double. D is returned 
%   as a matrix of class double.
%   
%   Example
%   -------
%       A = im2double(imread('rice.png'));
%       D = dstmtx(size(A,1));
%       dst = D*A*D';
%       figure, imshow(dst)
%
%   See also DST2.
%   Author: Said BOUREZG  
%   Electronics Engineer   Option:Communication .
%   Date: 07.02.2015
%   Adress:                        Said BOUREZG
%                               Elbassatine street
%                                 28038 Tarmount
%                               M'sila --- Algeria 
%   Email:  said.bourezg@yahoo.fr
%   Mobile: +213 796 018049 
%   If you can improve this code furtherly, please let me know. Thanks
%   Filename: dst2.m (Matlab)
%   Copyright 2015 Said BOUREZG.
%   I/O Spec
%      N - input must be double
%      D - output DST transform matrix is double
iptchecknargin(1,1,nargin,mfilename);
iptcheckinput(n,{'double'},{'integer' 'scalar'},mfilename,'n',1);
[cc,rr] = meshgrid(0:n-1);
c = sqrt(2 / n) * sin(pi * (2*cc + 1) .* rr / (2 * n));
c(1,:) = c(1,:) / sqrt(2);