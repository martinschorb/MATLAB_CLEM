function out = martin_combin(n,k)

% martin_combin(input,k)
% 
%  returns all possible combinations of the elements of the input vector taking k at the same time.
% 
%   out: is a k by n_over_k matrix


if n == k
   out=1:n;
elseif n == k + 1
   tmp =1:n;
   out = tmp(ones(n,1),:);
   out(1:n+1:n*n) = [];
   out = reshape(out,n,n-1);
   
else   
      a=2^n;
      sel=zeros(a,n);
      for ii=1:n
        binary=[repmat(0,[a/(2^ii) 1]);repmat(1,[a/(2^ii) 1])];
        sel(:,ii)=repmat(binary,[2^(ii-1) 1]);
      end           
end

   idx = sel(sum(sel,2) == k,:);
   nrows = size(idx,1);
   [rows,xxx] = find(idx');
   
   vec=1:n;
   out = reshape(vec(rows),k,nrows).';

