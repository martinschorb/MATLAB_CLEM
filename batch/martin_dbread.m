function [val db] = martin_dbread(search,colnum,dbname,separator)

% % [val db] = martin_dbread(search,colnum,seperator,dbname)
% 
%  version MS 100712
% 
%   returns content of database file (raw text)
%   output is either an array of the desired values or a cell array of
%   strings.
% 
%   the input search either indicates the rownumber or if it is a string
%   the proper row is the first one containing that string.
%   
%   dbname can be either the filename of the database or the database
%   itself returned as second output of this function being run previously.
%   Use of the database speeds up algorithm by the factor of approx. 4.

% tic

if ~isempty (nargchk(2,4,nargin))
    error('MATLAB:martin_dbread:Nargin',nargchk(2,4,nargin));
end

if nargin < 3
   dbname = '/struct/briggs2/schorb/_Endo-Data/Endocytosis_DATA.csv' ;
end

if nargin < 4
    separator = ';' ;
end

if isstr(dbname)
    db = readtext(dbname,separator);
else
    db = dbname;
end

val={};

if isnumeric(search)
    for i = 1:length(search)
        for j = 1:length(colnum)
            val{i,j} = db{floor(search(i)),colnum(j)};
        end
    end
  
elseif ~isstr(search)
    error('MATLAB:martin_dbread','search argument must be integer or string');
else
    
    if ~(isequal(search(end),'1') & isequal(search(end),'"'))
        search = ['"',search,'"'];
    end
    
    
    IDX = find(cellfun(@(x) isstr(x), db));     
    
    for k = 1:length(IDX)
        s{k} = db{IDX(k)};
    end
    
    indeces = find(ismember(s,search));
    [I,J]=ind2sub(size(db),IDX(indeces));
    
    for i = 1:length(I)
        for j = 1:length(colnum)
            val{i,j} = db{I(i),colnum(j)};
        end
    end
    
end



if isempty([find(cellfun(@(x) isstr(x), val)) find(cellfun('isempty',val))])&~isempty(val)
   for l = 1:size(val,1)
       for m = 1:size(val,2)
           temp(l,m) = val{l,m};
       end
   end
    val = temp;
end

% toc



