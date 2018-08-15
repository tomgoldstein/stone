%%  
%  This function performs the STO transform by calling the c implementation
%  of STO_fast.  This method can be called on either a vector or a matrix.
%  When called on a matrix, it performs a separate transform on each
%  column.


function [ out count ] = STO( in )

persistent acount;
if isempty(acount) 
    acount = 0;
end

assert(ndims(in)<3, 'Input must be a vector or matrix');

%  Handle case that 'in' is a row vector
if min(size(in))==1
    out = in;
    % This line forces matlab NOT to use a lazy copy.  This is important
    % because otherwise a call to "STO_fast" will overwrite the input
    out(1)=out(1);
    STO_fast(out);
else

    out = in;
    for c=1:size(in,2);
        slice = out(:,c);
        STO_fast(slice);
        out(:,c) = slice;
    end

end
%acount = 0;
acount = acount + 1;
count = acount;
