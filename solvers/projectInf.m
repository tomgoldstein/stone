function rval = projectInf( x,y )


norm = sqrt(max(x.*x+y.*y,1));

x = x./norm;
y = y./norm;

rval = [x;y];
end

