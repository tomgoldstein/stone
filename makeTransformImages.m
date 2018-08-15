% Create a visualization of a STOne transform matrix

close all;
N = 16;
s = STO(eye(N));
s = flipud(s);
imagesc(s);
axis off

colormap gray;

s = round(sqrt(N)*s);

for r = 1:N
    fprintf(' %d ',s(r,1));
    for c=2:N
         fprintf('& %d ',s(r,c));
    end
     fprintf(' \\\\ \n ');
end
        