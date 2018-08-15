%%  Write an image to a file in the 'figs' folder without compression
function saveImage(im, name, format)
   
    if(max(max(im)) >100)
        im = uint8(im);
    end

    if strcmp(format,'bmp') || strcmp(format,'jpg')
         imwrite(im, ['figs/' name  '.' format], format) %Save figure
         return;
    end
    imwrite(im, ['figs/' name  '.' format], format,'Compression','none') %Save figure
end