function makeGifFromFrames(frames, filename)
imdata = permute(frames,[1 2 4 3]);
imwrite(imdata,[filename '.gif'],'DelayTime',0,'LoopCount',inf);
