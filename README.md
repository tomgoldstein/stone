# Overview
This repository implements the stone transform, and accompanying demos, for compressive sensing and reconstruction of video.  The implementation is in Matlab, with compilable C-langauge "mex" components.  The STOne transform and its applications are described here.
https://arxiv.org/abs/1311.3405


To get started, run the `setupSTO.m` script.  Then, run a demo script like one of these...

- `testImageReconstruction.m`: compressed-sensing recovery of a synthetic image using TV prior  
- `reconstructPendCar.m`: compressed-sensing recovery of "PendCar" video frames 
- `testPreview.m`: create preview images from synthetic datasets using a direct method with no optimization
- `reconstructDMDVideosFromDatasets.m`:  recovery of videos from real emprical DMD measurements



# Some notes
The stone transform is implemented as a simple linear operator that acts on vectors with power-of-four length.  The transform itself is implemented in the compileable mex file `STO_fast.c`, and will need to be compiled before use.  This compilation should hopefully happen automatically with you call `setupSTO.m`.

To apply the transform to an image, you need to map the image into a vector.  First, call the `createOrderingData.m` method like so...  
```>> order = createOrderingData(256,'semi')```  
or   
```>> order = createOrderingData(256,'full')```   
This will create a struct called `order` that contrains information about the order in which image pixels should get mapped into a vector, assuming a 256x256 image size.  The 'semi' options creates a more regular emebedding patters that produces smoother previews, while the 'full' option creates a highly random embedding that is a little better for compessed sensing, but creates more noisy previews.

Once you have the `order` struct, you can use it to map images into vectors.  A stone transform on an image would look like this...  
```matlab
>> image = phantom(256);                      % Create 256x256 synthetic image . 
>> order = createOrderingData(256,'full');    % Create a struct with information on how to map image into vector  
>> vec = imageToNestedVector(image, order);   % Map the image into a vector using the specified pixel ordering    
>> transformed = STO(vec);                    % Apply the stone transform to the vector
```

To reverse this transform, one simply does this...
```matlab
>> untransformed = STO(transformed);          % The stone transform is self-adjoint, so it's its own inverse
>> image = nestedVectorToImage(untransformed, order); % Map the 1D vector of pixels back into a 2D image array
```
Note this works because the stone transform is unitary.
