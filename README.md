# Overview
This repository implements the stone transform, and accompanying demos, for compressive sensing and reconstruction of video.  The implementation is in Matlab, with compilable C-langauge "mex" components.  The STOne transform and its applications are described here.
https://arxiv.org/abs/1311.3405


To get started, run the `setupSTO.m` script.  Then, run a demo script like `testImagereconstruction.m` or `reconstructPendCar.m`.

# Some notes
The stone transform is implemented as a simple linear operator that acts on vectors with power-of-four length.  The transform itself is implemented in the compileable mex file `STO_fast.c`, and will need to be compiled before use.  This compilation should hopefully happen automatically with you call `setupSTO.m`.

To apply the transform to an image, you need to map the image into a vector.  First, call the `createOrderingData.m` method like so...  
```>> order = createOrderingData(256,'semi')```  
or   
```>> order = createOrderingData(256,'full')```   
This will create a struct called `order` that contrains information about the order in which image pixels should get mapped into a vector, assuming a 256x256 image size.  The 'semi' options creates a more regular emebedding patters that produces smoother previews, while the 'full' option creates a highly random embedding that is a little better for compessed sensing, but creates more noisy previews.

Once you have the `order` struct, you can use it to map images into vectors.  A stone transform on an image would look like this...  
```
>> order = createOrderingData(256,'full')  
>> vec = imageToNestedVector(<image>, order)  
>> transformed = STO(vec)  
```

To reverse this transform, one simply does this...
```
>> untransformed = STO(transformed)  
>> image = nestedVectorToImage(untransformed, order)
```
Note this works because the stone transform is unitary.
