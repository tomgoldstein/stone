# Overview
This repository implements the stone transform, and accompanying demos, for compressive sensing and reconstruction of video.  The STOne transform and its applications are described here.
https://arxiv.org/abs/1311.3405


To get started, run the `setupSTO.m` script.  Then, run a demo script like `testImagereconstruction.m` or `reconstructPendCar.m`.

# Some notes
The stone transform is implemented as a simple linear operator that acts on vectors with power-of-four length.  The transform itself is implemented in the max file `STO_fast.c`, and will need to be compiled before use.  This compilation should hopefully happen automatically with you call `setupSTO.m`.

To apply the transform to an image, you first need to map the image into a vector.  First, call the `createOrderingData.m` method like so...  
```>> order = createOrderingData(256,'semi')```  
or   
```>> order = createOrderingData(256,'full')```   
This will create a struct called `order` that contrains information about the order in which image pixels should get mapped into a vector.  The 'semi' options creates a more regular emebedding patters that produces smoother previews, which the 'full' option creates a highly random embedding that is a little better for compessed sensing, but does not preview as well.
