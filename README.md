# Fingerprint matching based on minutiae-triplets

This repository provides an implementation of following works in C# using EmguCV.

  1. ***Directionally Weighted Distance for Minutiae-Triplets Preservation on Elastic Deformation of Fingerprint Recognition***, please refer to [Rungchokanun2022](https://www.sciencedirect.com/science/article/abs/pii/S0167865522001878). This work is implemented as function in ```class program```, see function ```exampleMinutiaeTripletsFormation()``` for detail.
  
  2. ***Fingerprint minutiae matching based on the local and global structures***, please refer to [Jiang2000](https://ieeexplore.ieee.org/abstract/document/906252). This work is implemented as ```class MTriplet_JY2000``` and ```class MTriplet_JY2000V2```. The differences between these two versions are how the angle differences are calculated, the ```MTriplet_JY2000V2``` is prefer. See function ```exampleJYMatching()``` for detail.
  
  3. ***Improving Fingerprint Verification Using Minutiae Triplets***, please refer to [Medina-Pérez2012](https://www.mdpi.com/1424-8220/12/3/3418). This work is implemented as ```class M3gl_Perez2012```, See function ```exampleM3glMatching()``` for detail.
  
  4. ***Fingerprint Reference Point Detection Based on Local Axial Symmetry***, please refer to [Lui2006](https://ieeexplore.ieee.org/document/1699069). This work is implemented as ```class LAS_Liu2006```, See function ```exampleLAS()``` for detail.
  
  5. ***Localization of corresponding points in fingerprints by complex filtering***, please refer to [Nilsson2003](https://www.sciencedirect.com/science/article/abs/pii/S0167865503000837). This work is implemented as ```class ComplexFilter_Nilson2003```, See function ```exampleComplexFilter()``` for detail.
  
  6. ***Systematic Methods for the Computation of the Directional Fields and Singular Points of Fingerprints***, please refer to [Bazen2002](https://ieeexplore.ieee.org/abstract/document/1017618). The orientation field estimation part is implemented as ```class OF_Bazen2002``` and the singular point detection is ignored, See function ```exampleOrientationField()``` for detail.

This solution contains two projects which are DirectionallyWeightedDistance_Rungchokanun2022 and Graph.
The former contains the above implementation and the latter is used as a helper class for plotting matched minutiae. 

### Please note that the source code is implemented for the study #1. Directionally Weighted Distance for Minutiae-Triplets Preservation on Elastic Deformation of Fingerprint Recognition".

#### Usage
Unzip and rename the folder to make it shorter, otherwise you will not be able to download packages.
Open the solution with Microsoft Visual Studio and 'Rebuild' the solution to gather all required packages.
If you encounter an error from Emgu CV that is “Emgu CV is not able to detect the project type", you can unload the ```DirectionallyWeightedDistance_Rungchokanun2022``` project and reload it.

#### Requirement
EmguCV 4.2.0.3662, .NET Framework 4.7.2
