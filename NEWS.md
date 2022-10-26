# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.1.1                    
  
## NEW FEATURES:

* New struct.clustering() function to create a new volume in which voxels 
are clustered and labeled by Regions of Interest (RoI) defined in an 
rt-struct.

* New get.roi.connection() function to describe the interconnections between 
RoI, from an imaging volume of "cluster" modality.

* New fan.beam(), fan.planar(), fan.sphere() functions to create different 
ray patterns: fan of rays uniformly distributed by angle in a pyramid, or 
fan of rays passing through all pixels of a plane, or fan of rays 
regularly or randomly equidistributed in a sphere. A new espadon class object 
'fan' is created, described in espadon.class() help.

* New fan.to.voxel() function to computes the indices of voxels crossed by a 
fan.

* New nesting.bin() function to restrict a volume to the rectangular 
parallelepiped circumscribed to the selected voxels.

* New vol.median() function to apply a 3D median filter on a volume.

* New functions rt.gamma.index() and rt.chi.index() to compute respectively 
the local or global Gamma and Chi 2D/3D indices from a measurement and a 
reference.
  
* New get.obj.connection() function to describe the connections between the
DICOM objects of the patient.

## ENHANCEMENTS
  
* load.patient.from.dicom(), and hence dicom.to.Rdcm.converter() and
load.object.from.dicom(), increase the information coming from a DICOM 
rt-plan:, in particular, those of the beam control points, the position of 
the source and the direction of the beam for each control point. The new 
information are listed in espadon.class() function. They now also generate more 
useful warnings to inform the user of problems with their data, and avoid 
uninformative warnings.

* display_obj_links() now differentiates the origin (by arrows) and the 
links (by lines).
  
* load.Rdcm.raw.data() now detects espadon version changes.

* get.value.from.ijk() now calculates the value when the volume contains 
  only one plane.
  
* vol.regrid() now allows resampling when volume planes are missing.

* display.3D.contour() to allow the display of RoI (Region of Interest) of 
  type point.
  
* display.plane() now allows the display of volumes with only one plane, the 
display of RoI of type point and the display of overlays with different 
reference frames.

* toy.load.patient() now generates data faster when an rt-dose is requested.

* rt.indices.from.roi() and rt.indices.from.bin() now compute dosimetry index 
D.xcc

## BUG FIXES:

* display.plane() and display.kplane() now make the values outside 
the color palette transparent when sat.transp = TRUE and not when 
FALSE.
  
* load.T.MAT(), and hence load.patient.from.dicom() and 
dicom.to.Rdcm.converter(), no longer make errors even if the patient's source 
DICOM files do not give the same information about their PIN, birthday and sex. 
A simple warning about the uniqueness of the patient is generated.

* vol.gradient() now returns volume with border values at NA.
  
* grid.equal(), which compares the grids of two volumes, now returns TRUE even 
if grids are almost equal (difference due to rounding error).



