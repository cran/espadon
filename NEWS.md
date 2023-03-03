# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.3.1

## NEW FEATURES:

* study.deployment() to deploy DICOM files from multiple patients. This function
simplifies the analysis of multi-center or multi-expert studies in dosimetry 
challenges, contouring consensus searches, etc.

* vol.repair() to repair missing planes in volumes.

* get.value.from.mesh() to retrieve the values of an object of class "volume" at 
the desired depth of a surface described by a mesh.

* Rdcm.upgrade() to update Rdcm files that were created with a previous version.

* set.reference.obj() to add to an espadon object the information identifying the 
espadon object from which it derives

## ENHANCEMENTS

* The patient.name field, in addition to the PIN, is now handled by espadon 
objects. It impacts all functions that create or load espadon objects, like 
load.patient.from.dicom, load.patient.from.Rdcm(), vol.create(), 
dicom.to.Rdcm.converter(), load_Rdcm_raw_data() ...

* For all functions with the alias argument, when specified, the returned object 
will include the reference fields ref.object.alias and ref.object.info, which 
specify from which objects it was derived. 

* load.patient.from.Rdcm() and load.TMAT() now accept a vector of Rdcm files as argument. 

* New argument ignore.duplicates appears in dicom.to.Rdcm.converter(), 
load.patient.from.dicom() and load.patient.from.Rdcm(). It allows to ignore 
identical DICOM objects (with the same SOP UID), but which have been registered 
with different protocols.

* struct.in.new.ref() no longer recalculates the geometric information on the 
RoI, in case the new reference frame is identical to the old one.

* ref.cut.plane.add() now has new default value for the origin argument.

* get.extreme.pt() now returns the coordinates of the extreme points of the 
espadon objects of the mesh class, the struct class, the volume class, or from 
a binary selection .

* nesting.cube() function restricts the espadon objects of the mesh class, in 
addition to the espadon objects of the volume class, to the cube defined by its 
two extreme vertices.

* rt.gamma.index() and rt.chi.index() now correct the pixel size to match that 
found if the sensor was at the isocenter, in the case where the objects to be 
compared are of rtimage modality.

* get.value.from.xyz() now sends an error if the xyz argument is not a vector of 
length 3 or  a 3-column matrix or dataframe.

* histo.from.roi() now has the over.sampling.factor argument to oversample the 
volume voxels before treatment for a better result.

* the 'dcm' argument of dicom.parser() and dicom.viewer() can now be either a 
dicom file name, or a Rdcm file name or dicom raw data. 

## BUG FIXES:

* orientation.create() now provides orientation composed of unitary vectors.

* display.ob.links() gives error-free information even if the exclusion argument 
is non-NULL.

* save.to.Rdcm() now saves the espadon version.


# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.2.0 

## NEW FEATURES:

* New vol.oversampling() and vol.sub.sampling functions to oversample or 
subsample the grid of a volume class object.

## ENHANCEMENTS
* load.patient.from.dicom() now issues a warning when espadon cannot import 
one of several DICOM objects, instead of crashing.

* get.rigid.M() now returns diag(4) when source reference frame is equal to 
destination reference frame, regardless T.MAT. It therefore has an impact on  
all functions that use a T.MAT to change the reference frame.

* rt.indices.from.roi no longer uses all RoI containing "ptv" as target RoI 
by default.

* rt.gamma.index() and rt.chi.index() now details, in the gamma.info field 
or the chi.info field, the number of dose points, the number of evaluated dose 
points, the rate of evaluated dose points.

## BUG FIXES:

* load.patient.from.dicom() and load.patient.from.Rdcm() no longer fail when 
the dicom.dvh file has no region of interest as reference.

* get.obj.connection() now provides an error-free components matrix.

* display.obj.links() now have a better representation of dicom.dvh

* rt.incices.from.roi, and hence rt.indices.from.bin no longer crash when no 
volume information is requested.

* vol.create() now has a well calculated midpoint.


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




