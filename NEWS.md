# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.11.3

## BUG FIXES:

* Corrections in RERERENCES.bib formatting

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.11.1

## ENHANCEMENTS

* The qs package, used to create Rdcm files, will soon be deprecated. The espadon 
package now uses qs2, but allows loading of older *.Rdcm files, with load.T.MAT(),
load.obj.from.Rdcm(), load.patient.from.Rdcm(), load.Rdcm.raw.data(). 
Rdcm files can be updated with the Rdcm_upgrade() function.

* DICOM MR export is now possible with the export() function.

## BUG FIXES:

* The export() function now exports the instance number (0020,0013).

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.11.0

## NEW FEATURE:

* struct.update_roiinfo() updates the roi.info element in the espadon object 
of the struct class when point coordinates have been modified in the roi.data 
element.

* err.metrics.from.bin() is added to complement the err.metrics.from.roi() 
function. It is used to calculate errors on binary selections.

## ENHANCEMENTS

* Execution speed of dicom.raw.data.anonymizer() and dicom.patient.anonymiser() 
has been increased.

* Min and max values in err.metrics.from.roi() have been added.

* The dicom.viewer() cursor update has been modified.

## BUG FIXES:

* detection of Rdcm version in the load_Rcm_raw_data() function and all functions 
that use it, such as load.obj.from.Rdcm, load.patient.from.Rdcm

* The error in the load.obj.from.Rdcm() function when loading a single contour 
among others has been corrected.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.10.0

## NEW FEATURE:
* err.metrics.from.roi() calculates various metrics (ME, MAE, MSE, RMSE) to 
compare 2 “volume” class objects in the zones delimited by requested RoI.

* vol.error() and vol.abserror() provide the error and the absolute error volume 
between 2 volumes.

## ENHANCEMENTS
* the toy-patient calculated by toy.load.patient() includes a synthetic CT.

## BUG FIXES:
* nesting.bin() and nesting.cube() corrected for 2D use.

* plot.struct() and display.plane() have been corrected for RoI selection, when 
some RoIs do not contain data.

* get.value.from.ijk() and all functions that depend on it (vol.regrid(), 
bin.from.roi() for example) has been corrected in the counting of subvoxels 
contained in a new voxel (and condition for reaching these voxels) in 'average' 
mode.

* calculation of std value in rt.indices.from.roi() and rt.indices.from.nin() is 
corrected.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.9.0

## NEW FEATURE:

* struct.create  creates a struct object from a list of polygons, representing 
the contours of a shape.

## ENHANCEMENTS

* bin.cuboid(), bin.cylinder() and bin.ellipsoid() have an orientation argument, 
making it easier to tilt shapes. The run time of bin.ellipsoid() is improved.

## BUG FIXES:

* Voxel values calculated by vol.regrid() in “Average” mode are corrected.

* Correction of file name creation error in the export() function.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.8.0

## NEW FEATURES:

* bin.cuboid(), bin.cylinder() and bin.ellipsoid() create espadon objects of 
class “volume”, and of modality “binary” or “weight”, by selecting the voxels 
defining a rectangular cuboid, an elliptical cylinder or an ellipsoid.

* add.shape() adds the shape defined by espadon volume object of the modality 
"binary" or "weight" to a 3D volume.

## ENHANCEMENTS

* The execution time of the bin.opening(), bin.closing(), bin.dilation() and 
bin.erosion() functions has been improved.

* vol.create() now accepts random Gaussian values.

## BUG FIXES:
* nesting.bin() works as described for the volume objects.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.7.4

## BUG FIXES:
* export() now manages espadon object fields with the value NA. RT-struct export
writes the IUD from each CT slice to each ROI contour slice.

* nesting.bin() now takes the xyz.margin argument into account, for mesh class 
objects.

## ENHANCEMENTS
* bin.from.roi() now handles the “weight” modality, in which the value of each 
voxel is its volume fraction included in the ROI. In the “binary” modality, the 
voxel whose fraction is greater than 0.5 is set to TRUE. The emergence of this 
new modality impacts histo.from.roi(), histo.from.bin(), rt.indices.from.roi() 
and rt.indices.from.bin() functions. The use of the “weight” modality enables us 
to get closer to the ROI contours. This slightly modifies the results compared 
with the previous version. 

* histo.from.roi() and histo.from.bin() return a dataframe of dose-volume 
histograms (in percent) for ROIs if requested.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.7.2

## ENHANCEMENTS
* nesting.bin() now works on meshes.

* bin.closing(), bin.erosion(), bin.dilation(), bin.opening (), nesting.roi(), 
add.margin(),get.value.from.ijk now stop instead of warning and return NULL 

* enhancement of struct.clustering example.

## BUG FIXES:
* display.plane() now stops directly on error. It takes into account the argument 
'struct.dxyz'.

* vol.subsampling(): The xyz-ranges error has been corrected.

* plot.mesh() can now display opened meshes

* plot.volume(), in case of add=FALSE, has its map range rectified.

* load.obj.from.dicom() and load.patient.from.dicom() no longer uses the tag (
0020,0012) to separate DICOM objects

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.7.0

## NEW FEATURES:
* dicom.patient.anonymiser() anonymises all DICOM files in a patient's directory.

* export() function exports struct class objects and  volume class objects with 
CT or RTDOSE modality in DICOM format.

## ENHANCEMENTS
* load.obj.from.dicom() and load.obj.from.Rdcm now manage cardiac phase.

## BUG FIXES:
* bin.from.roi() : back to the more robust version 1.5.1.

* load.obj.from.dicom(), load.obj.from.Rdcm() and load.obj.data() modify their 
strategy for calculating the thickness between polygons in DICOM rt-struct files: 
all z-slices in each ROI are taken into account to calculate all thicknesses. 
The resulting thickness will be the median of all these thicknesses.

* dicom_set_tag_value(): correction of value acceptance for tags with "CS" VR.


# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.6.0

## NEW FEATURES:

* plot() function to plot in 2D espadon objects of class 'volume', 'struct' and
'mesh'.

* pal.rainbow() to produce a color palette suitable for displaying the dose in 
overlay.

## ENHANCEMENTS

* display.plane() runs faster in sagital and frontal view, and accepts ...
arguments. If the color palette is the RVV palette, breaks are automatically
calculated between -1000 and 1000, as with the display.kplane function.

* bin.from.roi() runs faster. This affects the functions histo.from.roi(), 
struct.clustering() and rt.indices.from.roi().

* sp.similarity.from.mesh() now computes surface Added Path Length. 

## BUG FIXES:

* load.patient.from.dicom(), load.obj.from.dicom() and dicom.to.Rdcm.converter()
now correctly handle the 'scan spot meterset weight' of rt ion plan modality.

* sp.similarity.from.mesh() now computes correctly sDSC metrics. 

* struct.from.mesh() now correctly calculates complex mesh contours, which 
depend on curvature.


# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.5.1

## NEW FEATURES:

* sp.similarity.from.bin() to calculate spatial similarity such as volumetric Dice 
similarity or Dice-Jaccard coefficient, from binary volumes

* sp.similarity.from.mesh() to calculate spatial similarity such as Hausdorff 
distances and surface Dice similarity coefficient, from meshes.

## BUG FIXES:

* display.kplane() and display.plane() now display images correctly with MAC OS.

* nesting.cube() : Rectification of points included in the cube and located at 
its edges.

* load.patient.from.dicom(), load.obj.from.dicom() and dicom.to.Rdcm.converter()
now handle special characters as soon as DICOM files are decrypted. The functions 
load.Rdcm.raw.data() and Rdcm.upgrade() functions have been updated accordingly.

* load.patient.from.dicom() and dicom.to.Rdcm.converter() now take into account 
all rtdose even those without tag (300C,0006). The import of tag content 
(3006,00A6) from rtstruct files has been corrected.

* dicom_set_tag_value() now accepts "IS" DICOM value representation with many 
integers.

* struct.from.bin() and struct.from.roi() create a espadon struct object with 
roi.obs$label equal to "".

* get.extreme.pt() now takes the max argument into account.

## ENHANCEMENTS
* display.plane() can now display the RoI legend with Roi pseudonyms or RoI names 
depending on the legend.roi.pseudo argument.

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.4.1

## NEW FEATURES:

* obj.create() to create an espadon object with the essential properties it must 
have.

* struct.from.mesh() to create a struct object with a unique RoI, defined by the 
contours of a mesh.

## ENHANCEMENTS

* display.plane() displays now point or open.planar contour, provided that all 
points of the contour are located in the same plane.

* mesh.from.bin() now has default tol argument equal to min(bin$dxyz)/2, to 
limit meshing errors.

* ref.srctodest.add() only warns and does not stop when the link between src.ref 
and dest.ref is already defined in the T.MAT object.

## BUG FIXES:

* vol.regrid() : The voxels of the volume edges are calculated.

* vol.oversampling() now makes the center of the volume invariant.

* get.extreme.pt() : management of NA value if obj is an rtstruct.

* struct.merge() : roi.obs field is reducted to selected RoI. The RoI number is 
incremented in roi.obs and roi.info fields
                   
## OTHERS
* package sp is not longer used
* new link to cite espadon package

# SIGNIFICANT USER-VISIBLE CHANGES IN espadon 1.3.2


## BUG FIXES:

* gamma.index() : the 1 voxel shift has been corrected

* rtplan loading : correction of the calculation of the returned value 
$beam.orientation for functions that can load an object of class rtplan, like 
load.data.from.patient(),load.patient.from.Rdcm(), load.obj.from.dicom(), 
load.obj.from.Rdcm() and load.obj.data().

* struct loading : All functions that can load an object of class struct, like 
load.data.from.patient(),load.patient.from.Rdcm(), load.obj.from.dicom(), 
load.obj.from.Rdcm() and load.obj.data() now calculate better the thickness 
between planes of the contours when the link with the reference imaging object 
has not been done. When this link is made, the thickness is now equal to the 
interplan of the reference object. This change impacts the calculation of volumes 
from contours, of DVH, and rt.indices.from.roi.

* vol.oversampling() and vol.subsampling() handles now the modality 'binary' 

* display.obj.links() : correction of the color management of the different
frames of reference.

## ENHANCEMENTS

* study.deployment() have now pid.prefix to add a prefix to the patient ID 
generated in the deployed DICOM files (tag(0010,0020)).


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




