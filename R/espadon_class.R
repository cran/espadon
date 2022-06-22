#' ESPADON class

#' @return Returns a vector of \pkg{espadon} class names. 
#' @note Each object of a class has specific features that are used to display or 
#' process that object. 
#' 
#' @note \emph{- the "patient" class includes :}
#' \itemize{
#' \item \code{$patient} : dataframe providing patient's information as PIN, 
#' birth date and gender.
#' \item \code{$pat.pseudo} : patient's pseudonym, initialized to the patient's 
#' PIN of \code{$patient} dataframe.
#' \item \code{$description} : dataframe describing the patient's DICOM objects: 
#' their modality (rtstruct, ct, mr, rtplan ...), the base name of the relevant 
#' source file in the patient's directory, the pseudonym of their frame of reference 
#' (ref1, ref2 ...), their number of sub-objects, their description if any, their 
#' numbers of slices/RoI for all sub-objects, their maximum voxels (for volume 
#' sub-objects), and finally the aliases of the sub-objects.
#' \item \code{$description.by.reg}: list of DICOM objects descriptions that are 
#' linked by a transfer matrix.
#' \item \code{$T.MAT} : list of class "t.mat" containing the information of the 
#' transfer matrices to move from one frame of reference to another. 
#' See \link[espadon]{load.T.MAT}.
#' \item \code{$ct} : list of CT, if any. They are named by their \code{$object.alias}
#' See \link[espadon]{load.obj.from.Rdcm}.
#' \item \code{$mr} : list of MRI, if any. They are formatted like the \code{$ct}.
#' \item \code{$rtdose} : list of dose matrices. They are formatted like the \code{$ct}.
#' \item \code{$rtstruct}: list of struct objects.
#' \item ...any DICOM objects other than the reg files, and those previously mentioned, or
#' any modalities created by \pkg{espadon}. 
#' \item \code{$dicom.dvh}: if any, list of DVH computed in rt-dose DICOM files.
#' }
#' 
#' @note \emph{- the "t.mat" class includes :}
#' \itemize{
#' \item \code{$ref.info}: dataframe giving the correspondence between the frame of
#' reference (column \code{$ref}) of the DICOM object (TAG (0020,0052) ) and its
#' pseudonym (column \code{$ref_pseudo}).
#' \item \code{$reg.info}: list of dataframes : the first one gives the PID, 
#' birthday, and sex of the patient, the second one gives the name of the source 
#' file of transfer matrices.
#' \item \code{$matrix.description}: dataframe giving the transfer matrix names
#' (column \code{$t}), its source frame of reference (column \code{$src}), the 
#' destination frame of reference (column \code{$dest}), and its type (\code{$type}). 
#' Note that only the RIGID type is supported.
#' \item \code{$matrix.list}: list of 4X4 transfer matrices. This list contains
#' at least as many Identity matrices as there are \code{ref.pseudo}.
#' }
#' 
#' @note  A \pkg{espadon} object  of class "dvh","histo","histo2D","mesh", "reg", 
#' "struct", "t.mat","undef","volume" is a list containing at least:
#' \itemize{
#' \item \code{$patient}: patient's PIN.
#' \item \code{$patient.bd}: patient's birthday.
#' \item \code{$patient.sex}: patient's sex
#' \item \code{$file.basename}: vector of .Rdcm or .dcm file basenames of the 
#' object, if it exists
#' \item \code{$file.dirname }: directory including the  .Rdcm or .dcm file, 
#' if it exists
#' \item \code{$object.name}: name of the object.  
#' \item \code{$object.alias}: alias of the object.
#' \item \code{$frame.of.reference}: value of TAG (0020,0052).
#' \item \code{$ref.pseudo}: pseudonym of the \code{$frame.of.reference}
#' \item \code{$modality}: modality of the object (e.g. ct, mr, bin, rtplan..)
#' \item \code{$description}: description of the object.
#' \item \code{$creation.date }: creation date of the object.
#'}
#' @note If the object was generated from a DICOM file, the list also contains:
#' \itemize{
#' \item \code{$object.info}: Information of the object. It includes:
#' \tabular{rl}{
#' -\tab the SOP ID (value of TAG (0008,0016)),\cr
#' -\tab the transfer syntax UID (value of TAG (0002,0010)),\cr
#' -\tab the SOP implementation ID (value of TAG (0002,0012)),\cr
#' -\tab the SOP type (value of TAG (0008,0008)),\cr
#' -\tab the study ID (value of TAG (0020,0010)),\cr
#' -\tab the study UID (value of TAG (0020,000D)),\cr
#' -\tab the serie UID (value of TAG (0020,000E)),\cr
#' -\tab the scanning sequence (value of TAG (0018,0020)),\cr
#' -\tab the list of SOP labels (values of TAG (0008,0018)),\cr
#' -\tab the dicom source files,\cr
#' -\tab the encoding of the content of text tags (values of TAG (0008,0005)) and\cr 
#' -\tab the number of sub-objects.
#' }
#' }
#' 
#' @note if the object is linked to another DICOM object, the list also contains:
#' \itemize{
#' \item \code{$ref.object.name}: Name of the reference object. Available only for 
#' rtdose.
#' \item \code{$ref.object.info}: Information of the reference object (not available 
#' for mr and ct). It includes:
#' \tabular{rl}{
#' -\tab the SOP ID of the reference object,\cr
#' -\tab the list of SOP names of the reference object.
#' }
#' }
#' 
#' 
#' @note \emph{- the "dvh" class also includes :}
#' \itemize{
#' \item \code{$nb.MC}: set to \code{histo$nb.MC}.
#' \item \code{$breaks}: vector breakpoints.
#' \item \code{$mids}: vector of cell centers.
#' \item \code{$mids.unit}: Character string, representing the unit of the abcissa
#' of the DVH. For instance, "Gy".
#' \item \code{$vol}: cumulative volume receiving at least the doses defined by \code{$mids}.
#' \item \code{$pcv}: percentage of the total volume receiving at least the doses defined by \code{$mids}.
#' \item if \code{$nb.MC} is different from 0, the arrays \code{MC.vol}, \code{MC.pcv} and 
#' \code{MC.dxyz} are added. See \link[espadon]{histo.DVH}.
#' } 
#'
#' @note \emph{- the "histo" class also includes :}
#' \itemize{
#' \item \code{$nb.MC}: number of Monte-Carlo simulations
#' \item \code{$breaks}: vector breakpoints
#' \item \code{$mids}: vector of cell centers.
#' \item \code{$mids.unit}: Character string, representing the unit of the abcissa
#' of the histogram. For instance, "Gy".
#' \item \code{counts}: count of voxels whose value is included in the limits 
#' defined by \code{$breaks}.
#' \item \code{dV_dx}: differential histogram, expressed in \mjeqn{cm^3}{ascii} by voxel units, 
#' at each \code{$mids}.
#' \item if \code{$nb.MC} is different from 0, the arrays \code{MC.counts}, \code{MC.dV_dx} and 
#' \code{MC.dxyz} are added. See \link[espadon]{histo.from.roi}.
#' }
#' 
#' @note \emph{- the "histo2D" class also includes :}
#' \itemize{
#' \item \code{$nb.pixels}: number of elements in the \code{density.map}.
#' \item \code{$x.file.src}: x label. See \link[espadon]{histo.2D}.
#' \item \code{$y.file.src}: y label. See \link[espadon]{histo.2D}.
#' \item \code{x.breaks}: vector of x-axis breakpoints.
#' \item \code{y.breaks}: vector of y-axis breakpoints.
#' \item \code{x.mids}: vector of x-axis cell centers.
#' \item \code{y.mids}: vector of y-axis cell centers.
#' \item \code{density.map}: array of densities.
#' \item \code{total.counts}: number of counted voxels.
#' }
#' @note \emph{- the "mesh" class also includes :} 
#' \itemize{
#' \item \code{$nb.faces}: set to the number of faces of the mesh.
#' \item \code{$mesh}: list of 3 elements defining the mesh :\code{$vb}, \code{$it}
#' and \code{$normals}. See \link[espadon]{mesh.from.bin}.
#' }
#' 
#' @note \emph{- the "reg" class also includes :}
#' \itemize{
#' \item \code{$nb.of.ref}: number of transfer matrices.
#' \item \code{$ref.data}: list including the lists of information on transfer
#' matrices, namely: the source frame of reference (\code{$src}), the matrix type
#' (\code{$type}, for example 'RIGID') and the transfer matrix (\code{$matrix}).
#' }
#' 
#' 
#' @note {- the "rtplan" class also includes :}
#' \itemize{
#' \item \code{$approval.status}: value of TAG (300E,0002).
#' \item \code{$number}: sub-object number.
#' \item \code{$plan.info}: dataframe containing, if they exist, 
#' \tabular{rl}{
#' \tab - \code{$label} the label for the treatment plan,\cr  
#' \tab - \code{$plan.name} the name for the treatment plan,\cr  
#' \tab - \code{$plan.description} description of treatment plan,\cr  
#' \tab - \code{$tt.protocol} the treatment protocol,\cr 
#' \tab - \code{$plan.intent} the intent of this plan,\cr  
#' \tab - \code{$tt.site} describing the anatomical treatment site,\cr  
#' \tab - \code{$geometry} describing whether RT Plan is based on patient or 
#' treatment device geometry.
#' }
#' 
#' \item \code{$presc.dose}: dataframe containing, if they exist,
#' \tabular{rl}{
#' \tab - \code{$ref.roi.nb} value of TAG (3006,0084),\cr 
#' \tab - \code{$dose.ref.nb} value of TAG (300A,0012),\cr
#' \tab - \code{$dose.ref.id} value of TAG (300A,0013),\cr 
#' \tab - \code{$struct.type} value of TAG (300A,0014),\cr
#' \tab - \code{$description} value of TAG (300A,0016),\cr 
#' \tab - \code{$pt.coord} value of TAG (300A,0018),\cr
#' \tab - \code{$nominal.prior.dose} value of TAG (300A,001A),\cr 
#' \tab - \code{$dose.type} value of TAG (300A,0020),\cr
#' \tab - \code{$constraint.weight} value of TAG (300A,0021),\cr 
#' \tab - \code{$deliv.warn.dose} value of TAG (300A,0022),\cr
#' \tab - \code{$deliv.max.dose} value of TAG (300A,0023),\cr 
#' \tab - \code{$targ.min.dose} value of TAG (300A,0025),\cr
#' \tab - \code{$targ.presc.dose} value of TAG (300A,0026),\cr 
#' \tab - \code{$targ.max.dose} value of TAG (300A,0027),\cr
#' \tab - \code{$targ.underdose.vol.frac} value of TAG (300A,0028),\cr 
#' \tab - \code{$org.risk.full.vol.dose} value of TAG (300A,002A),\cr
#' \tab - \code{$org.risk.lim.dose} value of TAG (300A,002B),\cr 
#' \tab - \code{$org.risk.max.dose} value of TAG (300A,002C),\cr
#' \tab - \code{$org.risk.overdose.vol.frac} value of TAG (300A,002D)
#' }

#' \item \code{$fraction.info}: dataframe containing, if they exist,
#' \tabular{rl}{
#' \tab - \code{$fraction.id} the id of the fraction group,\cr
#' \tab - \code{$description} its description,\cr
#' \tab - \code{$planned.frac.nb} the total number of treatments (Fractions) 
#' prescribed for current fraction group,\cr
#' \tab - \code{$frac.pattern.digit.per.day.nb} the number of digits in \code{$frac.pattern} 
#' used to represent one day,\cr
#' \tab - \code{$repeat.frac.cycle.le} the number of weeks needed to describe 
#' treatment pattern,\cr
#' \tab - \code{$frac.pattern} the value of TAG (300A,007B) describing treatment 
#' pattern every day,\cr
#' \tab - \code{$beam.nb} the number of beams in current fraction group,\cr
#' \tab - \code{$beam.dose.meaning} the value of TAG (300A,008B) indicating the 
#' meaning of Beam Dose,\cr
#' \tab - \code{$brachy.app.nb} the number of brachy application setups in current 
#' fraction group.
#' }

#' \item \code{$beam.info} (in case of beam treatment): dataframe containing, if 
#' they exist,
#' \tabular{rl}{
#' \tab - \code{$fraction.id},\cr
#' \tab - \code{$planned.frac.nb},\cr
#' \tab - \code{$beam.dose} the value of TAG (00A,0084),\cr
#' \tab - \code{$beam.specif.pt} the value of TAG (300A,0082),\cr
#' \tab - \code{$beam.meterset} the value of TAG (300A,0086),\cr
#' \tab - \code{$beam.type} the value of TAG (300A,0090,\cr
#' \tab - \code{$alt.dose} the value of TAG (300A,0091),\cr
#' \tab - \code{$alt.type} the value of TAG (300A,0092,\cr
#' \tab - \code{$duration.lim} the value of TAG (300A,00C5),\cr
#' \tab - \code{$beam.nb} the value of TAG (300C,0006) or (300A,00C0),\cr
#' \tab - \code{$beam.name} the value of TAG (300A,00C2),\cr
#' \tab - \code{$beam.description} the value of TAG (300A,00C3),\cr
#' \tab - \code{$beam.type} the value of TAG (300A,00C4),\cr 
#' \tab - \code{$radiation.type} the value of TAG (300A,00C6),\cr 
#' \tab - \code{$high.dose.technique.type} the value of TAG (300A,00C7),\cr 
#' \tab - \code{$treatment.machine.name} the value of TAG (300A,00B2),\cr 
#' \tab - \code{$device.serial.nb} the value of TAG (0018,1000),\cr                                 
#' \tab - \code{$primary.dosimeter.unit} the value of TAG (300A,00B3),\cr
#' \tab - \code{$referenced.tolerance.table.nb} the value of TAG (300C,00A0),\cr                   
#' \tab - \code{$src.axis.dist} the value of TAG (300A,00B4),\cr
#' \tab - \code{$referenced.patient.setup.nb} the value of TAG (300C,006A),\cr
#' \tab - \code{$treatment.delivery.type} the value of TAG (300A,00CE),\cr
#' \tab - \code{$wedges.nb} the value of TAG (300A,00D0),\cr                                     
#' \tab - \code{$compensators.nb} the value of TAG (300A,00E0),\cr
#' \tab - \code{$total.compensator.tray.factor} the value of TAG (300A,00E2),\cr                    
#' \tab - \code{$boli.nb} the value of TAG (300A,00ED),\cr      
#' \tab - \code{$blocks.nb} the value of TAG (300A,00F0),\cr       
#' \tab - \code{$total.block.tray.factor} the value of TAG (300A,00F2),\cr      
#' \tab - \code{$final.cumul.meterset.weight} the value of TAG (300A,010E),\cr                      
#' \tab - \code{$ctl.pts.nb} the value of TAG (300A,0110),\cr      
#' \tab - \code{$radiation.mass.nb} the value of TAG (300A,0302),\cr      
#' \tab - \code{$radiation.atomic.nb} the value of TAG (300A,0304),\cr      
#' \tab - \code{$radiation.charge.state} the value of TAG (300A,0306),\cr                           
#' \tab - \code{$scan.mode} the value of TAG (300A,0308),\cr      
#' \tab - \code{$modulated.scan.mode.type} the value of TAG (300A,0309),\cr       
#' \tab - \code{$virtual.src.axis.dist} the value of TAG (300A,030A),\cr      
#' \tab - \code{$total.wedge.tray.water.equ.thickness} the value of TAG (300A,00D7),\cr      
#' \tab - \code{$total.compensator.tray.water.equ.thickness} the value of TAG (300A,02E3),\cr      
#' \tab - \code{$total.block.tray.water.equ.thickness} the value of TAG (300A,00F3),\cr      
#' \tab - \code{$range.shifters.nb} the value of TAG (300A,0312),\cr      
#' \tab - \code{$lateral.spreading.devices.nb} the value of TAG (300A,0330),\cr 
#' \tab - \code{$range.modulators.nb} the value of TAG (300A,0340),\cr      
#' \tab - \code{$fixation.light.azimuthal.angle} the value of TAG (300A,0356),\cr                   
#' \tab - \code{$fixation.light.polar.angle} the value of TAG (300A,0358).
#' }
#' 
#' \item \code{$brachy.info} (in case of brachy treatment): dataframe containing, if they exist,
#' \tabular{rl}{
#' \tab - \code{$fraction.id}\cr
#' \tab - \code{$planned.frac.nb},\cr
#' \tab - \code{$brachy.dose} the value of TAG (300A,00A4),\cr
#' \tab - \code{$brachy.nb} the value of TAG (300C,000C),\cr
#' \tab - \code{$brachy.specif.pt} the value of TAG (300A,00A).
#' }
#' }
#' @note \emph{- the "struct" class also includes :}
#' \itemize{
#' \item \code{$nb.of.roi}: number of regions of interest (RoI).
#' \item \code{$thickness}: thickness between two consecutive planes of a contour.
#' \item \code{$ref.from.contour}: reference frame change matrix, from the contour 
#' reference frame to the ref.pseudo reference frame
#' \item \code{$roi.info}: dataframe. Information on RoI contours. See 
#' \link[espadon]{load.obj.from.Rdcm}
#' \item \code{$roi.data}: exists only if the data is loaded. Contains the list 
#' of contour coordinates. The RoI of list number i is that of line i of roi.info.
#' Each element of the list is a list giving the contour 
#' information for each plane, namely:
#' \tabular{rl}{
#' \tab - \code{$type}: value of TAG (3006,0042).\cr
#' \tab - \code{$pt}: dataframe of the coordinates of the contour points.\cr 
#'  \tab If the contour is closed (i.e.\code{$type = "CLOSED_PLANAR"}), \cr
#'  \tab then the first point is repeated at the end.\cr
#' \tab - \code{$level}: contour inclusion level. If this number is even,\cr
#' \tab the inside of the closed contour belongs to the RoI.\cr
#' \tab Otherwise, if odd, the inside of the closed contour is excluded from the RoI.
#' }
#' }
#' 
#' @note \emph{- the "undef" class :} is used for DICOM objects that will not be 
#' processed further by \pkg{espadon} functions. It can 
#' also include what the user wants.  
#' 
#' 
#' @note \emph{- the "volume" class also includes :}
#' \itemize{
#' \item \code{$number}: sub-object number.
#' \item \code{$n.ijk}: vector defining the number of indices i, j, k. The product
#'  \code{prod(...$n.ijk)} represents the number of voxels in the 3D volume.
#' \item \code{$slice.thickness}: thickness in mm of a plane.
#' \item \code{$min.pixel}: minimum value of voxels in the volume.
#' \item \code{$max pixel}: maximum value of voxels in the volume.
#' \item \code{$dxyz}: x, y, z steps in mm.
#' \item \code{$patient.orientation}: value of TAG (0020,0037).
#' Vector, comprising the vectors i and j defining the orientation of the patient
#' with respect to the volume planes.
#' \item \code{$patient.xyz0}: in the patient frame of reference, position of the
#' first voxel of each plane.
#' \item \code{$xyz.from.ijk}: transfer matrix of the voxels i, j, k indices to
#' the position x, y, z in mm in the patient's frame of reference.
#' \item \code{$k.idx}: index of planes in the 3D volume.
#' \item \code{$missing.k.idx}: Boolean, indicating if k is a continuous sequence of integers.
#' \item \code{$cube.idx}: 3D volume vertices indices.
#' \item \code{$vol3D.data}: exists only if the data is loaded. 3D array of the voxel
#'  values of the 3D volume.
#'  }
#' 
#' @seealso \link[espadon]{toy.load.patient}, \link[espadon]{load.patient.from.dicom},
#' \link[espadon]{load.patient.from.Rdcm}, \link[espadon]{load.T.MAT}
#' \link[espadon]{histo.DVH}, \link[espadon]{histo.vol}, 
#' \link[espadon]{histo.from.roi}, \link[espadon]{histo.from.bin},
#' \link[espadon]{histo.2D}, \link[espadon]{mesh.from.bin}, 
#' \link[espadon]{load.obj.from.Rdcm}
 
#' @examples
#' cat ("espadon class names are:", paste (espadon.class(), collapse = ", "))
#' @export
espadon.class <- function( ) {
return(c("dvh","histo","histo2D","mesh","patient", "reg", "rtplan",
         "struct", "t.mat","undef","volume"))
}