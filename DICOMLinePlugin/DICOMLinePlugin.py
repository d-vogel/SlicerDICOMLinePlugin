import glob, os, json
from datetime import datetime
import string
import vtk, qt, ctk, slicer
from DICOMLib import DICOMPlugin
from DICOMLib import DICOMLoadable
import logging
import numpy as np

#
# This is the plugin to handle translation of DICOM Surface SEG objects
#

class DICOMLinePluginClass(DICOMPlugin):

    def __init__(self,epsilon=0.01):
        super(DICOMLinePluginClass,self).__init__()
        self.loadType = "DICOMSurfaceSegmentation"

        self.surfaceSegmentationSOPClassUID = "1.2.840.10008.5.1.4.1.1.66.5"

        self.tags['seriesInstanceUID'] = "0020,000E"
        self.tags['seriesDescription'] = "0008,103E"
        self.tags['seriesNumber'] = "0020,0011"
        self.tags['modality'] = "0008,0060"
        self.tags['sopClassUID'] = "0008,0016"
        self.tags['instanceUID'] = "0008,0018"

    def examine(self,fileLists):
        """ Returns a list of DICOMLoadable instances
        corresponding to ways of interpreting the
        fileLists parameter.
        """
        loadables = []
        for files in fileLists:
            loadables += self.examineFiles(files)
        return loadables

    def examineFiles(self,files):

        """ Returns a list of DICOMLoadable instances
        corresponding to ways of interpreting the
        files parameter.
        """
        loadables = []

        # just read the modality type; need to go to reporting logic, since DCMTK
        #     is not wrapped ...

        for cFile in files:

            uid = slicer.dicomDatabase.fileValue(cFile, self.tags['instanceUID'])
            if uid == '':
                return []

            desc = slicer.dicomDatabase.fileValue(cFile, self.tags['seriesDescription'])
            if desc == '':
                desc = "Unknown"

            number = slicer.dicomDatabase.fileValue(cFile, self.tags['seriesNumber'])
            if number == '':
                number = "Unknown"

            isSurfaceSeg = (slicer.dicomDatabase.fileValue(cFile, self.tags['sopClassUID']) == self.surfaceSegmentationSOPClassUID)

            if isSurfaceSeg:
                loadable = DICOMLoadable()
                loadable.files = [cFile]
                loadable.name = desc + ' - as a DICOM Surface SEG object'
                loadable.tooltip = loadable.name
                loadable.selected = True
                loadable.confidence = 0.99
                loadable.uid = uid
                self.addReferences(loadable)
                refName = self.referencedSeriesName(loadable)
                if refName != "":
                    loadable.name = refName + " " + desc + " - SurfaceSegmentations"

                loadables.append(loadable)

                print('DICOM Surface SEG found')

        return loadables

    def referencedSeriesName(self,loadable):
        """Returns the default series name for the given loadable"""
        referencedName = "Unnamed Reference"
        if hasattr(loadable, "referencedSeriesUID"):
            referencedName = self.defaultSeriesNodeName(loadable.referencedSeriesUID)
        return referencedName

    def addReferences(self,loadable):
        """Puts a list of the referenced UID into the loadable for use
        in the node if this is loaded."""
        import pydicom as dicom
        dcm = dicom.read_file(loadable.files[0])

        if hasattr(dcm, "ReferencedSeriesSequence"):
            # look up all of the instances in the series, since segmentation frames
            #    may be non-contiguous
            if hasattr(dcm.ReferencedSeriesSequence[0], "SeriesInstanceUID"):
                loadable.referencedInstanceUIDs = []
                for f in slicer.dicomDatabase.filesForSeries(dcm.ReferencedSeriesSequence[0].SeriesInstanceUID):
                    refDCM = dicom.read_file(f)
                    # this is a hack that should probably fixed in Slicer core - not all
                    #    of those instances are truly referenced!
                    loadable.referencedInstanceUIDs.append(refDCM.SOPInstanceUID)
                    loadable.referencedSeriesUID = dcm.ReferencedSeriesSequence[0].SeriesInstanceUID

    def getValuesFromCodeSequence(self, segment, codeSequenceName, defaults=None):
        try:
            cs = segment[codeSequenceName]
            return cs["CodeValue"], cs["CodingSchemeDesignator"], cs["CodeMeaning"]
        except KeyError:
            return defaults if defaults else ['', '', '']

    def load(self,loadable):
        """ Load the DICOM SEG object
        """
        print('DICOM SEG load()')
        try:
            uid = loadable.uid
            print ('in load(): uid = ', uid)
        except AttributeError:
            return False

        import pydicom as dicom

        dataset = dicom.read_file(loadable.files[0])
        self.loadSurfaceSegmentationDataset(dataset)
        return True


    def loadSurfaceSegmentationDataset(self,dataset):
        """This is the functino actually called by slicer for the conversion.

        """
        def surfaceToPointsLPS(surface):
            lineName = dataset.SegmentSequence[0].SegmentLabel

            pointsNb = surface.SurfacePointsSequence[0].NumberOfSurfacePoints
            pointsCoordBin = surface.SurfacePointsSequence[0].PointCoordinatesData
            pointsIdx = surface.SurfaceMeshPrimitivesSequence[0].LineSequence[0].PrimitivePointIndexList
            print(f"surfaceToPoints:: {pointsNb} points found")
            points = [[0,0,0]]*pointsNb
            import struct
            for iPoint in range(pointsNb):
            
                thisPointIdx = struct.unpack('<H', pointsIdx[int(len(pointsIdx)/pointsNb) * iPoint : 
                                                            int(len(pointsIdx)/pointsNb) * (iPoint+1)]
                                            )[0] -1 # the point indexed start at 1
                thisPointCoordBin = pointsCoordBin[int(len(pointsCoordBin)/pointsNb) * iPoint:
                                                int(len(pointsCoordBin)/pointsNb) * (iPoint+1)]
                print(f"surfaceToPoints:: {iPoint} writing point at {thisPointIdx}")
                
                thisPoint = [struct.unpack('f', 
                                                        thisPointCoordBin[int(len(thisPointCoordBin)/3) * i:
                                                                            int(len(thisPointCoordBin)/3) * (i+1)]
                                                            )[0] for i in range(3)]
                # points are encoded RAS in DICOM, slicer expects LPS
                points[thisPointIdx] = [thisPoint[0]*-1,
                                        thisPoint[1]*-1,
                                        thisPoint[2],
                                        ]
            return (lineName, points)
        
        for surface in dataset.SurfaceSequence:
            points = surfaceToPointsLPS(surface)
            print(f"loadSurfaceSegmentationDataset:: points: {points}")
            lineNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsLineNode")
            slicer.util.updateMarkupsControlPointsFromArray(lineNode, np.array(points[1]))
            lineNode.SetName(points[0])

        return True

    def examineForExport(self, node):

        exportable = None

        if node.GetAssociatedNode() and node.GetAssociatedNode().IsA('vtkMRMLSegmentationNode'):

            # Check to make sure all referenced UIDs exist in the database.
            instanceUIDs = node.GetAttribute("DICOM.ReferencedInstanceUIDs").split()
            if instanceUIDs == "":
                    return []

            for instanceUID in instanceUIDs:
                inputDICOMImageFileName = slicer.dicomDatabase.fileForInstance(instanceUID)
                if inputDICOMImageFileName == "":
                    return []

            exportable = slicer.qSlicerDICOMExportable()
            exportable.confidence = 1.0
            exportable.setTag('Modality', 'SEG')

        if exportable is not None:
            exportable.name = self.loadType
            exportable.tooltip = "Create DICOM files from segmentation"
            exportable.nodeID = node.GetID()
            exportable.pluginClass = self.__module__
            # Define common required tags and default values
            exportable.setTag('SeriesDescription', 'No series description')
            exportable.setTag('SeriesNumber', '1')
            return [exportable]

        return []

    def export(self, exportables):

        exportablesCollection = vtk.vtkCollection()
        for exportable in exportables:
            vtkExportable = slicer.vtkSlicerDICOMExportable()
            exportable.copyToVtkExportable(vtkExportable)
            exportablesCollection.AddItem(vtkExportable)

        self.exportAsDICOMSEG(exportablesCollection)

    def exportAsDICOMSEG(self, exportablesCollection):
        """Export the given node to a segmentation object and load it in the
        DICOM database

        This function was copied and modified from the EditUtil.py function of the same name in Slicer.
        """

        import logging

        if hasattr(slicer.modules, 'segmentations'):

            exportable = exportablesCollection.GetItemAsObject(0)
            subjectHierarchyNode = slicer.mrmlScene.GetNodeByID(exportable.GetNodeID())

            instanceUIDs = subjectHierarchyNode.GetAttribute("DICOM.ReferencedInstanceUIDs").split()

            if instanceUIDs == "":
                raise Exception("Editor master node does not have DICOM information")

            # get the list of source DICOM files
            inputDICOMImageFileNames = ""
            for instanceUID in instanceUIDs:
                inputDICOMImageFileNames += slicer.dicomDatabase.fileForInstance(instanceUID) + ","
            inputDICOMImageFileNames = inputDICOMImageFileNames[:-1] # strip last comma

            # save the per-structure volumes in the temp directory
            inputSurfaceSegmentationsFileNames = ""

            import random # TODO: better way to generate temp file names?
            import vtkITK
            writer = vtkITK.vtkITKImageWriter()
            rasToIJKMatrix = vtk.vtkMatrix4x4()

            import vtkSegmentationCore
            import vtkSlicerSegmentationsModuleLogic
            logic = vtkSlicerSegmentationsModuleLogic.vtkSlicerSegmentationsModuleLogic()

            segmentationNode = subjectHierarchyNode.GetAssociatedNode()

            mergedSegmentationImageData = segmentationNode.GetImageData()
            mergedSegmentationLabelmapNode = slicer.vtkMRMLLabelMapVolumeNode()

            segmentationNode.GetRASToIJKMatrix(rasToIJKMatrix)
            mergedSegmentationLabelmapNode.SetRASToIJKMatrix(rasToIJKMatrix)
            mergedSegmentationLabelmapNode.SetAndObserveImageData(mergedSegmentationImageData)
            mergedSegmentationOrientedImageData = logic.CreateOrientedImageDataFromVolumeNode(mergedSegmentationLabelmapNode)

            segmentation = segmentationNode.GetSegmentation()

            segmentIDs = vtk.vtkStringArray()
            segmentation.GetSegmentIDs(segmentIDs)
            segmentationName = segmentationNode.GetName()

            for i in range(0, segmentIDs.GetNumberOfValues()):
                segmentID = segmentIDs.GetValue(i)
                segment = segmentation.GetSegment(segmentID)

                segmentName = segment.GetName()
                structureName = segmentName[len(segmentationName)+1:-1*len('-label')]

                structureFileName = structureName + str(random.randint(0,vtk.VTK_INT_MAX)) + ".nrrd"
                filePath = os.path.join(slicer.app.temporaryPath, structureFileName)
                writer.SetFileName(filePath)

                segmentImageData = segment.GetRepresentation(vtkSegmentationCore.vtkSegmentationConverter.GetSegmentationBinaryLabelmapRepresentationName())
                paddedImageData = vtkSegmentationCore.vtkOrientedImageData()
                vtkSegmentationCore.vtkOrientedImageDataResample.PadImageToContainImage(segmentImageData, mergedSegmentationOrientedImageData, paddedImageData)

                labelmapImageData = slicer.vtkMRMLLabelMapVolumeNode()
                logic.CreateLabelmapVolumeFromOrientedImageData(paddedImageData, labelmapImageData)

                writer.SetInputDataObject(labelmapImageData.GetImageData())

                labelmapImageData.GetRASToIJKMatrix(rasToIJKMatrix)
                writer.SetRasToIJKMatrix(rasToIJKMatrix)
                logging.debug("Saving to %s..." % filePath)
                writer.Write()
                inputSegmentationsFileNames += filePath + ","
            inputSegmentationsFileNames = inputSegmentationsFileNames[:-1] # strip last comma

            # save the per-structure volumes label attributes
            colorNode = segmentationNode.GetNodeReference('colorNodeID')

            terminologyName = colorNode.GetAttribute("TerminologyName")
            colorLogic = slicer.modules.colors.logic()
            if not terminologyName or not colorLogic:
                raise Exception("No terminology or color logic - cannot export")

            inputLabelAttributesFileNames = ""

            for i in range(0, segmentIDs.GetNumberOfValues()):
                segmentID = segmentIDs.GetValue(i)
                segment = segmentation.GetSegment(segmentID)

                segmentName = segment.GetName()
                structureName = segmentName[len(segmentationName)+1:-1*len('-label')]
                labelIndex = colorNode.GetColorIndexByName( structureName )

                rgbColor = [0,]*4
                colorNode.GetColor(labelIndex, rgbColor)
                rgbColor = map(lambda e: e*255., rgbColor)

                # get the attributes and convert to format CodeValue,CodeMeaning,CodingSchemeDesignator
                # or empty strings if not defined
                propertyCategoryWithColons = colorLogic.GetSegmentedPropertyCategory(labelIndex, terminologyName)
                if propertyCategoryWithColons == '':
                    logging.debug ('ERROR: no segmented property category found for label ',str(labelIndex))
                    # Try setting a default as this section is required
                    propertyCategory = "C94970,NCIt,Reference Region"
                else:
                    propertyCategory = propertyCategoryWithColons.replace(':',',')

                propertyTypeWithColons = colorLogic.GetSegmentedPropertyType(labelIndex, terminologyName)
                propertyType = propertyTypeWithColons.replace(':',',')

                propertyTypeModifierWithColons = colorLogic.GetSegmentedPropertyTypeModifier(labelIndex, terminologyName)
                propertyTypeModifier = propertyTypeModifierWithColons.replace(':',',')

                anatomicRegionWithColons = colorLogic.GetAnatomicRegion(labelIndex, terminologyName)
                anatomicRegion = anatomicRegionWithColons.replace(':',',')

                anatomicRegionModifierWithColons = colorLogic.GetAnatomicRegionModifier(labelIndex, terminologyName)
                anatomicRegionModifier = anatomicRegionModifierWithColons.replace(':',',')

                structureFileName = structureName + str(random.randint(0,vtk.VTK_INT_MAX)) + ".info"
                filePath = os.path.join(slicer.app.temporaryPath, structureFileName)

                # EncodeSEG is expecting a file of format:
                # labelNum;SegmentedPropertyCategory:codeValue,codeScheme,codeMeaning;SegmentedPropertyType:v,m,s etc
                attributes = "%d" % labelIndex
                attributes += ";SegmentedPropertyCategory:"+propertyCategory
                if propertyType != "":
                    attributes += ";SegmentedPropertyType:" + propertyType
                if propertyTypeModifier != "":
                    attributes += ";SegmentedPropertyTypeModifier:" + propertyTypeModifier
                if anatomicRegion != "":
                    attributes += ";AnatomicRegion:" + anatomicRegion
                if anatomicRegionModifier != "":
                    attributes += ";AnatomicRegionModifier:" + anatomicRegionModifier
                attributes += ";SegmentAlgorithmType:AUTOMATIC"
                attributes += ";SegmentAlgorithmName:SlicerSelfTest"
                attributes += ";RecommendedDisplayRGBValue:%g,%g,%g" % tuple(rgbColor[:-1])
                fp = open(filePath, "w")
                fp.write(attributes)
                fp.close()
                logging.debug ("filePath: %s", filePath)
                logging.debug ("attributes: %s", attributes)
                inputLabelAttributesFileNames += filePath + ","
            inputLabelAttributesFileNames = inputLabelAttributesFileNames[:-1] # strip last comma'''

            try:
                user = os.environ['USER']
            except KeyError:
                user = "Unspecified"
            segFileName = "editor_export.SEG" + str(random.randint(0,vtk.VTK_INT_MAX)) + ".dcm"
            segFilePath = os.path.join(slicer.app.temporaryPath, segFileName)
            # TODO: define a way to set parameters like description
            # TODO: determine a good series number automatically by looking in the database
            parameters = {
                "inputDICOMImageFileNames": inputDICOMImageFileNames,
                "inputSegmentationsFileNames": inputSegmentationsFileNames,
                "inputLabelAttributesFileNames": inputLabelAttributesFileNames,
                "readerId": user,
                "sessionId": "1",
                "timePointId": "1",
                "seriesDescription": "SlicerEditorSEGExport",
                "seriesNumber": "100",
                "instanceNumber": "1",
                "bodyPart": "HEAD",
                "algorithmDescriptionFileName": "Editor",
                "outputSEGFileName": segFilePath,
                "skipEmptySlices": False,
                "compress": False,
                }

            encodeSEG = slicer.modules.encodeseg
            cliNode = None

            cliNode = slicer.cli.run(encodeSEG, cliNode, parameters, delete_temporary_files=False)
            waitCount = 0
            while cliNode.IsBusy() and waitCount < 20:
                slicer.util.delayDisplay( "Running SEG Encoding... %d" % waitCount, 1000 )
                waitCount += 1

            if cliNode.GetStatusString() != 'Completed':
                raise Exception("encodeSEG CLI did not complete cleanly")

            logging.info("Added segmentation to DICOM database (%s)", segFilePath)
            slicer.dicomDatabase.insert(segFilePath)


#
# DICOMLinePlugin
#

class DICOMLinePlugin:
    """
    This class is the 'hook' for slicer to detect and recognize the plugin
    as a loadable scripted module
    """
    def __init__(self, parent):
        parent.title = "DICOM SurfaceSegmentation Object Import Plugin"
        parent.categories = ["Developer Tools.DICOM Plugins"]
        parent.contributors = ["Andrey Fedorov, BWH"]
        parent.helpText = """
        Plugin to the DICOM Module to parse and load DICOM SEG modality.
        No module interface here, only in the DICOM module
        """
        parent.dependencies = ['DICOM', 'Colors']
        parent.acknowledgementText = """
        This DICOM Plugin was developed by
        Andrey Fedorov, BWH.
        and was partially funded by NIH grant U01CA151261.
        """

        # Add this extension to the DICOM module's list for discovery when the module
        # is created.    Since this module may be discovered before DICOM itself,
        # create the list if it doesn't already exist.
        try:
            slicer.modules.dicomPlugins
        except AttributeError:
            slicer.modules.dicomPlugins = {}
        slicer.modules.dicomPlugins['DICOMLinePlugin'] = DICOMLinePluginClass
