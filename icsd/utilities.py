import os
import sys


class iCSD3d_Utilities():
    """ Connect to fwd modelling codes
    - read pymgimli formats
    - read resipy formats 
    """
    def __init__(self,dirName):
        self.dirName = dirName
       
        
    def DataImport(self,SimFile=None,ObsFile=None):
        """Data importer for common data files (Resipy and Gimli)
        
        Parameters
        ----------

        """
        
        if fileExt=='*.data':
            print('pygimli format import')
        if fileExt=='*.data':
            print('resipy format import') 
            
    def resipyDataImport(self,SimFile=None,ObsFile=None):

    def pygimliDataImport(self,SimFile=None,ObsFile=None):


    def saveInvData(self, outputdir):
        """Save inverted data
        
        Parameters
        ----------
        outputdir : str
            Path where the .csv files will be saved.
        """


