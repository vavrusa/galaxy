# -*- coding: utf-8 -*-

import data
import logging
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.sniff import *
import galaxy.datatypes.util.PyPDB.PyPDB as PyPDB
import galaxy.datatypes.util.Mol2 as Mol2
import galaxy.datatypes.util.sdf as sdf
import dataproviders

log = logging.getLogger(__name__)

@dataproviders.decorators.has_dataproviders
class Molecule ( data.Text ):

    MetadataElement(name="name", default="", desc="Entity", readonly=True, optional=True, visible=True, no_value="")
    MetadataElement(name="chains", default="", desc="Chains", readonly=True, optional=True, visible=True, no_value="")
    MetadataElement(name="models", default="", desc="Number of models", readonly=True, optional=True, visible=True, no_value=1)
    MetadataElement(name="residues", default="", desc="Number of residues", readonly=True, optional=True, visible=True, no_value=0)

class PDB( Molecule ):
    file_ext = "pdb"
    def parse( self, filename ):
        return PyPDB.PDB(os.path.abspath(filename))
    def sniff( self, filename ):
        pdb = self.parse(filename)
        return pdb.isPDB() and len(pdb) > 0
    def set_meta( self, dataset, **kwd ):
        pdb = self.parse(dataset.file_name)
        dataset.metadata.residues = len(pdb)
        dataset.metadata.models = pdb.nModels()
        dataset.metadata.chains = ', '.join(pdb.chnList())
        dataset.metadata.name = pdb.header()
        

class MOL2( Molecule ):
    file_ext = "mol2"
    def parse( self, filename ):
        return Mol2.mol2_set(os.path.abspath(filename))
    def sniff( self, filename ):
        return self.parse(filename).num_compounds > 0
        
class SDF( Molecule ):
    file_ext = "sdf"
    def parse( self, filename ):
        return sdf.SDF_file(os.path.abspath(filename))
    def sniff( self, filename ):
        return len(self.parse(filename).records) > 0

