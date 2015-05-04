#
# A parser for PDB files
#
# Written (2001-2003) by P. Tuffery, INSERM, France
# Contributions by R. Gautier, J. Maupetit, J. Herisson, F. Briand
# Version: 10.0 (2010 June)
#
# No warranty of any kind is provided
# This is free software. You can use it, modify it, distribute it
# but its origin (i.e. this present text) must remain clearly
# stated.
#
# Thanks for any feedback related to any bug fix, or any improvement.
#
#
# Main classe:
#   PDB (for one PDB file)
# Main function:
#   PDB() to obtain a PDB instance from a file
#   PDBList() to obtain a collection of PDB instance
#

"""
A PDB class to parse a PDB file

Written (2001-2003) by P. Tuffery, INSERM, France
Contributions by R. Gautier, J. Maupetit, J. Herisson, F. Briand
Version: 10.0 (2010 June)

The idea is to have an easy management of PDB files

Ex:

x = PDB("1tim")
y = PDB("1timA")

# 1st residue of x \n
x[0]

# Chain A of x \n
x["A"]

#1st atom of the 2nd residue of x\n
x[0][1]


Classes' Pattern:

PDBLine (a line of a PDB file)
atmLine (a line ATOM/HETATM of a PDB file)
atmList (a list of lines ATOM/HETATM of a PDB file (ex: a list of atoms lines))
residu (a residue of a file (atomic information))
PDB (a PDB)

PDBList (function, not a class for the moment) manages a collection of PDBs.

protein : Class to manage a protein type PDB
(it is not finish at present (ongoing adjustments), but it is partially functional)
"""

# NOTE:
# (The development has been incremental, and a facelift seems to be
# now necessary.)

import string
import sys
import os
import copy
import math
import gzip
import types
import urllib
import urllib2
import subprocess
import tempfile
import stat

## sys.path.append('/home/raid5/PyTools/Classes/')
## sys.path.append('/home/tuffery/proteineDBTools/PyScripts/')
sys.path.append("/data/PyTools/Classes/")

from FileBasics import *
## from Lines import *
from Geo3DUtils import *
from html2text import *

from Config import *

AA1 = "ACDEFGHIKLMNPQRSTVWY"
AA3 = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HSE", "HSD","HSP"]
AA1seq = "ACDEFGHIKLMNPQRSTVWYXXXSXWMCXWYMDPECXXYAMMSCMAHHH"
AA3STRICT = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

AA3new = ['PAQ', 'AGM', 'PR3', 'DOH', 'CCS', 'GSC', 'GHG', 'OAS', 'MIS', 'SIN', 'TPL', 'SAC', '4HT', 'FGP', 'HSO', 'LYZ', 'FGL', 'PRS', 'DCY', 'LYM', 'GPL', 'PYX', 'PCC', 'EHP', 'CHG', 'TPO', 'DAS', 'AYA', 'TYN', 'SVA', 'SCY', 'BNN', '5HP', 'HAR', 'IAS', 'SNC', 'AHB', 'PTR', 'PHI', 'NPH', 'PHL', 'SNN', 'A66', 'TYB', 'PHD', 'MAA', 'APN', 'TYY', 'TYT', 'TIH', 'TRG', 'CXM', 'DIV', 'TYS', 'DTH', 'MLE', 'CME', 'SHR', 'OCY', 'DTY', '2AS', 'AEI', 'DTR', 'OCS', 'CMT', 'BET', 'NLP', 'LLY', 'SCH', 'CEA', 'LLP', 'TRF', 'HMR', 'TYI', 'TRO', 'NLE', 'BMT', 'BUC', 'PEC', 'BUG', 'SCS', 'NLN', 'MHO', 'CSO', 'FTR', 'DLE', 'TRN', 'CSE', 'CSD', 'OMT', 'CSA', 'DSP', 'CSB', 'DSN', 'SHC', 'CSX', 'YCM', 'CSZ', 'TRQ', 'CSW', 'EFC', 'CSP', 'CSS', 'CSR', 'CZZ', 'MSO', 'BTR', 'HLU', 'MGN', 'HTI', 'TYQ', '4IN', 'M3L', 'C5C', 'HTR', 'MPQ', 'KCX', 'GLH', 'DIL', 'ACA', 'NEM', '5CS', 'LYX', 'DVA', 'ACL', 'GLX', 'MLZ', 'GLZ', 'SME', 'SMC', 'DLY', 'NEP', 'BCS', 'ASQ', 'SET', 'SEP', 'ASX', 'DGN', 'DGL', 'MHS', 'SEG', 'ASB', 'ASA', 'SEC', 'SEB', 'ASK', 'GGL', 'ASI', 'SEL', 'CGU', 'C6C', 'ASL', 'LTR', 'CLD', 'CLE', 'GMA', '1LU', 'CLB', 'MVA', 'S1H', 'DNP', 'SAR', 'FME', 'ALO', 'ALM', 'LEF', 'MEN', 'TPQ', 'NMC', 'SBD', 'ALY', 'MME', 'GL3', 'ALS', 'SBL', '2MR', 'CAY', '3AH', 'DPR', 'CAS', 'NC1', 'HYP', 'FLA', 'LCX', 'MSE', 'IYR', 'DPN', 'BAL', 'CAF', 'MSA', 'AIB', 'HIP', 'CYQ', 'PCA', 'DAL', 'BFD', 'DAH', 'HIC', 'CYG', 'DAR', 'CYD', 'IIL', 'CYM', 'CYL', 'CY3', 'CY1', 'HAC', '143', 'DHI', 'CY4', 'YOF', 'HPQ', 'SOC', 'DHA', '2LU', 'MLY', 'TRW', 'STY', 'MCL', 'BHD', 'NRQ', 'ARM', 'PRR', 'ARO', "5HP","ABA","PCA","FGL","BHD","HTR","MSE","CEA","ALS","TRO","TPQ","MHO","IAS","HYP","CGU","CSE","RON","3GA","TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HSE","HSD","HSP"]

dico_AA = {
 "ALA": 'A',
 "CYS": 'C',
 "ASP": "D",
 "GLU": 'E',
 "PHE": 'F',
 "GLY": 'G',
 "HIS": 'H',
 "HSE": 'H',
 "HSD": 'H',
 "HSP": 'H',
 "ILE": 'I',
 "LYS": 'K',
 "LEU": 'L',
 "MET": 'M',
 "ASN": 'N',
 "PRO": 'P',
 "GLN": 'Q',
 "ARG": 'R',
 "SER": 'S',
 "THR": 'T',
 "VAL": 'V',
 "TRP": 'W',
 "TYR": 'Y',
 'PAQ': 'Y',
 'AGM': 'R',
 'PR3': 'C',
 'DOH': 'D',
 'CCS': 'C',
 'GSC': 'G',
 'GHG': 'Q',
 'OAS': 'S',
 'MIS': 'S',
 'SIN': 'D',
 'TPL': 'W',
 'SAC': 'S',
 '4HT': 'W',
 'FGP': 'C',
 'HSO': 'H',
 'LYZ': 'K',
 'FGL': 'S',
 'PRS': 'P',
 'DCY': 'C',
 'LYM': 'K',
 'GPL': 'K',
 'PYX': 'C',
 'PCC': 'P',
 'EHP': 'F',
 'CHG': 'A',
 'TPO': 'T',
 'DAS': 'D',
 'AYA': 'A',
 'TYN': 'Y',
 'SVA': 'S',
 'SCY': 'C',
 'BNN': 'A',
 '5HP': 'E',
 'HAR': 'R',
 'IAS': 'D',
 'SNC': 'C',
 'AHB': 'N',
 'PTR': 'Y',
 'PHI': 'F',
 'NPH': 'C',
 'PHL': 'F',
 'SNN': 'D',
 'A66': 'A',
 'TYB': 'Y',
 'PHD': 'D',
 'MAA': 'A',
 'APN': 'A',
 'TYY': 'Y',
 'TYT': 'Y',
 'TIH': 'A',
 'TRG': 'K',
 'CXM': 'M',
 'DIV': 'V',
 'TYS': 'Y',
 'DTH': 'T',
 'MLE': 'L',
 'CME': 'C',
 'SHR': 'K',
 'OCY': 'C',
 'DTY': 'Y',
 '2AS': 'D',
 'AEI': 'T',
 'DTR': 'W',
 'OCS': 'C',
 'CMT': 'C',
 'BET': 'G',
 'NLP': 'L',
 'LLY': 'K',
 'SCH': 'C',
 'CEA': 'C',
 'LLP': 'K',
 'TRF': 'W',
 'HMR': 'R',
 'TYI': 'Y',
 'TRO': 'W',
 'NLE': 'L',
 'BMT': 'T',
 'BUC': 'C',
 'PEC': 'C',
 'BUG': 'L',
 'SCS': 'C',
 'NLN': 'L',
 'MHO': 'M',
 'CSO': 'C',
 'FTR': 'W',
 'DLE': 'L',
 'TRN': 'W',
 'CSE': 'C',
 'CSD': 'A',
 'OMT': 'M',
 'CSA': 'C',
 'DSP': 'D',
 'CSB': 'C',
 'DSN': 'S',
 'SHC': 'C',
 'CSX': 'C',
 'YCM': 'C',
 'CSZ': 'C',
 'TRQ': 'W',
 'CSW': 'C',
 'EFC': 'C',
 'CSP': 'C',
 'CSS': 'C',
 'CSR': 'C',
 'CZZ': 'C',
 'MSO': 'M',
 'BTR': 'W',
 'HLU': 'L',
 'MGN': 'Q',
 'HTI': 'C',
 'TYQ': 'Y',
 '4IN': 'W',
 'M3L': 'K',
 'C5C': 'C',
 'HTR': 'W',
 'MPQ': 'G',
 'KCX': 'K',
 'GLH': 'E',
 'DIL': 'I',
 'ACA': 'A',
 'NEM': 'H',
 '5CS': 'C',
 'LYX': 'K',
 'DVA': 'V',
 'ACL': 'R',
 'GLX': 'Z',
 'MLZ': 'K',
 'GLZ': 'G',
 'SME': 'M',
 'SMC': 'C',
 'DLY': 'K',
 'NEP': 'H',
 'BCS': 'C',
 'ASQ': 'D',
 'SET': 'S',
 'SEP': 'S',
 'ASX': 'B',
 'DGN': 'Q',
 'DGL': 'E',
 'MHS': 'H',
 'SEG': 'A',
 'ASB': 'D',
 'ASA': 'D',
 'SEC': 'C',
 'SEB': 'S',
 'ASK': 'D',
 'GGL': 'E',
 'ASI': 'N',
 'SEL': 'S',
 'CGU': 'E',
 'C6C': 'C',
 'ASL': 'D',
 'LTR': 'W',
 'CLD': 'S',
 'CLE': 'L',
 'GMA': 'E',
 '1LU': 'L',
 'CLB': 'S',
 'MVA': 'V',
 'S1H': 'S',
 'DNP': 'A',
 'SAR': 'G',
 'FME': 'M',
 'ALO': 'T',
 'ALM': 'A',
 'LEF': 'L',
 'MEN': 'N',
 'TPQ': 'Y',
 'NMC': 'G',
 'SBD': 'S',
 'ALY': 'K',
 'MME': 'M',
 'GL3': 'G',
 'ALS': 'C',
 'SBL': 'S',
 '2MR': 'R',
 'CAY': 'C',
 '3AH': 'H',
 'DPR': 'P',
 'CAS': 'C',
 'NC1': 'S',
 'HYP': 'P',
 'FLA': 'A',
 'LCX': 'K',
 'MSE': 'M',
 'IYR': 'Y',
 'DPN': 'F',
 'BAL': 'A',
 'CAF': 'C',
 'MSA': 'G',
 'AIB': 'A',
 'HIP': 'H',
 'CYQ': 'C',
 'PCA': 'E',
 'DAL': 'A',
 'BFD': 'D',
 'DAH': 'F',
 'HIC': 'H',
 'CYG': 'C',
 'DAR': 'R',
 'CYD': 'C',
 'IIL': 'I',
 'CYM': 'C',
 'CYL': 'C',
 'CY3': 'C',
 'CY1': 'C',
 'HAC': 'A',
 '143': 'C',
 'DHI': 'H',
 'CY4': 'C',
 'YOF': 'Y',
 'HPQ': 'F',
 'SOC': 'C',
 'DHA': 'A',
 '2LU': 'L',
 'MLY': 'K',
 'TRW': 'W',
 'STY': 'Y',
 'MCL': 'K',
 'BHD': 'D',
 'NRQ': 'Y',
 'ARM': 'R',
 'PRR': 'A',
 'ARO': 'R'
 }

RNA3 = ["U"]
DNA3 = ["A","T","G","C"]
SOLV = ["HOH","H2O","WAT","DOD"]

# BBATMS = ["N","CA","C","O","OXT"]
BBATMS = ["N","CA","C","O","OXT"]
SCATMS = ["-","N","CA","C","O","OXT"]
NCHIS  = [0,1,2,3,2,0,2,2,4,2,3,2,0,3,5,1,1,1,2,2]

CHIATMS = [ \
	[], \
	[["N","CA","CB","SG"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","ND1"]], \
	[["N","CA","CB","CG1"],["CA","CB","CG1","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","CE"],["CG","CD","CE","NZ"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","SD"], \
	 ["CB","CG","SD","CE"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","OD1"]], \
	[], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","OE1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD"], \
	 ["CB","CG","CD","NE"],["CG","CD","NE","CZ"], \
	 ["CD","NE","CZ","NH1"]], \
	[["N","CA","CB","OG"]], \
	[["N","CA","CB","OG1"]], \
	[["N","CA","CB","CG1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]], \
	[["N","CA","CB","CG"],["CA","CB","CG","CD1"]]]

AASC=[["CB"], \
      ["CB","SG"], \
      ["CB","CG","OD1","OD2"], \
      ["CB","CG","CD","OE1","OE2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ"], \
      [],["CB","CG","ND1","CD2","CE1","NE2"], \
      ["CB","CG1","CG2","CD1"], \
      ["CB","CG","CD","CE","NZ"], \
      ["CD","CG","CD1","CD2"], \
      ["CB","CG","SD","CE"], \
      ["CB","CG","OD1","ND2"], \
      ["CB","CG","CD"], \
      ["CB","CG","CD","OE1","NE2"], \
      ["CB","CG","CD","NE","CZ","NH1","NH2"], \
      ["CB","OG"], \
      ["CB","OG","OG1","CG2"], \
      ["CB","CG1","CG2"], \
      ["CB","CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"], \
      ["CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]]

AABB=["N","CA","C","O"]

#GBINPATH="/data/bin/"
GBINPATH="/data/bin/"
GBINPATH="/home/tintin/tuffery/bin/"
#GHMMPATH="/data/HMM/models/HMM1/"
GHMMPATH="/data/HMM/models/HMM1/"

# http://nmr.cmbi.ru.nl/~jd/Thesis/ChapterFive.html
normHNames = {
	"lIUPAC" : (
		{
			"CYS"	:	["HB2","HB3","HG"],
			"ASP"	:	["HB2","HB3","HD2"],
			"SER"	:	["HB2","HB3","HG"],
			"GLN"	:	["HB2","HB3","HG2","HG3","HE21","HE22"],
			"BCK"	:	["HN","H1","H2","H3","HA"],
			"ILE"	:	["HB","HG12","HG13","HG21","HG22","HG23","HD11","HD12","HD13"],
			"PRO"	:	["H2","H3","HB2","HB3","HG2","HG3","HD2","HD3"],
			"LYS"	:	["HB2","HB3","HG2","HG3","HD2","HD3","HE2","HE3","HZ1","HZ2","HZ3"],
			"THR"	:	["HB","HG1","HG21","HG22","HG23"],
			"PHE"	:	["HB2","HB3","HD1","HD2","HE1","HE2","HZ"],
			"ALA"	:	["HB1","HB2","HB3"],
			"GLY"	:	["HA2","HA3"],
			"HIS"	:	["HB2","HB3","HD1","HD2","HE1","HE2"],
			"GLU"	:	["HB2","HB3","HG2","HG3","HE2"],
			"LEU"	:	["HB2","HB3","HG","HD11","HD12","HD13","HD21","HD22","HD23"],
			"ARG"	:	["HB2","HB3","HG2","HG3","HD2","HD3","HE","HH11","HH12","HH21","HH22"],
			"TRP"	:	["HB2","HB3","HD1","HE1","HE3","HZ3","HH2","HZ2"],
			"VAL"	:	["HB","HG11","HG12","HG13","HG21","HG22","HG23"],
			"ASN"	:	["HB2","HB3","HD21","HD22"],
			"TYR"	:	["HB2","HB3","HD1","HD2","HE1","HE2","HH"],
			"MET"	:	["HB2","HB3","HG2","HG3","HE1","HE2","HE3"]
		}
		,
		{
			"BCK"	:	[("HN",0,0),("H",1,0),("H",2,0),("H",3,0),("HA",0,0)],
			"ALA"	:	[("HB",1,0),("HB",2,0),("HB",3,0)],
			"ARG"	:	[("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HD",2,0),("HD",3,0),("HE",0,0),("HH",1,1),("HH",1,2),("HH",2,1),("HH",2,2)],
			"ASN"	:	[("HB",2,0),("HB",3,0),("HD",2,1),("HD",2,2)],
			"ASP"	:	[("HB",2,0),("HB",3,0),("HD",2,0)],
			"CYS"	:	[("HB",2,0),("HB",3,0),("HG",0,0)],
			"GLN"	:	[("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HE",2,1),("HE",2,2)],
			"GLU"	:	[("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HE",2,0)],
			"GLY"	:	[("HA",2,0),("HA",3,0)],
			"HIS"	:	[("HB",2,0),("HB",3,0),("HD",1,0),("HD",2,0),("HE",1,0),("HE",2,0)],
			"ILE"	:	[("HB",0,0),("HG",1,2),("HG",1,3),("HG",2,1),("HG",2,2),("HG",2,3),("HD",1,1),("HD",1,2),("HD",1,3)],
			"LEU"	:	[("HB",2,0),("HB",3,0),("HG",0,0),("HD",1,1),("HD",1,2),("HD",1,3),("HD",2,1),("HD",2,2),("HD",2,3)],
			"LYS"	:	[("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HD",2,0),("HD",3,0),("HE",2,0),("HE",3,0),("HZ",1,0),("HZ",2,0),("HZ",3,0)],
			"MET"	:	[("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HE",1,0),("HE",2,0),("HE",3,0)],
			"PHE"	:	[("HB",2,0),("HB",3,0),("HD",1,0),("HD",2,0),("HE",1,0),("HE",2,0),("HZ",0,0)],
			"PRO"	:	[("H",2,0),("H",3,0),("HB",2,0),("HB",3,0),("HG",2,0),("HG",3,0),("HD",2,0),("HD",3,0)],
			"SER"	:	[("HB",2,0),("HB",3,0),("HG",0,0)],
			"THR"	:	[("HB",0,0),("HG",1,0),("HG",2,1),("HG",2,2),("HG",2,3)],
			"TRP"	:	[("HB",2,0),("HB",3,0),("HD",1,0),("HE",1,0),("HE",3,0),("HZ",3,0),("HH",2,0),("HZ",2,0),],
			"TYR"	:	[("HB",2,0),("HB",3,0),("HD",1,0),("HD",2,0),("HE",1,0),("HE",2,0),("HH",0,0)],
			"VAL"	:	[("HB",0,0),("HG",1,1),("HG",1,2),("HG",1,3),("HG",2,1),("HG",2,2),("HG",2,3)]
		}
	),
	"lPDB" : (
		{
			"CYS"	:	["1HB","2HB","HG"],
			"ASP"	:	["1HB","2HB","HD2"],
			"SER"	:	["1HB","2HB","HG"],
			"GLN"	:	["1HB","2HB","1HG","2HG","2HE2","1HE2"],
			"BCK"	:	["H","1H","2H","3H","HA"],
			"ILE"	:	["HB","1HG1","2HG1","1HG2","2HG2","3HG2","1HD1","2HD1","3HD1"],
			"PRO"	:	["H2","H1","1HB","2HB","1HG","2HG","1HD","2HD"],
			"LYS"	:	["1HB","2HB","1HG","2HG","1HD","2HD","1HE","2HE","1HZ","2HZ","3HZ"],
			"THR"	:	["HB","HG1","1HG2","2HG2","3HG2"],
			"PHE"	:	["1HB","2HB","HD1","HD2","HE1","HE2","HZ"],
			"ALA"	:	["1HB","2HB","3HB"],
			"GLY"	:	["1HA","2HA"],
			"HIS"	:	["1HB","2HB","HD1","HD2","HE1","HE2"],
			"GLU"	:	["1HB","2HB","1HG","2HG","HE2"],
			"LEU"	:	["1HB","2HB","HG","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2"],
			"ARG"	:	["1HB","2HB","1HG","2HG","1HD","2HD","HE","1HH1","2HH1","1HH2","2HH2"],
			"TRP"	:	["1HB","2HB","HD1","HE1","HE3","HZ3","HH2","HZ2"],
			"VAL"	:	["HB","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2"],
			"ASN"	:	["1HB","2HB","2HD2","1HD2"],
			"TYR"	:	["1HB","2HB","HD1","HD2","HE1","HE2","HH"],
			"MET"	:	["1HB","2HB","1HG","2HG","1HE","2HE","3HE"]
		}
		,
		{
			"BCK"	:	[("H",0,0),("H",1,0),("H",2,0),("H",3,0),("HA",0,0)],
			"ALA"	:	[("HB",1,0),("HB",2,0),("HB",3,0)],
			"ARG"	:	[("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HD",1,0),("HD",2,0),("HE",0,0),("HH",1,1),("HH",2,1),("HH",1,2),("HH",2,2)],
			"ASN"	:	[("HB",1,0),("HB",2,0),("HD",2,2),("HD",1,2)],
			"ASP"	:	[("HB",1,0),("HB",2,0),("HD",0,2)],
			"CYS"	:	[("HB",1,0),("HB",2,0),("HG",0,0)],
			"GLN"	:	[("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HE",2,2),("HE",1,2)],
			"GLU"	:	[("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HE",0,2)],
			"GLY"	:	[("HA",1,0),("HA",2,0)],
			"HIS"	:	[("HB",1,0),("HB",2,0),("HD",0,1),("HD",0,2),("HE",0,1),("HE",0,2)],
			"ILE"	:	[("HB",0,0),("HG",1,1),("HG",2,1),("HG",1,2),("HG",2,2),("HG",3,2),("HD",1,1),("HD",2,1),("HD",3,1)],
			"LEU"	:	[("HB",1,0),("HB",2,0),("HG",0,0),("HD",1,1),("HD",2,1),("HD",3,1),("HD",1,2),("HD",2,2),("HD",3,2)],
			"LYS"	:	[("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HD",1,0),("HD",2,0),("HE",1,0),("HE",2,0),("HZ",1,0),("HZ",2,0),("HZ",3,0)],
			"MET"	:	[("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HE",1,0),("HE",2,0),("HE",3,0)],
			"PHE"	:	[("HB",1,0),("HB",2,0),("HD",0,1),("HD",0,2),("HE",0,1),("HE",0,2),("HZ",0,0)],
			"PRO"	:	[("H",0,2),("H",0,1),("HB",1,0),("HB",2,0),("HG",1,0),("HG",2,0),("HD",1,0),("HD",2,0)],
			"SER"	:	[("HB",1,0),("HB",2,0),("HG",0,0)],
			"THR"	:	[("HB",0,0),("HG",0,1),("HG",1,2),("HG",2,2),("HG",3,2)],
			"TRP"	:	[("HB",1,0),("HB",2,0),("HD",0,1),("HE",0,1),("HE",0,3),("HZ",0,3),("HH",0,2),("HZ",0,2)],
			"TYR"	:	[("HB",1,0),("HB",2,0),("HD",0,1),("HD",0,2),("HE",0,1),("HE",0,2),("HH",0,0)],
			"VAL"	:	[("HB",0,0),("HG",1,1),("HG",2,1),("HG",3,1),("HG",1,2),("HG",2,2),("HG",3,2)]
		}
	)
}

"""
# OLD IUPACHNames
HNames = {
	"ALA" : ["HB1","HB2","HB3"],
	"CYS" : ["HB1","HB2","HG"],
	"ASP" : ["HB1","HB2"],
	"GLU" : ["HB1","HB2","HG1","HG2"],
	"PHE" : ["HB1","HB2","HD1","HE1","HZ","HE2","HD2"],
	"GLY" : ["HA2"],
	"HIS" : ["HB1","HB2","HD2","HE1","HD1"],
	"ILE" : ["HB","HG11","HG12","HD11","HD12","HD13","HG21","HG22","HG23"],
	"LYS" : ["HB1","HB2","HG1","HG2","HD1","HD2","HE1","HE2","HZ1","HZ2","HZ3"],
	"LEU" : ["HB1","HB2","HG","HD11","HD12","HD13","HD21","HD22","HD23"],
	"MET" : ["HB1","HB2","HG1","HG2","HE1","HE2","HE3"],
	"ASN" : ["HB1","HB2","HD21","HD22"],
	"PRO" : ["HB1","HB2","HG1","HG2","HD1","HD2"],
	"GLN" : ["HB1","HB2","HG1","HG2","HE21","HE22"],
	"ARG" : ["NH1","NH2","HB1","HB2","HG1","HG2","HD1","HD2","HE","HH11","HH12","HH21","HH22"],
	"SER" : ["HB1","HB2","HG"],
	"THR" : ["HB","HG1","HG21","HG22","HG23"],
	"VAL" : ["HB","HG11","HG12","HG13","HG21","HG22","HG23"],
	"TRP" : ["HB1","HB2","HD1","HE1","HZ2","HH2","HZ3","HE3"],
	"TYR" : ["HB1","HB2","HD1","HE1","HE2","HD2","HH"],
	"BCK" : ["HA","HN","HN1","HN2","HN3"]
	}

# OLD PDBHNames
PDBHNames = {
	"ALA" : ["1HB","2HB","3HB"],
	"CYS" : ["1HB","2HB"," HG"],
	"ASP" : ["1HB","2HB"],
	"GLU" : ["1HB","2HB","1HG","2HG"],
	"PHE" : ["1HB","2HB"," HD1"," HE1"," HZ"," HE2"," HD2"],
	"GLY" : ["2HA"],
	"HIS" : ["1HB","2HB"," HD2"," HE1"," HD1"],
	"ILE" : [" HB","1HG1","2HG1","1HD1","2HD1","3HD1","1HG2","2HG2","3HG2"],
	"LYS" : ["1HB","2HB","1HG","2HG","1HD","2HD","1HE","2HE","1HZ","2HZ","3HZ"],
	"LEU" : ["1HB","2HB"," HG","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2"],
	"MET" : ["1HB","2HB","1HG","2HG","1HE","2HE","3HE"],
	"ASN" : ["1HB","2HB","1HD2","2HD2"],
	"PRO" : ["1HB","2HB","1HG","2HG","1HD","2HD"],
	"GLN" : ["1HB","2HB","1HG","2HG","1HE2","2HE2"],
	"ARG" : ["1HB","2HB","1HG","2HG","1HD","2HD"," HE","1HH1","2HH1","1HH2","2HH2"],
	"SER" : ["1HB","2HB"," HG"],
	"THR" : [" HB"," HG1","1HG2","2HG2","3HG2"],
	"VAL" : [" HB","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2"],
	"TRP" : ["1HB","2HB"," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3"],
	"TYR" : ["1HB","2HB"," HD1"," HE1"," HE2"," HD2"," HH"],
	"BCK" : [" HA"," H","1H","2H","3H"]
	}
"""

def resType(aName):
	"""
	@author: P. Tuffery
	@param aName : a string of the sequence encoded using a 3 letters code
	@return: the position of the string
	"""
	if aName == "":
		raise ValueError,"resType expecting a 3 letters code for an amino acid"
	if string.count(AA3, aName) > 0:
		return string.index(AA3, aName)
	else:
		raise IndexError,"%s is an unknown amino acid" % aName

def aa3Type(aName):
	"""
	@author: P. Tuffery
	@param aName : a string of the sequence encoded using a 3 letters code
	@return: the position of the string
	"""
	if aName == "":
		raise ValueError,"aa3Type expecting a 3 letters code for an amino acid"
	if string.count(AA3, aName) > 0 :
		return string.index(AA3, aName)
	else:
		raise IndexError,"%s is an unknown amino acid" % aName

def aa1Type(aName):
	"""
	@author: P. Tuffery
	@param aName : a string of the sequence encoded using a 1 letter code
	@return: the position of the string
	"""
	if aName == "":
		raise ValueError,"aa1Type expecting a 1 letter code for an amino acid"
	if string.count(AA1, aName) > 0:
		return string.index(AA1, aName)
	else:
		raise IndexError,"%s is an unknown amino acid" % aName

# a series of AA3 separated by blanks into aa1 string
def SEQREStoAA1(seqres, verbose = 0):
	"""
	PDB SEQREStoAA1
	@param seqres: a sequence encoded using a 3 letters code (separated by blanks)
	@return: the sequence encoded using a 1 letter code (in a string)
	"""
	seq = ""
	aList = string.split(seqres)
	for aRes in aList:
		if verbose:
			print >> sys.stdout, "# SEQREStoAA1: %s" % aRes
#		if AA3.count(aRes) != 0:
#			if verbose:
#				print >> sys.stdout, "Found as AA3 %d\n" % AA3.index(aRes)
#			seq = seq + AA1seq[AA3.index(aRes)]
		if dico_AA.has_key(aRes):
			if verbose:
				print >> sys.stdout, "#   Found as AA3 %d" % AA3.index(aRes)
			seq = seq + dico_AA[aRes]
		else:
			seq = seq + "X"
	return seq

# This will convert an alignement into a selection mask.
# s1 and s2 must be 2 strings of identical lengths.
# gaps (indels) must be represented by '-'
def aln2mask(s1,s2):
	"""
	This will convert an alignement into a selection mask
	@param s1 : a string
	@param s2 : a string
	@note:
		- s1 and s2 must be 2 strings of identical lengths.\n
		- gaps (indels) must be represented by '-' \n
	"""
	res = ""
	if len(s1) != len(s2):
		return res
	for i in range(0,len(s1)):
		if s1[i] == '-':
			continue
		if s2[i] == '-':
			res = res + '-'
		else:
			res = res + s1[i]
	return res


# any PDB line
class PDBLine:
	"""
	PDBLine : basic management (mostly to access line type (ATOM, REMARK, etc) for one text line of a PDB datafile.
	"""
	def __init__(self, aLine = ""):
		"""
		PDBLine.__init__
		@param aLine: the text sets on the line
		@return: none
		"""
		self.txt = aLine

	def __getslice__(self,ffrom=0,tto=-1):
		"""
		PDBLine.__getslice__(ffrom = 0, tto = None)
		@param ffrom: the first residu considered
		@param tto: the last residu of the residues (excluded, as in python)
		@return: a string, a slice of residues
		"""
		return self.txt[ffrom:tto]

	def __repr__(self):
		"""
		PDBLine.__repr__
		@return: print the PDBLine
		"""
		return str(self.txt)

	def __getitem__(self,aPos):
		"""
		PDB.__getitem__(aPos)
		@param aPos: the position of one residue (PDB[aPos])
		@return: one residue
		"""
		return self.txt[aPos]

	def __len__(self):
		"""
		PDBLine.__len__
		@return: the length of the PDBLine
		"""
		return len(self.txt)

	# back to list of lines
	def flat(self):
		"""
		PDB.flat
		@return: the string of the current line
		"""
		return self.txt

	# header de la ligne
	def header(self):
		"""
		PDBLine header
		@return: the six first chars of the line (string)
		@note: could be one of:
			- HEADER
			- REMARK
			- ATOM
			- CONECT
			- etc
		"""
		#PDBLine header, i.e. its 6 first chars
		#one of:
		#HEADER
		#REMARK
		#ATOM
		#CONECT
		#etc
		try:
			return self.txt[0:6].split()[0]
		except:
			return ""


## ========================================
## a PDB ATOM (HETATM) line
## ========================================
class atmLine(PDBLine):
	"""
	class atmLine
	This models one PDB ATOM / HETATOM line
	and its accessors
	"""
	def __init__(self, aLine = ""):
		"""
		atmLine.__init__
		@param aLine: the text sets on the line
		@type aLine: could be atmLine, PDBLine or a string
		@return: none
		"""
		if isinstance(aLine,atmLine):
			## print "atmLine from atmLine"
			self.txt = aLine.txt
		elif isinstance(aLine,PDBLine):
			## print "atmLine from PDBLine"
			self.txt = aLine.txt
		elif isinstance(aLine,types.StringType):
			## print "atmLine from string"
			self.txt = aLine
		else:
			self.txt = aLine

	def header(self, hdr = ""):
		"""
		PDB line HEADER
		@param hdr: none, or a new header
		@return: ATOM or HETATM if the parameter hdr is not defined, otherwise it sets HETATM as header
		@type hdr: must be string
		"""
		if hdr != "":
			self.txt = "%-.6s%s" % (hdr,self.txt[6:])
			return hdr[:6]
		try:
			return self.txt[0:6].split()[0]
		except:
			return ""

	def atmNum(self, anum = ""):
		"""
		PDB line atom number
		@param anum: none, or an atom number
		@return: the atom number (string) if there is no parameter, otherwise it sets anum as atom number
		@type anum: anum may be int or string
		"""
		if anum != "":
			self.txt = "%s%5d%s" % (self.txt[:6],int(anum),self.txt[11:])
			return anum
		try:
			anum=self.txt[6:11].split()[0]
			return anum
		except ValueError:
			raise ValueError, "Incorrect ATOM atom number format for:\n%s" % self.txt

	def atmName(self, aname = ""):
		"""
		PDB line atom name
		@param aname: none, or an atom name
		@return: the atom name (string) if there is no parameter, otherwise it sets aname as atom number
		@type aname: aname must be string of size 4 at max, must begin by a blank if necessary !
		"""
		if aname != "":
			self.txt = "%s%-4s%s" % (self.txt[:12],aname,self.txt[16:])
			return aname
		try:
			rnum=self.txt[12:16].split()[0]
			return rnum
		except ValueError:
			raise ValueError, "Incorrect ATOM atom name format for:\n%s" % self.txt

	def atmType(self, atype = ""):
		"""
		PDB line atom type
		@param aname: none, or an atom type
		@return: the atom type (string) if there is no parameter, otherwise it sets atype as atom type
		@type atype: atype must be string of size 2 at max, right aligned !
		"""
		if atype != "":
			self.txt = "%s%-2s%s" % (self.txt[:76],atype,self.txt[78:])
			return atype
		try:
			atype=self.txt[76:78].strip()
			return atype
		except ValueError:
			raise ValueError, "Incorrect ATOM atom type format for:\n%s" % self.txt

	def atmBVal(self, abval=""):
		"""
		PDB line atom B value
		@param abfact: none, or an atom B value
		@return: the atom B value (string) if there is no parameter, otherwise it sets abfact as atom B value
		@type anum: abval may be float or string
		"""
		if abval != "":
			self.txt = "%s%6.2f%s" % (self.txt[:60],float(abval),self.txt[66:])
			return abval
		try:
			abval=self.txt[60:66].strip()
			return abval
		except ValueError:
			raise ValueError, "Incorrect ATOM atom B Value format for:\n%s" % self.txt

	def alt(self, acode = ""):
		"""
		PDB line alternate code
		@param acode: none, or an alternate code
		@return: the alternate code (1 char string) if there is no parameter, otherwise it sets acode as alternate code
		@type acode: acode must be one character
		"""
		if acode != "":
			self.txt = "%s%c%s" % (self.txt[:16],acode,self.txt[17:])
			return acode
		try:
			alt=self.txt[16]
			return alt
		except ValueError:
			raise valueError, "Incorrect ATOM alternate code format for:\n%s" % self.txt

	def resName(self, rName = ""):
		"""
		PDB line residue name
		@param rName: none, or a residue name
		@return: the residue name (string) if there is no parameter, otherwise it sets residue name
		@type rName: rName must be string (3 chars)
		"""
		if rName != "":
			self.txt = "%s%-.3s%s" % (self.txt[:17],rName,self.txt[20:])
			return rName[:3]
		try:
			rname=self.txt[17:20].split()[0]
			return rname
		except ValueError:
			raise ValueError, "Incorrect ATOM residue name format for:\n%s" % self.txt

	def chnLbl(self, lbl = ""):
		"""
		PDB line chain label
		@param lbl: none, or an atom chain label
		@return: the atom chain label (1 character string) if there is no parameter, otherwise it sets the atom chain label
		@type lbl:lbl must be 1 character
		"""
		if lbl != "":
			self.txt = "%s%c%s" % (self.txt[:21],lbl[0],self.txt[22:])
			return lbl
		try:
			lbl=self.txt[21]
			return lbl
		except ValueError:
			raise ValueError, "Incorrect ATOM chain label format for:\n%s" % self.txt

	def resNum(self, rnum = ""):
		"""
		PDB line residue number
		@param rnum: none, or a residue number rnum
		@return: the residue member (1 character string) if there is no parameter, otherwise it sets the atom residue number
		@type rnum: rnum may be string or int
		"""
		if rnum != "":
			self.txt = "%s%4d%s" % (self.txt[:22],int(rnum),self.txt[26:])
			return rnum
		try:
			rnum=self.txt[22:26].split()[0]
			return rnum
		except ValueError:
			raise ValueError, "Incorrect ATOM line format for:\n%s" % self.txt

	def resType(self,verbose = 0):
		"""
		PDB liste resType
		@return: the type of the residues of the list
		@note: it could be:
			- AMINO-ACID
			- RNA
			- DNA
			- SOLVENT
			- HETERO
		"""
		aName = self.resName()
		if AA3.count(aName) > 0:
			return "AMINO-ACID"
		elif RNA3.count(aName) > 0:
			return "RNA"
		elif DNA3.count(aName) > 0:
			return "DNA"
		elif SOLV.count(aName) > 0:
			return "SOLVENT"
		else:
			return "HETERO"

	def icode(self, thecode = ""):
		"""
		PDB line code number
		@param thecode: none, or a residue code
		@return: the residue code (1 character string) if there is no parameter, otherwise it sets the residue code
		@type thecode: the code must be 1 character
		"""
		if thecode != "":
			self.txt = "%s%c%s" % (self.txt[:26],thecode,self.txt[27:])
			return thecode
		try:
			icode=self.txt[26]
			return icode
		except ValueError:
			raise ValueError, "Incorrect ATOM line format for:\n%s" % self.txt

	def xyz(self):
		"""
		PDB line xyz
		@return: the atom's coordinated in 3D (x,y,z)
		@note: x, y, z returned are float
		"""
		try:
			x=string.split(self.txt[30:38])[0]
			y=string.split(self.txt[38:46])[0]
			z=string.split(self.txt[46:54])[0]
			return float(x), float(y), float(z)
		except ValueError:
			print >> sys.stderr, "Incorrect ATOM coordinates format for:\n%s" % self.txt
			return 0., 0., 0.

	def crds(self):
		"""
		PDB line crds
		@return: a string containing the atom's coordinated in 3D (x,y,z),
		or 0.,0.,0. if the crds are not found
		@note: x, y and z are contained in one string
		"""
		return self.txt[30:54]

	def setcrds(self,x,y,z):
		"""
		PDB line set crds
		@param x,y,z: the atom's coordinated in 3D
		@return: None
		@type x, y, z: must be a long float.
		"""
		self.txt = "%s%8.3lf%8.3lf%8.3lf%s" % (self.txt[:30], x, y, z, self.txt[54:])
		return

	def fpt(self, aQ=""):
		"""
		PDB line occupancy
		@param aQ: none, or a new occupancy for the atom
		@return: the occupancy read on the PDB line or the new occupancy set by the user
		@type aQ: must be a float
		"""
		if aQ != "":
			self.txt = string.replace(self.txt,"\n","") ### /!\
			self.txt = "%-60s %7.3f\n" % (self.txt[:61],aQ)
			return aQ
		try:
			occ=self.txt[54:61]
			return occ
		except ValueError:
			raise ValueError, "Incorrect ATOM occupancy format for:\n%s" % self.txt

	def q(self, aQ=""):
		"""
		PDB line occupancy
		@param aQ: none, or a new occupancy for the atom
		@return: the occupancy read on the PDB line or the new occupancy set by the user
		@type aQ: must be a float
		"""
		if aQ != "":
			self.txt = "%s%7.3f%s" % (self.txt[:54],aQ,self.txt[61:])
			return aQ
		try:
			occ=self.txt[54:61]
			return occ
		except ValueError:
			raise ValueError, "Incorrect ATOM occupancy format for:\n%s" % self.txt

	def occ(self, aOcc=""):
		"""
		PDB line occupancy
		@param aOcc: none, or a new occupancy for the atom
		@return: the occupancy read on the PDB line or the new occupancy set by the user
		@type aOcc: must be a float
		"""
		if aOcc != "":
			self.txt = "%s%6.2f%s" % (self.txt[:54],aOcc,self.txt[60:])
			return aOcc
		try:
			occ=self.txt[54:60]
			return occ
		except ValueError:
			raise ValueError, "Incorrect ATOM occupancy format for:\n%s" % self.txt

	def r(self, aR=""):
		"""
		PDB line B_iso_or_equiv
		@param aR: none, or a new B_iso_or_equiv for the atom
		@return: the B_iso_or_equiv read on the PDB line or the new B_iso_or_equiv set by the user
		@type aR: must be a float
		"""
		if aR != "":
			self.txt = "%-61s%7.3f%s" % (self.txt[:61],aR,self.txt[68:])
			return aR
		try:
			occ=self.txt[61:68]
			return occ
		except ValueError:
			raise ValueError, "Incorrect ATOM B iso format for:\n%s" % self.txt

	def tfac(self, tFac = ""):
		"""
		PBD line tfac
		@param tFac: none, or a new temperature factor for the atom
		@return: the temperature factor read on the PDB line or the new temperature factor set by the user
		@type tFac: must be a float
		"""
		if tFac != "":
			self.txt = "%-60s%6.2f%s\n" % (self.txt[:-1].ljust(60)[:60],float(tFac),self.txt[66:-1])
			return tFac
		try:
			tfac=self.txt[60:66]
			return tfac
		except ValueError:
			raise ValueError, "Incorrect ATOM temperature factor format for:\n%s" % self.txt

	def segId(self):
		"""
		PBD line segId
		@return: the segment Id
		"""
		try:
			segId=self.txt[72:76]
			return segId
		except ValueError:
			raise ValueError, "Incorrect ATOM segment ID format for:\n%s" % self.txt

	def ele(self):
		"""
		PDB line ele
		@return: the element symbol
		"""
		try:
			ele=self.txt[76:78]
			return ele
		except ValueError:
			raise ValueError, "Incorrect ATOM atom element symbol format for:\n%s" % self.txt

	def chrg(self):
		"""
		PDB line type symbol
		@return: the type symbol
		"""
		try:
			chrg=self.txt[78:80]
			return chrg
		except ValueError:
			raise ValueError, "Incorrect ATOM charge format for:\n%s" % self.txt


## ========================================
## A series of PDB ATOM (HETATM) lines
## Considered as a set of residues
## Tabulation of residues is achieved
##
## atmList always return atmList,
## EXCEPT for __getitem__ when requesting in 1 residue
## where it is desirable to return atmLine
##
## atom lines accessible as: x.atms
##
## With this class, we are simply manipulating text
## No semantics associated
## ========================================
class atmList(atmLine):
	"""
	class atmList
	This models a list of PDB ATOM / HETATM lines (atmLine)
	"""
	def __init__(self, data = "", chId = "", hetSkip = 0, verbose = 0):
		"""
		atmList.__init__ determine the type of data to initialize
		@param data: an instance
		@return: none
		"""
		# Order of parsing is important (inheritance)
		#
		# from PDB: just retain PDB.data field
		if isinstance(data,PDB):
			if verbose > 1:
				print >> sys.stdout, "# atmList from PDB"
			self.list = data.data
		# from residue: just retain residue.list field
		elif isinstance(data,residue):
			if verbose > 1:
				print >> sys.stdout, "# atmList from residue"
			self.list = data.atms
		# from atmList: just propagate
		elif isinstance(data,atmList):
			if verbose > 1:
				print >> sys.stdout, "# atmList from atmList"
			self.list = data.list
		# from atmLine: just wrap
		elif isinstance(data,atmLine):
			## We force one line as a residue
			if verbose > 1:
				print >> sys.stdout, "# atmList  from atmLine"
			self.list = []
			self.list.append(data)
		# from list: suppose a list of atomic lines
		elif isinstance(data,types.ListType):
			if verbose > 1:
				print >> sys.stdout, "# atmList from ListType"
			self.list = []
			for aLine in data:
				self.list.append(atmLine(aLine))
			## self.resTab(verbose)
		else:
			if verbose > 1:
				print >> sys.stdout, "# atmList from unknown"
			self.list  = []

	def __len__(self):
		"""
		atmList.__len__
		@return: the length of the atmList
		"""
		return len(self.list)

	def __add__(self,new):
		"""
		atmList.__add__
		@param new: the list added to the first one
		@return: the two lists concatenated
		"""
		print >> sys.stdout, "# __add__.atmList"
		return atmList(self.list[:] + new.list[:])

	def __getslice__(self,ffrom,tto):
		"""
		atmList.__getslice__(ffrom = 0, tto = None)
		@param ffrom: the first residu considered
		@param tto: the last residu of the residues (excluded, as in python)
		@return: a sub atmList, a slice of residues
		"""
		return atmList(self.list[ffrom:tto])

	def __getitem__(self,aPos):
		"""
		atmList.__getitem__(aPos)
		@param aPos: the number of one atom in an atmList (PDB[residue][aPos])
		@return: one atom line
		"""
		return self.list[aPos]

	def __delitem__(self,aPos):
		"""
		atmList.__delitem__
		@param aPos: the position of the atom to delete in the atmlist
		@return: none
		"""
		## delete old series
		#aDex = self.rt[aPos][0]
		###print "Removing atoms ",self.rt[aPos][0]," to ",self.rt[aPos+1][0]
		#for aAtm in range(self.rt[aPos][0],self.rt[aPos+1][0]):
		#	del self.atms[aDex]
		#self.resTab(0)
		self.list[aPos:aPos+1] = []

	def __setitem__(self,aPos, new):
		"""
		atmList.__setitem__
		@param aPos: the position of the atom to insert
		@param new: the new atom to insert in the position aPos
		@return: none
		"""
		# delete old
		del self[aPos]
		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.atms.insert(aDex,aAtm)		# /!\ aDex is not an index, method can't work
			aDex = aDex + 1
		self.resTab(0)

	def __repr__(self, altLbl = "", OXTSkip = 0, HSkip = 0):
		"""
		atmList.__repr__
		@param OXTSkip: does it skip OXT? (Yes:1 ; No:0)
		@param HSkip: does it skip hydrogen? (Yes:1 ; No:0)
		@param altLbl: alternate atom label
		@return: show atomic information of the atmList in a string
		"""
		res = ""
		for aAtm in self.list:
			if altLbl != "":
				alt = aAtm.alt()
				if alt != " " and alt != altLbl:
					continue
			if OXTSkip != 0:
				if aAtm.atmName() == "OXT":
					continue
			if HSkip != 0:
				atmName = aAtm.atmName()
				if atmName[0] == "H" or (atmName[0] in  "1234" and atmName[1] == "H"):
					continue
			res = res + str(aAtm)
		return res

	def flat(self, altLbl = "", OXTSkip = 0, PDBMac = 0, keepH = 1):
		"""
		PDB list flat
		@param altCare: does it take care about alternate atom? (Yes:1 ; No:0)
		@param PDBMac:
		@param altLbl: alternate atom label
		@param OXTSkip: does it skip OXT? (Yes:1 ; No:0)
		@param keepH: does it keep the hydrogen atoms?
		@return: a list of lines
		"""
		res = []
		for aAtm in self.list:
			if altLbl != "":
				alt = aAtm.alt()
				if alt != " " and alt != altLbl:
					continue
			if PDBMac and aAtm.atmName() == "O1":
				aAtm.atmName(" O")
			if PDBMac and aAtm.atmName() == "O2":
				aAtm.atmName(" OXT")
			if PDBMac and aAtm.atmName() == "OT1":
				aAtm.atmName(" O")
			if PDBMac and aAtm.atmName() == "OT2":
				aAtm.atmName(" OXT")
			if OXTSkip != 0:
				if aAtm.atmName() == "OXT":
					continue
			if not keepH:
				atmName = aAtm.atmName()
				if atmName[0] == "H":
					continue
				if atmName[0] in "1234":
					if atmName[1] == "H":
						continue
					if atmName[1] in "1234":
						if atmName[2] == "H":
							continue
			res.append(aAtm.flat())
		return res

	def insert(self,aPos, new):
		"""
		PDB list insert
		@param aPos: the position of new
		@param new: the list inserted
		@return: the list after insertion of new
		"""
		# insert new
		aDex = self.rt[aPos][0]
		for aAtm in new.atms:
			self.list.insert(aDex,aAtm)
			aDex = aDex + 1
		self.resTab(0)

	def crds(self, ffrom = 0, tto = -1): # /!\ same as xyz() ?
		"""
		PDB List coordinates (crds)
		@param ffrom: the first atom of the range that we need coordinates
		@param tto: the last one, -1 by default, in this case it will be equal to the last atom of the list
		@return: a list of all coordinates of the range (concatenated)
		"""
		if tto == -1:
			tto = len(self.list)
		res = []
		for aAtm in range(ffrom, tto):
			res.append(atmLine(self.list[aAtm]).crds())
		return res

	def xyz(self, ffrom = 0, tto = -1):
		"""
		PDB List xyz
		@param ffrom: the first atom of the range that we need coordinates
		@param tto: the last one, -1 by default: in this case it will be equal to the last atom of the list
		@return: a list of lists of coordinates x, y, z of the range
		"""
		if tto == -1:
			tto = len(self.list)
		if (ffrom == 0) and (len(self.list) == 1):
			return atmLine(self.list[ffrom]).xyz()
		else:
			res = []
			for aAtm in range(ffrom, tto):
				res.append(atmLine(self.list[aAtm]).xyz())
			return res

	def BC(self, ffrom = 0, tto = -1):
		"""
		PDB List BC give the center of geometry of a collection of atoms
		@param ffrom: the first atom of the range that we need the center of geometry
		@param tto: the last one, -1 by default: in this case it will be equal to the last atom of the list
		@return: a list of coordinates x, y, z of the center of geometry
		"""
		if tto == -1:
			tto = len(self)
		(x,y,z) = (0.,0.,0.)
		nAtm = 0.
		# for aAtm in self[ffrom:tto].atms:
		for aAtm in self[ffrom:tto]:
			(x1,y1,z1) = aAtm.xyz()
			x = x + x1
			y = y + y1
			z = z + z1
			nAtm = nAtm + 1.
		x = x / nAtm
		y = y / nAtm
		z = z / nAtm
		return (x,y,z)

	def radius(self, ffrom = 0, tto = -1): # /!\ compare to BC ? same ? or "help" error ?
		"""
		PDB List radius
		@param ffrom: the first atom of the range that we need the center of geometry
		@param tto: the last one, -1 by default: in this case it will be equal to the last atom of the list
		@return: a list of coordinates x, y, z of the center of geometry, the maximum one found
		"""
		if tto == -1:
			tto = len(self)
		(x,y,z) = self.BC(ffrom, tto)
		rs = 0.
		for aAtm in self[ffrom:tto]:
			(x1,y1,z1) = aAtm.xyz()
			r = distance(x,y,z,x1,y1,z1)
			if r > rs:
				rs = r
		return rs

	def gridSize(self, ffrom = 0, tto = -1): # /!\ compare to BC ? same ? or "help" error ?
		"""
		PDB List gridSize
		@param ffrom: the first atom of the range that we need the center of geometry
		@param tto: the last one, -1 by default: in this case it will be equal to the last atom of the list
		@return: a list of dimension x, y, z that correspond to the largest radii from BC on each axis.
		"""
		if tto == -1:
			tto = len(self)
		(x,y,z) = self.BC(ffrom, tto)
		rs = 0.
                dx = dy = dz = 0.
		for aAtm in self[ffrom:tto]:
			(x1,y1,z1) = aAtm.xyz()
			ldx, ldy, ldz = dxdydz(x,y,z,x1,y1,z1)
                        if abs(ldx) > dx:
                                dx = abs(ldx)
                        if abs(ldy) > dy:
                                dy = abs(ldy)
                        if abs(ldz) > dz:
                                dz = abs(ldz)
		return dx, dy, dz

	def oneChis(self):
		"""
		PDB List oneChis
		@return: the dihedral of the side chain
		"""
		resTpe = resType(self.list[0].resName())
		res = [AA3[resTpe]]
		for aChi in CHIATMS[resTpe]:
			aAtm = self.theAtm(aChi[0])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[1])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[2])
			if aAtm == []:
				return res
			aAtm = self.theAtm(aChi[3])
			if aAtm == []:
				return res
			a = self.theAtm(aChi[0]).xyz()
			b = self.theAtm(aChi[1]).xyz()
			c = self.theAtm(aChi[2]).xyz()
			d = self.theAtm(aChi[3]).xyz()
			res.append(apply(dihedral,a+b+c+d))
		return res

	def chis(self):
		"""
		PDB List chis
		@return: the dihedral of the side chain of the residues
		"""
		res = []
		if len(self) == 1:
			res.append(self.oneChis())
			return res
		for aRes in range(0,len(self)):
			res.append(self[aRes].oneChis()) ### /!\ self[aRes] == atmLine ... so ... no oneChis method ...
		return res

	def outChis(self):
		"""
		PDB List outChis
		@return: none, print the dihedral of the side chain of the residues
		"""
		chis = self.chis()
		for i in chis:
			print >> sys.stdout, i[0],
			for j in i[1:]:
				print >> sys.stdout, j,
			print >> sys.stdout

	def atmPos(self, aName):
		"""
		PDB List atmPos
		@param aName: name of the atom searched
		@return: aPos the position of the first atom named "aName"
		"""
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == aName:
				return aPos
		return None

	def Npos(self):
		"""
		PDB List Npos
		@return: aPos the position of the first atom N
		"""
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "N":
				return aPos
		return None

	def CApos(self):
		"""
		PDB List CApos
		@return: aPos the position of the firstatom CA
		"""
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "CA":
				return aPos
		return None

	def Cpos(self):
		"""
		PDB List Cpos
		@return: aPos the position of the first atom C
		"""
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "C":
				return aPos
		return None

	def Opos(self):
		"""
		PDB List Opos
		@return: aPos the position of the first atom O
		"""
		for aPos in range(0,len(self.list)):
			if self.list[aPos].atmName() == "O":
				return aPos
		return None

	def out(self): # /!\ ???
		"""
		PDB List out
		@return: none
		@note: just "pass"
		"""
		pass

	def resName(self):
		"""
		PDB list residue name
		@return: the residue name (string) if there is no parameter, otherwise it sets residue name
		"""
##		print self.list[0]
##		print self.list[0].__class__
##		print atmLine(self.list[0])
##		print "tutu"
##		print self.__class__, "resName",len(self.list)
##		print self.list[0].__class__, "resName",len(self.list[0])
##		return self.list[0].resName()
		return atmLine(self.list[0]).resName()

	def theAtm(self,atmName = ""):
		"""
		PDB list theAtm
		@param atmName: the name of the atom searched to know its line
		@return: the line of atmName
		"""
		for aLine in self.list:
			if atmLine(aLine).atmName() == atmName:
				return atmLine(aLine)
		return []

	def isPDB(self):
		"""
		PDB list isPDB:
		@return: 1
		"""
		return 1

	def write(self, outName = "", label="", hetSkip = 0,verbose = 0):
		"""
		PDB List write PDB or PDB chain(s) to file
		@param outName: the name of the file written, if no name sets, it will be write on the standard out
		@param label:
		@return: none
		@note: I{for example:}\n
			from PDB import * \n
			x = protein("/home/raid5/PDB/pdb1acc.ent.gz",hetSkip=1) \n
			x.frg(0).write()\n
		"""
		if outName == "":
			f = sys.stdout
		else:
			f = open(outName,"w")
		if verbose > 1:
			print >> sys.stderr, "# Writing in %s" % f.name
		f.write("HEADER %s (%d residues)\n" % (label, len(self)))
		for aAtm in self.list:
			f.write("%s" % aAtm)

	def oneHMMGeo(self, aCA):
		"""
		PDB List oneHMMGeo
		@param aCA: an atom
		@return: the seven descriptors of a fragment of one letter of the structural alphabet
		"""
		CA1x, CA1y, CA1z = self[aCA].xyz()
		CA2x, CA2y, CA2z = self[aCA+1].xyz()
		CA3x, CA3y, CA3z = self[aCA+2].xyz()
		CA4x, CA4y, CA4z = self[aCA+3].xyz()
		d1 = distance(CA1x, CA1y, CA1z, CA3x, CA3y, CA3z)
		d2 = distance(CA1x, CA1y, CA1z, CA4x, CA4y, CA4z)
		d3 = distance(CA2x, CA2y, CA2z, CA4x, CA4y, CA4z)
		x1, y1, z1 = vecteur(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		x2, y2, z2 = vecteur(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		x3, y3, z3 = vecteur(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		d4 = mixtproduct(x1, y1, z1, x2, y2, z2, x3, y3, z3)
		d5 = distance(CA1x, CA1y, CA1z, CA2x, CA2y, CA2z)
		d6 = distance(CA2x, CA2y, CA2z, CA3x, CA3y, CA3z)
		d7 = distance(CA3x, CA3y, CA3z, CA4x, CA4y, CA4z)
		return d1,d2,d3,d4,d5,d6,d7


## ========================================
## The ONE residue class
## ========================================
class residue(atmList):
	"""
	class residue
	This models a list of PDB ATOM / HETATOM lines with semantic significance
	"""
	def __init__(self,data="",verbose=0):
		"""
		residue.__init__ determine the type of data to initialize the residue
		@param data: an instance
		@return: none
		"""
		if data == "":
			self.atms = []
			self.type = None
			self.name = None
		else:
			if isinstance(data,residue): # residue instance
				## print "residue from residue"
				self.atms = data.atms
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data,atmList): # atmList instance
				## print "residue from atmList"
				self.atms = data
				self.type = self.rType()
				self.name = self.rName()
			elif isinstance(data,atmLine): # atmLine instance
				## print "residue from atmLine"
				self.atms = atmList(data)
				## self.type = self.rType()
				## self.name = self.rName()
			else:
				## print "residue from unknown",data.__class__
				self.atms = atmList(data)
				self.type = self.rType()
				self.name = self.rName()

	def __len__(self):
		"""
		residue.__len__
		@return: the length of the atmList
		"""
		return len(self.atms)

	def __repr__(self, altCare = 0, altLbl = "", OXTCare = 0, HSkip = 0):
		"""
		residue.__repr__
		@param altCare: does it take care about alternate atom? (Yes:1 ; No:0)
		@param OXTCare: does it take care about OXT? (Yes:1 ; No:0)
		@param HSkip: does it skip hydrogen? (Yes:1 ; No:0)
		@param altLbl: alternate atom label
		@return: show atomic information of the atmList in a string
		"""
		OXTSkip = 0
		if OXTCare != 0:
			if self.atms.atmPos("O") != None and self.atms.atmPos("OXT") != None :
				OXTSkip = 1
		if altCare != 0:
			altLbls = self.altLbls()
			if altLbl == "":
				if altLbls != "":
					altLbl = altLbls[0]
			else:
				if string.count(altLbls, altLbl):
					pass
				else:
					altLbl = altLbls[0]
		return self.atms.__repr__(altLbl, OXTSkip = OXTSkip, HSkip = HSkip)

	def __getslice__(self,ffrom,tto):
		"""
		residue.__getslice__(ffrom = 0, tto = None)
		@param ffrom: the first residu considered
		@param tto: the last residu of the residues (excluded, as in python)
		@return: a sub atmList, a slice of atoms
		"""
		if len(self.atms) == 1:
			return residue(self.atms)
		if tto > len(self.atms):
			tto = len(self.atms)
		return self.atms[ffrom:tto]

	def __getitem__(self,aPos):
		"""
		residue.__getitem__(aPos)
		@param aPos: the position of one residue (PDB[aPos]) or CHAINS (PDB["CHAINS"])
		@return: one residue or PDB instance of chains matching CHAINS (e.g. \"AB\")
		"""
		if isinstance(aPos,types.IntType):
			if aPos > len(self.atms):
				return None
			elif aPos < 0:
				if aPos + len(self.atms) < 0:
					return None
				else:
					return self.atms[aPos]
			## else, we return atmList
			return self.atms[aPos]
		elif isinstance(aPos,types.StringType):
			for iAtm in range(0,len(self.atms)):
				if self.atms[iAtm].atmName() == aPos:
					return self.atms[iAtm]

	def flat(self, altCare = 0, altLbl = "", OXTCare = 0, PDBMac = 0, keepH = 1):
		"""
		residue.flat
		@param altCare: does it take care about alternate atom? (Yes:1 ; No:0)
		@param altLbl: alternate atom label
		@param PDBMac:
		@param OXTCare: does it take care about OXT? (Yes:1 ; No:0)
		@return: a string of the list of lines
		"""
		OXTSkip = 0
		if OXTCare != 0:
			if self.atms.atmPos("O") != None and self.atms.atmPos("OXT") != None :
				OXTSkip = 1
		if altCare != 0:
			altLbls = self.altLbls()
			if altLbl == "":
				if altLbls != "":
					altLbl = altLbls[0]
			else:
				if string.count(altLbls, altLbl):
					pass
				else:
					altLbl = altLbls[0]
		return self.atms.flat(altLbl, OXTSkip = OXTSkip, PDBMac = PDBMac, keepH = keepH)

	def rName(self, name = "", verbose = 0):
		"""
		residue.rName give residue name
		@param name: none, or a residue name
		@return: the residue name (string) if there is no parameter, otherwise it sets residues name
		@type name: must be string (3 chars)
		"""
		if name == "":
			return self.atms[0].resName()
		else:
			for atm in self.atms:
				atm.resName(name)

	def rNum(self,aNum = "", verbose = 0):
		"""
		residue.rNum give residue number
		@param aNum: none, or a residue number aNum
		@return: the residue number (1 character string) if there is no parameter, otherwise it sets the atoms residue number
		@type aNum: aNum may be string or int
		"""
		if aNum == "":
			return self.atms[0].resNum()
		else:
			for atm in self.atms:
				atm.resNum(aNum)

	def tFac(self, tFac = "", verbose = 0):
		if tFac == "":
			return self.atms[0].tfac()
		else:
			for atm in self.atms:
				atm.tfac(tFac)

	def riCode(self,icode = "",verbose = 0):
		"""
		residue.riCode give code number
		@param icode: none, or a residue code
		@return: the residue code (1 character string) if there is no parameter, otherwise it sets the residues code
		@type icode: the code must be 1 character
		"""
		if icode == "":
			return self.atms[0].icode()
		else:
			for atm in self.atms:
				atm.icode(icode)

	def rType(self,verbose = 0):
		"""
		residue.rType give the type of residue
		@return: the type of the residues of the list
		@note: it could be:
			- AMINO-ACID
			- RNA
			- DNA
			- SOLVENT
			- HETERO
		"""
		aName = self.atms[0].resName()
		if AA3.count(aName) > 0:
			return "AMINO-ACID"
		elif RNA3.count(aName) > 0:
			return "RNA"
		elif DNA3.count(aName) > 0:
			return "DNA"
		elif SOLV.count(aName) > 0:
			return "SOLVENT"
		else:
			return "HETERO"

	def chnLbl(self,lbl = "", verbose = 0):
		"""
		residue.chnLbl give or set chain label
		@param lbl: none, or an atom chain label
		@return: the atom chain label (1 character string) if there is no parameter, otherwise it sets the atoms chain label
		@type lbl:lbl must be 1 character
		"""
		if lbl == "":
			return self.atms[0].chnLbl()
		else:
			for atm in self.atms:
				atm.chnLbl(lbl)

	def atmPos(self, aName):
		"""
		residue.atmPos give atmPos
		@param aName: name of the atom searched
		@return: aPos the position of the first atom named "aName"
		"""
		return self.atms.atmPos(aName)

	def hasAltAtms(self,verbose = 0):
		"""
		residue.hasAltAtms
		@return: Does the file has BBaltAtm or SCAltAtm? (Yes/No for each)
		"""
		BBAltAtm = False
		SCAltAtm = False
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
			alt = aAtm.alt()
			if alt != ' ':
				isAlt = 1
				if string.count(string.digits,aAtm.txt[12]):
					isAlt = 0
				if aAtm.txt[12] == ' ' and aAtm.txt[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = aAtm.atmName()
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "O":
					BBAltAtm = True
				else:
					SCAltAtm = True
		return BBAltAtm, SCAltAtm

	def altLbls(self,verbose = 0):
		"""
		residue.atlLbls
		@return: all the alternate atoms of the atom list
		"""
		rs = ""
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
			alt = aAtm.alt()
			if alt != ' ':
				if string.count(rs,alt) == 0:
					rs += alt
		return rs

	def select(self,awhat=[""]):
		"""
		residue.select
		@param awhat: what atoms
		@return: a selection of atoms, and atmList
		"""
		res = atmList()
		for iAtm in range(0,len(self.atms)):
			if awhat == [""]:
				res.list.append(atmLine(self.atms[iAtm].txt))
			else:
				if awhat[0] !=  "-":
					if awhat.count(self.atms[iAtm].atmName()) > 0:
						res.list.append(atmLine(self.atms[iAtm].txt))
				else:
					if awhat.count(self.atms[iAtm].atmName()) == 0:
						res.list.append(atmLine(self.atms[iAtm].txt))
		return res

	def delete(self,awhat=None):
		"""
		residue.delete
		@param awhat: none, or a list of atom names
		@return: none
		@note: This will remove atoms from the residue
		based on their names
		"""
		if awhat == None:
			return
		for iAtm in range(len(self.atms)-1, -1, -1):
			if self.atms[iAtm].atmName() in awhat:
				self.atms.list.remove(self.atms.list[iAtm])

	def BBAtmMiss(self, verbose = 0):
		"""
		residue.BBAtmMiss()
		@return: the positions of all BB atms missing in the atm list
		"""
		missp = []
		for atms in AABB:
			if self.atms.atmPos(atms) == None:
				missp.append(atms)
		if verbose:
			print >> sys.stdout, "# BB atoms missing :",missp
		return missp

	def findAtm(self, atmName = "CA", chId = None, rName = None, rNum = None, icode = None, verbose = 0):
		"""
		residue.findAtm identify an atom
		@param atmName: the atom name
		@param chId: the chain Id
		@param rName: the residu name
		@param rNum: the residu number
		@param icode: the line code number
		@Return: the atom instance
		"""
		if chId != "" and chId != None:
			if self.chnLbl() != chId:
				return None
		if rName != "" and rName != None:
			if self.rName() != rName:
				return None
		if rNum != "" and rNum != None:
			if self.rNum() != rNum:
				return None
		if icode != "" and icode != None:
			if self.riCode() != icode:
				return None
		for aAtm in self.atms:
			if verbose:
				print >> sys.stdout, "# aAtm loop", aAtm.atmName()
			if aAtm.atmName() == atmName:
				return aAtm
		if verbose:
			print >> sys.stdout, "# findAtm: ", atmName,": None"
		return None

	def setBBOrder(self, verbose = 0):
		"""
		residue.setBBOrder set backbone's atoms of a residue in the "right" order
		@author: F.Briand
		@return: none
		"""
		resBB = [None,None,None,None,None]
		resOth = []
		if verbose:
			print >> sys.stdout, "#",
		for atom in self:
			if verbose:
				print >> sys.stdout, atom.atmName(),"",
			try:
				iAtm = BBATMS.index(atom.atmName())
				resBB[iAtm] = atom
			except ValueError:
				resOth.append(atom)
		self.atms = atmList([atom for atom in resBB if atom]+resOth)
		if verbose:
			print >> sys.stdout, " => ",
			for atom in self:
				print >> sys.stdout, atom.atmName(),"",
			print >> sys.stdout

        def isResCode(self, rString, verbose = 0):
                """
                @return:
                True or False, depending on the fact a residue matches rString
                @note: the rString is in the format BBatmMiss 
                @note: For each residue, the string (BBatmMiss) consists of RName_ChLbl_RNum_RIcode, where RName is the name of the residue (3 letters), ChLbl is the chain label, RNum is the number of the residue (as in the PDB file) and RIcode the PDB insertion code of the residue.

                """
                resName = self.rName()
                resNum = self.rNum()
                icode  = self.riCode()
                lbl  = self.chnLbl()
                if icode == ' ':
                        icode = ''
                if lbl == ' ':
                        lbl = ''
                rCode=resName+"_"+lbl+"_"+str(resNum)+"_"+icode
                if verbose:
                        sys.stderr.write("isResCode: comparing \"%s\" and \"%s\"\n" % (rString, rCode))
                if rString == rCode:
                        return True
                return False


## ========================================
## The PDB file parser
## ========================================
class PDB(PDBLine,residue):
	"""
	PDB class
	Manage a PDB
	"""
	def __init__(self, fname = "", chId = "", model = 1, hetSkip = 0, altCare = 0, OXTCare = 0, PDBMac = 0, keepH = 1, id = None, isXML = False, verbose = 0):
		"""
		PDB.__init__ determine the type of fname to initialize
		@param fname: the file name
		@param chId: a chain Id
		@param model: the number of the model you want to set as working model (number 1 by default)
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@param altCare: does it take care about alternate atom? (Yes:1 ; No:0)
		@param OXTCare: does it take care about OXT? (Yes:1 ; No:0)
		@param PDBMac:
		@param keepH: does it keep the hydrogen atoms?
		@param id:
		@param isXML:
		@param altLbl: alternate atom label
		@return: none
		"""
		if fname != "":
			if fname == None:
				return None
			elif isinstance(fname,PDB):		 # already a PDB instance
				if verbose > 1:
					print >> sys.stdout, "# PDB from PDBType. hetSkip : ", hetSkip
				self.info  = fname.info
				self.remark350 = []
				for i in fname.remark350:
					self.remark350.append(i)
				self.id    = fname.id
				self.repository = fname.repository
				self.mdls  = fname.mdls
				# self.atms  = fname.atms
				iPDB = PDB(fname.flat(altCare = altCare, OXTCare = OXTCare, PDBMac = PDBMac, keepH = keepH), model = model, hetSkip = hetSkip, id = fname.id, verbose = verbose)
				self.atms  = iPDB.atms
				# self.data  = fname.data
				self.data  = iPDB.data
				self.mdls  = iPDB.mdls
				self.seq   = iPDB.seq
				self.seq3D = iPDB.seq3D
				self.ss    = iPDB.ss
				self.s2    = iPDB.s2
				self.mtrixn = iPDB.mtrixn
				self.nModel = iPDB.nModel
				self.dbref = iPDB.dbref
				self.chns  = iPDB.chns
				self.conect= iPDB.conect
				self.setModel(model, verbose)
				self.resTab(verbose)
				if altCare :
					self.atms = PDB(self.flat(altCare = 1), verbose = verbose).atms
					self.resTab(verbose)
				if OXTCare :
					# print "altCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)
				if PDBMac :
					# print "PDBMac set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)
			# a flat series of text lines
			elif isinstance(fname,types.ListType):    # a list of atoms
				if verbose > 1:
					print >> sys.stdout, "# PDB from ListType. hetSkip : ", hetSkip
				self.id = "unkwn"
				self.repository = "unkwn"
				self = self.parse(fname, "", chId, hetSkip, verbose)
				if id != None:
					self.id = id
				self.setModel(model, verbose)
				self.resTab(verbose)
				if altCare :
					# print "altCare set"
					self.atms = PDB(self.flat(altCare = 1)).atms
					self.resTab(verbose)
				if PDBMac :
					# print "PDBMac set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)
				if OXTCare :
					# print "OXTCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)
				if not keepH :
					# print "keepH set"
					self.atms = PDB(self.flat(keepH = 0)).atms
					self.resTab(verbose)
			# from disk file
			elif isinstance(fname,types.StringType):  # read file from disk
				if verbose > 1:
					print >> sys.stdout, "# PDB from StringType : "
				if isXML:
					self.loadXML(fname, chId, hetSkip, PDBDIR=GDFLTPDBDIR, verbose = verbose)
				else:
					self.load(fname, chId, hetSkip, PDBDIR=GDFLTPDBDIR, verbose = verbose)
				# if fname[0] in "0123456789" and len(fname < 6):
				#	self.id = fname
				if id != None:
					self.id = id
				self.setModel(model, verbose)
				self.resTab(verbose)
				if PDBMac :
					# print "PDBMac set"
					self.atms = PDB(self.flat(PDBMac = 1)).atms
					self.resTab(verbose)
				if OXTCare :
					# print "OXTCare set"
					self.atms = PDB(self.flat(OXTCare = 1)).atms
					self.resTab(verbose)
				if altCare :
					# print "altCare set"
					self.atms = PDB(self.flat(altCare = 1)).atms
					self.resTab(verbose)
				if not keepH :
					# print "keepH set"
					self.atms = PDB(self.flat(keepH = 0)).atms
					self.resTab(verbose)

	def __getslice__(self, ffrom = 0, tto = None):
		"""
		PDB.__getslice__(ffrom = 0, tto = None)
		@param ffrom: the first residu considered
		@param tto: the last residu of the residues (excluded, as in python)
		@return: a PDB instance, a slice of residues
		"""
		if tto == None:
			tto = len(self)
		res = self[ffrom].flat()
		for i in range(ffrom+1,tto):
			res= res + self[i].flat()
		return PDB(res)

	def __getitem__(self,aPos):
		"""
		PDB.__getitem__(aPos)
		@param aPos: the position of one residue (PDB[aPos]) or CHAINS (PDB["CHAINS"])
		@return: one residue or PDB instance of chains matching CHAINS (e.g. \"AB\")
		"""
		if isinstance(aPos,types.IntType):
			return self.rt[aPos]
		elif isinstance(aPos,types.StringType):
			return self.chn(aPos)

	def __len__(self):
		"""
		PDB.__len__()
		@return: number of residues of the PDB instance
		"""
		return len(self.rt)

	def __add__(self,new):
		"""
		PDB.__add__()
		@param new: a PDB to concatenate to a first one
		@return: the 2 PDBs concatenated
		I{in this example, it will return z: } \n
		x = PDB()\n
		new = PDB()\n
		z = x + new\n
		"""
		addPDB = PDB(self)
		addPDB.atms = self.atms + new.atms
		addPDB.resTab(0)
		return addPDB

	def __delitem__(self,aPos):
		"""
		del x[i]
		@param aPos: position of the residu
		@note: preserves anything in comments: modify atms then performs resTab() MIGHT BE BUGGY (P. Tuffery, 2007)
		"""
		try:
			rName  = self[aPos].rName()
			rNum   = self[aPos].rNum()
			riCode = self[aPos].riCode()
			rChn   = self[aPos].chnLbl()
		except:
			return
		for iAtm in range(len(self.atms)-1, -1, -1):
			aAtm = self.atms[iAtm]
			resName = aAtm.resName()
			resNum  = aAtm.resNum()
			iCode   = aAtm.icode()
			chn     = aAtm.chnLbl()
			if (resNum == rNum) and (resName == rName) and (iCode == riCode) and (chn == rChn):
				del self.atms[iAtm]
		self.resTab(0)

	def __repr__(self, altCare = 0, altLbl = "", OXTCare = 0, HSkip = 0):
		"""
		PDB.__repr__
		@param altCare: does it take care about alternate atoms? (Yes:1 ; No:0)
		@param OXTCare: does it take care about OXT? (Yes:1 ; No:0)
		@param HSkip: does it skip hydrogen? (Yes:1 ; No:0)
		@param altLbl: alternate atom label
		@return: show atomic information of PDB
		(print PDB ATOM lines)
		"""
		res = ""
		i = 0
		curChn = ""
		for aRes in self:
			if (i > 0) and (curChn != aRes.chnLbl()):
				res = res + "TER \n"
			curChn = aRes.chnLbl()
			res = res + aRes.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip)
			i += 1
		return res

	# back to a string
	# an internal vital function to transit from PDB to atmList etc
	def flat(self, altCare = 0, altLbl = "", OXTCare = 0, PDBMac = 0, keepH = 1):
		"""
		PDB list flat
		@param altCare: does it take care about alternate atom? (Yes:1 ; No:0)
		@param altLbl: alternate atom label
		@param OXTCare: does it take care about OXT? (Yes:1 ; No:0)
		@return: an atmList
		"""
		res = []
		for i in self:
			res = res + i.flat(altCare, altLbl, OXTCare = OXTCare, PDBMac = PDBMac, keepH = keepH)
		return res

	def out(self, outName = "", chainId = "", altCare = 0, altLbl = "", OXTCare = 0, hetSkip = 0, fmode = "w", header = 1, ter = 1, model = 0, end = 0, info = 0, HSkip = 0, allModels = 0, atmFrom = None, verbose = 0):
		"""
		PDB.out write a PDB
		@param outName: the name of the file written, if none it sets the standard out
		@param chainId:
		@param altCare: does the file contains the alternate atoms (Yes:1 by default)
		@param altLbl: does the file contains the alternate label (Yes:1 by default)
		@param OXTCare: does the file contains the information of each line (Yes:1 by default)
		@param hetSkip: does the file skip the het (No:0 by default)
		@param fmode: the opening file method, "w" writing by default
		@param header: does the file contains the header (Yes:1 by default)
		@param ter: does the file keep TER lines (Yes:1 by default)
		@param model: the number of the model that you want to write in the file
		@param end: does the file keep END line (No:0 by default)
		@param info: does the file contains the information of each line (Yes:1 by default)
		@param HSkip: does the file skip the hydrogen (No:0 by default)
		@param allModels: does the file contains all the models existing (No:0 by default)
		@return: none, it outputs PDB content to file
		"""
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName,fmode)
			except:
				raise IOError,"Failed to write to %s\n" % outName
		if verbose > 1:
			print >> sys.stderr, "# Writing in %s" % f.name
		# print "Opened :",outName,fmode
		# f = sys.stdout
		if header and not info:
			f.write("HEADER                                                        %s\n" % self.id)
		if info:
			rmk350 = False
			rmk = False
			formul = False
			for aLine in self.info:
				if not rmk and aLine[:6] == "REMARK":
					rmk = True
				try:
					if not rmk350 and aLine[:6] == "REMARK" and int(aLine[6:10]) > 350:
						for rmkLine in self.divRemark350:
							f.write("%s" % rmkLine)
						for rmk in self.remark350:
							for rmkLine in rmk.txt:
								f.write("%s" % rmkLine)
						rmk350 = True
				except:
					pass
				if rmk and aLine[:6] != "REMARK":
					for iLine in self.dbref:
						f.write(iLine)
					for iLine in self.seq:
						f.write(iLine)
					rmk = False
				if not formul and aLine[:6] == "FORMUL":
					formul = True
				if formul and aLine[:6] != "FORMUL":
					for iLine in self.s2:
						f.write(iLine)
					formul = False
				f.write("%s" % aLine)
		if allModels:
			for aModel in range(1, self.nModel + 1):
				self.setModel(model = aModel)
				if atmFrom is not None:
					self.atmsRenumber2(atmFrom)
				f.write("MODEL %d \n" % (aModel))
				res = self.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip )
				f.write("%s" % res)
				f.write("ENDMDL\n")
		else:
			if model:
				f.write("MODEL %d \n" % (model))
			res = self.__repr__(altCare, altLbl, OXTCare = OXTCare, HSkip = HSkip )
			f.write("%s" % res)
			if model:
				f.write("ENDMDL\n")
		if info:
			for aLine in self.ss:
				f.write(aLine)
			for aLine in self.s2:
				f.write(aLine)
			for aLine in self.conect:
				f.write(aLine)
		if end:
			f.write("END\n")
		f.flush()
		if f != sys.stdout:
			f.close()
		else:
			sys.stdout.flush()

	def xyzout(self, outName = "", chainId = "", hetSkip = 0, fmode = "w", verbose = 0):
		"""
		xyzout writes atoms coordinates in a file
		@param outName: the name of the file written, if none it sets the standard out
		@param fmode: the opening file method, "w" writing by default
		@return: none, the file "outName" is written
		"""
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
		if outName == "":
			f = sys.stdout
		else:
			try:
				f = open(outName,fmode)
			except:
				raise IOError, "Failed to write to %s" % outName
		if verbose > 1:
			print >> sys.stderr, "# Writing in %s" % f.name
		for aCrd in res:
			f.write("%s\n" % aCrd)
		if f != sys.stdout:
			f.close()
		else:
			sys.stdout.flush()

	def xyz(self, outName = "", chainId = "", hetSkip = 0, fmode = "w", verbose = 0):
		"""
		xyz
		@param outName: the name of the file written, if none it sets the standard out
		@param chainId:
		@param hetSkip:
		@param fmode: the opening file method, "w" writing by default
		@return: the coordinates of all the atoms, they are concatenated
		"""
		res = []
		for aRes in self:
			res = res + aRes.atms.crds()
		return res

	def loadXML(self,fname, chainId = "", hetSkip = 0, PDBDIR = GDFLTPDBDIR, CATHDIR = GDFLTCATHDIR, SCOPDIR = GDFLTSCOPDIR, verbose = 0, model = 1):
		"""
		loadXML currently converts atoms informations in flat and parses them. You must add "isXML=True" to load a XML file
		@param fname: the name of the file
		@param chainId: the Chain id
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@param PDBDIR: equal to GDFLTPDBDIR
		@param CATHDIR: equal to GDFLTCATHDIR
		@param SCOPDIR: equal to GDFLTSCOPDIR
		@return: the atoms informations parsed
		@note: it does not support/try cath, scop and astral
		"""
		from Ft.Xml.XPath import Evaluate
		from Ft.Xml.XPath.Context import Context
		from Ft.Xml.Domlette import NonvalidatingReader
		from Ft.Lib import Uri

		file_uri = Uri.OsPathToUri(fname)
		doc = NonvalidatingReader.parseUri(file_uri)

		nsBindings = {'PDBx': 'http://pdbml.pdb.org/schema/pdbx-v32.xsd'}

		context = Context(None, processorNss=nsBindings)
		path="/PDBx:datablock/PDBx:atom_siteCategory/PDBx:atom_site"
		nodes=Evaluate(path,doc,context)

		L=[]
		for n in nodes:
			atom={}
			atom["recname"] = n.xpath('PDBx:group_PDB/text()')[0].data
			atom["id"] = n.xpath('@id')[0].nodeValue
			atom["name"] = n.xpath('PDBx:auth_atom_id/text()')[0].data
			atom["res"] = n.xpath('PDBx:auth_comp_id/text()')[0].data
			atom["asym"] = n.xpath('PDBx:auth_asym_id/text()')[0].data
			atom["resnum"] = n.xpath('PDBx:auth_seq_id/text()')[0].data
			atom["Cartn_x"] = n.xpath('PDBx:Cartn_x/text()')[0].data
			atom["Cartn_y"] = n.xpath('PDBx:Cartn_y/text()')[0].data
			atom["Cartn_z"] = n.xpath('PDBx:Cartn_z/text()')[0].data
			atom["occupancy"] = n.xpath('PDBx:occupancy/text()')[0].data
			atom["B_iso"] = n.xpath('PDBx:B_iso_or_equiv/text()')[0].data
			atom["tfac"] = " "
			atom["type_symbol"] = n.xpath('PDBx:type_symbol/text()')[0].data
			try:
				atom["charge"] = n.xpath('PDBx:pdbx_formal_charge/text()')[0].data
			except:
				atom["charge"] = " "
			#try:
			#	atom["label_alt_id"] = n.xpath('PDBx:label_alt_id/text()')[0].data
			#except:
			#	atom["label_alt_id"] = " "
			#atom["label_asym"] = n.xpath('PDBx:label_asym_id/text()')[0].data
			#atom["label_atom_id"] = n.xpath('PDBx:label_atom_id/text()')[0].data
			#atom["label_comp_id"] = n.xpath('PDBx:label_comp_id/text()')[0].data
			#atom["label_entity_id"] = n.xpath('PDBx:label_entity_id/text()')[0].data
			#atom["label_seq_id"] = n.xpath('PDBx:label_seq_id/text()')[0].data
			#atom["model_num"] = n.xpath('PDBx:pdbx_PDB_model_num/text()')[0].data
			#print atom
			L.append(atom)

		D=[]
		for i in L:
			stringD = "%-6s%5s  %-4s%3s %1s%4s%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s%2s\n"%(i['recname'],i['id'],i['name'],i['res'],i["asym"],i["resnum"],"",float(i["Cartn_x"]),float(i["Cartn_y"]),float(i["Cartn_z"]),float(i["occupancy"]),float(i["B_iso"]),i["type_symbol"],i["charge"])
			D.append(stringD.encode('Utf8'))
#			print "%-6s%5s  %-4s%3s %1s%4s%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s%2s"%(i['recname'],i['id'],i['name'],i['res'],i["asym"],i["resnum"],"",float(i["Cartn_x"]),float(i["Cartn_y"]),float(i["Cartn_z"]),float(i["occupancy"]),float(i["B_iso"]),i["type_symbol"],i["charge"])
			#print "%-6s%5s  %-4s%3s %1s%4s%1s   %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2s%2s  %4s%-4s%3s %1s%4s%4s %3s"%(i['recname'],i['id'],i['name'],i['res'],i["asym"],i["resnum"],"",float(i["Cartn_x"]),float(i["Cartn_y"]),float(i["Cartn_z"]),float(i["occupancy"]),float(i["B_iso"]),i["type_symbol"],i["charge"],i["label_alt_id"],i["label_atom_id"],i["label_comp_id"],i["label_asym"],i["label_seq_id"],i["label_seq_id"],i["model_num"])
		#print D

		# Organize series of lines
		idName = fname
		if string.find(fname,"/") != -1:
			idName = fname[string.rindex(fname,"/")+1:]
		(prefix, tail) = os.path.split(fname)
		(Id, tail) =  os.path.splitext(tail)
			## self.parse(allPDB, idName[:40]+"_"+chainId, chainId, hetSkip, verbose, model)
		self.parse(D, Id[:40], chainId, hetSkip, verbose, model)
		return self

	def load(self,fname, chainId = "", hetSkip = 0, PDBDIR = GDFLTPDBDIR, CATHDIR = GDFLTCATHDIR, SCOPDIR = GDFLTSCOPDIR, verbose = 0, model = 1):
		"""
		load parse a PDB file
		@param fname: the name of the file
		@param chainId: the Chain id
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@param PDBDIR: equal to GDFLTPDBDIR
		@param CATHDIR: equal to GDFLTCATHDIR
		@param SCOPDIR: equal to GDFLTSCOPDIR
		@return: the atoms informations parsed
		"""
		try:
			if verbose > 1:
				print >> sys.stdout, "# Trying: %s \n" % fname
			allPDB=gsimpleload(fname, verbose=verbose)
			self.repository = "file"
			if verbose > 1:
				print >> sys.stdout, "# Succeeded!\n"
		except IOError, UnboundLocalError:
			if verbose:
				print >> sys.stderr, "# %s: not a local file ... or some problem occurred ...\n" % fname
			if fname[0] in "0123456789" and (len(fname) == 7) and (fname[-1] in "0123456789") and (fname[-2] in "0123456789"):
				try:
					if verbose > 1:
						print >> sys.stdout, "# Trying local CATH entry: %s/%s\n" % (CATHDIR,fname )
					allPDB=gsimpleload("%s/%s1" % (CATHDIR,fname),0)
					self.repository = "cath"
				except :
					try:
						if verbose > 1:
							print >> sys.stdout, "# Trying remote CATH entry: http://data.cathdb.info/v3_2_0/pdb/%s\n" % (fname )
						from urllib import urlretrieve
						file, log = urlretrieve("http://data.cathdb.info/v3_2_0/pdb/%s" % fname)
						allPDB=simpleload(file,0)
						if verbose > 2:
							print >> sys.stdout, "# urlretrieve at %s" % file
						del urlretrieve
						self.repository = "cath"
						# print len(allPDB)
						if len(allPDB) < 1:
							raise IOError
					except:
						raise UnboundLocalError, 'Sorry: PDB entry %s not found' % pdbEntry
			elif fname[0] in "0123456789":
				pdbEntry = fname[:4]
				if chainId == "":
					chainId = fname[4:]
				try:
				## Experimental structure
					if verbose > 1:
						print >> sys.stdout, "# Trying local PDB copy: %s/all/pdb/pdb%s.ent.Z\n" % (PDBDIR, pdbEntry)
					try:
						# Old PDB is .Z
						allPDB=gsimpleload(PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.Z",0)
						self.repository = "PDB"
					except:
						# wwPDB is .gz
						allPDB=gsimpleload(PDBDIR+"/all/pdb/pdb"+pdbEntry+".ent.gz",0)
						self.repository = "PDB"
				except IOError:
					if verbose > 1:
						print >> sys.stderr, "# Failed"
				## Model structure
					try:
						if verbose > 1:
							print >> sys.stdout, "# Trying local PDB model: %s/models/current/pdb/%s/pdb%s.ent.Z\n" % (PDBDIR,pdbEntry[1:3],pdbEntry )
						try:
							allPDB=gsimpleload(PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.Z",0)
							self.repository = "PDB"
						except:
							allPDB=gsimpleload(PDBDIR+"/models/current/pdb/"+pdbEntry[1:3]+"/pdb"+pdbEntry+".ent.gz",0)
							self.repository = "PDB"
					except IOError:
						try:
							if verbose > 1:
								print >> sys.stdout, "# Attempting: PDB from wwpdb.org for %s\n" % pdbEntry
							from urllib import urlretrieve
							# file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(pdbEntry,pdbEntry))
							# We query the RCSB
							# file, log = urlretrieve("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" %(pdbEntry))
							# We query the wwpdb
							file, log = urlretrieve("ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb%s.ent.gz" % pdbEntry.lower())
							subprocess.call("chmod 777 %s" % file, shell=True)
							allPDB=gsimpleload(file,verbose)
							if verbose > 1:
								print >> sys.stdout, "# urlretrieve at %s" % file
							del urlretrieve
							self.repository = "PDB"
							# return
						except IOError:
							error_txt = '# Sorry: PDB entry %s not found' % pdbEntry
                                                        sys.stderr.write("%s\n" % error_txt)
							raise UnboundLocalError # , error_txt
			elif fname[0] in "deg":
				## Astral / scop
				try:
					subdir  = fname[2:4]
					lname = fname
					if fname[0] != "d":
						lname = "d"+fname[1:]
					if verbose > 1:
						print >> sys.stdout, "# Trying Astral/Scop local copy: %s/%s/%s.ent\n" % (SCOPDIR, subdir, lname)
					allPDB=gsimpleload("%s/%s/%s.ent" % (SCOPDIR, subdir, lname),0)
					self.repository = "Astral"
				except:
					try:
						# We try for astral on the net.
						if verbose > 1:
							print >> sys.stdout, "# Attempting: Astral/Scop entry from astral.berkeley.edu %s\n" % lname

						from urllib import urlretrieve
						# file, log = urlretrieve("http://astral.berkeley.edu/pdbstyle.cgi?id=%s&output=text" % lname)
						file, log = urlretrieve("http://scop.berkeley.edu/astral/pdbstyle/ver=2.03&id=%s&output=text" % lname)
						allPDB=simpleload(file,0)
						del urlretrieve
						self.repository = "Astral"
						if verbose > 2:
							print >> sys.stdout, "# Astral entry %s retrieved via the net\n" % lname
					except:
						# raise UnboundLocalError, '# Sorry: Astral/SCOP entry % not found' % lname
						print >> sys.stderr, '# Sorry: Astral/SCOP entry % not found' % lname
						### /!\ No return ?
			else:
				print >> sys.stderr, "# Sorry: %s does not sound as PDB or Astral/SCOP entry!" % fname
				return self

		# Organize series of lines
		idName = fname
		if string.find(fname,"/") != -1:
			idName = fname[string.rindex(fname,"/")+1:]
		(prefix, tail) = os.path.split(fname)
		(Id, tail) =  os.path.splitext(tail)

		self.parse(allPDB, Id[:40], chainId, hetSkip, verbose, model)
		return self

	def parse(self, allPDB, id="", chainId = "",  hetSkip = 0, verbose = 0, model = 1):
		"""
		PDB.parse(self, allPDB, id="", chainId = "",  hetSkip = 0, verbose = 0, model = 1)
		@param allPDB: the flat lines to format
		@param id: the id of PDB
		@param chainId: the Chain id
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@param model: the model sets as working set (1 by default)
		@return: a PDB format
		"""
		if id == "":
			id = "unkwn"
		self.info  = []
		self.remark350 = []
		self.divRemark350 = []
		self.id    = string.replace(id," ","")
		self.data  = []   # All ATOM DATA, N MODELS
		self.mdls  = []
		self.atms  = []
		self.seq   = []
		self.seq3D = []
		self.ss    = []
		self.s2    = []
		self.mtrixn = []
		self.nModel = 0
		self.mdls.append(0)
		self.dbref = []
		self.chns  = ""
		self.conect = []
		for curLine in allPDB:
			aLine = PDBLine(curLine)
			header = aLine.header()
			if header == "ATOM" or header == "HETATM":
				aLine = atmLine(aLine)
				OK = 0
				if chainId == "":
					OK = 1
				elif chainId[0] != '-':
					if string.count(chainId, aLine.chnLbl()):
						OK = 1
				else:
					if string.count(chainId, aLine.chnLbl()) == 0:
						OK = 1
				if OK:
					if hetSkip:
						if AA3.count(aLine.resName()) > 0:
							self.data.append(aLine)
						elif hetSkip == 2:
							if SOLV.count(aLine.resName()) == 0:
								self.data.append(aLine)
					else:
						self.data.append(aLine)
			elif header == "TER":
				## self.data.append(curLine)
				pass
			elif header == "HEADER":
				self.info.append(curLine)
				try:
					if curLine[62:66] != "    ":
						self.id = string.split(curLine[62:])[0]
				except:
					pass
			elif header == "COMPND":
				self.info.append(curLine)
			elif header == "SOURCE":
				self.info.append(curLine)
			elif header == "REMARK":
				try:
					if int(curLine[7:10]) == 350:
						self.remark350.append(curLine)
					else:
						self.info.append(curLine)
				except:
					self.info.append(curLine)
			elif header == "SEQRES":
				self.seq.append(curLine)
			elif header == "HELIX" or header == "SHEET" or header == "TURN":
				## self.s2.append(allPDB[aLine])
				self.s2.append(curLine)
			elif header == "SSBOND":
				## self.ss.append(curLine)
				self.ss.append(curLine)
			elif header == "DBREF":
				## self.dbref = allPDB[aLine]
				self.dbref.append(curLine)
			elif header[:-1] == "MTRIX":
				self.mtrixn.append(curLine)
			elif header == "ENDMDL":
				## self.mdls.append(len(self.data))
##				## self.nModel = self.nModel+1
				self.mdls.append(len(self.data))
				self.nModel = self.nModel+1
			elif header == "CONECT":
				self.conect.append(curLine)
			else:
				if header not in ["MASTER", "END", "MODEL"]:
					self.info.append(curLine)
		if self.remark350:
			curRemark350 = []
			remark350Lines = self.remark350
			self.remark350 = []
			for line in remark350Lines:
				if line[11:23] == "BIOMOLECULE:":
					if curRemark350:
						self.remark350.append(remark350(curRemark350))
						curRemark350 = []
					curRemark350.append(line)
				else:
					if curRemark350:
						curRemark350.append(line)
					else:
						self.divRemark350.append(line)
			self.remark350.append(remark350(curRemark350))
		if self.nModel == 0:
			self.nModel = self.nModel+1
		self.mdls.append(len(self.data))
		return self

	def resTab(self, verbose):
		"""
		PDB.resTab reformats the data for an easier access
		@return: none
		@note: update PDB.rt atoms from PDB.atms atoms
		"""
		start   = 1
		self.rt = []
		curResNum = "-1000"
		curResName = "XXX"
		curICode  = ""
		curChn = ""
		atmFrom = 0
		if len(self.atms) == 0:
			if verbose:
				print >> sys.stderr, "# Empty PDB instance\n"
			return
		for iAtm in range(0,len(self.atms)):
			aAtm = self.atms[iAtm]
			resName = aAtm.resName()
			resNum  = aAtm.resNum()
			iCode   = aAtm.icode()
			chn     = aAtm.chnLbl()
			if resNum != curResNum or resName != curResName or iCode != curICode or chn != curChn:
				curResNum = resNum
				curResName = resName
				curICode = iCode
				curChn = chn
				if start:
					start = 0
				else:
					self.rt.append(residue(self.atms[atmFrom:iAtm]))
					atmFrom = iAtm
		self.rt.append(residue(self.atms[atmFrom:iAtm+1]))
		if verbose:
			print >> sys.stdout, "# Found %s residues" % len(self.rt)

	def atmTab(self):
		"""
		PDB.atmTab reformats the data for an easier access
		@author: F.Briand
		@return: none
		@note: update PDB.atms atoms from PDB.rt atoms
		"""
		atoms = []
		for res in self.rt:
			for atom in res:
				atoms.append(atom)
		self.atms = atoms

	def nModels(self):
		"""
		PDB.nmodels()
		@return: number of models of PDB instance (as defined by MODEL / ENDMDL lines)
		"""
		return self.nModel

	# Install current model
	# x = PDB("/home/raid5/PDB/pdb1g25.ent.gz", hetSkip = 1)
	# /!\ "delete" modifications make on the current model and set "originals" atoms (from PDB.data). Maybe make a dataTab (like resTab, atmTab ...) ?
	def setModel(self,model = 1,verbose = 0):
		"""
		This will install the model of rank specified by "model" as the current working set.
		PDB.setModel(model = 1,verbose = 0)
		@param model: the number of the model you want to set as working model
		@return: none
		"""
		#set model \# model as working model
		if model > self.nModels():
			print >> sys.stderr, "# Sorry: no model number %d (Total of %d)" % (model,self.nModels())
			return
		self.atms = []
		#print "before: self.atoms",self.atms
		if verbose:
			print >> sys.stdout, "# Installing model %d (atoms %d - %d)" % (model, self.mdls[model-1], self.mdls[model])
		for aLine in range(self.mdls[model-1],self.mdls[model]):
			self.atms.append(atmLine(self.data[aLine]))
		self.curMdl = model
		#print "after: self.atoms",self.atms
		self.resTab(0)
		return

	# "update" the PDB.chns parameter
	def chnList(self):
		"""
		PDB.chnList()
		@return: a string, the concatenated chain identifiers present in the PDB (including blank).
		"""
		curChn = ""
		self.chns = ""
		for aLine in range(0,len(self.atms)):
			if string.count(self.chns,self.atms[aLine][21]) == 0:
				curChn = self.atms[aLine][21]
				self.chns = self.chns + curChn
		return self.chns

	def nChn(self):
		"""
		PDB.nChn()
		@return: a number, the number of chain identifiers present in the PDB (including blank).
		"""
		if self.chns == "":
			return len(self.chnList())
		else:
			return len(self.chns)

	def hasChn(self, chnId):
		"""
		PDB.hasChn check if chain of id chainId is present in the PDB instance.
		@param chnId: the id of a chain
		@return: number of chnId in the list of chain
		"""
		if self.chns == "":
			return string.count(self.chnList(),chnId)
		else:
			return string.count(self.chns,chnId)

	# DEBUG VERSION (Florian Briand - Mars 2010)
	# => Goal : Keep all informations (about header, chains ...) with selection by chn()
	#
	# extract particular chain(s) passed in string chainId
	# the default is to return all the chains
	#
	def chn(self,chainId=None, hetSkip = 0, inplace = 0, verbose = 0):
		"""
		PDB.chn(chnId, hetSkip = 0)
		@author: F.Briand
		@param chainId: the id of a chain
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@param inplace: does the chn function will replace 'self' by the truncated PDB (No:0 by default)
		@return: a PDB instance of chains of the PDB.
		@type chainId: might contain several chain Ids (e.g. AB)
		@note: if chainId starts by \"-\" (minus) returns all but the chains specified.
		"""
		if not chainId:
			if verbose:
				print >> sys.stdout, "# Selecting all chains"
			return PDB(self, hetSkip=hetSkip)
		if verbose:
			print >> sys.stdout, "# Selecting %s " % chainId,
		if inplace:
			if verbose:
				print >> sys.stdout, "in place"
			chain = self
		else:
			if verbose:
				print >> sys.stdout, "in a new instance"
			chain = PDB(self)
		chainId = chainId.split("-")
		if len(chainId)==1:
			chainId.append("")
		if chainId[0] == "":
			chainId[0] = chain.chnList()
		for i in range(0,len(chain)):
			for iAtm in range(len(chain[i].atms)-1,-1,-1):
				if (chain[i].atms[iAtm].chnLbl() not in chainId[0]) or (chain[i].atms[iAtm].chnLbl() in chainId[1]):
					chain[i].atms.list.remove(chain[i].atms.list[iAtm])
		return PDB(chain, hetSkip=hetSkip)

	def chnRename(self, pattern='', verbose=0):
		"""
		PDB.chnRename Rename chains
		@author: F.Briand
		@param pattern: "chain1chain2:newChain1newChain2"
		@return: none
		@note: edit the PDB class "in place", ":newChain" give a name to a "one chain PDB"
		"""
		if pattern=='':
			print >> sys.stderr, "# PDB.chnRename() : pattern missing. USAGE : \"chain1chain2...:newChain1newChain2...\""
			return
		chainIn, chainOut = pattern.split(":")
		if len(chainIn) != len(chainOut):
			if chainIn:
				print >> sys.stderr, "# PDB.chnRename() : not the same number of chain and newChain"
				return
			elif len(self.chnList()) > 0:
				print >> sys.stderr, "# PDB.chnRename() : PDB contain more than one chain, please make a specific pattern"
				return
		for atom in self.atms:
			if chainIn != "":
				try:
					atom.chnLbl(chainOut[chainIn.index(atom.chnLbl())])
				except ValueError:
					pass
			else:
				atom.chnLbl(lbl=chainOut[0])
		if verbose:
			print >> sys.stdout, "# chnRename : %s => %s" % (chainIn, chainOut)
		self.chnList() # update self.chns
		self.resTab(0)

	def renameHydrogens(self, norm, verbose = 0):
		"""
		PDB.renameHydrogens
		@author: F.Briand
		@param norm: IUPAC to set Hydrogens to IUPAC norm, PDB to set Hydrogens to PDB norm
		@param verbose: set verbose mode
		@return: none
		@note: It only works for standard amino-acids.
		"""
		# Setting opposite norm to make conversion.
		# /!\ Maybe set opposite norm with PDB.hNorm() or use 2 params in input
		if norm == 'lPDB':
			oppositeNorm = 'lIUPAC'
		elif norm == 'lIUPAC':
			oppositeNorm = 'lPDB'
		else:
			raise ValueError, "%s is not a recognized norm" % norm
		selfnorm = self.hNorm()
		if selfnorm == norm:
			if verbose:
				print >> sys.stdout, "# Hydrogens already in %s norm" % norm
			return
		else:
			for atom in self.atms:
				# For each atom, check if this is an Hydrogen
				if atom.atmType() == 'H':
					rName = atom.resName()
					aName = atom.atmName()
					# Check if the selected hydrogen is an atom not in the wished norm
					if aName not in normHNames[norm][0][rName] and aName not in normHNames[norm][0]['BCK']:
						try:
							# Take the "type of hydrogen" in oppositeNorm and apply the corresponding hydrogen in norm
							hIndex = normHNames[oppositeNorm][0][rName].index(aName)
							if verbose:
								print >> sys.stdout, "# %s %s %s => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
							atom.atmName(makeHName(rName,hIndex,norm))
						except ValueError:
							try:
								# If this type of hydrogen not in "lateral chain" hydrogens, it can be in "backbone" hydrogens
								hIndex = normHNames[oppositeNorm][0]['BCK'].index(aName)
								if verbose:
									print >> sys.stdout, "# %s %s %s-BCK => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
								atom.atmName(makeHName('BCK',hIndex,norm))
							except ValueError:
								# If not in existing oppositeNorm, maybe a "modified" atom name, not totally standard => "factorization"
								nameSplit = re.split('([A-Z]+)', aName)
								hName = ''
								hNum = []
								for item in nameSplit:
									if item:
										try:
											if len(item)>1:
												hNum.append(int(item[0]))
												hNum.append(int(item[1]))
											else:
												if not hName or len(hNum):
													hNum.append(int(item))
												else:
													hNum.append(int(item))
													hNum[0:0] = [0]
										except ValueError:
											hName = item
									hTup = (hName,hNum[0],hNum[1])

								if hTup in normHNames[norm][1][rName]:
									hIndex = normHNames[norm][1][rName].index(hTup)
									if verbose:
										print >> sys.stdout, "# %s %s %s => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
									atom.atmName(makeHName(rName,hIndex,norm))

								elif hTup in normHNames[norm][1]['BCK']:
									hIndex = normHNames[norm][1]['BCK'].index(hTup)
									if verbose:
										print >> sys.stdout, "# %s %s %s-BCK => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
									atom.atmName(makeHName('BCK',hIndex,norm))

								elif hTup in normHNames[oppositeNorm][1][rName]:
									hIndex = normHNames[oppositeNorm][1][rName].index(hTup)
									if verbose:
										print >> sys.stdout, "# %s %s %s => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
									atom.atmName(makeHName(rName,hIndex,norm))

								elif hTup in normHNames[oppositeNorm][1]['BCK']:
									hIndex = normHNames[oppositeNorm][1]['BCK'].index(hTup)
									if verbose:
										print >> sys.stdout, "# %s %s %s-BCK => %s" % (aName, atom.atmNum(), rName, makeHName(rName,hIndex,norm))
									atom.atmName(makeHName('BCK',hIndex,norm))

								else:
									if verbose:
										print >> sys.stderr, "# Error in HName change for :\n#", atom
		self.resTab(0)

	def hNorm(self):
		"""
		PDB.hNorm give the hydrogen naming convention of the PDB.
		@author: F.Briand
		@return : the norm in {lIUPAC, lPDB, lpIUPAC, lpPDB, None}
		@note : "l" means legacy, "p" means "partial", None if norms are "indifferentiable"
		"""
		lIUPAC = 0
		lPDB = 0
		UNK = 0
		for atom in self.atms:
			# for each hydrogen, look the norm containing its name
			if atom.atmType() == 'H':
				aName = atom.atmName()
				rName = atom.resName()
				ilPDB, ilIUPAC, iUNK = False,False,False
				if aName in normHNames['lPDB'][0][rName] or aName in normHNames['lPDB'][0]['BCK']:
					ilPDB = True
				if aName in normHNames['lIUPAC'][0][rName] or aName in normHNames['lIUPAC'][0]['BCK']:
					ilIUPAC = True
				if not ilIUPAC^ilPDB:
					if not ilIUPAC:
						UNK += 1
				else:
					lIUPAC, lPDB = lIUPAC+int(ilIUPAC),lPDB+int(ilPDB)
		lIUPAC += UNK
		lPDB += UNK
		# if more than 1 norm have been found, it's a "partial" norm which his returned
		if not bool(lIUPAC)^bool(lPDB):
			if lIUPAC > lPDB:
				return "plIUPAC"
			elif lIUPAC < lPDB:
				return "plPDB"
			else:
				return None
		else:
			if lIUPAC:
				return "lIUPAC"
			else:
				return "lPDB"

	def chnType(self, chainId = "", verbose = 0):
		"""
		PDB.chnType give the molecular type of chain(s)
		@param chainId: the chain id
		@return: Heuristic detection if the chain type is one of:
			- Protein
			- RNA
			- DNA
			- SOLVENT
			- HETERO
		"""
		if chainId == "":
			chainId = self.chnList()
		res = []
		unres = []
		for aChain in chainId:
			theChain = self.chn(aChain)
			nAA  = 0
			nRNA  = 0
			nDNA  = 0
			nHET = 0
			nH2O = 0
			for i in range(0,len(theChain)):
				resName = theChain[i].rName()
				if AA3.count(resName) > 0:
					nAA = nAA +1
				elif RNA3.count(string.split(resName)[0]) > 0:
					nRNA = nRNA +1
				elif DNA3.count(string.split(resName)[0]) > 0:
					nDNA = nDNA +1
				elif SOLV.count(string.split(resName)[0]) > 0:
					nH2O = nH2O +1
				else:
					nHET = nHET + 1
					if verbose:
						if unres.count(resName) == 0:
							unres.append(resName)
							print >> sys.stderr, "# Unknown residue type : %s" % resName
			if verbose:
				print >> sys.stdout, "# nAA : ",nAA," nNA : ",nDNA + nRNA," nHET : ",nHET
			nOTHER = nHET + nDNA + nRNA
			if nOTHER < nAA:
				res.append("Protein")
			elif nAA > nDNA + nRNA:
				res.append("Protein")
			# elif nRNA + nDNA > nHET:
			elif nRNA + nDNA > 0:
				if nRNA > 0:
					res.append("RNA")
				else:
					res.append("DNA")
			else:
				if nH2O > nHET:
					res.append("SOLVENT")
				else:
					res.append("HETERO")
		if unres:
			print >> sys.stderr, "# Unknown residues : %s" % unres
		if len(chainId) == 1:
			return res[0]
		return res

	def resTypes(self, what = "all", types = "all", solvent = True):
		"""
		PDB.resTypes return a list of the residue names
		@param what: could be "all" or "aminoacid" or "nucleicacid"
		@param types: could be "all" or "none" or "std" or "lstd" or "het" or "shet"
		@note:
			- "all" : (default) all residue types.
			- "none": no residue types (only solvent mask is effective).
			- "std" : all standard residue types (standard amino acids, standard nucleotides).
			- "lstd": all amino acid types in addition to std.
			- "het" : all non standard residues (includes non standard amino-acids)
			- "shet": true heteros groups (does not include non standard amino-acids)
		@param solvent: consider solvent (true by default)
		@return: a list of residue names (3 characters), combining mask values of what and solvent
		"""
		rs = []
		for i in self:
			rName = i.rName()
			if (solvent == False) and (rName in SOLV):
				continue
			if (types == "none") and (rName not in SOLV):
				continue
			if (types == "std") and ((rName not in AA3STRICT) and (rName not in DNA3) and (rName not in RNA3)):
				continue
			if (types == "lstd") and ((rName not in AA3) and (rName not in DNA3) and (rName not in RNA3)):
				continue
			if (types == "het") and ((rName in AA3STRICT) or (rName in DNA3) or (rName in RNA3)):
				continue
			if (types == "shet") and ((rName in AA3) or (rName in DNA3) or (rName in RNA3)):
				continue
			if rName not in rs:
				rs.append(rName)
		return rs

	def select(self,rwhat=[""],awhat=[""]):
		"""
		PDB.select return a selection of (sub) residues
		@param rwhat: none, or what residues
		@param awhat: none, or what atoms
		@return: a selection of (sub) residues
		"""
		res = []
		for i in self:
			if rwhat == [""]:
				res = res + i.select(awhat).flat()
			elif rwhat[0] !=  "-":
				if rwhat.count(i.rName()) > 0:
					res = res + i.select(awhat).flat()
			else:
				if rwhat.count(i.rName()) == 0:
					res = res + i.select(awhat).flat()
		if res == []:
			return None
		return PDB(res, id = self.id)

	def mask(self,ffrom=0,tto=-1,mask=""):
		"""
		PDB.mask return a selection of (sub) residues for a structure
		@param ffrom: the first (sub) residues of the selection
		@param tto: the last (sub) residues of the selection
		@param mask: (if specified) is a string of length to-from, positions corresponding to '-' will be discarded
		@return: a selection of (sub) residues for a structure
		"""
		res = []
		aPos = 0
		if tto == -1:
			tto = len(self)
		if (mask != "") and (len(mask) < tto-ffrom):
			tto = ffrom + len(mask)
		for i in range(ffrom,tto):
			if mask == "" or ((mask != "") and (mask[aPos] != '-')):
				res = res + self[i].flat()
			aPos = aPos + 1
		if res == []:
			return None
		return PDB(res, id = self.id)


	def zonesOverlap(self,  theZones ):
	    """
	    We test if zones overlaps.
	    """

	    for aZone1 in range( 0, len(theZones)-1, 1 ):
		for aZone2 in range( aZone1+1, len(theZones), 1 ):
		    if ( theZones[ aZone1 ][0] >= theZones[ aZone2 ][0] and theZones[ aZone1 ][0] <= theZones[ aZone2 ][1] ) or ( theZones[ aZone1 ][1] >= theZones[ aZone2 ][0] and theZones[ aZone1 ][1] <= theZones[ aZone2 ][1] ):
			return True

	    return False

	def parseZone(self, znDef, verbose = False):
	    """
	    Parse user defined paired Zones.

	    syntax is : 
		 zone1|zone2|...|zoneN
	    with zoneX    = from:to
	    exemple:
		 12:18|24:30

	    @param znDef   : a string formatted zone
	    @return        : a list of list of 2 indexes defining the ranges, or None is any error.
	    """
	    if verbose:
		    print >> sys.stderr, "parseZone: will parse user zone definition: %s " % znDef
    
	    # Parse the string defining the zones
	    zones = []
	    a = znDef.split("|") # fragmented zone: each subzone
	    for b in a:
		    c = b.split(":")  # reference, mobi sub zones
		    tmp2 = []
		    for d in c:
			    try:
				    tmp2.append( int(d) )
			    except:
				    print >> sys.stderr, "Error: Invalid pattern for zones !"
				    return None
		    zones.append(tmp2)

	    # check lists overlapping
	    if self.zonesOverlap( zones ):
		    print >> sys.stderr, "Error: zones definition must not overlap ! Reference:", zones
		    return None
	    
	    return zones

	def subPDB(self, zoneStr = None, zones = None, verbose = False):
		"""
		PDB.subPDB return a selection of (sub) residues
		@param zoneStr: None, or what residues on the form x:y|a:b where x,y corerspond to one fragment,
		                a,b to another. Numbers from 0, included.
				If specified, zoneStr will supersed the zones.
		@param zones  : a list of list of intervals or residues to consider.
		@return: a PDB, a selection of residues

		Warning: This function does not consider chains or whatever.
		"""
		if zoneStr is not None:
			zones = self.parseZone(zoneStr, verbose = verbose)
		# Select zones if any
		aPDB = self
		if zones is not None:
			oPDB = None
			for aZone in zones:
			    # negative index compatibility
			    if aZone[0] < 0:
				aZone[0] = len(self) + aZone[0]
			    if aZone[1] < 0:
				aZone[1] = len(self) + aZone[1]
			    if oPDB is not None:
				oPDB += self[aZone[0]:aZone[1]+1]
			    else:
				oPDB = self[aZone[0]:aZone[1]+1]
			aPDB = oPDB
			
		return aPDB

        def subPDBFromResList(self, rList, verbose = 0):
		"""
		PDB.subPDBFromResList return a selection of residues in rList.
                @param rList :a list of residues matching the findRes function
                specified on the form _chId_rName_rNum_icode_
                blank icode and blank icode must be encoded as ""

                import PyPDB.PyPDB as PDB
                x = PDB.PDB("test/1crnA.pdb")
                rList = ["_A_THR_1__", "_A_THR_2__", "_A_ILE_7__" ]
                y = x.subPDBFromResList(rList, verbose = True)
                
                """
                if rList == None:
                        return self
                if len(rList) == 0:
                        return self
                res = []
                for r in rList:
                        it = r.split("_")
                        chId = it[1]
                        rName = it[2]
                        rNum  = it[3]
                        icode = it[4]
                        # findRes(self,chId,rName,rNum, icode, what = None, verbose = 0):
                        if verbose:
                                sys.stderr.write("Looking for %s \n"% r)
                        aRes = self.findRes(chId, rName, rNum, icode)
                        if aRes is None:
                                return None
                        res += aRes
                return PDB(res)
                
        def BC(self, rList = None, verbose = 0):
		"""
                PDB.BC(): return barycenter of the PDB or of the PDB subset if rList is not None
                @param rList: a list of residues, specified as for the subPDBFromResList()
                @return : a tuple of either coordinate or None if some residue of the rList was missing
		"""

                if rList is not None:
                        sub = self.subPDBFromResList(rList, verbose = verbose)
                        if not sub:
                                return None, None, None
                else:
                        sub = self
                atomList = atmList(sub.atms)
                return atomList.BC()

        def gridSize(self, rList = None, verbose = 0):
		"""
                PDB.gridSize(): return barycenterthe size of the grid centered on the PDB or of the PDB subset if rList is not None
                @param rList: a list of residues, specified as for the subPDBFromResList()
                @return : a tuple of either coordinate or None if some residue of the rList was missing
		"""
                if rList is not None:
                        sub = self.subPDBFromResList(rList, verbose = verbose)
                        if not sub:
                                return None, None, None
                else:
                        sub = self
                atomList = atmList(sub.atms)
                return atomList.gridSize()

	def header(self):
		"""
		PDB.HEADER
		@return: the title of the file
		"""
		title=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items[1:]:
					if string.count(aItem,"-") == 2:
						break
					if title != '':
						title = title + " "
					title = title + aItem
		return title

	def compound(self):
		"""
		PDB.compound
		@return: the nature of the file (string)
		"""
		title=''
		for Line in self.info:
			if Line[:6]=='COMPND':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title

	def source(self):
		"""
		PDB.source
		@return: where the molecule come from (string)
		"""
		title=''
		for Line in self.info:
			if Line[:6]=='SOURCE':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title

	def author(self):
		"""
		PDB.author
		@return: the author of the structure
		"""
		title=''
		for Line in self.info:
			if Line[:6]=='AUTHOR':
				items = string.split(Line[6:70])
				for aItem in items:
					if title != '':
						title = title + " "
					title = title + aItem
		return title

	def keywords(self):
		"""
		PDB.keywords
		@return: a list of keywords
		"""
		keylist = ''
		for Line in self.info:
			if Line[:6]=='KEYWDS':
				keylist=keylist+Line[10:-1]
		aPos = 0
		OK = 1
		while string.find(keylist,'\'',aPos) != -1:
			aPos = string.find(keylist,'\'',aPos)
			afunc = keylist[0:aPos]+"\\"+keylist[aPos:]
			keylist = afunc
			aPos = aPos + 1
		return keylist

	def date(self):
		"""
		PDB.date
		@return: creation date of the file (a string)
		"""
		date=''
		for Line in self.info:
			if Line[:6]=='HEADER':
				items = string.split(Line)
				for aItem in items:
					if string.count(aItem,"-") == 2:
						date = aItem
				break
		if date != '':
			return date
		# If no creation date, try revision date
		return self.revdate()

	def revdate(self):
		"""
		PDB.revdate (supposes last revision is first REVDAT)
		@return: a string coding for the last revision date.
		"""
		date=''
		for Line in self.info:
			if Line[:6]=='REVDAT':
				date=string.split(Line[13:22])[0]
				break
		return date

	def expmethod(self, verbose = 0):
		"""
		PDB.expmethod
		@return: method by which crds were generated, it corresponds to the values of the EXPDTA field:'X-RAY DIFFRACTION', 'NMR', 'ELECTRON DIFFRACTION', etc.
		"""
		for Line in self.info:
			if Line[:6]=='EXPDTA':
				if string.find(Line,'X-RAY DIFFRACTION')!=-1:
					return 'X-RAY DIFFRACTION'
				if string.find(Line,'X-RAY POWDER DIFFRACTION')!=-1:
					return 'X-RAY POWDER DIFFRACTION'
				elif string.find(Line,'NMR')!=-1:
					return 'NMR'
				elif string.find(Line,'ELECTRON DIFFRACTION')!=-1:
					return 'ELECTRON DIFFRACTION'
				elif string.find(Line,'FIBER DIFFRACTION')!=-1:
					return 'FIBER DIFFRACTION'
				elif string.find(Line,'FLUORESCENCE TRANSFER')!=-1:
					return 'FLUORESCENCE TRANSFER'
				elif string.find(Line,'NEUTRON DIFFRACTION')!=-1:
					return 'NEUTRON DIFFRACTION'
				elif string.find(Line,'THEORETICAL MODEL')!=-1:
					return 'THEORETICAL MODEL'
				elif string.find(Line,'SYNCHROTRON')!=-1:
					return 'SYNCHROTRON'
				elif string.find(Line,'ELECTRON MICROSCOPY')!=-1:
					return 'ELECTRON MICROSCOPY'
				elif string.find(Line,'INFRARED SPECTROSCOPY')!=-1:
					return 'INFRARED SPECTROSCOPY'
                                        
# /!\ attention, not totally pertinent. NMR can have resolution
#				else:
#					# Suppose if resolution set: Xray
#					if self.resolution() != -1.:
#						return 'X-RAY DIFFRACTION'
#					return ''
#		# Suppose if resolution set: Xray
#		if self.resolution() != -1.:
#			return 'X-RAY DIFFRACTION'
		return ''

	def resolution(self, verbose = 0):
		"""
		PDB.resolution()
		@return: the Resolution of the file, if specified somewhere (return -1 if not found). (data mining in the REMARK lines)
		@note: This is only available for files determined using Xray.
		"""
		resol = -1.
		for Line in self.info:
			if string.find(Line,'REMARK   2 RESOLUTION')!=-1:
				posMax=string.find(Line,'ANGSTROM')-1
				posMin=string.find(Line,'RESOLUTION')+11
				if posMax!=-1:
					try:
						resol=float(Line[posMin:posMax])
					except ValueError:
						pass
		return resol

	def rvalue(self, verbose = 0):
		"""
		PDB.rvalue()
		A corresponding method is defined for the free R value:
		@return: the RValue of the file (float), if specified somewhere (-1 if not). (data mining in the REMARK lines)
		"""
		R_VALUE = "NULL"
		checkRValue = 0
		for Line in self.info:
			if string.find(Line,'REMARK   3') != -1:
				# Case where it is on the next line !!
				if R_VALUE == "NULL" and ((checkRValue == 1) or (checkRValue == 2)):
					if checkRValue == 1:
						if string.find(Line,'.') != -1:
							pos=string.find(Line,'.')-1
							checkRValue == 0
							try:
								R_VALUE=float(Line[pos:pos+5])
							except ValueError:
								R_VALUE   = "NULL"
					elif checkRValue == 2:
						startPos = string.find(Line,'VALUE')
						if string.find(Line,'.', startPos) != -1:
							pos=string.find(Line,'.', startPos)-1
							toPos = pos+5
							# check for cases such as: 0.20.
							if string.count(Line,'.', pos,toPos) > 1:
								toPos = string.find(Line,'.', pos+2)
							try:
								R_VALUE=float(Line[pos:toPos])
							except ValueError:
								R_VALUE   = "NULL"
					checkRValue = 0
				# On one line ?
				# 2009 mod to match WORKING + TEST SET
				# if R_VALUE == "NULL" and (string.find(Line,' R ') != -1 or string.find(Line,'R VALUE') != -1 or string.find(Line,'R-VALUE') != -1 or string.find(Line,'R-FACTOR') != -1) and string.find(Line,'TEST') == -1 and string.find(Line,'FREE') == -1 and string.find(Line,'ESTIMATE') == -1 and string.find(Line,'BIN') == -1 and  string.find(Line,'ERROR') == -1:
				if R_VALUE == "NULL" and (string.find(Line,' R ') != -1 or string.find(Line,'R VALUE') != -1 or string.find(Line,'R-VALUE') != -1 or string.find(Line,'R-FACTOR') != -1) and string.find(Line,'FREE') == -1 and string.find(Line,'ESTIMATE') == -1 and string.find(Line,'BIN') == -1 and  string.find(Line,'ERROR') == -1:
					startPos = string.find(Line,'R VALUE')
					if startPos == -1:
						startPos = string.find(Line,'R-VALUE')
					if startPos == -1:
						startPos = string.find(Line,'R-FACTOR')
					if startPos == -1:
						if string.find(Line,' R '):
							checkRValue = 2
					if verbose:
						print >> sys.stdout, "#", Line[:-1]
						print >> sys.stdout, "#", Line[startPos:-1]
					if string.find(Line,'.', startPos) != -1:
						pos=string.find(Line,'.', startPos)-1
						toPos = pos+5
						# check for cases such as: 0.20.
						if string.count(Line,'.', pos,toPos) > 1:
							toPos = string.find(Line,'.', pos+2)
						try:
							R_VALUE=float(Line[pos:toPos])
						except ValueError:
							if Line[pos] == 'O':
								try:
									R_VALUE=float(Line[pos+1:toPos])
								except ValueError:
									R_VALUE   = "NULL"
							else:
								R_VALUE   = "NULL"
					else:
						checkRValue = 1
		return R_VALUE

	def freervalue(self):
		"""
		PDB.freervalue() look for the Rvalue
		@return: the rvalue or NULL if not found
		"""
		FREE_R_VALUE   = "NULL"
		for Line in self.info:
			if string.find(Line,'FREE R VALUE') != -1 and string.find(Line,'TEST') == -1 and string.find(Line,'ESTIMATE') == -1 and string.find(Line,'BIN') == -1 and  string.find(Line,'ERROR') == -1:
				if string.find(Line,'.') != -1:
					pos=string.find(Line,'.')-1
					try:
						FREE_R_VALUE=float(Line[pos:pos+5])
					except ValueError:
						FREE_R_VALUE   = "NULL"
		return FREE_R_VALUE


        def pH(self, verbose = 0):
                """
                PDB.pH()
                @return: PH at which structure was solved (if specified)
                """
                rs = "NULL"
		for Line in self.info:
			if string.find(Line,' PH        ') != -1:
                                it = Line.split()
                                rs = it[-1]
                                break
                return rs
                                
        def isMembrane(self, verbose = 0):
                """
                PDB.isMembrane()
                @return: True if the information contains anything indicating the protein is related to membrane
                """
                rs = False
		for Line in self.info:
			if string.find(Line.upper(),'MEMBRAN') != -1:
                                rs = True
                                break
                return rs
                                
	def seqresaa3(self, chIds=None, verbose = 0):
		"""
		PDB.seqresaa3()
		@param chIds: the chains Ids
		@return: the sequence of the PDB as specified in the SEQRES lines.
		"""
		if not chIds:
			chIds = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES':
					if string.count(chIds,Line[11]) == 0:
						chIds = chIds + Line[11]
		if len(chIds) > 1:
			rs = []
		for chId in chIds:
			aseqres = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES' and Line[11] == chId:
					aseqres = aseqres + Line[19:70]+' '
			if len(chIds) > 1:
				rs.append(aseqres.split())
			else:
				rs = aseqres.split()
		return rs

	def seqres(self,chIds=None, verbose = 0):
		"""
		PDB.seqres()
		@param chIds: It is possible to specify chain Id
		@return:  a string of the sequence in the SEQRES lines. If several chains exist: a list of all the sequences is returned.
		@note: I{example:}
		x.seqres(``ABD'') will return a list of the three sequences (if exist) corresponding to the chains A, B and D.
		"""
		if not chIds:
			chIds = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES':
					if string.count(chIds,Line[11]) == 0:
						chIds = chIds + Line[11]
		if len(chIds) > 1:
			rs = []
		for chId in chIds:
			aseqres = ''
			for Line in self.seq:
				if Line[:6]=='SEQRES' and Line[11] == chId:
					aseqres = aseqres + Line[19:70]+' '
			type = self.chnType(chId)
			if type == 'Protein':
				if verbose:
					sys.stderr.write("seqres for %s\n" % aseqres)
				aa1seq=SEQREStoAA1(aseqres, verbose = verbose)
			elif type == "DNA" or type == "RNA":
				curseqres = string.split(aseqres)
				aa1seq = ""
				for i in curseqres:
					aa1seq = aa1seq + i[0]
			else:
				aa1seq = aseqres
			if len(chIds) > 1:
				rs.append(aa1seq)
			else:
				rs = aa1seq
		return rs

	def CAonly(self,verbose=0):
		"""
		PDB.CAonly()
		@return: Does the file contain only CAs ? (Yes/No)
		"""
		res=True
		for aLine in self.data:
			if string.find(aLine[12:15],"CA")==-1:
				res=False
				if verbose:
					print >> sys.stdout, '# PDB do not contain only CA atoms'
				return res
				break
		if verbose:
			print >> sys.stdout, '# PDB contain only CA atoms'
		return res

	def SCatmMiss(self, verbose = 0):
		"""
		PDB.SCAtmMiss()
		@return:
			- nSCMiss: the number of amino-acid residues having some side chain missing atom
			- SCatmMiss:  a string containing the information about all the residues with missing side chain atoms? (Yes/No)
		@note: For each residue, the string (SCatmMiss) consists of RName_ChLbl_RNum_RIcode, where RName is the name of the residue (3 letters), ChLbl is the chain label, RNum is the number of the residue (as in the PDB file) and RIcode the PDB insertion code of the residue.
		@note: For residues having at least one atomic coordinate present
		"""
		SCatmMiss=""
		nSCMiss = 0
		for i in range(0,len(self)):
			resName = self[i].rName()
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue
			aaTpe = AA3STRICT.index(resName)
			chaine = ""
			for atm in self[i].atms:
				chaine=chaine+atm.atmName()+' '
			if verbose:
				print >> sys.stdout, "#", chaine
			missp = 0
			for atms in AASC[aaTpe]:
				if string.find(chaine,atms)==-1:
					missp = 1
					break
			if missp:
				nSCMiss = nSCMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				SCatmMiss = SCatmMiss+Res
		return nSCMiss, SCatmMiss

	def BBatmMiss(self, verbose = 0):
		"""
		PDB.BBAtmMiss()
		@return:
			- nBBmiss: the number of amino-acid residues having some backbone missing atom (one of N, CA, C, O)
			- BBatmMiss: nd a string concatening the information about all the residues with missing peptidic chain atoms
		@note: For each residue, the string (BBatmMiss) consists of RName_ChLbl_RNum_RIcode, where RName is the name of the residue (3 letters), ChLbl is the chain label, RNum is the number of the residue (as in the PDB file) and RIcode the PDB insertion code of the residue.
		@note: For residues having at least one atomic coordinate present
		"""
		BBatmMiss=""
		status = False
		nBBMiss = 0
		for i in range(0,len(self)):
			resName = self[i].rName()
			if AA3STRICT.count(resName) == 0:
			## Suppose nonstandard amino-acids are OK
				continue
			aaTpe = AA3STRICT.index(resName)
			theCheck = self[i].BBAtmMiss()
			missp = 0
			if theCheck != []:
				missp = 1
			if (missp == 1) and (theCheck[0] == "O") and (self[i].atmPos("OXT") != None):
				missp = 0
##			chaine = ""
##			for atm in self[i].atms:
##				## chaine=chaine+string.split(str(atm))[2]+' '
##				chaine=chaine+atm.atmName()+' '
##			if verbose:
##				print chaine
##			missp = 0
##			for atms in AABB:
##				if string.find(chaine,atms)==-1:
##					missp = 1
##					break
			if missp:
				status = True
				if i == 0:
					status = False
				if i == len(self) -1 and not status:
					status = False
				if i > 0 and i < len(self) -1:
					status = True
				nBBMiss = nBBMiss+1
				resName = self[i].rName()
				resNum = self[i].rNum()
				icode  = self[i].riCode()
				lbl  = self[i].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				BBatmMiss = BBatmMiss+Res
		return nBBMiss, BBatmMiss


        def BBclean(self, verbose = 0):
                """
import PyPDB.PyPDB as PDB
x = PDB.PDB("1crnA.pdb")
rs = x.BBclean(verbose = True)
                """

                nBBMiss, BBatmMiss = self.BBatmMiss(verbose = verbose)
                if verbose:
                        sys.stderr.write("BBclean: residues with missing BB atoms: %s\n" % BBatmMiss)
                if not nBBMiss:
                        return self
                oPDB = None
                for i, aRes in enumerate(self):
                        isOK = True
                        for rString in BBatmMiss.split():
                                if aRes.isResCode(rString, verbose = verbose):
                                   isOK = False
                        if isOK:
                                if not oPDB:
                                        oPDB = self[i:i+1]
                                else:
                                        oPDB += self[i:i+1]
                return oPDB


	def hasAltAtms(self,verbose = 0):
		"""
		PDB.hasAltAtms
		@return: Does the file has BBaltAtm or SCAltAtm? (Yes/No for each)
		"""
		BBAltAtm = False
		SCAltAtm = False
		for i in self:
			BB, SC = i.hasAltAtms()
			if BB:
				BBAltAtm = True
			if SC:
				SCAltAtm = True
		return BBAltAtm, SCAltAtm

	def altAtmsResList(self,verbose = 0):
		"""
		PDB.altAtmsResList
		@return: nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm:
			- nBBAltAtm: number of Back Bones in alternate atoms
			- BBAltAtm: the Back Bones in alternate atoms
			- nSCAltAtm: number of side chain in alternate atoms
			- SCAltAtm: side chain in alternate atoms
		"""
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		for i in self:
			BB, SC = i.hasAltAtms()
			if not BB and not SC:
				continue
			resName = i.rName()
			resNum = i.rNum()
			icode  = i.riCode()
			lbl  = i.chnLbl()
			if icode == ' ':
				icode = ''
			if lbl == ' ':
				lbl = ''
			resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
			if BB:
				nBBAltAtm = nBBAltAtm + 1
				BBAltAtm = BBAltAtm + resLabel
			if SC:
				nSCAltAtm = nSCAltAtm + 1
				SCAltAtm = SCAltAtm + resLabel
		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm

	def geomCheck(self,verbose=0):
		"""
		PDB.geomCheck() This will scan and check that the peptidic bonds geometry is rather correct. It is based on the value of the peptidic bond.
		@return: Is the BB peptidic geometry (distance) correct? (OK/Poor/Bad)
		@note: THIS WILL NOT DETECT FRAGMENTS. IF MANY, THE GAPS ARE IGNORED AND DO NOT RESULT IN "Bad" RETURN. \n
		This allows to scan that all the fragments are correct at once.
		"""
		aN = None
		aC = None
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		for aRes in self:
			if AA3.count(aRes.rName()) == 0:
				continue
			aN = aRes.atmPos("N")
			if aN != None:
				Nx, Ny, Nz = aRes[aN].xyz()
				theN = aRes[aN]
			if aC != None:
				if theN.chnLbl() == theC.chnLbl():
					aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
					if aDist > 1.50 and aDist < 3.:
						if verbose:
							print >> sys.stdout, "# Poor peptidic bond of ",aDist," for ", theC.resName(), theC.resNum(), theN.resName(), theN.resNum()
						if BBGeoOK == "Ok":
							BBGeoOK = "Poor"
					elif aDist > 3.:
						if verbose:
							print >> sys.stdout, "# Bad peptidic bond  of ",aDist," for :", theC.resName(), theC.resNum(), theN.resName(), theN.resNum()
						BBGeoOK = "Bad"
			aC  = aRes.atmPos("C")
			if aC != None:
				Cx, Cy, Cz = aRes[aC].xyz()
				theC = aRes[aC]
		return BBGeoOK

	def traceCheck(self,hetSkip = 0, maxCADist = 4.2, verbose = 0):
		"""
		PDB.traceCheck check if BB peptidic geometry is correct (distance)
		@param maxCADist: the maximum distance between 2 CA consecutive, 4.2 angstrom by default
		@return: traceOK (OK/bad), tracePB (residues with bad geometry), nCISPRO (number of cis prolines), CISPRO (cis prolines), nCISPep (number of cis peptides), CISPep (cis peptides)
		"""
		theTrace = self.select(awhat=["CA"])
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Ok"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		tracePB = ""
		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1][0].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:\n# %s" % theTrace[aRes - 1]
				return CisWarning,"No"
			try:
				x2, y2, z2 = theTrace[aRes][0].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:\n# %s" % theTrace[aRes]
				return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			if aDist < 3.60: # CIS peptide
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if resName != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + resLabel
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + resLabel
			if aDist > maxCADist: # bad geometry
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					print >> sys.stdout, "# Bad Trace for ",theTrace[aRes-1]," dist = ",aDist," / ",maxCADist
		return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep

	def traceCheck2(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		"""
		PDB.traceCheck2 check if BB peptidic geometry is correct (distance)
		@param minCADist: the minimum distance between 2 CA consecutive, 3.7 angstrom by default
		@param maxCADist: the maximum distance between 2 CA consecutive, 3.9 angstrom by default
		@return: traceOK (OK/bad), tracePB (residues with bad geometry), nCISPRO (number of cis prolines), CISPRO (cis prolines), nCISPep (number of cis peptides), CISPep (cis peptides)
		"""
		theTrace = self.select(awhat=["CA"])
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Ok"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		tracePB = ""
		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1][0].xyz()
				raise ValueError
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:\n# %s" % theTrace[aRes - 1]
				return CisWarning,"No"
			try:
				x2, y2, z2 = theTrace[aRes][0].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:\n# %s" % theTrace[aRes]
				return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			if aDist < 3.3: # CIS peptide
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if resName != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + resLabel
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + resLabel
				if verbose:
					print >> sys.stdout, "# %s : CIS Trace for %s dist = %f" % (self.id, theTrace[aRes-1],aDist)
			if (aDist < minCADist) and (aDist >= 3.3): # bad geometry
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					print >> sys.stdout, "# %s : Bad Trace for %s dist = %f / %f" % (self.id, theTrace[aRes-1],aDist,minCADist)
			if aDist > maxCADist: # bad geometry
				resName = theTrace[aRes].rName()
				resNum = theTrace[aRes].rNum()
				icode  = theTrace[aRes].riCode()
				lbl  = theTrace[aRes].chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				tracePB  = tracePB + resLabel
				traceOK = "Bad"
				if verbose:
					print >> sys.stdout, "# %s : Bad Trace for %s dist = %f / %f" % (self.id, theTrace[aRes-1],aDist,maxCADist)
		return traceOK, tracePB, nCISPRO, CISPRO, nCISPep, CISPep

	def resLabel(self,aRes):
		"""
		PDB.resLabel()
		@param aRes: a residue
		@return: the label of the residue
		"""
		resName = self[aRes].rName()
		resNum = self[aRes].rNum()
		icode  = self[aRes].riCode()
		lbl  = self[aRes].chnLbl()
		if icode == ' ':
			icode = ''
		if lbl == ' ':
			lbl = ''
		resLabel = resName+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
		return resLabel

	def traceCheck3(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		"""
		PDB.traceCheck3 check if BB peptidic geometry is correct (distance)
		@param minCADist: the minimum distance between 2 CA consecutive, 3.7 angstrom by default
		@param maxCADist: the maximum distance between 2 CA consecutive, 3.9 angstrom by default
		@return: traceOK (OK/bad), tracePB (residues with bad geometry), nCISPRO (number of cis prolines), nCISPep (number of cis peptides)
		"""
		traceOK   = "Ok"
		tracePB   = ""
		nCISPRO   = 0
		nCISPEP   = 0
		for aRes in range(0,len(self)-1):
			CAPos = self[aRes].atms.CApos()
			if CAPos == None:
				traceOK = "No"
				tracePB = tracePB + self.resLabel(aRes)
				continue
			try:
				x1, y1, z1 = self[aRes][CAPos].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:", self[aRes].rName(),self[aRes].rNum()
				traceOK = "No"
				tracePB = tracePB + self.resLabel(aRes)
				continue
			CAPos1 = self[aRes+1].atms.CApos()
			if CAPos1 == None:
				traceOK = "No"
				continue
			try:
				x2, y2, z2 = self[aRes+1][CAPos1].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: fname incorrect ATOM format for:", self[aRes+1].rName(),self[aRes+1].rNum()
				traceOK = "No"
				continue
			aDist = distance(x1, y1, z1, x2, y2, z2)
			if aDist > maxCADist: # bad geometry
				resLabel = self.resLabel(aRes)
				tracePB  = tracePB + resLabel
				if traceOK != "No":
					traceOK = "Bad"
				if verbose:
					print >> sys.stdout, "# %s : Bad Trace for %s dist = %f / %f\n" % (self.id, resLabel,aDist,maxCADist)
			if aDist < minCADist: # bad geometry
				# We check the case for CIS
				try:
					CA  = self[aRes].atms.theAtm("CA")
					C   = self[aRes].atms.theAtm("C")
					x1,y1,z1 = CA.xyz()
					x2,y2,z2 = C.xyz()
				except:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
				try:
					N   = self[aRes+1].atms.theAtm("N")
					CA2 = self[aRes+1].atms.theAtm("CA")
					x3,y4,z3 = N.xyz()
					x4,y3,z4 = CA2.xyz()
				except:
					resLabel = self.resLabel(aRes+1)
					tracePB  = tracePB + resLabel
				ome = dihedral(x1,y1,z1,x2,y2,z2,x3,y4,z3,x4,y3,z4)
				if abs(ome) < 20:
					resName = self[aRes+1].rName()
					resLabel = self.resLabel(aRes)
					if resName != "PRO":
						nCISPEP += 1
						if verbose:
							print >> sys.stdout, "# %s : CISPEP conformation for %s ome = %f" % (self.id, resLabel, ome)
					else:
						print >> sys.stdout, "# %s : CISPRO conformation for %s ome = %f" % (self.id, resLabel, ome)
						nCISPRO += 1
				else:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
					if traceOK != "No":
						traceOK = "Bad"
					if verbose:
						print >> sys.stdout, "# %s : Bad Trace for %s dist = %f / %f" % (self.id, resLabel,aDist,maxCADist)
		return traceOK, tracePB, nCISPEP, nCISPRO

	def CISSeq(self,hetSkip = 0, minCADist = 3.7, maxCADist = 3.9, verbose = 0):
		"""
		CISSeq
		@param minCADist: the minimum distance between 2 CA consecutive, 3.7 angstrom by default
		@param maxCADist: the maximum distance between 2 CA consecutive, 3.9 angstrom by default
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@return: the prolines and other peptides found to be in the cis conformation.
		"""
		oSeq = ""
		for aRes in range(0,len(self)-1):
			CAPos = self[aRes].atms.CApos()
			if CAPos == None:
				continue
				#raise ValueError, "Incorrect ATOM line format for: %s %s" % (self[aRes].rName(),self[aRes].rNum())
			try:
				x1, y1, z1 = self[aRes][CAPos].xyz()
			except ValueError:
				continue
				#raise ValueError, "Incorrect ATOM line format for: %s %s" % (self[aRes].rName(),self[aRes].rNum())

			CAPos1 = self[aRes+1].atms.CApos()
			if CAPos1 == None:
				continue
				#raise ValueError, "Incorrect ATOM line format for: %s %s" % (self[aRes+1].rName(),self[aRes+1].rNum())
			try:
				x2, y2, z2 = self[aRes+1][CAPos1].xyz()
			except ValueError,AttributeError:
				continue
				#raise ValueError, "Incorrect ATOM line format for: %s %s" % (self[aRes+1].rName(),self[aRes+1].rNum())
			aDist = distance(x1, y1, z1, x2, y2, z2)
			if aDist < minCADist: # bad geometry
				# We check the case for CIS
				try:
					CA  = self[aRes].atms.theAtm("CA")
					C   = self[aRes].atms.theAtm("C")
					x1,y1,z1 = CA.xyz()
					x2,y2,z2 = C.xyz()
				except:
					resLabel = self.resLabel(aRes)
					tracePB  = tracePB + resLabel
				try:
					N   = self[aRes+1].atms.theAtm("N")
					CA2 = self[aRes+1].atms.theAtm("CA")
					x3,y4,z3 = N.xyz()
					x4,y3,z4 = CA2.xyz()
				except:
					resLabel = self.resLabel(aRes+1)
					tracePB  = tracePB + resLabel
				ome = dihedral(x1,y1,z1,x2,y2,z2,x3,y4,z3,x4,y3,z4)
				if abs(ome) < 20:
					resName = self[aRes+1].rName()
					resLabel = self.resLabel(aRes)
					if resName != "PRO":
						oSeq += "2"
						if verbose:
							print >> sys.stderr, "# %s : CISPEP conformation for %s ome = %f\n" % (self.id, resLabel, ome)
					else:
						oSeq += "1"
				else:
					oSeq += "0"
			else:
				oSeq += "0"
		return oSeq

	def chnCAFrgList(self, chId = "", maxDist = 4.10):
		"""
		PDB.chnCAFrgList determine fragments based on alpha carbon inter-atomic distance alone
		@param chId: a chain ID
		@param maxDist: the  maximal distance between two consecutive AC to be in the same fragment (default = 4.10)
		@return: the chain with fragments separated and the number of fragments.
		"""
		if chId == "" and len(self.chnList()) > 1:
			print >> sys.stderr, "# PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\""
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0
		for aRes in range(1,len(self)):
			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
				# skip non amino acid
				continue
			aC = self[aRes-1].atmPos("CA")
			if aC == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1
			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			aN = self[aRes].atmPos("CA")
			if aN == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def asOneChn(self,chnId = ' '):
		"""
		PDB.asOneChn
		@param chnId: a chain ID
		@return: a PDB instance with residues renumbered as if there were only one chain.
		"""
		for aRes in range(0,len(self)):
			self[aRes].chnLbl(chnId)
			self[aRes].rNum(aRes+1)
		return PDB(self.flat())

	def renumber(self, ffrom, verbose = 0):
		"""
		PDB.renumber renumber residues
		@author: P. Tuffery
		@param ffrom: renumber residues from
		@return: none
		"""
		for res in self.rt:
			res.rNum(ffrom)
			ffrom += 1
		self.atmTab()

	def resRenumber(self,pattern='',verbose=0):
		"""
		PDB.resRenumber Renumber residus
		@author: F.Briand
		@param pattern: "chnName1(string):ffrom1(integer):tto1(integer):index1(integer) chnName2(string):ffrom2(integer):tto2(integer):index2(integer) ..."
		@return: none
		@note: edit the PDB class "in place"
		"""
		# pattern parsing
		if pattern=='':
			raise ValueError, "PDB.resRenumber() : pattern missing. USAGE : \"chnName1(string):ffrom1(integer):tto1(integer):index1(integer) chnName2(string):ffrom2(integer):tto2(integer):index2(integer) ...\""
		for subPattern in pattern.split():
			try:
				chnName, ffrom, tto, index = subPattern.split(":")
				chnName, ffrom, tto, index = str(chnName), int(ffrom), int(tto), int(index)
			except:
				raise ValueError, "PDB.resRenumber() : wrong pattern. USAGE : \"chnName1(string):ffrom1(integer):tto1(integer):index1(integer) chnName2(string):ffrom2(integer):tto2(integer):index2(integer) ...\""
			for iChn in chnName:
				if iChn not in self.chnList():
					raise ValueError, "PDB.resRenumber() : \"%s\" is not a chain. Chains : %s" % (iChn, self.chnList())
				if verbose:
					print >> sys.stdout, "# Renumbering residues from %d to %d of chain %s with id starting to %d" % (ffrom, tto, chnName, index)
				# renumbering
				delta = None
				for res in self.rt:
					resNum = int(res.rNum())
					if res.chnLbl() == iChn and resNum >= ffrom and resNum < tto:
						if not delta:
							# setting "delta" between old number and new number
							delta = index - resNum
						# setting new res number, based on delta
						newindex = resNum + delta
						try:
							# check the "format" of the new number (<10000)
							if not int(math.log10(newindex))+1 > 4:
								res.rNum(newindex)
							else:
								print >> sys.stderr, "# PDB.resRenumber() : indexation reach max value for residue sequence number. Error on residue\""+res.rName(),res.rNum()+"\"("+str(newindex),"> 9999)"
								return
						except ValueError:
							print >> sys.stderr, "# PDB.resRenumber() : indexation under zero. Error on residue\"",res.rName(),res.rNum(),"(",res.chnLbl(),")\" (HETATM atom could have \"out of range\" number.). Residue ignored."
							pass
		self.atmTab()

	def atmsRenumber(self, pattern='', verbose = 0):
		"""
		PDB.atmsRenumber Renumber atoms
		@author: F.Briand
		@param pattern: index(integer):ffrom(integer):tto(integer)
		@return: True, False if interrupt (indexation problem)
		@note: edit the PDB class "in place"
		"""
		# pattern parsing
		if not pattern:
			raise ValueError, "PDB.atmsRenumber() : pattern missing. USAGE : \"index(integer):ffrom(integer):tto(integer)\""
		index, ffrom, tto = pattern.split(":")
		index, ffrom, tto = int(index), int(ffrom), int(tto)
		if verbose:
			print >> sys.stdout, "# Renumbering atoms from",ffrom,"to",tto,"with id starting to",index
		delta = index - ffrom
		if verbose:
			print >> sys.stdout, "# Renumbering atoms ..."
		# ATOM renumbering
		for atom in self.atms:
			if int(atom.atmNum()) >= ffrom and int(atom.atmNum()) < tto:
				index = int(atom.atmNum()) + delta
				# check the "format" of the new number (<100000)
				if not int(math.log10(index))+1 > 5:
					atom.atmNum(index)
				else:
					print >> sys.stderr, "# PDB.atmsRenumber() : indexation reach max value for atom serial number. Error on residue \""+atom.resName(),atom.resNum()+"\" ("+str(index),"> 99999)"
					return
		# CONECT renumbering
		if verbose:
			print >> sys.stdout, "# Renumbering CONECT fields ..."
		for iConect in range(0,len(self.conect)):
			if verbose:
				print >> sys.stdout, self.conect[iConect],
			for iAtom in range(0,14):
				try:
					# parsing conect fields
					atom = int(self.conect[iConect][6+(iAtom*5):11+(iAtom*5)])
					if atom >= ffrom and atom <= tto:
						index = atom + delta
						# check the "format" of the new number (<100000)
						if not int(math.log10(index))+1 > 5:
							self.conect[iConect] = "%s%5d%s" % (self.conect[iConect][:6+(iAtom*5)],index,self.conect[iConect][11+(iAtom*5):])
						else:
							print >> sys.stderr, "# PDB.atmsRenumber() : indexation reach max value for atom serial number. Error on connect\""+self.conect[iConect].split()[1]+"\"("+str(index),"> 99999)"
							return
				except:
					pass
			if verbose:
				print >> sys.stdout, "#", self.conect[iConect]
		self.resTab(0)

	def atmsRenumber2(self, ffrom=1, verbose = 0):
		"""
		PDB.atmsRenumber2 Renumber atoms
		@author: P. Tuffery
		@param ffrom: number of the first atom
		@return: True, False if interrupt (indexation problem)
		@note: edit the PDB class "in place"
		"""
		# pattern parsing
		if verbose:
			print >> sys.stdout, "# Renumbering atoms from",ffrom
		if verbose:
			print >> sys.stdout, "# Renumbering atoms ..."
		# ATOM renumbering
		index = ffrom
		for atom in self.atms:
			atom.atmNum(index)
			index += 1
		# CONECT renumbering
		if verbose:
			print >> sys.stdout, "# Renumbering CONECT fields ..."
		for iConect in range(0,len(self.conect)):
			if verbose:
				print >> sys.stdout, self.conect[iConect],
			for iAtom in range(0,14):
				try:
					# parsing conect fields
					atom = int(self.conect[iConect][6+(iAtom*5):11+(iAtom*5)])
					if atom >= ffrom and atom <= tto:
						index = atom + delta
						# check the "format" of the new number (<100000)
						if not int(math.log10(index))+1 > 5:
							self.conect[iConect] = "%s%5d%s" % (self.conect[iConect][:6+(iAtom*5)],index,self.conect[iConect][11+(iAtom*5):])
						else:
							print >> sys.stderr, "# PDB.atmsRenumber() : indexation reach max value for atom serial number. Error on connect\""+self.conect[iConect].split()[1]+"\"("+str(index),"> 99999)"
							return
				except:
					pass
			if verbose:
				print >> sys.stdout, "#", self.conect[iConect]
		self.resTab(0)

	def resMeanBVal(self, pattern='', verbose = 0):
		"""
		PDB.resMeanBVal Renumber the B Value with the mean on each residue
		@author: F.Briand
		@param pattern: chn1:ffrom1:tto1_chn2:ffrom2:tto2_...
		@return: Nothing, edit "in place"
		"""
		if not pattern:
			pattern = ""
			chns = self.chnList()
			for aChn in chns:
				ffrom =  self.rt[0].rNum()
				tto = self.rt[-1].rNum()
				pattern += aChn+":"+ffrom+":"+tto+"_"
				
		print "THE PATTERN", pattern
		
		for iPattern in pattern.split("_")[:-1]:
			chns, ffrom, tto = iPattern.split(":")
			ffrom, tto = int(ffrom), int(tto)
			if verbose:
				print >> sys.stdout, "# Average B Value on chain(s) %s, on residues from %d to %d" % (chns,ffrom,tto)
			bValDict = {}
			for res in self.rt:
				if res.chnLbl() in chns and int(res.rNum()) >= ffrom and int(res.rNum()) <= tto:
					meanBVal = 0.
					for atom in res:
						meanBVal += float(atom.atmBVal())
					if verbose:
						print >> sys.stdout, "# Residue : %s => B Value : %6.2f" % (res.chnLbl() + res.rNum(), meanBVal / len(res))
					meanBVal = meanBVal / len(res)
					for atom in res:
						try:
							atom.atmBVal(meanBVal)
						except KeyError:
							pass
		self.atmTab()

	def atmsBValRenumber(self, input, verbose = 0):
		"""
		PDB.atmsBValRenumber renumber B Values according to the pattern in input
		@author: F.Briand
		@param input: "newBValue" for all the atoms or pairs "AtomNumber newBValue". 1 BValue per line, separated by blank (tabs, spaces ...)
		@return: Nothing, edit "in place"
		"""
		try:
			# opening input file
			if verbose:
				print >> sys.stdout, "# Trying: %s" % input
			input_file=simpleload(input, verbose=verbose)
			if verbose:
				print >> sys.stdout, "# succeeded!"
		except IOError:
			raise IOError, "%s : not a local file ... or some problem occurred ..." % input
		atmsList = {}
		# check the "type" of file : 2 column(if) or 1 column(else)
		if len(input_file[0].split())>1:
			for line in input_file:
				try:
					result = line.split()
					atmsList[result[0]] = float(result[1])
				except (ValueError, AttributeError):
					print >> sys.stderr, "# Wrong line format :\n# %s" % line
			if verbose:
				print >> sys.stdout, "#", atmsList
			for atom in self.atms:
				try:
					atom.atmBVal(atmsList[atom.atmNum()])
					if verbose:
						print >> sys.stdout, "# %s => BValue : %6.2f" % (atom.atmNum(),atmsList[atom.atmNum()])
				except KeyError:
					pass
		else:
			i = 0
			imax = len(input_file)
			for atom in self.atms:
				try:
					atom.atmBVal(float(input_file[i]))
					if verbose:
						print >> sys.stdout, "# %s => BValue : %6.2f" % (atom.atmNum(),float(atom.atmBVal()))
				except ValueError:
					raise ValueError, "Wrong line format :\n# %s" % input_file[line]
				i+=1
				if i>=imax:
					break
		self.resTab(0)

	def contactMap(self, what = "CA", dist = 11.):
		"""
		PDB.contactMap: return an internal  contact map based on
		CA distance
		@param what : the atom name (CA, CB)
		@param dist : the threshold to detect contacts
		@return     : an array of the contactMap
		"""
		nCtct = 0
		rs = []
		for i in range(len(self)):
			rs.append([0] * len(self))
		for i in range(len(self)):
			iwhat = what
			if (iwhat == "CB") and (self[i].rName() == "GLY"):
				iwhat = "CA"
			x1, y1, z1 = self[i][iwhat].xyz()
			for j in range(i+1, len(self)):
				jwhat = what
				if (what == "CB") and (self[j].rName() == "GLY"):
					jwhat = "CA"
				x2, y2, z2 = self[j][jwhat].xyz()
				# print i, j, self[i].rName(), iwhat, self[j].rName(), jwhat
				aDist = distance(x1, y1, z1, x2, y2, z2)
				if aDist < 12:
					rs[i][j] = 1
					rs[j][i] = 1
					nCtct += 1
		return rs, nCtct

	def chnFrgList(self, chId = "", maxDist = 1.7):
		"""
		PDB.chnFrgList determine fragments based on inter-atomic distance C'-N alone
		@param chId: a chain ID
		@param maxDist: the  maximal C'-N distance to be in the same fragment
		@return: the chain with fragments separated and the number of fragments.
		@note: 1.70 is default threshold
		"""
		if chId == "" and len(self.chnList()) > 1:
			print >> sys.stderr, "# PDB.chnFrgList() : cannot handle several chains as \""+self.chnList()+"\""
			return []
		res = []
		oriRes = 0
		lastRes = 0
		nFrg = 0
		for aRes in range(1,len(self)):
			try:
				aaTpe = AA3.index(self[aRes-1].rName())
			except ValueError:
			# skip non amino acid
				continue
			aC = self[aRes-1].atmPos("C")
			if aC == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = self[aRes-1][aC].xyz()
			CchnLbl = self[aRes-1].chnLbl()
			lastRes = aRes-1
			try:
				aaTpe = AA3.index(self[aRes].rName())
			except ValueError:
				continue
			aN = self[aRes].atmPos("N")
			if aN == None:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx, Ny, Nz = self[aRes][aN].xyz()
			NchnLbl = self[aRes].chnLbl()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			lastRes = aRes
			if aDist > maxDist or (CchnLbl != NchnLbl):
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		# lRes.append(len(self) - 1)
		lRes.append(lastRes)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def frgList(self, maxNCDist = 1.7, maxCADist = 4.1, verbose = 0):
		"""
		PDB.frgList will return the number of fragments and their boundaries
		@param maxNCDist: the maximum distance between N and C to be in the same fragment
		@param maxCADist: the maximum distance between 2 consecutive AC to be in the same fragment
		@return: a list of the fragments of the PDB if some geometric inconsistencies are
		detected and the numbers of fragments
		"""
		# (3, [[0, 326], [327, 655], [656, 860]])
		res = []
		oriRes = 0
		nFrg = 0
		chnIds = self.chnList()
		curDex = 0
		for chId in chnIds:
			curChn = self.chn(chId)
			if self.chnType(chId) != "Protein":
				curDex = curDex + len(curChn)
				continue
			CAonly = curChn.CAonly()
			if not CAonly:
				curNFrg, curFrgList = curChn.chnFrgList(maxDist=maxNCDist)
			else:
				curNFrg, curFrgList = curChn.chnCAFrgList(maxDist=maxCADist)
				# curNFrg = 1
				# curFrgList = [[0,len(curChn)-1]]
			for i in range(0,len(curFrgList)):
				curFrgList[i][0] = curFrgList[i][0] + curDex
				curFrgList[i][1] = curFrgList[i][1] + curDex
				res.append(curFrgList[i])
			nFrg = nFrg + curNFrg
			curDex = curDex + len(curChn)
		return nFrg, res

	def nFrg(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		"""
		nFrg
		@param maxNCDist: the maximum distance between N and C to be in the same fragment
		@param maxCADist: the maximum distance between 2 consecutive AC to be in the same fragment
		@return: the numbers of fragments detected based on NC and CA distances
		"""
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)
		return nFrg

	# This will not check for fragments
	def aaseq_ori(self, verbose = 0):
		"""
		PDB.aaseq_ori()
		@return: the sequence of residues present in the PDB file, having coordinates.
		@note: Converts non standard amino-acids to equivalent standard amino-acid.
		"""
		res = ""
		unres = []
		for aRes in self:
			if AA3STRICT.count(aRes.rName()):
				res = res + AA1[AA3STRICT.index(aRes.rName())]
			elif AA3.count(aRes.rName()):
				rName = aRes.rName()
				if verbose:
					print >> sys.stdout, "# Unfrequent residue type : %s" % rName
				if rName == "MSE":   # seleno MET
					res = res+"M"
				elif rName == "CSE": # seleno CYS
					res = res+"C"
				elif rName == "FGL": # amino propane dioique
					res = res+"S"
				elif rName == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif rName == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif rName == "TRO": # HYDROXY TRYPTOPHANE
					res = res+"W"
				elif rName == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif rName == "MHO": # Hydroxy-MET
					res = res+"M"
				elif rName == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif rName == "HYP": # HYDROXY PRO
					res = res+"P"
				elif rName == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				elif rName == "AYA": # acetyl ALA
					res = res+"A"
				elif rName == "SEG": # hydroxy ALA
					res = res+"A"
				elif rName == "MME": # N methyl MET
					res = res+"M"
##				elif rName == "BET": # 3methyl GLY
##					res = res+"G"
##				elif rName == "DAR": # D ARG
##					res = res+"R"
				elif rName == "FME": # formyl MET
					res = res+"M"
				elif rName == "CXM": # carboxy MET
					res = res+"M"
				elif rName == "SAC": # acetyl SER
					res = res+"S"
				elif rName == "CSO": # -HYDROXYCYSTEINE
					res = res+"C"
				elif rName == "HTR": # BETA-HYDROXYTRYPTOPHANE
					res = res+"W"
				else:
					res = res+'X'
			else:
				if unres.count(aRes.rName()) == 0:
					unres.append(aRes.rName())
				# res = res+'X'
		if verbose:
			print >> sys.stderr, "# Unknown residue type (2):",
			for iRes in unres:
				print >> sys.stderr, " %s " % iRes,
				print >> sys.stderr
		return res

	# This will not check for fragments
	def aaseq(self, matchAtms = None, hetSubst = True, verbose = 0):
		"""
		PDB.aaseq()
		@param matchAtms: a list of atoms searched
		@param hetSubst: converts non standards to their standard equivalent
		@return: a list of atoms matching
		"""
		res = ""
		unres = []
		for aRes in self:
			if matchAtms != None:
				lOK = 1
				for aAtm in matchAtms:
					if verbose:
						print >> sys.stdout, "# aaseq checking atmName %s" % aAtm
					if not aRes.findAtm(atmName = aAtm, verbose = verbose):
						lOK = 0
						break
				if not lOK:
					continue
			if AA3STRICT.count(aRes.rName()):
				res = res + AA1[AA3STRICT.index(aRes.rName())]
			# elif AA3.count(aRes.rName()):
			elif AA3new.count(aRes.rName()):
				rName = aRes.rName()
				if verbose:
					print >> sys.stderr, "Unfrequent residue type : %s" % rName
					
				# dico instead of too much elif
				if hetSubst and (dico_AA.has_key(rName)):
					res = res+dico_AA[rName]
				else:
					res = res+'X'
			else:
				if unres.count(aRes.rName()) == 0:
					unres.append(aRes.rName())
                                res = res+'X' # TUFFERY 2013
		if verbose and len(unres):
			print >> sys.stderr, "# Unknown residue type (2):",
			for iRes in unres:
				print >> sys.stderr, " %s " % iRes,
				print >> sys.stderr
		return res

	def frgseq(self, maxNCDist = 1.7, maxCADist = 4.1, verbose=0):
		"""
		PDB.frgseq fragments the sequence according to maxNCDist and maxCADist
		@param maxNCDist: the maximum distance between N and C to be in the same fragment
		@param maxCADist: the maximum distance between 2 consecutive AC to be in the same fragment
		@return: the fragments of the PDB if some geometric inconsistencies are
		detected and the numbers of fragments
		"""
		res = []
		nFrg, frgList = self.frgList(maxNCDist,maxCADist,verbose)
		if verbose:
			print >> sys.stdout, "# frgSeq :", nFrg, frgList
		for i in frgList:
			res.append( self[i[0]:i[1]+1].aaseq())
		return res

	def SGList(self):
		"""
		PDB.SGList
		@return: a list of all the coordinates of the gamma sulfur
		"""
		SGList = []
		for aRes in self:
			if aRes.rName() == "CYS":
				lSGList = []
				for aAtm in aRes.atms:
					if aAtm.atmName() == "SG":
						lSGList.append(aAtm.xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList

	def nSSIntra(self, maxDist = 2.35):
		"""
		nSSIntra()
		@return: nSSbond the number of SSbond in a protein
		"""
		nSSBond = 0
		aSGList = self.SGList()
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < maxDist:
							nSSBond = nSSBond + 1
							break
		return nSSBond

	def SSIntra(self, maxDist = 2.35):
		"""
		SSIntra()
		@return: the position of the cysteins involed in a SSbond
		"""
		SSBonds = []
		for aRes1 in range(0,len(self)):
			if self[aRes1].rName() != "CYS":
				continue
			for aAtm in self[aRes1].atms:
				if aAtm.atmName() == "SG":
					xyz1 = aAtm.xyz()
			for aRes2 in range(aRes1+1,len(self)):
				if self[aRes2].rName() != "CYS":
					continue
				for aAtm in self[aRes2].atms:
					if aAtm.atmName() == "SG":
						xyz2 = aAtm.xyz()
				if apply(distance, xyz1+xyz2) < maxDist:
					SSBonds.append([aRes1, aRes2])
		return SSBonds

	def isHalfCys(self, aRes, maxDist = 2.35):
		"""
		isHalfCys(aRes) checks if the distance between the sulfur of two cysteins allows a disulfide bond
		@param aRes: the number of the residue, must be a "CYS"
		@return: the position of the second cys and the distance separating them
		"""
		if self[aRes].rName() != "CYS":
			return 0,0,0
		x = 0.
		y = 0.
		z = 0.
		isSet = 0
		for aAtm in range(0,len(self[aRes].atms)):
			if self[aRes].atms[aAtm].atmName() == "SG":
				x,y,z = self[aRes].atms[aAtm].xyz()
				isSet = 1
		if isSet == 0:
			return 0,0,0
		for aPos in range(0,len(self)):
			if self[aPos].rName() != "CYS":
				continue
			if aPos == aRes:
				continue
			for aAtm in range(0,len(self[aPos].atms)):
				if self[aPos].atms[aAtm].atmName() == "SG":
					x1,y1,z1 = self[aPos].atms[aAtm].xyz()
					if distance(x,y,z,x1,y1,z1) < maxDist:
						return 1, aPos, distance(x,y,z,x1,y1,z1)
		return 0,0,0

	def findRes(self,chId,rName,rNum, icode, what = None, verbose = 0):
		"""
		PDB.findRes To identify a residue given its chain Id, name, PDB number, insertion code
		@return:
			- either the residue if what == None
			- or the residue rank (from 0) if what != None
		"""
		aPos = -1
		for aRes in self:
			aPos = aPos + 1
			if verbose:
				print >> sys.stdout, "# \"%s\" \"%s\" \"%s\"" % (aRes.chnLbl(),aRes.rName(), aRes.rNum())
			if chId != "" and chId != None:
				if aRes.chnLbl() != chId:
					continue
			if rName != "" and rName != None:
				if aRes.rName() != rName:
					continue
			if rNum != "" and rNum != None:
				if aRes.rNum() != rNum:
					continue
			if icode != "" and icode != None:
				if aRes.riCode() != icode:
					continue
			if what != None:
				return aPos
			return aRes
		return None

	def findAtm(self,chId,rName,rNum, icode, atmName = "CA", verbose = 0):
		"""
		PDB.findAtm To identify an atom given residue chain Id, name, PDB number, insertion code and atom Name
		@param chId: the chain Id
		@param rName: the residu name
		@param rNum: the residu number
		@param icode: the line code number
		@param atmName: the atom name
		@param verbose: if verbose=1 it will print the aAtm loop or None if not found
		@return: either the atom instance
		"""
		res = self.findRes(chId,rName,rNum, icode)
		if res != None:
			for aAtm in res.atms:
				if aAtm.atmName() == atmName:
					return aAtm
		return None

	#########################################################################
	#																		#
	# clean_ori : Modification de la structure PDB :						#
	#																		#
	#	+ suppression des residus PCA										#
	#	+ transformation des residus MSE en MET								#
	#	- transformation des atomes SE en S									#
	#	+ transformation des residus CSE en CYS								#
	#	- transformation des atomes SE en S									#
	#	+ transformation des residus CEA en CYS								#
	#	- suppression des atomes O1 et HO1									#
	#	+ transformation des residus CGU en GLU								#
	#	- suppression des atomes CD2, OE3, OE4, HE4							#
	#	+ transformation des residus HTR en TRP								#
	#	- suppression des atomes O et OH									#
	#	+ transformation des residus TPQ en PHE								#
	#	- suppression des atomes O2, O4, O5 et HO4							#
	#																		#
	#########################################################################
	def clean_ori(self, whatRes = ["PCA", "5HP", "FGL", "MSE", "CSE", "CEA", "CGU", "HTR", "MHO", "IAS", "TPQ", "TYS", "AYA", "FME", "CXM", "SAC", "CSO", "MME", "SEG", "HYP", "TRO"], verbose = 0 ):
		"""
		clean(whatRes = ["5HP", "PCA","MSE", "CSE", "CEA", "CGU", "HTR", "TPQ"])
		cleanup PDB files by converting some non standard residue into standard ones.
		@param whatRes: it is the list of residues that may be affected.
		@type whatRes: the default is all. whatRes could be only part of the default list.
		@note:
		+ transformation des residus PCA en GLU
		+ transformation des residus MHO en MET
		+ transformation des residus IAS en ASP
		+ transformation des residus HYP en PRO
		+ transformation des residus TPQ en TYR
		+ transformation des residus TRO en TRP
		+ transformation des residus TYS en TYR
		+ transformation des residus MSE en MET
		   - transformation des atomes SE en S
		+ transformation des residus CSE en CYS
		   - transformation des atomes SE en S
		+ transformation des residus CEA en CYS
		   - suppression des atomes O1 et HO1
		+ transformation des residus CGU en GLU
		   - suppression des atomes CD2, OE3, OE4, HE4
		+ transformation des residus HTR en TRP
		   - suppression des atomes O et OH
		+ transformation des residus TPQ en PHE
		   - suppression des atomes O2, O4, O5 et HO4
		+ transformation des residus FGL en SER
		   - suppression des atomes OG1, renomme OG1 en OG
		+ transformation des residus AYA (acetyl ALA) en ALA
		   - suppression des atomes CT, OT, CM
		+ transformation des residus FME (formyl MET) en MET
		   - suppression des atomes OF, CF, (also CN, O1, HCN)
		+ transformation des residus CXM (carboxy MET) en MET
		   - suppression des atomes CN, O1, O2, HO1, HO2
		+  transformation des residus SAC (acetyl SER) en SER
		   - suppression des atomes C1A, C2A, OAC, 1H2A, 2H2A, 3H2A
		+  transformation des residus CSO (s-hydroxycysteine) en CYS
		   - suppression des atomes OD, HD
		+  transformation des residus BET (3methyl GLY) en GLY (NOT BY DEFAULT)
		   - suppression des atomes C1, C2, C3, 1H1, 1H2, 1H3, 2H1, 2H2, 2H3, 3H1, 3H2, 3H3
		+  transformation des residus DAR (D ARG) en ARG (NOT BY DEFAULT)
		+  transformation des residus MME (n-methyl MET) en MET
		   - suppression des atomes CM, 1HM, 2HM, 3HM
		+  transformation des residus SEG (HYDROXYALANINE) en ALA
		   - suppression des atomes OH, HOD
		"""
		#   To do:
		#   FGL->CYS ??, ERROR !
		#   TPQ->TYR *done* (not PHE),
		#   MHO->MET *done*,
		#   IAS->ASP *done: BUT requires atom addition!
		#   CGU->GLU *done*,
		#   TYS->TYR *done*
		# for iR in range(len(self)-1,-1,-1):
		# residu = self[iR]
		if verbose:
			print >> sys.stdout, "# Searching non standard residues ..."
		for residu in self:
			if residu.rName() not in whatRes:
				continue
			if verbose:
				print >> sys.stdout, "#", residu.rName()
			if residu.rName()=="PCA":
				# self.__delitem__(iR)
				residu.rName("GLU")
				residu.delete(["CD","OE"])
			if residu.rName()=="5HP":
				residu.rName("GLU")
				residu.delete(["CD", "OD"])
			if residu.rName()=="FGL":
				residu.rName("SER")
				residu.delete(["OG2"])
				for atom in residu:
					if atom.atmName()=="OG1":
						atom.atmName("OG")
			if residu.rName()=="MSE":
				residu.rName("MET")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CSE":
				residu.rName("CYS")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CEA":
				residu.rName("CYS")
				residu.delete(["O1","HO1"])
			if residu.rName()=="CGU":
				residu.rName("GLU")
				residu.delete(["CD2","OE3","OE4","HE4"])
			if residu.rName()=="HTR":
				residu.rName("TRP")
				residu.delete(["O","HO"])
			if residu.rName()=="MHO":
				residu.rName("MET")
				residu.delete(["OD1"])
			if residu.rName()=="AYA":
				residu.rName("ALA")
				residu.delete(["CT", "OT", "CM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="HYP":
				residu.rName("PRO")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="SEG":
				residu.rName("ALA")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="DAR":
				residu.rName("ARG")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="FME":
				residu.rName("MET")
				residu.delete(["OF", "CF", "CN", "O1", "HCN"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="MME":
				residu.rName("MET")
				residu.delete(["CM", "1HM", "2HM", "3HM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CXM":
				residu.rName("MET")
				residu.delete(["CN", "O1", "O2", "HO1", "HO2"])
				for atom in residu:
					atom.header("ATOM  ")
					atom.header("ATOM  ")
			if residu.rName()=="SAC":
				residu.rName("SER")
				residu.delete(["C1A", "C2A", "OAC", "1H2A", "2H2A", "3H2A"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CSO":
				residu.rName("CYS")
				residu.delete(["OD", "HD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="BET":
				residu.rName("GLY")
				residu.delete(["C1", "C2", "C3", "1H1", "1H2", "1H3", "2H1", "2H2", "2H3", "3H1", "3H2", "3H3"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="TRO":
				residu.rName("TRP")
				residu.delete(["OD1", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="IAS":
				residu.rName("ASP")
			if residu.rName()=="TPQ":
				residu.rName("TYR")
				residu.delete(["O2","O5"])
				for atom in residu:
					if atom.atmName()=="O4":
						atom.atmName("OH")
					if atom.atmName()=="HO4":
						atom.atmName("HH")
			if residu.rName()=="TYS":
				residu.rName("TYR")
				residu.delete(["S","O1","O2","O3","HO3"])
			for atom in residu:
				atom.header("ATOM  ") # 6 chars for header
		return None

	def clean(self, whatRes = AA3new, verbose = 0 ):
		"""
		clean(whatRes = ["5HP", "PCA","MSE", "CSE", "CEA", "CGU", "HTR", "TPQ"])
		@param whatRes: it is the list of residues that may be affected.
		@type whatRes: the default is all. whatRes could be only part of the default list.
		@note:
		cleanup PDB files by converting some non standard residue into standard ones.
		+ transformation des residus PCA en GLU
		+ transformation des residus MHO en MET
		+ transformation des residus IAS en ASP
		+ transformation des residus HYP en PRO
		+ transformation des residus TPQ en TYR
		+ transformation des residus TRO en TRP
		+ transformation des residus TYS en TYR
		+ transformation des residus MSE en MET
		   - transformation des atomes SE en S
		+ transformation des residus CSE en CYS
		   - transformation des atomes SE en S
		+ transformation des residus CEA en CYS
		   - suppression des atomes O1 et HO1
		+ transformation des residus CGU en GLU
		   - suppression des atomes CD2, OE3, OE4, HE4
		+ transformation des residus HTR en TRP
		   - suppression des atomes O et OH
		+ transformation des residus TPQ en PHE
		   - suppression des atomes O2, O4, O5 et HO4
		+ transformation des residus FGL en SER
		   - suppression des atomes OG1, renomme OG1 en OG
		+ transformation des residus AYA (acetyl ALA) en ALA
		   - suppression des atomes CT, OT, CM
		+ transformation des residus FME (formyl MET) en MET
		   - suppression des atomes OF, CF, (also CN, O1, HCN)
		+ transformation des residus CXM (carboxy MET) en MET
		   - suppression des atomes CN, O1, O2, HO1, HO2
		+  transformation des residus SAC (acetyl SER) en SER
		   - suppression des atomes C1A, C2A, OAC, 1H2A, 2H2A, 3H2A
		+  transformation des residus CSO (s-hydroxycysteine) en CYS
		   - suppression des atomes OD, HD
		+  transformation des residus BET (3methyl GLY) en GLY (NOT BY DEFAULT)
		   - suppression des atomes C1, C2, C3, 1H1, 1H2, 1H3, 2H1, 2H2, 2H3, 3H1, 3H2, 3H3
		+  transformation des residus DAR (D ARG) en ARG (NOT BY DEFAULT)
		+  transformation des residus MME (n-methyl MET) en MET
		   - suppression des atomes CM, 1HM, 2HM, 3HM
		+  transformation des residus SEG (HYDROXYALANINE) en ALA
		   - suppression des atomes OH, HOD
		"""
		#   To do:
		#   FGL->CYS ??, ERROR !
		#   TPQ->TYR *done* (not PHE),
		#   MHO->MET *done*,
		#   IAS->ASP *done: BUT requires atom addition!
		#   CGU->GLU *done*,
		#   TYS->TYR *done*
		# for iR in range(len(self)-1,-1,-1):
		# residu = self[iR]
		if verbose:
			print >> sys.stdout, "# Searching non standard residues"
		for residu in self:
			if residu.rName() not in whatRes:
				continue
			if verbose:
				print >> sys.stdout, "#", residu.rName()
			if residu.rName()=="HSE":
				# self.__delitem__(iR)
				residu.rName("HIS")
			if residu.rName()=="HSD":
				# self.__delitem__(iR)
				residu.rName("HIS")
			if residu.rName()=="HSP":
				# self.__delitem__(iR)
				residu.rName("HIS")
			if residu.rName()=="PCA":
				# self.__delitem__(iR)
				residu.rName("GLU")
				residu.delete(["CD","OE"])
			if residu.rName()=="5HP":
				residu.rName("GLU")
				residu.delete(["CD", "OD"])
			if residu.rName()=="FGL":
				residu.rName("SER")
				residu.delete(["OG2"])
				for atom in residu:
					if atom.atmName()=="OG1":
						atom.atmName("OG")
			if residu.rName()=="MSE":
				residu.rName("MET")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CSE":
				residu.rName("CYS")
				for atom in residu:
					if atom.atmName()=="SE":
						atom.atmName("S")
			if residu.rName()=="CEA":
				residu.rName("CYS")
				residu.delete(["O1","HO1"])
			if residu.rName()=="CGU":
				residu.rName("GLU")
				residu.delete(["CD2","OE3","OE4","HE4"])
			if residu.rName()=="HTR":
				residu.rName("TRP")
				residu.delete(["O","HO"])
			if residu.rName()=="MHO":
				residu.rName("MET")
				residu.delete(["OD1"])
			if residu.rName()=="AYA":
				residu.rName("ALA")
				residu.delete(["CT", "OT", "CM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="HYP":
				residu.rName("PRO")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="SEG":
				residu.rName("ALA")
				residu.delete(["OD", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="DAR":
				residu.rName("ARG")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="FME":
				residu.rName("MET")
				residu.delete(["OF", "CF", "CN", "O1", "HCN"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="MME":
				residu.rName("MET")
				residu.delete(["CM", "1HM", "2HM", "3HM"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CXM":
				residu.rName("MET")
				residu.delete(["CN", "O1", "O2", "HO1", "HO2"])
				for atom in residu:
					atom.header("ATOM  ")
					atom.header("ATOM  ")
			if residu.rName()=="SAC":
				residu.rName("SER")
				residu.delete(["C1A", "C2A", "OAC", "1H2A", "2H2A", "3H2A"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="CSO":
				residu.rName("CYS")
				residu.delete(["OD", "HD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="BET":
				residu.rName("GLY")
				residu.delete(["C1", "C2", "C3", "1H1", "1H2", "1H3", "2H1", "2H2", "2H3", "3H1", "3H2", "3H3"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="TRO":
				residu.rName("TRP")
				residu.delete(["OD1", "HOD"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="IAS":
				residu.rName("ASP")
			if residu.rName()=="TPQ":
				residu.rName("TYR")
				residu.delete(["O2","O5"])
				for atom in residu:
					if atom.atmName()=="O4":
						atom.atmName("OH")
					if atom.atmName()=="HO4":
						atom.atmName("HH")
			if residu.rName()=="TYS":
				residu.rName("TYR")
				residu.delete(["S","O1","O2","O3","HO3"])
				for atom in residu:
					atom.header("ATOM  ")
			# add other residus here :
			if residu.rName()=="OMT":
				residu.rName("MET")
				residu.delete(["OD1","OD2"])
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="ACL":
				residu.rName("ARG")
				residu.delete(["CM","1HM","2HM"])
				for atom in residu:
					if atom.atmName()=="CL":
						atom.atmName("3HM")
					atom.header("ATOM  ")
			if residu.rName()=="AGM":
				residu.rName("ARG")
				residu.delete(["1HE2","2HE2","3HE2"])
				for atom in residu:
					if atom.atmName()=="CE2":
						atom.atmName("HD")
					atom.header("ATOM  ")
			if residu.rName()=="ARM":
				residu.rName("ARG")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					if atom.atmName()=="CM": #replace CM by O
						atom.atmName("O")
					if atom.atmName()=="1HM": #replace CL by H
						atom.atmName("H")
					atom.header("ATOM  ")
			if residu.rName()=="HAR":
				residu.rName("ARG")
				for atom in residu:
					if atom.atmName()=="OH1":
						atom.atmName("HH")
					atom.header("ATOM  ")
			if residu.rName()=="HMR":
				residu.rName("ARG")
				residu.delete(["2HC","OXT","HXT","2HC"])
				for atom in residu:
					if atom.atmName()=="1HC":
						atom.atmName("O")
					if atom.atmName()=="C":
						atom.atmName("OH")
					if atom.atmName()=="O":
						atom.atmName("HO")
					atom.header("ATOM  ")
			if residu.rName()=="AIB":
				residu.rName("ALA")
				residu.delete(["1HB2","2HB2","3HB2"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="CB2":
						atom.atmName("HA")
			if residu.rName()=="ALM":
				residu.rName("ALA")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="O":
						atom.atmName("O1")
					if atom.atmName()=="CM":
						atom.atmName("O2")
					if atom.atmName()=="1HM":
						atom.atmName("HO2")
			# HETNAM PHE mais suppose ALA --> a verifier
			if residu.rName()=="BNN":
				residu.rName("ALA")
				residu.delete(["O1","CH","N16","N17","C1","C2","C3","C4","C5","C6",
				"C15","1HH1","2HH1","3HH1","1H16","2H16","H17","H2","H3","H5","H6"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C7":
						atom.atmName("CB")
					if atom.atmName()=="C1":
						atom.atmName("3HB")
					if atom.atmName()=="2H7":
						atom.atmName("1HB")
					if atom.atmName()=="1H7":
						atom.atmName("2HB")
					if atom.atmName()=="C11":
						atom.atmName("2HN")
					if atom.atmName()=="H":
						atom.atmName("1HN")
			# EN GLY plutot que ALA
			if residu.rName()=="CHG":
				residu.rName("GLY")
				residu.delete(["H1","C2","C3","C4","C5","C6",
				"1H2","2H2","1H3","2H3","1H4","2H4","1H5","2H5","1H6","2H6"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C7":
						atom.atmName("CA")
					if atom.atmName()=="H7":
						atom.atmName("1HA")
					if atom.atmName()=="N7":
						atom.atmName("N")
					if atom.atmName()=="1HN7":
						atom.atmName("1HN")
					if atom.atmName()=="2HN7":
						atom.atmName("2HN")
					if atom.atmName()=="C1":
						atom.atmName("2HA")
			if residu.rName()=="CSD":
				residu.rName("ALA")
				residu.delete(["OD1","OD2","HD2"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="SG":
						atom.atmName("3HB")
			if residu.rName()=="DAL":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
			# 2 H missing
			if residu.rName()=="DHA":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="DNP":
				residu.rName("ALA")
				residu.delete(["1HG","2HG","3HG"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="NG":
						atom.atmName("3HB")
			if residu.rName()=="FLA":
				residu.rName("ALA")
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="F1":
						atom.atmName("1HB")
					if atom.atmName()=="F2":
						atom.atmName("2HB")
					if atom.atmName()=="F3":
						atom.atmName("3HB")
			if residu.rName()=="HAC":
				residu.rName("ALA")
				residu.delete(["C2","C3","C4","C5","C6","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","1H6","2H6","H1"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C1":
						atom.atmName("3HB")
			if residu.rName()=="MAA":
				residu.rName("ALA")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="CM":
						atom.atmName("2HN")
			if residu.rName()=="PRR":
				residu.rName("ALA")
				residu.delete(["C4","H4","C3","C9","H9","N1","1H10","2H10","3H10",
				"C10","C2","H2","1H3"])
				for atom in residu:
					atom.header("ATOM  ")
					if atom.atmName()=="C9":
						atom.atmName("3H5")
			# one H missing on OG1
			if residu.rName()=="ALO":
				residu.rName("THR")
				for atom in residu:
					atom.header("ATOM  ")
			if residu.rName()=="BMT":
				residu.rName("THR")
				residu.delete(["1HH","2HH","3HH","CH","CZ","HZ","HE","CE",
				"1HD2","2HD2","2HD1","1HD1","3HD1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD2":
						atom.atmName("2HG2")
					if atom.atmName()=="CD1":
						atom.atmName("3HG2")
			if residu.rName()=="DTH":
				residu.rName("THR")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="TPO":
				residu.rName("THR")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HG1")
			if residu.rName()=="BCS":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","CE","CZ1","CZ2","HZ1","HZ2",
				"CT1","CT2","HT1","HT2","CH"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
			if residu.rName()=="BUC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3","1H4","2H4","3H4",
				"C4","c3","C2","C1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			if residu.rName()=="C5C":
				residu.rName("CYS")
				residu.delete(["H1","1H2","2H2","1H3","2H3","1H4","2H4",
				"1H5","2H5","C1","C2","C3","C4","C5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			if residu.rName()=="C6C":
				residu.rName("CYS")
				residu.delete(["H1","1H2","2H2","1H3","2H3","1H4","2H4",
				"1H5","2H5","1H6","2H6","C1","C2","C3","C4","C5","C6"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			if residu.rName()=="CCS":
				residu.rName("CYS")
				residu.delete(["HOZ","OZ1","OZ2","CE","1HD","2HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
			if residu.rName()=="CME":
				residu.rName("CYS")
				residu.delete(["HO","OH","CZ","1HZ","2HZ","CE","1HE","2HE"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			# Verifier la valence du residu csp
			if residu.rName()=="CSP":
				residu.rName("CYS")
				residu.delete(["O1P","O2P","O3P","PHO2","PHO3","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("H")
			if residu.rName()=="CSS":
				residu.rName("CYS")
				residu.delete(["HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			if residu.rName()=="CSW":
				residu.rName("CYS")
				residu.delete(["OD1","OD2"])
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="CSX":
				residu.rName("CYS")
				residu.delete(["OD"])
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="CY3":
				residu.rName("CYS")
				residu.delete(["2H1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1H1":
						atom.atmName("HXT")
					if atom.atmName()=="N1":
						atom.atmName("OXT")
			if residu.rName()=="CYG":
				residu.rName("CYS")
				residu.delete(["OE2","CG1","1HG1","2HG1","CB1","1HB1","2HB1",
				"CA1","HA1","N1","1HN1","2HN1","C1","O1","O2","HO2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD1":
						atom.atmName("H")
			if residu.rName()=="CYM":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","3HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
			if residu.rName()=="DCY":
				residu.rName("CYS")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="EFC":
				residu.rName("CYS")
				residu.delete(["F2","1H2","2H2","C2","C1","1H1","2H1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			if residu.rName()=="OCS":
				residu.rName("CYS")
				residu.delete(["HD2","OD1","OD3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OD2":
						atom.atmName("H")
			if residu.rName()=="PEC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","3H5","C1","C2","C3","C4","C5",])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("H")
			# one H missing on C-ter
			if residu.rName()=="PEC":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","1HZ","2HZ","1HH","2HH","3HH",
				"CE","CZ","CH"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="H":
						atom.atmName("O")
					if atom.atmName()=="SD":
						atom.atmName("HS")
			if residu.rName()=="SCH":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("HS")
			# one H missing on C-ter
			if residu.rName()=="SCS":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE","1HZ","2HZ","CZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SD":
						atom.atmName("HS")
					if atom.atmName()=="H":
						atom.atmName("O")
			if residu.rName()=="SCY":
				residu.rName("CYS")
				residu.delete(["1HE","2HE","3HE","CE","OCD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("HS")
			if residu.rName()=="SHC":
				residu.rName("CYS")
				residu.delete(["1H1","2H1","1H2","2H2","1H3","2H3",
				"1H4","2H4","1H5","2H5","1H6","2H6","3H6","C6","C2","C3",
				"C4","C5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1":
						atom.atmName("H")
			if residu.rName()=="SMC":
				residu.rName("CYS")
				residu.delete(["1HCS","2HCS","3HCS"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CS":
						atom.atmName("H")
			if residu.rName()=="SOC":
				residu.rName("CYS")
				residu.delete(["OD2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="SE":
						atom.atmName("SG")
					if atom.atmName()=="OD1":
						atom.atmName("H")
			if residu.rName()=="ALY":
				residu.rName("LYS")
				residu.delete(["1HH3","2HH3","3HH3","CH3","HO"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("H")
			if residu.rName()=="DLY":
				residu.rName("LYS")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="LLP":
				residu.rName("LYS")
				residu.delete(["1H4A","2H4A","C4","C5","C5A","1H5A","2H5A",
				"04P","P","O1P","O2P","O3P","2HOP","3HOP","C6","N1","H6","C2",
				"C2A","1H2A","2H2A","3H2A","C3","O3","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C4A":
						atom.atmName("H")
			if residu.rName()=="LLY":
				residu.rName("LYS")
				residu.delete(["C1","O1","O2","HO2","C2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("H")
			if residu.rName()=="LYM":
				residu.rName("LYS")
				residu.delete(["3HM","2HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH":
						atom.atmName("OXT")
					if atom.atmName()=="1HM":
						atom.atmName("HXT")
			if residu.rName()=="LYZ":
				residu.rName("LYS")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OH":
						atom.atmName("2HD")
			# Possiblite d'un ASP ?
			if residu.rName()=="SHR":
				residu.rName("LYS")
				residu.delete(["C1","C2","C3","C5","O1","O2","HO1",
				"1H2","2H2","1H3","2H3","O3","O4","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C4":
						atom.atmName("2HN")
			if residu.rName()=="TRG":
				residu.rName("LYS")
				residu.delete(["1HH1","2HH1","3HH1","1HH2","2HH2","3HH2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CH1":
						atom.atmName("1HZ")
					if atom.atmName()=="CH2":
						atom.atmName("2HZ")
			#Valine plutot que Leucine
			if residu.rName()=="BUG":
				residu.rName("VAL")
				residu.delete(["2HN2","1HG3","2HG3","3HG3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="NA":
						atom.atmName("OXT")
					if atom.atmName()=="1HN2":
						atom.atmName("OXT")
					if atom.atmName()=="CG3":
						atom.atmName("H")
			if residu.rName()=="CLE":
				residu.rName("LEU")
				residu.delete(["2H2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="N2":
						atom.atmName("OXT")
					if atom.atmName()=="1H2":
						atom.atmName("OXT")
			if residu.rName()=="DLE":
				residu.rName("LEU")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="MLE":
				residu.rName("LEU")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
			# Transforme en MET plutot que LEU
			if residu.rName()=="NLE":
				residu.rName("MET")
				residu.delete(["1HD","2HD"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
			# Transforme en MET plutot que LEU
			if residu.rName()=="NLN":
				residu.rName("MET")
				residu.delete(["1HD","2HD","2HH2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
					if atom.atmName()=="NH2":
						atom.atmName("OXT")
					if atom.atmName()=="1HH2":
						atom.atmName("HXT")
			# Transforme en MET plutot que LEU
			if residu.rName()=="NLP":
				residu.rName("MET")
				residu.delete(["1HD","2HD","O3","HO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("S")
					if atom.atmName()=="P":
						atom.atmName("C")
			if residu.rName()=="CYQ":
				residu.rName("CYS")
				residu.delete(["1HD","2HD","P","O1P","O2P","O3P","2HOP","3HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CD":
						atom.atmName("H")
			if residu.rName()=="DVA":
				residu.rName("VAL")
				for atom in residu:
					atom.header("ATOM   ")
			# Trasnforme en ALA plutot que VAL
			if residu.rName()=="DIV":
				residu.rName("VAL")
				residu.delete(["1HB2","2HB2","3HB2","1HG1","2HG1","3HG1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CB2":
						atom.atmName("HA")
					if atom.atmName()=="CG1":
						atom.atmName("3HB1")
			if residu.rName()=="MVA":
				residu.rName("VAL")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
			if residu.rName()=="2AS":
				residu.rName("ASP")
				residu.delete(["1HBB","2HBB","3HBB"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CBB":
						atom.atmName("2HB")
			# one H missing on C-ter
			if residu.rName()=="ASA":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HXT":
						atom.atmName("OXT")
			if residu.rName()=="ASB":
				residu.rName("ASP")
				residu.delete(["O1","O2","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1":
						atom.atmName("HOD")
			if residu.rName()=="ASK":
				residu.rName("ASP")
				residu.delete(["1HM","2HM","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("OXT")
					if atom.atmName()=="3HM":
						atom.atmName("HXT")
			if residu.rName()=="ASL":
				residu.rName("ASP")
				residu.delete(["1HC3","2HC3","3HC3","C3","HC2","C1","O1","O2","HO1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C2":
						atom.atmName("HOD")
			if residu.rName()=="ASQ":
				residu.rName("ASP")
				residu.delete(["O1P","O2P","O3P","2HOP","3HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOD")
			if residu.rName()=="BHD":
				residu.rName("ASP")
				residu.delete(["HOB"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OB":
						atom.atmName("2HB")
			if residu.rName()=="DAS":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DSP":
				residu.rName("ASP")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DSN":
				residu.rName("SER")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="MIS":
				residu.rName("SER")
				residu.delete(["C1","C2","C3","1H2","2H2","3H2","H1",
				"1H3","2H3","3H3","1HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOG")
			if residu.rName()=="OAS":
				residu.rName("SER")
				residu.delete(["1HC2","2HC2","3HC2","C2A","OAC"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C1A":
						atom.atmName("HOG")
			# one H missing on C-ter
			if residu.rName()=="SEL":
				residu.rName("SER")
				residu.delete(["2HB2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1HB2":
						atom.atmName("OXT")
			if residu.rName()=="SEP":
				residu.rName("SER")
				residu.delete(["2HOP","3HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HOG")
			if residu.rName()=="SET":
				residu.rName("SER")
				residu.delete(["2HNT"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="1HNT":
						atom.atmName("HXT")
					if atom.atmName()=="NT":
						atom.atmName("OXT")
			if residu.rName()=="SVA":
				residu.rName("SER")
				residu.delete(["O1","O2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="V":
						atom.atmName("HOG")
			if residu.rName()=="DGL":
				residu.rName("GLU")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="GGL":
				residu.rName("GLU")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="GMA":
				residu.rName("GLU")
				residu.delete(["2HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="N2":
						atom.atmName("OXT")
					if atom.atmName()=="1HN":
						atom.atmName("HXT")
			if residu.rName()=="DIL":
				residu.rName("ILE")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="IIL":
				residu.rName("ILE")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="GL3":
				residu.rName("GLY")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("OXT")
					if atom.atmName()=="HS":
						atom.atmName("HXT")
			# one H missing on C-ter
			if residu.rName()=="GLZ":
				residu.rName("GLY")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HXT":
						atom.atmName("OXT")
			if residu.rName()=="GSC":
				residu.rName("GLY")
				residu.delete(["1H1","2H1","1H2","2H2","3H2","C1","C2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("2HA")
			if residu.rName()=="MPQ":
				residu.rName("GLY")
				residu.delete(["1HM","2HM","3HM","CD1","CD2","CE1","CE2","CZ",
				"1HD1","1HE1","1HE2","1HD2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("2HN")
					if atom.atmName()=="CG":
						atom.atmName("2HA")
			if residu.rName()=="MSA":
				residu.rName("GLY")
				residu.delete(["1HN","3HN","1HG","2HG","3HG"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
					if atom.atmName()=="SB":
						atom.atmName("2HA")
			if residu.rName()=="NMC":
				residu.rName("GLY")
				residu.delete(["1HCN","2HCN","CX1","CX2","CX3","1HC2",
				"2HC2","1HC3","2HC3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
			if residu.rName()=="SAR":
				residu.rName("GLY")
				residu.delete(["1HN","2HN","3HN"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CN":
						atom.atmName("2HN")
			if residu.rName()=="DGL":
				residu.rName("GLN")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DHI":
				residu.rName("HIS")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DPN":
				residu.rName("PHE")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DPR":
				residu.rName("PRO")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DTR":
				residu.rName("TRP")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="DTY":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="TYY":
				residu.rName("TYR")
				residu.delete(["HN5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="O2":
						atom.atmName("HD1")
					if atom.atmName()=="N5":
						atom.atmName("HE2")
			if residu.rName()=="TYQ":
				residu.rName("TYR")
				residu.delete(["HN51","HN52","HOZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OZ":
						atom.atmName("HD1")
					if atom.atmName()=="N5":
						atom.atmName("HE2")
			# one H missing on C-ter
			if residu.rName()=="TYB":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="HC":
						atom.atmName("OXT")
			if residu.rName()=="STY":
				residu.rName("TYR")
				residu.delete(["O2","O3","O4","HO4"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="S":
						atom.atmName("HH")
			if residu.rName()=="PTR":
				residu.rName("TYR")
				residu.delete(["O1P","O2P","O3P","PHO2","PHO3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HH")
			if residu.rName()=="PAQ":
				residu.rName("TYR")
				residu.delete(["HN1","N2","NH2","C1","C2","C3","C4","C5",
				"H2","H3","H4","H5"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="O2":
						atom.atmName("HD1")
					if atom.atmName()=="N1":
						atom.atmName("HE2")
			if residu.rName()=="IYR":
				residu.rName("TYR")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="IE":
						atom.atmName("HE")
			if residu.rName()=="PHL":
				residu.rName("PHE")
				residu.delete(["H2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="H1":
						atom.atmName("O2")
			if residu.rName()=="PHI":
				residu.rName("PHE")
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="I":
						atom.atmName("HH")
			if residu.rName()=="MEN":
				residu.rName("ASN")
				residu.delete(["1HE2","2HE2","3HE2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CE2":
						atom.atmName("2HD2")
			if residu.rName()=="KCX":
				residu.rName("LYS")
				residu.delete(["HX2","OX2","OX1"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CX":
						atom.atmName("2H2")
			# one H missing on ND1
			if residu.rName()=="3AH":
				residu.rName("HIS")
				residu.delete(["N1","HN1","N2","C3","N4","N3A","1HN3","2HN3"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="C5":
						atom.atmName("HNE2")
			if residu.rName()=="DAH":
				residu.rName("TYR")
				residu.delete(["HOE"])
				residu.delete(["HE2"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="OE2":
						atom.atmName(" HE2")
					if atom.atmName()=="OZ":
						atom.atmName(" OH")
					if atom.atmName()=="HZ":
						atom.atmName(" HH")
			# one H missing on ND1
			if residu.rName()=="HIC":
				residu.rName("HIS")
				residu.delete(["1HZ","2HZ","3HZ"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CZ":
						atom.atmName("HE2")
			if residu.rName()=="HIP":
				residu.rName("HIS")
				residu.delete(["2HOP","3HOP","O1P","O2P","O3P"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HD1")
			#Nomenclature atom 31HN ?
			if residu.rName()=="HPQ":
				residu.rName("PHE")
				residu.delete(["2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("OXT")
					if atom.atmName()=="1HM":
						atom.atmName("HXT")
			if residu.rName()=="LTR":
				residu.rName("TRP")
				for atom in residu:
					atom.header("ATOM   ")
			if residu.rName()=="TPL":
				residu.rName("TRP")
				residu.delete(["1HC"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="2HC":
						atom.atmName("O")
			if residu.rName()=="MHS":
				residu.rName("HIS")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("HD1")
			if residu.rName()=="NEM":
				residu.rName("HIS")
				residu.delete(["1HM","2HM","3HM"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="CM":
						atom.atmName("HE2")
			if residu.rName()=="NEP":
				residu.rName("HIS")
				residu.delete(["O1P","O2P","O3P","1HOP","2HOP"])
				for atom in residu:
					atom.header("ATOM   ")
					if atom.atmName()=="P":
						atom.atmName("HE2")
# P. Tuffery, september 2008
#			if residu.rName()=="CSD":
#				residu.rName("CYS")
#				residu.delete(["O1P","O2P","O3P","1HOP","2HOP"])
#				for atom in residu:
#					atom.header("ATOM   ")
#					if atom.atmName()=="P":
#						atom.atmName("HE2")
			for atom in residu:
				atom.header("ATOM  ") # 6 chars for header

	def site(self, verbose = 0):
		"""
		Parse the info lines and check for a site description
		according to the PDB format
		@return: a dictionnary or None
		"""
		rs = None
		for Line in self.info:
			if Line[:6]=='SITE  ':
				if rs == None:
					rs = {}
				try:
					siteName = Line[11:14]
				except:
					return rs
				try:
					siteNRes = Line[15:17]
				except:
					return rs
				try:
					resName = Line[18:21]
				except:
					return rs
				try:
					resChn = Line[22]
				except:
					return rs
				try:
					resNum = Line[23:27]
				except:
					return rs
				try:
					resIcode = Line[27]
				except:
					return rs
				if rs.has_key(siteName):
					rs[siteName].append([resName,resChn,resNum,resIcode])
				else:
					rs[siteName] = [[resName,resChn,resNum,resIcode]]
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[29:32]
				except:
					pass
				try:
					resChn = Line[33]
				except:
					pass
				try:
					resNum = Line[34:38]
				except:
					pass
				try:
					resIcode = Line[38]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[40:43]
				except:
					pass
				try:
					resChn = Line[44]
				except:
					pass
				try:
					resNum = Line[45:49]
				except:
					pass
				try:
					resIcode = Line[49]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
				resName = resChn = resNum = resIcode = None
				try:
					resName = Line[51:54]
				except:
					pass
				try:
					resChn = Line[55]
				except:
					pass
				try:
					resNum = Line[56:60]
				except:
					pass
				try:
					resIcode = Line[60]
				except:
					pass
				if resName != None and resName != "   ":
					rs[siteName].append([resName,resChn,resNum,resIcode])
		return rs

	def CSAsite(self, id = None, verbose = 0):
		"""
		PDB.CSAsite
		@param id: a PDB id
		@return: a list of dictionnaries
		@note:
			- Attempt to retrieve site from Catalytic Site Atlas:
		at http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=2lzm
		returns either None or a dictionnary of the sites
			- CSA does not consider chains. Hence, one must parse it later.
			- Status: "Literature" or "PsiBLAST"
			- Referer: "PsiBLAST" match
			- Comment: Comment on psiBlast and EC.
			- Site:    Atoms involved
		"""
		#Attempt to retrieve site from Catalytic Site Atlas:
		#at http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=2lzm
		#returns either None or a dictionnary of the sites
		#Id must be a PDB id. CSA does not consider chains. Hence, one must parse it later.
		#Return a list of dictionnaries.
		#Status: "Literature" or "PsiBLAST"
		#Referer: "PsiBLAST" match
		#Comment: Comment on psiBlast and EC.
		#Site:    Atoms involved
		from urllib import urlretrieve
		chIds = None
		if id == None:
			cmd = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=%s" % self.id
		else:
			if len(id) > 4:
				chIds = id[4:]
				id    = id[:4]
			cmd = "http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/CSA/CSA_Site_Wrapper.pl?pdb=%s" % id
		file, log = urlretrieve(cmd)
		if verbose:
			print >> sys.stdout, "# urlretrieve at %s" % file
		CSA=simpleload(file,0)
		del urlretrieve
		grs = []
		rs = {}
		on = 0
		EC = ""
		Title = ""
		Compound = ""
		StatusOn = False
		TitleCount = 0
		CompoundCount = 0
		for i in CSA:
			# Site status
			if StatusOn:
				rs["EC"] = EC
				rs["Title"] = Title
				rs["Compound"] = Compound
				if string.count(i,"Literature reference"):
					rs["Status"] = "Literature"
				if string.count(i,"PsiBLAST"):
					# PsiBLAST alignment on <a href="CSA_Site_Wrapper.pl?pdb=1ds2">1ds2</a><br><font color=red>1ds2 has EC code 3.4.21.81 whereas 1ssx has EC code 0....<br>The difference in function suggests that the transfer of annotation from 1ds2 to 1ssx may be incorrect.</font>
					rs["Status"]  = "PsiBLAST"
					s = string.split(i,">")[1]
					r = string.split(s,"<")[0]
					rs["Referer"] = r
					try:
						aPos = string.index(i,">")
						aPos = string.index(i,">", aPos+1)
						aPos = string.index(i,">", aPos+1)
						aPos = string.index(i,">", aPos+1)
						rs["Comment"] = i[aPos+1:]
						rs["Comment"] = string.replace(rs["Comment"],"<br>"," ")
						rs["Comment"] = string.replace(rs["Comment"],"</font>\n","")
					except:
						rs["Comment"] = None
				StatusOn = False
			# New Site
			if string.count(i,"Found by:"):
				if rs != {}:
					if rs["Site"] != []:
						grs.append(rs)
				rs = {}
				rs["Site"] = []
				StatusOn = True
			if string.count(i,"http://www.ebi.ac.uk/intenz/query?cmd=SearchEC"):
				# <pre><a href="http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&amp;ec=1.1.1.1" target="_top">1.1.1.1</a>
				if EC == "":
					it = string.split(i,">")[2]
					lEC = string.replace(it,"</a","")
					if lEC != "\n":
						EC = lEC
			if TitleCount > 0:
				TitleCount -= 1
				if TitleCount == 0:
					Title = string.replace(i,"</div>","")
					Title = string.replace(Title,"\n","")
			if CompoundCount > 0:
				CompoundCount -= 1
				if CompoundCount == 0:
					Compound = string.replace(i,"</div>","")
					Compound = string.replace(Compound,"\n","")
			if string.count(i,"Title:"):
				TitleCount = 4
			if string.count(i,"Compound:"):
				CompoundCount = 4
			if string.count(i,"Residue"):
				on = 1
				continue
			if on:
				if string.count(i,"<td>"):
					l = string.replace(i,"<td>"," ")
					l = string.replace(l,"</td>\n"," ")
					resName = string.split(l)[0]
					resNum  = string.split(l[6:])[0]
					resChn  = l[5]
					what    = string.split(l[6:])[-1]
					if (chIds != None):
						if (resChn in chIds):
							rs["Site"].append( [resName,resNum,resChn,what])
					else:
						rs["Site"].append( [resName,resNum,resChn,what])
				if string.count(i,"</table>"):
					on = 0
		if rs != {}:
			if rs["Site"] != []:
				grs.append(rs)
		return grs

	def exposedAminoAcids(self, rH2O = "1.4", ASALimit = "0.25", what = "E", verbose = 0):
		"""
		PDB.exposedAminoAcids: determine if an amino acid is exposed to the solvant
		@param what:
		@param rH2O: the radius of the H2O molecule
		@param ASALimit: exposure threshold beyond which the residue is considered as exposed
		@return: a PDB instance with all the residues exposed to the solvant
		"""
		import ASA
		lines = ASA.ASA2(self, rH2O = rH2O, ASALimit = ASALimit, verbose = verbose)
		del ASA
		aaseq  = ""
		be     = ""
		toKeep = []
		for i in range(1,len(lines)-3):
			it = string.split(lines[i])
			# if it[-1] == "E":
			if it[-1] == what: # "E" or "B"
				toKeep.append([it[0],it[1],it[2]])
				be += it[-1]
			else:
				be += "-"
		res = []
		for i in self:
			if [i.rName(),i.chnLbl(),i.rNum()] in toKeep:
				res = res + i.flat()
			elif AA3STRICT.count(i.rName()) == 0:
				res = res + i.flat()
		return PDB(res)

	def BB(self):
		"""
		PDB.BB()
		@return: a PDB of the backbone atoms only
		"""
		theBB = self.select(awhat=BBATMS)
		return theBB

	def SC(self):
		"""
		PDB.SC()
		@return: a PDB of the side-chain atoms only
		"""
		theSC = self.select(awhat=SCATMS)
		return theSC

	def around(self, elt, aPos, dist = 3., verbose = 0):
		"""
		PDB.around search all the atoms or residus around aPos to a distance dist
		@param elt: "a" or "r", respectively for atoms or residus
		@param dist: distance in angstrom
		@param aPos: atom position
		@return: none, print the name and number of atoms around the atom position
		"""
		dist = float(dist)
		if elt == 'r':
			chain, entity = aPos.split("_")
			rs = []
			residu = self.findRes(chId=chain,rName=None,rNum=entity,icode=None,verbose=0)
			if not residu:
				print >> sys.stderr, "# No residu %s in chain %s" % (entity,chain)
				return
			for aAtm in residu:
				x, y, z = aAtm.xyz()
				for aRes in self:
					if aRes.rNum() == entity:
						continue
					aCa = aRes.findAtm("CA")
					if aCa == None:
						continue
					CAx, CAy, CAz = aCa.xyz()
					d =  distance(x, y, z, CAx, CAy, CAz)
					if d > 15.:
						continue
					if d < dist:
						if aRes not in rs:
							rs.append(aRes)
						continue
					for aAtm2 in aRes:
						CAx, CAy, CAz = aAtm2.xyz()
						d =  distance(x, y, z, CAx, CAy, CAz)
						if d < dist:
							if aRes not in rs:
								rs.append(aRes)
							break
					continue
			rs.sort()
			if verbose:
				print >> sys.stdout, "# Residu around \"%s %s %s\" :" % (residu.chnLbl(), residu.rNum(), residu.rName())
			for aRes in rs:
				print >> sys.stdout, aRes.chnLbl(), aRes.rNum(), aRes.rName()
		elif elt == 'a':
			if len(str(aPos).split("_")) > 1:
				aPos = int(aPos.split("_")[1])
			else:
				aPos = int(aPos)
			# have to be improved to have a "cutoff, like for residus
			atms = []
			atom = None
			for atm in self.atms:
				if int(atm.atmNum()) == aPos:
					atom = atm
					x,y,z = atom.xyz()
					break
			if not atom:
				print >> sys.stderr, "# No atom %d" % aPos
				return
			for atm in self.atms:
				if atm == atom:
					continue
				xAtm,yAtm,zAtm = atm.xyz()
				d = distance(x,y,z, xAtm,yAtm,zAtm)
				if d < dist:
					if atm not in atms:
						atms.append(atm)
			atms.sort()
			if verbose:
				print >> sys.stdout, "# Atoms around \"%s %s %s - %s %s\" :" % (atom.chnLbl(), atom.resNum(), atom.resName(), atom.atmNum(), atom.atmName())
			for atm in atms:
				print >> sys.stdout, "%s %s %s - %s %s" % (atm.chnLbl(), atm.resNum(), atm.resName(), atm.atmNum(), atm.atmName())

	def addHydrogens(self, algorithm = "HAAD", HSkip = True, norm = None, verbose = 0):
		"""
		PDB.addHydrogens add Hydrogens to the current structure using HAAD (Li,
		    et al.(2009) "HAAD: A Quick Algorithm for Accurate Prediction of
		    Hydrogen Atoms in Protein Structures" PLoS One, 4: e6701.) or Reduce
		    (Word, et al.(1999) "Asparagine and glutamine: using hydrogen atom
		    contacts in the choice of sidechain amide orientation" J. Mol. Biol.
		    285, 1735-1747) methods
		@author: F.Briand
		@param method: Method used (HAAD, Reduce)
		@param HSkip: Does the hydrogens already in place are skipped or not ? (True/False)
		@param norm: Which norm have to be applied ? None (keep norm of the algorithm), lIUPAC, lPDB
		@param verbose: if > 0, step avancement is printed; if > 1, CONECT fields changes are printed
		@return: None
		@note: Hydrogens are added "in place"
		"""
		inLineID = 0
		stderrPipe = " 2> /dev/null"
		tempDir = tempfile.mkdtemp()
		os.chmod(tempDir,stat.S_IRWXU)
		for chain in self.chnList():
			if verbose:
				print >> sys.stdout, "# Adding hydrogen to chain %s ..." % chain
			beforeOutfile = tempDir + "/pdb_" + chain + "_before_add_H.pdb"
			afterOutfile = tempDir + "/pdb_" + chain + "_after_add_H.pdb"
			if method == "HAAD":
				if verbose:
					print >> sys.stdout, "# Running HAAD ..."
				self.renameHydrogens(norm = "lPDB")
				self[chain].out(beforeOutfile, header = 0, HSkip = int(HSkip))
				rt = subprocess.call(HAADBIN + " " + beforeOutfile, shell=True)
				os.rename(beforeOutfile+".h",afterOutfile)
				chainHydrogen = PDB(afterOutfile, hetSkip = int(HSkip))
			elif method == "Reduce":
				if verbose:
					print >> sys.stdout, "# Running Reduce ... (log file : %sReduce_Add_H_Info.log )" % os.getcwd()
				self.renameHydrogens(norm = "lIUPAC")
				self[chain].out(beforeOutfile, header = 0, HSkip = int(HSkip))
				if verbose:
					stderrPipe = " 2> "+ os.getcwd() + "Reduce_Add_H_Info.log"
				rt = subprocess.call(REDUCEBIN + " " + beforeOutfile + " 1> " + afterOutfile + stderrPipe, shell=True)
				chainHydrogen = PDB(afterOutfile, hetSkip = True)
			for line in chainHydrogen.atms:
				try:
					if self.atms[inLineID].txt[:27] != line.txt[:27]:
						line.chnLbl(chain)
						if re.match("[A-Z]",line.txt[12]):
							atmType = line.txt[12]
						else:
							atmType = line.txt[13]
						if method == "HAAD":
							endline = "  1.00  0.00           %s  \n" % atmType
						elif method == "Reduce":
							endline = line.txt[54:80]+"\n"
						line.txt = "%s%s" % (line.txt[:54], endline)
						self.atms[inLineID:inLineID] = [line]
					inLineID += 1
				except IndexError:
					break
			os.remove(beforeOutfile)
			os.remove(afterOutfile)
		os.rmdir(tempDir)
		self.atmsForceRenumber(verbose=verbose)
		if norm:
			if verbose:
				print >> sys.stdout, "# Renaming hydrogens in norm %s" % norm
			self.renameHydrogens(norm=norm)

	def atmsForceRenumber(self, verbose = 0):
		"""
		PDB.atmsForceRenumber renumber all the atoms of a PDB structure, starting to the first atom with index equal to its atom number. Also increase the index by 1 on TER lines.
		@author: F.Briand
		@param verbose: if > 0, step avancement is printed; if > 1, CONECT fields changes are printed
		@return: None
		@note: renumbering "in place"
		"""
		if verbose:
			print >> sys.stdout, "# Renumbering new atoms ..."
		index = int(self.atms[0].atmNum())
		chain = self[0][0].chnLbl()
		atomsDict = {}
		for line in self.atms:
			if line.chnLbl() != chain:
				index += 1
				chain = line.chnLbl()
			atomsDict[line.atmNum()] = index
			if verbose > 1:
				print >> sys.stdout, "# " + str(line.atmNum()) + " => " + str(index)
			line.atmNum(index)
			index += 1
		if verbose:
			print >> sys.stdout, "# Renumbering CONECT fields ..."
		for iConect in range(0,len(self.conect)):
			if verbose > 1:
				print >> sys.stdout, "# " + self.conect[iConect],
			for iAtom in range(0,14):
				try:
					atom = int(self.conect[iConect][6+(iAtom*5):11+(iAtom*5)])
					if atom:
						self.conect[iConect] = "%s%5d%s" % (self.conect[iConect][:6+(iAtom*5)],atomsDict[str(atom)],self.conect[iConect][11+(iAtom*5):])
				except ValueError:
					pass
				except KeyError: #!# "Keep" CONECT fields not matching an existing atom. Have to be improved to "delete" these fields.
					pass
			if verbose > 1:
				print >> sys.stdout, "# " + self.conect[iConect]
		self.resTab(0)

	def delHydrogens(self, verbose = 0):
		"""
		PDB.delHydrogens delete Hydrogens of a PDB structure then renumber atoms
		@author: F.Briand
		@param verbose: if > 0, step avancement is printed; if > 1, CONECT fields changes are printed
		@return: None
		@note: Hydrogens are deleted "in place"
		"""
		if verbose:
			print >> sys.stdout, "# Deleting Hydrogens ..."
		self.atms = PDB(self, keepH = 0).atms
		self.atmsForceRenumber(verbose=verbose)

	def trace(self, hetSkip = 0, verbose = 0):
		"""
		PDB.trace:
		@param hetSkip: does hetero atom are skipped.0 skip nothing, 1 skip no-amino acids atoms, 2 skip no-amino acids & "special" amino acids.
		@return: a PDB of the CA
		"""
		atmsTraceList = PDB(self)
		atmId = 0
		while atmId < len(atmsTraceList.atms):
			if (atmsTraceList.atms[atmId].atmName() != "CA") or (hetSkip == 1 and AA3.count(atmsTraceList.atms[atmId].resName()) == 0) or (hetSkip > 1 and AA3STRICT.count(atmsTraceList.atms[atmId].resName()) == 0):
				atmsTraceList.atms[atmId:atmId+1] = []
				atmId -= 1
			atmId += 1
		atmsTraceList.resTab(0)
		return atmsTraceList

	def getAtoms(self, pattern='', byID = False, verbose = 0):
		"""
		PDB.getAtoms extract a list of atoms from a PDB structure
		@author: F.Briand
		@param pattern: \"ffrom:tto\", tto not included
		@param byID: does the function select atoms by python list index (False), or by atom number (True).
		@return: a PDB structure
		"""
		try:
			ffrom, tto = (int(x) for x in pattern.split(":"))
		except ValueError:
			raise ValueError, "PDB.getAtoms() : pattern error. USAGE : \"ffrom(integer):tto(integer)\""
		subPDB = PDB(self)
		if byID:
			iAtom = 0
			while iAtom < len(subPDB.atms):
				if int(subPDB.atms[iAtom].atmNum()) < ffrom or int(subPDB.atms[iAtom].atmNum()) >= tto:
					subPDB.atms[iAtom:iAtom+1] = []
					iAtom -= 1
				iAtom += 1
		else:
			subPDB.atms = self.atms[ffrom:tto]
		subPDB.resTab(0)
		return subPDB

	def getMultiAtoms(self, pattern='', verbose = 0):
		"""
		PDB.getMultiAtoms extract lists of atoms from a PDB structure
		@author: F.Briand
		@param pattern: \"from1:to1 from2:to2 ...\", to not included
		@return a PDB structure
		"""
		for patterns in pattern.split():
			if verbose:
				print >> sys.stdout, "# Extracting %s atoms" % patterns
			try:
				newPDB += self.getAtoms(patterns, byID = 1, verbose = verbose)
			except NameError:
				newPDB = self.getAtoms(patterns, byID = 1, verbose = verbose)
		return newPDB

	def addOXT(self, chainId = "", verbose = 0):
		"""
		PDB.addOXT add terminal oxygens
		@author: F.Briand
		@param chainId: set chains which have their OXT to be added. "" = all the chains
		"""
		chnList = self.chnList()
		if not chainId:
			chainId = chnList
		else:
			for chain in chainId:
				if chain not in chnList:
					chainId = chainId.replace(chain, "")
		if verbose:
			print >> sys.stdout, "# Add of OXT in chains %s ..." % chainId
		for chain in chainId:
			OXT = False
			CA, CP, O = None, None, None
			for atom in self.atms:
				if atom.chnLbl() == chain and atom.resType() == 'AMINO-ACID':
					aName = atom.atmName()
					if aName == "CA":
						CA = atom
					elif aName == "C":
						CP = atom
					elif aName == "O":
						O = atom
					elif aName == "OXT":
						OXT = True
					lastAtom = atom
			if not OXT:
				OXT = xyzOXT(ca=CA, cp=CP, o=O, verbose=verbose)
				insertKey = self.atms.index(lastAtom)
				lineOXT = atmLine(lastAtom)
				lineOXT.atmName(" OXT")
				lineOXT.atmNum(0)
				lineOXT.atmBVal(0)
				lineOXT.atmType('O')
				lineOXT.occ(1)
				lineOXT.setcrds(OXT[0],OXT[1],OXT[2])
				self.atms.insert(insertKey+1, lineOXT)
				if verbose:
					print >> sys.stdout, "#", self.atms[insertKey+1],
		self.atmsForceRenumber()
		self.resTab(0)

	def rmk350apply(self, biomolecule = None, verbose = 0):
		"""
		PDB.rmk350apply apply the translation/rotation matrix which is in REMARK 350 fields
		@author: F.Briand
		@param: specify which biomolecule (model) is targeted by the transformation
		@return: a PDB with the matrix applied
		"""
		toDel = []
		newAllAtms = PDB(self)
		newAllAtms.atms = []
		for rmkId in range(0,len(self.remark350)):
			if str(biomolecule) in self.remark350[rmkId].biomolecule or not biomolecule:
				for chain in self.remark350[rmkId].chains:
					if verbose:
						print >> sys.stdout, "# Operating on biomolecule "+str(biomolecule)+", chain "+chain
					for key in self.remark350[rmkId]['biomt'].keys():
						for chnLetter in string.ascii_uppercase+string.join(map(str,range(0,10)),""):
							if chnLetter not in newAllAtms.chnList():
								break
						rotation = [list(self.remark350[rmkId]['biomt'][key][0][:-1]),list(self.remark350[rmkId]['biomt'][key][1][:-1]),list(self.remark350[rmkId]['biomt'][key][2][:-1])]
						translation = [self.remark350[rmkId]['biomt'][key][0][-1], self.remark350[rmkId]['biomt'][key][1][-1], self.remark350[rmkId]['biomt'][key][2][-1]]
						if verbose:
							print >> sys.stdout, "# Rotation Matrix :", rotation
							print >> sys.stdout, "# Translation Matrix :", translation
						newAtms = PDB(self[chain].atms)
						for atm in newAtms.atms:
							x1,y1,z1 = atm.xyz()
							x2 = float(rotation[0][0])*x1 + float(rotation[0][1])*y1 + float(rotation[0][2])*z1
							y2 = float(rotation[1][0])*x1 + float(rotation[1][1])*y1 + float(rotation[1][2])*z1
							z2 = float(rotation[2][0])*x1 + float(rotation[2][1])*y1 + float(rotation[2][2])*z1
							x3, y3, z3 = x2+float(translation[0]), y2+float(translation[1]), z2+float(translation[2])
							atm.setcrds(x3,y3,z3)
							atm.chnLbl(chnLetter)
						newAtms.resTab(0)
						newAllAtms += newAtms
				toDel.append(rmkId)
		toDel.sort(reverse=True)
		for id in toDel:
			del(newAllAtms.remark350[id])
		newAllAtms.atmsForceRenumber()
		return newAllAtms

	def checkBBOrder(self, verbose = 0):
		"""
		PDB.checkBBOrder check that backbone atoms order is : N, CA, C, O, (OXT) and restore correct order if it wasn't.
		@author: F.Briand
		"""
		for res in self.rt:
			BB = 0
			for atom in res:
				if atom.atmName() in BBATMS:
					if atom.atmName() == BBATMS[BB]:
						BB += 1
					else:
						if len(res) != 1:
							if verbose:
								print >> sys.stderr, "# Error in residue " + atom.resName() + atom.resNum()
							res.setBBOrder(verbose = verbose)
							break
		self.atmTab()
		self.atmsForceRenumber()


def xyzOXT(ca, cp, o, verbose = 0):
	"""
	xyzOXT return terminal oxygen position
	@author: F.Briand
	@param ca: atmLine of the Carbon Alpha
	@param cp: atmLine of the backbone oxygene (prime)
	@param o: atmLine of the backbone oxygene
	@return: tuple (x,y,z) with OXT coordinates
	"""
	ca, cp, o = ca.xyz(), cp.xyz(), o.xyz()
	# Vectors coordinates
	CpCa = (ca[0]-cp[0], ca[1]-cp[1], ca[2]-cp[2]) # (ax, ay, az)
	CpO = (o[0]-cp[0], o[1]-cp[1], o[2]-cp[2]) # (bx, by, bz)
	# First vectorial product : define Z axe
	#    cx = ay * bz - az * by
	#    cy = az * bx - ax * bz
	#    cz = ax * by - ay * bx
	CpZ = (CpCa[1]*CpO[2]-CpCa[2]*CpO[1], CpCa[2]*CpO[0]-CpCa[0]*CpO[2], CpCa[0]*CpO[1]-CpCa[1]*CpO[0])
	# Compute vectors norm
	CpZ_norm = math.sqrt((CpZ[0]*CpZ[0]) + (CpZ[1]*CpZ[1]) + (CpZ[2]*CpZ[2]))
	# Normalize k vector director (Z)
	CpZ = (CpZ[0]/CpZ_norm, CpZ[1]/CpZ_norm, CpZ[2]/CpZ_norm)
	# Second vectorial product : define Y axe (X is along CpCa)
	CpY = (CpZ[1]*CpCa[2]-CpZ[2]*CpCa[1], CpZ[2]*CpCa[0]-CpZ[0]*CpCa[2], CpZ[0]*CpCa[1]-CpZ[1]*CpCa[0])
	# Compute vectors norm
	CpY_norm = math.sqrt((CpY[0]*CpY[0]) + (CpY[1]*CpY[1]) + (CpY[2]*CpY[2]))
	# Normalize j vector director (Y)
	CpY = (CpY[0]/CpY_norm, CpY[1]/CpY_norm, CpY[2]/CpY_norm)
	# Compute vectors norm
	CpCa_norm = math.sqrt((CpCa[0]*CpCa[0]) + (CpCa[1]*CpCa[1]) + (CpCa[2]*CpCa[2]))
	# Normalize i vector director (X)
	CpX = (CpCa[0]/CpCa_norm, CpCa[1]/CpCa_norm, CpCa[2]/CpCa_norm)
	# Compute O coordinates in local repair
	Otr = (o[0]-cp[0], o[1]-cp[1], o[2]-cp[2])
	Oloc = (Otr[0]*CpX[0]+Otr[1]*CpX[1]+Otr[2]*CpX[2], Otr[0]*CpY[0]+Otr[1]*CpY[1]+Otr[2]*CpY[2], Otr[0]*CpZ[0]+Otr[1]*CpZ[1]+Otr[2]*CpZ[2])
	# Transform O local coordinates ( axial symetry
	# along x ) to have local OXT coordinates
	OxtLoc = (Oloc[0], -Oloc[1],  Oloc[2])
	# Return O world coodinates
	oxt = (OxtLoc[0]*CpX[0]+OxtLoc[1]*CpY[0]+OxtLoc[2]*CpZ[0]+cp[0], OxtLoc[0]*CpX[1]+OxtLoc[1]*CpY[1]+OxtLoc[2]*CpZ[1]+cp[1], OxtLoc[0]*CpX[2]+OxtLoc[1]*CpY[2]+OxtLoc[2]*CpZ[2]+cp[2])
	return oxt

def PDBBiologicalUnit(PDBid = None, verbose = 0):
	"""
	PDBBiologicaUnit:
	@param PDBid: none, or a PDBid
	@return: the PDB entry biological unit at PQS server (EBI: http://pqs.ebi.ac.uk/pqs-doc/macmol/2lzm.mmol)
	"""
	from urllib import urlretrieve
	if PDBid == None:
		return []
	cmd = "http://pqs.ebi.ac.uk/pqs-doc/macmol/%s.mmol" % (string.lower(PDBid))
	file, log = urlretrieve(cmd)
	if verbose:
		print >> sys.stdout, "# urlretrieve at %s" % file
	list=simpleload(file,0)
	del urlretrieve
	x  =PDB(list)
	return x

def PDBEntries(what = None, isauthor = "no", verbose = 0):
	"""
	PDBEntries
	@param what: a word searched on PDB.org
	@param isauthor: none, or an author
	@return: a list of entries matching a word on http://www.pdb.org/pdb/navbarsearch.do?newSearch=yes&isAuthorSearch=no&radioset=All &inputQuickSearch=calpain&image.x=0&image.y=0&image=Search
	"""
	from urllib import urlretrieve
	if what == None:
		return []
	cmd = "http://www.pdb.org/pdb/navbarsearch.do?newSearch=yes&isAuthorSearch=%s&radioset=All&inputQuickSearch=%s&image.x=0&image.y=0&image=Search" % (isauthor, what)
	file, log = urlretrieve(cmd)
	if verbose:
		print >> sys.stdout, "urlretrieve at %s" % file
	list=simpleload(file,0)
	del urlretrieve
##	f = open("calpain","w")
##	for i in list:
##		f.write("%s" % i)
##	f.close()
	rs = []
	for i in list:
		if string.count(i,"<input type=\"checkbox\" name="):
			j = string.split(string.replace(i,"<input type=\"checkbox\" name=",""))[0]
			k = string.replace(j,"\"","")
			rs.append(k)
	return rs

def CSASite2Escan(PDBid = None, patterns = "strict", purge = 0, verbose = 0):
	"""
	CSASite2Escan extract Catalytic Site atoms for a PDB.
	It gets the information directly at CSA.
	@param PDBid: a PDB file
	@param purge: does it keep only compatible sites? (No: 0 by default)
	@return: a list of atoms involved.
	@note: Will not take into account psiblast sites.
	@param patterns: one of "strict", "medium", "light"
	@type patterns: "strict": exact atom name match (by default); "medium": atom name match using atom class compatible pattern!; "light" : atom name match using light maks (atomic type)
	"""
	# extracted from Catalytic Site Atlas
	catalyticAtoms = {
		"ASP" : [["CG",  "OD1", "OD2"],["CG",  "OD1", "OD2"]],
		"GLU" : [["CD",  "OE1", "OE2"],["CG",  "CD",  "OE1", "OE2"]],
		"ASN" : [["CG",  "OD1", "ND2"],["CG",  "OD1", "ND2"]],
		"GLN" : [["CD",  "OE1", "NE2"],["CG",  "CD",  "OE1", "NE2"]],
		"HIS" : [["ND1", "NE2"],["CG",  "ND1", "NE2"]],
		"ARG" : [["NH1", "NH2"],["NE",  "NH1", "NH2"]],
		"LYS" : [["CE",  "NZ"],["CD",  "CE",  "NZ"]],
		"SER" : [["CB",  "OG"],["CB",  "OG"]],
		"THR" : [["CB",  "OG", "OG1"], ["CB",  "OG", "OG1"]],
		"CYS" : [["CB",  "SG"],["CB",  "SG"]],
		"ALA" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"GLY" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"LEU" : [["CG",  "CD1", "CD2"],["CG",  "CD1", "CD2"]],
		"PHE" : [["CE1", "CE2", "CZ"],["CE1", "CE2", "CZ"]],
		"TRP" : [["NE1"],["NE1", "CZ2", "CH2"]],
		"TYR" : [["CZ",  "OH"],["CE1", "CZ",  "OH"]],
		"MET" : [["SD"],["SD",  "CE"]],
		}
	# Remark: Some hetero atoms such as ZN/NAD (1qlh) could interfer
	# Patterns to match
	# The first is for similar class, the second pattern is light
	matchPatterns = {
		"OD" : ["OD|OE", "O.*"],
		"OE" : ["OD|OE", "O.*"],
		"OG" : ["OG|OH", "SG|OG|OH"],
		"OH" : ["OG|OH", "SG|OG|OH"],
		"NE" : ["ND.*|NE.*", "N.*"],
		"ND2": ["ND2|NE2", "ND|NE"],
		"ND1": ["ND1|NE1", "ND|NE"],
		"NZ" : ["NZ", "N.*"],
		"NH" : ["NH.*", "N.*"],
		"SD" : ["SD", "S.*"],
		"SG" : ["SG", "SG|OG|OH"],
		"CG" : ["CG|CD|CE", "C.*"],
		"CD" : ["CG|CD|CE", "C.*"],
		"CE" : ["CG|CD|CE", "C.*"],
		"CZ" : ["CZ", "C.*"],
		"CH" : ["CH", "C.*"],
		"CB" : ["CB|CG|CD|CE", "C.*"],
		"N"  : ["N","N"],
		"CA" : ["CA","CA"],
		"C"  : ["C","C"],
		}
	rs = None
	if PDBid == None:
		return None
	x = PDB(PDBid)
	if (x == None) or (len(x) == 0):
		return rs
	sites = x.CSAsite()
	if purge:
		sites = purgeSites(sites)
	if verbose > 1:
		print >> sys.stdout, "#", sites
	if len(sites) == 0:
		return None
	chnList = x.chnList()
	if patterns == "light":
		rank = 1
	elif patterns == "medium":
		rank = 0
	grs = []
	if verbose:
		print >> sys.stdout, "# %s : Found %d site(s)" % (PDBid, len(sites))
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				print >> sys.stdout, "# %s: PsiBlast referer : %s" % (PDBid, site["Referer"])
			continue
		rs = []
		cmpLine = "COMPND    %-70s\n" % site["Compound"]
		ECLine  = "REMARK    EC: %-60s\n" % site["EC"]
		rs.append(cmpLine)
		rs.append(ECLine)
		# Check for Hetero Groups
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			if AA3.count(aItem[0]) == 0:
				HTLine = "REMARK    HETGRP: %s\n" % aItem[0]
				rs.append(HTLine)
		# if len(site["Site"]) == 1:
		# Check if monoresidue output !
		count = 0
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			count += 1
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					count -= 1
					break
		if count < 2:
			rs.append("%-70s\n" % "REMARK    MONORESIDUE SITE")
		# Here, we format the site
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					if verbose:
						print >> sys.stderr, "%s: unreferenced catalytic atom %s %s" % (PDBid, aAtm.resName(), aAtm.atmName())
					continue
				if aAtm.atmName() not in catalyticAtoms[aItem[0]][0]:
					continue
				rs.append(aAtm.flat())
				if patterns != "strict":
					try:
						rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
					except:
						try:
							rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
						except:
							pass
				else:
					rs.append("REMARK     RESMATCH. %s\n" % aItem[0])
		if len(rs) > 2:
			grs.append(atmList(rs))
	return grs

def purgeSites(sites, hetCheck = 0, verbose = 0):
	"""
	purgeSites only keep compatible sites
	@param sites: a list of sites
	@param hetCheck: if 1, also check heteros
	@return: the list of compatible sites.
	"""
	rs = []
	for aSite in range(0,len(sites)):
		OK = 1
		for aSite2 in range(aSite+1, len(sites)):
			id = identicalp(sites[aSite],sites[aSite2], hetCheck = hetCheck, verbose =verbose)
			if id == 1:
				 OK = 0
				 break
		if OK:
			rs.append(sites[aSite])
	return rs


def identicalp(site1, site2, hetCheck = 0, verbose = 0):
	"""
	Check if site1 is compatible with site2	on the basis of residue names, residue number and which part
	(sidechain, backbone)
	@param hetCheck: if 1, also check heteros
	@return: the compatibility of the both sites (0: no Matches, 1: Matches)
	"""
	if verbose:
		print >> sys.stdout, "#",site1
		print >> sys.stdout, "#",site2
	if site1["EC"] != site2["EC"]:
		return 0
	m = []
	h1 = []
	h2 = []
	for aRes2 in site1["Site"]:
		if AA3.count(aRes2[0]) == 0:
			h1.append(aRes2[0])
	for aRes2 in site2["Site"]:
		m.append(0)
		if AA3.count(aRes2[0]) == 0:
			h2.append(aRes2[0])
	if len(site1["Site"]) - len(h1) != len(site2["Site"]) - len(h2):
		return 0
	for aRes in  site1["Site"]:
		OK = 0
		count = -1
		for aRes2 in site2["Site"]:
			count += 1
			if m[count]: # this residue already assigned
				continue
			if (aRes[0] == aRes2[0]) and (aRes[1] == aRes2[1]) and (aRes[3] == aRes2[3]):
				OK = 1
				# we flag aRes2 as matched
				if verbose:
					print >> sys.stdout, "#", aRes, "matches", aRes2
				break
		if (OK == 0):
			if AA3.count(aRes[0]) == 0:
				if hetCheck == 1:
					return 0
			else:
				if verbose:
					print >> sys.stdout, "# No match for",aRes
				return 0
	if verbose:
		print >> sys.stdout, "# samep" ### /!\
	# amino acids are OK
	# what about het groups
	if len(h1) and len(h2):
		if len(h1) == len(h2):
			OK = 1
			for i in h1:
				if i not in h2:
					OK = 0
					break
			return OK
		else:
			return 0
	elif len(h1) or len(h2):
		if hetCheck == 0:
			if len(h1) > len(h2):
				return 2
			else:
				return 1
	return 1

def samep(site1, site2, hetCheck = 0, verbose = 0):
	"""
	Check if site1 is compatible with site2	on the basis of residue names, residue number and which part
	(sidechain, backbone)
	@param hetCheck: if 1, also check heteros
	@return: the compatibility of the both sites (0: no Matches, 1: Matches)
	"""
	if verbose:
		print >> sys.stdout, "#", site1
		print >> sys.stdout, "#", site2
##	if site1["EC"] != site2["EC"]:
##		return 0
	h1 = []
	h2 = []
	m = []
	for aRes2 in site1["Site"]:
		if AA3.count(aRes2[0]) == 0:
			h1.append(aRes2[0])
	for aRes2 in site2["Site"]:
		m.append(0)
		if AA3.count(aRes2[0]) == 0:
			h2.append(aRes2[0])
	if len(site1["Site"]) - len(h1) != len(site2["Site"]) - len(h2):
		return 0
	for aRes in  site1["Site"]:
		OK = 0
		count = -1
		for aRes2 in site2["Site"]:
			count += 1
			if m[count]: # this residue already assigned
				continue
			if (aRes[0] == aRes2[0]) and (aRes[3] == aRes2[3]):
				OK = 1
				# we flag aRes2 as matched
				if verbose:
					print >> sys.stdout, "#", aRes, "matches", aRes2
				break
		if (OK == 0):
			if AA3.count(aRes[0]) == 0:
				if hetCheck == 1:
					if verbose:
						print >> sys.stdout, "# Het inconsistency"
					return 0
			else:
				if verbose:
					print >> sys.stdout, "# No match for",aRes
				return 0
	if verbose:
		print >> sys.stdout, "# samep"
	return 1

def EscanCASSites(PDBid = None, patterns = "strict", purge = 0, verbose = 0):
	"""
	extract Catalytic Site atoms for a PDB.
	recurse to validated catalytic sites (follow referer if required).
	Then gets the information directly at CSA.
	@param PDBid: a PDB file
	@param purge: does it keep only compatible sites? (No: 0 by default)
	@return: a list of atoms involved.
	@param patterns: one of "strict", "medium", "light"
	@type patterns: "strict": exact atom name match; "medium": atom name match using atom class compatible pattern; "light" : atom name match using light maks (atomic type)
	"""
	# extracted from Catalytic Site Atlas
	catalyticAtoms = {
		"ASP" : [["CG",  "OD1", "OD2"],["CG",  "OD1", "OD2"]],
		"GLU" : [["CD",  "OE1", "OE2"],["CG",  "CD",  "OE1", "OE2"]],
		"ASN" : [["CG",  "OD1", "ND2"],["CG",  "OD1", "ND2"]],
		"GLN" : [["CD",  "OE1", "NE2"],["CG",  "CD",  "OE1", "NE2"]],
		"HIS" : [["ND1", "NE2"],["CG",  "ND1", "NE2"]],
		"ARG" : [["NH1", "NH2"],["NE",  "NH1", "NH2"]],
		"LYS" : [["CE",  "NZ"],["CD",  "CE",  "NZ"]],
		"SER" : [["CB",  "OG"],["CB",  "OG"]],
		"THR" : [["CB",  "OG", "OG1"], ["CB",  "OG", "OG1"]],
		"CYS" : [["CB",  "SG"],["CB",  "SG"]],
		"ALA" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"GLY" : [["N",   "CA",  "C"],["N",   "CA",  "C"]],
		"LEU" : [["CG",  "CD1", "CD2"],["CG",  "CD1", "CD2"]],
		"PHE" : [["CE1", "CE2", "CZ"],["CE1", "CE2", "CZ"]],
		"TRP" : [["NE1"],["NE1", "CZ2", "CH2"]],
		"TYR" : [["CZ",  "OH"],["CE1", "CZ",  "OH"]],
		"MET" : [["SD"],["SD",  "CE"]],
		}
	# Remark: Some hetero atomes such as ZN/NAD (1qlh) peuvent intervenir
	# Patterns to match
	# The first is for similar class, the second pattern is light
	matchPatterns = {
		"OD" : ["OD|OE", "O.*"],
		"OE" : ["OD|OE", "O.*"],
		"OG" : ["OG|OH", "SG|OG|OH"],
		"OH" : ["OG|OH", "SG|OG|OH"],
		"NE" : ["ND.*|NE.*", "N.*"],
		"ND2": ["ND2|NE2", "ND|NE"],
		"ND1": ["ND1|NE1", "ND|NE"],
		"NZ" : ["NZ", "N.*"],
		"NH" : ["NH.*", "N.*"],
		"SD" : ["SD", "S.*"],
		"SG" : ["SG", "SG|OG|OH"],
		"CG" : ["CG|CD|CE", "C.*"],
		"CD" : ["CG|CD|CE", "C.*"],
		"CE" : ["CG|CD|CE", "C.*"],
		"CZ" : ["CZ", "C.*"],
		"CH" : ["CH", "C.*"],
		"CB" : ["CB|CG|CD|CE", "C.*"],
		"N"  : ["N","N"],
		"CA" : ["CA","CA"],
		"C"  : ["C","C"],
		}
	rs = None
	if PDBid == None:
		return None
	x = PDB(PDBid)
	if (x == None) or (len(x) == 0):
		return rs
	sites = x.CSAsite()
	if purge:
		sites = purgeSites(sites)
	if verbose > 1:
		print >> sys.stdout, "#", sites
	if len(sites) == 0:
		return None
	chnList = x.chnList()
	if patterns == "light":
		rank = 1
	elif patterns == "medium":
		rank = 0
	grs = []
	if verbose:
		print >> sys.stdout, "# %s : Found %d site(s)" % (PDBid, len(sites))
	# 1. We check that these are validated sites
	tsites = []
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				print >> sys.stdout, "# %s: PsiBlast referer : %s" % (PDBid, site["Referer"])
			y = PDB(site["Referer"])
			ysites = y.CSAsite()
			if purge:
				ysites = purgeSites(ysites)
			for ysite in ysites:
				if ysite["Status"] != "Literature":
					continue
				# Now we check on residue names and what part (sidechain, etc)
				if samep(site,ysite, verbose = 0):
					ysite["Id"] = site["Referer"]
					tsites.append(ysite)
					break
			continue
		else:
			site["Id"] = PDBid
			tsites.append(site)
	sites = tsites
	if verbose:
		print >> sys.stdout, "# Before purge:",
		print >> sys.stdout, sites
	if purge:
		sites = purgeSites(sites)
	if verbose:
		print >> sys.stdout, "# After purge:",
		print >> sys.stdout, sites
	orix = x
	oriChnList = chnList
	for site in sites:
		if site["Status"] != "Literature":
			if verbose:
				print >> sys.stdout, "# %s: PsiBlast referer : %s" % (PDBid, site["Referer"])
			continue
		if site["Id"] == PDBid:
			x = orix
			chnList = oriChnList
		else:
			if verbose:
				print >> sys.stdout, "# Loading %s" % site["Id"]
			x = PDB(site["Id"])
			chnList = x.chnList()
		rs = []
		if site["Id"] == PDBid:
			# PDB HEADER ID is 63-66
			# headLine = "HEADER   %s CSA literature site %s\n" % (x.id,x.id)
			headLine = "HEADER    %s CSA literature site                            %s\n" % (x.id,x.id)
		else:
			# headLine = "HEADER   %s CSA psiblast site  %s\n" % (x.id,x.id)
			headLine = "HEADER    %s CSA psiblast %s site                         %s\n" % (PDBid, x.id,x.id)
		cmpLine = "COMPND    %-70s\n" % site["Compound"]
		ECLine  = "REMARK    EC: %-60s\n" % site["EC"]
		rs.append(headLine)
		rs.append(cmpLine)
		rs.append(ECLine)
		# Check for Hetero Groups
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			if AA3.count(aItem[0]) == 0:
				HTLine = "REMARK    HETGRP: %s\n" % aItem[0]
				rs.append(HTLine)
		# if len(site["Site"]) == 1:
		# Check if monoresidue output !
		count = 0
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			count += 1
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					count -= 1
					break
		if count < 2:
			rs.append("%-70s\n" % "REMARK    MONORESIDUE SITE")
		# Here, we format the site
		for aItem in site["Site"]:
			if aItem[2] not in chnList:
				continue
			# rs.append(aItem)
			try:
				aRes  = x.findRes(aItem[2],aItem[0],aItem[1], icode= "", what = None, verbose = 0)
			except:
				continue
			for aAtm in aRes:
				if (aItem[-1] == "Sidechain") and (aAtm.atmName() in ["N","CA","C","O","OXT","CT"]):
					continue
				try:
					isAtm = aAtm.atmName() not in catalyticAtoms[aItem[0]][0]
				except:
					if verbose:
						print >> sys.stderr, "# %s: unreferenced catalytic atom %s %s" % (PDBid, aAtm.resName(), aAtm.atmName())
					continue
				if aAtm.atmName() not in catalyticAtoms[aItem[0]][0]:
					continue
				rs.append(aAtm.flat())
				if patterns != "strict":
					try:
						rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
					except:
						try:
							rs.append("REMARK     MATCH. %s\n" % (matchPatterns[aAtm.atmName()[:2]][rank]))
						except:
							pass
				else:
					rs.append("REMARK     RESMATCH. %s\n" % aItem[0])
		if len(rs) > 2:
			rs.append("END    \n")
			grs.append(atmList(rs))
	return grs

def makeHName(rName, index, norm):
	"""
	makeHName() return the hydrogen name for the H in position "index" of the
	entry "rName" of globals dictionnary relativ to the specified "norm"
	@author: F.Briand
	@param rName: residue name
	@param index: index of the hydrogen in the "rName" line of "normHNames
				  dictionnaries
	@param norm: formating norm, "IUPAC" or "PDB"
	@return: a 4 characters string corresponding to the hydrogen name formated
			 with "norm" rules
	"""
	hName = normHNames[norm][1][rName][index]
	if norm == 'lIUPAC':
		if len(hName[0]) == 2:
			if not hName[1] or not hName[2]:
				return " %s" % hName[0] + str(hName[1]).replace("0"," ")
			else:
				return "%s" % hName[0] + "".join([str(hName[1]), str(hName[2])])
		if len(hName[0]) == 1:
			return " %s" % hName[0] + (str(hName[1])+str(hName[2])).replace("0"," ")
	elif norm == 'lPDB':
		if len(hName[0]) == 1:
			return str(hName[1]).replace("0"," ")+hName[0]+str(hName[2]).replace("0"," ")+' '
		else:
			return str(hName[1]).replace("0"," ")+hName[0]+str(hName[2]).replace("0"," ")


## ========================================
## Protein specific tools
## y = protein(x.chn("A"))
## ========================================
class protein(PDB):
	"""
	class protein
	This models protein
	"""
	def __init__(self, data, chId = "", model = 1, hetSkip = 0, verbose = 0):
		"""
		PDB.__init__ determine the type of data to initialize
		@param data: an instance
		@param chId: a chain Id
		@param model: the number of the model you want to set as working model (number 1 by default)
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@return: none
		"""
		if data != "":
			if isinstance(data,PDB):
				if verbose:
					print >> sys.stdout, "# protein init from PDB"
				self.atms = data.atms
				self.rt   = data.rt
				self.nFrg    = 0
				self.resTypes(verbose)
				self.frgs = []
				self.chns  = data.chns

			elif isinstance(data,types.ListType):
				if verbose:
					print >> sys.stdout, "# protein init from listType"
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
##				self.atms = []
##				for aLine in data:
##					self.atms.append(aLine)
##				self.atms = data
				self.resTab(verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = []
			elif isinstance(data,atmList):
				if verbose:
					print >> sys.stdout, "# protein init from atmList"
				self.atms = []
				for aLine in data.atms:
					self.atms.append(aLine)
##				self.atms = data
				PDB.resTab(self,verbose)
				self.resTypes(verbose)
				self.nFrg    = 0
				self.frgs = []
				self.chns  = ""
				self.id    = ""
				self.dbref = []
			elif isinstance(data,types.StringType):
				if verbose:
					print >> sys.stdout, "# protein init from string"
				self.atms  = []
				self.info  = []
				self.seq   = []
				self.seq3D = []
				self.ss    = []
				self.s2    = []
				self.id    = ""
				self.dbref = []
				self.chns  = ""
				self.nFrg    = 0
				self.frgs = []
				self.nModel = 0
				PDB.__init__(self,data, chId, model, hetSkip, verbose)
				## self.load(data, chId, hetSkip, verbose)
				## self.setModel(model, verbose)
				self.resTab(verbose)
				self.resTypes(verbose)

	def resTypes(self, verbose = 0):
		"""
		PDB.resTypes
		@return: a list of restypes of the protein
		"""
		self.tpe = []
		unres = []
		for aRes in range(0,len(self.rt) -1):
			aAtm = self.rt[aRes][0]
			if AA3.count(aAtm.resName()) != 0:
				idex = AA3.index(aAtm.resName())
				self.tpe.append(idex)
			else:
				if unres.count(aAtm.resName()) == 0:
					if verbose:
						print >> sys.stderr, "# Unknown residue type (3): ",aAtm.resName()
					unres.append(aAtm.resName())
				self.tpe.append(-1)

	def frgList(self):
		"""
		PDB.frgList()
		@return: a list of the fragments of the PDB if some geometric inconsistencies are
		detected and the numbers of fragments
		@note: the detection of fragments is based on the NC distance
		"""
		res = []
		theBB = self.BB()
		oriRes = 0
		nFrg = 0
		for aRes in range(1,len(theBB)):
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		self.frgs = res
		return nFrg, res

	def nFrgs(self):
		"""
		nFrgs
		@return: the number of fragments
		@note: the detection of fragments is based on the NC distance
		"""
		res = []
		theBB = self.BB()
		oriRes = 0
		nFrg = 0
		for aRes in range(1,len(theBB)):
			if self.tpe[aRes-1] == -1 and atmList(theBB[aRes-1].atms).theAtm("C") == []:
				continue
			if atmList(theBB[aRes-1].atms).theAtm("C") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
				continue
			Cx,Cy,Cz = atmList(theBB[aRes-1].atms).theAtm("C").xyz()
			if self.tpe[aRes] == -1 and atmList(theBB[aRes].atms).theAtm("N") == []:
				continue
			if atmList(theBB[aRes].atms).theAtm("N") == []:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes)
				oriRes = aRes+1
				res.append(lRes)
				continue
			Nx,Ny,Nz = atmList(theBB[aRes].atms).theAtm("N").xyz()
			aDist = distance(Cx,Cy,Cz,Nx,Ny,Nz)
			if aDist > 1.7:
				nFrg = nFrg + 1
				lRes = []
				lRes.append(oriRes)
				lRes.append(aRes-1)
				oriRes = aRes
				res.append(lRes)
		lRes = []
		lRes.append(oriRes)
		lRes.append(len(theBB) - 1)
		res.append(lRes)
		nFrg = nFrg + 1
		self.nFrg = nFrg
		return nFrg

	def trace(self,fname = "", chId = "", hetSkip = 0, altSel = " "): ### /!\ Dont work
		"""
		PDB.Trace:
		@param fname: the file name.
		@param hetSkip: does the file skip the hetero atoms (No:0 by default)
		@return: an atom list of the CA
		"""
		self.resTab(0)
		res = []
		for aRes in range(0,len(self.rt) -1):
			for aAtm in range(self.rt[aRes][0], self.rt[aRes+1][0]):
				aLine = self.atms[aAtm]
				if (hetSkip == 2) and (AA3STRICT.count(atmLine(aLine).resName()) == 0):
					#print "HETPEP : ", atmLine(aLine).resName()
					break
				if hetSkip and (AA3.count(atmLine(aLine).resName()) == 0):
					#print "HET : ", atmLine(aLine).resName()
					break
				#print atmLine(aLine).resName(), "Checking for CA"
				if atmLine(aLine).atmName() == "CA":
					res.append(aLine)
					break
##		for aLine in self.atms:
##			if hetSkip and AA3.count(atmLine(aLine).resName()) == 0:
##				continue
##			if atmLine(aLine).atmName() == "CA":
##				res.append(aLine)
##		print "trace :", res.__class__
		return atmList(res)

##	def chis(self):
##		"""
##		protein.chis()
##		"""
##		res = []
##		for aRes in range(0,len(self)):
##			res.append(self[aRes].chis())
##		return res

	def outSeq(self, label, hetSkip = 0, verbose = 0): ### /!\ Dont work
		"""
		protein.outseq()
		@param hetSkip: does the file skip the het (No:0 by default)
		@return: none, it prints the trace of the protein sequence
		"""
		#print self
		#print hetSkip
		theTrace = self.trace("","",hetSkip, verbose)
		#print theTrace
		#sys.exit(0)
		seq = protein(theTrace).aaseq()
		print >> sys.stdout, "> ",label,len(seq)
		while len(seq) > 0:
			print >> sys.stdout, seq[:80]
			seq = seq[80:]

	def outRawSeq(self, hetSkip = 0, verbose = 0):
		"""
		protein.outrawseq()
		@param hetSkip: does the file skip the het (No:0 by default)
		@return: none, it prints the trace of the protein sequence
		"""
		#print self
		#print hetSkip
		theTrace = self.trace("","",hetSkip, verbose)
		#print theTrace
		#sys.exit(0)
		seq = protein(theTrace).aaseq()
		print >> sys.stdout, seq

	def aaseq(self, verbose = 0):
		"""
		PDB.aaseq()
		@return: the sequence of residues present in the PDB file, having coordinates.
		@note: Converts non standard amino-acids to equivalent standard amino-acid.
		"""
		res = ""
		unres = []
		for aRes in self:
			if AA3STRICT.count(aRes[0].resName()):
				res = res + AA1[AA3STRICT.index(aRes[0].resName())]
			elif AA3.count(aRes[0].resName()):
				if verbose:
					print >> sys.stdout, "# Unfrequent residue type: ",aRes[0].resName()
				if aRes[0].resName() == "MSE": # seleno MET
					res = res+"M"
				elif aRes[0].resName() == "CSE": # seleno CYS
					res = res+"C"
				elif aRes[0].resName() == "FGL": # Formyl GLY
					res = res+"C"
				elif aRes[0].resName() == "CEA": # SHYDROXY-CYS
					res = res+"C"
				elif aRes[0].resName() == "TPQ": # 2,4,5-TRIHYDROXYPHE
					res = res+"Y"
				elif aRes[0].resName() == "CGU": # GAMMA-CARBOXY-GLU
					res = res+"E"
				elif aRes[0].resName() == "MHO": # Hydroxy-MET
					res = res+"M"
				elif aRes[0].resName() == "IAS": # BETA-CARBOXY ASP
					res = res+"D"
				elif aRes[0].resName() == "TYS": # SULFONATED TYROSINE
					res = res+"Y"
				else:
					res = res+'X'
			else:
				if unres.count(aRes[0].resName()) == 0:
					unres.append(aRes[0].resName())
		if verbose:
			print >> sys.stderr, "# Unknown residue type (2):",
			for iRes in unres:
				print >> sys.stderr, " %s " % iRes,
				print >> sys.stderr

		return res

	def frg(self,whatFrg, frgs = []):
		"""
		PDB.frg
		@param whatFrg: the fragment of the protein we need
		@param frgs: a list of fragments
		@return: the atoms of the fragment "whatFrg"
		"""
		if frgs == [] and self.frgs == []:
			self.nFrg, self.frgs = self.frgList()
		return protein(self[self.frgs[whatFrg][0]:self.frgs[whatFrg][1]+1])

	def hasAltAtms(self,verbose):
		"""
		PDB.hasAltAtms
		This will return 2 values consisting of ``Yes'' or ``No''. The first answers the question: does some x amino-acid backbone atoms have alternate coordinates (as specified in PDB files). The second anwsers the corresponding question for side chains.
		@return: Does the file has BBaltAtm or SCAltAtm? (Yes/No for each)
		"""
		BBAltAtm = False
		SCAltAtm = False
		for aLine in self.atms:
			if aLine[16] != ' ':
				isAlt = 1
				if string.count(string.digits,aLine[12]):
					isAlt = 0
				if aLine[12] == ' ' and aLine[13] == 'H':
					isAlt = 0
				if isAlt == 0:
					continue
				theAtmTpe = string.split(aLine[12:15])[0]
				if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
					BBAltAtm = True
				else:
					SCAltAtm = True
		return BBAltAtm, SCAltAtm

	def altAtmsResList(self,verbose):
		"""
		PDB.altAtmsResList
		This function is related to ``hasAltAtms''. It will, for backbone and side chains return the number of residues having backbone alt coordinates followed by a string containing the information about these residues (in a format similar to that described for SCatmMiss. The two first values concern backbone, the two next side chains.
		@return: nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm:
			- nBBAltAtm: number of Back Bones in alternate atoms
			- BBAltAtm: the Back Bones in alternate atoms
			- nSCAltAtm: number of side chain in alternate atoms
			- SCAltAtm: side chain in alternate atoms
		"""
		nBBAltAtm = 0
		nSCAltAtm = 0
		BBAltAtm  = ""
		SCAltAtm  = ""
		for aPos in range(0,len(self)):
			curRes = self[aPos]
			for aLine in curRes.atms:
				if aLine[16] != ' ':
					isAlt = 1
					if string.count(string.digits,aLine[12]):
						isAlt = 0
					if aLine[12] == ' ' and aLine[13] == 'H':
						isAlt = 0
					if isAlt == 0:
						continue
					theAtmTpe = string.split(aLine[12:15])[0]
					res    = aLine.resName()
					resNum = aLine.resNum()
					icode  = aLine.icode()
					lbl    = aLine.chnLbl()
					if icode == ' ':
						icode = ''
					if lbl == ' ':
						lbl = ''
					Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
					if theAtmTpe == "CA" or theAtmTpe == "N" or theAtmTpe == "C" or theAtmTpe == "C":
						nBBAltAtm = nBBAltAtm + 1
						BBAltAtm = BBAltAtm + Res
						break
					else:
						nSCAltAtm = nSCAltAtm + 1
						SCAltAtm = SCAltAtm + Res
						break
		return nBBAltAtm, BBAltAtm, nSCAltAtm, SCAltAtm

	def hasAllBBAtms(self,verbose):
		"""
		hasAllBBatms checks if all BB atoms are present
		@return: the position of the BB atoms missing
		"""
		CAWarning = 0
		CWarning  = 0
		OWarning  = 0
		NWarning  = 0
		residuNameMissing=[]
		cp=0
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			if aRes.Npos() == None:
				if aPos == 0:
					NWarning  = 1
				elif aPos == len(self) - 1:
					if NWarning < 1:
						NWarning  = 1
				else:
					NWarning  = 2
					cpt=1
					residuNameMissing.append(aPos)
			if aRes.CApos() == None:
				if aPos == 0:
					CAWarning  = 1
				elif aPos == len(self) - 1:
					if CAWarning < 1:
						CAWarning  = 1
				else:
					CAWarning  = 2
					if cp==0:
						cp=1
						residuNameMissing.append(aPos)
			if aRes.Cpos() == None:
				if aPos == 0:
					CWarning  = 1
				elif aPos == len(self) - 1:
					if CWarning < 1:
						CWarning  = 1
				else:
					CWarning  = 2
					if cp==0:
						cp=1
						residuNameMissing.append(aPos)
			if aRes.Opos() == None:
				if aPos == 0:
					OWarning  = 1
				elif aPos == len(self) - 1:
					if OWarning < 1:
						OWarning  = 1
				else:
					OWarning  = 2
					if cp==0:
						cp=1
						residuNameMissing.append(aPos)
			cp=0
		### /!\ Yes/Ext/No pertinence ?
		if OWarning == 2 or NWarning == 2 or CAWarning == 2 or CWarning == 2:
			BBAtmMiss = "Yes"
		elif OWarning == 1 or NWarning == 1 or CAWarning == 1 or CWarning == 1:
			BBAtmMiss = "Ext"
		else:
			BBAtmMiss = "No"
		return BBAtmMiss,residuNameMissing

	# Check if BB peptidic geometry is correct (distance)
	def geomCheck(self,verbose):
		"""
		PDB.geomCheck()
		@return: Is the BB peptidic geometry (distance) correct? (OK/Poor/Bad)
		@note: THIS WILL NOT DETECT FRAGMENTS. IF MANY, THE GAPS ARE IGNORED AND DO NOT RESULT IN "Bad" RETURN. \n
		This allows to scan that all the fragments are correct at once.
		"""
		aN = None
		aC = None
		Cx, Cy, Cz = 0., 0., 0.
		BBGeoOK = "Ok"
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			aN = aRes.Npos()
			if aN != None:
				# Nx, Ny, Nz = atmLine.atmCrds(aRes[aN])
				Nx, Ny, Nz = aRes[aN].xyz()
			if aC != None:
				aDist = distance(Nx, Ny, Nz, Cx, Cy, Cz)
				if aDist > 1.50 and aDist < 3.:
					if verbose:
						print >> sys.stdout, "# Poor peptidic bond of ",aDist," for ", resName(theChain[aC]), resNum(theChain[aC]), resName(theChain[aN]), resNum(theChain[aN])
					if BBGeoOK == "Ok":
						BBGeoOK = "Poor"
				elif aDist > 3.:
					if verbose:
						print >> sys.stdout, "# Bad peptidic bond  of ",aDist," for :", resName(theChain[aC]), resNum(theChain[aC]), resName(theChain[aN]), resNum(theChain[aN])
					BBGeoOK = "Bad"
			aC  = aRes.Cpos()
			if aC != None:
				# Cx, Cy, Cz =atmLine.atmCrds(aRes[aC])
				Cx, Cy, Cz = aRes[aC].xyz()
		return BBGeoOK

	# Check if BB peptidic geometry is correct (distance)
	# /!\ bad description ...
	def traceCheck(self,hetSkip = 0, verbose = 0):
		"""
		PDB.traceCheck check if BB peptidic geometry is correct (distance)
		@param hetSkip: does the file skip the het (No:0 by default)
		@return: traceOK (OK/bad), tracePB (residues with bad geometry), nCISPRO (number of cis prolines), CISPRO (cis prolines), nCISPep (number of cis peptides), CISPep (cis peptides),CisWarning (CisPRO/CisPEP), hasCisPRO(Yes/No), hasCisPEP(Yes,No)
		"""
		theTrace = self.trace("","",hetSkip, verbose)
		CisWarning = "None"
		hasCisPRO = "No"
		hasCisPEP = "No"
		traceOK = "Yes"
		nCISPRO = 0
		nCISPep = 0
		CISPRO  = ""
		CISPep  = ""
		for aRes in range(1,len(theTrace)):
			try:
				x1, y1, z1 = theTrace[aRes - 1].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:", theTrace[aRes - 1]
				return CisWarning,"No"
			try:
				x2, y2, z2 = theTrace[aRes].xyz()
			except ValueError:
				if verbose:
					print >> sys.stderr, "# Sorry: incorrect ATOM line format for:", theTrace[aRes]
				return CisWarning,"No"
			aDist = distance(x1, y1, z1, x2, y2, z2)
			if aDist < 3.60: # CIS peptide
				res    = atmLine(self[aRes].atms[0]).resName()
				resNum = atmLine(self[aRes].atms[0]).resNum()
				icode  = atmLine(self[aRes].atms[0]).icode()
				lbl  = atmLine(self[aRes].atms[0]).chnLbl()
				if icode == ' ':
					icode = ''
				if lbl == ' ':
					lbl = ''
				Res=res+"_"+lbl+"_"+str(resNum)+"_"+icode+" "
				if CisWarning == "None":
					CisWarning = "CISPRO"
				if theTrace[aRes][17:20] != "PRO": # CIS PROLINES
					CisWarning = "CISPEP"
					hasCisPEP  = "Yes"
					nCISPep = nCISPep + 1
					CISPep  = CISPep + Res
				else:
					hasCisPRO  = "Yes"
					nCISPRO = nCISPRO + 1
					CISPRO  = CISPRO + Res
			if aDist > 4.10: # bad geometry
				traceOK = "No"
				if verbose:
					print >> sys.stdout, "# Bad Trace for ",theTrace[aRes-1]
		return CisWarning, hasCisPRO, hasCisPEP, traceOK, nCISPRO, CISPRO, nCISPep, CISPep

	def BBAngles(self,aRes = -1000):
		"""
		PDB BBangles calculate phi psi ome of 2 consecutives residues
		@param aRes: a residue number, -1000 by default in this case it will calculate all phi psi ome of BB angles.
		@return: angles phi psi ome, or a list of them if aRes is not set by user.
		"""
		res = []
		if aRes == -1000:
			rFrom = 0
			rTo = len(self)
		else:
			rFrom = aRes
			rTo = aRes+1
		for aPos in range(rFrom,rTo):
			phi = -1000.
			psi = -1000.
			ome = -1000.
			if aPos > 0:
				OK = 1
				aAtm = self[aPos-1].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos-1].theAtm("C").xyz()
					b = self[aPos].theAtm("N").xyz()
					c = self[aPos].theAtm("CA").xyz()
					d = self[aPos].theAtm("C").xyz()
					phi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("N").xyz()
					b = self[aPos].theAtm("CA").xyz()
					c = self[aPos].theAtm("C").xyz()
					d = self[aPos+1].theAtm("N").xyz()
					psi = apply(dihedral,a+b+c+d)
			if aPos < len(self) - 1:
				OK = 1
				aAtm = self[aPos].theAtm("CA")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos].theAtm("C")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("N")
				if aAtm == []:
					OK = 0
				aAtm = self[aPos+1].theAtm("CA")
				if aAtm == []:
					OK = 0
				if OK:
					a = self[aPos].theAtm("CA").xyz()
					b = self[aPos].theAtm("C").xyz()
					c = self[aPos+1].theAtm("N").xyz()
					d = self[aPos+1].theAtm("CA").xyz()
					ome = apply(dihedral,a+b+c+d)
			res.append([phi,psi,ome])
		return res

	def SGList(self):
		"""
		PDB.SGList
		@return: a list of all the coordinates of the gamma sulfur
		"""
		SGList = []
		for aPos in range(0,len(self)):
			aRes = self[aPos]
			if aRes[0].resName() == "CYS":
				lSGList = []
				for aAtm in range(0,len(aRes.atms)):
					if atmLine(aRes.atms[aAtm]).atmName() == "SG":
						lSGList.append(atmLine(aRes.atms[aAtm]).xyz())
				if lSGList != []:
					SGList.append(lSGList)
		return SGList

	def nSSIntra(self):
		"""
		nSSIntra()
		@return: the number of SSbonds in the PDB instance
		"""
		nSSBond = 0
		aSGList = self.SGList()
		for aRes1 in range(0,len(aSGList)):
			for aSG1 in range(0,len(aSGList[aRes1])):
				for aRes2 in range(aRes1+1,len(aSGList)):
					for aSG2 in range(0,len(aSGList[aRes2])):
						if apply(distance, aSGList[aRes1][aSG1]+aSGList[aRes2][aSG2]) < 2.35:
							nSSBond = nSSBond + 1
							break
		return nSSBond

	def BB(self):
		"""
		protein.BB()
		@return: a protein with the backbone atoms only
		"""
		res = []
		for aLine in self.atms:
			theName = aLine.atmName()
			if BBATMS.count(theName) > 0:
				res.append(aLine)
		return  PDB(res)

	def SC(self):
		"""
		protein.SC()
		@return: a protein with the side-chains atoms only
		"""
		res = []
		for aLine in self.atms:
			theName = atmLine(aLine).atmName()
			if BBATMS.count(theName) == 0:
				res.append(aLine)
		return  PDB(res)


def PDBList(input = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	This is to organize the iterative treatment of PDB instances.
	@param input: a list (of lines), or a file (list of lines).
	@param altCare: does the file contains the alternate atoms (No:0 by default)
	@param OXTCare: does the file contains the information of each line (No:0 by default)
	@param hetSkip: does the file skip the hetero atoms (No:0 by default)
	@type input: each line can specify:
		 a local file, it may contain a multi PDB separated with HEADER / END lines
		 aPDB Id compatible with the PDB class
		 an url
	@return: a list of PDB instances.
	@note: The list length is the number of lines of the input
	"""
	from urllib import urlretrieve
	rs = []
	if input == None:
		return rs
	# Try as if input already specifies the data itself
	# try reading PDB(s) as local PDB file
	# String can be:
	#    file
	#    PDBid
	#    url
	# file can be:
	#    PDBList
	#    id/url list
	# url can be:
	#    PDBList
	#    id/url list
	if isinstance(input,types.StringType):
		if verbose:
			print >> sys.stdout, "# PDBList: Input is string: %s" % input
		try:
			if verbose:
				print >> sys.stdout, "# PDBList: Attempting local multiPDB file"
			open(input).close()
			rs = fileInput(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if rs != []:
				try:
					id = rs[0].rt
					return rs
				except:
					pass
		except:
			pass
		# direct PDBid
		if verbose:
			print >> sys.stderr, "# PDBList: Not a local multiPDB file"
		try:
			if verbose:
				print >> sys.stdout, "# PDBList: Attempting PDBid", input
			x = PDB(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (x != None) and (len(x)):
				return [x]
		except:
			if verbose:
				print >> sys.stderr, "# PDBList: Not a PDB"
			try:
				if verbose:
					print >> sys.stdout, "# PDBList: Attempting PDBid at pdb.org"
				file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(input[:4],input[:4]))
				x = PDB(file, chId = string.replace(input[4:]," ",""), hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
				if (x != None) and (len(x)):
					x.id = input
					return [x]
			except:
				pass
			try:
				if verbose:
					print >> sys.stdout, "# PDBList: Attempting URL to direct multiPDB or PDB"
				file, log = urlretrieve(input)
				if verbose:
					print >> sys.stdout, "#", file
				x = fileInput(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
				if len(x) == 0:
					if verbose:
						print >> sys.stderr, "# PDBList: Not a direct multiPDB URL"
				elif x != None:
					return x
			except:
				pass
##			try: # here we try explicit url if pdb.org url did not returned exception
##				if verbose:
##					print "Attempting URL"
##				file, log = urlretrieve(input)
##				x = PDBList(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
##				if x != None:
##					return x
##			except:
##				pass
	# a list of atoms
	elif isinstance(input,types.ListType):
		if verbose:
			print >> sys.stdout, "# PDBList: Input is list"
		try:
			rs = parseInput(input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if rs != []:
				return rs
		except:
			pass
	if verbose:
		print >> sys.stderr, "# PDBList: Input might be complex input"
	# From here, we consider we have indirections to the data
	# We expect 1 dataset each valid line
	rs = []
	# try file of Ids or URLs
	if isinstance(input,types.StringType):
		try:
			input = open(input).readlines()
		except:
			try:
				if verbose:
					print >> sys.stdout, "# PDBList: Attempting URL to list of proteins"
				file, log = urlretrieve(input)
				input = open(file).readlines()
				if verbose:
					print >> sys.stdout, "#", input ### /!\
			except:
				pass
	# try list of Ids or URLs
	if isinstance(input,types.ListType):
		for aInput in input:
			x = None
			# remove \n, \r if any
			try:
				aInput = string.replace(aInput,"\n","")
				aInput = string.replace(aInput,"\r","")
				if aInput == "":
					continue
			except:
				pass
			if aInput[0] == "#":
				continue
			if isinstance(aInput,PDB):		 # already a PDB instance
				rs.append(aInput)
				continue
			elif isinstance(aInput,types.StringType):  # read file from disk
				if verbose:
					print >> sys.stdout, "# PDBList: Considering: ",aInput
				try:
					if verbose:
						print >> sys.stdout, "# PDBList: Trying local file"
					lines = open(aInput).readlines()
					x = parseInput(lines, Id = aInput, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
				except:
					try:
						if verbose:
							print >> sys.stdout, "# PDBList: Trying PDBId"
						x = PDB(aInput, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
						x.id = aInput
						if (x != None) and len(x):
							x = [x]
						else:
							x = None
					except:
						pass
					# if len(aInput) == 4:
					# if 1:
					if x == None:
						if verbose:
							print >> sys.stdout, "# PDBList: Trying PDB entry (at PDB), ChnIds: \"%s\"" % aInput[4:]
						try:
							file, log = urlretrieve("http://www.rcsb.org/pdb/cgi/export.cgi/%s.pdb?format=PDB&pdbId=%s&compression=None" %(aInput[:4],aInput[:4]))
							x = PDB(file, chId = string.replace(aInput[4:]," ",""), hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
							if verbose:
								print >> sys.stdout, "# PDBList: Completed PDB at pdb.org"
							if (x != None) and len(x):
								x.id = aInput
								x = [x]
							else:
								x = None
						except:
							pass
					if x == None and (string.lower(aInput[:4]) == "http"):
						try:
							if verbose:
								print >> sys.stdout, "# PDBList: Trying URL file"
							file, log = urlretrieve(aInput)
							x = PDBList(file, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
							if (x == []):
								x = None
						except:
							pass
			if x != None:
				rs += x
	del urlretrieve
	return rs

def fileInput(input = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	fileInput: to read a multiPDB file from disk.
	@param input: the file to read
	@param altCare: does the file contains the alternate atoms (Yes:1 by default)
	@param OXTCare: does the file contains the information of each line (No:0 by default)
	@param hetSkip: does the file skip the hetero atoms (No:0 by default)
	@return: a parsed PDB instance of the inputs contents
	"""
	rs = []
	try:
		inputs = open(input).readlines()
	except:
		return rs
	return parseInput(inputs, Id = input, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)

def parseInput(inputs, Id = None, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	parseInput:
	We have a list of PDBlines. We parse them and return a list of PDBs.
	@param inputs: the PDB instance to parse
	@param Id: the Id of the protein
	@param altCare: does the file contains the alternate atoms (Yes:1 by default)
	@param OXTCare: does the file contains the information of each line (No:0 by default)
	@param hetSkip: does the file skip the hetero atoms (No:0 by default)
	@return: a parsed PDB instance of the inputs contents
	"""
	rs = []
	count = 0
	hcount = 0
	ecount = 0
	for i in inputs:
		if string.count(i[:6],"HEADER"):
			hcount +=1
		if string.count(i[:3],"END"):
			ecount +=1
	if verbose:
		print >> sys.stdout, "# parseInput: detected %d HEADER lines" % hcount
	if hcount < 2:
		if verbose:
			print >> sys.stdout, "# Got",hcount,"HEADER"
		try:
			x = PDB(inputs, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (Id != None) and (x.id == "unkwn"):
				x.id = Id
		except:
			return []
		if (x != None) and len(x):
			return [x]
		return []
	else:
		hcount = 0
		for i in range(0,len(inputs)):
			if string.count(inputs[i][:6],"HEADER"):
				if hcount != 0:
					try:
						x = PDB(inputs[hcount-1:i], hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
						if (x != None) and len(x):
							count += 1
							if x.id == "unkwn":
								if (Id != None):
									x.id = "%s_%d" % (Id,count)
								else:
									x.id = "%s_%d" % ("unkwn",count)
							rs.append(x)
					except:
						pass
				hcount = i+1 # +1 to avoid 0 again
		try:
			x = PDB(inputs[hcount-1:i], hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
			if (x != None) and len(x):
				count += 1
				if x.id == "unkwn":
					if (Id != None):
						x.id = "%s_%d" % (Id,count)
					else:
						x.id = "%s_%d" % ("unkwn",count)
				rs.append(x)
		except:
			pass
	return rs

def outPDBList(pdbList, outName = "", initMode =  "w", altCare = 0, altLbl = "", OXTCare = 0, hetSkip = 0, fmode = "w", header = 1, ter = 1, end = 1, info = 0, verbose = 0):
	"""
	outPDBList
	@param pdbList: a list of PDB
	@param initMode: mode of creating file, "w" write by default
	@param outName: the name of file where the PDB list will be written
	@param fmode: mode of opening file, "w" write by default
	@return: none, it writes the PDB list in a file
	"""
	for i in range(0,len(pdbList)):
		if i == 0:
			pdbList[i].out(outName = outName, fmode=initMode, altCare = altCare, altLbl = altLbl, OXTCare = OXTCare, hetSkip = hetSkip, header = header, ter = ter, end = end, info = info, verbose = verbose)
		else:
			pdbList[i].out(outName = outName, fmode="a", altCare = altCare, altLbl = altLbl, OXTCare = OXTCare, hetSkip = hetSkip, header = header, ter = ter, end = end, info = info, verbose = verbose)

def PDBSumHeaders(what):
	"""
	PDBSumHeaders
	@param what: the word that you want to search on EBI website
	@return: the PDB ids corresponding to your research on what
	"""
	url = 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/SearchHeaders.pl'
	values = {'string' : 'hydrolase', 'toolbar' : 'biobar'}
	values['all'] = "TRUE"
	values['string'] = what
	data = urllib.urlencode(values)
	req = urllib2.Request(url, data)
	response = urllib2.urlopen(req)
	the_page = response.read()
	# return the_page
	the_lines =  string.split(html2text(the_page),"\n")
	# print len(the_lines)
	# sys.exit(0)
	rs = []
	for i in the_lines:
		if string.count(i,"/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode="):
			rs.append(string.split(i,"=")[1].encode())
	return rs

def PDBListFromPDBSum(what, hetSkip = 0, altCare = 0, OXTCare = 0, verbose = 0):
	"""
	PDBListFromPDBSum
	@param what: the word that you want to search on EBI website
	@param altCare: does the file contains the alternate atoms (No:0 by default)
	@param OXTCare: does the file contains the information of each line (No:0 by default)
	@param hetSkip: does the file skip the hetero atoms (No:0 by default)
	@return: a list of all the PDB ids corresponding to your research on what
	"""
	rs = PDBSumHeaders(what)
	if verbose:
		print >> sys.stdout, "#", what
		print >> sys.stdout, "#", rs
	pdbrs = PDBList(rs, hetSkip = hetSkip, altCare = altCare, OXTCare = OXTCare, verbose = verbose)
	return pdbrs

def subPDB(pdb, seedSeq):
	"""
	Given a PDB instance, return a part corresponding to the seed sequence
	@author: P. Tuffery
	@param pdb : a PPDB instance
	@param seedSeq: the sequence to fetch (a string)
	@return: a PDB instance corresponding to the seeSeq
		 or None if the seedSeq is not found
	"""
	aas = pdb.aaseq()
	if aas.count(seedSeq):
		pos = aas.index(seedSeq)
		return pdb[pos: pos+len(seedSeq)]
	return None


class remark350():
	"""
	class remark350
	this models the content of REMARK 350 fields
	"""
	def __init__(self, input = "", verbose = 0):
		"""
		remark350.__init__ initialize datas
		@author: F.Briand
		@param input: list of REMARK 350 lines
		"""
		self.biomolecule = []
		self.chains = []
		self.biomt = {}
		self.txt = input
		for line in input:
			if line[11:22] == "BIOMOLECULE":
				m = re.split("\D+",line[10:])
				for biomol in m:
					if biomol and biomol not in self.biomolecule:
						self.biomolecule.append(biomol)
			elif line[34:40] == "CHAINS":
				self.chains = list(set(re.split("\W+",line[41:])))
				self.chains.sort()
				self.chains = "".join(self.chains).strip()
			elif line[13:18] == "BIOMT":
				m = re.match("^REMARK 350   BIOMT[1-3] *[0-9]+ *([0-9\-]+\.[0-9]+) *([0-9\-]+\.[0-9]+) *([0-9\-]+\.[0-9]+) *([0-9\-]+\.[0-9]+).*$",line)
				try:
					self.biomt[line[22]].append(m.groups())
				except KeyError:
					self.biomt[line[22]] = [m.groups()]

	def __repr__(self):
		"""
		remark350.__repr__
		@author: F.Briand
		@return: return the content of instance (list of lines REMARK 350)
		"""
		txt = ""
		for line in self.txt:
			txt += str(line)
		return txt

	def __len__(self):
		"""
		remark350.__len__
		@author: F.Briand
		@return: number of REMARK 350 lines
		"""
		return len(self.txt)

	def __getitem__(self, rmkPos):
		"""
		remark350.__getitem__ return information calling self[rmkPos]
		@author: F.Briand
		@param rmkPos:
			- "biomol"	: list of biomolecules
			- "biomt"	: translation/rotation matrix
			- "chains"	: list of chains
			- n (int)	: n'ieme line of REMARK 350 fields
		@return:
		"""
		if rmkPos == "biomol":
			return self.biomolecule
		elif rmkPos == "biomt":
			return self.biomt
		elif rmkPos == "chains":
			return self.chains
		else:
			try:
				return self.txt[int(rmkPos)]
			except ValueError:
				raise ValueError, "List indices must be integers."

	def __getslice__(self, ffrom = 0, tto = None):
		"""
		remark350.__getslice__
		@author: F.Briand
		@return: a slice of REMARK 350 lines : self[ffrom:tto]
		"""
		return self.txt[ffrom:tto]

def isPDB(fname):
        """
        check if fname corresponds to a valid PDB
        """
        try:
                x = PDB(fname)
        except:
                msg = "Sorry: file %s does not seem a valid PDB." % fname
                return False, msg
        if not len(x):
                msg = "Sorry: file %s does not seem a valid PDB." % fname
                return False, msg
        return True, ""
