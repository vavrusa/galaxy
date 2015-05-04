
import string
import sys
import os
import copy
import math
import gzip
## import bz2

def simpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = open(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	return lines

def load(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = open(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def gzsimpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = gzip.GzipFile(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	return lines

def gzload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = gzip.GzipFile(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def bz2simpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = bz2.BZ2File(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	return lines

def bz2load(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	lines = bz2.BZ2File(fname,'r').readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def zsimpleload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	cmd = "gunzip -c "+fname
	lines = os.popen(cmd).readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	return lines

def zload(fname, verbose = 0):
	try:
		os.stat(fname)
	except:
		raise IOError, "%s : failed to stat file !!" % fname
	cmd = "gunzip -c "+fname
	lines = os.popen(cmd).readlines()
	if verbose > 1:
		print >> sys.stdout, "# %s : Read %d lines" % (fname, len(lines))
	for i in range(0,len(lines)):
		lines[i] = string.split(lines[i])
	return lines

def gsimpleload(fname, verbose = 0):
	if fname[-3:] == ".gz":
		if verbose:
			print >> sys.stdout, "# Trying: %s as gz ..." % fname
		lines = gzsimpleload(fname, verbose = verbose)
##	elif fname[-2:] == ".bz2":
##		lines = bz2simpleload(fname, 0)
	elif fname[-2:] == ".Z":
		if verbose:
			print >> sys.stdout, "# Trying: %s as .Z ..." % fname
		lines = zsimpleload(fname, verbose = verbose)
	else:
		if verbose:
			print >> sys.stdout, "# Trying: %s as flat ..." % fname
		lines = simpleload(fname, verbose = verbose)
	return lines

def gload(fname, verbose = 0):
	if fname[-3:] == ".gz":
		lines = gzload(fname, 0)
##	elif fname[-2:] == ".bz2":
##		lines = bz2load(fname, 0)
	elif fname[-2:] == ".Z":
		lines = zload(fname, 0)
	else:
		lines = load(fname, 0)
	return lines

