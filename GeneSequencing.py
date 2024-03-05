#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		if not banded:
			editDist = self.NW(self, seq1, seq2)
		else:
			editDist = self.NWB(seq1, seq2)

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		score = random.random()*100;
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def NW(self, x, y):
		E = [[0] * len(y) for _ in range(len(x))]
		for i in range(len(x)):
			E[i][0] = i * INDEL
		for j in range(len(y)):
			E[0][j] = j * INDEL
		for i in range(1,len(x)):
			for j in range(1,len(y)):
				E[i][j] = min(self.diff(i,j,x,y) + E[i-1][j-1], E[i][j-1] + INDEL, E[i-1][j] + INDEL)
		return E[len(x)][len(y)]

	def NWB(self, x, y):
		E = [[0] * len(y) for _ in range(len(x))]

		return E[len(x),len(y)]

	def diff(self,i,j,x,y):
		if x[i] == y[j]:
			return MATCH
		else:
			return SUB