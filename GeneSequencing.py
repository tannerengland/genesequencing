#!/usr/bin/python3
import math

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

		if len(seq1) > align_length:
			seq1 = seq1[:align_length]
		if len(seq2) > align_length:
			seq2 = seq2[:align_length]

		if len(seq1) > len(seq2):
			seq1, seq2 = seq2, seq1

		if not banded:
			score = self.NW(seq1, seq2)
		else:
			score = self.NWB(seq1, seq2)

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100;
		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}

	def NW(self, x, y):
		E = [[0] * (len(y)+1) for _ in range(len(x)+1)]
		for i in range(0, len(x)+1):
			E[i][0] = i * INDEL
		for j in range(0, len(y)+1):
			E[0][j] = j * INDEL
		for i in range(1, len(x)+1):
			for j in range(1, len(y)+1):
				E[i][j] = min(self.diff(i,j,x,y) + E[i-1][j-1], E[i][j-1] + INDEL, E[i-1][j] + INDEL)
		return E[len(x)][len(y)]

		# northwest, west, north

	def NWB(self, x, y):
		bandwidth = MAXINDELS * 2 + 1
		E = [[float('inf')] * (bandwidth) for _ in range(len(x)+1)]

		for col in range(0, MAXINDELS + 1):
			E[col][MAXINDELS - col] = col * INDEL
		for row in range(0, MAXINDELS + 1):
			E[0][MAXINDELS + row] = row * INDEL
		for row in range(1, len(x)+1):
			for col in range(0, bandwidth):
				# if (j <= (i + MAXINDELS)) and (j >= (i - MAXINDELS)):
				currIndex = row - MAXINDELS + col
				diffIndex = col - math.ceil(bandwidth/2)
				if currIndex >= 0:
					# if row - 1 >= 0 and col - 1 >= 0 and col + 1 <= bandwidth-1:
					west = float('inf')
					north = float('inf')
					northwest = float('inf')

					if row - 1 >= 0:
						northwest = self.diff(row,diffIndex,x,y) + E[row-1][col]

					if col - 1 >= 0:
						west = E[row][col - 1] + INDEL

					if row-1 >= 0 and col + 1 <= bandwidth-1:
						north = E[row-1][col+1] + INDEL

					# E[row][col] = min(self.diff(row,currIndex,x,y) + E[row-1][col], E[row][col-1] + INDEL, E[row-1][col+1] + INDEL)
					E[row][col] = min(west,north,northwest)
		return E[len(x)][math.floor(bandwidth/2)]


	def diff(self,i,j,x,y):
		if x[i-1] == y[j-1]:
			return MATCH
		else:
			return SUB