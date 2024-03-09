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

		swapped = False
		if len(seq1) > len(seq2):
			seq1, seq2 = seq2, seq1
			swapped = True

		if (len(seq2) - len(seq1)) > (MAXINDELS * 2 + 1):
			score = float('inf')
			alignment1 = "No Alignment Possible."
			alignment2 = "No Alignment Possible."
		else:
			if not banded:
				score,alignment1,alignment2 = self.NW(seq1, seq2)
			else:
				score,alignment1,alignment2 = self.NWB(seq1, seq2)

		if swapped == True:
			alignment1, alignment2 = alignment2, alignment1

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
# 		score = random.random()*100;
# 		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq1), align_length, ',BANDED' if banded else '')
# 		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
# 			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################
		# print(alignment1)
		# print(alignment2)

		return {'align_cost':score, 'seqi_first100':alignment1[:100], 'seqj_first100':alignment2[:100]}

	def NW(self, x, y):
		Pointers =  [[0] * (len(y)+1) for _ in range(len(x)+1)]
		E = [[0] * (len(y)+1) for _ in range(len(x)+1)]
		for i in range(0, len(x)+1):
			E[i][0] = i * INDEL
			Pointers[i][0] = "n"
		for j in range(0, len(y)+1):
			E[0][j] = j * INDEL
			Pointers[0][j] = "w"
		for i in range(1, len(x)+1):
			for j in range(1, len(y)+1):
				west = E[i][j - 1] + INDEL
				north = E[i - 1][j] + INDEL
				northwest = self.diff(i, j, x, y) + E[i - 1][j - 1]
				E[i][j] = min(northwest, west, north)
# i indicates insert, d indicates delete, o indicates either match or sub
				pointer = E[i][j]
				if (pointer == west):
					pointer = "w"
				elif (pointer == north):
					pointer = "n"
				else:
					pointer = "nw"

				Pointers[i][j] = pointer
		# horizontal is insert by dash in x, diagnol match or sub by value, vertical delete dash in y

		seq1,seq2 = self.align_characters(Pointers,x,y)
		# E[i][j] = min(self.diff(i,j,x,y) + E[i-1][j-1], E[i][j-1] + INDEL, E[i-1][j] + INDEL)
		return E[len(x)][len(y)],seq1,seq2

	def align_characters(self,Pointers,x,y):
		seqx = ""
		seqy = ""

		i = len(x)
		j = len(y)
		while i > 0 or j > 0:
			currPointer = Pointers[i][j]

			if currPointer == "w":
				seqx += "-"
				seqy += y[j-1]
				j = j - 1
			elif currPointer == "n":
				seqx += x[i-1]
				seqy += "-"
				i = i - 1
			elif currPointer == "nw":
				seqx += x[i-1]
				seqy += y[j-1]
				i = i - 1
				j = j - 1

		seqx = seqx[::-1]
		seqy = seqy[::-1]
		return seqx,seqy

	# horizontal is insert by dash in x, diagnol match or sub by value, vertical delete dash in y

	# northwest, west, north

	def NWB(self, x, y):
		bandwidth = MAXINDELS * 2 + 1

		Pointers =  [[0] * (bandwidth) for _ in range(len(x)+1)]
		E = [[float('inf')] * (bandwidth) for _ in range(len(x)+1)]

		for col in range(0, MAXINDELS + 1):
			E[col][MAXINDELS - col] = col * INDEL
			Pointers[col][MAXINDELS - col] = "n"
		for row in range(0, MAXINDELS + 1):
			E[0][MAXINDELS + row] = row * INDEL
			Pointers[0][MAXINDELS + row] = "w"
		for row in range(1, len(x)+1):
			for col in range(0, bandwidth):
				# if (j <= (i + MAXINDELS)) and (j >= (i - MAXINDELS)):
				currIndex = row - MAXINDELS + col
				diffIndex = col - MAXINDELS + row
				if currIndex >= 0 and currIndex <= len(y):
					# if row - 1 >= 0 and col - 1 >= 0 and col + 1 <= bandwidth-1:
					west = float('inf')
					north = float('inf')
					northwest = float('inf')

					if row - 1 >= 0 and diffIndex <= len(y):
						northwest = self.diff(row,diffIndex,x,y) + E[row-1][col]

					if col - 1 >= 0:
						west = E[row][col - 1] + INDEL

					if row-1 >= 0 and col + 1 <= bandwidth-1:
						north = E[row-1][col+1] + INDEL

					# E[row][col] = min(self.diff(row,currIndex,x,y) + E[row-1][col], E[row][col-1] + INDEL, E[row-1][col+1] + INDEL)
					E[row][col] = min(west,north,northwest)

					pointer = E[row][col]
					if (pointer == west):
						pointer = "w"
					elif (pointer == north):
						pointer = "n"
					else:
						pointer = "nw"

					Pointers[row][col] = pointer

		last_non_inf_index = None
		for i in range(len(E[len(x)]) - 1, -1, -1):
			if E[len(x)][i] != float('inf'):
				last_non_inf_index = i
				break

		seq1,seq2 = self.align_characters_banded(Pointers,x,y,last_non_inf_index)
		return E[len(x)][last_non_inf_index],seq1,seq2

	def align_characters_banded(self,Pointers,x,y,last_non_inf_index):
		seqx = ""
		seqy = ""
		i = len(x)
		j = last_non_inf_index
		currIndex = i - MAXINDELS + j
		while currIndex > 1 and currIndex <= len(y):
			currIndex = i - MAXINDELS + j
			diffIndex = j - MAXINDELS + i
			# if currIndex >= 0 and currIndex <= len(y):
			currPointer = Pointers[i][j]

			if currPointer == "w":
				seqx += "-"
				seqy += y[diffIndex-1]
				j = j - 1
			elif currPointer == "n":
				seqx += x[i-1]
				seqy += "-"
				i = i - 1
				j = j + 1
			elif currPointer == "nw":
				seqx += x[i-1]
				seqy += y[diffIndex-1]
				i = i - 1

		seqx = seqx[::-1]
		seqy = seqy[::-1]
		return seqx,seqy


	def diff(self,i,j,x,y):
		if x[i-1] == y[j-1]:
			return MATCH
		else:
			return SUB