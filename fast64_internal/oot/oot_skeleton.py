import shutil, copy, mathutils, bpy, math, os, re
from bpy.utils import register_class, unregister_class

from ..utility import *
from .oot_utility import *
from .oot_constants import *
from .oot_model_classes import *
from ..panels import OOT_Panel

from ..f3d.f3d_parser import *
from ..f3d.f3d_writer import *
from ..f3d.f3d_material import TextureProperty, tmemUsageUI
from .oot_f3d_writer import *

# for debugging only
#from pprint import pprint

ootEnumBoneType = [
	("Default", "Default", "Default"),
	("Custom DL", "Custom DL", "Custom DL"),
	("Ignore", "Ignore", "Ignore"),
]

class OOTSkeleton():
	def __init__(self, name):
		self.name = name
		self.segmentID = None
		self.limbRoot = None
		self.hasLOD = False

	def createLimbList(self):
		if self.limbRoot is None:
			return []

		limbList = []
		self.limbRoot.getList(limbList)
		self.limbRoot.setLinks()
		return limbList

	def getNumDLs(self):
		if self.limbRoot is not None:
			return self.limbRoot.getNumDLs()
		else:
			return 0

	def getNumLimbs(self):
		if self.limbRoot is not None:
			return self.limbRoot.getNumLimbs()
		else:
			return 0

	def isFlexSkeleton(self):
		if self.limbRoot is not None:
			return self.limbRoot.isFlexSkeleton()
		else:
			return False

	def limbsName(self):
		return self.name + "Limbs"

	def toC(self):
		limbData = CData()
		data = CData()

		if self.limbRoot is None:
			return data

		limbList = self.createLimbList()
		isFlex = self.isFlexSkeleton()

		data.source += "void* " + self.limbsName() + "[" + str(self.getNumLimbs()) + "] = {\n"
		for limb in limbList:
			limbData.source += limb.toC(self.hasLOD)
			data.source += '\t&' + limb.name() + ',\n'
		limbData.source += '\n'
		data.source += "};\n\n"

		if isFlex:
			data.source += "FlexSkeletonHeader " + self.name + " = { " + self.limbsName() + ", " +\
				str(self.getNumLimbs()) + ", " + str(self.getNumDLs()) + " };\n\n" 
			data.header = "extern FlexSkeletonHeader " + self.name + ";\n"
		else:
			data.source += "SkeletonHeader " + self.name + " = { " + self.limbsName() + ", " +\
				str(self.getNumLimbs()) + " };\n\n" 
			data.header = "extern SkeletonHeader " + self.name + ";\n"

		limbData.append(data)

		return limbData

class OOTLimb():
	def __init__(self, skeletonName, index, translation, DL, lodDL):
		self.skeletonName = skeletonName
		self.translation = translation
		self.firstChildIndex = 0xFF
		self.nextSiblingIndex = 0xFF
		self.DL = DL
		self.lodDL = lodDL

		self.isFlex = False
		self.index = index
		self.children = []
		self.inverseRotation = None

	def toC(self, isLOD):
		if not isLOD:
			data = "StandardLimb "
		else:
			data = "LodLimb "

		data += self.name() + " = { " +\
			"{ " + str(int(round(self.translation[0]))) + ", " + \
			str(int(round(self.translation[1]))) + ", " + \
			str(int(round(self.translation[2]))) + " }, " + \
			str(self.firstChildIndex) + ", " +\
			str(self.nextSiblingIndex) + ", "

		if not isLOD:
			data += (self.DL.name if self.DL is not None else "NULL")
		else:
			data += "{ " + (self.DL.name if self.DL is not None else "NULL") + ", " +\
				(self.lodDL.name if self.lodDL is not None else "NULL") + " }"

		data += " };\n" 

		return data

	def name(self):
		return self.skeletonName + "Limb_" + format(self.index, '03')

	def getNumLimbs(self):
		numLimbs = 1
		for child in self.children:
			numLimbs += child.getNumLimbs()
		return numLimbs

	def getNumDLs(self):
		numDLs = 0
		if self.DL is not None or self.lodDL is not None:
			numDLs += 1

		for child in self.children:
			numDLs += child.getNumDLs()

		return numDLs

	def isFlexSkeleton(self):
		if self.isFlex:
			return True
		else:
			for child in self.children:
				if child.isFlexSkeleton():
					return True
			return False

	# should be same order as ootProcessBone
	def getList(self, limbList):
		limbList.append(self)
		for child in self.children:
			child.getList(limbList)

	def setLinks(self):
		if len(self.children) > 0:
			self.firstChildIndex = self.children[0].index
		for i in range(len(self.children)):
			if i < len(self.children) - 1:
				self.children[i].nextSiblingIndex = self.children[i + 1].index
			self.children[i].setLinks()
		# self -> child -> sibling

def setArmatureToNonRotatedPose(armatureObj):
	restPoseRotations = {}
	poseBoneName = getStartBone(armatureObj)
	setBoneNonRotated(armatureObj, poseBoneName, restPoseRotations)
	return restPoseRotations

def setBoneNonRotated(armatureObj, boneName, restPoseRotations):
	bone = armatureObj.data.bones[boneName]
	poseBone = armatureObj.pose.bones[boneName]

	while len(poseBone.constraints) > 0:
		poseBone.constraints.remove(poseBone.constraints[0])

	rotation = bone.matrix_local.inverted().decompose()[1]
	armatureObj.pose.bones[boneName].rotation_mode = "QUATERNION"
	armatureObj.pose.bones[boneName].rotation_quaternion = rotation

	restPoseRotations[boneName] = rotation

	for child in bone.children:
		setBoneNonRotated(armatureObj, child.name, restPoseRotations)

def getGroupIndices(meshInfo, armatureObj, meshObj, rootGroupIndex):
	meshInfo.vertexGroupInfo = OOTVertexGroupInfo()
	for vertex in meshObj.data.vertices:
		meshInfo.vertexGroupInfo.vertexToGroup[vertex.index] = getGroupIndexOfVert(vertex, armatureObj, meshObj, rootGroupIndex)
		#print(vertex.index, meshInfo.vertexGroupInfo.vertexToGroup[vertex.index], vertex.co)

def getGroupIndexOfVert(vert, armatureObj, obj, rootGroupIndex):
	actualGroups = []
	nonBoneGroups = []
	for group in vert.groups:
		groupName = getGroupNameFromIndex(obj, group.group)
		if groupName is not None:
			if groupName in armatureObj.data.bones:
				actualGroups.append(group)
			else:
				nonBoneGroups.append(groupName)

	if len(actualGroups) == 0:
		return rootGroupIndex
		#highlightWeightErrors(obj, [vert], "VERT")
		#raise VertexWeightError("All vertices must be part of a vertex group that corresponds to a bone in the armature.\n" +\
		#	"Groups of the bad vert that don't correspond to a bone: " + str(nonBoneGroups) + '. If a vert is supposed to belong to this group then either a bone is missing or you have the wrong group.')
	
	vertGroup = actualGroups[0]
	for group in actualGroups:
		if group.weight > vertGroup.weight:
			vertGroup = group
	#if vertGroup not in actualGroups:
	#raise VertexWeightError("A vertex was found that was primarily weighted to a group that does not correspond to a bone in #the armature. (" + getGroupNameFromIndex(obj, vertGroup.group) + ') Either decrease the weights of this vertex group or remove it. If you think this group should correspond to a bone, make sure to check your spelling.')
	return vertGroup.group

class LimbInfo:
	# Info about a limb, for the optimizer
	def __init__(self, meshObj, meshInfo, boneName, groupIndex):
		self.boneName = boneName
		self.groupIndex = groupIndex
		self.vertIndices = getVertIndices(meshObj, meshInfo, groupIndex)
		self.matlsFaces, _ = getMatlsFaces(meshInfo, self.vertIndices, boneName, groupIndex)
		self.matlsSplitVertTris = {}
		if isinstance(self.matlsFaces, dict):
			for matlIndex, faceList in self.matlsFaces.items():
				self.matlsSplitVertTris[matlIndex] = getSplitVertTris(meshInfo, faceList, meshObj.data)
		# SVIs = Split Vertex Indices
		self.allPastSVIs = set()
		self.allOwnOnlySVIs = set()
		self.allOwnAndFutureSVIs = set()
		self.allFutureOnlySVIs = set()
		# These four: past or future limbIndex : list of SVIndices
		# (They're lists because they get sorted later)
		self.pastUsedSVIsByLimb = {}
		self.pastUnusedSVIsByLimb = {}
		self.ownAndFutureSVIsByLimb = {}
		self.futureOnlySVIsByLimb = {}

class SplitVert:
	# This represents one RCP vtx, which is different from both a Blender vertex
	# and a Blender loop
	def __init__(self, groupIndex, f3dvert):
		self.groupIndex = groupIndex
		self.limbIndex = None
		self.f3d = f3dvert

def getSplitVerts(meshInfo, mesh):
	meshInfo.splitVerts = []
	for loopIndex, f3dvert in meshInfo.f3dVert.items():
		vertIndex = mesh.loops[loopIndex].vertex_index
		groupIndex = meshInfo.vertexGroupInfo.vertexToGroup[vertIndex]
		# See if equivalent f3d already exists. If in same group, they can
		# be merged.
		for sv in meshInfo.splitVerts:
			if sv.f3d == f3dvert and sv.groupIndex == groupIndex:
				break
		else:
			meshInfo.splitVerts.append(SplitVert(groupIndex, f3dvert))
		
def getSplitVertTris(meshInfo, faces, mesh):
	ret = []
	for face in faces:
		faceSplitVertIndices = []
		for loopIndex in face.loops:
			f3dvert = meshInfo.f3dVert[loopIndex]
			vertIndex = mesh.loops[loopIndex].vertex_index
			groupIndex = meshInfo.vertexGroupInfo.vertexToGroup[vertIndex]
			for i, sv in enumerate(meshInfo.splitVerts):
				if sv.f3d == f3dvert and sv.groupIndex == groupIndex:
					faceSplitVertIndices.append(i)
					break
			else:
				raise PluginError("Internal error in getSplitVertTris")
		ret.append(faceSplitVertIndices)
	return ret

def ootOptimFinalSetup(meshInfo):
	for sv in meshInfo.splitVerts:
		sv.limbIndex = meshInfo.vertexGroupInfo.vertexGroupToLimb[sv.groupIndex]
	# First arrange into past, own, and future ignoring the "only" part, then
	# afterwards separate those in both.
	for limbIndex, limbInfo in enumerate(meshInfo.limbs):
		for matlIndex, faceList in limbInfo.matlsSplitVertTris.items():
			for face in faceList:
				# Tri must belong to max limb of its three verts
				assert limbIndex == max([meshInfo.splitVerts[svIdx].limbIndex for svIdx in face])
				for svIdx in face:
					sv = meshInfo.splitVerts[svIdx]
					if limbIndex == sv.limbIndex:
						limbInfo.allOwnOnlySVIs.add(svIdx)
					else:
						limbInfo.allPastSVIs.add(svIdx)
						if sv.limbIndex not in limbInfo.pastUsedSVIsByLimb:
							limbInfo.pastUsedSVIsByLimb[sv.limbIndex] = set()
						limbInfo.pastUsedSVIsByLimb[sv.limbIndex].add(svIdx)
						meshInfo.limbs[sv.limbIndex].allFutureOnlySVIs.add(svIdx)
	for li in meshInfo.limbs:
		li.allOwnAndFutureSVIs = li.allOwnOnlySVIs & li.allFutureOnlySVIs
		li.allOwnOnlySVIs = li.allOwnOnlySVIs - li.allOwnAndFutureSVIs
		li.allFutureOnlySVIs = li.allFutureOnlySVIs - li.allOwnAndFutureSVIs
	# Separate SV indices which were used in tris of the past limb, from those
	# which have the transform of the past limb but not used in any of its tris.
	for limbInfo in meshInfo.limbs:
		# Sort in descending order of limbs (closest to present first)
		limbInfo.pastUsedSVIsByLimb = dict(sorted(limbInfo.pastUsedSVIsByLimb.items(), reverse=True))
		for pastLimbIndex in limbInfo.pastUsedSVIsByLimb.keys():
			newUsed = []
			newUnused = []
			for svIdx in limbInfo.pastUsedSVIsByLimb[pastLimbIndex]:
				if svIdx in meshInfo.limbs[pastLimbIndex].allOwnAndFutureSVIs:
					newUsed.append(svIdx)
				elif svIdx in meshInfo.limbs[pastLimbIndex].allFutureOnlySVIs:
					newUnused.append(svIdx)
				else:
					raise RuntimeError('Internal error in ootOptimFinalSetup')
			if len(newUsed) > 0:
				limbInfo.pastUsedSVIsByLimb[pastLimbIndex] = newUsed
			else:
				del limbInfo.pastUsedSVIsByLimb[pastLimbIndex]
			# These are added in order so should be sorted
			if len(newUnused) > 0:
				limbInfo.pastUnusedSVIsByLimb[pastLimbIndex] = newUnused
	# Set up ownAndFutureSVIsByLimb and futureOnlySVIsByLimb as the symmetric
	# data structures to pastUsedSVIsByLimb and pastUnusedSVIsByLimb.
	for limbIndex, limbInfo in enumerate(meshInfo.limbs):
		for pastLimbIndex, SVIs in limbInfo.pastUsedSVIsByLimb.items():
			meshInfo.limbs[pastLimbIndex].ownAndFutureSVIsByLimb[limbIndex] = SVIs
		for pastLimbIndex, SVIs in limbInfo.pastUnusedSVIsByLimb.items():
			meshInfo.limbs[pastLimbIndex].futureOnlySVIsByLimb[limbIndex] = SVIs
	for limbInfo in meshInfo.limbs:
		# Sort in ascending order of limbs (closest to present first)
		limbInfo.ownAndFutureSVIsByLimb = dict(sorted(limbInfo.ownAndFutureSVIsByLimb.items()))
		limbInfo.futureOnlySVIsByLimb = dict(sorted(limbInfo.futureOnlySVIsByLimb.items()))
	'''
	for limbIndex, limbInfo in enumerate(meshInfo.limbs):
		print('Limb {}'.format(limbIndex))
		print('pastUsedSVIsByLimb', limbInfo.pastUsedSVIsByLimb)
		print('pastUnusedSVIsByLimb', limbInfo.pastUnusedSVIsByLimb)
		print('ownAndFutureSVIsByLimb', limbInfo.ownAndFutureSVIsByLimb)
		print('futureOnlySVIsByLimb', limbInfo.futureOnlySVIsByLimb)
	'''

def isInMaterial(meshInfo, svIdx, matlIndex, limbIndex):
	# Does this split vertex belong to any tris of the given limb and given material?
	# This is done this way rather than by having each split vertex having a
	# single material index, is because a split vertex may be shared by multiple
	# materials. For example a single Blender vertex with tris which have different
	# solid color materials, but all of them have the same UVs of 0, the same
	# lighting setting, no sharp, etc.
	if not isinstance(meshInfo.limbs[limbIndex].matlsFaces, dict):
		return False
	for tri in meshInfo.limbs[limbIndex].matlsSplitVertTris[matlIndex]:
		if svIdx in tri:
			return True
	return False

def getBestStartMaterial(meshInfo, matlChoices, limbIndex):
	if limbIndex >= len(meshInfo.limbs):
		return None
	if not isinstance(meshInfo.limbs[limbIndex].matlsFaces, dict):
		return None
	ownMatls = meshInfo.limbs[limbIndex].matlsFaces.keys()
	okayStartChoices = []
	for m in matlChoices:
		if m in ownMatls:
			okayStartChoices.append(m)
	if len(okayStartChoices) == 0:
		return None
	if len(okayStartChoices) == 1:
		return okayStartChoices[0]
	bestEndChoice = getBestStartMaterial(meshInfo, ownMatls, limbIndex+1)
	if bestEndChoice in okayStartChoices:
		# Want to end with this material, so don't start with it
		okayStartChoices.remove(bestEndChoice)
	# Sort by how many verts would be carried over
	carriedOverVerts = {}
	for o in okayStartChoices:
		c = 0
		if limbIndex-1 in meshInfo.limbs[limbIndex].pastUsedSVIsByLimb:
			sharedSVs = meshInfo.limbs[limbIndex].pastUsedSVIsByLimb[limbIndex-1]
			for svIdx in sharedSVs:
				if isInMaterial(meshInfo, svIdx, o, limbIndex-1):
					c += 1
		carriedOverVerts[o] = c
	return max(carriedOverVerts, key=carriedOverVerts.get)

def ootOptimSeqMaterials(opSeq, meshInfo):
	# Compute an attempted optimal order of material loads across all limbs.
	curMatl = None
	for limbIndex in range(len(meshInfo.limbs)):
		opSeq.append(('limbstart', limbIndex))
		if not isinstance(meshInfo.limbs[limbIndex].matlsFaces, dict):
			opSeq.append(('limbend', limbIndex))
			continue
		ownMatls = list(meshInfo.limbs[limbIndex].matlsSplitVertTris.keys())
		endMatl = getBestStartMaterial(meshInfo, ownMatls, limbIndex+1)
		while len(ownMatls) > 0:
			if curMatl not in ownMatls:
				for m in ownMatls:
					if m != endMatl or len(ownMatls) == 1:
						curMatl = m
						break
				opSeq.append(('matl', curMatl))
			ownMatls.remove(curMatl)
			opSeq.append(('trilist', meshInfo.limbs[limbIndex].matlsSplitVertTris[curMatl]))
		opSeq.append(('limbend', limbIndex))
	#pprint(opSeq)

def ootOptimCombineSVBeginEnd(opSeq):
	# Combine beginnings and endings between each pair of tri commands.
	# This is currently unused.
	ptr = 0
	beginnings = set()
	endings = set()
	while ptr < len(opSeq):
		cmd = opSeq[ptr]
		if cmd[0] == 'svbegin':
			beginnings.add(cmd[1])
			del opSeq[ptr]
		elif cmd[0] == 'svend':
			endings.add(cmd[1])
			del opSeq[ptr]
		else:
			if len(endings) > 0:
				opSeq.insert(ptr, ('svendset', endings))
				endings = set()
				ptr += 1
			if len(beginnings) > 0:
				opSeq.insert(ptr, ('svbeginset', beginnings))
				beginnings = set()
				ptr += 1
			ptr += 1
			
def ootOptimPrintSeq(opSeq, meshInfo):
	# Print opSeq as sv lifetimes and uses
	# Sort svs by end
	svOrder = []
	for cmd in opSeq:
		# if cmd[0] == 'tri':
		# 	for sv in cmd[1]:
		# 		if sv not in svOrder:
		# 			svOrder.append(sv)
		if cmd[0] == 'svend' and cmd[1] not in svOrder:
			svOrder.append(cmd[1])
		elif cmd[0] == 'svendset':
			for sv in cmd[1]:
				if sv not in svOrder:
					svOrder.append(sv)
	print("\nopSeq\n")
	for i in range(len(meshInfo.splitVerts)):
		svIdx = svOrder[i]
		print("{:4d}".format(svIdx), end="")
		loaded = False
		for cmd in opSeq:
			ch = '=' if loaded else ' '
			if cmd[0] == 'limbstart':
				print("L{}".format(cmd[1]), end="")
			elif cmd[0] == 'limbend':
				print("E{}".format(cmd[1]), end="")
			elif cmd[0] == 'matl':
				print("M{}".format(cmd[1]), end="")
			elif cmd[0] == 'matrix':
				print("X{}".format(cmd[1]), end="")
			elif cmd[0] == 'svbeginset':
				if svIdx in cmd[1]:
					print("[", end="")
					loaded = True
				else:
					print(ch, end="")
			elif cmd[0] == 'svendset':
				if svIdx in cmd[1]:
					print("]", end="")
					loaded = False
				else:
					print(ch, end="")
			elif cmd[0] == 'svbegin':
				if cmd[1] == svIdx:
					print("[", end="")
					loaded = True
				else:
					print(ch, end="")
			elif cmd[0] == 'svend':
				if cmd[1] == svIdx:
					print("]", end="")
					loaded = False
				else:
					print(ch, end="")
			elif cmd[0] == 'tri':
				idxs = [svOrder.index(v) for v in cmd[1]]
				if svIdx in cmd[1]:
					print("O", end="")
				elif i > min(idxs) and i < max(idxs):
					print("|", end="")
				else:
					print(ch, end="")
			else:
				print("?", end="")
		print("")
	print("\n")

def ootOptimTriOrder(opSeq, limbIndex, meshInfo, 
	allStartLive, allOwnAndFutureEndLive, allFutureOnlyEndLive,
	pastSVIsPriority, ownAndFutureSVIsPriority):
	# Given a starting and ending state and priority list, compute the "optimal"
	# order of drawing tris in the limb. opSeq should be composed of trilists.
	limb = meshInfo.limbs[limbIndex]
	# Get bounds of sequence
	startPtr = 0
	while opSeq[startPtr] != ('limbstart', limbIndex):
		startPtr += 1
	startPtr += 1
	if opSeq[startPtr][0] == 'matl':
		startPtr += 1
	endPtr = len(opSeq) - 1
	while opSeq[endPtr] != ('limbend', limbIndex):
		endPtr -= 1
	origStartPtr = startPtr
	origEndPtr = endPtr
	# Load future only verts
	for svIdx in allFutureOnlyEndLive:
		opSeq.insert(endPtr, ('svbegin', svIdx))
	# More setup
	nsv = len(meshInfo.splitVerts)
	allTrisLeft = [tri for triList in limb.matlsSplitVertTris.values() for tri in triList]
	pastSVsLeft = allStartLive.copy()
	futureSVsLeft = allOwnAndFutureEndLive.copy()
	startLiveState = [svIdx in allStartLive for svIdx in range(nsv)]
	endLiveState = [svIdx in allOwnAndFutureEndLive for svIdx in range(nsv)]
	pastSVsNotStartAlive = limb.allPastSVIs - allStartLive
	startMatrixState = limbIndex
	startEnabled = True
	endEnabled = True
	wereEndVerts = False
	# Loop
	while True:
		limbDone = False
		# Loop as long as there was a new trilist
		while True:
			assert opSeq[startPtr][0] == 'trilist'
			startTris = opSeq[startPtr][1]
			assert opSeq[endPtr-1][0] == 'trilist'
			endTris = opSeq[endPtr-1][1]
			# Insert all tris possible with current verts
			for tri in reversed(startTris):
				if all(startLiveState[tri[i]] for i in range(3)):
					#print("At start, adding tri ", tri)
					opSeq.insert(startPtr, ('tri', tri))
					startPtr += 1
					endPtr += 1
					startTris.remove(tri)
					allTrisLeft.remove(tri)
			for tri in reversed(endTris):
				if all(endLiveState[tri[i]] for i in range(3)):
					#print("At end, adding tri ", tri)
					opSeq.insert(endPtr, ('tri', tri))
					endTris.remove(tri)
					allTrisLeft.remove(tri)
			#print("{} tris left".format(len(allTrisLeft)))
			#print("startLiveState", [i for i in range(nsv) if startLiveState[i]])
			# End or begin lifetimes for all verts which are no longer needed
			for svIdx in range(nsv):
				if startLiveState[svIdx] == endLiveState[svIdx]:
					# If alive on neither end, not relevant; if alive on both
					# ends, let it stay alive throughout
					continue
				for tri in allTrisLeft:
					if svIdx in tri:
						break
				else:
					if startLiveState[svIdx]:
						#print("Done with sv {} at start".format(svIdx))
						startLiveState[svIdx] = False
						pastSVsLeft.discard(svIdx)
						opSeq.insert(startPtr, ('svend', svIdx))
						startPtr += 1
						endPtr += 1
					else:
						assert endLiveState[svIdx]
						#print("Done with sv {} at end".format(svIdx))
						endLiveState[svIdx] = False
						futureSVsLeft.discard(svIdx)
						if svIdx in ownAndFutureSVIsPriority:
							ownAndFutureSVIsPriority.remove(svIdx)
						opSeq.insert(endPtr, ('svbegin', svIdx))
			# Only traverse in both directions if there's past or future SVs,
			# otherwise continue traversing in the last direction.
			if startEnabled and endEnabled:
				if len(pastSVsLeft) == 0:
					#print("Disabling start traversal")
					startEnabled = False
				elif len(futureSVsLeft) == 0:
					#print("Disabling end traversal")
					endEnabled = False
			# Update pointers and possibly exit
			gotNewTriLists = False
			if len(startTris) == 0:
				#print("Done with trilist at start")
				gotNewTriLists = True
				del opSeq[startPtr]
				endPtr -= 1
				if opSeq[startPtr][0] == 'matl':
					startPtr += 1
				if startPtr >= endPtr:
					limbDone = True
					break
			if len(endTris) == 0:
				#print("Done with trilist at end")
				gotNewTriLists = True
				del opSeq[endPtr-1]
				endPtr -= 1
				if opSeq[endPtr-1][0] == 'matl':
					endPtr -= 1
				if startPtr >= endPtr:
					limbDone = True
					break
			if not gotNewTriLists:
				break
		# end loop while got new trilists
		if limbDone:
			break
		def getPriorityWeighting(idx):
			return 40 + 10000 // ((2 ** idx) * (idx + 1))
		def getStartSVWeighting(svIdx):
			if pastSVIsPriority is not None:
				return getPriorityWeighting(pastSVIsPriority.index(svIdx))
			# Weight more if used closer to the present, and also
			# weight more if actually used by the past/present limb
			wgt = 40
			for SVIs in reversed(list(limb.pastUnusedSVIsByLimb.values())):
				if svIdx in SVIs:
					return wgt
				wgt += 50
			wgt = 80
			for SVIs in reversed(list(limb.pastUsedSVIsByLimb.values())):
				if svIdx in SVIs:
					return wgt
					break
				wgt += 90
			raise RuntimeError('Internal error in getStartSVWeighting')
		def getEndSVWeighting(svIdx):
			if ownAndFutureSVIsPriority is not None:
				return getPriorityWeighting(ownAndFutureSVIsPriority.index(svIdx))
			# Weight more than 2x more if used closer to the present, and
			# weight more for SVs whose ends are sooner in the future limb
			wgt = 40
			for SVIs in reversed(list(limb.ownAndFutureSVIsByLimb.values())):
				for s in reversed(SVIs):
					if s == svIdx:
						return wgt
					wgt += 2
				wgt = wgt * 2 + 1
			raise RuntimeError('Internal error in getEndSVWeighting')
		# Weight allTrisLeft based on:
		def getTriWeighting(tri, alive, allSet, wtgFunc, isStart):
			wgt = 1 # bonus for all verts in a tri so not 0
			for svIdx in tri:
				isLoaded = alive[svIdx]
				isOnlyUse = False # svIdx is only used in this tri (of the remaining ones)
				for tri2 in allTrisLeft:
					if tri2 != tri and svIdx in tri2:
						break
				else:
					isOnlyUse = True
				if svIdx in allSet:
					wgt += (2 if isOnlyUse else 1) * wtgFunc(svIdx)
				elif isLoaded:
					wgt += 15 if isOnlyUse else 5
				elif isStart and meshInfo.splitVerts[svIdx].limbIndex != startMatrixState:
					wgt -= 100
				else:
					wgt += 1 if isOnlyUse else 0
			return wgt
		startTrisWgt = [getTriWeighting(tri, startLiveState, allStartLive, 
			getStartSVWeighting, True) for tri in startTris]
		endTrisWgt = [getTriWeighting(tri, endLiveState, allOwnAndFutureEndLive,
			getEndSVWeighting, False) for tri in endTris]
		# Weight all SVs by the total score of all remaining tris they
		# contribute to
		def getSVWeighting(svIdx, tris, triswgt, alive, otheralive, enabled):
			wgt = sum([triswgt[t] for t in range(len(tris)) if svIdx in tris[t] ])
			# If the vert is not used by any tris, or already loaded,
			# or this dir not enabled, really don't load it
			if wgt == 0 or alive[svIdx] or not enabled: wgt = -100000
			# Avoid loading verts which belong to the other end
			if otheralive[svIdx]: wgt -= 1000
			return wgt
		svsStartWgt = [getSVWeighting(svIdx, startTris, startTrisWgt, 
			startLiveState, endLiveState, startEnabled) for svIdx in range(nsv)]
		svsEndWgt = [getSVWeighting(svIdx, endTris, endTrisWgt,
			endLiveState, startLiveState, endEnabled) for svIdx in range(nsv)]
		# Find the best vert to load, with ties to end
		maxWgt = max(max(svsStartWgt), max(svsEndWgt))
		assert maxWgt > -100000
		if maxWgt in svsEndWgt:
			svIdx = svsEndWgt.index(maxWgt)
			assert not endLiveState[svIdx]
			assert svIdx not in allEndLive
			#print("Adding sv {} to end with wgt {}".format(svIdx, maxWgt))
			opSeq.insert(endPtr, ('svend', svIdx))
			endLiveState[svIdx] = True
			wereEndVerts = True
		else:
			svIdx = svsStartWgt.index(maxWgt)
			assert not startLiveState[svIdx]
			assert svIdx not in allStartLive
			#print("Adding sv {} to start with wgt {}".format(svIdx, maxWgt))
			if meshInfo.splitVerts[svIdx].limbIndex != startMatrixState:
				startMatrixState = meshInfo.splitVerts[svIdx].limbIndex
				opSeq.insert(startPtr, ('matrix', startMatrixState))
			opSeq.insert(startPtr, ('svbegin', svIdx))
			startPtr += 1
			endPtr += 1
			startLiveState[svIdx] = True
	if wereEndVerts and startMatrixState != limbIndex:
		opSeq.insert(endPtr, ('matrix', limbIndex))

def ootOptimGetSVEndOrder(opSeq, limbIndex):
	# Returns a dict svIdx : endOrderIndex, where endOrderIndex is an integer
	# representing the order in which the SV lifetimes end. 0 is the first SV
	# to end, some larger number is the last. The dict is sorted by endOrderIndex.
	ptr = 0
	while opSeq[ptr] != ('limbstart', limbIndex):
		ptr += 1
	ptr += 1
	order = 0
	ret = {}
	while opSeq[ptr] != ('limbend', limbIndex):
		if opSeq[ptr][0] == 'svend':
			svIdx = opSeq[ptr][1]
			assert svIdx not in ret
			ret[svIdx] = order
			order += 1
		elif opSeq[ptr][0] == 'svendset':
			raise RuntimeError('Internal error, should not have SVs in sets at this point')
		ptr += 1
	return ret

def ootOptimSeqTrisAndSVLifetimes(opSeq, meshInfo):
	# First pass of computing triangle draw order. The only thing this is used
	# for is computing ownAndFutureSVIsPriority for each limb, which is used to
	# determine which verts to carry over between limbs. The actual output
	# sequence is discarded.
	opSeq2 = copy.deepcopy(opSeq)
	for limbIndex, limb in reversed(list(enumerate(meshInfo.limbs))):
		print("Limb {}".format(limbIndex))
		if not isinstance(limb.matlsFaces, dict):
			continue
		# No priority order for the past.
		# Priority order for the future is the order the verts end--verts which
		# end sooner are higher priority, because they can be gotten out of DMEM
		# faster.
		limb.ownAndFutureSVIsPriority = []
		for futureLimbIndex in limb.ownAndFutureSVIsByLimb.keys():
			l = limb.ownAndFutureSVIsByLimb[futureLimbIndex]
			d = ootOptimGetSVEndOrder(opSeq2, futureLimbIndex)
			l = sorted(l, key = lambda svIdx: d[svIdx])
			print("Limb {} -> {}: ".format(limbIndex, futureLimbIndex), l)
			#limb.ownAndFutureSVIsByLimb[futureLimbIndex] = l
			limb.ownAndFutureSVIsPriority.extend(l)
		# Assume all past/future SVs are alive are both ends.
		ootOptimTriOrder(opSeq2, limbIndex, meshInfo,
			limb.allPastSVIs, limb.allOwnAndFutureSVIs, limb.allFutureOnlySVIs,
			None, limb.ownAndFutureSVIsPriority)
	ootOptimPrintSeq(opSeq2, meshInfo)

def ootOptimLimbAssignVertsToSlots(opSeq, limbIndex, meshInfo, 
	nSlots, minSlot, startList, endState):
	# Given a start and end state and a certain range of DMEM vertex slots,
	# compute the "optimal" assignment of verts to vertex slots.
	# Returns a new, updated sequence, and the desired initial state of the
	# vertex slots before starting this limb.
	# startList: A list of SVs, from highest priority / shortest lifetime to
	# lowest / longest, to have resident before the limb begins. The returned
	# initial state is similar to this, but may start at either end of the DMEM
	# range.
	assert len(startList) <= nSlots
	# endState: A list of the end DMEM vertex state this limb must leave, of
	# length nSlots, where each entry is a SV index or None for don't care.
	assert len(endState) == nSlots
	limb = meshInfo.limbs[limbIndex]
	# Now that we know the actual start and end verts, recompute the tri
	# ordering.
	opSeq2 = copy.deepcopy(opSeq)
	endSet = set([svIdx in endState if svIdx is not None])
	endPriority = [svIdx in limb.ownAndFutureSVIsPriority if svIdx in endSet]
	ootOptimTriOrder(opSeq2, limbIndex, meshInfo,
		set(startList), endSet & limb.allOwnAndFutureSVIs, endSet & limb.allFutureOnlySVIs,
		startList, endPriority)
	# TODO actual slot assignment


def ootOptimComputeSeqCost(opSeq):
	# Return the estimated cost in cycles to run this sequence on the RSP.
	# Tris are considered free, since the same number of tris have to be drawn
	# regardless of which path is taken. Also materials are ignored since it is
	# assumed the material order has already been decided. So basically what is
	# counted is vertex loads (including the DMA and the transformations), and
	# matrix loads. The matrix load done by the SkelAnime system at the 
	# beginning of each limb is considered free.
	# 
	# Costs
	# 
	# The average number of cycles from having called through returning from
	# a RSP dma_read_write call, minus the actual cycles consumed transferring data. 
	# For now, completely guessed.
	dmaOverheadCost = 100
	# The number of cycles to DMA transfer the data of one vertex. RDRAM transfers
	# 8 bytes per cycle.
	dmaVertexCost = 2
	# The number of cycles to DMA transfer the data of a (modelview) matrix.
	# RDRAM transfers 8 bytes per cycle, and a matrix is 0x40 bytes.
	dmaMatrixCost = 8
	# The number of cycles to execute RSP code (not counting dma_read_write or
	# the DMA itself) for gsSPMatrix on the codepath of G_MTX_LOAD.
	# Estimated from reading code and accounting for pipeline latencies; probably
	# quite close but not exact.
	gsSPMatrixCodeCost = 16 + 8
	# The number of cycles to update the MVP matrix. This is marked dirty after
	# each matrix load and updated on the next vertex load, presumably so that
	# if both modelview and projection matrices are to be loaded, it doesn't
	# bother recomputing between those two loads. Since we only modify modelview,
	# we can just count this as another cost to loading a matrix.
	# Estimated from reading code and accounting for pipeline latencies; probably
	# quite close but not exact.
	calculateMVPCost = 6 + TODO
	
	

def ootDuplicateArmature(originalArmatureObj):
	# Duplicate objects to apply scale / modifiers / linked data
	bpy.ops.object.select_all(action = 'DESELECT')
	
	for originalMeshObj in [obj for obj in originalArmatureObj.children if isinstance(obj.data, bpy.types.Mesh)]:
		originalMeshObj.select_set(True)
		originalMeshObj.original_name = originalMeshObj.name

	originalArmatureObj.select_set(True)
	originalArmatureObj.original_name = originalArmatureObj.name
	bpy.context.view_layer.objects.active = originalArmatureObj
	bpy.ops.object.duplicate()

	armatureObj = bpy.context.view_layer.objects.active
	meshObjs = [obj for obj in bpy.context.selected_objects if obj is not armatureObj]

	try:
		for obj in meshObjs:
			setOrigin(armatureObj, obj)

		bpy.ops.object.select_all(action = "DESELECT")
		armatureObj.select_set(True)
		bpy.context.view_layer.objects.active = armatureObj
		bpy.ops.object.transform_apply(location = False, rotation = False,
			scale = True, properties = False)

		# convert blender to n64 space, then set all bones to be non-rotated
		applyRotation([armatureObj], math.radians(90), 'X')
		restPoseRotations = setArmatureToNonRotatedPose(armatureObj)
			
		# Apply modifiers/data to mesh objs
		bpy.ops.object.select_all(action = 'DESELECT')
		for obj in meshObjs:
			obj.select_set(True)
			bpy.context.view_layer.objects.active = obj

		bpy.ops.object.make_single_user(obdata = True)
		bpy.ops.object.transform_apply(location = False, 
			rotation = True, scale = True, properties =  False)
		for selectedObj in meshObjs:
			bpy.ops.object.select_all(action = 'DESELECT')
			selectedObj.select_set(True)
			bpy.context.view_layer.objects.active = selectedObj

			for modifier in selectedObj.modifiers:
				attemptModifierApply(modifier)
		
		# Apply new armature rest pose
		bpy.ops.object.select_all(action = "DESELECT")
		bpy.context.view_layer.objects.active = armatureObj
		bpy.ops.object.mode_set(mode = "POSE")
		bpy.ops.pose.armature_apply()
		bpy.ops.object.mode_set(mode = "OBJECT")

		return armatureObj, meshObjs, restPoseRotations
	except Exception as e:
		cleanupDuplicatedObjects(meshObjs + [armatureObj])
		originalArmatureObj.select_set(True)
		bpy.context.view_layer.objects.active = originalArmatureObj
		raise Exception(str(e))

def ootConvertArmatureToSkeletonWithoutMesh(originalArmatureObj, convertTransformMatrix, name):
	skeleton, fModel, restPoseRotations = ootConvertArmatureToSkeleton(originalArmatureObj, convertTransformMatrix, 
		None, name, False, True, "Opaque")
	return skeleton, restPoseRotations

def ootConvertArmatureToSkeletonWithMesh(originalArmatureObj, convertTransformMatrix, fModel, name, convertTextureData, drawLayer):
	
	skeleton, fModel, restPoseRotations = ootConvertArmatureToSkeleton(originalArmatureObj, convertTransformMatrix, 
		fModel, name, convertTextureData, False, drawLayer)
	return skeleton, fModel

def ootConvertArmatureToSkeleton(originalArmatureObj, convertTransformMatrix, 
	fModel, name, convertTextureData, skeletonOnly, drawLayer):
	checkEmptyName(name)

	armatureObj, meshObjs, restPoseRotations = ootDuplicateArmature(originalArmatureObj)
	
	try:
		skeleton = OOTSkeleton(name)

		if len(armatureObj.children) == 0:
			raise PluginError("No mesh parented to armature.")

		#startBoneNames = sorted([bone.name for bone in armatureObj.data.bones if bone.parent is None])
		#startBoneName = startBoneNames[0]
		checkForStartBone(armatureObj)
		startBoneName = getStartBone(armatureObj)
		meshObj = meshObjs[0]

		meshInfo = getInfoDict(meshObj)
		getGroupIndices(meshInfo, armatureObj, meshObj, getGroupIndexFromname(meshObj, startBoneName))

		convertTransformMatrix = convertTransformMatrix @ \
			mathutils.Matrix.Diagonal(armatureObj.scale).to_4x4()

		#for i in range(len(startBoneNames)):
		#	startBoneName = startBoneNames[i]
		
		if bpy.context.scene.ootSkeletonExportOptimize:
			getSplitVerts(meshInfo, meshObj.data)
			meshInfo.limbs = []
			ootBoneOptimSetup(0, startBoneName, armatureObj, meshObj, meshInfo)
			ootOptimFinalSetup(meshInfo)
			opSeq = []
			ootOptimSeqMaterials(opSeq, meshInfo)
			ootOptimSeqTrisAndSVLifetimes(opSeq, meshInfo)
			raise NotImplementedError("Export not done")
		else:
			ootProcessBone(fModel, startBoneName, skeleton, 0,
				meshObj, armatureObj, convertTransformMatrix, meshInfo, convertTextureData, 
				name, skeletonOnly, drawLayer, None)

		cleanupDuplicatedObjects(meshObjs + [armatureObj])
		originalArmatureObj.select_set(True)
		bpy.context.view_layer.objects.active = originalArmatureObj

		return skeleton, fModel, restPoseRotations
	except Exception as e:
		cleanupDuplicatedObjects(meshObjs + [armatureObj])
		originalArmatureObj.select_set(True)
		bpy.context.view_layer.objects.active = originalArmatureObj
		raise Exception(str(e))

def ootBoneOptimSetup(limbIndex, boneName, armatureObj, meshObj, meshInfo):
	groupIndex = getGroupIndexFromname(meshObj, boneName)
	meshInfo.vertexGroupInfo.vertexGroupToLimb[groupIndex] = limbIndex
	limb = LimbInfo(meshObj, meshInfo, boneName, groupIndex)
	meshInfo.limbs.append(limb)
	limbIndex += 1
	bone = armatureObj.data.bones[boneName]
	childrenNames = getSortedChildren(armatureObj, bone)
	for childName in childrenNames:
		limbIndex = ootBoneOptimSetup(limbIndex, childName, armatureObj, meshObj, meshInfo)
	return limbIndex

def ootProcessBone(fModel, boneName, parentLimb, nextIndex, meshObj, armatureObj, 
	convertTransformMatrix, meshInfo, convertTextureData, namePrefix, skeletonOnly,
	drawLayer, lastMaterialName):
	bone = armatureObj.data.bones[boneName]
	if bone.parent is not None:
		transform = convertTransformMatrix @ bone.parent.matrix_local.inverted() @ bone.matrix_local
	else:
		transform = convertTransformMatrix @ bone.matrix_local

	translate, rotate, scale = transform.decompose()

	groupIndex = getGroupIndexFromname(meshObj, boneName)

	meshInfo.vertexGroupInfo.vertexGroupToLimb[groupIndex] = nextIndex

	if skeletonOnly:
		mesh = None
		hasSkinnedFaces = None
	else:
		mesh, hasSkinnedFaces, lastMaterialName = ootProcessVertexGroup(
			fModel, meshObj, boneName, 
			convertTransformMatrix, armatureObj, namePrefix,
			meshInfo, drawLayer, convertTextureData, lastMaterialName)

	if bone.ootBoneType == "Custom DL":
		if mesh is not None:
			raise PluginError(bone.name + " is set to use a custom DL but still has geometry assigned to it. Remove this geometry from this bone.")
		else:
			# Dummy data, only used so that name is set correctly
			mesh = FMesh(bone.ootCustomDLName, DLFormat.Static)

	DL = None
	if mesh is not None:
		if not bone.use_deform:
			raise PluginError(bone.name + " has vertices in its vertex group but is not set to deformable. Make sure to enable deform on this bone.")
		DL = mesh.draw
		
	if isinstance(parentLimb, OOTSkeleton):
		skeleton = parentLimb
		limb = OOTLimb(skeleton.name, nextIndex, translate, DL, None)
		skeleton.limbRoot = limb
	else:
		limb = OOTLimb(parentLimb.skeletonName, nextIndex, translate, DL, None)
		parentLimb.children.append(limb)

	limb.isFlex = hasSkinnedFaces
	nextIndex += 1

	childrenNames = getSortedChildren(armatureObj, bone)
	for childName in childrenNames:
		nextIndex, lastMaterialName = ootProcessBone(fModel, childName, limb, nextIndex, meshObj, 
			armatureObj, convertTransformMatrix, meshInfo, convertTextureData, 
			namePrefix, skeletonOnly, drawLayer, lastMaterialName)

	return nextIndex, lastMaterialName

def ootConvertArmatureToC(originalArmatureObj, convertTransformMatrix, 
	f3dType, isHWv1, skeletonName, folderName, DLFormat, savePNG, exportPath, isCustomExport, drawLayer, removeVanillaData):
	skeletonName = toAlnum(skeletonName)

	fModel = OOTModel(f3dType, isHWv1, skeletonName, DLFormat, drawLayer)
	skeleton, fModel = ootConvertArmatureToSkeletonWithMesh(originalArmatureObj, convertTransformMatrix, 
		fModel, skeletonName, not savePNG, drawLayer)

	if originalArmatureObj.ootFarLOD is not None:
		lodSkeleton, fModel = ootConvertArmatureToSkeletonWithMesh(originalArmatureObj.ootFarLOD, convertTransformMatrix, 
			fModel, skeletonName + "_lod", not savePNG, drawLayer)
	else:
		lodSkeleton = None

	if lodSkeleton is not None:
		skeleton.hasLOD = True
		limbList = skeleton.createLimbList()
		lodLimbList = lodSkeleton.createLimbList()

		if len(limbList) != len(lodLimbList):
			raise PluginError(originalArmatureObj.name + " cannot use " + originalArmatureObj.ootFarLOD.name + \
				"as LOD because they do not have the same bone structure.")

		for i in range(len(limbList)):
			limbList[i].lodDL = lodLimbList[i].DL
			limbList[i].isFlex |= lodLimbList[i].isFlex
	

	data = CData()
	data.source += '#include "ultra64.h"\n#include "global.h"\n'
	if not isCustomExport:
		data.source += '#include "' + folderName + '.h"\n\n'
	else:
		data.source += '\n'

	exportData = fModel.to_c(
		TextureExportSettings(False, savePNG, "test"), OOTGfxFormatter(ScrollMethod.Vertex))
	skeletonC = skeleton.toC()

	data.append(exportData.all())
	data.append(skeletonC)

	path = ootGetPath(exportPath, isCustomExport, 'assets/objects/', folderName, False, False)
	writeCData(data, 
		os.path.join(path, skeletonName + '.h'),
		os.path.join(path, skeletonName + '.c'))

	if not isCustomExport:
		addIncludeFiles(folderName, path, skeletonName)
		if removeVanillaData:
			ootRemoveSkeleton(path, folderName, skeletonName)

class OOTDLEntry:
	def __init__(self, dlName, limbIndex):
		self.dlName = dlName
		self.limbIndex = limbIndex

def ootGetSkeleton(skeletonData, skeletonName, continueOnError):
	# TODO: Does this handle non flex skeleton?
	matchResult = re.search("(Flex)?SkeletonHeader\s*" + re.escape(skeletonName) + \
		"\s*=\s*\{\s*\{?\s*([^,\s]*)\s*,\s*([^,\s\}]*)\s*\}?\s*(,\s*([^,\s]*))?\s*\}\s*;\s*", skeletonData)
	if matchResult is None:
		if continueOnError:
			return None
		else:
			raise PluginError("Cannot find skeleton named " + skeletonName)
	return matchResult

def ootGetLimbs(skeletonData, limbsName, continueOnError):
	matchResult = re.search("(static\s*)?void\s*\*\s*" + re.escape(limbsName) + \
		"\s*\[\s*[0-9]*\s*\]\s*=\s*\{([^\}]*)\}\s*;\s*", skeletonData, re.DOTALL)
	if matchResult is None:
		if continueOnError:
			return None
		else:
			raise PluginError("Cannot find skeleton limbs named " + limbsName)
	return matchResult

def ootGetLimb(skeletonData, limbName, continueOnError):
	matchResult = re.search("([A-Za-z0-9\_]*)Limb\s*" + re.escape(limbName), skeletonData)

	if matchResult is None:
		if continueOnError:
			return None
		else:
			raise PluginError("Cannot find skeleton limb named " + limbName)

	limbType = matchResult.group(1)
	if limbType == "Lod":
		dlRegex = "\{\s*([^,\s]*)\s*,\s*([^,\s]*)\s*\}"
	else:
		dlRegex = "([^,\s]*)"

	matchResult = re.search("[A-Za-z0-9\_]*Limb\s*" + re.escape(limbName) + \
		"\s*=\s*\{\s*\{\s*([^,\s]*)\s*,\s*([^,\s]*)\s*,\s*([^,\s]*)\s*\},\s*([^, ]*)\s*,\s*([^, ]*)\s*,\s*" +\
		dlRegex + "\s*\}\s*;\s*", skeletonData, re.DOTALL)

	if matchResult is None:
		if continueOnError:
			return None
		else:
			raise PluginError("Cannot handle skeleton limb named " + limbName + " of type " + limbType)
	return matchResult

def ootImportSkeletonC(filepaths, skeletonName, actorScale, removeDoubles, importNormals, basePath, drawLayer):
	skeletonData = getImportData(filepaths)

	matchResult = ootGetSkeleton(skeletonData, skeletonName, False)
	limbsName = matchResult.group(2)

	matchResult = ootGetLimbs(skeletonData, limbsName, False)
	limbsData = matchResult.group(2)
	limbList = [entry.strip()[1:] for entry in limbsData.split(',')]

	#print(limbList)
	isLOD, armatureObj = ootBuildSkeleton(skeletonName, skeletonData, limbList, 
		actorScale, removeDoubles, importNormals, False, basePath, drawLayer)
	if isLOD:
		isLOD, LODArmatureObj = ootBuildSkeleton(skeletonName, skeletonData, limbList, 
			actorScale, removeDoubles, importNormals, True, basePath, drawLayer)
		armatureObj.ootFarLOD = LODArmatureObj
	

def ootBuildSkeleton(skeletonName, skeletonData, limbList, actorScale, removeDoubles, 
	importNormals, useFarLOD, basePath, drawLayer):
	lodString = "_lod" if useFarLOD else ""

	# Create new skinned mesh
	mesh = bpy.data.meshes.new(skeletonName + '_mesh' + lodString)
	obj = bpy.data.objects.new(skeletonName + '_mesh' + lodString, mesh)
	bpy.context.scene.collection.objects.link(obj)

	# Create new armature
	armature = bpy.data.armatures.new(skeletonName + lodString)
	armatureObj = bpy.data.objects.new(skeletonName + lodString, armature)
	armatureObj.show_in_front = True
	armatureObj.ootDrawLayer = drawLayer
	#armature.show_names = True

	bpy.context.scene.collection.objects.link(armatureObj)
	bpy.context.view_layer.objects.active = armatureObj
	#bpy.ops.object.mode_set(mode = 'EDIT')

	f3dContext = OOTF3DContext(F3D("F3DEX2/LX2", False), limbList, basePath)
	f3dContext.mat().draw_layer.oot = armatureObj.ootDrawLayer
	transformMatrix = mathutils.Matrix.Scale(1 / actorScale, 4)
	isLOD = ootAddLimbRecursively(0, skeletonData, obj, armatureObj, transformMatrix, None, f3dContext, useFarLOD)
	for dlEntry in f3dContext.dlList:
		limbName = f3dContext.getLimbName(dlEntry.limbIndex)
		boneName = f3dContext.getBoneName(dlEntry.limbIndex)
		parseF3D(skeletonData, dlEntry.dlName, obj, f3dContext.matrixData[limbName], 
			limbName, boneName, "oot", drawLayer, f3dContext)
		if f3dContext.isBillboard:
			armatureObj.data.bones[boneName].ootDynamicTransform.billboard = True
		f3dContext.clearMaterial() # THIS IS IMPORTANT
	f3dContext.createMesh(obj, removeDoubles, importNormals)
	armatureObj.location = bpy.context.scene.cursor.location

	# Set bone rotation mode.
	bpy.ops.object.select_all(action = "DESELECT")
	armatureObj.select_set(True)
	bpy.context.view_layer.objects.active = armatureObj
	bpy.ops.object.mode_set(mode = 'POSE')
	for bone in armatureObj.pose.bones:
		bone.rotation_mode = 'XYZ'

	# Apply mesh to armature.	
	if bpy.context.mode != 'OBJECT':
		bpy.ops.object.mode_set(mode = 'OBJECT')
	bpy.ops.object.select_all(action = 'DESELECT')
	obj.select_set(True)
	armatureObj.select_set(True)
	bpy.context.view_layer.objects.active = armatureObj
	bpy.ops.object.parent_set(type = "ARMATURE")

	applyRotation([armatureObj], math.radians(-90), 'X')

	return isLOD, armatureObj

def ootAddBone(armatureObj, boneName, parentBoneName, currentTransform, loadDL):
	if bpy.context.mode != 'OBJECT':
		bpy.ops.object.mode_set(mode="OBJECT")
	bpy.ops.object.select_all(action = 'DESELECT')
	bpy.context.view_layer.objects.active = armatureObj
	bpy.ops.object.mode_set(mode="EDIT")
	bone = armatureObj.data.edit_bones.new(boneName)
	bone.use_connect = False
	if parentBoneName is not None:
		bone.parent = armatureObj.data.edit_bones[parentBoneName]
	bone.head = currentTransform @ mathutils.Vector((0,0,0))
	bone.tail = bone.head + (currentTransform.to_quaternion() @ \
		mathutils.Vector((0,0.3,0)))
	
	# Connect bone to parent if it is possible without changing parent direction.
	
	if parentBoneName is not None:
		nodeOffsetVector = mathutils.Vector(bone.head - bone.parent.head)
		# set fallback to nonzero to avoid creating zero length bones
		if(nodeOffsetVector.angle(bone.parent.tail - bone.parent.head, 1) \
			< 0.0001 and loadDL):
			for child in bone.parent.children:
				if child != bone:
					child.use_connect = False
			bone.parent.tail = bone.head
			bone.use_connect = True
		elif bone.head == bone.parent.head and bone.tail == bone.parent.tail:
			bone.tail += currentTransform.to_quaternion() @ mathutils.Vector((0,0.2,0))

	if bpy.context.mode != 'OBJECT':
		bpy.ops.object.mode_set(mode="OBJECT")

def ootAddLimbRecursively(limbIndex, skeletonData, obj, armatureObj, parentTransform, parentBoneName, f3dContext, useFarLOD):

	limbName = f3dContext.getLimbName(limbIndex)
	boneName = f3dContext.getBoneName(limbIndex)
	matchResult = ootGetLimb(skeletonData, limbName, False)

	isLOD = matchResult.lastindex > 6

	if isLOD and useFarLOD:
		dlName = matchResult.group(7)
	else:
		dlName = matchResult.group(6)

	# Animations override the root translation, so we just ignore importing them as well.
	if limbIndex == 0:
		translation = [0,0,0]
	else:
		translation = [hexOrDecInt(matchResult.group(1)), hexOrDecInt(matchResult.group(2)), hexOrDecInt(matchResult.group(3))]
	nextChildIndex = hexOrDecInt(matchResult.group(4))
	nextSiblingIndex = hexOrDecInt(matchResult.group(5))

	#str(limbIndex) + " " + str(translation) + " " + str(nextChildIndex) + " " + \
	#	str(nextSiblingIndex) + " " + str(dlName))

	currentTransform = parentTransform @ mathutils.Matrix.Translation(mathutils.Vector(translation))
	f3dContext.matrixData[limbName] = currentTransform
	loadDL = dlName != "NULL"

	ootAddBone(armatureObj, boneName, parentBoneName, currentTransform, loadDL)

	# DLs can access bone transforms not yet processed.
	# Therefore were delay F3D parsing until after skeleton is processed.
	if loadDL:
		f3dContext.dlList.append(OOTDLEntry(dlName, limbIndex))
		#parseF3D(skeletonData, dlName, obj, transformMatrix, boneName, f3dContext)

	if nextChildIndex != 255:
		isLOD |= ootAddLimbRecursively(nextChildIndex, skeletonData, obj, armatureObj, currentTransform, boneName, f3dContext, useFarLOD)
	
	if nextSiblingIndex != 255:
		isLOD |= ootAddLimbRecursively(nextSiblingIndex, skeletonData, obj, armatureObj, parentTransform, parentBoneName, f3dContext, useFarLOD)
	
	return isLOD

def ootRemoveSkeleton(filepath, objectName, skeletonName):
	headerPath = os.path.join(filepath, objectName + '.h')
	sourcePath = os.path.join(filepath, objectName + '.c')

	skeletonDataC = readFile(sourcePath)
	originalDataC = skeletonDataC

	skeletonDataH = readFile(headerPath)
	originalDataH = skeletonDataH

	matchResult = ootGetSkeleton(skeletonDataC, skeletonName, True)
	if matchResult is None:
		return
	skeletonDataC = skeletonDataC[:matchResult.start(0)] + skeletonDataC[matchResult.end(0):]
	limbsName = matchResult.group(2)

	headerMatch = getDeclaration(skeletonDataH, skeletonName)
	if headerMatch is not None:
		skeletonDataH = skeletonDataH[:headerMatch.start(0)] + skeletonDataH[headerMatch.end(0):] 

	matchResult = ootGetLimbs(skeletonDataC, limbsName, True)
	if matchResult is None:
		return
	skeletonDataC = skeletonDataC[:matchResult.start(0)] + skeletonDataC[matchResult.end(0):]
	limbsData = matchResult.group(2)
	limbList = [entry.strip()[1:] for entry in limbsData.split(',')]

	headerMatch = getDeclaration(skeletonDataH, limbsName)
	if headerMatch is not None:
		skeletonDataH = skeletonDataH[:headerMatch.start(0)] + skeletonDataH[headerMatch.end(0):] 

	for limb in limbList:
		matchResult = ootGetLimb(skeletonDataC, limb, True)
		if matchResult is not None:
			skeletonDataC = skeletonDataC[:matchResult.start(0)] + skeletonDataC[matchResult.end(0):]
		headerMatch = getDeclaration(skeletonDataH, limb)
		if headerMatch is not None:
			skeletonDataH = skeletonDataH[:headerMatch.start(0)] + skeletonDataH[headerMatch.end(0):] 

	if skeletonDataC != originalDataC:
		writeFile(sourcePath, skeletonDataC)

	if skeletonDataH != originalDataH:
		writeFile(headerPath, skeletonDataH)

class OOT_ImportSkeleton(bpy.types.Operator):
	# set bl_ properties
	bl_idname = 'object.oot_import_skeleton'
	bl_label = "Import Skeleton"
	bl_options = {'REGISTER', 'UNDO', 'PRESET'}

	# Called on demand (i.e. button press, menu item)
	# Can also be called from operator search menu (Spacebar)
	def execute(self, context):
		armatureObj = None
		if context.mode != 'OBJECT':
			bpy.ops.object.mode_set(mode = "OBJECT")

		try:
			importPath = bpy.path.abspath(context.scene.ootSkeletonImportCustomPath)
			isCustomImport = context.scene.ootSkeletonImportUseCustomPath
			folderName = context.scene.ootSkeletonImportFolderName
			skeletonName = context.scene.ootSkeletonImportName
			scale = context.scene.ootActorBlenderScale
			removeDoubles = context.scene.ootActorRemoveDoubles
			importNormals = context.scene.ootActorImportNormals
			decompPath = bpy.path.abspath(bpy.context.scene.ootDecompPath)
			drawLayer = bpy.context.scene.ootActorImportDrawLayer

			filepaths = [ootGetObjectPath(isCustomImport, importPath, folderName)]
			if not isCustomImport:
				filepaths.append(os.path.join(bpy.context.scene.ootDecompPath, "assets/objects/gameplay_keep/gameplay_keep.c"))

			ootImportSkeletonC(filepaths, skeletonName, scale, removeDoubles, importNormals, decompPath, drawLayer)

			self.report({'INFO'}, 'Success!')		
			return {'FINISHED'}

		except Exception as e:
			if context.mode != 'OBJECT':
				bpy.ops.object.mode_set(mode = 'OBJECT')
			raisePluginError(self, e)
			return {'CANCELLED'} # must return a set

class OOT_ExportSkeleton(bpy.types.Operator):
	# set bl_ properties
	bl_idname = 'object.oot_export_skeleton'
	bl_label = "Export Skeleton"
	bl_options = {'REGISTER', 'UNDO', 'PRESET'}

	# Called on demand (i.e. button press, menu item)
	# Can also be called from operator search menu (Spacebar)
	def execute(self, context):
		armatureObj = None
		if context.mode != 'OBJECT':
			bpy.ops.object.mode_set(mode = "OBJECT")
		if len(context.selected_objects) == 0:
			raise PluginError("Armature not selected.")
		armatureObj = context.active_object
		if type(armatureObj.data) is not bpy.types.Armature:
			raise PluginError("Armature not selected.")

		if len(armatureObj.children) == 0 or \
			not isinstance(armatureObj.children[0].data, bpy.types.Mesh):
			raise PluginError("Armature does not have any mesh children, or " +\
				'has a non-mesh child.')

		obj = armatureObj.children[0]
		finalTransform = mathutils.Matrix.Scale(context.scene.ootActorBlenderScale, 4)

		try:
			#exportPath, levelName = getPathAndLevel(context.scene.geoCustomExport, 
			#	context.scene.geoExportPath, context.scene.geoLevelName, 
			#	context.scene.geoLevelOption)

			saveTextures = bpy.context.scene.saveTextures or bpy.context.scene.ignoreTextureRestrictions
			isHWv1 = context.scene.isHWv1
			f3dType = context.scene.f3d_type
			skeletonName = context.scene.ootSkeletonExportName
			folderName = context.scene.ootSkeletonExportFolderName
			exportPath = bpy.path.abspath(context.scene.ootSkeletonExportCustomPath)
			isCustomExport = context.scene.ootSkeletonExportUseCustomPath
			drawLayer = armatureObj.ootDrawLayer
			removeVanillaData = context.scene.ootSkeletonRemoveVanillaData

			ootConvertArmatureToC(armatureObj, finalTransform, 
				f3dType, isHWv1, skeletonName, folderName, DLFormat.Static, saveTextures,
				exportPath, isCustomExport, drawLayer, removeVanillaData)

			self.report({'INFO'}, 'Success!')		
			return {'FINISHED'}

		except Exception as e:
			if context.mode != 'OBJECT':
				bpy.ops.object.mode_set(mode = 'OBJECT')
			raisePluginError(self, e)
			return {'CANCELLED'} # must return a set

class OOT_ExportSkeletonPanel(OOT_Panel):
	bl_idname = "OOT_PT_export_skeleton"
	bl_label = "OOT Skeleton Exporter"

	# called every frame
	def draw(self, context):
		col = self.layout.column()
		col.operator(OOT_ExportSkeleton.bl_idname)
		
		prop_split(col, context.scene, 'ootSkeletonExportName', "Skeleton")
		if context.scene.ootSkeletonExportUseCustomPath:
			prop_split(col, context.scene, 'ootSkeletonExportCustomPath', "Folder")
		else:		
			prop_split(col, context.scene, 'ootSkeletonExportFolderName', "Object")
		col.prop(context.scene, "ootSkeletonExportUseCustomPath")
		col.prop(context.scene, "ootSkeletonExportOptimize")
		if context.scene.ootSkeletonExportOptimize:
			b = col.box().column()
			b.label(icon = 'LIBRARY_DATA_BROKEN', text = "Do not draw anything in SkelAnime")
			b.label(text = "callbacks or cull limbs, will be corrupted.")

		col.operator(OOT_ImportSkeleton.bl_idname)

		prop_split(col, context.scene, 'ootSkeletonImportName', "Skeleton")
		if context.scene.ootSkeletonImportUseCustomPath:
			prop_split(col, context.scene, 'ootSkeletonImportCustomPath', "File")
		else:		
			prop_split(col, context.scene, 'ootSkeletonImportFolderName', "Object")
		prop_split(col, context.scene, "ootActorImportDrawLayer", "Import Draw Layer")

		col.prop(context.scene, "ootSkeletonImportUseCustomPath")
		col.prop(context.scene, "ootActorRemoveDoubles")
		col.prop(context.scene, "ootActorImportNormals")
		col.prop(context.scene, "ootSkeletonRemoveVanillaData")
		

class OOT_SkeletonPanel(bpy.types.Panel):
	bl_idname = "OOT_PT_skeleton"
	bl_label = "OOT Skeleton Properties"
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "object"
	bl_options = {'HIDE_HEADER'} 

	@classmethod
	def poll(cls, context):
		return context.scene.gameEditorMode == "OOT" and \
			hasattr(context, "object") and context.object is not None and isinstance(context.object.data, bpy.types.Armature)

	# called every frame
	def draw(self, context):
		col = self.layout.box().column()
		col.box().label(text = "OOT Skeleton Inspector")
		prop_split(col, context.object, "ootDrawLayer", "Draw Layer")
		prop_split(col, context.object, "ootFarLOD", "LOD Skeleton")
		if context.object.ootFarLOD is not None:
			col.label(text = "Make sure LOD has same bone structure.", icon = 'BONE_DATA')

class OOT_BonePanel(bpy.types.Panel):
	bl_idname = "OOT_PT_bone"
	bl_label = "OOT Bone Properties"
	bl_space_type = 'PROPERTIES'
	bl_region_type = 'WINDOW'
	bl_context = "bone"
	bl_options = {'HIDE_HEADER'} 

	@classmethod
	def poll(cls, context):
		return context.scene.gameEditorMode == "OOT" and context.bone is not None

	# called every frame
	def draw(self, context):
		col = self.layout.box().column()
		col.box().label(text = "OOT Bone Inspector")
		prop_split(col, context.bone, "ootBoneType", "Bone Type")
		if context.bone.ootBoneType == "Custom DL":
			prop_split(col, context.bone, "ootCustomDLName", "DL Name")
		if context.bone.ootBoneType == "Custom DL" or\
			context.bone.ootBoneType == "Ignore":
			col.label(text = "Make sure no geometry is skinned to this bone.", icon = 'BONE_DATA')
		
		if context.bone.ootBoneType != "Ignore":
			col.prop(context.bone.ootDynamicTransform, 'billboard')

def pollArmature(self, obj):
	return isinstance(obj.data, bpy.types.Armature)

oot_skeleton_classes = (
	OOT_ExportSkeleton,
	OOT_ImportSkeleton,
)

oot_skeleton_panels = (
	OOT_ExportSkeletonPanel,
	OOT_SkeletonPanel,
	OOT_BonePanel,
)

def oot_skeleton_panel_register():
	for cls in oot_skeleton_panels:
		register_class(cls)

def oot_skeleton_panel_unregister():
	for cls in oot_skeleton_panels:
		unregister_class(cls)

def oot_skeleton_register():
	bpy.types.Scene.ootSkeletonExportName = bpy.props.StringProperty(
		name = "Skeleton Name", default = "gGerudoRedSkel")
	bpy.types.Scene.ootSkeletonExportFolderName = bpy.props.StringProperty(
		name = "Skeleton Folder", default = "object_geldb")
	bpy.types.Scene.ootSkeletonExportCustomPath = bpy.props.StringProperty(
		name ='Custom Skeleton Path', subtype = 'FILE_PATH')
	bpy.types.Scene.ootSkeletonExportUseCustomPath = bpy.props.BoolProperty(
		name = "Use Custom Path")
	bpy.types.Scene.ootSkeletonExportOptimize = bpy.props.BoolProperty(
		name = "Optimize",
		description = "Applies various optimizations between the limbs in a skeleton. "
			+ "If enabled, the skeleton limbs must be drawn in their normal order, "
			+ "with nothing in between and no culling, otherwise the mesh will be corrupted.")

	bpy.types.Scene.ootSkeletonImportName = bpy.props.StringProperty(
		name = "Skeleton Name", default = "gGerudoRedSkel")
	bpy.types.Scene.ootSkeletonImportFolderName = bpy.props.StringProperty(
		name = "Skeleton Folder", default = "object_geldb")
	bpy.types.Scene.ootSkeletonImportCustomPath = bpy.props.StringProperty(
		name ='Custom Skeleton Path', subtype = 'FILE_PATH')
	bpy.types.Scene.ootSkeletonImportUseCustomPath = bpy.props.BoolProperty(
		name = "Use Custom Path")

	bpy.types.Scene.ootActorRemoveDoubles = bpy.props.BoolProperty(name = "Remove Doubles On Import", default = True)
	bpy.types.Scene.ootActorImportNormals = bpy.props.BoolProperty(name = "Import Normals", default = True)
	bpy.types.Scene.ootSkeletonRemoveVanillaData = bpy.props.BoolProperty(name = "Replace Vanilla Headers On Export", default = True)
	bpy.types.Scene.ootActorImportDrawLayer = bpy.props.EnumProperty(name = "Import Draw Layer", items = ootEnumDrawLayers)

	bpy.types.Object.ootFarLOD = bpy.props.PointerProperty(type = bpy.types.Object, poll = pollArmature)

	bpy.types.Bone.ootBoneType = bpy.props.EnumProperty(name = 'Bone Type', items = ootEnumBoneType)
	bpy.types.Bone.ootDynamicTransform = bpy.props.PointerProperty(type = OOTDynamicTransformProperty)
	bpy.types.Bone.ootCustomDLName = bpy.props.StringProperty(name = 'Custom DL', default = "gEmptyDL")

	for cls in oot_skeleton_classes:
		register_class(cls)

def oot_skeleton_unregister():
	del bpy.types.Scene.ootSkeletonExportName
	del bpy.types.Scene.ootSkeletonExportFolderName
	del bpy.types.Scene.ootSkeletonExportCustomPath
	del bpy.types.Scene.ootSkeletonExportUseCustomPath
	del bpy.types.Scene.ootSkeletonExportOptimize

	del bpy.types.Scene.ootSkeletonImportName
	del bpy.types.Scene.ootSkeletonImportFolderName
	del bpy.types.Scene.ootSkeletonImportCustomPath
	del bpy.types.Scene.ootSkeletonImportUseCustomPath

	del bpy.types.Scene.ootActorRemoveDoubles
	del bpy.types.Scene.ootActorImportNormals
	del bpy.types.Scene.ootSkeletonRemoveVanillaData
	del bpy.types.Scene.ootActorImportDrawLayer

	del bpy.types.Object.ootFarLOD
	
	del bpy.types.Bone.ootBoneType
	del bpy.types.Bone.ootDynamicTransform
	del bpy.types.Bone.ootCustomDLName

	for cls in reversed(oot_skeleton_classes):
		unregister_class(cls)
