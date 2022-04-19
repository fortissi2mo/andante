import sys
import re
import glob
#dag version 041 date 210414

input_file=sys.argv[1]
path_file = sys.argv[2]
gr_mut_file = sys.argv[3]
gr_se_file = sys.argv[4]

output_path=sys.argv[5]

cnt_interval = 10
X_list=[201912, 202001, 202002, 202003, 202004, 202005, 202006, 202007, 202008, 202009, 202010, 202011, 202012, 202101, 202102, 202103, 202104, 202105, 202106, 202107]

print '#################################################################'

INPUT= open(input_file, 'r')
PATH = open(path_file, 'r')
GR_MUT = open(gr_mut_file, 'r')
GR_SE = open(gr_se_file, 'r')

MERGEINFO = open(output_path + '_merge_node.tsv' , 'w')


class Mutation :
	def __init__(self, name, date, geo, X_list, before=None, after=None, end=None) :
		self.name = name
		self.date = date
		self.geo = geo.split(',')
		self.group = group
		self.x = X_list.index(int(date))
		self.y = None
		if before != None :
			self.before = [before]
		else :
			self.before = []
		if after != None :
			self.after = [after]
		else :
			self.after = []
		if end != None :
			self.end = {end : ''}
		else :
			self.end = {}
class MutationInfo :
	def __init__(self) :
		self.info = {}
		self.matrix = {} ## y : {x : node}
		##self.matrix_x = {}  ## x : {y : node}
		self.wlist = []
	def insert(self, line, X_list) :
		header = line.rstrip('\n').split('\t')
		before = header[0]
		before_date = header[2]
		before_geo = header[4]
		after = header[1]
		after_date = header[3]
		after_geo = header[5]
		self.wlist.append([before, after])
		if before not in self.info :
			self.info[before] = Mutation(before, before_date, before_geo, X_list, None, after)
		else :
			self.info[before].after.append(after)
		if after not in self.info :
			self.info[after] = Mutation(after, after_date, after_geo, X_list, before, None)
		else :
			self.info[after].before.append(before)
	def positioning(self, maxroot, pathinfo, maxkey_list = None, y_flag = False, endlist = None, turnback = True, start_stop = True) :
		global cnt_y
		global cnt_interval
		global complete_root
		before_x = -1
		if maxroot in complete_root :
			return
		if maxkey_list == None :
			maxkey_list = ['\t'.join(maxroot)]
		else :
			maxkey_list.append('\t'.join(maxroot))
		for idx, node in enumerate(maxroot) :
			if idx == 0 :
				start_node = node
			if self.info[node].y == None :
				if y_flag :
					cnt_y = cnt_y - cnt_interval
					y_flag = False
				self.info[node].y = cnt_y
#				print maxroot
#				print 'Positioning with', node
##				print 'add y'
#			else :
#				print maxroot
#				print maxroot in complete_root
#				print node, 'was not positioned'
			if self.info[node].y not in self.matrix.keys() :
				self.matrix[self.info[node].y] = {}
			if int(before_x) == int(self.info[node].x) :
				binning = pathinfo.max_binning[int(self.info[node].date)]
				self.info[node].x = float(before_x) + float(1/binning)
			if self.info[node].x not in self.matrix[self.info[node].y].keys() :
##				print self.info[node].y, self.info[node].x
				self.matrix[self.info[node].y][self.info[node].x] = node
			else :
				if self.matrix[self.info[node].y][self.info[node].x] != node :
					print 'Error!! In this node ', node, ',', self.info[node].x, self.info[node].y, self.matrix[self.info[node].y][self.info[node].x], 'is already exist!!'
			end_node = node
			before_x = self.info[node].x
		other_root_list = list(pathinfo.pathlist[start_node][end_node].keys())
##		print 'other root list', other_root_list
		sort_dict = {}
		for mk in maxkey_list :
			if mk in other_root_list :
				del other_root_list[other_root_list.index(mk)]
		if len(other_root_list) > 0 :
			for other_root in other_root_list :
				node_cnt = pathinfo.pathlist[start_node][end_node][other_root]['node']
				if node_cnt not in sort_dict.keys() :
					sort_dict[node_cnt] = [other_root]
				else :
					sort_dict[node_cnt].append(other_root)
			max_cnt = max(list(sort_dict.keys()))
			inputpath = sort_dict[max_cnt][0].split('\t')
			if inputpath not in complete_root :
				self.positioning(inputpath, pathinfo, maxkey_list, True, endlist, False)
				complete_root.append(inputpath)
		if endlist == None :
##			print 'end node node node'
			endlist = [end_node]
		else :
			endlist.append(end_node)
##		print 'end list', endlist

	## next node ...
		if turnback :
			for j in range(len(maxroot) - 2, -1, -1) :
				curr_node = maxroot[j]
				if curr_node == start_node :
					if start_stop :
						continue
##				print 'current node', curr_node
				end_node_list = []
				sort_dict = {}  ##length : {node_cnt : endnode}
				total_end_node_list = list(self.info[curr_node].end.keys())
##				print 'total end node list', total_end_node_list
				for tend in total_end_node_list :
					if endlist :
#					if endlist == None :
#						endlist = []
						if tend not in endlist :
							end_node_list.append(tend)
				for en in end_node_list :
					length = self.info[en].x - self.info[curr_node].x
					node_cnt = self.info[curr_node].end[en]
					if length not in sort_dict.keys() :
						sort_dict[length] = {node_cnt : [en]}
					else :
						if node_cnt not in sort_dict[length].keys() :
							sort_dict[length][node_cnt] = [en]
						else :
							sort_dict[length][node_cnt].append(en)
				lengthsort = list(sort_dict.keys())
				lengthsort.sort()
				lengthsort.reverse()
				for ls in lengthsort :
					nodesort = list(sort_dict[ls].keys())
					nodesort.sort()
					nodesort.reverse()
					for ns in nodesort :
						for enode in sort_dict[ls][ns] :
							endlist.append(enode)  ### add 210828
							maxroot2 = self.findroot(pathinfo,start_node, enode, curr_node)
							if maxroot not in complete_root :
								self.positioning(maxroot2, pathinfo, None, True, endlist, True)
								complete_root.append(maxroot2)
	def findroot(self, pathinfo, start_node, end_node, searchnode) :
		path_list = pathinfo.pathlist[start_node][end_node].keys()
		for path in path_list :
			nodelist = path.split('\t')
			if searchnode in nodelist :
				return nodelist






class Pathinfo :
	def __init__(self) :
		self.start = {} # start : [end]
		self.pathlist = {} ## start : { end : {'\t'.join(pathlist1) : 'length' : int, 'node' : int}, [pathlist2] .. ]}
		self.maxpath = {} ## start_node : {'length' : int, 'node' : int, 'path' : list}
		self.max_binning = {} ## group {date1 : int, date2 : int ...}
	def insert(self, readline, mutinfo) :
		pathline = readline.rstrip('\n').split('\t')
		plist = []
		path_date = {}
		node_cnt = len(pathline)
		start_node = pathline[0].split(':')[0]
		end_node = pathline[-1].split(':')[0]
		if start_node not in self.pathlist :
			self.pathlist[start_node] = {}
		for idx, node_date in enumerate(pathline) :
			header = node_date.split(':')
			node = header[0]
			date = int(header[1])
			if date not in path_date :
				path_date[date] = 1
			else :
				path_date[date] += 1
			plist.append(node)
			if idx == 0 :
				start_date = date
				mutinfo.info[node].end[end_node] = node_cnt - 1 - idx
				if node not in plist :
					self.pathlist[node] = {}
					self.start = []
			elif idx == node_cnt - 1 :
				end_date = date
			else :
				mutinfo.info[node].end[end_node] = node_cnt - 1 - idx
		pl_key = '\t'.join(plist)
		length = end_date - start_date
		if end_node not in self.pathlist[start_node] :
			self.pathlist[start_node][end_node] = {pl_key : {'length' : length, 'node' : node_cnt}  }
			self.start[start_node] = [end_node]
		else :
			self.pathlist[start_node][end_node][pl_key] = {'length' : length, 'node' : node_cnt}
			self.start[start_node].append(end_node)
		## max binning
		for d in path_date.keys() :
			if d not in self.max_binning :
				self.max_binning[d] = path_date[d]
			else  :
				if path_date[d] > self.max_binning[d] :
					self.max_binning[d] = path_date[d]

		## max length check
		if start_node not in self.maxpath :
			self.maxpath[start_node] = {'length' : end_date - start_date , 'node' : node_cnt, 'path' : plist}
		else :
			if self.maxpath[start_node]['length'] < end_date - start_date :
				self.maxpath[start_node] = {'length' : end_date - start_date , 'node' : node_cnt, 'path' : plist}
			elif self.maxpath[start_node]['length'] == end_date - start_date :
				if self.maxpath[start_node]['node'] < node_cnt :
					self.maxpath[start_node] = {'length' : end_date - start_date , 'node' : node_cnt, 'path' : plist}
##		if start_node not


	def maxroot(self, exc_str) :
		start_node = self.start.keys()
		max_len = 0
		max_node = 0
		max_path = ''

		for start in start_node :
			chg_flag = False
			if start in exc_str : continue
			if self.maxpath[start]['length'] > max_len :
				chg_flag = True
			elif self.maxpath[start]['length'] == max_len :
				if self.maxpath[start]['node'] > max_node :
					chg_flag = True
			if chg_flag :
				max_len = self.maxpath[start]['length']
				max_node = self.maxpath[start]['node']
				max_path = self.maxpath[start]['path']
		output_path = []
		for p in max_path :
			output_path.append(p)
		return output_path



######################## Positioning ######################




lines = GR_MUT.readlines()
groups = {} # {group1 : [mut1, mut2], group2 : [mut3, mut4, mut5] ... }
for idx, line in enumerate(lines) :
	if idx == 0 : continue
	header=line.rstrip('\n').split('\t')

	group = header[0]
	component = header[1].split(':')
	groups[group] = component



##insert mutation information
lines = INPUT.readlines()
mutinfo = {}
for idx, line in enumerate(lines) :
	if idx == 0 : continue
	gg = ''
	sample_node = line.split('\t')[0]
	for gr in groups.keys() :
		if sample_node in groups[gr] :
			gg = gr
			break
	if gg not in mutinfo  :
		mutinfo[gg] = MutationInfo()
		mutinfo[gg].insert(line, X_list)
	else :
		mutinfo[gg].insert(line, X_list)



##insert path infomation
path_info = {}  ###Pathinfo()
path_lines = PATH.readlines()
for idx, line in enumerate(path_lines) :
	if idx == 0 : continue
	ggg=''
	sample_node = line.split('\t')[0].split(':')[0]
	for gr in groups.keys() :
		if sample_node in groups[gr] :
			ggg = gr
			break
	if ggg not in path_info :
		path_info[ggg] = Pathinfo()
		path_info[ggg].insert(line, mutinfo[ggg])
	else :
		path_info[ggg].insert(line, mutinfo[ggg])

##insert coordinate position
gps = groups.keys()
gps.sort()
merge_number = 0
for g in gps :
	cnt_y = 0
	complete_root = []
	start_list = path_info[g].start
	coord_flag = True
	start_node = []
	total_start = len(start_list)
	print 'total start', total_start
	while coord_flag :
		max_root = path_info[g].maxroot(start_node)
		print '\nmax root', max_root
		if __name__=='__main__' :
			sys.setrecursionlimit(1000000)
			mutinfo[g].positioning(max_root, path_info[g], None, False, None, True, False)
			complete_root.append(max_root)
		start_node.append(max_root[0])
		cnt_y = cnt_y - cnt_interval
		if total_start == len(start_node) :
			coord_flag = False
##############################
##Adjusting position ver 2.0
	node_list = mutinfo[g].info.keys()
	mutinfo[g].matrix = {}
	coord_x = {} ## {x : [node1, node2, ...]}
	for nl in node_list : 
		nx = float(mutinfo[g].info[nl].x)
		ny = float(mutinfo[g].info[nl].y)
		if nx not in coord_x.keys() : 
			# coord_x[nx] = [nl]
			coord_x[nx] = {ny : nl}
		else : 
			if ny not in coord_x[nx] : 
				coord_x[nx][ny] = nl
			else : 
				print('Error!!', 'Node', nl, 'in' ,nx, ny, 'is exist!!', coord_x[nx][ny] )
			##coord_x[nx].append(nl)
	nxkeys = list(coord_x.keys())
	for nxk in nxkeys : 
		ncnt = len(coord_x[nxk]) + 1
		term = cnt_y / ncnt
		nykeys = coord_x[nxk].keys()
		nykeys.sort()
		nykeys.reverse()
		for idx, nyk in enumerate(nykeys) : 
			edit_y = term * (idx + 1)
			nn = coord_x[nxk][nyk]
			if edit_y not in mutinfo[g].matrix : 
				mutinfo[g].matrix[edit_y] = {nxk : nn}
			else : 
				mutinfo[g].matrix[edit_y][nxk] = nn
			mutinfo[g].info[nn].x = nxk
			mutinfo[g].info[nn].y = edit_y




###########^^^          VER2.0           ^^^^##########



###### Writing Results ######
##	merge_number = 0

	merge_dict = {}
	output = open(output_path + '_' + g + '.R.input.tsv' , 'w')
	output_table = open(output_path + '_' + g + '.R.table.tsv' , 'w')
	tribble = 'name\tlabel\tx\ty\tprotein\tcoding\tgeo'
	table_wlist = [tribble]
	y_list = list(mutinfo[g].matrix.keys())
	y_list.sort()
	y_list.reverse()
	output_wline = ''
	for y in y_list :
		print 'y',y
		x_list = list(mutinfo[g].matrix[y].keys())
		x_list.sort()
		for x in x_list :
			node = mutinfo[g].matrix[y][x]
			print x, node
			if node[0:5] == 'merge' :
				merge_number += 1
				merge_dict[node] = 'M' + str(merge_number)
				MERGEINFO.write('M' + str(merge_number) + '\t' + node + '\n')
	for pair in mutinfo[g].wlist :
		before = pair[0]
		after = pair[1]
		if before[0:5] == 'merge' :
			mergeid = merge_dict[before]
			before_form = 'Merge_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + '.'.join(mutinfo[g].info[before].geo)
			before_table = before_form + '\t' + mergeid + '\t' + str(mutinfo[g].info[before].x) + '\t' + str(mutinfo[g].info[before].y) + '\tMerge\tMerge\t' + '.'.join(mutinfo[g].info[before].geo)
		else :
			before_form = before.replace('>','_') + '_' + '.'.join(mutinfo[g].info[before].geo)
			before_table = before_form + '\t' + before.split('_')[2].replace('>','_') + '\t' + str(mutinfo[g].info[before].x) + '\t' + str(mutinfo[g].info[before].y) + '\t' + before.split('_')[0].replace('>','_') + '\t' + before.split('_')[1].replace('>','_') + '\t' + '.'.join(mutinfo[g].info[before].geo)
		if after[0:5] == 'merge' :
			mergeid = merge_dict[after]
			after_form = 'Merge_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + mergeid + '_' + '.'.join(mutinfo[g].info[after].geo)
			after_table = after_form + '\t' + mergeid + '\t' + str(mutinfo[g].info[after].x) + '\t' + str(mutinfo[g].info[after].y) + '\tMerge\tMerge\t' + '.'.join(mutinfo[g].info[after].geo)
		else :
			after_form = after.replace('>','_') + '_' + '.'.join(mutinfo[g].info[after].geo)
			after_table = after_form + '\t' + after.split('_')[2].replace('>','_') + '\t' + str(mutinfo[g].info[after].x) + '\t' + str(mutinfo[g].info[after].y) + '\t' + after.split('_')[0].replace('>','_') + '\t' + after.split('_')[1].replace('>','_') + '\t' + '.'.join(mutinfo[g].info[after].geo)
		if before_table not in table_wlist :
			table_wlist.append(before_table)
		if after_table not in table_wlist :
			table_wlist.append(after_table)

		output_wline = output_wline + ',\n' + after_form + ' ~ ' + before_form

	output.write(output_wline.lstrip(',\n'))
	for tw in table_wlist :
		output_table.write(tw + '\n')


	output.close()
	output_table.close()


INPUT.close()
PATH.close()
GR_MUT.close()
GR_SE.close()
MERGEINFO.close()
