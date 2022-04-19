#####i######
##Version 7. After merging node, cut pathway were fixed
##Version 7 ISSUE
###issue :  When after merging node, similar mutation (only changed nucleotide) were changed before and after node
###			This issue maybe occured from node overwriting
##Version 8

import sys
import copy

input_file = sys.argv[1]
path_file = sys.argv[2]
gr_mut_file = sys.argv[3]
gr_se_file = sys.argv[4]
output_file = sys.argv[5]


INPUT= open(input_file, 'r')
PATH = open(path_file, 'r')
OUTPUT = open(output_file, 'w')

#################################
class Node :
	def __init__(self, key, data=None) :
		self.key = key
		self.data = data
		self.children = {}

class Path :
	def __init__(self) :
		self.head = Node(None)
		self.allpath = []
	def insert(self, readline) :
		pathline = readline.rstrip('\n').split('\t')
		curr_node = self.head
		for ch in pathline :
			if ch not in curr_node.children :
				curr_node.children[ch] = Node(ch)
			curr_node = curr_node.children[ch]
		curr_node.data = pathline
		self.allpath.append(readline)
	def search(self, pathline) :
		##pathline = readline.rstrip('\n').split('\t')
		curr_node = self.head
		for ch in pathline :
			if ch in curr_node.children :
				curr_node = curr_node.children[ch]
			else :
				return False
		if curr_node.data != None:
			return True
		##curr_node.data = pathline

	def merge(self, readline, mergepath) :
		pathline = readline.rstrip('\n').split('\t')
		curr_node = self.head
		if readline in mergepath :
			print(readline,'is merged')
			return 'No'
		ass_line = [pathline]
		merge_idx = {} ## 1 : [node1, node2, ...]
		merge_node = {}
		for idx, ch in enumerate(pathline) :
			if idx == 0 :
				curr_node = curr_node.children[ch]
				continue
			elif idx == len(pathline) -1  : 
				continue
			merge_node_list = []
			if len(curr_node.children.keys()) > 1 :
				for ch_child in curr_node.children.keys() :
					search_path = copy.deepcopy(pathline)
					search_path[idx] = ch_child
					if self.search(search_path) :
						if ch_child not in merge_node_list :
							merge_node_list.append(ch_child)
				if len(merge_node_list) > 1 :
					merge_node_list.sort()
					m_node = 'merge'
					for m in merge_node_list :
						m_node = m_node + '|' + m.split(':')[0]
					merge_idx[idx] = merge_node_list
					merge_node[idx] = m_node

		for i in merge_idx.keys() :
			ass_tmp_line = copy.deepcopy(ass_line)
			ass_line = []
			for atl in ass_tmp_line :
				addline = copy.deepcopy(atl)
				for mn in merge_idx[i] :
					tmp_add = copy.deepcopy(addline)
					tmp_add[i] = mn
					if tmp_add not in ass_line :
						ass_line.append(tmp_add)
		for li in ass_line :
			flagpath = '\t'.join(li) 
			if flagpath not in mergepath :
				mergepath.append(flagpath)

		for j in merge_node.keys() :
			pathline[j] = merge_node[j]

		return pathline




print('Inserting path ...')
lines = PATH.readlines()
path_lines = lines
path_dict = Path()
for idx, line in enumerate(lines) :
	path_dict.insert(line)





## lines = INPUT.readlines()
print('Making mutation information')
lines = INPUT.readlines()
input_lines = lines
var_dict = {} # var : date
##vec_dict = {}
geo_dict = {} # var : list(=geo_list)
for idx, line in enumerate(lines) :
	if line[0] == '#' :
		OUTPUT.write(line)
		continue
	header = line.rstrip('\n').split('\t')
	before = header[0]
	after = header[1]
	b_date = int(header[2])
	a_date = int(header[3])
	b_geo = header[4].split(',')
	b_geo.sort()
	a_geo = header[5].split(',')
	a_geo.sort()

	if before not in var_dict.keys() :
		var_dict[before] = b_date
	if after not in var_dict.keys() :
		var_dict[after] = a_date
	if before not in geo_dict.keys() :
		geo_dict[before] = b_geo
	if after not in geo_dict.keys() :
		geo_dict[after] = a_geo



##merge stage
print('merged line')
merged_line = []
w_mline = []
tcnt = len(path_lines)
print('Total Path :', tcnt)
for idx, line in enumerate(path_lines) :
	if idx % 100 == 0 :
		print(round(float(idx) / float(tcnt) * 100, 2), '%' )
	m_line = path_dict.merge(line.rstrip('\n'), merged_line)
	if m_line == 'No' : continue
	if m_line not in w_mline :
		w_mline.append(m_line)
print('Writing..')
w_log = []
w_line = []
for ll in w_mline :
	header = ll   ###.split('\t')
	length = len(header)
	for i in range(0, length -1) :
		before = header[i].split(':')[0]
		after = header[i+1].split(':')[0]
		if before + '\t' + after in w_log :
			continue
		else :
			w_log.append(before + '\t' + after)

		if before[0:5] == 'merge' :
			node_list = before.split('|')[1:]
			m_node_geo = []
			da = 0
			for node in node_list :
				m_geo = geo_dict[node]
				for mg in m_geo :
					if mg not in m_node_geo :
						m_node_geo.append(mg)
				m_date = var_dict[node]
				if da == 0 :
					da = m_date
				if da > m_date :
					da = m_date
			before_date = str(da)
			m_node_geo.sort()
			before_geo = ','.join(m_node_geo)
		else :
			before_date = str(var_dict[before])
			before_geo = ','.join(geo_dict[before])

		if after[0:5] == 'merge' :
			node_list = after.split('|')[1:]
			m_node_geo = []
			da = 0
			for node in node_list :
				m_geo = geo_dict[node]
				for mg in m_geo :
					if mg not in m_node_geo :
						m_node_geo.append(mg)
				m_date = var_dict[node]
				if da == 0 :
					da = m_date
				if da > m_date :
					da = m_date
			after_date = str(da)
			m_node_geo.sort()
			after_geo = ','.join(m_node_geo)
		else :
			after_date = str(var_dict[after])
			after_geo = ','.join(geo_dict[after])
		OUTPUT.write(before + '\t' + after + '\t' + before_date + '\t' + after_date + '\t' + before_geo + '\t' + after_geo + '\n')
