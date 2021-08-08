#####i######
##Version 7. After merging node, cut pathway were fixed
##Version 7 ISSUE
###issue :  When after merging node, similar mutation (only changed nucleotide) were changed before and after node 
###			This issue maybe occured from node overwriting
##Version 8

import sys

input_file = sys.argv[1]
path_file = sys.argv[2]
gr_mut_file = sys.argv[3]
gr_se_file = sys.argv[4]
output_file = sys.argv[5]


INPUT= open(input_file, 'r')
PATH = open(path_file, 'r')
GR_MUT = open(gr_mut_file, 'r')
GR_SE = open(gr_se_file, 'r')


## lines = INPUT.readlines()

lines = INPUT.readlines()
input_lines = lines
var_dict = {} # var : date
##vec_dict = {}
geo_dict = {} # var : list(=geo_list)
for idx, line in enumerate(lines) :
	if line[0] == '#' : continue
#		output.write(line)
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

##	vec_dict[before] = after

lines = GR_MUT.readlines()

#grouping list
group_var_dict = {} # group1 : [var1, var2 ...]
for idx, line in enumerate(lines) :
	if idx == 0 : continue
	header = line.rstrip('\n').split('\t')
	g_key = header[0]
	v_list = header[1].split(':')
	group_var_dict[g_key] = v_list

single_vec = {}


lines = PATH.readlines()
before_dict = {} # group_1 : { variant : [var1, var2, var3 ...], variant2 : [var1, var2, ...]}
after_dict = {} # group_1 : { variant : [var1, var2, var3 ...], variant2 : [var1, var2, ...]}
start_node = []
end_node = []
for g in group_var_dict.keys() :
	if g not in before_dict.keys() :
		before_dict[g] = {}
		after_dict[g] = {}
	for idx, line in enumerate(lines) :
		if idx == 0 : continue
		header = line.rstrip('\t\n').split('\t')
		start = header[0].split(':')[0]
		if start not in group_var_dict[g] : continue
		path_length = len(header)
		if path_length <= 2 :
			start2 = header[0].split(':')[0]
			end2 = header[1].split(':')[0]
			if start2 not in after_dict[g].keys() :
				after_dict[g][start2] = [end2]
			elif start2 in after_dict[g].keys() :
				if end2 in after_dict[g][start2] :
					pass
				elif end2 not in after_dict[g][start2] :
					after_dict[g][start2].append(end2)
				else :
					print('Error!!, check your after dictionary!!')
			else :
				print('Error!!, check your path length!!')
			if end2 not in before_dict[g].keys() :
				before_dict[g][end2] = [start2]
			elif end2 in before_dict[g].keys() :
				if start2 in before_dict[g][end2] :
					pass
				elif start2 not in before_dict[g][end2] :
					before_dict[g][end2].append(start2)
				else :
					print('Error!!, check your after dictionary!!')
			else :
				print('Error!!, check your path length!!')
			### add to single_vec
			if start2 not in single_vec.keys() : 
				single_vec[start2] = [end2]
			else : 
				if end2 not in single_vec[start2] : 
					single_vec[start2].append(end2)
			if start2 not in start_node : 
				start_node.append(start2)
			if end2 not in start_node : 
				start_node.append(end2)
			continue
		for i in range(1, path_length - 1) :
			mid = header[i].split(':')[0]
			bef = header[i-1].split(':')[0]
			aft = header[i+1].split(':')[0]
#			print i,aft
			if mid in group_var_dict[g] : 
				if mid not in before_dict[g].keys() : 
					before_dict[g][mid]=[bef]
				else : 
					if bef not in before_dict[g][mid] : 
						before_dict[g][mid].append(bef)	

				if aft not in before_dict[g].keys() : 
					before_dict[g][aft]=[mid]
				else : 
					if mid not in before_dict[g][aft] : 
						before_dict[g][aft].append(mid)

				if mid not in after_dict[g].keys() : 
					after_dict[g][mid]=[aft]
				else : 
					if aft not in after_dict[g][mid] : 
						after_dict[g][mid].append(aft)	

				if bef not in after_dict[g].keys() : 
					after_dict[g][bef]=[mid]
				else : 
					if mid not in after_dict[g][bef] : 
						after_dict[g][bef].append(mid)
			if i == 1 : 
				if bef not in start_node : 
					start_node.append(bef)
			if i == path_length -2 : 
				if aft not in end_node : 
					end_node.append(aft)



merge_list = []
merge_true_dict = {}
edit_d = {}
##before_del = []
after_del = []
merge_dict = {} # {var1|var2 : [var1, var2]
merge_date = {}
a_to_m = {}

all_merge = []
b_to_m = {}

#######MERGE STAGE#####
for g in group_var_dict.keys() :
	print(g)
	mid_list = []
	mid_seed = {}

	for var in before_dict[g].keys() : 
		mid_list.append(var)
	for var in after_dict[g].keys() :
		mid_list.append(var)
	for var in mid_list : 
		if var not in mid_seed.keys() : 
			var_cnt = mid_list.count(var)
			if int(var_cnt) == 2 : 
				mid_seed[var] = int(var_cnt)
			elif int(var_cnt) == 1 : 
				continue
			else : 
				print( 'Error!!', var, 'count is ', str(var_cnt))

	##Merge node##
	for seed in mid_seed.keys() :

		seed_before_list = before_dict[g][seed]
		seed_after_list = after_dict[g][seed]

	##	pt_line = {} # start : [end1, end2]
		
		## Front to Back Method : before to after, 
		for sb in seed_before_list : 
			mfsb_list = after_dict[g][sb]
			candi_merge_list = []

			for sa in seed_after_list : 
				
				mfsa_list = before_dict[g][sa]
				intersection = list(set(mfsb_list) & set(mfsa_list))
				if seed not in intersection	: 
					print('Error!! Check your seed and intersection list!')
					print(seed, sb, sa)
					continue
				if len(intersection) == 1 : 
#######################################################################################
					if sb not in single_vec.keys() : 
						single_vec[sb] = [seed]
					else : 
						if seed not in single_vec[sb] :
							single_vec[sb].append(seed)
					if seed not in single_vec.keys() : 
						single_vec[seed] = [sa]
					else : 
						if sa not in single_vec[seed] :
							single_vec[seed].append(sa)
					continue
				intersection.sort()
				m_node = 'merge'
				m_node_geo = []
				da = 0
				for k in intersection : 
					if k not in all_merge : 
						all_merge.append(k)
					m_node = m_node + '|' + k
					for k_geo in geo_dict[k] : 
						m_node_geo.append(k_geo)
					if da == 0 : 
						da = var_dict[k]
					if var_dict[k] < da : 
						da = var_dict[k]
				if m_node not in merge_dict.keys() : 
					merge_dict[m_node] = [sa]
				else : 
					if sa not in merge_dict[m_node] : 
						merge_dict[m_node].append(sa)
				if sa not in a_to_m.keys() :
					a_to_m[sa] = [m_node] 
				else : 
					if m_node not in a_to_m[sa] : 
						a_to_m[sa].append(m_node)
				if sb not in b_to_m.keys() :
					b_to_m[sb] = [m_node] 
				else : 
					if m_node not in b_to_m[sb] : 
						b_to_m[sb].append(m_node)
				if m_node not in merge_date.keys() : 
					merge_date[m_node] = str(da)
				if m_node not in geo_dict.keys() : 
					geo_dict[m_node] = m_node_geo
	


##single vec check

for sin in single_vec.keys() : 
	##group check
	for grp in group_mut.keys() : 
		if sin in group_mut[grp] : 
			sin_gr = grp
			break
	for sinext in







output = open(output_file , 'w')

w_list = []
## Node overwriting
print('single vec')
for sing in single_vec.keys() : 
	print('single start :',sing)
	print(single_vec[sing])
	print('-----------------------')

print('merge path')
for mg in merge_dict.keys() : 
	print('merge node :', mg)
	print(merge_dict[mg])
	print('-----------------------')

print('before to merge')
for bbb in b_to_m.keys() : 
	print('before node :', bbb)
	print(b_to_m[bbb])
	print('-----------------------')

print('after to merge')
for aaa in a_to_m.keys() : 
	print('after node :', aaa)
	print(a_to_m[aaa])
	print('-----------------------')




for idx, line in enumerate(input_lines) :
	if line[0] == '#' :
		output.write(line)
		continue
	header = line.rstrip('\n').split('\t')
	before = header[0]
	after = header[1]
	b_date = header[2]
	a_date = header[3]
	b_geo = header[4]
	a_geo = header[5]
	ori_before = before
	ori_after = after
	if ori_before in single_vec.keys() : 
		if ori_after in single_vec[ori_before] : 
			b_geo_w = b_geo
			a_geo_w = a_geo 
			
			w_line = ori_before + '\t' + ori_after + '\t' + b_date + '\t' + a_date + '\t' + b_geo_w + '\t' + a_geo_w + '\n'
			if w_line not in w_list : 
				w_list.append(w_line)

		if ori_after in all_merge : 
			if ori_before in b_to_m.keys() : 
				for mer in b_to_m[ori_before] : 
					mer_n_list = mer.split('|')
					if ori_after in mer_n_list : 
						b_geo_w = b_geo
						a_geo_w = ','.join(geo_dict[mer])
						w_line = ori_before + '\t' + mer + '\t' + b_date + '\t' + merge_date[mer] + '\t' + b_geo_w + '\t' + a_geo_w + '\n'
						if w_line not in w_list : 
							w_list.append(w_line)


	if ori_before in all_merge : 
		if ori_after not in a_to_m.keys() : 
			print('Error!!', ori_after, 'is not in a_to_m key!!')
			continue
		if ori_after in a_to_m.keys() : 
			for mer in a_to_m[ori_after] : 
				mer_n_list = mer.split('|')
				if ori_before == 'envelope_nonsyn_49V>L_145G>T' : print(mer_n_list, '\nflag')
				if ori_before in mer_n_list : 								
					if ori_after in single_vec.keys() : 
						b_geo_w = ','.join(geo_dict[mer])
						a_geo_w = a_geo 
						w_line = mer + '\t' + ori_after + '\t' + merge_date[mer] + '\t' + a_date + '\t' + b_geo_w + '\t' + a_geo_w + '\n'
						if w_line not in w_list : 
							w_list.append(w_line)
					if ori_after in all_merge : 
						if ori_before in b_to_m.keys() : 
							for amer in b_to_m[ori_before] : 
								mer_n_list = amer.split('|')
								if ori_after in mer_n_list : 
									b_geo_w = ','.join(geo_dict[mer])
									a_geo_w = ','.join(geo_dict[amer])
									w_line = mer + '\t' + amer + '\t' + merge_date[mer] + '\t' + merge_date[amer] + '\t' + b_geo_w + '\t' + a_geo_w + '\n'
									if w_line not in w_list : 
										w_list.append(w_line)
					if ori_after in end_node : 
						b_geo_w = ','.join(geo_dict[mer])
						a_geo_w = a_geo
						w_line = mer + '\t' + ori_after + '\t' + merge_date[mer] + '\t' + a_date + '\t' + b_geo_w + '\t' + a_geo_w + '\n'
						if w_line not in w_list : 
							w_list.append(w_line)
						



for w in w_list :
	output.write(w)

output.close()


