import os
import sys
import re
import numpy as np

def file_open(fname, line = False): #open the file line = True readlines line = false read entire file as string 
	with open (fname, "r") as file:
		if line:
			each_line = file.readlines()
			return each_line
		all_line = file.read()
		return all_line

def open_values(seq_list, foldn = "rmsf-msl"): 
	cwd = os.getcwd()
	path = cwd +"\\"+ foldn
	#print(path)
	list_files = os.listdir(path)
	#print(a)
	_ret_val = {}
	for i in seq_list:
		for j in list_files:
			if (i[4:]) in j:
				t = file_open(path+"\\"+j, line = True)
				_tmp = {}
				for k in range(1,len(t)+1):
					_tmp[k]=t[k-1][8:].strip("\n")
				_ret_val[i] = _tmp
	return(_ret_val)
	
def read2dat(m_input, START = 26837, LENGTH = 618, END = 30):
	_tmp = {}
	_tmplist = []
	_tmpt, _tmpl = sequ(m_input)
	for i in range(0, 43):
		_tmplist.append(_tmpl[i])
		_tmp[_tmplist[i]] = _tmpt[_tmplist[i]]
	for i in range(0,END):
		_tmp[m_input[START+6*(i+1)+i*LENGTH:START+6*(i+1)+i*LENGTH+11]] = parse(m_input[START+6*(i+1)+i*LENGTH+12:START+6*(i+1)+(i+1)*LENGTH])
		_tmplist.append(m_input[START+6*(i+1)+i*LENGTH:START+6*(i+1)+i*LENGTH+11])
	
	return _tmp, _tmplist
	
"""def read_dat(x, START = 26843, LENGTH = 615, END = 30):
	_tmp = {}
	_tmplist = []
	b = 0
	c = 0
	for i in range(0,END):
		t = x[START+LENGTH*i+6*i+b+c:START+LENGTH*i+6*i+12+b+c]
		a = re.search("TYR_(.*)/",t)
		tmp = a.group()
		b = len(tmp)-8	
		if b > 0:
			c += b
			b = 0
		_tmp[tmp[:-1]] = x[START+LENGTH*i+6*i+b+c+len(tmp)+c:START+LENGTH*(1+i)+6*i+b+c]
		_tmplist.append(tmp[:-1])
	return _tmp, _tmplist """
		
"""def filt_dat(): #filter the data in two sequence and written the position
	for i in range(len(_t)):
		if _t[i] != "-" and _t[i+1] != "-":
			if _r[i] != "-" and _r[i+1] != "-":
				print(_d.index(i))
"""			
def chec_numseq(seq): # to check the number of sequences (1byg)in the file for eg; pdb:1byg 
	count = 0
	for i in range(len(seq)):
		if seq[i] != "-":
			count+=1
	return count 

def print_numseq(seq, pos): #print position and sequence for verification purpose
	for i in range(1, len(pos)+1):
		print(i, seq[pos[i-1]-1])

def final_match(seq, pos = False): # works only for 1byg needs configuration for other 43 sequences
	rfile = reduce2single(seq)
	nfile = crea_numseq(seq)
	mfile, seqlist = sequ(seq)
	vrmsf = open_values(seqlist)
	_rmsfmatrix = {}
	_rmsfposition = {}
	for i in nfile:
		_tmplist1 = [] #postion
		_tmplist2 = [] #value matching rmsf-msl
		for j in rfile:
			_tmplist1.append(nfile[i][j])
		for j in _tmplist1:
			_tmplist2.append(vrmsf[i][j])
		_rmsfmatrix[i] = _tmplist2
		if pos:
			_rmsfposition[i] = _tmplist1
	if pos:
		return _rmsfmatrix, _rmsfposition
	return _rmsfmatrix

def make_matrix(seq, l):
	arr = []
	for i in l:
		tmp = np.zeros(len(seq[i]))
		for j in range(len(seq[i])):
			tmp[j] = float(seq[i][j])
		arr.append(tmp)
	return(np.array(arr))

def save_matrix(fname, matrix, type="%0.4f"):
	np.savetxt(fname, matrix, type, delimiter = ",")
	return
	
def crea_numseq(seq): # do not disturb to create the sequence index in the file for eg: pdb:1byg
	_seq = {}
	mfile, seqlist = sequ(seq)
	for i in seqlist:
		_ret_val = {}
		var = init(mfile[i]) #file_name
		for j in range(1, len(var)+1):
			_ret_val[var[j-1]] = j
		_seq[i] = _ret_val
	return _seq

def init(x): # to remove the "-" in the file for eg: remove "-" in pdb:1byg
	_tmp = []
	for i in range(1,len(x)):
		if x[i] != "-":
			_tmp.append(i+1)
	return _tmp

def parse(seq, chart = "\n"): # to remove the "\n" from the file
    _tmp = []
    for i in seq:
        if i != chart:
            _tmp.append(i)
    return "".join(_tmp)

def reduce2single(x):
	a = _redu2single(x)
	b = _filter(a)
	return b[0]
	
def _redu2single(x):        # 1st level
	_ret_val = {}
	t, a = read2dat(x)
	for i in range(0,len(t)-1, 2):
		_ret_val[i] = compare(t[a[i]],t[a[i+1]])
	_ret_val[i+2] = compare(t[a[0]],t[a[i+2]])
	return _ret_val
	
"""def _redu2single(x):        # 1st level
	_tmp = {}
	t, a = sequ(x)
	r, b = read2dat(x)
	for i in range(0,len(t)-1, 2):
		_tmp[i] = compare(t[a[i]],t[a[i+1]])
	_tmp[i+2] = compare(t[a[0]],t[a[i+2]])
	return _tmp
"""	
def _filter(fn):   			#  2nd level to one value 
	_dat = None
	for i in range(0, 5):    # it will iterate 6 times to reduce the value from 37 to 1
		_dat = filter_dat(fn, 2)
		fn = _dat
	return _dat
	
def filter_dat(fn, node):
	_tmp = {}
	for i in range(0, (len(fn)-1), node):
		_dat = [k for k in fn[node*i] if k in fn[node*i+node]]
		_tmp[i] = _dat
	return _tmp

def compare(seq_1, seq_2):
	_tmp, count=_compare_id(seq_1,seq_2)
	return _tmp

def sequ(x, t = 43, LENGTH = 595, START = 158): # collecting data from sequence file.. do not disturb
	_tmp = {}
	_tmplist = []
	a = 0
	for i in range(0,t):
		if i > 8:
			_dat = parse(x[START+(i)*LENGTH+a+(i-8):START+(i+1)*LENGTH+a+10+(i-8)])			
			_val = x[START-10+(i)*LENGTH+(i)*25+(i-8):START-10+8+(i)*LENGTH+(i)*25+(i-8)] 
			_tmp[_val] = _dat
			a+=25
		#	print(x[149+i*LENGTH+i*25+(i-8):149+8+i*LENGTH+i*25+(i-8)])
		else:
			_dat = parse(x[START+(i)*LENGTH+a:START+(i+1)*LENGTH+a+10])
			_val = x[START-10+(i)*LENGTH+(i)*25:START-10+8+(i)*LENGTH+(i)*25]
			_tmp[_val] = _dat
			a+=25
		#	print(x[149+i*LENGTH+i*25:149+8+i*LENGTH+i*25])
		#sprint(check_numseq(_dat))
		_tmplist.append(_val)
	return _tmp, _tmplist

def _compare_id(seq_1, seq_2): #for testing sequence and the length
	count = 0
	_tmp = []
	for i in range(1, len(seq_1)+ 1):
		if seq_1[i-1] != "-":
			if seq_2[i-1] != "-":
				_tmp.append(i)
				count+=1
	return _tmp, count
