import sys
import os

import json
import numpy as np
from scipy import stats
from scipy.signal import find_peaks_cwt

if len(sys.argv)<>2:
	sys.exit("python {} config".format(sys.argv[0]))

class ECDF:
	def __init__(self, observations):
		self.observations = observations
		self._p = self._calp

	def _calp(self, x):
		counter = 0.0
		for obs in self.observations:
			if obs <= x:
				counter += 1
		return counter / len(self.observations)
	
	def __call__(self,x):
		return self._p(x)

def pick_peak(dat,step):
	dat.sort()
	maxv = np.max(dat)
	intervar = np.arange(0,maxv,step)
	k = 0
	density = []
	meand = []
	for i in range(1,len(intervar)):
		n = 0
		t = []
		for j in range(0,len(dat)):
			if dat[j]>intervar[i-1] and dat[j]<=intervar[i]:
				t.append(dat[j])
				n+=1
		if n==0:
			meand.append(intervar[i])
		else:
			meand.append(np.mean(t))
		density.append(n)
	ind = density.index(np.max(density))
	return meand[ind]


def config_reader(configfile):
	with open(configfile) as f:
		param_d = json.load(f)
	return param_d


def target_hapid_reader(targetfile):
	target_hapid = []
	with open(targetfile) as f:
		for line in f.readlines():
			target_hapid.append(line.strip())
	return target_hapid


def raw_ab_reader(abfile):
	ab = {}
	window = {}
	with open(abfile) as f:
		w = None
		for line in f.readlines():
			if "sites" in line:
				continue
			ii = line.strip().split()
			hapid = "{}_{}".format(ii[0],ii[1])
			if "mixed" in ii[1]:
				w = ii[0]
				window[w]={}
			window[w][hapid]=None
			a = float(ii[2])
			b = float(ii[5])
			if b==0:
				b=0.0001037452
			ab[hapid]=[a,b]
	return ab,window


def fix_a(indir,samplefile,strand,ab):
	hap_ratio = {}
	with open(samplefile) as f:
		for line in f.readlines():
			ii = line.strip().split()
			haplotype_filename = "{}/{}_{}.readbed.gz_haplotype_{}.txt".format(indir,ii[0],strand,strand)
			with open(haplotype_filename) as h:
				for line_h in h.readlines():
					jj = line_h.strip().split()
					if jj[-1] == "nan":
						continue
					hapid = "{}_{}".format(jj[0],jj[1])
					if hapid in hap_ratio:
						hap_ratio[hapid].append(float(jj[-1]))
					else:
						hap_ratio[hapid]=[]
						hap_ratio[hapid].append(float(jj[-1]))
	fix_a = {}
	for h in hap_ratio:
		m=np.median(hap_ratio[h])
		if m>ab[h][0]:
			fix_a[h]=m
		else:
			fix_a[h]=ab[h][0]
	return fix_a, hap_ratio



def select_one_haplotype_by_window(hap_ratio,fix_a,window):
	h_fix_a = {}
	'''
	for w in window:
		cv = []
		hh = []
		for h in window[w]:
			if h in hap_ratio:
				r = hap_ratio[h]
				m = np.mean(r)
				s = np.std(r)
				cv.append(m/s)
				hh.append(h)
		if hh:
			minh = hh[cv.index(min(cv))]
			h_fix_a[minh] = fix_a[minh]
	'''
	for h in fix_a:
		if "mixed" in h:
			h_fix_a[h]=fix_a[h]
	return h_fix_a


def fetus_distribution(target_hapid,sample,indir,strand,h_fix_a,ab): # get the weight of each window
	sample_filename = "{}/{}_{}.readbed.gz_haplotype_{}.txt".format(indir,sample,strand,strand)
	non_cnv_fetus_ratio = {}
	fetus_ratio = []
	cnv_hap_count = {}
	cnv_hap_noncount = {}
	with open(sample_filename) as f:
		for line in f.readlines():
			ii = line.strip().split()
			chrid = ii[0].split(";")[0]
			hapid = "{}_{}".format(ii[0],ii[1])
			if hapid in target_hapid:
				cnv_hap_count[hapid] = int(ii[-2])
				cnv_hap_noncount[hapid] = int(ii[-3])-int(ii[-2])
				continue
			if hapid in h_fix_a:
				f = (float(ii[-1])-ab[hapid][1])/(h_fix_a[hapid]-ab[hapid][1])
				if f<0:
					continue
				if f>=1:
					continue
				non_cnv_fetus_ratio[hapid] = f
				fetus_ratio.append(f)
	ecdf = ECDF(fetus_ratio)
	fetus_ratio = np.array(fetus_ratio)
	#peakind = find_peaks_cwt(fetus_ratio, np.arange(9,10))
	#peak = np.mean(fetus_ratio[peakind])
	peak = pick_peak(fetus_ratio,0.01)
	non_cnv_fetus_ratio_p = {}
	allp = 0
	alln = 0
	tmpr = []
	for h in non_cnv_fetus_ratio:
		tmpr = non_cnv_fetus_ratio[h]
		non_cnv_fetus_ratio_p[h]=1-abs(ecdf(tmpr)-ecdf(peak))
		allp+=1-abs(ecdf(tmpr)-ecdf(peak))
		alln+=1
	non_cnv_fetus_ratio_w = {}
	for h in non_cnv_fetus_ratio:
		non_cnv_fetus_ratio_w[h]=non_cnv_fetus_ratio_p[h]/allp
	
	return non_cnv_fetus_ratio_p,non_cnv_fetus_ratio_w,non_cnv_fetus_ratio,cnv_hap_count,cnv_hap_noncount


def genotyping(cnv_fetus_ratio_p,cnv_fetus_ratio,cnv_hap_count,cnv_hap_noncount,h_fix_a,ab):
	N = np.arange(0,2.1,0.1)
	P = {}
	for cnvh in cnv_hap_count:
		for n in N:
			p = 0
			a = h_fix_a[cnvh]
			b = ab[cnvh][1]
			cnvcount = cnv_hap_count[cnvh]
			cnvnoncount = cnv_hap_noncount[cnvh]
			for noncnvh in cnv_fetus_ratio_p:
				w = cnv_fetus_ratio_p[noncnvh]
				f = cnv_fetus_ratio[noncnvh]
				F = n*f/(1-f)
				ph = a*F/(F+1)+b*(1-F/(F+1))
				pnh = 1-ph
				p += w*(ph**cnvcount)*(pnh**cnvnoncount)
			if cnvh in P:
				P[cnvh][n]=p
			else:
				P[cnvh]={}
				P[cnvh][n]=p
	return P


if __name__ == '__main__':
	PARA = config_reader(sys.argv[1])
	target_hapid = target_hapid_reader(PARA["target_hapid"])
	(ab_plus,window_plus) = raw_ab_reader(PARA["raw_hap_plus"])
	(ab_minus,window_minus) = raw_ab_reader(PARA["raw_hap_minus"])
	(fix_a_plus, hap_ratio_plus) = fix_a(PARA["villus_input"],PARA["villus_sample"],"plus",ab_plus)
	(fix_a_minus, hap_ratio_minus) = fix_a(PARA["villus_input"],PARA["villus_sample"],"minus",ab_minus)
	h_fix_a_plus = select_one_haplotype_by_window(hap_ratio_plus,fix_a_plus,window_plus)
	h_fix_a_minus = select_one_haplotype_by_window(hap_ratio_minus,fix_a_minus,window_minus)
	L = open(PARA["output"],"w")
	N = np.arange(0,2.1,0.1)
	L.write("sample\thapid")
	for nn in N:
		L.write("\t"+str(nn))
	L.write("\tgenotype\n")
	with open(PARA["pregnant_plasma_sample"]) as f:
		for line in f.readlines():
			ii = line.strip().split()
			sample = ii[0]
			(cnv_fetus_ratio_p_plus,cnv_fetus_ratio_w_plus,cnv_fetus_ratio_plus,cnv_hap_count_plus,cnv_hap_noncount_plus) = fetus_distribution(target_hapid,sample,PARA["pregnant_plasma_input"],"plus",h_fix_a_plus,ab_plus)
			(cnv_fetus_ratio_p_minus,cnv_fetus_ratio_w_minus,cnv_fetus_ratio_minus,cnv_hap_count_minus,cnv_hap_noncount_minus) = fetus_distribution(target_hapid,sample,PARA["pregnant_plasma_input"],"minus",h_fix_a_minus,ab_minus)
			P_plus = genotyping(cnv_fetus_ratio_p_plus,cnv_fetus_ratio_plus,cnv_hap_count_plus,cnv_hap_noncount_plus,h_fix_a_plus,ab_plus)
			P_minus = genotyping(cnv_fetus_ratio_p_minus,cnv_fetus_ratio_minus,cnv_hap_count_minus,cnv_hap_noncount_minus,h_fix_a_minus,ab_minus)
			
			for h in P_plus:
				L.write("{}\t{}".format(sample,h))
				pp = []
				for n in N:
					L.write("\t{}".format(str(P_plus[h][n])))
					pp.append(P_plus[h][n])
				L.write("\t"+str(N[pp.index(np.max(pp))]))
				L.write("\n")
