#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import math
import os
import sys
import getopt
from joblib import Parallel, delayed
from datetime import datetime

def read_chr_list(chr_list_file):
    chrs=[]
    f = open(chr_list_file)
    for line in f:
        line = line.rstrip('\n').split('\t')
        chr = line[0]
        if len(chr) != 0:
            chrs.append(chr)
    f.close()
    return chrs

def make_frag_space(fragment_path, chrs_list, window_frag, step_frag):
    frag_space = {x: [] for x in chrs_list}
    f = open(fragment_path)
    for line in f:
        line = line.rstrip('\n').split('\t')
        if line[0] in frag_space:
            chr, start, end = line[0], int(line[1]), int(line[2]) - 1
            frag_space[chr].append([start, end])
    f.close()
    frag_space_slid = {}
    for chr in frag_space:
        fs = frag_space[chr]
        fs = sorted(fs, key = lambda x: (x[0], x[1]))
        n = len(fs)
        fs_one_chr = []
        i = 0
        start_idx = 0
        end_idx = window_frag - 1
        while end_idx < n:
            fs_one_chr.append([fs[start_idx][0], fs[end_idx][1]])
            i += 1
            start_idx = i * step_frag
            end_idx = i* step_frag + window_frag -1
        fs_one_chr.append([fs[start_idx][0], fs[n - 1][1]])
        frag_space_slid[chr] = fs_one_chr
    return frag_space, frag_space_slid

def BinarySearch(point, region_list):
    low = 0
    height = len(region_list)-1
    while low <= height:
        mid = (low+height)/2
        if region_list[mid][1] < point:
            low = mid + 1
        elif region_list[mid][0] > point:
            height = mid - 1
        else:
            return mid
    return -1

def map_pets_to_frag_intra(one_chr_pets_file, fs_one_chr, fs_slid_one_chr, window_frag, step_frag):
    slid_over_windows = window_frag / step_frag
    bin_num = len(fs_slid_one_chr)
    frag_matrix = np.zeros((bin_num, bin_num), dtype='int32')
    f = open(one_chr_pets_file)
    for line in f:
        line = line.rstrip('\n').split('\t')
        pos1 = int(line[1])
        pos2 = int(line[3])
        idx1 = BinarySearch(pos1, fs_one_chr)
        idx2 = BinarySearch(pos2, fs_one_chr)
        idx1_end = idx1 / step_frag
        idx1_end = min(idx1_end, bin_num -1)
        idx1_start = idx1_end - slid_over_windows + 1
        idx1_start = max(0, idx1_start)
        idx2_end = idx2 / step_frag
        idx2_end = min(idx2_end, bin_num -1)
        idx2_start = idx2_end - slid_over_windows + 1
        idx2_start = max(0, idx2_start)
        idxes = [[y, x] for x in range(idx1_start, idx1_end + 1) for y in range(idx2_start, idx2_end + 1) if y>=x]\
                + [[x,y] for x in range(idx1_start, idx1_end + 1) for y in range(idx2_start, idx2_end + 1) if y<x]
        idxes_row = [x[0] for x in idxes]
        idxes_col = [x[1] for x in idxes]
        idxes = (np.array(idxes_row), np.array(idxes_col))
        frag_matrix[idxes] += 1
    f.close()
    return frag_matrix

###step2 gpava
def weighted_fractile(y, w, p):
    #### y and w is 1d np.array
    a = 1 - p
    b = p
    ox = np.argsort(y)
    y = y[ox]
    w = w[ox]
    low = np.cumsum(np.hstack((np.array([0]), w)))
    up = np.sum(w) - low
    df = a * low - b * up
    for k in range(len(df)):
        if df[k] < 0:
            pass
        elif df[k] == 0:
            return (float(w[k]*y[k]+w[k-1]*y[k-1]))/(w[k]+w[k-1])
        else:
            return y[k-1]

def gpava(y, p):
    idx_start = 0
    n = y.size
    u, indices = np.unique(y, return_inverse=True)
    y1 = np.zeros(n)
    while idx_start < n:
        sub_indices = indices[idx_start:n]
        min_inice = min(sub_indices)
        idx_num = np.where(sub_indices == min_inice)[0][-1] + 1
        idx_end = idx_start + idx_num
        a = weighted_fractile(y[idx_start:idx_end], np.ones(idx_num), p)
        y1[idx_start:idx_end] = a
        idx_start = idx_end
    return y1

def matrix_gpava(matrix, p, pet_cut, ratio_cut):
    ratio_select = []
    back_select = []
    idxes_select = []
    pets_select = []
    for i in xrange(len(matrix)):
        y = matrix[i]
        y1 = y[:i]
        b1 = gpava(y1, p)
        b1_ratio = np.zeros(b1.size)
        idx_zero_no = np.where(b1 != 0)[0]
        b1_ratio[idx_zero_no] = y[idx_zero_no]/b1[idx_zero_no]
        pets_cut_idx_zero = np.where(y1 < pet_cut)
        b1_ratio[pets_cut_idx_zero] = 0
        select_idxes = np.where(b1_ratio >= ratio_cut)[0]
        idxes_select += [[j,i] for j in select_idxes]
        ratio_select += b1_ratio[select_idxes].tolist()
        back_select += b1[select_idxes].tolist()
        pets_select += y[select_idxes].tolist()
    ###sort by idx
    idxes_select = [idxes_select[x] + [x] for x in xrange(len(idxes_select))]
    idxes_select = sorted(idxes_select, key=lambda x: (x[0], x[1]))
    idx_sorted = [x[2] for x in idxes_select]
    idx_sorted = np.array(idx_sorted)
    pets_select = np.array(pets_select)[idx_sorted]
    ratio_select = np.array(ratio_select)[idx_sorted]
    back_select = np.array(back_select)[idx_sorted]
    idxes_select = [idxes_select[x][:2] for x in xrange(len(idxes_select))]
    return pets_select, ratio_select, back_select, idxes_select

def write_output(output, pets_select, ratio_select, back_select, idxes_select, fs_one_chr, chr):
    f = open(output, 'w')
    for i in xrange(len(pets_select)):
        idx1 = idxes_select[i][0]
        idx2 = idxes_select[i][1]
        line = [chr, fs_one_chr[idx1][0], fs_one_chr[idx1][1], fs_one_chr[idx2][0], fs_one_chr[idx2][1],
                chr, pets_select[i], ratio_select[i], back_select[i]] + idxes_select[i]
        f.write('\t'.join([str(x) for x in line]) + '\n')
    f.close()

### step3 select overlap rectangle
def merge_windows(regions_list,overlap_windows,pets_select):
    idx_merge = [regions_list[0]]
    for i in xrange(1, len(regions_list)):
        idx_last = regions_list[i-1]
        idx = regions_list[i]
        idx_max = idx_merge[-1]
        if idx[0] - idx_last[0] <= overlap_windows and idx[1] - idx_last[1] <= overlap_windows:
            if pets_select[idx_max[2]] < pets_select[idx[2]]:
                idx_merge[-1] = idx
        else:
            idx_merge.append(idx)
    return idx_merge

def select_recetangle(idxes_select, window_frag, step_frag, pets_select):
    overlap_windows = math.ceil(window_frag / step_frag) - 1
    idxes_select = [idxes_select[x] + [x] for x in range(len(idxes_select))]
    idx_merge = [x for x in idxes_select if x[1]-x[0]>=50]
    while 1:
        idx_merge = sorted(idx_merge, key=lambda x : (x[0], x[1]))
        num1 = len(idx_merge)
        idx_merge = merge_windows(idx_merge,overlap_windows, pets_select)
        num2 = len(idx_merge)
        if num1 == num2:
            break
    return idx_merge

def write_recetangle_select(pets_select, ratio_select, back_select, idx_merge, fs_slid_one_chr, chr, output):
    f = open(output, 'w')
    f.write("chr\tstart\tend\tchr\tstart\tend\tpets_num\tratio\tback_pets\n")
    for i in range(len(idx_merge)):
        idx1 = idx_merge[i][0]
        idx2 = idx_merge[i][1]
        idx_merge_select = idx_merge[i][2]
        region1 = fs_slid_one_chr[idx1]
        region2 = fs_slid_one_chr[idx2]
        pet_num = pets_select[idx_merge_select]
        ratio = ratio_select[idx_merge_select]
        back = back_select[idx_merge_select]
        line = [chr] + region1 + [chr] + region2 + [pet_num, round(ratio, 4), back]
        f.write('\t'.join([str(x) for x in line]) + '\n')
    f.close()

def step123_main(pets_file, fs_one_chr, fs_slid_one_chr, window_frag, step_frag, p, pet_cut, ratio_cut,output,chr):
    print datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Processing %s: map pets to matrix!" % chr
    matrix = map_pets_to_frag_intra(pets_file, fs_one_chr, fs_slid_one_chr, window_frag, step_frag)
    print datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Processing %s: Gpava!" % chr
    pets_select, ratio_select, back_select, idxes_select = matrix_gpava(matrix, p, pet_cut, ratio_cut)
    print datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Processing %s: rectangle overlap!" % chr
    idx_merge = select_recetangle(idxes_select, window_frag, step_frag, pets_select)
    write_recetangle_select(pets_select, ratio_select, back_select, idx_merge, fs_slid_one_chr, chr, output)
    print datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Finished %s!" % chr

def step2_main(pets_select, ratio_select, back_select, idxes_select, fs_slid_one_chr, chr, output, overlap_windows):
    idx_merge = select_recetangle(idxes_select, overlap_windows, pets_select)
    write_recetangle_select(pets_select, ratio_select, back_select, idx_merge, fs_slid_one_chr, chr, output)

def usage():
    print 'Usage: HiC loop calling first step: map pets to fragment space.'
    print '     -i          The input pets file dir.Pets file name format:\"chr1.vs.chr2.pets\".'
    print '                 Pets format:chr1, pos1, chr2, pos2.'
    print '     -f          The REsite fragment bed file.'
    print '     -o          The output dir.'
    print '     -c          The chromosome list file path to select to process chrs.'
    print '     -w          The fragment to make window. Default: 50.'
    print '     -s          The step fragment of window.Default:5.'
    print '     -p          The cpu used. Default:1. Warning: no more than5!'

if __name__ == '__main__':
    start_time = datetime.now()
    cpu = 1
    window_frag = 50
    step_frag = 5
    p = 0.25
    pet_cut = 30
    ratio_cut =1.2
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:f:o:c:w:s:p:')
        for o, a in opts:
            if o == '-h':
                usage()
                sys.exit()
            if o == '-i':
                input_dir = a
            if o == '-f':
                fragment_path = a
            if o == '-o':
                output_dir = a
            if o == '-c':
                chr_list_path = a
            if o == '-w':
                window_frag = int(a)
            if o == '-s':
                step_frag = int(a)
            if o == '-p':
                cpu =int(a)
    except getopt.GetoptError:
        print 'Error in getting parametres!'
        sys.exit()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_files = os.listdir(input_dir)
    pets_names = [file for file in all_files if '.vs.' in file and file.endswith('.pets')]
    input_pets_pathes = []
    output_loop_pathes = []
    chrs_list = read_chr_list(chr_list_path)
    for chr in chrs_list:
        file = [f for f in pets_names if chr + '.' in f][0]
        if len(file) == 0:
            print "Waring: %s/%s.vs.%s.matrix does not exist!" %(input_dir, chr, chr)
        input_pets_pathes.append(input_dir + os.sep + file)
        output_loop_pathes.append(output_dir + os.sep + chr + '_loop.txt')

    print datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Strat make fragment space for all chr!"
    fs_all, fs_slid_all = make_frag_space(fragment_path, chrs_list, window_frag, step_frag)

    Parallel(n_jobs=cpu)(delayed(step123_main)(input_pets_pathes[x], fs_all[chrs_list[x]], fs_slid_all[chrs_list[x]],
                                    window_frag, step_frag, p, pet_cut, ratio_cut,output_loop_pathes[x], chrs_list[x])
                         for x in range(len(chrs_list)))

    time_caused = datetime.now() - start_time
    print "The time caused is : ", time_caused, '\n'
