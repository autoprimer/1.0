#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 17:23:31 2017

@author: HCC604
"""
import atexit
import hashlib
import math
import json
import re
import sys
import time

import urllib.parse
import urllib.request

debug = False
version = '1.0g build 20180502'
message = 'Updated batch submission mode'

def print_disclaimer():
    """Print disclaimer at exit"""
    print('# AutoPrimer program by Tom C.C. Ho (c) 2017', file=sys.stderr)
    print('# Version ' + version, file=sys.stderr)
    print('# Release note: ' + message, file=sys.stderr)
    if debug:
        print('# Debug mode is ON')

atexit.register(print_disclaimer)

def check_url(url):
    """Check server status using URL. Return HTTP status code."""
    print('# Now checking status of -', file=sys.stderr)
    print(url, file=sys.stderr)
    code = None
    try:
        code = urllib.request.urlopen(url).getcode()
    except:
        code = False
    return code

def jget(url):
    """Return a dictionary of parsed JSON info from Ensembl"""
    d = None
    with urllib.request.urlopen(url) as response:
        json_data = response.read().decode('ascii')
        d = json.loads(json_data)
    return d

def wget(url):
    """Return plain text response from a URL"""
    d = None
    passed = False
    while not passed:
        try:
            with urllib.request.urlopen(url) as response:
                d = response.read().decode('ascii')
            passed = True
        except:
            print('# Waiting...', file=sys.stderr)
            time.sleep(5)
    return d

def autoSNPfree(pair_name, fp, rp, server, max_wait=1200, batch_mode=False):
    """Perform online SNPfree checking and return an objective score"""
    # prepare the query
    if server == 'HKU':
        base_url = 'http://grass.cgs.hku.hk/SNPfree/SNPfree.php?primers='
    if server == 'PYN':
        base_url = 'http://ebensee.hopto.org:8082/pyneh/SNPfree/SNPfree.php?primers='
    if server == 'CYK':
        base_url = 'http://vonbismarck.hopto.org:8080/SNPfree/SNPfree.php?primers='
    if batch_mode:
        query = ''
        for i in range(len(pair_name)):
            query += '\t'.join([pair_name[i], fp[i], rp[i]]) + '\n'
            query = urllib.parse.quote_plus(query)
    else:
        query = urllib.parse.quote_plus(pair_name + '\t' + fp + '\t' + rp)

    # submit the query and retrieve the jobid
    jobid_page = wget(base_url + query)
    jobid_found = re.search('jobid=([0-9a-f]{40})', jobid_page)
    if jobid_found:
        print('# Online SNPfree job successfully created.', file=sys.stderr)
        jobid = jobid_found.group(1)
        if server == 'HKU':
            result_url = 'http://grass.cgs.hku.hk/SNPfree/view.php?jobid='
        if server == 'PYN':
            result_url = 'http://ebensee.hopto.org:8082/pyneh/SNPfree/view.php?jobid='
        if server == 'CYK':
            result_url = 'http://vonbismarck.hopto.org:8080/SNPfree/view.php?jobid='
        ready = False
        # wait for the result
        waited = 0
        while not ready:
            result_page = wget(result_url + jobid)
            if 'You will be redirected' in result_page:
                if waited > max_wait:
                    print('# Maximum waiting time exceeded. Aborting SNPfree checking.', file=sys.stderr)
                    return None
                else:
                    time.sleep(5)
                    waited += 5
            else:
                ready = True
        return result_page

def sp_penalty(ratio, mism, gaps):
    """Return the penalized score given a non-specific product"""
    if ratio <= 1:
        deduction = math.log10(1 + ratio * 10)
        if deduction > 1:
            deduction = 1
    else:
        deduction = math.exp(-ratio)/math.exp(-1)
    if mism > 20:
        mism = 20
    deduction = deduction * (1 - mism/20)
    if gaps > 10:
        gaps = 10
    deduction = deduction * (1 - gaps/10)
    return deduction

def snp_penalty(snp_count, threeprime=False):
    """Return the penalized score given a SNP count"""
    if threeprime:
        if snp_count > 10:
            snp_count = 10
        deduction = snp_count / 10
    else:
        if snp_count > 20:
            snp_count = 20
        deduction = snp_count / 20
    return deduction

def SNPfree_to_score(autoSNPfree_result, get_size=False):
    """Parse the online SNPfree result"""
    autoSNPfree_result = autoSNPfree_result.split('\n')
    primer_name = None
    # parameters in calculating primer score
    product_size = list()
    product_mism = list()
    product_gaps = list()
    fp_total_snps = None
    fp_3prime_snps = None
    rp_total_snps = None
    rp_3prime_snps = None
    # exclude QC results
    query_result = list()
    append_to_result = False
    for line in autoSNPfree_result:
        primer_name_found = re.search('^Primer name: (.*)$', line)
        if primer_name_found and primer_name_found.group(1) != 'QC':
            append_to_result = True
            primer_name = primer_name_found.group(1)
        terminator_string_found = re.search('^-----$', line)
        if terminator_string_found:
            append_to_result = False
        if append_to_result:
            query_result.append(line)
    for line in query_result:
        hits_tag_found = re.search('^<H:([^>]+)>$', line)
        if hits_tag_found:
            hits = hits_tag_found.group(1).rstrip('|').split('|')
            print(hits, file=sys.stderr)
            for i in range(len(hits)):
                size, mism, gaps = hits[i].split(',')
                product_size.append(int(size))
                product_mism.append(int(mism))
                product_gaps.append(int(gaps))
        fpsnp_tag_found = re.search('^<SF:([^>]+)>$', line)
        if fpsnp_tag_found:
            fpsnp_data = fpsnp_tag_found.group(1).split(',')
            fp_total_snps, fp_3prime_snps = int(fpsnp_data[0]), int(fpsnp_data[1])
        rpsnp_tag_found = re.search('^<SR:([^>]+)>$', line)
        if rpsnp_tag_found:
            rpsnp_data = rpsnp_tag_found.group(1).split(',')
            rp_total_snps, rp_3prime_snps = int(rpsnp_data[0]), int(rpsnp_data[1])
    # calculate the objective SNPfree score
    score = 100
    # non-specific priming penalty
    pri_prod_size = product_size[0]

    # special function to return primary product size only
    if get_size:
        return pri_prod_size
    if len(product_size) > 1:
        print('#', len(product_size) - 1, 'non-specific products found by SNPfree')
        for i in range(1, len(product_size)):
            sec_prod_size = product_size[i]
            sec_prod_mism = product_mism[i]
            sec_prod_gaps = product_gaps[i]
            size_ratio = sec_prod_size / pri_prod_size
            score = score * (1 - sp_penalty(size_ratio, sec_prod_mism, sec_prod_gaps))
    # SNP penalty
    if fp_total_snps > 0 or  rp_total_snps > 0:
        print('#', fp_total_snps, 'and', rp_total_snps, 'SNPs on FP and RP respectively')
    score = score * (1 - snp_penalty(fp_total_snps + rp_total_snps, threeprime=False))
    if fp_3prime_snps > 0 or  rp_3prime_snps > 0:
        print('#', fp_3prime_snps, 'and', rp_3prime_snps, "SNPs within 5 bp form 3' end of FP and RP respectively")
    score = score * (1 - snp_penalty(fp_3prime_snps + rp_3prime_snps, threeprime=True))
    score = round(score, 1)

    return (primer_name, score)



def get_sequence(region, mask='soft'):
    """Return the sequence as specified by region"""
    chromosome = region[0]
    start = region[1]
    end = region[2]
    strand = region[3]
    base_url = 'http://rest.ensembl.org/sequence/region/human/'
    url = base_url + str(chromosome) + ':' + str(start) + '..' + str(end) + ':' + str(strand)
    url += '?content-type=text/plain;'
    if mask == 'soft':
        url += 'mask=soft'
    elif mask == 'hard':
        url += 'mask=hard'
    sequence = wget(url)
    return sequence

def get_flanking_regions(exon_region, flank=400):
    """Return the 5' and 3' flanking regions"""
    chromosome = exon_region[0]
    start = exon_region[1]
    end = exon_region[2]
    strand = exon_region[3]
    if strand == 1:
        flank_1_start = start - flank
        flank_1_end = start - 1
        flank_2_start = end + 1
        flank_2_end = end + flank
    if strand == -1:
        flank_1_start = end + 1
        flank_1_end = end + flank
        flank_2_start = start - flank
        flank_2_end = start - 1
    return [[chromosome, flank_1_start, flank_1_end, strand], [chromosome, flank_2_start, flank_2_end, strand]]

def vcf_to_dict(vcf_data):
    """Return a prased VCF dictionary for variant query"""
    snp_dict = dict()
    vcf_data = vcf_data.split('\n')
    for line in vcf_data:
        if not line.startswith('#'):
            line = line.rstrip().split('\t')
            if len(line) == 8:
                pos = line[1]
                maf = 'NA'
                m = re.search('CAF=([^;]+);', line[7])
                if m:
                    maf = m.group(1)
                m = re.search('COMMON=1', line[7])
                if m:
                    maf = 'COMMON'
                snp_dict[int(pos)] = maf
    return snp_dict

def webix_get_snps(chrom, start, end, server):
    """Return a parsed SNP dictionary from 'webix'"""
    if server == 'HKU':
        base_url = 'http://grass.cgs.hku.hk/SNPfree/webix.php?'
    if server == 'PYN':
        base_url = 'http://ebensee.hopto.org:8082/pyneh/SNPfree/webix.php?'
    if server == 'CYK':
        base_url = 'http://vonbismarck.hopto.org:8080/SNPfree/webix.php?'
    url = base_url + 'chr=' + str(chrom) + '&'
    url += 'start=' + str(start) + '&'
    url += 'end=' + str(end) + '&'
    query = str(chrom) + ':' + str(start) + '-' + str(end)
    secret_key = '29901946'
    key = hashlib.sha1(query.encode('ascii')).hexdigest()
    key = hashlib.sha1((key + secret_key).encode('ascii')).hexdigest()
    url += 'key=' + str(key)
    vcf_data = wget(url)
    snps = vcf_to_dict(vcf_data)
    return snps

def wrimer3(seq_id, seq, target, server, min_prod=100, max_prod=1200):
    """Return the Primer3 result from 'wrimer3'"""
    if server == 'HKU':
        base_url = 'http://grass.cgs.hku.hk/SNPfree/wrimer3.php?'
    if server == 'PYN':
        base_url = 'http://ebensee.hopto.org:8082/pyneh/AutoPrimer/wrimer3.php?'
    if server == 'CYK':
        base_url = 'http://vonbismarck.hopto.org:8080/AutoPrimer/wrimer3.php?'
    url = base_url + 'sequence_id=' + str(seq_id) + '&'
    url += 'sequence_template=' + str(seq) + '&'
    url += 'sequence_target=' + str(target) + '&'
    url += 'min_prod=' + str(min_prod) + '&'
    url += 'max_prod=' + str(max_prod) + '&'
    file_content = 'SEQUENCE_ID=' + seq_id + '\n'
    file_content += 'SEQUENCE_TEMPLATE=' + seq + '\n'
    file_content += 'SEQUENCE_TARGET=' + target + '\n'
    file_content += 'PRIMER_PRODUCT_SIZE_RANGE=' + str(min_prod) + '-' + str(max_prod) + '\n'
    file_content += 'PRIMER_EXPLAIN_FLAG=1\n'
    file_content += 'PRIMER_NUM_RETURN=10\n'
    file_content += '=\n'
    secret_key = '29901946'
    key = hashlib.sha1(file_content.encode('ascii')).hexdigest()
    key = hashlib.sha1((key + secret_key).encode('ascii')).hexdigest()
    url += 'key=' + str(key)
    print('# Submitting request for ' + seq_id, file=sys.stderr)
    print('# Key =', key, file=sys.stderr)
    print('#', url, file=sys.stderr)
    primer3_result = wget(url)
    return primer3_result

def wrimer3_to_primers(wrimer3_result):
    """Return the primer pairs returned by 'wrimer3'"""
    primers = dict()
    wrimer3_result = wrimer3_result.split('\n')
    for line in wrimer3_result:
        lp = re.search('PRIMER_LEFT_([0-9])_SEQUENCE=([ATCGatcg]+)', line)
        if lp:
            if not lp.group(1) in primers:
                primers[lp.group(1)] = dict()
            primers[lp.group(1)]['FP'] = lp.group(2)
        rp = re.search('PRIMER_RIGHT_([0-9])_SEQUENCE=([ATCGatcg]+)', line)
        if rp:
            primers[rp.group(1)]['RP'] = rp.group(2)
    print('#', len(primers), 'primer pairs returned', file=sys.stderr)
    if len(primers) == 0:
        print(wrimer3_result, file=sys.stderr)
    primer_list = list()
    for i in range(len(primers)):
        primer_list.append([primers[str(i)]['FP'], primers[str(i)]['RP']])
    return primer_list

def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_masked_sequence(region, webix_server, threshold='COMMON', repeat_mask='soft'):
    """Return the masked sequence"""
    chromosome = region[0]
    start = region[1]
    end = region[2]
    sequence = get_sequence(region, mask=repeat_mask)
    snp_dict = webix_get_snps(chromosome, start, end, webix_server)
    sequence_output = ''
    this_pos = start
    for s in sequence:
        if threshold == 'COMMON':
            if this_pos not in snp_dict:
                sequence_output += s
            elif snp_dict[this_pos] == 'COMMON':
                sequence_output += 'N'
            else:
                sequence_output += s
        if threshold == 'ALL':
            if this_pos in snp_dict:
                sequence_output += 'N'
            else:
                sequence_output += s
        if is_numeric(threshold):
            if this_pos not in snp_dict:
                sequence_output += s
            else:
                maf = snp_dict[this_pos]
                if maf == 'COMMON':
                    sequence_output += 'N'
                elif maf == 'NA':
                    sequence_output += s
                else:
                    maf = maf.split(',')
                    aaf = maf[1:]
                    flag = False
                    for af in aaf:
                        if is_numeric(af) and float(af) >= threshold:
                            flag = True
                    if flag:
                        sequence_output += 'n'
                    else:
                        sequence_output += s

        this_pos += 1
    return sequence_output

def is_compliant(translated_regions):
    """Check for compliance before splitting / joining"""
    # check for directionality and compliance to special data format
    direction = None
    for tr in translated_regions:
        if tr[3] and direction == None:
            direction = tr[3]
        if tr[3] != direction:
            print('Exons do not locate on the same strand!', file=sys.stderr)
            return False
        if tr[1] and tr[2] and tr[1] > tr[2]:
            print('Non-compliant data format!', file=sys.stderr)
            return False
    return True

def region_size_stat(translated_regions):
    """Return mean, min and max size of translated regions"""
    sizes = list()
    if is_compliant(translated_regions):
        for tr in translated_regions:
            if tr[1] and tr[2]:
                sizes.append(tr[2] - tr[1] + 1)
        return [sum(sizes)/len(sizes), min(sizes), max(sizes)]
    else:
        return None


def split_big_regions(translated_regions, translated_region_names, chunk=700):
    """Split big regions into smaller chunks for easier sequencing"""
    out_regions = list()
    out_names = list()
    sub_region_suffix = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

    if is_compliant(translated_regions) == False:
        raise ValueError('Non-compliant regions cannot be split.')

    # split region if size > chunk
    for tr, trn in zip(translated_regions, translated_region_names):
        if tr[1] and tr[2]:
            size = tr[2] - tr[1] + 1
            if size > chunk:
                print('#', trn, 'is > ', str(chunk), 'bp. Splitting...', file=sys.stderr)
                remaining_region = size
                if tr[3] == 1:
                    start_coordinate = tr[1]
                    end_coordinate = tr[2]
                    this_chunk = chunk
                    add_to_start = this_chunk - 1
                if tr[3] == -1:
                    start_coordinate = tr[2]
                    end_coordinate = tr[1]
                    this_chunk = -chunk
                    add_to_start = this_chunk + 1
                sub_region_index = 0
                this_start_coordinate = start_coordinate
                while remaining_region > 0:
                    this_end_coordinate = this_start_coordinate + add_to_start
                    if abs(this_end_coordinate - this_start_coordinate) > abs(end_coordinate - this_start_coordinate):
                        this_end_coordinate = end_coordinate
                    if tr[3] == 1:
                        out_regions.append([tr[0], this_start_coordinate, this_end_coordinate, tr[3], tr[4]])
                    if tr[3] == -1:
                        out_regions.append([tr[0], this_end_coordinate, this_start_coordinate, tr[3], tr[4]])
                    out_names.append(trn + sub_region_suffix[sub_region_index])
                    this_start_coordinate += this_chunk
                    remaining_region -= chunk
                    sub_region_index += 1
            else:
                out_regions.append(tr)
                out_names.append(trn)
        else:
            print('#', trn, 'skipped', file=sys.stderr)
            out_regions.append(tr)
            out_names.append(trn)
    return out_regions, out_names

def combine_adj_regions(translated_regions, translated_region_names, chunk=800):
    """Combine adjacent regions for efficient use of primers"""
    out_regions = list()
    out_names = list()
    processed_regions = dict()

    if is_compliant(translated_regions) == False:
        raise ValueError('Non-compliant regions cannot be combined.')

    for i, tr in enumerate(translated_regions):
        if translated_region_names[i] in processed_regions:
            continue
        if (i + 1 < len(translated_regions) and
            translated_regions[i][1] and translated_regions[i][2] and
            translated_regions[i+1][1] and translated_regions[i+1][2]):
            this_region = translated_regions[i]
            next_region = translated_regions[i+1]
            if this_region[3] == 1:
                    new_start = min(this_region[1], this_region[2])
                    new_end = max(next_region[1], next_region[2])
            if this_region[3] == -1:
                new_start = max(this_region[1], this_region[2])
                new_end = min(next_region[1], next_region[2])
            combined_size = abs(new_end - new_start) + 1
            if combined_size <= chunk:
                # combine the contiguous regions
                new_name = translated_region_names[i] + translated_region_names[i+1]
                new_transcript_info = this_region[4] + ';' + next_region[4]
                if this_region[0] != next_region[0] or this_region[3] != next_region[3]:
                    raise ValueError('Regions on different chromosomes / strands cannot be combined.')
                if new_start > new_end:
                    new_start, new_end = new_end, new_start
                out_regions.append([this_region[0], new_start, new_end, this_region[3], new_transcript_info])
                out_names.append(new_name)
                processed_regions[translated_region_names[i]] = True
                processed_regions[translated_region_names[i+1]] = True
            else:
                # combined size is bigger than chunk size
                out_regions.append(tr)
                out_names.append(translated_region_names[i])
        else:
            print('#', translated_region_names[i], 'skipped', file=sys.stderr)
            out_regions.append(tr)
            out_names.append(translated_region_names[i])
    return out_regions, out_names

class Gene:
    """Base container for genetic information from Ensembl"""
    def __init__(self, name):
        self.name = str(name).upper()
        self.get_gene_info()
        self.set_transcript(self.name)

    def get_gene_info(self):
        """Get gene info from Ensembl"""
        base_url = 'http://rest.ensembl.org/lookup/symbol/homo_sapiens/'
        url = base_url + str(self.name).upper() + '?content-type=application/json;expand=1'
        gene_info = jget(url)
        self.gene_info= gene_info
        return 'Info for gene ' + self.name + ' retrieved'

    def list_transcripts(self):
        """Return the list of available transcripts"""
        transcripts = list()
        for t in self.gene_info['Transcript']:
            transcripts.append(t['display_name'])
        return transcripts

    def set_transcript(self, transcript_name):
        """Set the preferred transcript"""
        self.transcript_info = None
        self.transcript_name = None
        for t in self.gene_info['Transcript']:
            if t['display_name'] == transcript_name:
                self.transcript_info = t
                self.transcript_name = transcript_name
                self.transcript_id = t['id']
        return transcript_name + ' set as working transcript'

    def list_exon_regions(self):
        """Return the exonic regions in self.transcript_info"""
        regions = list()
        for e in self.transcript_info['Exon']:
            regions.append([e['seq_region_name'], e['start'], e['end'], e['strand'], e['id']])
        return regions

    def list_translation_boundaries(self):
        """Return the list of translation boundaires"""
        if self.transcript_info['strand'] == 1:
            return [self.transcript_info['Translation']['start'], self.transcript_info['Translation']['end']]
        if self.transcript_info['strand'] == -1:
            return [self.transcript_info['Translation']['end'], self.transcript_info['Translation']['start']]

    def exon_to_translated(self, exon_region):
        """Return the translated region of the exon"""
        t_start, t_end = self.list_translation_boundaries()
        e_start, e_end = exon_region[1], exon_region[2]
        new_start, new_end = e_start, e_end
        if self.transcript_info['strand'] == 1:
            if not ((t_start >= e_start and t_start <= e_end) or
                    (t_end >= e_start and t_end <= e_end)):
                if t_start > e_end or t_end < e_start:
                    return [exon_region[0], None, None, exon_region[3], exon_region[4]]
                return exon_region
            else:
                if (t_start >= e_start and t_start <= e_end):
                    new_start = t_start
                if (t_end >= e_start and t_end <= e_end):
                    new_end = t_end
                return [exon_region[0], new_start, new_end, exon_region[3], exon_region[4]]
        if self.transcript_info['strand'] == -1:
            e_start, e_end = e_end, e_start
            new_start, new_end = new_end, new_start
            if not ((t_start <= e_start and t_start >= e_end) or
                    (t_end <= e_start and t_end >= e_end)):
                if t_start < e_end or t_end > e_start:
                    return [exon_region[0], None, None, exon_region[3], exon_region[4]]
                return exon_region
            else:
                if (t_start <= e_start and t_start >= e_end):
                    new_start = t_start
                if (t_end <= e_start and t_end >= e_end):
                    new_end = t_end
                new_start, new_end = new_end, new_start
                return [exon_region[0], new_start, new_end, exon_region[3], exon_region[4]]

class PrimerPool:
    """Container for generated primers"""
    def __init__(self):
        self.primers = list()

    def add(self, amplicon, primer_id, fp, rp, score, size):
        self.primers.append(
                {'amplicon': amplicon,
                 'primer_id': primer_id,
                 'fp': fp,
                 'rp': rp,
                 'score': float(score),
                 'pp_size': int(size)
                }
                            )
    def getbest(self):
        # get the complete list of amplicons
        amplicons = list()
        for p in self.primers:
            if p['amplicon'] not in amplicons:
                amplicons.append(p['amplicon'])
        # get the best pair for each amplicon
        amplicon_scores = list()
        for a in amplicons:
            score = 0
            for p in self.primers:
                if p['amplicon'] == a and p['score'] >= score:
                    score = p['score']
            amplicon_scores.append([a, score])
        # output the best pair for each amplicon
        output_list = list()
        for amp_score in amplicon_scores:
            a, s = amp_score
            for p in self.primers:
                if p['amplicon'] == a and p['score'] == s:
                    output_list.append([p['amplicon'], p['primer_id'],
                                   p['fp'], p['rp'], p['score'], p['pp_size']])
                    break
        return output_list

def check_servers():
    ''' Return a dictionary of server status '''
    hku_status = check_url('http://grass.cgs.hku.hk/SNPfree/')
    pyn_status = check_url('http://ebensee.hopto.org:8082/pyneh/SNPfree/')
    cyk_status = check_url('http://vonbismarck.hopto.org:8080/SNPfree/')
    rest_status = check_url('https://rest.ensembl.org/')

    status_dict = {'hku_status': hku_status,
                   'pyn_status': pyn_status,
                   'cyk_status': cyk_status,
                   'rest_status': rest_status
                   }

    return status_dict

def select_server():
    server_status = check_servers()
    hku_status = server_status['hku_status']
    pyn_status = server_status['pyn_status']
    cyk_status = server_status['cyk_status']
    rest_status = server_status['rest_status']
    print('Server status -', file=sys.stderr)
    print('HKU', hku_status, file=sys.stderr)
    print('PYNEH', pyn_status, file=sys.stderr)
    print('CYK', cyk_status, file=sys.stderr)
    print('ENSEMBL', rest_status, file=sys.stderr)

    if rest_status == False:
        print('Essential webserver (ENSEMBL) down!', file=sys.stderr)
        return None

    if hku_status != False:
        print('Choosing HKU-CGS server...', file=sys.stderr)
        server = 'HKU'
    else:
        if cyk_status != False:
            print('Choosing CYK (private) server...', file=sys.stderr)
            server = 'CYK'
        else:
            if pyn_status != False:
                print('Choosing PYNEH server...', file=sys.stderr)
                server = 'PYN'
            else:
                print('All SNPfree servers down! Exiting...', file=sys.stderr)
                return None
    return server

def main(args):

    if len(args) < 2:
        print('No gene specified.', file=sys.stderr)
        print('USAGE: autoprimer.py [GENE] [TRANSCRIPT] [MODE] [AMPLICON_SIZE]')
        sys.exit()

    gene = str(args[1]).upper()

    if len(args) < 3:
        g = Gene(gene)
        print('No transcript specified. Choose from below:', file=sys.stderr)
        print(g.list_transcripts(), file=sys.stderr)
        print('USAGE: autoprimer.py [GENE] [TRANSCRIPT] [MODE] [AMPLICON_SIZE]')
        sys.exit()

    transcript = str(args[2]).upper()

    if len(sys.argv) < 4:
        print('Please specify operation mode. Either "manual" or "auto".')
        print('USAGE: autoprimer.py [GENE] [TRANSCRIPT] [MODE] [AMPLICON_SIZE]')
        sys.exit()

    mode = str(args[3]).upper()

    if len(sys.argv) < 5:
        print('Please specify amplicon size. Either "standard" or "small".')
        print('USAGE: autoprimer.py [GENE] [TRANSCRIPT] [MODE] [AMPLICON_SIZE]')
        sys.exit()

    amplicon = str(args[4]).upper()

    ###########################################################################
    ''' SNP and repeat masking settings '''
    routine_SNP_threshold = 0.0001
    retry_SNP_threshold = 'COMMON'
    routine_repeat_masking = 'hard'
    retry_repeat_masking = 'soft'
    ''' Amplicon size, exon-cutting and exon-combining settings '''
    routine_total_flank = 500
    retry_total_flank = 1000
    routine_target_flank = 100
    retry_target_flank = 80
    if amplicon == 'SMALL':
        routine_amplicon_max = 600
        retry_amplicon_max = 800
    elif amplicon == 'STANDARD':
        routine_amplicon_max = 1200
        retry_amplicon_max = 1500
    else:
        raise ValueError('Amplicon size setting cannot be', amplicon)

    ''' Intelligent assay design settings '''
    auto_retry = True
    auto_retry_level = 3
    auto_cut_exons = True
    auto_combine_exons = True

    if amplicon == 'SMALL':
        auto_cut_chunk = 400
        auto_combine_chunk = 500
    elif amplicon == 'STANDARD':
        auto_cut_chunk = 800
        auto_combine_chunk = 900
    else:
        raise ValueError('Amplicon size setting cannot be', amplicon)

    ###########################################################################

    default_server = select_server()
    print('Default server:', default_server, file=sys.stderr)

    if default_server == None:
        print('Please check internet connection!', file=sys.stderr)
        sys.exit()


    g = Gene(gene)
    print(g.list_transcripts(), file=sys.stderr)
    print(g.set_transcript(transcript), file=sys.stderr)
    print(g.transcript_id)

    print('# Now retrieving exon information...', file=sys.stderr)
    exon_regions = g.list_exon_regions()
    print(region_size_stat(exon_regions), file=sys.stderr)
    if debug:
        for i, j in zip(['E' + str(i) for i in range(1, len(exon_regions) + 1)], exon_regions):
            if j[1] and j[2]:
                print('## ORIG', i, j[1], j[2], j[4], sep='\t')
    print('# Now trimming down exons to coding sequences...', file=sys.stderr)
    translated_regions = [g.exon_to_translated(er) for er in exon_regions]
    translated_region_names = ['E' + str(i+1) for i in range(len(translated_regions))]
    print(region_size_stat(translated_regions), file=sys.stderr)
    if debug:
        for i, j in zip(translated_region_names, translated_regions):
            if j[1] and j[2]:
                print('## TRIMMED', i, j[1], j[2], j[4], sep='\t')
    print('# Now scanning for big coding regions...', file=sys.stderr)
    translated_regions, translated_region_names = split_big_regions(translated_regions, translated_region_names, chunk=auto_cut_chunk)
    if debug:
        for i, j in zip(translated_region_names, translated_regions):
            if j[1] and j[2]:
                print('## SPLIT', i, j[1], j[2], j[4], sep='\t')
    print(region_size_stat(translated_regions), file=sys.stderr)
    print('# Now iteratively combining adjacent coding regions...', file=sys.stderr)
    print('#', len(translated_regions), 'regions before combining adjacent coding regions', file=sys.stderr)
    reduction = 1
    while reduction > 0:
        c_translated_regions, c_translated_region_names = combine_adj_regions(translated_regions, translated_region_names, chunk=auto_combine_chunk)
        reduction = len(translated_regions) - len(c_translated_regions)
        if reduction > 0:
            print('#', reduction, 'regions combined...', file=sys.stderr)
            translated_regions, translated_region_names = c_translated_regions, c_translated_region_names
        if reduction < 0:
            raise RuntimeError('An error occured. The number of combined regions is greated than input. Aborting.')
    print('#', len(translated_regions), 'regions after combining adjacent coding regions', file=sys.stderr)
    if debug:
        for i, j in zip(translated_region_names, translated_regions):
            if j[1] and j[2]:
                print('## COMBINED', i, j[1], j[2], j[4], sep='\t')
    print(region_size_stat(translated_regions), file=sys.stderr)

    pp = PrimerPool()

    for translated_region, translated_region_name in zip(translated_regions, translated_region_names):
        if translated_region[1] and translated_region[2]:
            upstream, downstream = get_flanking_regions(translated_region,
                                                        flank=routine_total_flank)
            upseq = get_masked_sequence(upstream,
                                        default_server,
                                        threshold=routine_SNP_threshold,
                                        repeat_mask=routine_repeat_masking)
            coreseq = get_sequence(translated_region)
            downseq = get_masked_sequence(downstream,
                                          default_server,
                                          threshold=routine_SNP_threshold,
                                          repeat_mask=routine_repeat_masking)
            # prepare wrimer3 input
            if mode == 'AUTO':
                sequence_id = g.name + '_' + translated_region_name
                sequence = upseq + coreseq + downseq
                target_start = len(upseq) - routine_target_flank
                target_length = len(coreseq) + 2 * routine_target_flank
                print('# Total length =', len(sequence), file=sys.stderr)
                print('#', target_start, routine_target_flank, len(coreseq), routine_target_flank, routine_total_flank - routine_target_flank)
                target = str(target_start) + ',' + str(target_length)
                print('#', target, file=sys.stderr)
                print('#', sequence_id)
                wrimer3_result = wrimer3(sequence_id, sequence, target, default_server, max_prod=routine_amplicon_max)
                primers = wrimer3_to_primers(wrimer3_result)
                # if no primers found, retry using level 1 settings
                if len(primers) == 0 and auto_retry == True and auto_retry_level >= 1:
                    print('# Intelligent-retry with relaxed settings, Level 1: relaxing amplicon size to ' +
                          str(retry_amplicon_max) + ' bp', file=sys.stderr)
                    wrimer3_result = wrimer3(sequence_id, sequence, target, default_server, max_prod=retry_amplicon_max)
                    primers = wrimer3_to_primers(wrimer3_result)
                    # if no primers found, retry using level 2 settings
                    if len(primers) == 0 and auto_retry == True and auto_retry_level >= 2:
                        print('# Intelligent-retry with relaxed settings, Level 1: relaxing amplicon size to ' +
                              str(retry_amplicon_max) + ' bp', file=sys.stderr)
                        print('# Intelligent-retry with relaxed settings, Level 2: increasing flanking primer search region to ' +
                              str(retry_total_flank) + ' bp', file=sys.stderr)
                        print('# Intelligent-retry with relaxed settings, Level 2: decreasing obligatory flanking region to ' +
                              str(retry_target_flank) + ' bp', file=sys.stderr)
                        upstream, downstream = get_flanking_regions(translated_region,
                                                        flank=retry_total_flank)
                        upseq = get_masked_sequence(upstream,
                                                    default_server,
                                                    threshold=routine_SNP_threshold,
                                                    repeat_mask=routine_repeat_masking)
                        coreseq = get_sequence(translated_region)
                        downseq = get_masked_sequence(downstream,
                                                      default_server,
                                                      threshold=routine_SNP_threshold,
                                                      repeat_mask=routine_repeat_masking)
                        sequence = upseq + coreseq + downseq
                        target_start = len(upseq) - retry_target_flank
                        target_length = len(coreseq) + 2 * retry_target_flank
                        print('# Total length =', len(sequence), file=sys.stderr)
                        print('#', target_start, retry_target_flank, len(coreseq), retry_target_flank, retry_total_flank - retry_target_flank)
                        target = str(target_start) + ',' + str(target_length)
                        print('#', target, file=sys.stderr)
                        print('#', sequence_id)
                        wrimer3_result = wrimer3(sequence_id, sequence, target, default_server, max_prod=retry_amplicon_max)
                        primers = wrimer3_to_primers(wrimer3_result)
                        # if no primers found, retry using level 3 settings
                        if len(primers) == 0 and auto_retry == True and auto_retry_level >= 3:
                            print('# Intelligent-retry with relaxed settings, Level 1: relaxing amplicon size to ' +
                                  str(retry_amplicon_max) + ' bp', file=sys.stderr)
                            print('# Intelligent-retry with relaxed settings, Level 2: increasing flanking primer search region to ' +
                                  str(retry_total_flank) + ' bp', file=sys.stderr)
                            print('# Intelligent-retry with relaxed settings, Level 2: decreasing obligatory flanking region to ' +
                                  str(retry_target_flank) + ' bp', file=sys.stderr)
                            print('# Intelligent-retry with relaxed settings, Level 3: relaxing SNP masking setting to ' +
                                  str(retry_SNP_threshold), file=sys.stderr)
                            print('# Intelligent-retry with relaxed settings, Level 3: relaxing repeat masking setting to ' +
                                  str(retry_repeat_masking), file=sys.stderr)
                            upstream, downstream = get_flanking_regions(translated_region,
                                                            flank=retry_total_flank)
                            upseq = get_masked_sequence(upstream,
                                                        default_server,
                                                        threshold=retry_SNP_threshold,
                                                        repeat_mask=retry_repeat_masking)
                            coreseq = get_sequence(translated_region)
                            downseq = get_masked_sequence(downstream,
                                                          default_server,
                                                          threshold=retry_SNP_threshold,
                                                          repeat_mask=retry_repeat_masking)
                            sequence = upseq + coreseq + downseq
                            target_start = len(upseq) - retry_target_flank
                            target_length = len(coreseq) + 2 * retry_target_flank
                            print('# Total length =', len(sequence), file=sys.stderr)
                            print('#', target_start, retry_target_flank, len(coreseq), retry_target_flank, retry_total_flank - retry_target_flank)
                            target = str(target_start) + ',' + str(target_length)
                            print('#', target, file=sys.stderr)
                            print('#', sequence_id)
                            wrimer3_result = wrimer3(sequence_id, sequence, target, default_server, max_prod=retry_amplicon_max)
                            primers = wrimer3_to_primers(wrimer3_result)
                            if len(primers) == 0:
                                print('# No primers found after 3 levels of intelligent-retry. Proceeding to next exon...', file=sys.stderr)
                primer_count = 1
                #####
                # Batch submission for improved performance of autoSNPfree
                pair_name, fp, rp = [], [], []
                for p in primers:
                    pair_name.append(sequence_id + '_' + str(primer_count))
                    fp.append(p[0])
                    rp.append(p[1])
                autoSNPfree(pair_name, fp, rp, default_server, batch_mode=True)
                pair_name, fp, rp = None, None, None                
                #####
                for p in primers:
                    pair_name = sequence_id + '_' + str(primer_count)
                    fp = p[0]
                    rp = p[1]
                    print(pair_name, fp, rp, sep='\t')
                    autoSNPfree_result = autoSNPfree(pair_name, fp, rp, default_server)
                    primer_name, primer_score = SNPfree_to_score(autoSNPfree_result)
                    print('#', primer_name, 'SNPfree score:', primer_score)
                    pp.add(sequence_id, primer_name, fp, rp,
                           primer_score, SNPfree_to_score(autoSNPfree_result, get_size=True))
                    primer_count += 1

            if mode == 'MANUAL':
                print('>' + g.name + '_' + translated_region_name)
                print(upseq[0:len(upseq)-routine_target_flank])
                print('[')
                print(upseq[len(upseq)-routine_target_flank:])
                print(coreseq)
                print(downseq[0:routine_target_flank])
                print(']')
                print(downseq[routine_target_flank:])
        else:
            print('Skipping non-translated exon', translated_region_name, file=sys.stderr)

    print('##############################')
    print('# Below is the list of selected primers for the target regions:')

    final_list = pp.getbest()
    for pair in final_list:
        print('#', pair[1], 'chosen for', pair[0], 'with SNPfree score', pair[4])
        print(pair[0], pair[2], pair[3], sep='\t')

    print('##############################')
    print('# Below is the list for ordering:')
    for pair in final_list:
        p_gene, p_region = str(pair[0]).split('_')
        print(pair[2], len(pair[2]), 'No', p_gene, p_region, 'F', pair[5], sep='\t')
        print(pair[3], len(pair[3]), 'No', p_gene, p_region, 'R', pair[5], sep='\t')



if __name__ == '__main__':
    main(sys.argv)



