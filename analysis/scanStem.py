import argparse
import csv
import os
from Bio import SeqIO


class ParsingHmmResult(object):
    def __init__(self):
        self.seq_name = ''
        self.sequence = ''
        self.motif_58s = ''
        self.its2 = ''
        self.motif_28s = ''
        self.its2_len = -1
        self.its2_start = -1
        self.its2_end = -1
        self.remark = ''
        self.hits = Hit()

    def hmmQuery(self, query_context):
        query_lines = query_context.split('\n')
        i = -1
        start_hit = False
        list_hit = []
        for query_line in query_lines:
            i = i+1
            query_line = query_line.lstrip()
            if query_line.startswith('>>'):
                start_hit = True
            if query_line.endswith(' PP'):
                start_hit = False
            if start_hit:
                if query_line.startswith('Internal pipeline statistics summary'):
                    break
                if query_line.startswith('Alignment:'):
                    hit = Hit()
                    hit.score = float(query_lines[i+1].lstrip().split(' ')[1])
                    if i >= 2:
                        line_reduce2_info = ' '.join(
                            query_lines[i-2].strip().split()).split(' ')
                        hit.reliable = True if line_reduce2_info[0] == "!" else False
                    line2_info = ' '.join(
                        query_lines[i+2].strip().split()).split(' ')
                    hit.type = line2_info[0]
                    hit.hmm_start = int(line2_info[1])
                    hit.hmm_end = int(line2_info[3])
                    hit.hmm = line2_info[2]
                    line3_info = query_lines[i+3].strip()
                    hit.align = line3_info
                    line4_info = ' '.join(
                        query_lines[i+4].strip().split()).split(' ')
                    self.seq_name = line4_info[0]
                    hit.motif = line4_info[2]
                    hit.motif_start = line4_info[1]
                    hit.motif_end = line4_info[3]
                    hit.hmm_start = int(line2_info[1])
                    hit.hmm_end = int(line2_info[3])
                    hit.hmm = line2_info[2]
                    if (hit.type == "5.8s" and hit.hmm_end < 24) or (hit.type == "28s" and hit.hmm_start > 1):
                        hit.reliable = False
                    list_hit.append(hit)
        self.hits = list_hit


class Hit():
    def __init__(self):
        self.score = -1
        self.type = ''
        self.hmm_start = -1
        self.hmm_end = -1
        self.hmm = ''
        self.motif_start = -1
        self.motif_end = -1
        self.motif = ''
        self.align = ''
        self.reliable = False


def main():
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--seqfile_path', type=str, default=None,
                        help="The sequence file path in fasta format.")
    parser.add_argument('--hmmfile_path', type=str,
                        default=None, help="hmm database file path.")
    parser.add_argument('--hmm_result_path', type=str,
                        default=None, help="hmm run results save path.")
    parser.add_argument('--its2_minlength', type=int,
                        default=140, help="ITS2 Mininum size; default 140")
    parser.add_argument('--e_value', type=float, default=1,
                        help="e-value; default 1")
    args = parser.parse_args()
    its2_minlength = args.its2_minlength
    e_value = args.e_value
    seqfile_path = args.seqfile_path
    hmmfile_path = args.hmmfile_path
    hmm_result_path = args.hmm_result_path
    if seqfile_path == None or hmmfile_path == None or hmm_result_path == None:
        print(
            'usage: [-h]\n--seqfile_path, --hmmfile_path, --hmm_result_path is a must')
        return
    print('begin run hmmer')
    process = os.popen('nhmmscan -E '+str(e_value)+' --watson ' +
                       hmmfile_path + ' ' + seqfile_path + ' > ' + hmm_result_path)
    process.close()
    print('begin parsing hmmer result')
    result_handle = open(hmm_result_path)
    result_context = result_handle.read()
    query_records = result_context.split('Query:')
    list_its2 = []
    seq_records = SeqIO.parse(seqfile_path, "fasta")
    seq_records = list(seq_records)
    error_58s = '5.8s motif possible errors;'
    error_28s = '28s motif possible errors;'
    for query_record in query_records:
        parsingHmm = ParsingHmmResult()
        if query_record.startswith('# nhmmscan'):
            continue
        parsingHmm.hmmQuery(query_record.lstrip())
        hits = parsingHmm.hits
        motif_28s_list = list(filter(lambda x: x.type == '28s', hits))
        motif_58s_list = list(filter(lambda x: x.type == '5.8s', hits))
        score_temp = -1
        its2_len_temp = -1
        if len(motif_28s_list) == 1 and len(motif_58s_list) == 0:
            parsingHmm.motif_28s = motif_28s_list[0].motif
            parsingHmm.its2_end = int(motif_28s_list[0].motif_start) - 1
            parsingHmm.remark += '' if motif_28s_list[0].reliable else error_28s
        if len(motif_58s_list) == 1 and len(motif_28s_list) == 0:
            parsingHmm.motif_58s = motif_58s_list[0].motif
            parsingHmm.its2_start = int(motif_58s_list[0].motif_end) + 1
            parsingHmm.remark = '' if motif_58s_list[0].reliable else error_58s
        if len(motif_28s_list) == 1 and len(motif_58s_list) == 1:
            parsingHmm.motif_58s = motif_58s_list[0].motif
            parsingHmm.its2_start = int(motif_58s_list[0].motif_end) + 1
            parsingHmm.motif_28s = motif_28s_list[0].motif
            parsingHmm.its2_end = int(motif_28s_list[0].motif_start) - 1
            parsingHmm.its2_len = parsingHmm.its2_end-parsingHmm.its2_start+1
            parsingHmm.remark = '' if motif_58s_list[0].reliable else error_58s
            parsingHmm.remark += '' if motif_28s_list[0].reliable else error_28s
        if len(motif_28s_list) > 1 and len(motif_58s_list) > 1:
            for motif_28s in motif_28s_list:
                for motif_58s in motif_58s_list:
                    its2_len_temp = int(motif_28s.motif_start) - \
                        int(motif_58s.motif_end)-1
                    if its2_len_temp >= its2_minlength and ((its2_len_temp < parsingHmm.its2_len or parsingHmm.its2_len == -1) or
                                                            (its2_len_temp == parsingHmm.its2_len and
                                                                score_temp < int(motif_28s.score) + int(motif_58s.score))):
                        score_temp = int(motif_28s.score) + \
                            int(motif_58s.score)
                        parsingHmm.its2_len = its2_len_temp
                        if motif_58s.reliable or parsingHmm.its2_start == -1 or \
                                (motif_58s.reliable == False and parsingHmm.remark.find('5.8s') > -1):
                            parsingHmm.its2_start = int(
                                motif_58s.motif_end) + 1
                            parsingHmm.motif_58s = motif_58s.motif
                        if motif_28s.reliable or parsingHmm.its2_end == -1 or \
                                (motif_28s.reliable == False and parsingHmm.remark.find('28s') > -1):
                            parsingHmm.its2_end = int(
                                motif_28s.motif_start) - 1
                            parsingHmm.motif_28s = motif_28s.motif
                        parsingHmm.remark = '' if motif_58s.reliable else error_58s
                        parsingHmm.remark += '' if motif_28s.reliable else error_28s
        if len(motif_28s_list) > 1 and len(motif_58s_list) == 1:
            parsingHmm.its2_start = int(motif_58s_list[0].motif_end)+1
            parsingHmm.motif_58s = motif_58s_list[0].motif
            for motif_28s in motif_28s_list:
                its2_len_temp = int(motif_28s.motif_start) - \
                    int(motif_58s_list[0].motif_end)-1
                if its2_len_temp >= its2_minlength and ((its2_len_temp < parsingHmm.its2_len or parsingHmm.its2_len == -1) or
                                                        (its2_len_temp == parsingHmm.its2_len and score_temp < int(motif_28s.score))) and \
                        (motif_28s.reliable or parsingHmm.its2_end == -1 or
                            (motif_28s.reliable == False and parsingHmm.remark.find('28s') > -1)):
                    score_temp = int(motif_28s.score)
                    parsingHmm.its2_len = its2_len_temp
                    parsingHmm.its2_end = int(motif_28s.motif_start) - 1
                    parsingHmm.motif_28s = motif_28s.motif
                    parsingHmm.remark = '' if motif_58s_list[0].reliable else error_58s
                    parsingHmm.remark += '' if motif_28s.reliable else error_28s
        if len(motif_28s_list) > 1 and len(motif_58s_list) == 0:
            for motif_28s in motif_28s_list:
                if score_temp == -1 or int(motif_28s.score) > score_temp and \
                    (motif_28s.reliable or parsingHmm.its2_end == -1 or
                        (motif_28s.reliable == False and parsingHmm.remark.find('28s') > -1)):
                    score_temp = int(motif_28s.score)
                    parsingHmm.its2_end = int(motif_28s_list[0].motif_start)-1
                    parsingHmm.motif_28s = motif_28s.motif
                    parsingHmm.remark = '' if motif_28s.reliable else error_28s
        if len(motif_58s_list) > 1 and len(motif_28s_list) == 1:
            parsingHmm.its2_end = int(motif_28s_list[0].motif_start)-1
            parsingHmm.motif_28s = motif_28s_list[0].motif
            for motif_58s in motif_58s_list:
                its2_len_temp = int(
                    motif_28s_list[0].motif_start) - int(motif_58s.motif_end)-1
                if its2_len_temp >= its2_minlength and \
                    ((its2_len_temp < parsingHmm.its2_len or parsingHmm.its2_len == -1) or
                        (its2_len_temp == parsingHmm.its2_len and score_temp < int(motif_58s.score))) and \
                    (motif_58s.reliable or parsingHmm.its2_start == -1 or
                     (motif_58s.reliable == False and parsingHmm.remark.find('5.8') > -1)):
                    score_temp = int(motif_58s.score)
                    parsingHmm.its2_len = its2_len_temp
                    parsingHmm.its2_start = int(motif_58s.motif_end) + 1
                    parsingHmm.motif_58s = motif_58s.motif
                    parsingHmm.remark = '' if motif_58s.reliable else error_58s
                    parsingHmm.remark += '' if motif_28s_list[0].reliable else error_28s
        if len(motif_58s_list) > 1 and len(motif_28s_list) == 0:
            for motif_58s in motif_58s_list:
                if (score_temp == -1 or int(motif_58s.score) > score_temp) and \
                    (motif_58s.reliable or parsingHmm.its2_start == -1 or
                        (motif_58s.reliable == False and parsingHmm.remark.find('5.8') > -1)):
                    score_temp = int(motif_58s.score)
                    parsingHmm.its2_end = int(motif_58s_list[0].motif_start)-1
                    parsingHmm.motif_58s = motif_58s.motif
                    parsingHmm.remark = '' if motif_58s.reliable else error_58s
        seq_filter = list(
            filter(lambda x: x.id == parsingHmm.seq_name, seq_records))
        if len(seq_filter)>0:                
            parsingHmm.sequence = str(seq_filter[0].seq)
        if len(motif_58s_list) > 0 and len(motif_28s_list) > 0:
            parsingHmm.its2 = parsingHmm.sequence[parsingHmm.its2_start -
                                                  1:parsingHmm.its2_end]
        if parsingHmm.motif_28s == '':
            parsingHmm.remark += 'No 28s motif was found.'
        if parsingHmm.motif_58s == '':
            parsingHmm.remark += 'No 5.8s motif was found.'
        list_its2.append(parsingHmm)
    with open("scanStemResult.csv", "w", newline="") as datacsv:
        csvwriter = csv.writer(datacsv, dialect=("excel"))
        csvwriter.writerow(["seq_name", "sequence", "5.8s_motif", "its2",
                            "28s_motif", "its2_len", "its2_start", "its2_end", "remark"])
        for its2 in list_its2:
            csvwriter.writerow([its2.seq_name, its2.sequence, its2.motif_58s, its2.its2,
                                its2.motif_28s, its2.its2_len, its2.its2_start, its2.its2_end, its2.remark])
    print('Success, results are saved in scanStemResult.csv.\nend.')


if __name__ == "__main__":
    main()
