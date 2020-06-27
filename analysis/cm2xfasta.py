import argparse
import csv
import os
import re
from Bio import SeqIO


def buildDotStr(dot_num):
    str_list = []
    num = 0
    while num < dot_num:
        str_list.append('.')
        num = num+1
    return "".join(str_list)


def checkPair(s):
    symbols = {'}': '{', ']': '[', ')': '('}
    symbols_l, symbols_r = symbols.values(), symbols.keys()
    arr = []
    for c in s:
        if c in symbols_l:
            arr.append(c)
        elif c in symbols_r:
            if arr and arr[-1] == symbols[c]:
                arr.pop()
            else:
                return False
    return not arr


class ParsingCmResult(object):
    def __init__(self):
        self.seq_name = ''
        self.sequence = ''
        self.structure = ''
        self.part_seq = ''
        self.part_struct = ''
        self.part_start = -1
        self.part_end = -1
        self.remark = ''

    def CmQuery(self, query_context):
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
                if query_line.startswith('Internal CM pipeline statistics summary'):
                    break
                if query_line.startswith('>>'):
                    line3_info = ' '.join(
                        query_lines[i+3].strip().split()).split(' ')
                    if line3_info[1] != '!':
                        break
                    line6_info = ' '.join(
                        query_lines[i+6].strip().split()).split(' ')
                    self.part_struct = line6_info[0]
                    line9_info = ' '.join(
                        query_lines[i+9].strip().split()).split(' ')
                    self.seq_name = line9_info[0]
                    i = 2
                    while i < len(line9_info)-1:
                        if(i < len(line9_info)-1) and i > 2:
                            self.part_seq += ' '
                        self.part_seq += line9_info[i]
                        i += 1
                    self.part_start = int(line9_info[1])
                    self.part_end = int(line9_info[len(line9_info)-1])
                    break


class XfastaFormat():
    def __init__(self):
        self.name = ''
        self.sequence = ''
        self.structure = ''


def main():
    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--seqfile_path', type=str, default=None,
                        help="The sequence file path in fasta format.")
    parser.add_argument('--cm_result_path', type=str,
                        default=None, help="The infernal result file path.")
    args = parser.parse_args()
    cm_result_path = args.cm_result_path
    seqfile_path = args.seqfile_path
    if seqfile_path == None or cm_result_path == None:
        print('usage: [-h]\n--seqfile_path, --cm_result_path is a must')
        return
    result_handle = open(cm_result_path)
    result_context = result_handle.read()
    query_records = result_context.split('Query:')
    list_its2 = []
    seq_records = SeqIO.parse(seqfile_path, "fasta")
    seq_records = list(seq_records)
    xfasta = []
    for query_record in query_records:
        parsingCm = ParsingCmResult()
        if query_record.startswith('# cmscan'):
            continue
        parsingCm.CmQuery(query_record.lstrip())
        seq_filter = list(
            filter(lambda x: x.id == parsingCm.seq_name, seq_records))
        if len(seq_filter) <= 0:
            print(parsingCm.seq_name)
            continue
        seq_filter = seq_filter[0]
        sequence = str(seq_filter.seq)
        parsingCm.sequence = sequence[0:parsingCm.part_start-1] + \
            parsingCm.part_seq+sequence[parsingCm.part_end:]

        parsingCm.structure = buildDotStr(
            parsingCm.part_start-1)+parsingCm.part_struct+buildDotStr(len(sequence)-parsingCm.part_end)
        parsingCm.structure = re.sub('[^<>]', '.', parsingCm.structure)
        parsingCm.structure = parsingCm.structure.replace(
            '<', '(').replace('>', ')')
        # omitting_seqs = re.findall("(?<=\*\[).*?(?=\]\*)", parsingCm.sequence)
        omitting_seqs = re.findall("\*\[.*?\]\*", parsingCm.sequence)
        for omitting_seq in omitting_seqs:
            omitting_num = int(omitting_seq.replace(
                '*[', '').replace(']*', ''))
            sequence_temp = parsingCm.sequence.replace('-', '')
            index_temp = sequence_temp.find(omitting_seq)
            index = parsingCm.sequence.find(omitting_seq)
            parsingCm.sequence = parsingCm.sequence[0:index]+sequence[index_temp:index_temp +
                                                                      omitting_num]+parsingCm.sequence[index+len(omitting_seq):]
            parsingCm.structure = parsingCm.structure[0:index]+buildDotStr(
                omitting_num)+parsingCm.structure[index+len(omitting_seq):]
        char_finds = re.finditer('-', parsingCm.sequence)
        for char_f in char_finds:
            char_position = char_f.start()
            if parsingCm.structure[char_position:char_position+1] == '.':
                parsingCm.sequence = parsingCm.sequence[:char_position] + \
                    ','+parsingCm.sequence[char_position+1:]
                parsingCm.structure = parsingCm.structure[:char_position] + \
                    ','+parsingCm.structure[char_position+1:]
        parsingCm.sequence = parsingCm.sequence.replace(',', '').upper()
        parsingCm.structure = parsingCm.structure.replace(',', '')
        parsingCm.structure = parsingCm.structure.replace('()', '..')
        parsingCm.structure = parsingCm.structure.replace('(.)', '...')
        parsingCm.structure = parsingCm.structure.replace('(..)', '....')
        parsingCm.sequence = parsingCm.sequence.replace('T', 'U')
        if parsingCm.sequence.find('<[0]*') > -1:
            parsingCm.sequence = parsingCm.sequence.replace('<[0]*', '')
            parsingCm.structure = parsingCm.structure[5:]
        if parsingCm.sequence.find('<[ 0]*') > -1:
            parsingCm.sequence = parsingCm.sequence.replace('<[ 0]*', '')
            parsingCm.structure = parsingCm.structure[6:]
        if parsingCm.sequence.find('*[0]>') > -1:
            parsingCm.sequence = parsingCm.sequence.replace('*[0]>', '')
            parsingCm.structure = parsingCm.structure[:-5]
        if parsingCm.sequence.find('*[ 0]>') > -1:
            parsingCm.sequence = parsingCm.sequence.replace('*[ 0]>', '')
            parsingCm.structure = parsingCm.structure[:-6]

        xfastaFormat = XfastaFormat()
        xfastaFormat.name = parsingCm.seq_name
        xfastaFormat.sequence = parsingCm.sequence
        xfastaFormat.structure = parsingCm.structure
        xfasta.append(xfastaFormat)
    filename = "core set secondary structures.xfasta"
    out_handle = open(filename, "w")
    for x in xfasta:
        if checkPair(x.structure) == False:
            print(x.name+' structural error, please note')
        out_handle.write('>'+x.name+'\n')
        out_handle.write(x.sequence+'\n')
        out_handle.write(x.structure+'\n')
    out_handle.close()
    print('end.')


if __name__ == "__main__":
    main()