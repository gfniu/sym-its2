import argparse
import csv
import os
import re
from Bio import SeqIO


class SegmentSeq(object):
    def __init__(self):
        self.seq_name = ''
        self.start_gap = 0
        self.len_18s = 0
        self.before_gap_its1 = 0
        self.len_its1 = 0
        self.before_gap_58s = 0
        self.len_58s = 0
        self.before_gap_its2 = 0
        self.len_its2 = 0
        self.before_gap_28s = 0
        self.len_28s = 0
        self.total = 0


def main():

    parser = argparse.ArgumentParser(description='manual to this script')
    parser.add_argument('--genbank_file_path', type=str, default=None,
                        help="The genbank file path in genbank format.")
    parser.add_argument('--alignment_file_path', type=str, default=None,
                        help="The alignment file path in fasta format.")
    args = parser.parse_args()
    genbank_file_path = args.genbank_file_path
    alignment_file_path = args.alignment_file_path
    # genbank_file_path = '/Users/niuguigui/资料/hmm/gb_its2.gb'
    # alignment_file_path='/Users/niuguigui/资料/hmm/gb_its2_anno_aln.fasta'
    if genbank_file_path == None:
        print(
            'usage: [-h]\n--genbank_file_path is a must')
        return
    print('waiting...')
    aln_records = SeqIO.parse(alignment_file_path, "fasta")
    aln_records = list(aln_records)
    records = SeqIO.parse(genbank_file_path, "genbank")
    list_segment = []
    for record in records:
        segmentSeq = SegmentSeq()
        segmentSeq.seq_name = record.id
        end_18s = 0
        end_its1 = 0
        end_58s = 0
        end_its2 = 0
        end_28s = 0
        for f in record.features:
            if f.type != 'source':
                liststr = list(f.qualifiers.values())[0][0]
                if liststr.find('contains') > -1:
                    continue
                if liststr.find('18S') > -1:
                    end_18s = int(f.location.end)
                if liststr.find('ITS1') > -1 or liststr.find('spacer 1') > -1:
                    end_its1 = int(f.location.end)
                if liststr.find('5.8S') > -1:
                    end_58s = int(f.location.end)
                if liststr.find('ITS2') > -1 or liststr.find('spacer 2') > -1:
                    end_its2 = int(f.location.end)
                if liststr.find('28S') > -1:
                    end_28s = int(f.location.end)
        aln_filter = list(
            filter(lambda x: x.id == segmentSeq.seq_name, aln_records))
        if len(aln_filter) > 0:
            aln_char_list = list(aln_filter[0])
            i = 0
            position = 0
            base_start = 0
            segment_end = False
            gap = 0
            for aln_char in aln_char_list:
                if segment_end and aln_char == '-':
                    gap = gap+1
                if aln_char != '-':
                    segment_end = False
                    i = i+1
                    if base_start == 0:
                        base_start = position-1
                    if i > 0:
                        if i == end_18s:
                            segmentSeq.len_18s = position-base_start
                            segment_end = True
                            gap = 0
                        if i == end_its1:
                            segmentSeq.before_gap_its1 = gap
                            segmentSeq.len_its1 = position-segmentSeq.len_18s - \
                                base_start-segmentSeq.before_gap_its1
                            segment_end = True
                            gap = 0
                        if i == end_58s:
                            segmentSeq.before_gap_58s = gap
                            segmentSeq.len_58s = position-segmentSeq.len_its1-segmentSeq.len_18s - \
                                base_start-segmentSeq.before_gap_58s-segmentSeq.before_gap_its1
                            segment_end = True
                            gap = 0
                        if i == end_its2:
                            segmentSeq.before_gap_its2 = gap
                            segmentSeq.len_its2 = position-segmentSeq.len_58s - \
                                segmentSeq.len_its1-segmentSeq.len_18s-base_start-segmentSeq.before_gap_its2 - \
                                segmentSeq.before_gap_58s-segmentSeq.before_gap_its1
                            segment_end = True
                            gap = 0
                        if i == end_28s:
                            segmentSeq.before_gap_28s = gap
                            segmentSeq.len_28s = position-segmentSeq.len_its2-segmentSeq.len_58s - \
                                segmentSeq.len_its1-segmentSeq.len_18s-base_start-segmentSeq.before_gap_28s - \
                                segmentSeq.before_gap_its2-segmentSeq.before_gap_58s-segmentSeq.before_gap_its1
                            segment_end = True
                            gap = 0
                position = position+1
            segmentSeq.start_gap = base_start+1
            segmentSeq.total = sum([segmentSeq.start_gap, segmentSeq.len_18s, segmentSeq.before_gap_its1, segmentSeq.len_its1, segmentSeq.before_gap_58s,
                                    segmentSeq.len_58s, segmentSeq.before_gap_its2, segmentSeq.len_its2, segmentSeq.before_gap_28s, segmentSeq.len_28s])
        if segmentSeq.len_its2 != 0:
            list_segment.append(segmentSeq)
    with open("gb_its2_anno.csv", "w", newline="") as datacsv:
        csvwriter = csv.writer(datacsv, dialect=("excel"))
        csvwriter.writerow(["seq_name", "start_gap", "len_18s", "before_gap_its1",
                            "len_its1", "before_gap_5.8s", "len_5.8s", "before_gap_its2", "len_its2", "before_gap_28s", "len_28s", "total"])
        for segment in list_segment:
            csvwriter.writerow([segment.seq_name, segment.start_gap, segment.len_18s, segment.before_gap_its1,
                                segment.len_its1, segment.before_gap_58s, segment.len_58s, segment.before_gap_its2, segment.len_its2,
                                segment.before_gap_28s, segment.len_28s, segment.total])
    print('Success, results are saved in gb_its2_anno.csv.\nend.')


if __name__ == "__main__":
    main()