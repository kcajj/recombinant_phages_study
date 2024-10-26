import pysam


def create_coordinate_conversion_map(bam_file):
    map = {'gaps': []} #one key and value with gaps of the read on the reference and the rest with the conversion map
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not (read.is_secondary):
                alignment = read.get_aligned_pairs()
                for read_base, ref_base in alignment:
                    if read_base == None:
                        map['gaps'].append(ref_base)
                    if read_base in map.keys(): #we are dealing with a summplementary alignment
                        if map[read_base] == None:
                            map[read_base] = ref_base #it could be that the supplementary alignment brings more info
                        else:
                            continue #but we don't want to overwrite the primary info
                    else:
                        map[read_base] = ref_base
    return map
