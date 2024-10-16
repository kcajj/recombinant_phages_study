import pysam

def get_hyb_ref_map(bam_file):
    map_hyb_ref={}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not(read.is_secondary):
                alignment=read.get_aligned_pairs()
                for hyb,ref in alignment:
                    if hyb==None:
                        continue
                    if hyb in map_hyb_ref.keys():
                        if map_hyb_ref[hyb]==None:
                            map_hyb_ref[hyb]=ref
                        else:
                            continue
                    else:
                        map_hyb_ref[hyb]=ref
    return map_hyb_ref