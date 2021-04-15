import sys
from pysam import VariantFile

def check(with_gt = True):
    vcf_path = sys.argv[1]
    vcf_path_1 = sys.argv[2]
    vcf_path_2 = sys.argv[3]

    vcf = VariantFile(vcf_path, 'r', drop_samples = True)
    callable = set()
    for record in vcf:
        pos = record.pos
        callable.add(pos)

    vcf_1 = VariantFile(vcf_path_1, 'r', drop_samples = False)
    #ref_real = set()
    alt_real = set()
    for record in vcf_1:
        pos = record.pos
        if pos < 51 or pos > 29852:
            continue
        gt = record.samples[0]['GT'][0]
        #if gt == 0:
        #    ref_real.add(pos)
        #else:
        #if gt > 0:
        if pos in callable:
            alt_real.add(pos)

    alt_called = set()
    try:
        vcf_2 = VariantFile(vcf_path_2, 'r', drop_samples = False)
        for record in vcf_2:
            pos = record.pos
            if pos < 51 or pos > 29852:
                continue
            if pos not in callable:
                continue
            if with_gt:
                gt = record.samples[0]['GT'][0]
                if gt > 0:
                    alt_called.add(pos)
            else:
                if record.qual > 1000:
                    alt_called.add(pos)
    except:
        pass
    print("Truth", alt_real)
    print("Called", alt_called)
    #refTPs = len(ref_called & ref_real)
    #refFPs = len(ref_called - ref_real)
    #refFNs = len(ref_real - ref_called)
    #refP = round(refTPs/(refTPs+refFPs), 3) if refTPs+refFPs != 0 else 0
    #refR = round(refTPs/(refTPs+refFNs), 3) if refTPs+refFNs != 0 else 0

    altTPs = len(alt_called & alt_real)
    altFPs = len(alt_called - alt_real)
    altFNs = len(alt_real - alt_called)
    altP = round(altTPs/(altTPs+altFPs), 3) if altTPs+altFPs != 0 else 0
    altR = round(altTPs/(altTPs+altFNs), 3) if altTPs+altFNs != 0 else 0

    print("TPs", "FPs", "FNs", "P", "R", sep = ',')
    print(#refTPs, refFPs, refFNs, refP, refR,
          altTPs, altFPs, altFNs, altP, altR, sep=',')


if __name__ == "__main__":
    mode = sys.argv.pop(1)
    if mode == "malva" or mode == "bcftools":
        check()
    elif mode == "lofreq":
        check(False)
