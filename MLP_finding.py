# #######################
# 31/10/2022, wangxm
#
# 1）genotyping_format_transfer(). This function transfer the format of the iuput haplotype file
# 2) get_weighted_matrix(). This function generate the matrix which define the weight coefficient to modify the
# penalty score in haplotype comparisons for each haplotype window based on the haplotype number in the compared window.
# 3) hap_smooth (). This function smooths the mismathed windows basen the penalty socres and the flanking matched windows
# 4) contri_matrix_counter(). This function find the minimum landrace path based on the compared haplotypes and write the
# related results.

# #######################################
import argparse
import gzip
import math
import os
import re
import numpy as np

def getchrlength(chr_lengthFile, pangenome, windowsize):
    lengs = {}

    for li in open(chr_lengthFile, "r+"):
        refer, chr, start, end = li.strip().split("\t")
        chr = chr.split("_")[0]

        if refer == pangenome and chr != "chrUn":
            leng = int(math.ceil(int(end)/windowsize))
            lengs[chr] = leng

    return lengs


def genotyping_format_transfer(querysamples, targetsamples, genotypingFile, outputpath):
    #fullfilepath = path + genotypingFile
    #infile = open(genotypingFile, "r+")

    ourputFile = open(os.path.join(outputpath, "standerd_genotype.txt"), "w")

    tempdata = {}
    for position in gzip.open(genotypingFile, "rt"):
    #for position in open(genotypingFile, "r"):
        if position.startswith("chr\t"):
            continue

        (chr, start, end, query, dmp_max) = position.strip().split("\t")
        chr = chr.split("_")[0]
        start = int(start)

        if query in querysamples or query in targetsamples:

            if chr not in tempdata.keys():
                tempdata[chr] = {}

            if query not in tempdata[chr].keys():
                tempdata[chr][query] = {}

            if start not in tempdata[chr][query].keys():
                tempdata[chr][query][start] = dmp_max

    writtensamples = []
    for chr in sorted(tempdata.keys()):
        for query in sorted(tempdata[chr].keys()):
            flag = 1
            #if query in querysamples or query in targetsamples:
            ourputFile.write(chr + "\t" + query + "\t")
            for i in sorted(tempdata[chr][query].keys()):

                if int(i) < flag:   ### making sure the haolotype were write in order
                    print(chr + " " + query + " in" + start + " has a error in genotyping transform")
                else:
                    ourputFile.write(tempdata[chr][query][i] + ";")
                flag = int(i)

            ourputFile.write("\n")
            writtensamples.append(query)

    for sample in querysamples:
        if sample not in writtensamples:
            print(sample + " have no genotyping data")

    for sample in targetsamples:
        if sample not in writtensamples:
            print(sample + " have no genotyping data")

    ourputFile.close()

def readsamples(sampleFile):
    samples = []

    for li in open(sampleFile, "r+"):
        sample = li.strip().split()[0]
        samples.append(sample)

    return samples

def get_weighted_matrix(hap_interval, outputpath):
    genotypingFile = os.path.join(outputpath, "standerd_genotype.txt")
    hapnumwriter = open(os.path.join(outputpath, "hap_num_window.txt"), "w")
    hapnumwriter.write("\t".join(["chr", "posi", "hapnum"]) + "\n")

    scorewriter = open(os.path.join(outputpath, "weighted_matrix.txt"), "w")

    genotypes = {}
    for li in open(genotypingFile, "r+"):
        ichr, iquery, ihaps = li.strip().split("\t")

        if ichr not in genotypes.keys():
            genotypes[ichr] = {}

        samplegenotype = ihaps[:-1].split(";")
        location = 0
        for windowhap in samplegenotype:
            location += 1

            if location not in genotypes[ichr].keys():
                genotypes[ichr][location] = []

            if windowhap not in genotypes[ichr][location]:
                genotypes[ichr][location].append(windowhap)

    ######recording the haplotype number for each window
    hapnums = []
    for ichr in genotypes.keys():
        for location in genotypes[ichr].keys():
            hapnum = len(genotypes[ichr][location])
            hapnumwriter.write("\t".join([ichr, str(location), str(hapnum)]) + "\n")
            hapnums.append(hapnum)

    hapnums = sorted(hapnums)

    a = np.array(hapnums)
    thresholds = [np.percentile(a, x) for x in range(100 // hap_interval, 100, 100 // hap_interval)]

    ####define the weighted coefficient for each window
    weights = {}
    for ichr in genotypes.keys():
        weights[ichr] = {}
        for location in genotypes[ichr].keys():
            hapnum = len(genotypes[ichr][location])

            if hapnum <= thresholds[0]:
                weights[ichr][location] = 1
            elif hapnum > thresholds[-1]:
                weights[ichr][location] = 1 - len(thresholds) * 0.1
            else:

                for i in range(0, (len(thresholds) - 1)):
                    if hapnum > thresholds[i] and hapnum <= thresholds[i + 1]:
                        weights[ichr][location] = 1 - (i + 1) * 0.1
                        break

    ###writing the results
    for chr in weights.keys():
        for location in weights[chr]:
            scorewriter.write(chr + "\t" + str(location) + "\t" + str(weights[chr][location]) + "\n")

    scorewriter.close()
    hapnumwriter.close()

def read_weighted_matrix(path):
    weighted_matrix = {}

    for li in open(os.path.join(path, "weighted_matrix.txt"), "r+"):
        chr, startposi, coefficient = li.strip().split("\t")

        if chr not in weighted_matrix.keys():
            weighted_matrix[chr] = {}

        weighted_matrix[chr][int(startposi) -1] = float(coefficient)


    return weighted_matrix

def read_Trans_Genotype(querysample, targetsamples, hapnum_thresh, path):
    genotypingFile = os.path.join(path, "standerd_genotype.txt")
    hap_num_File = os.path.join(path, "hap_num_window.txt")

    hap_nums = {}
    for li in open(hap_num_File, "r+"):
        ichr, posi, num = li.strip().split("\t")

        if ichr == "chr":
            continue

        if ichr not in hap_nums.keys():
            hap_nums[ichr] = {}

        hap_nums[ichr][posi] = int(num)


    genotypes = {}
    queryGenotype = {}
    for li in open(genotypingFile, "r+"):
        ichr, iquery, ihaps = li.strip().split("\t")

        if iquery == querysample:
            queryGenotype[ichr] = ihaps[:-1].split(";")

        if ichr not in genotypes.keys():
            genotypes[ichr] = {}
        if iquery in targetsamples:
            genotypes[ichr][iquery] = ihaps[:-1].split(";")

    transferedGeno = {}
    for ichr in genotypes.keys():
        if ichr not in transferedGeno.keys():
            transferedGeno[ichr] = {}
        for targetsample in genotypes[ichr].keys():
            transferedGeno[ichr][targetsample] = ""

            if len(queryGenotype[ichr]) != len(genotypes[ichr][targetsample]):
                print(targetsample + "at " + ichr + " has different hap number" + "\n")
            else:
                if querysample == targetsample:
                    for i in range(0, len(queryGenotype[ichr])):
                        transferedGeno[ichr][targetsample] += "0"

                else:
                    for i in range(0, len(queryGenotype[ichr])):
                        if queryGenotype[ichr][i] == genotypes[ichr][targetsample][i]:
                            transferedGeno[ichr][targetsample] += "1"
                        elif hap_nums[ichr][str(i + 1)] > hapnum_thresh:
                            transferedGeno[ichr][targetsample] += "1"
                        else:
                            transferedGeno[ichr][targetsample] += "0"

    return (transferedGeno, queryGenotype, genotypes)

def hap_smooth(type, transferredGenotype, weighted_matrix, max_smooth_length, mismatch_penalty_coefficient,
               gap_intro_penalty_coefficient, match_coefficient, revise_coefficient, chrslength):
    if type == "0to1":
        mismath_pattern = r"10{1," + str(max_smooth_length) + "}1"
    elif type == "1to0":
        mismath_pattern = r"01{1," + str(max_smooth_length) + "}0"

    for chr in transferredGenotype.keys():
        for sample in transferredGenotype[chr].keys():

            ### find the mismathed win
            for substr in re.finditer(mismath_pattern, transferredGenotype[chr][sample]):
                startposi = substr.start() + 1
                mismath_len = len(substr.group()) - 2

                penaltyscore = 0

                if mismath_len == 1:  ### single mismatched window
                    penaltyscore = weighted_matrix[chr][startposi] * mismatch_penalty_coefficient
                elif mismath_len > 1:
                    penaltyscore = weighted_matrix[chr][startposi] * gap_intro_penalty_coefficient
                    for l in range(0, mismath_len):
                        posi = startposi + l
                        penaltyscore += weighted_matrix[chr][posi] * mismatch_penalty_coefficient

                flankingscore = 0
                leftposi = startposi - 1
                rightposi = startposi + mismath_len

                revised_threshold = penaltyscore * revise_coefficient * 2

                while flankingscore < revised_threshold:

                    if leftposi < 0 or rightposi >= chrslength[chr]:
                        break

                    if type == "0to1":
                        if transferredGenotype[chr][sample][leftposi] == "1" and transferredGenotype[chr][sample][
                            rightposi] == "1":
                            leftscore = weighted_matrix[chr][leftposi] * match_coefficient
                            rightscore = weighted_matrix[chr][rightposi] * match_coefficient
                            flankingscore = flankingscore + leftscore + rightscore

                        elif transferredGenotype[chr][sample][leftposi] == "0" or transferredGenotype[chr][sample][
                            rightposi] == "0":
                            break

                    if type == "1to0":
                        if transferredGenotype[chr][sample][leftposi] == "0" and transferredGenotype[chr][sample][
                            rightposi] == "0":
                            leftscore = weighted_matrix[chr][leftposi] * match_coefficient
                            rightscore = weighted_matrix[chr][rightposi] * match_coefficient
                            flankingscore = flankingscore + leftscore + rightscore

                        elif transferredGenotype[chr][sample][leftposi] == "1" or transferredGenotype[chr][sample][
                            rightposi] == "1":
                            break

                    leftposi -= 1
                    rightposi += 1


                if flankingscore >= revised_threshold:
                    seqs = list(transferredGenotype[chr][sample])
                    if mismath_len == 1:
                        if type == "0to1":
                            seqs[startposi] = "1"
                        elif type == "1to0":
                            seqs[startposi] = "0"

                        #transferredGenotype[chr][sample][startposi] = "1"
                    elif mismath_len > 1:
                        for l in range(0, mismath_len):
                            reposi = startposi + l
                            if type == "0to1":
                                seqs[reposi] = "1"
                            elif type == "1to0":
                                seqs[reposi] = "0"

                            #transferredGenotype[chr][sample][reposi] = "1"

                    transferredGenotype[chr][sample] = "".join(seqs)

def slidingwidowCheck(endposi, matchedblock, slidingwindowsize, maxmismatchednuminslidingwidow):
    newendposi = endposi
    if len(matchedblock) <= slidingwindowsize:
        if matchedblock.count("0") > maxmismatchednuminslidingwidow or matchedblock.count("0") > (len(matchedblock) / 2):
            newendposi = matchedblock.find("0")
    else: ## using sliding window to check the density of the mismatched windows
        for i in range(0, (len(matchedblock) - slidingwindowsize)):
            checkblock = matchedblock[i:(i+slidingwindowsize)]
            if checkblock.count("0") > maxmismatchednuminslidingwidow:
                newendposi = checkblock.find("0") + i
                break

    return newendposi

def extendblock(sample, chr, targetSamplesGenotype, startposi, linkerlength, slidingwindowsize, maxmismatchednuminslidingwidow):

    genotyping = targetSamplesGenotype[chr][sample][(startposi-1):]
    mismatchnum = 0
    endposi = 0

    if not genotyping.startswith("0"):

        linker = "0" * (linkerlength + 1)
        firstposi = genotyping.find(linker)

        if firstposi == -1:
            endposi = len(genotyping)
        else:
            endposi = firstposi

        matchedblock = genotyping[0:endposi]

        endposi = slidingwidowCheck(endposi, matchedblock, slidingwindowsize, maxmismatchednuminslidingwidow)

        mismatchnum = genotyping[0:endposi].count("0")

    if endposi == 0:
        endposi = startposi + endposi
        mismatchnum = 1
    else:
        endposi = startposi + endposi - 1


    return endposi, mismatchnum

def MTPfinding(targetSamplesGenotype, chrs, chrslength, linkerlength, slidingwindowsize, maxmismatchednuminslidingwidow):
    # finding the longest matched window for each starting window
    window_blocks = {}
    MLPinfo = {}

    for chr in chrs:
        window_blocks[chr] = {}
        maxblocklength_formerblock = 0

        for startposi in range(1, chrslength[chr]):
            maxblocklength_presentblock = 0
            tempdata = {}  # recorde the information of the current window

            for sample in targetSamplesGenotype[chr].keys():

                extendendposi, mismatchnum = extendblock(sample, chr, targetSamplesGenotype, startposi, linkerlength,
                                                          slidingwindowsize, maxmismatchednuminslidingwidow)

                if extendendposi >= maxblocklength_presentblock:
                    maxblocklength_presentblock = extendendposi

                    tempdata[sample] = {}
                    tempdata[sample]["e"] = extendendposi
                    tempdata[sample]["mis"] = mismatchnum

            # finding the mostly matched block if there are more than one matched block and they are equal in length
            if maxblocklength_presentblock > maxblocklength_formerblock:
                window_blocks[chr][startposi] = {}
                minmismatchnum = chrslength[chr]

                for sample in tempdata.keys():
                    if tempdata[sample]["e"] == maxblocklength_presentblock:
                        window_blocks[chr][startposi][sample] = {}

                        window_blocks[chr][startposi][sample]["end"] = tempdata[sample]["e"]
                        window_blocks[chr][startposi][sample]["mis"] = tempdata[sample]["mis"]

                        if tempdata[sample]["mis"] < minmismatchnum:
                            minmismatchnum = tempdata[sample]["mis"]

                        ###### recored mapping locations of the MLP varieties
                        if sample not in MLPinfo.keys():
                            MLPinfo[sample] = {}
                        if chr not in MLPinfo[sample].keys():
                            MLPinfo[sample][chr] = {}

                        MLPinfo[sample][chr][startposi] = tempdata[sample]["e"]

                if len(window_blocks[chr][startposi].keys()) > 1:
                    for sample in list(window_blocks[chr][startposi].keys()):
                        if window_blocks[chr][startposi][sample]["mis"] > minmismatchnum:
                            window_blocks[chr][startposi].pop(sample)

                            MLPinfo[sample][chr].pop(startposi)

                maxblocklength_formerblock = maxblocklength_presentblock

    return window_blocks, MLPinfo

def contri_matrix_counter(outfilewriter4, outfilewriter5, querysample, new_ordered_contri, mlpinfo, chrs, chrslength, total_chr_len):
    fill_chrs = {}
    cumulative_contris =[]
    cumulative_samples = []

    for chr in chrs:
        fill_chrs[chr] = "0" * chrslength[chr]

    for i in range(0, len(new_ordered_contri)):
        sample = new_ordered_contri[i][0]
        cumulative_samples.append(sample)

        cumulative_contri = 0
        if sample in mlpinfo.keys():
            for chr in mlpinfo[sample].keys():
                for startposi in mlpinfo[sample][chr].keys():
                    endposi = mlpinfo[sample][chr][startposi]

                    if (endposi - startposi + 1) > 2:
                        genotype = list(fill_chrs[chr])

                        for matchedposi in range(startposi, (endposi + 1)):
                            if genotype[matchedposi - 1] == "0":
                                genotype[matchedposi - 1] = "1"
                                cumulative_contri += 1

                        fill_chrs[chr] = "".join(genotype)

        cumulative_contri = round(cumulative_contri / total_chr_len * 100, 2)
        cumulative_contris.append(cumulative_contri)

    outfilewriter4.write(querysample)
    for sample in cumulative_samples:
        outfilewriter4.write("\t" + str(sample))

    outfilewriter4.write("\n")


    outfilewriter5.write(querysample)
    for contri in cumulative_contris:
        outfilewriter5.write("\t" + str(contri))

    outfilewriter5.write("\n")


def mtp_similarity_counter(outfilewriter3, outfilewriter4, outfilewriter5, outfilewriter6,
                           querysample, targetsamples, transferredGenotype, chrs, chrslength, linkerlength, total_chr_len,
                           slidingwindowsize, maxmismatchednuminslidingwidow):

    similirity_result = {}
    ordered_contri = {}
    for sample in targetsamples:
        similirity_result[sample] = 0

    mtp_variety, mlpinfo = MTPfinding(transferredGenotype, chrs, chrslength, linkerlength, slidingwindowsize, maxmismatchednuminslidingwidow)

    for chr in mtp_variety.keys():
        for startposi in mtp_variety[chr].keys():
            for sample in mtp_variety[chr][startposi].keys():
                endposi = mtp_variety[chr][startposi][sample]["end"]

                #genotype = transferredGenotype[chr][sample]
                #genotype = genotype[(startposi-1):endposi]
                #similirity_result[sample] += genotype.count("1")
                if (endposi - startposi + 1) > 2:
                    similirity_result[sample] += (endposi - startposi + 1)

                    outfilewriter6.write("\t".join([querysample, sample, str(chr), str(startposi), str(endposi)]) + "\n")

    outfilewriter3.write(querysample)
    for sample in sorted(targetsamples):
        simipercent = round(similirity_result[sample] / total_chr_len * 100, 2)
        ordered_contri[sample] = simipercent

        outfilewriter3.write("\t" + str(simipercent))

    outfilewriter3.write("\n")

    #################################### caculating the cumulative matrix
    new_ordered_contri = sorted(ordered_contri.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)

    contri_matrix_counter(outfilewriter4, outfilewriter5, querysample, new_ordered_contri, mlpinfo, chrs, chrslength, total_chr_len)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--haplptype_file')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-q', '--query_file')
    parser.add_argument('-t', '--target_file')
    parser.add_argument('-c', '--chr_length')
    parser.add_argument('-w', '--window_size', type=int)
    parser.add_argument('-o', '--output_path')

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    genotypingFile = args.haplptype_file
    pangenome = args.reference
    querysamplesFile = args.query_file
    targetsamplesFile = args.target_file
    chr_lengthFile = args.chr_length
    windowsize = args.window_size
    outputpath = args.output_path


    chrs = ["chr1A", "chr2A", "chr3A", "chr4A", "chr5A", "chr6A", "chr7A",
            "chr1B", "chr2B", "chr3B", "chr4B", "chr5B", "chr6B", "chr7B",
            "chr1D", "chr2D", "chr3D", "chr4D", "chr5D", "chr6D", "chr7D"]

    chrslength = getchrlength(chr_lengthFile, pangenome, windowsize)

    total_chr_len = 0
    for ichr in chrslength.keys():
        total_chr_len += chrslength[ichr]


    linkerlength = 3 ## It means less than three continuous mismathed windows were allowed in the haplotype comparison
    slidingwindowsize = 10  ###
    maxmismatchednuminslidingwidow = 4  ### maximum mismatched window in a given sliding window

    hap_interval = 5  ## It define the interval number to delimite the haplotype number
    max_smooth_length = 2   # The maximum continous windows which could be smoothed.
    mismatch_penalty_coefficient = 1
    gap_intro_penalty_coefficient = 1
    match_coefficient = 1
    revise_coefficient = 5  ### It is related with how many flanking matched windows were needed to smooth the mismatched window

    threshold = 50  ## The window where there are more than 50 haplotypes is directly regarded as "match" in the haplotype comparison


    querysamples = readsamples(querysamplesFile)
    targetsamples = readsamples(targetsamplesFile)
    querysamples = sorted(querysamples)
    targetsamples = sorted(targetsamples)

    genotyping_format_transfer(querysamples, targetsamples, genotypingFile, outputpath)

    get_weighted_matrix(hap_interval, outputpath)

    weighted_matrix = read_weighted_matrix(outputpath)

    outfile3 = open(os.path.join(outputpath, "matrix_filtered_similirity.txt"), "w+")
    outfile4 = open(os.path.join(outputpath, "matrix_cumulative_sample.txt"), "w+")
    outfile5 = open(os.path.join(outputpath, "matrix_cumulative_percent.txt"), "w+")
    outfile6 = open(os.path.join(outputpath, "matrix_inheritage_blocks.txt"), "w+")

    outfile3.write("query" + "\t" + "\t".join(targetsamples) + "\n")
    outfile4.write("query" + "\t" + "\t".join([str(i) for i in range(1, (len(targetsamples) + 1))]) + "\n")
    outfile5.write("query" + "\t" + "\t".join([str(i) for i in range(1, (len(targetsamples) + 1))])+ "\n")

    for querysample in querysamples:
        print(querysample + " is running")

        transferredGenotype, queryGnotype, targetGenotype = read_Trans_Genotype(querysample, targetsamples, threshold, outputpath)

        ################################################################################
        hap_smooth("0to1", transferredGenotype, weighted_matrix, max_smooth_length, mismatch_penalty_coefficient,
                       gap_intro_penalty_coefficient, match_coefficient, revise_coefficient, chrslength)

        hap_smooth("1to0", transferredGenotype, weighted_matrix, max_smooth_length, mismatch_penalty_coefficient,
                       gap_intro_penalty_coefficient, match_coefficient, revise_coefficient, chrslength)


        ################################################################################
        mtp_similarity_counter(outfile3, outfile4, outfile5, outfile6, querysample,
                               targetsamples, transferredGenotype, chrs, chrslength, linkerlength, total_chr_len,
                               slidingwindowsize, maxmismatchednuminslidingwidow)

    outfile3.close()
    outfile4.close()
    outfile5.close()
    outfile6.close()


if __name__ == '__main__':
    main()

