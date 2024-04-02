#!/usr/bin/env python

import sys
import gzip
import re
from io import BufferedReader
from subprocess import check_output
from os.path import exists


class VcfRecord:
    def __init__(self, inline):
        tokens = inline.strip().split('\t')

        self.chr = tokens[0]
        self.pos = int(tokens[1])
        self.vid = tokens[2]
        self.ref = tokens[3]
        self.alt = tokens[4]
        self.qual = tokens[5]
        self.filter = tokens[6]
        self.info = tokens[7].split(';')
        self.others = "\t".join(tokens[8:])

        # Create a dictionary for INFO
        self.infoDict = {}
        for infoItem in self.info:
            items = infoItem.split('=')
            if len(items) == 1:
                self.infoDict[items[0]] = True
            elif len(items) > 1:
                self.infoDict[items[0]] = items[1]

        self.mateChr = self.chr
        if "END" in self.infoDict:
            self.matePos = int(self.infoDict["END"])
        else:
            self.matePos = -1

        if self.infoDict["SVTYPE"] == "BND":
            items = self.alt.replace("[", " ").replace("]", " ").split(" ")
            [self.mateChr, matePos] = items[1].split(':')
            self.matePos = int(matePos)
            self.mateChr = self.mateChr.replace("CHR", "chr")

    def makeLine(self):
        infoStr = ";".join(self.info)
        self.line = "\t".join((self.chr, str(self.pos), self.vid, self.ref, self.alt, self.qual, self.filter, infoStr,
                               self.others)) + "\n"


def scanVcf(vcfFile):
    bnd_mate_dict = {}
    chr_dict = {}

    if vcfFile.endswith('gz'):
        gzfp = gzip.open(vcfFile, 'rb')
        fpVcf = BufferedReader(gzfp)
    else:
        fpVcf = open(vcfFile, 'rb')

    for line in fpVcf.readlines():
        line = line.decode()
        if line[0] == '#':
            if "##contig" in line:
                m = re.match("##contig=<ID=(.+),length=(\d+)>", line)
                chr_dict[m.group(1)] = int(m.group(2))
            continue

        vcfRec = VcfRecord(line)
        if "BND" in vcfRec.infoDict["SVTYPE"]:
            if vcfRec.vid in bnd_mate_dict:
                # update mate INFO
                bnd_mate_dict[vcfRec.vid] = vcfRec.infoDict
            else:
                mateId = vcfRec.infoDict["MATEID"]
                bnd_mate_dict[mateId] = ""

    return bnd_mate_dict, chr_dict


def writeLines(lines):
    for line in lines:
        sys.stdout.write(line)


def convertBND(vcfFile, bndMateDict, chrDict):
    isHeaderInfoAdded = False
    lineBuffer = []
    bufferedChr = ""
    bufferedPos = -1

    if vcfFile.endswith('gz'):
        gzfp = gzip.open(vcfFile, 'rb')
        fpVcf = BufferedReader(gzfp)
    else:
        fpVcf = open(vcfFile, 'rb')

    for line in fpVcf:
        line = line.decode()
        if line.startswith('#'):
            if (not isHeaderInfoAdded) and line.startswith("##FORMAT="):
                sys.stdout.write(
                    "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">\n")
                sys.stdout.write(
                    "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
                sys.stdout.write(
                    "##INFO=<ID=CHREND,Number=1,Type=String,Description=\"The short distance to the chromosome end\">\n")
                isHeaderInfoAdded = True

            sys.stdout.write(line)
            continue

        vcfRec = VcfRecord(line)
        vcfRec.info.append("CHR2=%s" % vcfRec.mateChr)
        vcfRec.info.append("POS2=%d" % int(vcfRec.matePos))
        vcfRec.info.append("CHREND=%d,%d" % (min(vcfRec.pos - 1, chrDict[vcfRec.chr] - vcfRec.pos), min(vcfRec.matePos - 1, chrDict[vcfRec.mateChr] - vcfRec.matePos)))

        # skip mate record
        if vcfRec.vid in bndMateDict:
            continue

        if "SVLEN" not in vcfRec.infoDict:
            vcfRec.info.append("SVLEN=-1")

        if vcfRec.infoDict["SVTYPE"] == "BND":
            # update manta ID
            vidSuffix = vcfRec.vid.split("DRAGEN:BND")[1]
            idx = vidSuffix.rfind(':')
            vcfRec.vid = "DRAGEN:BND%s" % vidSuffix[:idx]

            newInfo = []
            for infoItem in vcfRec.info:
                # skip BND-specific tags
                if infoItem.startswith("MATEID") or infoItem.startswith("BND_DEPTH") or infoItem.startswith("MATE_BND_DEPTH"):
                    continue

                # update event ID
                elif infoItem.startswith("EVENT"):
                    eidSuffix = vcfRec.infoDict["EVENT"].split("DRAGEN:BND")[1]
                    idx = vidSuffix.rfind(':')
                    infoEventStr = "EVENT=DRAGEN:BND%s" % eidSuffix[:idx]
                    newInfo.append(infoEventStr)

                # apply all other tags
                else:
                    newInfo.append(infoItem)

            vcfRec.info = newInfo

        vcfRec.makeLine()

        # make sure the vcf is sorted in genomic order
        if (not vcfRec.chr == bufferedChr) or (vcfRec.pos > bufferedPos):
            if lineBuffer:
                writeLines(lineBuffer)

            lineBuffer = [vcfRec.line]
            bufferedChr = vcfRec.chr
            bufferedPos = vcfRec.pos
        elif vcfRec.pos < bufferedPos:
            lineBuffer.insert(0, vcfRec.line)
        else:
            lineBuffer.append(vcfRec.line)

    if lineBuffer:
        writeLines(lineBuffer)


if __name__ == '__main__':

    usage = "convertBND.py <vcf file>\n"
    if len(sys.argv) <= 1:
        sys.stderr.write(usage)
        sys.exit(1)

    vcfFile = sys.argv[1]
    bndMateDict, chrDict = scanVcf(vcfFile)
    convertBND(vcfFile, bndMateDict, chrDict)
