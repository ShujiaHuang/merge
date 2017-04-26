"""
A general module for merging files.

Author: Shujia Huang
Date: 2017-04-26
"""

import sys
import heapq
import gzip


class FileForQueueing(object):
    """
    """
    def __init__(self, the_file, line):
        """
        Store the file, and initialise the current value
        """
        self.the_file = the_file
        self.finishedReadingFile = False
        self.heap = []

        line = line
        cols = line.strip().split("\t")
        chrom = cols[0]

        # Where possible, convert chromosome names into
        # integers for sorting. If not possible, use
        # original names.
        try:
            chrom = int(chrom.upper().strip("CHR"))
        except Exception:
            pass

        pos = int(cols[1])
        heapq.heappush(self.heap, (chrom, pos, line))

        while not self.finishedReadingFile and len(self.heap) < 100:

            try:
                line = self.the_file.next()
                cols = line.strip().split("\t")
                chrom = cols[0]

                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
            except StopIteration:
                self.finishedReadingFile = True
                break

            heapq.heappush(self.heap, (chrom, pos, line))

        # Now take the top line
        self.chrom, self.pos, self.line = heapq.heappop(self.heap)

    def __cmp__(self, other):
        """
        Comparison function. Utilises the comparison function defined in
        the AlignedRead class.
        """
        return cmp(self.chrom, other.chrom) or cmp(self.pos, other.pos)

    def __del__(self):
        """
        Destructor
        """
        self.the_file.close()

    def next(self):
        """
        Increment the iterator and yield the new value. Also, store the
        current value for use in the comparison function.
        """
        if not self.finishedReadingFile:

            try:
                line = self.the_file.next()
                #cols = line.strip().split('\t')
                cols = line.strip().split()
                chrom = cols[0]

                # Where possible, convert chromosome names into
                # integers for sorting. If not possible, use
                # original names.
                try:
                    chrom = int(chrom.upper().strip("CHR"))
                except Exception:
                    pass

                pos = int(cols[1])
                heapq.heappush(self.heap, (chrom, pos, line))

            except StopIteration:
                self.finishedReadingFile = True

        if len(self.heap) != 0:
            # Now take the top line
            self.chrom, self.pos, self.line = heapq.heappop(self.heap)
        else:
            raise StopIteration


def expandedOpen(path, mode):
    try:
        return open(path, mode)
    except IOError:
        return open(os.path.expanduser(path), mode)


def Open(fileName, mode, compressLevel=9):
    """
    Function that allows transparent usage of dictzip, gzip and
    ordinary files
    """
    if fileName.endswith(".gz") or fileName.endswith(".GZ"):
        fileDir = os.path.dirname(fileName)
        if os.path.exists(fileDir):
            return gzip.GzipFile(fileName, mode, compressLevel)
        else:
            return gzip.GzipFile(os.path.expanduser(fileName), mode,
                                 compressLevel)
    else:
        return expandedOpen(fileName, mode)


def merge_files(temp_file_names, final_file_name):
    """
    Merging output VCF/CVG files into a final file
    """

    # Final output file
    if final_file_name == "-":
        output_file = sys.stdout
    else:
        output_file = Open(final_file_name, 'wb')
    the_heap = []

    # Initialise queue
    for index, file_name in enumerate(temp_file_names):
        the_file = Open(file_name, 'rb')

        for line in the_file:

            # End of this file
            if line[0] == "#":
                if index == 0:
                    output_file.write(line)
            else:
                the_file_for_queueing = FileForQueueing(the_file, line)
                heapq.heappush(the_heap, the_file_for_queueing)
                break

        # If there are no calls in the temp file, we still want to
        # remove it.
        else:
            the_file.close()

    # Merge-sort the output using a priority queue
    while len(the_heap) != 0:

        # Get file from heap in right order
        next_file = heapq.heappop(the_heap)
        output_file.write(next_file.line)

        # Put file back on heap
        try:
            next_file.next()
            heapq.heappush(the_heap, next_file)
        except StopIteration:
            continue

    # Close final output file
    if final_file_name != "-":
        output_file.close()

    #log.info("Finished merging %s file(s)"%final_file_name)
    return
