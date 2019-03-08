#!/usr/bin/env python3

from collections import namedtuple
from os.path import splitext
from os.path import basename
from io import BytesIO
from shutil import copyfileobj
from copy import deepcopy
from collections.abc import Iterator

import sys
import argparse


IndexRow = namedtuple("IndexRow", ["name", "start", "size"])


class FFIndex(object):

    def __init__(self, index=None):
        if index is None:
            index = []
        else:
            assert all(isinstance(r, IndexRow) for r in index)

        self.index = sorted(index, key=lambda x: x.start)
        self.lookup = dict()

        for i, idx in enumerate(self.index):
            # A well formed ffindex should never have duplicate names.
            assert idx.name not in self.lookup
            self.lookup[idx.name] = i

        return

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return self.index[key]
        elif isinstance(key, (str, bytes)):
            return self.index[self.lookup[key]]
        else:
            raise ValueError("Expected either a string, an int, or a slice.")

    def __contains__(self, key):
        return key in self.lookup

    def __iter__(self):
        for row in self.index:
            yield row
        return

    def __len__(self):
        return len(self.index)

    @classmethod
    def from_file(cls, handle):
        indices = cls._parse_ffindex(handle)
        return cls(index=indices)

    @staticmethod
    def _parse_ffindex(handle):
        ffindex = list()
        for line in handle:
            name, start, size = line.strip().split()
            row = IndexRow(name, int(start), int(size))  # size - 1?
            ffindex.append(row)

        return ffindex

    def append(self, value):
        assert isinstance(value, IndexRow)

        name, start, size = value

        assert name not in self

        if len(self.index) > 0:
            last_row = self.index[-1]
            last_end = last_row.start + last_row.size
            start = last_end
        else:
            start = 0

        self.index.append(IndexRow(name, start, size))

        self.lookup[name] = len(self.index)
        return

    def extend(self, values):
        for value in values:
            self.append(value)
        return len(values)

    def write_to(self, handle):
        length = 0
        for ind in sorted(self.index, key=lambda x: x.name):
            line = "{}\t{}\t{}\n".format(
                ind.name.decode("utf-8"),
                ind.start,
                ind.size
            )
            length += handle.write(line.encode())

        return length

    def bump_starts(self, by=0):
        if by == 0:
            return deepcopy(self)

        new_index = []
        for name, start, size in self.index:
            new_index.append(IndexRow(name, start + by, size))

        # Construct a new object
        return self.__class__(new_index)


class FFData(object):

    def __init__(self, handle):
        self.handle = handle

    def __getitem__(self, key):
        if isinstance(key, IndexRow):
            name, start, size = key
            self.handle.seek(start)
            return self.handle.read(size)
        elif isinstance(key, list):
            records = []
            for name, start, size in key:
                self.handle.seek(start)
                records.append(self.handle.read(size))
            return records
        else:
            raise ValueError("Must be an IndexRow or a list of IndexRows")

    def append(self, b):
        self.handle.seek(0, 2)  # Go to end of file.
        assert b[-1:] == b"\0"
        return self.handle.write(b)

    def write_to(self, handle):
        self.handle.seek(0)
        return copyfileobj(self.handle, handle)

    def write_sized(self, start, size, handle):
        self.handle.seek(start)
        return handle.write(self.handle.read(size))


class FFDB(object):

    def __init__(self, data, index):
        self.data = data
        self.index = index
        return

    @classmethod
    def from_file(cls, data_handle, index_handle):
        data = FFData(data_handle)
        index = FFIndex.from_file(index_handle)
        return cls(data, index)

    @classmethod
    def new(cls, data_handle=None):
        if data_handle is None:
            data_handle = BytesIO()

        data = FFData(data_handle)
        index = FFIndex()
        return cls(data, index)

    def __getitem__(self, key):
        indices = self.index[key]
        return self.data[indices]

    def __contains__(self, key):
        return key in self.index

    def __len__(self):
        return len(self.index)

    def append(self, data, key):
        if isinstance(key, IndexRow):
            this_key = key
        else:
            this_key = data.index[key]

        self.index.append(this_key)
        to_write = data.data[this_key]
        return self.data.append(to_write)

    def extend(self, data, keys):
        if isinstance(keys, slice):
            keys = data.index[keys]
        elif keys is None:
            keys = data.index

        length = 0
        for key in keys:
            length += self.append(data, key)
        return length

    def write_to(self, data_handle, index_handle):
        assert data_handle.tell() == 0

        l1 = self.index.write_to(index_handle)
        l2 = self.data.write_to(data_handle)
        return l1, l2

    def concat(self, dbs):
        # Go to end of file
        self.data.handle.seek(0, 2)
        for db in dbs:
            self.index.extend(db.index.index)
            db.data.write_to(self.data.handle)
        return

    def partition(self, name, template="{name}_{index}.{ext}", n=10000):
        start_pos = 0
        pindices = []
        partition = 1

        for i, p in enumerate(self.index, 1):
            if i % n == 0:
                self._write_partition(
                    start_pos,
                    p.start,
                    template,
                    name,
                    pindices,
                    partition
                )

                pindices = []
                partition += 1
                start_pos = p.start

            pindices.append(p)

        if len(pindices) > 1:
            end = pindices[-1].start + pindices[-1].size
            self._write_partition(
                start_pos,
                end,
                template,
                name,
                pindices,
                partition
            )

        return partition

    def _write_partition(self, start, end, template, name, indices, partition):
        size = (end - start)

        ffindex_name = template.format(
            name=name,
            index=partition,
            ext="ffindex"
        )

        ffdata_name = template.format(
            name=name,
            index=partition,
            ext="ffdata"
        )

        partition_index = FFIndex(indices).bump_starts(by=(-1 * start))

        with open(ffindex_name, "wb") as handle:
            partition_index.write_to(handle)

        with open(ffdata_name, "wb") as handle:
            self.data.write_sized(start, size, handle)

        return


class Seq(object):

    def __init__(self, id, desc, seq):
        """ Construct a new Seq object.
        Keyword arguments:
        id -- The sequence id <str>.
        desc -- A short description of the sequence <str>.
        seq -- The biological sequence <str>.
        """
        self.id = id
        self.desc = desc
        self.seq = seq
        return

    def __str__(self):
        """ Returns a FASTA string from the object """
        line_length = 60

        if self.desc is None:
            lines = [">{}".format(self.id)]
        else:
            lines = [">{} {}".format(self.id, self.desc)]

        for i in range(0, len(self), line_length):
            lines.append(self.seq[i:i+line_length].decode("utf-8"))

        return "\n".join(lines) + "\n"

    def __repr__(self):
        """ Returns a simple string representation of the object. """
        cls = self.__class__.__name__
        return "{}(id='{}', desc='{}', seq='{}')".format(cls, self.id,
                                                         self.desc, self.seq)

    def __getitem__(self, key):
        """ Allow us to access indices from the seq directly. """
        seq = self.seq[key]
        return self.__class__(self.id, self.desc, seq)

    def __eq__(self, other):
        """ Allows us to compare two Seq objects directly using '==' .
        NB. python internally implements != based on this too.
        """
        if isinstance(other, self.__class__):
            return self.seq == other.seq
        elif isinstance(other, str):
            return self.seq == other
        else:
            raise ValueError((
                "Equality comparisons not implemented between {} and {}."
                ).format(type(self), type(other)))
        return

    def __len__(self):
        """ The length of a Seq object should be the length of the seq. """
        return len(self.seq)

    @classmethod
    def read(cls, handle):
        """ Read a single FASTA record.
        Parses a single FASTA record into a Seq object.
        Assumes that the first line will always contain the id line,
        and that there is a single FASTA sequence.
        Keyword arguments:
        handle -- An iterable containing lines (newline split) of the file.
        Returns:
        A Seq object.
        """
        if not isinstance(handle, Iterator):
            handle = iter(handle)

        try:
            id_, desc = cls._split_id_line(next(handle).strip())
        except ValueError as e:
            raise ValueError("Fasta parsing failed. " + str(e))

        # tuple comprehensions are generators so we're still doing lazy eval
        seq = "".join((l.strip() for l in handle))
        return Seq(id_, desc, seq.encode())

    @classmethod
    def parse(cls, handle):
        """ Parse multiple fasta records.
        Parses a multi-fasta formatted file-like object.
        Keyword arguments:
        handle -- A file-like object or any iterable over the fasta file lines.
        Yields:
        Seq objects.
        """

        # Store the initial state to avoid outputting empty record.
        first = True
        # Store lines for this block here.
        current_record = []

        for line in handle:
            if line.startswith(">"):
                if not first:
                    # Yield makes this function a generator.
                    # NB we reuse the read method to avoid repetition.
                    # It's also easier to test.
                    yield cls.read(current_record)
                else:
                    # Once we've passed the first sequence this passes.
                    first = False

                # Start a new block
                current_record = [line]
            else:
                # For lines containing sequences we simply append the sequence.
                current_record.append(line)

        # The last sequence in the file won't have a ">" following it.
        # so we yield the last block too.
        yield cls.read(current_record)
        return

    @classmethod
    def parse_many(cls, handles):
        for handle in handles:
            for record in cls.parse(handle):
                yield record
        return

    @staticmethod
    def _split_id_line(line):
        """ Parse the FASTA header line into id and description components.
        NB expects the '>' character to be present at start of line.
        Keyword arguments:
        line -- A string containing the header.
        Returns:
        Tuple -- id and description strings.
        """

        if not line.startswith(">"):
            raise ValueError(("Encountered malformed fasta header. "
                              "Offending line is '{}'").format(line))
        # Strip the ">" character and split at most 1 time on spaces.
        sline = line[1:].split(" ", 1)

        if len(sline) == 1:
            return sline[0], None
        else:
            return sline[0], sline[1]

        # We should never reach this point.
        return

    def checksum(self):
        from hashlib import sha1
        from base64 import b64encode
        hash_ = sha1(self.seq).digest()
        return b64encode(hash_).rstrip(b"=").decode("utf-8")


def cli(prog, args):
    """ Process command line arguments.
    Often this is simply put in the main function or the
    if __name__ == "__main__" block, but keeping it as a function allows
    testing if we had more complex command line interfaces.
    Keyword arguments:
    prog -- The name of the program. Usually this will be the first element of
        the sys.argv list.
    args -- The command line arguments for the program. Usually this will be
        a slice from i.e sys.argv[1:].
    Returns:
    An argparse Args object containing the parsed args.
    Raises:
    Standard argparse exceptions if the CLI conditions are not met.
    Required parameters, type constraints etc.
    """

    parser = argparse.ArgumentParser(
        prog=prog,
        description=""
    )

    subparsers = parser.add_subparsers(dest='subparser_name')

    split_subparser = subparsers.add_parser(
        "split",
        help="Split an ffindex database into n partitions."
    )

    # Add split parameters with sep function
    cli_split(split_subparser)

    combine_subparser = subparsers.add_parser(
        "combine",
        help="Collects many ffindex databases into a single one."
    )

    cli_combine(combine_subparser)

    fasta_subparser = subparsers.add_parser(
        "fasta",
        help=("Creates an ffindex database from a multifasta, "
              "with many sequences per document.")
    )

    cli_fasta(fasta_subparser)

    collect_subparser = subparsers.add_parser(
        "collect",
        help=("Collects all ffdata documents into a single file."
              "Essentially it just filters out any null bytes and makes "
              "sure there is a newline between documents.")
    )

    cli_collect(collect_subparser)

    parsed = parser.parse_args(args)

    # Validate arguments passed to combine
    if parsed.subparser_name in ("combine", "collect"):
        files = []
        files.extend(parsed.ffdata)
        files.extend(parsed.ffindex)

        if len(files) % 2 != 0:
            parser.error((
                "There should be the same number of ffindex and "
                "ffdata files provided to `combine`."
            ))

        print(files)
        parsed.ffdata = files[:len(files)//2]
        parsed.ffindex = files[len(files)//2:]

    return parsed


def cli_collect(parser):
    parser.add_argument(
        "-t", "--trim",
        type=int,
        default=None,
        help=("Trim this many lines from the start of each document. "
              "Useful for headers in csv documents."),
    )

    parser.add_argument(
        "ffdata",
        metavar="FFDATA",
        nargs="+",
        type=argparse.FileType('rb'),
        help="The ffindex .ffdata file.",
    )

    parser.add_argument(
        "ffindex",
        metavar="FFINDEX",
        nargs="+",
        type=argparse.FileType('rb'),
        help="The ffindex .ffindex file.",
    )


def cli_combine(parser):
    parser.add_argument(
        "-d", "--data",
        required=True,
        type=argparse.FileType('wb'),
        help="The path to write the ffdata file to.",
    )

    parser.add_argument(
        "-i", "--index",
        required=True,
        type=argparse.FileType('wb'),
        help="The path to write the ffindex file to.",
    )

    parser.add_argument(
        "ffdata",
        metavar="FFDATA",
        nargs="+",
        type=argparse.FileType('rb'),
        help="The ffindex .ffdata file.",
    )

    parser.add_argument(
        "ffindex",
        metavar="FFINDEX",
        nargs="+",
        type=argparse.FileType('rb'),
        help="The ffindex .ffindex file.",
    )

    return


def cli_fasta(parser):
    parser.add_argument(
        "-d", "--data",
        required=True,
        type=argparse.FileType('wb'),
        help="The path to write the ffdata file to.",
    )

    parser.add_argument(
        "-i", "--index",
        required=True,
        type=argparse.FileType('wb'),
        help="The path to write the ffindex file to.",
    )

    parser.add_argument(
        "-c", "--checksum",
        default=False,
        action="store_true",
        help="Replace the sequence ids with checksums and remove duplicates",
    )

    parser.add_argument(
        "-n", "--size",
        type=int,
        default=1,
        help="The number of fasta records to use per document.",
    )

    parser.add_argument(
        "-f", "--filter",
        type=argparse.FileType('r'),
        default=None,
        help="A newline delimited file of ids or checksums to filter with."
    )

    parser.add_argument(
        "--mapping",
        type=argparse.FileType('w'),
        default=None,
        help="Write the mapping of seqids to checksums to file.",
    )

    parser.add_argument(
        "fasta",
        metavar="FASTA",
        nargs="+",
        type=argparse.FileType('r'),
        help="The fasta files to pull in.",
    )

    return


def cli_split(parser):
    parser.add_argument(
        "-n", "--size",
        type=int,
        default=100000,
        help="The number of records for each partition to have.",
    )

    parser.add_argument(
        "-b", "--basename",
        type=str,
        default="{name}_{index}.{ext}",
        help=(
            "The output database partition names. "
            "Can use python format syntax. "
            "Some values are available, `name` will be the basename of the "
            "input database (no extensions), `index` will be the 1-based "
            "partition number, and ext will be ffindex or ffdata as "
            "appropriate."
        )
    )

    parser.add_argument(
        "ffdata",
        metavar="FFDATA_FILE",
        type=argparse.FileType('rb'),
        help="The ffindex .ffdata files.",
    )

    parser.add_argument(
        "ffindex",
        metavar="FFINDEX_FILE",
        type=argparse.FileType('rb'),
        help="The ffindex .ffindex files.",
    )

    return


def simplename(path):
    splitext(basename(path))[0]
    return


def ffsplit(args):
    ffdb = FFDB.from_file(args.ffdata, args.ffindex)

    file_basename = simplename(args.ffdata.name)
    ffdb.partition(
        name=file_basename,
        template=args.basename,
        n=args.size,
    )
    return


def ffcombine(args):
    outdb = FFDB.new(args.data)

    indbs = []
    for (data, index) in zip(args.ffdata, args.ffindex):
        indb = FFDB.from_file(data, index)
        indbs.append(indb)

    # Writes to ffdata since new was provided handle.
    outdb.concat(indbs)
    outdb.index.write_to(args.index)
    return


def from_fasta(args):
    outdb = FFDB.new(args.data)

    if args.filter is not None:
        filter_out = {
            f.strip()
            for f
            in args.filter
        }

    seen = set()

    chunk_data = bytearray()
    chunk_name = None
    chunk_size = 1

    for record in Seq.parse_many(args.fasta):
        if args.checksum:
            checksum = record.checksum()

            if args.filter:
                if checksum in filter_out:
                    continue

            print("{}\t{}".format(record.id, checksum), file=args.mapping)

            if checksum in seen:
                continue
            else:
                record.id = checksum
                record.desc = None
                seen.add(checksum)

        elif args.filter:
            if record.id in filter_out:
                continue

        chunk_data.extend(str(record).encode())

        # Handles first case after write, or just first case.
        if chunk_name is None:
            chunk_name = record.id

        if chunk_size % args.size != 0:
            chunk_size += 1
            continue

        chunk_data.extend(b'\0')
        index = IndexRow(chunk_name.encode(), 0, len(chunk_data))

        outdb.data.append(chunk_data)
        outdb.index.append(index)

        chunk_data = bytearray()
        chunk_name = None
        chunk_size = 1

    if chunk_name is not None:
        chunk_data.extend(b'\0')
        index = IndexRow(chunk_name.encode(), 0, len(chunk_data))

        outdb.data.append(chunk_data)
        outdb.index.append(index)

    outdb.index.write_to(args.index)
    return


def collect(args):
    outfile = sys.stdout.buffer

    for (data, index) in zip(args.ffdata, args.ffindex):
        db = FFDB.from_file(data, index)
        for index in db.index:

            # Take up to :-1 to strip the null byte
            document = db.data[index][:-1]

            if args.trim is not None:
                sdocument = document.split(b'\n')[args.trim:]
                if len(sdocument) == 0:
                    continue
                if sdocument[-1] != b'':
                    # Adds a newline when we join.
                    sdocument.append(b'')

                outfile.write(b'\n'.join(sdocument))
            else:
                outfile.write(document)
                if not document.endswith(b'\n'):
                    outfile.write(b'\n')
    return


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    try:
        if args.subparser_name == "split":
            ffsplit(args)
        elif args.subparser_name == "combine":
            ffcombine(args)
        elif args.subparser_name == "fasta":
            from_fasta(args)
        elif args.subparser_name == "collect":
            collect(args)
        else:
            raise ValueError("I shouldn't reach this point ever")
    except EnvironmentError as e:
        print((
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again."
        ), file=sys.stderr)
        raise e
    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print((
            "I'm so sorry, but we've encountered an unexpected error.\n"
            "This shouldn't happen, so please file a bug report with the "
            "authors.\nWe will be extremely grateful!\n\n"
        ), file=sys.stderr)
        raise e
    return


if __name__ == "__main__":
    main()
