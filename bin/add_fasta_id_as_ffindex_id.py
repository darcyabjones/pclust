#!/usr/bin/env python3

import re
import sys
import argparse
import mmap

from typing import List, Sequence, Iterator
from typing import Optional
from typing import BinaryIO

from ffdb.ffindex import FFDB, FFIndex, IndexRow
from ffdb.exceptions import FFError, EmptySequenceError
from ffdb.seq import Seq


def cli(prog: str, args: Sequence[str]) -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        prog=prog,
        description="Add first fasta id as ffindex index ids."
    )

    parser.add_argument(
        "-i", "--index",
        required=True,
        type=argparse.FileType('wb'),
        help="The path to write the ffindex file to.",
    )

    parser.add_argument(
        "-s", "--strip-consensus",
        action="store_true",
        default=False,
        help=(
            "Strip the suffix '_consensus' from the string before adding it "
            "to the id."
        )
    )

    parser.add_argument(
        "--mmap",
        action="store_true",
        default=False,
        help=("Memory map the input hhdata file before reading chunks. "
              "This will significantly reduce IO overhead when doing balanced "
              "or sorted chunks, but requires enough memory to store the "
              "entire ffdata file."),
    )

    parser.add_argument(
        "ffdata",
        metavar="FFDATA_FILE",
        type=argparse.FileType('r+b'),
        help="The ffindex .ffdata files.",
    )

    parser.add_argument(
        "ffindex",
        metavar="FFINDEX_FILE",
        type=argparse.FileType('rb'),
        help="The ffindex .ffindex files.",
    )

    return parser.parse_args(args)


class Timeout(Exception):
    pass


def filter_null(handle: Iterator[str]) -> Iterator[str]:
    for line in handle:
        sline = line.strip()
        if sline != b'\0':
            yield sline
    return


def run(
    ffdata: BinaryIO,
    ffindex: BinaryIO,
    index: BinaryIO,
    should_mmap: bool,
    strip_consensus: bool,
):
    try:
        if should_mmap:
            mm: Optional[mmap.mmap] = mmap.mmap(ffdata.fileno(), 0)
            ffdb = FFDB.from_file(mm, ffindex)
        else:
            mm = None
            ffdb = FFDB.from_file(ffdata, ffindex)

        # Store the new ffindex records with new names here.
        new_indices: List[IndexRow] = []

        # Loop through each of the original indices.
        for old_index in ffdb.index:

            # Get the contents of the data for that index
            data = ffdb.data[old_index]

            # Parse the contents as a fasta file, and get the id of the first
            # record.
            seqs = Seq.parse(filter_null(data.split(b'\n')))
            try:
                new_name = next(seqs).id
            except StopIteration:
                EmptySequenceError((
                    "Encountered an empty record at ffindex "
                    f"document with name {old_index.name}."
                ))

            if strip_consensus:
                new_name = re.sub(r"_consensus$", '', new_name)

            # Create a new index record with the new name.
            new_index = IndexRow(
                new_name.encode('utf-8'),
                old_index.start,
                old_index.size
            )

            new_indices.append(new_index)

        # Write the updated ffindex to a file.
        # Note that we haven't touched the ffdata file, so there's no need to
        # rewrite it.
        new_ffindex = FFIndex(new_indices)
        new_ffindex.write_to(index)

    finally:
        if mm is not None:
            mm.close()
    return


def main():  # noqa
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    try:
        run(
            ffdata=args.ffdata,
            ffindex=args.ffindex,
            index=args.index,
            should_mmap=args.mmap,
            strip_consensus=args.strip_consensus,
        )
    except FFError as e:
        print(f"Error: {e.msg}")
        sys.exit(e.ecode)

    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(130)

    except BrokenPipeError:
        # Pipes get closed and that's normal
        sys.exit(0)

    except EnvironmentError as e:
        print((
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again."
        ), file=sys.stderr)
        raise e

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
