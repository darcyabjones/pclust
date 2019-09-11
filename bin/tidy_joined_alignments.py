#!/usr/bin/env python3

import sys
from ffdb.seq import Seq


def main():

    first_length = None
    first_end_was_x = False

    for seq in Seq.parse(sys.stdin):
        seq.seq = seq.seq.replace(b"*", b"X")

        if first_length is None:
            first_length = len(seq)
            first_end_was_x = seq.seq[-1:].upper() == b"X"
        elif len(seq) == first_length - 1 and first_end_was_x:
            seq.seq = seq.seq + b"-"
        else:
            raise ValueError((
                "Encountered MSA with abnormal number of columns, "
                f"or with first seq missing X. Seq was {seq.id}."
            ))

        print(seq, file=sys.stdout)

    return


if __name__ == "__main__":
    main()
