#!/usr/bin/env python3

import sys
import signal
from ffdb.seq import Seq


class Timeout(Exception):
    pass


def raiser(err):
    raise err


def stringify(handle):
    for line in handle:
        yield line.strip().decode("utf-8")
    return


def main():

    first_id = None
    first_length = None
    first_end_was_x = False
    output = []

    signal.signal(signal.SIGALRM, lambda i, j: raiser(Timeout()))

    try:
        signal.alarm(60 * 15)
        infile = sys.stdin.buffer.read()
        for seq in Seq.parse(stringify(infile.split(b"\n"))):
            seq.seq = seq.seq.replace(b"*", b"X")

            if first_length is None:
                first_id = seq.id
                first_length = len(seq)
                first_end_was_x = seq.seq[-1:].upper() == b"X"
            elif len(seq) == first_length:
                pass
            elif (len(seq) == first_length - 1) and first_end_was_x:
                seq.seq = seq.seq + b"-"
            else:
                raise ValueError((
                    "Encountered MSA with abnormal number of columns, "
                    f"or with first seq missing X. Seq was {seq.id}."
                ))

            output.append(bytes(str(seq), encoding="utf-8"))
        signal.alarm(0) 
    except Timeout:
        print(
            "Process ", first_id, "timed out. Got through", len(output),
            file=sys.stderr
        )
    finally:
        sys.stdout.buffer.write(b"\n".join(output))

    return


if __name__ == "__main__":
    main()
