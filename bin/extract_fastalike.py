#!/usr/bin/env python3

import sys


def check_last_was(block):
    if len(block) == 0:
        return False
    else:
        return block[-1].startswith(">")


def get_block_filename(block):
    header = block[0]
    basename = header.lstrip("> ").split(" ", 1)[0]
    return basename + ".fasta"


def main():

    if len(sys.argv) != 2:
        print("Expected a file to process")
        sys.exit(1)

    fastalike = sys.argv[1]

    current = []
    with open(fastalike, "r") as handle:
        _ = next(handle)  # Skip first line
        for line in handle:
            line = line.strip()
            if line.startswith(">") and check_last_was(current):
                assert len(current) > 1

                filename = get_block_filename(current)
                with open(filename, "w") as out_handle:
                    print("\n".join(current[:-1]), file=out_handle)
                current = []

            current.append(line)

        assert len(current) > 1

        filename = get_block_filename(current)
        with open(filename, "w") as out_handle:
            print("\n".join(current[:-1]), file=out_handle)


if __name__ == "__main__":
    main()
