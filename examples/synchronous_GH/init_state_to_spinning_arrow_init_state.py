import surface_crns.readers.manifest_readers as manifest_readers
import sys

'''
Converts the initial state from a surface CRN manifest file into an equivalent
initial state for a spinning-arrow synchronous emulation surface CRN.

Intended for use with the Greenberg-Hastings examples. No guarantee of success
in other examples.
'''

def main():
    input_filename = sys.argv[1]
    manifest_options = manifest_readers.read_manifest(input_filename)
    init_state = manifest_options["init_state"].transpose()

    output_filename = input_filename[:input_filename.rfind(".")] \
                        + "_spinning_arrow_init.txt"
    with open(output_filename, 'w') as outfile:
        for j in range(init_state.shape[1] + 2):
            outfile.write("Edge_D ")
        for i in range(init_state.shape[0]):
            outfile.write("\nEdge_R ")
            for j in range(init_state.shape[1]):
                orig_state = init_state[i, j]
                position   = str(j%3 + 3*(i%3) + 1)
                orientation = "D" if (i+j)%2 == 0 else "U"
                outfile.write(f"{orig_state}_{orientation}_{position}_None ")
            outfile.write("Edge_L")
        outfile.write("\n")
        for j in range(init_state.shape[1] + 2):
            outfile.write("Edge_U ")


if __name__ == "__main__":
    main()