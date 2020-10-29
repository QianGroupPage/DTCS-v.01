def read_line_blocks(file_path, encoding=None):
    with open(file_path, "r+", encoding=encoding) as f:
        lines = f.readlines()

        # Split into lists of lines by newline characters
        size = len(lines)
        idx_list = [idx for idx, val in
                    enumerate(lines) if val == "\n"]

        splitted_lines = [lines[i: j] for i, j in
                          zip([0] + [idx + 1 for idx in idx_list], idx_list +
                              ([size] if idx_list[-1] != size else [])) if i != j]
        return splitted_lines
