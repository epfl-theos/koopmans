from shutil import copy


def extract(filename_in, filename_out, start=0, end=None, heading=None):
    with open(filename_in) as fd:
        flines = fd.readlines()

    # If heading is provided, find the line number of the heading
    if heading:
        for i, line in enumerate(flines):
            if heading in line:
                start = i
                break
        else:
            raise ValueError(f'Heading {heading} not found in file {filename_in}')

        indent = len(flines[start]) - len(flines[start].lstrip())

        for i in range(start + 1, len(flines)):
            line = flines[i]
            if not line.startswith(' ' * (indent + 1)):
                end = i
                break

        if end is None:
            raise ValueError(f'Could not find the end of the {heading} block')

    flines = flines[start:end]

    # Find the shortest leading whitespace and strip this from all lines (as otherwise the markdown will be rendered as a code block)
    min_indent = min(len(line) - len(line.lstrip()) for line in flines)
    flines = [line[min_indent:] for line in flines]

    # Manual conversion of double spaces to /
    flines = [line[:-2] + '\ \n' if (line.endswith('  \n')
                                     and not line.strip().startswith('-')) else line for line in flines]

    with open(filename_out, 'w') as fd:
        fd.writelines(flines)


if __name__ == '__main__':

    # Tutorial 3
    extract('../../../tutorials/tutorial_3/01-ki/zno.md',
            'tutorial_3/md_excerpts/zno_wannierize_section.md', heading='Wannierize')
    extract('../../../tutorials/tutorial_3/01-ki/zno.md', 'tutorial_3/md_excerpts/zno_w2kc.md', -4, -3)
    extract('../../../tutorials/tutorial_3/01-ki/zno.md', 'tutorial_3/md_excerpts/zno_ham.md', -3, -2)
    copy('../../../tutorials/tutorial_3/01-ki/01-koopmans-dfpt/Koopmans_DFPT_bandstructure.png', 'tutorial_3/')
