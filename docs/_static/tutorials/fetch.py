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
        else:
            raise ValueError(f'Could not find the end of the {heading} block')

    flines = flines[start:end]

    # Find the shortest leading whitespace and strip this from all lines (as otherwise the markdown will be rendered as a code block)
    min_indent = min(len(line) - len(line.lstrip()) for line in flines if line.strip())
    flines = [line if line == "\n" else line[min_indent:]for line in flines]

    # Manual conversion of double spaces to /
    flines = [line[:-2] + '\ \n' if (line.endswith('  \n')
                                     and not line.strip().startswith('-')) else line for line in flines]

    with open(filename_out, 'w') as fd:
        fd.writelines(flines)


if __name__ == '__main__':
    # Tutorial 1
    extract('../../../tutorials/tutorial_1/ozone.md', 'tutorial_1/md_excerpts/ozone_init.md', heading='Initialization')
    extract('../../../tutorials/tutorial_1/ozone.md', 'tutorial_1/md_excerpts/ozone_alpha.md', 24, 33)
    extract('../../../tutorials/tutorial_1/ozone.md', 'tutorial_1/md_excerpts/ozone_alpha_10.md', 45, 49)
    extract('../../../tutorials/tutorial_1/ozone.md', 'tutorial_1/md_excerpts/ozone_tables.md', 50, 62)
    extract('../../../tutorials/tutorial_1/ozone.md', 'tutorial_1/md_excerpts/ozone_final.md', -4)

    # Tutorial 2
    extract('../../../tutorials/tutorial_2/si_wannierize.md', 'tutorial_2/md_excerpts/si_wannierize.md', start=18)
    extract('../../../tutorials/tutorial_2/si_ki.md', 'tutorial_2/md_excerpts/si_ki_wannierize.md', heading="Wannierize")
    extract('../../../tutorials/tutorial_2/si_ki.md', 'tutorial_2/md_excerpts/si_ki_fold.md', heading="Fold To Supercell")
    extract('../../../tutorials/tutorial_2/si_ki.md', 'tutorial_2/md_excerpts/si_ki_screening.md', 51, 61)
    extract('../../../tutorials/tutorial_2/si_ki.md',
            'tutorial_2/md_excerpts/si_ki_postproc.md', heading="Unfold And Interpolate")

    # Tutorial 3
    extract('../../../tutorials/tutorial_3/01-ki/zno.md',
            'tutorial_3/md_excerpts/zno_wannierize_section.md', heading='Wannierize')
    extract('../../../tutorials/tutorial_3/01-ki/zno.md', 'tutorial_3/md_excerpts/zno_w2kc.md', -4, -3)
    extract('../../../tutorials/tutorial_3/01-ki/zno.md', 'tutorial_3/md_excerpts/zno_ham.md', -3, -2)
    copy('../../../tutorials/tutorial_3/01-ki/01-koopmans-dfpt/Koopmans_DFPT_bandstructure.png', 'tutorial_3/')

    # Tutorial 5
    extract('../../../tutorials/tutorial_5/01-train/h2o_train.md', 'tutorial_5/md_excerpts/train.md', 45, 92)
    extract('../../../tutorials/tutorial_5/02-predict/h2o_predict.md',
            'tutorial_5/md_excerpts/predict.md', heading='Calculate Screening Via DSCF')
    extract('../../../tutorials/tutorial_5/03-test/h2o_test.md', 'tutorial_5/md_excerpts/test.md', 45, 96)
