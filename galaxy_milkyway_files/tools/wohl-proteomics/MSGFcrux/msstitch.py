import argparse
import shutil
import tempfile
import os
import sys
import subprocess

import sample_pools


class MSStitch(object):
    def __init__(self):
        self.temp_wd = None
        self.outdir = os.getcwd()
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('--command', dest='command', required=True)
        self.parser.add_argument('--infile', dest='infile')
        self.parser.add_argument('--workdir', dest='workdir')
        self.parser.add_argument('--dbfile', dest='dbfile', default=None)
        self.parser.add_argument('--outdst', dest='outdst', required=True)

    def run(self):
        self.define_parse_options()
        self.parse_args()
        self.prepare_input()
        self.cmd = [self.command, '-c', self.args.command, '-i',
                    self.args.infile, '-d', self.outdir]
        self.build_command()
        self.run_command()
        self.move_output()
        self.remove_temp_workdir()

    def parse_args(self):
        self.args = self.parser.parse_args()

    def define_parse_options(self):
        raise NotImplementedError

    def prepare_input(self):
        pass

    def build_command(self, excluded=None):
        to_exclude = ['command', 'workdir', 'infile', 'outdst', 'specext']
        if excluded:
            to_exclude.extend(excluded)
        listargs = [('--{0}'.format(x), vars(self.args)[x])
                    for x in vars(self.args)
                    if x not in to_exclude
                    and vars(self.args)[x] is not None]
        self.cmd.extend([x for y in listargs for x in y])

    def move_output(self):
        outfile = self.get_msstitch_outfile_path()
        shutil.move(outfile, self.args.outdst)

    def run_command(self):
        subprocess.call(self.cmd)

    def get_msstitch_outfile_path(self):
        outfn = os.path.basename(self.args.infile) +\
            self.suffixes[self.args.command]
        return os.path.join(self.outdir, outfn)

    def copy_input_db(self):
        if self.args.dbfile is not None:
            copied_db = os.path.join(self.outdir,
                                     os.path.basename(self.args.dbfile)
                                     + '_copy')
            shutil.copy(self.args.dbfile, copied_db)
            self.args.dbfile = copied_db

    def set_tmp_workdir(self):
        self.temp_wd = tempfile.mkdtemp(dir=self.args.workdir)
        self.outdir = self.temp_wd

    def remove_temp_workdir(self):
        if self.temp_wd is not None:
            shutil.rmtree(self.temp_wd)


class MSLookup(MSStitch):
    command = 'mslookup.py'

    def prepare_input(self):
        self.set_tmp_workdir()
        self.copy_input_db()

    def set_mslookup_command_stub(self):
        self.cmd = [self.command, '-c', self.args.command,
                    '-d', self.outdir, '-i']

    def build_command(self, excluded=None):
        if excluded is None:
            excluded = []
        self.set_mslookup_command_stub()
        self.cmd.extend(self.args.infile)
        if self.args.spectra is not None:
            self.cmd.append('--spectra')
            self.cmd.extend(self.args.spectra)
            excluded.append('spectra')
        super(MSLookup, self).build_command(excluded=excluded)

    def define_parse_options(self):
        self.parser.add_argument('--spectra', dest='spectra', default=None)
        self.parser.add_argument('--setnames', dest='setnames', default=None)
        self.parser.add_argument('--spectracol', dest='spectracol',
                                 default=None)
        self.parser.add_argument('--isobquantcolpattern', dest='isobquantcolpattern',
                                 default=None)
        self.parser.add_argument('--ms1quantcolpattern', dest='ms1quantcolpattern',
                                 default=None)
        self.parser.add_argument('--psmnrcolpattern', dest='psmnrcolpattern',
                                 default=None)
        self.parser.add_argument('--probcolpattern', dest='probcolpattern',
                                 default=None)
        self.parser.add_argument('--quanttype', dest='quanttype',
                                 default=None)
        self.parser.add_argument('--protcol', dest='protcol')
        self.parser.add_argument('--fasta', dest='fasta')
        self.parser.add_argument('--confidence-col', dest='confidence-col')
        self.parser.add_argument('--confidence-lvl', dest='confidence-lvl')
        self.parser.add_argument('--confidence-better',
                                 dest='confidence-better')
        self.parser.add_argument('--mztol', dest='mztol')
        self.parser.add_argument('--mztoltype', dest='mztoltype')
        self.parser.add_argument('--rttol', dest='rttol')

    def parse_args(self):
        """Gets multifile input from dataset dir"""
        args = self.parser.parse_args()
        for argname in ['infile', 'spectra']:
            content = getattr(args, argname)
            if content is None:
                continue
            if os.path.isdir(content):
                parsedcontent = [os.path.join(content, x)
                                 for x in sorted(os.listdir(content))]
            elif os.path.isfile(content):
                parsedcontent = [content]
            setattr(args, argname, parsedcontent)

        self.args = args

    def get_msstitch_outfile_path(self):
        if self.args.dbfile is not None:
            return self.args.dbfile
        else:
            return os.path.join(self.outdir,
                                'msstitcher_lookup.sqlite')


class MSLookupBiosets(MSLookup):

    def parse_args(self):
        self.args = self.parser.parse_args()

    def build_command(self):
        files_sets = self.get_files_and_sets()
        self.args.infile = [fset[0] for fset in files_sets]
        super(MSLookupBiosets, self).build_command(['setnames'])
        self.cmd.append('--setnames')
        self.cmd.extend(['"{0}"'.format(fset[1]) for fset in files_sets])
    

class MSLookupBiosetsSpectra(MSLookupBiosets):
    command = 'mslookup.py'

    def get_files_and_sets(self):
        setmap = sample_pools.get_all_sets(self.args.setnames)
        files_sets = []
        for fn in os.listdir(self.args.infile):
            setname = sample_pools.map_set_fn(setmap, fn)
            files_sets.append((os.path.join(self.args.infile, fn), setname))
        return files_sets


class MSLookupBiosetsProtquant(MSLookupBiosets):
    command = 'mslookup.py'

    def get_files_and_sets(self):
        setmap = sample_pools.get_all_sets(self.args.setnames)
        files_sets = []
        for fn in os.listdir(self.args.infile):
            tasknr = fn.split('_')[1]
            setname = sample_pools.get_setname_from_tasknr_percolator_batch(tasknr, setmap)
            files_sets.append((os.path.join(self.args.infile, fn), setname))
        return files_sets


class MzidPlus(MSStitch):
    command = 'mzidplus.py'
    suffixes = {
        'spectratsv': '_spectradata.tsv',
        'percotsv': '_percolated.tsv',
        'quanttsv': '_quant.tsv',
        'proteingroup': '_protgroups.txt',
        'mergetsv': '_concat.tsv',
        'splittsv': '_split.tsv',
        'peptable': '_peptable.tsv',
    }

    def define_parse_options(self):
        self.parser.add_argument('--mzid', dest='mzid')
        self.parser.add_argument('--fasta', dest='fasta')
        self.parser.add_argument('--isobaric', dest='isobaric',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--spectracol', dest='spectracol',
                                 default=None)
        self.parser.add_argument('--confidence-col', dest='confidence-col')
        self.parser.add_argument('--confidence-lvl', dest='confidence-lvl')
        self.parser.add_argument('--confidence-better',
                                 dest='confidence-better')

        self.parser.add_argument('--precursor', dest='precursor',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--fncol', dest='fncol')
        self.parser.add_argument('--scorecol', dest='scorecol')

    def build_command(self, excluded=None):
        flags = ['isobaric', 'precursor']
        to_exclude = flags[:]
        if excluded is not None:
            to_exclude.extend(excluded)
        super(MzidPlus, self).build_command(excluded=to_exclude)
        for flag in flags:
            if hasattr(self.args, flag) and getattr(self.args, flag):
                self.cmd.append('--{0}'.format(flag))

    def prepare_input(self):
        self.set_tmp_workdir()
        self.copy_input_db()


class MzidPlusSplit(MzidPlus):
    def define_parse_options(self):
        self.parser.add_argument('--bioset', dest='bioset',
                                 action="store_const", const=True,
                                 default=False)
        self.parser.add_argument('--splitcol', dest='splitcol')
        self.parser.add_argument('--rename-cols', dest='rename-cols',
                                 default=None)
        self.parser.add_argument('--rename-col-startswith',
                                 dest='rename-col-startswith')

    def prepare_input(self):
        self.splitoutdir = 'mzidplussplitout'
        os.mkdir(self.splitoutdir)
        self.outdir = os.path.join(os.getcwd(), self.splitoutdir)

    def parse_args(self):
        super(MzidPlusSplit, self).parse_args()
        args = vars(self.args)
        if args['rename-cols'] is None:
            return
        parsed_colnrs = []
        cols = args['rename-cols'].split(',')
        for col in cols:
            try:
                parsed_colnrs.append(int(col) - 1)
            except ValueError:
                col = col.split('-')
                parsed_colnrs.extend([x - 1 for x in range(col[0],
                                                           col[1] + 1)])
        self.args['rename-cols'] = parsed_colnrs[:]

    def build_command(self, excluded=None):
        to_exclude = ['bioset']
        if excluded is not None:
            to_exclude.extend(excluded)
        super(MzidPlusSplit, self).build_command(excluded=to_exclude)
        if hasattr(self.args, 'bioset') and self.args.bioset:
            self.cmd.append('--bioset')

    def move_output(self):
        outfiles = os.listdir(self.splitoutdir)
        outdir = '{0}_files'.format(os.path.splitext(self.args.outdst)[0])
        os.mkdir(outdir)
        for task_nr, fn in enumerate(sorted(outfiles)):
            outfile = 'task_{0}_{1}'.format(task_nr,
                                            os.path.basename(self.args.outdst))
            infile = os.path.join(os.getcwd(), self.splitoutdir, fn)
            shutil.move(infile, os.path.join(outdir, outfile))


class MzidPlusMerge(MzidPlus):
    def define_parse_options(self):
        self.parser.add_argument('--multifiles', dest='multifiles', nargs='+')

    def parse_args(self):
        in_args = self.parser.parse_args()
        infiles = [os.path.join(in_args.infile, x) for x in
                   sorted(os.listdir(in_args.infile))]
        args_to_parse = ['--command', in_args.command,
                         '--infile', infiles[0],
                         '--multifiles', ]
        args_to_parse.extend(infiles[1:])
        args_to_parse.extend(['--outdst', in_args.outdst])
        self.args = self.parser.parse_args(args_to_parse)

    def build_command(self, excluded=None):
        to_exclude = excluded
        if to_exclude is None:
            to_exclude = []
        if 'multifiles' not in to_exclude:
            to_exclude.append('multifiles')
        super(MzidPlusMerge, self).build_command(to_exclude)
        self.cmd.append('--multifiles')
        self.cmd.extend(self.args.multifiles)


class Prottable(MSStitch):
    command = 'prottable.py'
    suffixes = {
        'addprotdata': '_proteindata.txt',
        'addms1quant': '_ms1q.txt',
        'addprob': '_protprob.txt',
        'buildquant': '',
        'createlabelfree': '_labelfree.txt',
        'qvality': '_protqvality.txt',
    }

    def define_parse_options(self):
        self.parser.add_argument('--fasta', dest='fasta')
        self.parser.add_argument('--psmtable', dest='psmtable')
        self.parser.add_argument('--peptable', dest='peptable')
        self.parser.add_argument('--decoy', dest='fasta')
        self.parser.add_argument('--proteindata', dest='proteindata',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--isobaric', dest='isobaric',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--precursor', dest='precursor',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--probability', dest='probability',
                                 action='store_const',
                                 default=False, const=True)
        self.parser.add_argument('--setname', dest='setname', default=None)

    def build_command(self, excluded=None):
        self.cmd = [self.command, '-c', self.args.command,
                    '-d', self.outdir]
        
        flags = ['proteindata', 'isobaric', 'precursor', 'probability']
        to_exclude = flags[:] + ['setname']
        if excluded is not None:
            to_exclude.extend(excluded)
        super(Prottable, self).build_command(excluded=to_exclude)
        for flag in flags:
            if hasattr(self.args, flag) and getattr(self.args, flag):
                self.cmd.append('--{0}'.format(flag))
        if hasattr(self.args, 'setname') and getattr(self.args, 'setname'):
            self.cmd.extend(['--setname', sample_pools.get_setname_from_setfile(self.args.setname)])
        if hasattr(self.args, 'infile') and self.args.infile:
            self.cmd.extend(['-i', self.args.infile])
        else:
            self.args.infile = 'built_protein_table.txt'
    
    def prepare_input(self):
        self.set_tmp_workdir()
        self.copy_input_db()


class Pycolator(MSStitch):
    command = 'pycolator.py'
    suffixes = {}


def main():
    actions = {
        'percotsv': MzidPlus,
        'peptable': MzidPlus,
        'quanttsv': MzidPlus,
        'spectratsv': MzidPlus,
        'proteingroup': MzidPlus,
        'mergetsv': MzidPlusMerge,
        'splittsv': MzidPlusSplit,
        'spectra': MSLookupBiosetsSpectra,
        'psms': MSLookup,
        'isoquant': MSLookup,
        'ms1quant': MSLookup,
        'proteingrouplookup': MSLookup,
        'protquant': MSLookupBiosetsProtquant,
        'addprotdata': Prottable,
        'addms1quant': Prottable,
        'addprob': Prottable,
        'buildquant': Prottable,
        'createlabelfree': Prottable,
        'qvality': Prottable,
    }
    app = actions[sys.argv[2]]()
    app.run()

if __name__ == '__main__':
    main()
